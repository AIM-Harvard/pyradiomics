from collections import OrderedDict
import csv
from datetime import datetime
import logging.config
import os
import threading

import SimpleITK as sitk
import six

import radiomics.featureextractor

caseLogger = logging.getLogger('radiomics.script')
_parallel_extraction_configured = False


def extractSegment(case_idx, case, extractor, **kwargs):
  global caseLogger

  out_dir = kwargs.get('out_dir', None)

  if out_dir is None:
    return _extractFeatures(case_idx, case, extractor)

  filename = os.path.join(out_dir, 'features_%s.csv' % case_idx)
  if os.path.isfile(filename):
    # Output already generated, load result (prevents re-extraction in case of interrupted process)
    with open(filename, 'r') as outputFile:
      reader = csv.reader(outputFile)
      headers = six.next(reader)
      values = six.next(reader)
      feature_vector = OrderedDict(zip(headers, values))

    caseLogger.info('Patient %s already processed, reading results...', case_idx)
  else:
    # Extract the set of features. Set parallel_config flag to None, as any logging initialization is already handled.
    feature_vector = _extractFeatures(case_idx, case, extractor)

    # Store results in temporary separate files to prevent write conflicts
    # This allows for the extraction to be interrupted. Upon restarting, already processed cases are found in the
    # TEMP_DIR directory and loaded instead of re-extracted
    with open(filename, 'w') as outputFile:
      writer = csv.DictWriter(outputFile, fieldnames=list(feature_vector.keys()), lineterminator='\n')
      writer.writeheader()
      writer.writerow(feature_vector)

  return feature_vector


def _extractFeatures(case_idx, case, extractor):
  global caseLogger

  # Instantiate the output
  feature_vector = OrderedDict(case)

  try:
    caseLogger.info('Processing case %s', case_idx)
    t = datetime.now()

    imageFilepath = case['Image']  # Required
    maskFilepath = case['Mask']  # Required
    label = case.get('Label', None)  # Optional
    if isinstance(label, six.string_types):
      label = int(label)
    label_channel = case.get('Label_channel', None)  # Optional
    if isinstance(label_channel, six.string_types):
      label_channel = int(label_channel)

    # Extract features
    feature_vector.update(extractor.execute(imageFilepath, maskFilepath, label, label_channel))

    # Display message
    delta_t = datetime.now() - t
    caseLogger.info('Case %s processed in %s', case_idx, delta_t)

  except (KeyboardInterrupt, SystemExit):  # Cancel extraction by forwarding this 'error'
    raise
  except SystemError:
    # Occurs when Keyboard Interrupt is caught while the thread is processing a SimpleITK call
    raise KeyboardInterrupt()
  except Exception:
    caseLogger.error('Feature extraction failed!', exc_info=True)

  return feature_vector


def extractSegment_parallel(args, logging_config=None, **kwargs):
  try:
    # set thread name to patient name
    threading.current_thread().name = 'case %s' % args[0]  # args[0] = case_idx

    if logging_config is not None:
      _configureParallelExtraction(logging_config)

    return extractSegment(*args, **kwargs)
  except (KeyboardInterrupt, SystemExit):
    # Catch the error here, as this represents the interrupt of the child process.
    # The main process is also interrupted, and cancellation is further handled there
    return None


def _configureParallelExtraction(logging_config, add_info_filter=True):
  """
  Initialize logging for parallel extraction. This needs to be done here, as it needs to be done for each thread that is
  created.
  """
  global _parallel_extraction_configured
  if _parallel_extraction_configured:
    return

  # Configure logging
  ###################

  logging.config.dictConfig(logging_config)

  if add_info_filter:
    # Define filter that allows messages from specified filter and level INFO and up, and level WARNING and up from
    # other loggers.
    class info_filter(logging.Filter):
      def __init__(self, name):
        super(info_filter, self).__init__(name)
        self.level = logging.WARNING

      def filter(self, record):
        if record.levelno >= self.level:
          return True
        if record.name == self.name and record.levelno >= logging.INFO:
          return True
        return False

    # Adding the filter to the first handler of the radiomics logger limits the info messages on the output to just
    # those from radiomics.script, but warnings and errors from the entire library are also printed to the output.
    # This does not affect the amount of logging stored in the log file.
    outputhandler = radiomics.logger.handlers[0]  # Handler printing to the output
    outputhandler.addFilter(info_filter('radiomics.script'))

  # Ensure the entire extraction for each cases is handled on 1 thread
  ####################################################################

  sitk.ProcessObject_SetGlobalDefaultNumberOfThreads(1)

  _parallel_extraction_configured = True
  radiomics.logger.debug('parallel extraction configured')
