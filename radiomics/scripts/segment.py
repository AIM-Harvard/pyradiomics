from collections import OrderedDict
import csv
from datetime import datetime
from functools import partial
import json
import logging.config
import os
import threading

import numpy
import SimpleITK as sitk
import six

import radiomics.featureextractor

caseLogger = logging.getLogger('radiomics.script')
_parallel_extraction_configured = False


def extractSegment(case_idx, case, config, config_override):
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

    # Instantiate Radiomics Feature extractor
    extractor = radiomics.featureextractor.RadiomicsFeaturesExtractor(config, **config_override)

    # Extract features
    feature_vector.update(extractor.execute(imageFilepath, maskFilepath, label))

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


def extractSegment_parallel(args, logging_config=None):
  try:
    if logging_config is not None:
      # set thread name to patient name
      threading.current_thread().name = 'case %s' % args[0]  # args[0] = case_idx
      _configureParallelExtraction(logging_config)
    return extractSegment(*args)
  except (KeyboardInterrupt, SystemExit):
    # Catch the error here, as this represents the interrupt of the child process.
    # The main process is also interrupted, and cancellation is further handled there
    return None


def extractSegmentWithTempFiles(case_idx, case, config, config_override, temp_dir):
  global caseLogger

  filename = os.path.join(temp_dir, 'features_%s.csv' % case_idx)
  if os.path.isfile(filename):
    # Output already generated, load result (prevents re-extraction in case of interrupted process)
    with open(filename, 'w') as outputFile:
      reader = csv.reader(outputFile)
      headers = reader.rows[0]
      values = reader.rows[1]
      feature_vector = OrderedDict(zip(headers, values))

    caseLogger.info('Patient %s already processed, reading results...', case_idx)
  else:
    # Extract the set of features. Set parallel_config flag to None, as any logging initialization is already handled.
    feature_vector = extractSegment(case_idx, case, config, config_override)

    # Store results in temporary separate files to prevent write conflicts
    # This allows for the extraction to be interrupted. Upon restarting, already processed cases are found in the
    # TEMP_DIR directory and loaded instead of re-extracted
    with open(filename, 'w') as outputFile:
      writer = csv.DictWriter(outputFile, fieldnames=list(feature_vector.keys()), lineterminator='\n')
      writer.writeheader()
      writer.writerow(feature_vector)

  return feature_vector


def extractSegmentWithTempFiles_parallel(args, logging_config=None):
  try:
    if logging_config is not None:
      # set thread name to patient name
      threading.current_thread().name = 'case %s' % args[0]  # args[0] = case_idx
      _configureParallelExtraction(logging_config)
    return extractSegmentWithTempFiles(*args)
  except (KeyboardInterrupt, SystemExit):
    # Catch the error here, as this represents the interrupt of the child process.
    # The main process is also interrupted, and cancellation is further handled there
    return None


def processOutput(results,
                  outStream,
                  skip_nans=False,
                  format_output='csv',
                  format_path='absolute',
                  relative_path_start=''):
  global caseLogger
  caseLogger.info('Processing results...')

  # Store the header of all calculated features
  headers = results[0].keys()

  # Set the formatting rule for image and mask paths
  if format_path == 'absolute':
    pathFormatter = os.path.abspath
  elif format_path == 'relative':
    pathFormatter = partial(os.path.relpath, start=relative_path_start)
  elif format_path == 'basename':
    pathFormatter = os.path.basename
  else:
    caseLogger.warning('Unrecognized format for paths (%s), reverting to default ("absolute")', format_path)
    pathFormatter = os.path.abspath

  for case_idx, case in enumerate(results, start=1):
    # if specified, skip NaN values
    if skip_nans:
      for key in list(case.keys()):
        if isinstance(case[key], float) and numpy.isnan(case[key]):
          caseLogger.debug('Case %d, feature %s computed NaN, removing from results', case_idx, key)
          del case[key]

    # Format paths of image and mask files
    case['Image'] = pathFormatter(case['Image'])
    case['Mask'] = pathFormatter(case['Mask'])

    # Write out results
    if format_output not in ('csv', 'json', 'txt'):
      caseLogger.warning('Unrecognized format for output (%s), reverting to default ("csv")', format_output)
      format_output = 'csv'

    if format_output == 'csv':
      writer = csv.DictWriter(outStream, headers, lineterminator='\n')
      if case_idx == 1:
        writer.writeheader()
      writer.writerow(case)  # if skip_nans is enabled, nan-values are written as empty strings
    elif format_output == 'json':
      json.dump(case, outStream)
      outStream.write('\n')
    else:  # txt
      for k, v in six.iteritems(case):
        outStream.write('Case-%d_%s: %s\n' % (case_idx, k, v))


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
