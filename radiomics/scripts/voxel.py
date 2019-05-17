from collections import OrderedDict
from datetime import datetime
import logging.config
import os
import threading

import SimpleITK as sitk
import six

import radiomics.featureextractor

caseLogger = logging.getLogger('radiomics.script')
_parallel_extraction_configured = False


def extractVoxel(case_idx, case, extractor, **kwargs):
  global caseLogger

  out_dir = kwargs.get('out_dir', None)
  unix_path = kwargs.get('unix_path', False)

  # Instantiate the output
  feature_vector = OrderedDict(case)

  try:
    if out_dir is None:
      out_dir = '.'
    elif not os.path.isdir(out_dir):
      caseLogger.debug('Creating output directory at %s' % out_dir)
      os.makedirs(out_dir)

    caseLogger.info('Processing case %s', case_idx)
    t = datetime.now()

    imageFilepath = case['Image']  # Required
    maskFilepath = case['Mask']  # Required
    label = case.get('Label', None)  # Optional
    if isinstance(label, six.string_types):
      label = int(label)
    label_channel = case.get('Label_channel', None)  # Optional
    if isinstance(label_channel, six.string_types):
      label_channel = int(label)

    # Extract features
    result = extractor.execute(imageFilepath, maskFilepath, label, label_channel, voxelBased=True)

    for k in result:
      if isinstance(result[k], sitk.Image):
        target = os.path.join(out_dir, 'Case-%i_%s.nrrd' % (case_idx, k))
        sitk.WriteImage(result[k], target, True)
        if unix_path and os.path.sep != '/':
          target = target.replace(os.path.sep, '/')
        feature_vector[k] = target
      else:
        feature_vector[k] = result[k]

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


def extractVoxel_parallel(args, logging_config=None, **kwargs):
  try:
    # set thread name to patient name
    threading.current_thread().name = 'case %s' % args[0]  # args[0] = case_idx

    if logging_config is not None:
      _configureParallelExtraction(logging_config)

    return extractVoxel(*args, **kwargs)
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
