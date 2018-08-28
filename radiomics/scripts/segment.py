from collections import OrderedDict
import csv
from datetime import datetime
from functools import partial
import json
import logging
import os
import threading

import numpy
import SimpleITK as sitk
import six

from radiomics import featureextractor, setVerbosity

caseLogger = logging.getLogger('radiomics.script')


def extractSegment(case_idx, case, config, config_override):
  global caseLogger

  # Instantiate the output
  feature_vector = OrderedDict(case)

  try:
    t = datetime.now()

    imageFilepath = case['Image']  # Required
    maskFilepath = case['Mask']  # Required
    label = case.get('Label', None)  # Optional
    if isinstance(label, six.string_types):
      label = int(label)

    # Instantiate Radiomics Feature extractor
    extractor = featureextractor.RadiomicsFeaturesExtractor(config, **config_override)

    # Extract features
    feature_vector.update(extractor.execute(imageFilepath, maskFilepath, label))

    # Display message
    delta_t = datetime.now() - t
    caseLogger.info('Patient %s processed in %s', case_idx, delta_t)

  except Exception:
    caseLogger.error('Feature extraction failed!', exc_info=True)

  return feature_vector


def extractSegment_parallel(args, parallel_config=None):
  if parallel_config is not None:
    _configurParallelExtraction(parallel_config)
    # set thread name to patient name
    threading.current_thread().name = 'case %s' % args[0]  # args[0] = case_idx
  return extractSegment(*args)


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


def extractSegmentWithTempFiles_parallel(args, parallel_config=None):
  if parallel_config is not None:
    _configurParallelExtraction(parallel_config)
    # set thread name to patient name
    threading.current_thread().name = 'case %s' % args[0]  # args[0] = case_idx
  return extractSegmentWithTempFiles(*args)


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


def _configurParallelExtraction(parallel_config):
  """
  Initialize logging for parallel extraction. This needs to be done here, as it needs to be done for each thread that is
  created.
  """
  # Configure logging
  ###################

  rLogger = logging.getLogger('radiomics')

  # Add logging to file is specified
  logFile = parallel_config.get('logFile', None)
  if logFile is not None:
    logHandler = logging.FileHandler(filename=logFile, mode='a')
    logHandler.setLevel(parallel_config.get('logLevel', logging.INFO))
    rLogger.addHandler(logHandler)

  # Include thread name in Log-message output for all handlers.
  parallelFormatter = logging.Formatter('[%(asctime)-.19s] %(levelname)-.1s: (%(threadName)s) %(name)s: %(message)s')
  for h in rLogger.handlers:
    h.setFormatter(parallelFormatter)

  if parallel_config.get('addFilter', True):
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
    outputhandler = rLogger.handlers[0]  # Handler printing to the output
    outputhandler.addFilter(info_filter('radiomics.script'))

  # Ensures that log messages are being passed to the filter with the specified level
  setVerbosity(parallel_config.get('verbosity', logging.INFO))

  # Ensure the entire extraction for each cases is handled on 1 thread
  ####################################################################

  sitk.ProcessObject_SetGlobalDefaultNumberOfThreads(1)
