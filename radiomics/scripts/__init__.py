#!/usr/bin/env python
import argparse
import csv
from functools import partial
import json
import logging.config
import logging.handlers
from multiprocessing import cpu_count, Manager, Pool
import os
import sys
import threading

import numpy
from pykwalify.compat import yaml
import pykwalify.core
import six.moves

import radiomics.featureextractor
from . import segment, voxel


class PyRadiomicsCommandLine:

  def __init__(self, custom_arguments=None):
    self.logger = logging.getLogger('radiomics.script')  # holds logger for script events
    self.relative_path_start = os.getcwd()
    self.args = self.getParser().parse_args(args=custom_arguments)  # Exits with code 2 if parsing fails

    self.logging_config, self.queue_listener = self._configureLogging()

    if self.args.mode == 'segment':
      self.serial_func = segment.extractSegment
      self.parallel_func = segment.extractSegment_parallel
    else:
      self.serial_func = voxel.extractVoxel
      self.parallel_func = voxel.extractVoxel_parallel

    self.case_count = 0
    self.num_workers = 0

  @classmethod
  def getParser(cls):
    parser = argparse.ArgumentParser(usage='%(prog)s image|batch [mask] [Options]',
                                     formatter_class=argparse.RawTextHelpFormatter)

    inputGroup = parser.add_argument_group(title='Input',
                                           description='Input files and arguments defining the extraction:\n'
                                                       '- image and mask files (single mode) '
                                                       'or CSV-file specifying them (batch mode)\n'
                                                       '- Parameter file (.yml/.yaml or .json)\n'
                                                       '- Overrides for customization type 3 ("settings")\n'
                                                       '- Multi-threaded batch processing')
    inputGroup.add_argument('input', metavar='{Image,Batch}FILE',
                            help='Image file (single mode) or CSV batch file (batch mode)')
    inputGroup.add_argument('mask', nargs='?', metavar='MaskFILE', default=None,
                            help='Mask file identifying the ROI in the Image. \n'
                                 'Only required when in single mode, ignored otherwise.')
    inputGroup.add_argument('--param', '-p', metavar='FILE', default=None,
                            help='Parameter file containing the settings to be used in extraction')
    inputGroup.add_argument('--setting', '-s', metavar='"SETTING_NAME:VALUE"', action='append', default=[], type=str,
                            help='Additional parameters which will override those in the\n'
                                 'parameter file and/or the default settings. Multiple\n'
                                 'settings possible. N.B. Only works for customization\n'
                                 'type 3 ("setting").')
    inputGroup.add_argument('--jobs', '-j', metavar='N', type=int, default=1,
                            choices=six.moves.range(1, cpu_count() + 1),
                            help='(Batch mode only) Specifies the number of threads to use for\n'
                                 'parallel processing. This is applied at the case level;\n'
                                 'i.e. 1 thread per case. Actual number of workers used is\n'
                                 'min(cases, jobs).')
    inputGroup.add_argument('--validate', action='store_true',
                            help='If specified, check if input is valid and check if file locations point to exisiting '
                                 'files')

    outputGroup = parser.add_argument_group(title='Output', description='Arguments controlling output redirection and '
                                                                        'the formatting of calculated results.')
    outputGroup.add_argument('--out', '-o', metavar='FILE', type=argparse.FileType('a'), default=sys.stdout,
                             help='File to append output to.')
    outputGroup.add_argument('--out-dir', '-od', type=str, default=None,
                             help='Directory to store output. If specified in segment mode, this writes csv output for '
                                  'each processed case. In voxel mode, this directory is used to store the featuremaps.'
                                  ' If not specified in voxel mode, the current working directory is used instead.')
    outputGroup.add_argument('--mode', '-m', choices=['segment', 'voxel'], default='segment',
                             help='Extraction mode for PyRadiomics.')
    outputGroup.add_argument('--skip-nans', action='store_true',
                             help='Add this argument to skip returning features that have an\n'
                                  'invalid result (NaN)')
    outputGroup.add_argument('--format', '-f', choices=['csv', 'json', 'txt'], default='txt',
                             help='Format for the output.\n'
                                  '"txt" (Default): one feature per line in format "case-N_name:value"\n'
                                  '"json": Features are written in a JSON format dictionary\n'
                                  '(1 dictionary per case, 1 case per line) "{name:value}"\n'
                                  '"csv": one row of feature names, followed by one row of\n'
                                  'feature values per case.')
    outputGroup.add_argument('--format-path', choices=['absolute', 'relative', 'basename'], default='absolute',
                             help='Controls input image and mask path formatting in the output.\n'
                                  '"absolute" (Default): Absolute file paths.\n'
                                  '"relative": File paths relative to current working directory.\n'
                                  '"basename": Only stores filename.')
    outputGroup.add_argument('--unix-path', '-up', action='store_true',
                             help='If specified, ensures that all paths in the output\n'
                                  'use unix-style path separators ("/").')

    loggingGroup = parser.add_argument_group(title='Logging',
                                             description='Controls the (amount of) logging output to the '
                                                         'console and the (optional) log-file.')
    loggingGroup.add_argument('--logging-level', metavar='LEVEL',
                              choices=['NOTSET', 'DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL'],
                              default='WARNING', help='Set capture level for logging')
    loggingGroup.add_argument('--log-file', metavar='FILE', default=None, help='File to append logger output to')
    loggingGroup.add_argument('--verbosity', '-v', action='store', nargs='?', default=3, const=4, type=int,
                              choices=[1, 2, 3, 4, 5],
                              help='Regulate output to stderr. By default [3], level\n'
                                   'WARNING and up are printed. By specifying this\n'
                                   'argument without a value, level INFO [4] is assumed.\n'
                                   'A higher value results in more verbose output.')
    parser.add_argument('--label', '-l', metavar='N', default=None, type=int,
                        help='(DEPRECATED) Value of label in mask to use for\n'
                             'feature extraction.')

    parser.add_argument('--version', action='version', help='Print version and exit',
                        version='%(prog)s ' + radiomics.__version__)
    return parser

  def run(self):
    # Run the extraction
    try:
      self.logger.info('Starting PyRadiomics (version: %s)', radiomics.__version__)
      caseGenerator = self._processInput()
      if caseGenerator is not None:
        if self.args.validate:
          self._validateCases(caseGenerator)
        else:
          results = self._processCases(caseGenerator)
          self._processOutput(results)
          self.logger.info('Finished %s-based extraction successfully...', self.args.mode)
      else:
        return 1  # Feature extraction error
    except (KeyboardInterrupt, SystemExit):
      self.logger.info('Cancelling Extraction')
      return -1
    except Exception:
      self.logger.error('Error extracting features!', exc_info=True)
      return 3  # Unknown error
    finally:
      if self.queue_listener is not None:
        self.queue_listener.stop()
    return 0  # success

  def _processInput(self):
    self.logger.info('Processing input...')

    self.case_count = 1
    self.num_workers = 1

    # Check if input represents a batch file
    if self.args.input.endswith('.csv'):
      self.logger.debug('Loading batch file "%s"', self.args.input)
      self.relative_path_start = os.path.dirname(self.args.input)
      with open(self.args.input, mode='r') as batchFile:
        cr = csv.DictReader(batchFile, lineterminator='\n')

        # Check if required Image and Mask columns are present
        if 'Image' not in cr.fieldnames:
          self.logger.error('Required column "Image" not present in input, unable to extract features...')
          return None
        if 'Mask' not in cr.fieldnames:
          self.logger.error('Required column "Mask" not present in input, unable to extract features...')
          return None

        cases = []
        for row_idx, row in enumerate(cr, start=2):
          if row['Image'] is None or row['Mask'] is None:
            self.logger.warning('Batch L%d: Missing required Image or Mask, skipping this case...', row_idx)
            continue
          imPath = row['Image']
          maPath = row['Mask']
          if not os.path.isabs(imPath):
            imPath = os.path.abspath(os.path.join(self.relative_path_start, imPath))
            self.logger.debug('Considered image filepath to be relative to input CSV. Absolute path: %s', imPath)
          if not os.path.isabs(maPath):
            maPath = os.path.abspath(os.path.join(self.relative_path_start, maPath))
            self.logger.debug('Considered mask filepath to be relative to input CSV. Absolute path: %s', maPath)
          cases.append(row)
          cases[-1]['Image'] = imPath
          cases[-1]['Mask'] = maPath

          self.case_count = len(cases)
        caseGenerator = enumerate(cases, start=1)
        self.num_workers = min(self.case_count, self.args.jobs)
    elif self.args.mask is not None:
      caseGenerator = [(1, {'Image': self.args.input, 'Mask': self.args.mask})]
    else:
      self.logger.error('Input is not recognized as batch, no mask specified, cannot compute result!')
      return None

    return caseGenerator

  def _validateCases(self, case_generator):
    self.logger.info('Validating input for %i cases', self.case_count)
    errored_cases = 0
    for case_idx, case in case_generator:
      if case_idx == 1 and self.args.param is not None:
        if not os.path.isfile(self.args.param):
          self.logger.error('Path for specified parameter file does not exist!')
        else:
          schemaFile, schemaFuncs = radiomics.getParameterValidationFiles()

          c = pykwalify.core.Core(source_file=self.args.param, schema_files=[schemaFile], extensions=[schemaFuncs])
          try:
            c.validate()
          except (KeyboardInterrupt, SystemExit):
            raise
          except Exception:
            self.logger.error('Parameter validation failed!', exc_info=True)
            self.logger.debug("Validating case (%i/%i): %s", case_idx, self.case_count, case)

      case_error = False
      if not os.path.isfile(case['Image']):
        case_error = True
        self.logger.error('Image path for case (%i/%i) does not exist!', case_idx, self.case_count)
      if not os.path.isfile(case['Mask']):
        case_error = True
        self.logger.error('Mask path for case (%i/%i) does not exist!', case_idx, self.case_count)

      if case_error:
        errored_cases += 1

    self.logger.info('Validation complete, errors found in %i case(s)', errored_cases)

  def _processCases(self, case_generator):
    setting_overrides = self._parseOverrides()

    extractor = radiomics.featureextractor.RadiomicsFeatureExtractor(self.args.param, **setting_overrides)

    if self.args.out_dir is not None and not os.path.isdir(self.args.out_dir):
      os.makedirs(self.args.out_dir)

    if self.num_workers > 1:  # multiple cases, parallel processing enabled
      self.logger.info('Input valid, starting parallel extraction from %d cases with %d workers...',
                       self.case_count, self.num_workers)
      pool = Pool(self.num_workers)
      try:
        task = pool.map_async(partial(self.parallel_func,
                                      extractor=extractor,
                                      out_dir=self.args.out_dir,
                                      logging_config=self.logging_config,
                                      unix_path=self.args.unix_path),
                              case_generator,
                              chunksize=min(10, self.case_count))
        # Wait for the results to be done. task.get() without timeout performs a blocking call, which prevents
        # the program from processing the KeyboardInterrupt if it occurs
        while not task.ready():
          pass
        results = task.get()
        pool.close()
      except (KeyboardInterrupt, SystemExit):
        pool.terminate()
        raise
      finally:
        pool.join()
    elif self.num_workers == 1:  # single case or sequential batch processing
      self.logger.info('Input valid, starting sequential extraction from %d case(s)...',
                       self.case_count)
      results = []
      for case in case_generator:
        results.append(self.serial_func(*case,
                                        extractor=extractor,
                                        out_dir=self.args.out_dir,
                                        unix_path=self.args.unix_path))
    else:
      # No cases defined in the batch
      self.logger.error('No cases to process...')
      results = None
    return results

  def _processOutput(self, results):
    self.logger.info('Processing results...')

    # Store the header of all calculated features
    # By checking all headers of cases > 1, and subtracting those already in case 1, original ordering of case 1 is
    # preserved. Additional headers are by definition generated by pyradiomics, and therefore can be appended at
    # the end.
    additional_headers = set()
    for case in results[1:]:
      additional_headers.update(set(case.keys()))
      additional_headers -= set(results[0].keys())  # Subtract all headers found in the first case

    headers = list(results[0].keys()) + sorted(additional_headers)

    # Set the formatting rule for image and mask paths
    if self.args.format_path == 'absolute':
      pathFormatter = os.path.abspath
    elif self.args.format_path == 'relative':
      pathFormatter = partial(os.path.relpath, start=self.relative_path_start)
    elif self.args.format_path == 'basename':
      pathFormatter = os.path.basename
    else:
      self.logger.warning('Unrecognized format for paths (%s), reverting to default ("absolute")',
                          self.args.format_path)
      pathFormatter = os.path.abspath

    for case_idx, case in enumerate(results, start=1):
      # if specified, skip NaN values
      if self.args.skip_nans:
        for key in list(case.keys()):
          if isinstance(case[key], float) and numpy.isnan(case[key]):
            self.logger.debug('Case %d, feature %s computed NaN, removing from results', case_idx, key)
            del case[key]

      # Format paths of image and mask files
      case['Image'] = pathFormatter(case['Image'])
      case['Mask'] = pathFormatter(case['Mask'])

      if self.args.unix_path and os.path.sep != '/':
        case['Image'] = case['Image'].replace(os.path.sep, '/')
        case['Mask'] = case['Mask'].replace(os.path.sep, '/')

      # Write out results per case if format is 'csv' or 'txt', handle 'json' outside of this loop (issue #483)
      if self.args.format == 'csv':
        writer = csv.DictWriter(self.args.out, headers, lineterminator='\n', extrasaction='ignore')
        if case_idx == 1:
          writer.writeheader()
        writer.writerow(case)  # if skip_nans is enabled, nan-values are written as empty strings
      elif self.args.format == 'txt':
        for k, v in six.iteritems(case):
          self.args.out.write('Case-%d_%s: %s\n' % (case_idx, k, v))

    # JSON dump of cases is handled outside of the loop, otherwise the resultant document would be invalid.
    if self.args.format == 'json':
      # JSON cannot serialize numpy arrays, even when that array represents a scalar value (PyRadiomics Feature values)
      # Therefore, use this encoder, which first casts numpy arrays to python lists, which are JSON serializable
      class NumpyEncoder(json.JSONEncoder):
        def default(self, obj):
          if isinstance(obj, numpy.ndarray):
            return obj.tolist()
          return json.JSONEncoder.default(self, obj)

      json.dump(results, self.args.out, cls=NumpyEncoder, indent=2)

  def _parseOverrides(self):
    setting_overrides = {}

    # parse overrides
    if len(self.args.setting) == 0:
      self.logger.debug('No overrides found')
      return setting_overrides

    self.logger.debug('Reading parameter schema')
    schemaFile, schemaFuncs = radiomics.getParameterValidationFiles()
    with open(schemaFile) as schema:
      settingsSchema = yaml.safe_load(schema)['mapping']['setting']['mapping']

    # parse single value function
    def parse_value(value, value_type):
      if value_type == 'str':
        return value  # no conversion
      elif value_type == 'int':
        return int(value)
      elif value_type == 'float':
        return float(value)
      elif value_type == 'bool':
        return value == '1' or value.lower() == 'true'
      else:
        raise ValueError('Cannot understand value_type "%s"' % value_type)

    for setting in self.args.setting:  # setting = "setting_key:setting_value"
      if ':' not in setting:
        self.logger.warning('Incorrect format for override setting "%s", missing ":"', setting)
        continue
      # split into key and value
      setting_key, setting_value = setting.split(':', 2)

      # Check if it is a valid PyRadiomics Setting
      if setting_key not in settingsSchema:
        self.logger.warning('Did not recognize override "%s", skipping...', setting_key)
        continue

      # Try to parse the value by looking up its type in the settingsSchema
      try:
        setting_def = settingsSchema[setting_key]
        setting_type = 'str'  # If type is omitted in the schema, treat it as string (no conversion)
        if 'seq' in setting_def:
          # Multivalued setting
          if len(setting_def['seq']) > 0 and 'type' in setting_def['seq'][0]:
            setting_type = setting_def['seq'][0]['type']

          setting_overrides[setting_key] = [parse_value(val, setting_type) for val in setting_value.split(',')]
          self.logger.debug('Parsed "%s" as list (element type "%s"); value: %s',
                            setting_key, setting_type, setting_overrides[setting_key])
        else:
          if 'type' in setting_def:
            setting_type = setting_def['type']
          setting_overrides[setting_key] = parse_value(setting_value, setting_type)
          self.logger.debug('Parsed "%s" as type "%s"; value: %s', setting_key, setting_type,
                            setting_overrides[setting_key])

      except (KeyboardInterrupt, SystemExit):
        raise
      except Exception:
        self.logger.warning('Could not parse value "%s" for setting "%s", skipping...', setting_value, setting_key)

    # Section for deprecated argument label
    if self.args.label is not None:
      self.logger.warning(
        'Argument "label" is deprecated. To specify a custom label, use argument "setting" as follows:'
        '"--setting=label:N", where N is the a label value.')
      setting_overrides['label'] = self.args.label
    # End deprecated section

    return setting_overrides

  def _configureLogging(self):
    # Listener to process log messages from child processes in case of multiprocessing
    queue_listener = None

    logfileLevel = getattr(logging, self.args.logging_level)
    verboseLevel = (6 - self.args.verbosity) * 10  # convert to python logging level
    logger_level = min(logfileLevel, verboseLevel)

    logging_config = {
      'version': 1,
      'disable_existing_loggers': False,
      'formatters': {
        'default': {
          'format': '[%(asctime)s] %(levelname)-.1s: %(name)s: %(message)s',
          'datefmt': '%Y-%m-%d %H:%M:%S'
        }
      },
      'handlers': {
        'console': {
          'class': 'logging.StreamHandler',
          'level': verboseLevel,
          'formatter': 'default'
        }
      },
      'loggers': {
        'radiomics': {
          'level': logger_level,
          'handlers': ['console']
        }
      }
    }

    if self.args.jobs > 1:
      # Update the logger format to include the threadname if multiprocessing
      # is enabled
      logging_config['formatters']['default']['format'] = \
        '[%(asctime)s] %(levelname)-.1s: (%(threadName)s) %(name)s: %(message)s'

    # Set up optional logging to file
    if self.args.log_file is not None:
      py_version = (sys.version_info.major, sys.version_info.minor)
      if self.args.jobs > 1 and py_version >= (3, 2):
        # Multiprocessing! Use a QueueHandler, FileHandler and QueueListener
        # to implement thread-safe logging.

        # However, QueueHandler and Listener were added in python 3.2.
        # Therefore, only use this if the python version > 3.2
        q = Manager().Queue(-1)
        threading.current_thread().setName('Main')

        logging_config['handlers']['logfile'] = {
          'class': 'logging.handlers.QueueHandler',
          'queue': q,
          'level': logfileLevel,
          'formatter': 'default'
        }

        file_handler = logging.FileHandler(filename=self.args.log_file, mode='a')
        file_handler.setFormatter(logging.Formatter(fmt=logging_config['formatters']['default'].get('format'),
                                                    datefmt=logging_config['formatters']['default'].get('datefmt')))

        queue_listener = logging.handlers.QueueListener(q, file_handler)
        queue_listener.start()
      else:
        logging_config['handlers']['logfile'] = {
          'class': 'logging.FileHandler',
          'filename': self.args.log_file,
          'mode': 'a',
          'level': logfileLevel,
          'formatter': 'default'
        }
      logging_config['loggers']['radiomics']['handlers'].append('logfile')

    logging.config.dictConfig(logging_config)

    self.logger.debug('Logging initialized')
    return logging_config, queue_listener


def parse_args():
  try:
    return PyRadiomicsCommandLine().run()
  except Exception as e:
    logging.getLogger().error("Error executing PyRadiomics command line!", exc_info=True)
    print("Error executing PyRadiomics command line!\n%s" % e)
    return 4
