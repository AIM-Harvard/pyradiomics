#!/usr/bin/env python
import argparse
import csv
from functools import partial
import logging.config
import logging.handlers
from multiprocessing import cpu_count, Manager, Pool
import os
import sys
import threading

from pykwalify.compat import yaml
import pykwalify.core
import six.moves

import radiomics
from . import segment


scriptlogger = logging.getLogger('radiomics.script')  # holds logger for script events
relative_path_start = os.getcwd()


def parse_args(custom_arguments=None):
  global relative_path_start
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
  inputGroup.add_argument('--jobs', '-j', metavar='N', type=int, default=1, choices=six.moves.range(1, cpu_count() + 1),
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
                           help='File to append output to')
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

  args = parser.parse_args(args=custom_arguments)  # Exits with code 2 if parsing fails

  # variable to hold the listener needed for processing parallel log records
  queue_listener = None

  # Run the extraction
  try:
    logging_config, queue_listener = _configureLogging(args)
    scriptlogger.info('Starting PyRadiomics (version: %s)', radiomics.__version__)
    input_tuple = _processInput(args)
    if input_tuple is not None:
      if args.validate:
        _validateCases(*input_tuple)
      else:
        results = _extractSegment(*input_tuple, logging_config=logging_config)
        segment.processOutput(results, args.out, args.skip_nans, args.format, args.format_path, relative_path_start)
        scriptlogger.info('Finished extraction successfully...')
    else:
      return 1  # Feature extraction error
  except (KeyboardInterrupt, SystemExit):
    scriptlogger.info('Cancelling Extraction')
    return -1
  except Exception:
    scriptlogger.error('Error extracting features!', exc_info=True)
    return 3  # Unknown error
  finally:
    if queue_listener is not None:
      queue_listener.stop()
  return 0  # success


def _processInput(args):
  global relative_path_start, scriptlogger
  scriptlogger.info('Processing input...')

  caseCount = 1
  num_workers = 1

  # Check if input represents a batch file
  if args.input.endswith('.csv'):
    scriptlogger.debug('Loading batch file "%s"', args.input)
    relative_path_start = os.path.dirname(args.input)
    with open(args.input, mode='r') as batchFile:
      cr = csv.DictReader(batchFile, lineterminator='\n')

      # Check if required Image and Mask columns are present
      if 'Image' not in cr.fieldnames:
        scriptlogger.error('Required column "Image" not present in input, unable to extract features...')
        return None
      if 'Mask' not in cr.fieldnames:
        scriptlogger.error('Required column "Mask" not present in input, unable to extract features...')
        return None

      cases = []
      for row_idx, row in enumerate(cr, start=2):
        if row['Image'] is None or row['Mask'] is None:
          scriptlogger.warning('Batch L%d: Missing required Image or Mask, skipping this case...', row_idx)
          continue
        imPath = row['Image']
        maPath = row['Mask']
        if not os.path.isabs(imPath):
          imPath = os.path.abspath(os.path.join(relative_path_start, imPath))
          scriptlogger.debug('Updated relative image filepath to be relative to input CSV: %s', imPath)
        if not os.path.isabs(maPath):
          maPath = os.path.abspath(os.path.join(relative_path_start, maPath))
          scriptlogger.debug('Updated relative mask filepath to be relative to input CSV: %s', maPath)
        cases.append(row)
        cases[-1]['Image'] = imPath
        cases[-1]['Mask'] = maPath

      caseCount = len(cases)
      caseGenerator = _buildGenerator(args, cases)
      num_workers = min(caseCount, args.jobs)
  elif args.mask is not None:
    caseGenerator = _buildGenerator(args, [{'Image': args.input, 'Mask': args.mask}])
  else:
    scriptlogger.error('Input is not recognized as batch, no mask specified, cannot compute result!')
    return None

  return caseGenerator, caseCount, num_workers


def _extractSegment(case_generator, case_count, num_workers, logging_config):
  if num_workers > 1:  # multiple cases, parallel processing enabled
    scriptlogger.info('Input valid, starting parallel extraction from %d cases with %d workers...',
                      case_count, num_workers)
    pool = Pool(num_workers)
    try:
      task = pool.map_async(partial(segment.extractSegment_parallel, logging_config=logging_config),
                            case_generator,
                            chunksize=min(10, case_count))
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
  elif num_workers == 1:  # single case or sequential batch processing
    scriptlogger.info('Input valid, starting sequential extraction from %d case(s)...',
                      case_count)
    results = []
    for case in case_generator:
      results.append(segment.extractSegment(*case))
  else:
    # No cases defined in the batch
    scriptlogger.error('No cases to process...')
    results = None
  return results


def _validateCases(case_generator, case_count, num_workers):
  global scriptlogger
  scriptlogger.info('Validating input for %i cases', case_count)
  errored_cases = 0
  for case_idx, case, param, setting_overrides in case_generator:
    if case_idx == 1 and param is not None:
      if not os.path.isfile(param):
        scriptlogger.error('Path for specified parameter file does not exist!')
      else:
        schemaFile, schemaFuncs = radiomics.getParameterValidationFiles()

        c = pykwalify.core.Core(source_file=param, schema_files=[schemaFile], extensions=[schemaFuncs])
        try:
          c.validate()
        except (KeyboardInterrupt, SystemExit):
          raise
        except Exception as e:
          scriptlogger.error('Parameter validation failed!\n%s' % e.message)
    scriptlogger.debug("Validating case (%i/%i): %s", case_idx, case_count, case)

    case_error = False
    if not os.path.isfile(case['Image']):
      case_error = True
      scriptlogger.error('Image path for case (%i/%i) does not exist!', case_idx, case_count)
    if not os.path.isfile(case['Mask']):
      case_error = True
      scriptlogger.error('Mask path for case (%i/%i) does not exist!', case_idx, case_count)

    if case_error:
      errored_cases += 1

  scriptlogger.info('Validation complete, errors found in %i case(s)', errored_cases)


def _buildGenerator(args, cases):
  global scriptlogger
  setting_overrides = _parseOverrides(args.setting)

  # Section for deprecated argument label
  if args.label is not None:
    scriptlogger.warning('Argument "label" is deprecated. To specify a custom label, use argument "setting" as follows:'
                         '"--setting=label:N", where N is the a label value.')
    setting_overrides['label'] = args.label
  # End deprecated section

  for case_idx, case in enumerate(cases, start=1):
    yield case_idx, case, args.param, setting_overrides


def _parseOverrides(overrides):
  global scriptlogger
  setting_overrides = {}

  # parse overrides
  if len(overrides) == 0:
    scriptlogger.debug('No overrides found')
    return setting_overrides

  scriptlogger.debug('Reading parameter schema')
  schemaFile, schemaFuncs = radiomics.getParameterValidationFiles()
  with open(schemaFile) as schema:
    settingsSchema = yaml.load(schema)['mapping']['setting']['mapping']

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

  for setting in overrides:  # setting = "setting_key:setting_value"
    if ':' not in setting:
      scriptlogger.warning('Incorrect format for override setting "%s", missing ":"', setting)
      continue
    # split into key and value
    setting_key, setting_value = setting.split(':', 2)

    # Check if it is a valid PyRadiomics Setting
    if setting_key not in settingsSchema:
      scriptlogger.warning('Did not recognize override "%s", skipping...', setting_key)
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
        scriptlogger.debug('Parsed "%s" as list (element type "%s"); value: %s',
                           setting_key, setting_type, setting_overrides[setting_key])
      else:
        if 'type' in setting_def:
          setting_type = setting_def['type']
        setting_overrides[setting_key] = parse_value(setting_value, setting_type)
        scriptlogger.debug('Parsed "%s" as type "%s"; value: %s', setting_key, setting_type, setting_overrides[setting_key])

    except (KeyboardInterrupt, SystemExit):
      raise
    except Exception:
      scriptlogger.warning('Could not parse value "%s" for setting "%s", skipping...', setting_value, setting_key)

  return setting_overrides


def _configureLogging(args):
  global scriptlogger

  # Listener to process log messages from child processes in case of multiprocessing
  queue_listener = None

  logfileLevel = getattr(logging, args.logging_level)
  verboseLevel = (6 - args.verbosity) * 10  # convert to python logging level
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

  if args.jobs > 1:
    # Update the logger format to include the threadname if multiprocessing
    # is enabled
    logging_config['formatters']['default']['format'] = \
      '[%(asctime)s] %(levelname)-.1s: (%(threadName)s) %(name)s: %(message)s'

  # Set up optional logging to file
  if args.log_file is not None:
    py_version = (sys.version_info.major, sys.version_info.minor)
    if args.jobs > 1 and py_version >= (3, 2):
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

      file_handler = logging.FileHandler(filename=args.log_file, mode='a')
      file_handler.setFormatter(logging.Formatter(fmt=logging_config['formatters']['default'].get('format'),
                                                  datefmt=logging_config['formatters']['default'].get('datefmt')))

      queue_listener = logging.handlers.QueueListener(q, file_handler)
      queue_listener.start()
    else:
      logging_config['handlers']['logfile'] = {
        'class': 'logging.FileHandler',
        'filename': args.log_file,
        'mode': 'a',
        'level': logfileLevel,
        'formatter': 'default'
      }
    logging_config['loggers']['radiomics']['handlers'].append('logfile')

  logging.config.dictConfig(logging_config)

  scriptlogger.debug('Logging initialized')
  return logging_config, queue_listener
