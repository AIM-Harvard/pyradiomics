#!/usr/bin/env python
import argparse
import csv
from functools import partial
import logging
from multiprocessing import cpu_count, Pool
import os
import sys

from pykwalify.compat import yaml
import six.moves

import radiomics
from . import segment


scriptlogger = logging.getLogger('radiomics.script')  # holds logger for script events
logging_config = {}
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

  # Run the extraction
  try:
    _configureLogging(args)
    scriptlogger.info('Starting PyRadiomics (version: %s)', radiomics.__version__)
    results = _processInput(args)
    if results is not None:
      segment.processOutput(results, args.out, args.skip_nans, args.format, args.format_path, relative_path_start)
      scriptlogger.info('Finished extraction successfully...')
    else:
      return 1  # Feature extraction error
  except Exception:
    scriptlogger.error('Error extracting features!', exc_info=True)
    return 3  # Unknown error
  return 0  # success


def _processInput(args):
  global logging_config, relative_path_start, scriptlogger
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

  from radiomics.scripts import segment

  if num_workers > 1:  # multiple cases, parallel processing enabled
    scriptlogger.info('Input valid, starting parallel extraction from %d cases with %d workers...',
                      caseCount, num_workers)
    pool = Pool(num_workers)
    results = pool.map(partial(segment.extractSegment_parallel, parallel_config=logging_config), caseGenerator)
  elif num_workers == 1:  # single case or sequential batch processing
    scriptlogger.info('Input valid, starting sequential extraction from %d case(s)...',
                      caseCount)
    results = []
    for case in caseGenerator:
      results.append(segment.extractSegment(*case))
  else:
    # No cases defined in the batch
    scriptlogger.error('No cases to process...')
    return None
  return results


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
      raise ValueError('Cannot understand value_type %s' % value_type)

  for setting in overrides:  # setting = "setting_key:setting_value"
    if ':' not in setting:
      scriptlogger.warning('Incorrect format for override setting "%s", missing ":"', setting)
    # split into key and value
    setting_key, setting_value = setting.split(':', 2)

    # Check if it is a valid PyRadiomics Setting
    if setting_key not in settingsSchema:
      scriptlogger.warning('Did not recognize override %s, skipping...', setting_key)
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

    except Exception:
      scriptlogger.warning('Could not parse value %s for setting %s, skipping...', setting_value, setting_key)

  return setting_overrides


def _configureLogging(args):
  global scriptlogger, logging_config

  # Initialize Logging
  logLevel = getattr(logging, args.logging_level)
  rLogger = radiomics.logger
  logging_config['logLevel'] = logLevel

  # Set up optional logging to file
  if args.log_file is not None:
    rLogger.setLevel(logLevel)
    handler = logging.FileHandler(filename=args.log_file, mode='a')
    handler.setFormatter(logging.Formatter("%(levelname)s:%(name)s: %(message)s"))
    rLogger.addHandler(handler)
    logging_config['logFile'] = args.log_file

  # Set verbosity of output (stderr)
  verboseLevel = (6 - args.verbosity) * 10  # convert to python logging level
  radiomics.setVerbosity(verboseLevel)
  logging_config['verbosity'] = verboseLevel

  scriptlogger.debug('Logging initialized')
