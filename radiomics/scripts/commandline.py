#!/usr/bin/env python

import argparse
import collections
import csv
import json
import logging
import os.path
import sys

import numpy
import six

import radiomics
from radiomics import featureextractor

parser = argparse.ArgumentParser(usage='%(prog)s image mask [Options]')
parser.add_argument('image', metavar='Image',
                    help='Features are extracted from the Region Of Interest (ROI) in the image')
parser.add_argument('mask', metavar='Mask', help='Mask identifying the ROI in the Image')

parser.add_argument('--out', '-o', metavar='FILE', nargs='?', type=argparse.FileType('w'), default=sys.stdout,
                    help='File to append output to')
parser.add_argument('--format', '-f', choices=['txt', 'csv', 'json'], default='txt',
                    help='Format for the output. '
                         'Default is "txt": one feature per line in format "name:value". For "csv": one row of feature '
                         'names, followed by one row of feature values. For "json": Features are written in a JSON '
                         'format dictionary "{name:value}"')
parser.add_argument('--skip-nans', action='store_true',
                    help='Add this argument to skip returning features that have an invalid result (NaN)')
parser.add_argument('--param', '-p', metavar='FILE', type=str, default=None,
                    help='Parameter file containing the settings to be used in extraction')
parser.add_argument('--label', '-l', metavar='N', default=None, type=int,
                    help='Value of label in mask to use for feature extraction')
parser.add_argument('--logging-level', metavar='LEVEL',
                    choices=['NOTSET', 'DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL'],
                    default='WARNING', help='Set capture level for logging')
parser.add_argument('--log-file', metavar='FILE', type=argparse.FileType('w'), default=None,
                    help='File to append logger output to')
parser.add_argument('--verbosity', '-v', action='store', nargs='?', default=3, const=4, type=int, choices=range(0, 6),
                    help='Regulate output to stderr. By default [3], level WARNING and up are printed. By specifying '
                    'this argument without a value, level INFO [4] is assumed. A higher value results in more verbose '
                    'output.')
parser.add_argument('--shorten-path', dest='shorten', action='store_true',
                    help='specify this argument to image and mask path to just the file names')
parser.add_argument('--version', action='version', help='Print version and exit',
                    version='%(prog)s ' + radiomics.__version__)


def main():
  args = parser.parse_args()

  # Initialize Logging
  logLevel = getattr(logging, args.logging_level)
  rLogger = radiomics.logger

  # Set up optional logging to file
  if args.log_file is not None:
    rLogger.setLevel(logLevel)
    handler = logging.StreamHandler(args.log_file)
    handler.setFormatter(logging.Formatter("%(levelname)s:%(name)s: %(message)s"))
    rLogger.addHandler(handler)

  # Set verbosity of output (stderr)
  verboseLevel = (6 - args.verbosity) * 10  # convert to python logging level
  radiomics.setVerbosity(verboseLevel)

  # Initialize logging for script log messages
  logger = rLogger.getChild('script')

  # Initialize extractor
  try:
    if args.param is not None:
      extractor = featureextractor.RadiomicsFeaturesExtractor(args.param)
    else:
      extractor = featureextractor.RadiomicsFeaturesExtractor()
    logger.info('Extracting features with kwarg settings: %s\n\tImage: %s\n\tMask: %s',
                str(extractor.settings), os.path.abspath(args.image), os.path.abspath(args.mask))
    featureVector = collections.OrderedDict()
    if args.shorten:
      featureVector['Image'] = os.path.basename(args.image)
      featureVector['Mask'] = os.path.basename(args.mask)

    featureVector.update(extractor.execute(args.image, args.mask, args.label))

    # if specified, skip NaN values
    if args.skip_nans:
      for key in list(featureVector.keys()):
        if isinstance(featureVector[key], float) and numpy.isnan(featureVector[key]):
          logger.debug('Feature %s computed NaN, removing from results', key)
          del featureVector[key]

    if args.format == 'csv':
      writer = csv.writer(args.out, lineterminator='\n')
      writer.writerow(list(featureVector.keys()))
      writer.writerow(list(featureVector.values()))
    elif args.format == 'json':
      json.dump(featureVector, args.out)
      args.out.write('\n')
    else:
      for k, v in six.iteritems(featureVector):
        args.out.write('%s: %s\n' % (k, v))
  except Exception:
    logger.error('FEATURE EXTRACTION FAILED', exc_info=True)

  args.out.close()
  if args.log_file is not None:
    args.log_file.close()


if __name__ == "__main__":
  main()
