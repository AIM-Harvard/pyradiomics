#!/usr/bin/env python

import argparse
import collections
import csv
import json
import logging
import os.path
import sys
import traceback

import six

from radiomics import featureextractor

parser = argparse.ArgumentParser()
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
parser.add_argument('--param', '-p', metavar='FILE', nargs=1, type=str, default=None,
                    help='Parameter file containing the settings to be used in extraction')
parser.add_argument('--label', '-l', metavar='N', nargs=1, default=None, type=int,
                    help='Value of label in mask to use for feature extraction')
parser.add_argument('--logging-level', metavar='LEVEL',
                    choices=['NOTSET', 'DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL'],
                    default='WARNING', help='Set capture level for logging')
parser.add_argument('--log-file', metavar='FILE', nargs='?', type=argparse.FileType('w'), default=sys.stderr,
                    help='File to append logger output to')


def main():
  args = parser.parse_args()

  # Initialize Logging
  logLevel = eval('logging.' + args.logging_level)
  rLogger = logging.getLogger('radiomics')
  rLogger.handlers = []
  rLogger.setLevel(logLevel)

  logger = logging.getLogger()
  logger.setLevel(logLevel)
  handler = logging.StreamHandler(args.log_file)
  handler.setLevel(logLevel)
  handler.setFormatter(logging.Formatter("%(levelname)s:%(name)s: %(message)s"))
  logger.addHandler(handler)

  # Initialize extractor
  try:
    if args.param is not None:
      extractor = featureextractor.RadiomicsFeaturesExtractor(args.param[0])
    else:
      extractor = featureextractor.RadiomicsFeaturesExtractor()
    logging.info('Extracting features with kwarg settings: %s\n\tImage:%s\n\tMask:%s',
                 str(extractor.kwargs), os.path.abspath(args.image), os.path.abspath(args.mask))
    featureVector = collections.OrderedDict()
    featureVector['image'] = os.path.basename(args.image)
    featureVector['mask'] = os.path.basename(args.mask)

    featureVector.update(extractor.execute(args.image, args.mask, args.label))

    if args.format == 'csv':
      writer = csv.writer(args.out, lineterminator='\n')
      writer.writerow(featureVector.keys())
      writer.writerow(featureVector.values())
    elif args.format == 'json':
      json.dump(featureVector, args.out)
      args.out.write('\n')
    else:
      for k, v in six.iteritems(featureVector):
        args.out.write('%s: %s\n' % (k, v))
  except Exception:
    logging.error('FEATURE EXTRACTION FAILED:\n%s', traceback.format_exc())

  args.out.close()
  args.log_file.close()


if __name__ == "__main__":
  main()
