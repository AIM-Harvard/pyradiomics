#!/usr/bin/env python

import argparse
import sys
import os.path
import logging
import collections
import traceback
import csv
import json
from radiomics import featureextractor

parser = argparse.ArgumentParser()
parser.add_argument('inFile', metavar='In', type=argparse.FileType('r'),
                    help='CSV file containing combinations of image and mask. Each row represents one combination with '
                         'the following elements: (1) Patient Name, (2) Image type, (3) Reader, (4) Image location and '
                         '(5) Mask location')
parser.add_argument('outFile', metavar='Out', type=argparse.FileType('w'), help='File to write results to')
parser.add_argument('--format', '-f', choices=['csv', 'json'], default='csv', help='Format for the output. '
                    'Default is "csv": one row of feature names, followed by one row of feature values for each '
                    'image-mask combination. For "json": Features are written in a JSON format dictionary '
                    '"{name:value}", one line per image-mask combination')
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

  # Load patient list
  flists = []
  try:
    cr = csv.reader(args.inFile, lineterminator='\n')
    flists = [row for row in cr]
    args.inFile.close()
  except Exception:
    logging.error('CSV READ FAILED:\n%s', traceback.format_exc())
    args.inFile.close()
    args.outFile.close()
    args.log_file.close()
    exit(-1)

  # Initialize extractor
  try:
    if args.param is not None:
      extractor = featureextractor.RadiomicsFeaturesExtractor(args.param[0])
    else:
      extractor = featureextractor.RadiomicsFeaturesExtractor()
  except Exception:
    logging.error('EXTRACTOR INITIALIZATION FAILED:\n%s', traceback.format_exc())
    args.outFile.close()
    args.log_file.close()
    exit(-1)

  # Extract features
  logging.info('Extracting features with kwarg settings: %s', str(extractor.kwargs))

  headers = False
  for idx, entry in enumerate(flists, start=1):

    logging.info("(%d/%d) Processing Patient: %s, Study: %s, Reader: %s", idx, len(flists), entry[0], entry[1],
                 entry[2])

    imageFilepath = entry[3]
    maskFilepath = entry[4]

    if (imageFilepath is not None) and (maskFilepath is not None):
      featureVector = collections.OrderedDict()
      featureVector['PatientID'] = entry[0]
      featureVector['Study'] = entry[1]
      featureVector['Reader'] = entry[2]
      featureVector['image'] = os.path.basename(imageFilepath)
      featureVector['mask'] = os.path.basename(maskFilepath)

      try:
        featureVector.update(extractor.execute(imageFilepath, maskFilepath, args.label))

        if args.format == 'csv':
          writer = csv.writer(args.outFile, lineterminator='\n')
          if not headers:
            writer.writerow(featureVector.keys())
            headers = True
          writer.writerow(featureVector.values())
        elif args.format == 'json':
          json.dump(featureVector, args.out)
          args.out.write('\n')
      except Exception:
        logging.error('FEATURE EXTRACTION FAILED:\n%s', traceback.format_exc())

  args.outFile.close()
  args.log_file.close()


if __name__ == "__main__":
  main()
