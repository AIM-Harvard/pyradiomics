#!/usr/bin/env python

import argparse
import logging
import os.path

import pandas

import radiomics
from radiomics import featureextractor

parser = argparse.ArgumentParser(usage='%(prog)s In Out [Options]')
parser.add_argument('inFile', metavar='In', type=argparse.FileType('r'),
                    help='CSV file containing combinations of image and mask. First row should contain the headers, '
                         'where "Image" and "Mask" must be present and identify the image and mask locations, '
                         'respectively. All columns present in CSV file are copied to the output, this enables '
                         'specification of additional identifiers, such as patient ID, Study ID and Reader.')
parser.add_argument('outFile', metavar='Out', type=argparse.FileType('w'),
                    help='File to write results to')
parser.add_argument('--format', '-f', choices=['csv', 'json'], default='csv', help='Format for the output. '
                    'Default is "csv": one row of feature names, followed by one row of feature values for each '
                    'image-mask combination. For "json": Features are written in a JSON format dictionary '
                    '"{name:value}", one line per image-mask combination')
parser.add_argument('--param', '-p', metavar='FILE', type=str, default=None,
                    help='Parameter file containing the settings to be used in extraction')
parser.add_argument('--label', '-l', metavar='N', default=None, type=int,
                    help='Value of label in mask to use for feature extraction')
parser.add_argument('--logging-level', metavar='LEVEL',
                    choices=['NOTSET', 'DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL'],
                    default='WARNING', help='Set capture level for logging')
parser.add_argument('--log-file', metavar='FILE', type=argparse.FileType('w'), default=None,
                    help='File to append logger output to')
parser.add_argument('--shorten-path', dest='shorten', action='store_true',
                    help='specify this argument to image and mask path to just the file names')
parser.add_argument('--verbosity', '-v', action='store', nargs='?', default=3, const=4, type=int, choices=range(0, 6),
                    help='Regulate output to stderr. By default [3], level WARNING and up are printed. By specifying '
                    'this argument without a value, level INFO [4] is assumed. A higher value results in more verbose '
                    'output.')
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
  logger = rLogger.getChild('batch')

  # Load patient list
  try:
    logger.info('Loading patients')
    flists = pandas.read_csv(args.inFile).T  # Transpose data so that each column represents one test case
    logger.info('Loading succesful, found %d patients', len(flists))
  except Exception:
    logging.error('CSV READ FAILED', exc_info=True)
    args.inFile.close()
    args.outFile.close()
    if args.log_file is not None:
      args.log_file.close()
    exit(-1)

  # Initialize extractor
  try:
    logger.debug("Initializing extractor")
    if args.param is not None:
      extractor = featureextractor.RadiomicsFeaturesExtractor(args.param)
    else:
      extractor = featureextractor.RadiomicsFeaturesExtractor()
  except Exception:
    logger.error('EXTRACTOR INITIALIZATION FAILED', exc_info=True)
    args.outFile.close()
    args.log_file.close()
    exit(-1)

  # Extract features
  logger.info('Extracting features with kwarg settings: %s', str(extractor.kwargs))

  results = pandas.DataFrame()
  for entry in flists:

    logger.info("(%d/%d) Processing Patient (Image: %s, Mask: %s)",
                entry + 1,
                len(flists),
                flists[entry]['Image'],
                flists[entry]['Mask'])

    imageFilepath = flists[entry]['Image']
    maskFilepath = flists[entry]['Mask']

    if (imageFilepath is not None) and (maskFilepath is not None):
      featureVector = flists[entry]
      if args.shorten:
        featureVector['Image'] = os.path.basename(imageFilepath)
        featureVector['Mask'] = os.path.basename(maskFilepath)

      try:
        featureVector = featureVector.append(extractor.execute(imageFilepath, maskFilepath, args.label))
      except Exception:
        logger.error('FEATURE EXTRACTION FAILED', exc_info=True)

      featureVector.name = entry
      results = results.join(featureVector, how='outer')

  try:
    logger.info('Extraction complete, writing output (format %s)', args.format)
    if args.format == 'csv':
      results.T.to_csv(args.outFile, index=False, na_rep='NaN')
    elif args.format == 'json':
      results.T.to_json(args.outFile, format='records')
    logger.info('writing complete')
  except Exception:
    logger.error('Unable to write to output file', exc_info=True)

  args.outFile.close()
  if args.log_file is not None:
    args.log_file.close()


if __name__ == "__main__":
  main()
