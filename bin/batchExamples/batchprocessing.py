
from __future__ import print_function

import collections
import csv
import logging
import os

import SimpleITK as sitk

import radiomics
from radiomics import featureextractor


def main():
  outPath = r'E:\Git-Repos\PyRadiomicsNKI\pyradiomics-API'

  inputCSV = outPath + os.path.sep + "TestCases.csv"
  outputFilepath = outPath + os.path.sep + "radiomics_features.csv"
  progress_filename = outPath + os.path.sep + "pyrad_log.txt"

  # Configure logging
  rLogger = logging.getLogger('radiomics')

  # Set logging level
  # rLogger.setLevel(logging.INFO)  # Not needed, default log level of logger is INFO

  # Create handler for writing to log file
  handler = logging.FileHandler(filename=progress_filename, mode='w')
  handler.setFormatter(logging.Formatter("%(levelname)s:%(name)s: %(message)s"))
  rLogger.addHandler(handler)

  # Initialize logging for batch log messages
  logger = rLogger.getChild('batch')

  # Set verbosity level for output to stderr (default level = WARNING)
  radiomics.setVerbosity(logging.INFO)

  logger.info('Loading CSV')

  flists = []
  try:
    with open(inputCSV, 'r') as inFile:
      cr = csv.DictReader(inFile, lineterminator='\n')
      flists = [row for row in cr]
  except Exception:
    logging.error('CSV READ FAILED', exc_info=True)

  logger.info('Loading Done')
  logger.info('Patients: %d', len(flists))

  kwargs = {}
  kwargs['binWidth'] = 25
  kwargs['resampledPixelSpacing'] = None  # [3,3,3]
  kwargs['interpolator'] = sitk.sitkBSpline
  kwargs['enableCExtensions'] = False

  logger.info('pyradiomics version: %s', radiomics.__version__)
  logger.info('Extracting features with kwarg settings: %s', str(kwargs))

  extractor = featureextractor.RadiomicsFeaturesExtractor(**kwargs)
  extractor.enableInputImages(Original={})
  # extractor.enableInputImages(wavelet= {'level': 2})

  headers = None

  for idx, entry in enumerate(flists, start=1):

    logger.info("(%d/%d) Processing Patient (Image: %s, Mask: %s)", idx, len(flists), entry['Image'], entry['Mask'])

    imageFilepath = entry['Image']
    maskFilepath = entry['Mask']

    if (imageFilepath is not None) and (maskFilepath is not None):
      featureVector = collections.OrderedDict(entry)
      featureVector['Image'] = os.path.basename(imageFilepath)
      featureVector['Mask'] = os.path.basename(maskFilepath)

      try:
        featureVector.update(extractor.execute(imageFilepath, maskFilepath))

        with open(outputFilepath, 'a') as outputFile:
          writer = csv.writer(outputFile, lineterminator='\n')
          if headers is None:
            headers = list(featureVector.keys())
            writer.writerow(headers)

          row = []
          for h in headers:
            row.append(featureVector.get(h, "N/A"))
          writer.writerow(row)
      except Exception:
        logger.error('FEATURE EXTRACTION FAILED', exc_info=True)

if __name__ == '__main__':
  main()
