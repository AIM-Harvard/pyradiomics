
from __future__ import print_function

import collections
import csv
import logging
import os
import traceback

import SimpleITK as sitk

import radiomics
from radiomics import featureextractor


def main():
  outPath = r'E:\Git-Repos\PyRadiomicsNKI\pyradiomics-API'

  inputCSV = outPath + os.path.sep + "TestCases.csv"
  outputFilepath = outPath + os.path.sep + "radiomics_features.csv"
  progress_filename = outPath + os.path.sep + "pyrad_log.txt"

  # Enable writing out the log using radiomics logger
  logLevel = logging.INFO
  rLogger = logging.getLogger('radiomics')
  rLogger.setLevel(logLevel)
  handler = logging.FileHandler(filename=progress_filename, mode='w')
  handler.setLevel(logLevel)
  handler.setFormatter(logging.Formatter("%(levelname)s:%(name)s: %(message)s"))
  rLogger.addHandler(handler)

  # Prevent radiomics logger from printing out log entries with level < WARNING to the console
  radiomics.logger.handlers[0].setLevel(logging.WARNING)

  # Initialize logging for batch log messages
  logger = logging.getLogger('radiomics.batch')

  logger.info('Loading CSV')
  print("Loading CSV")

  flists = []
  try:
    with open(inputCSV, 'r') as inFile:
      cr = csv.reader(inFile, lineterminator='\n')
      flists = [row for row in cr]
  except Exception:
    logging.error('CSV READ FAILED:\n%s', traceback.format_exc())

  print("Loading Done")
  print("Patients: " + str(len(flists)))

  kwargs = {}
  kwargs['binWidth'] = 25
  kwargs['resampledPixelSpacing'] = None  # [3,3,3]
  kwargs['interpolator'] = sitk.sitkBSpline
  kwargs['verbose'] = True
  kwargs['enableCExtensions'] = False

  logger.info('pyradiomics version: %s', radiomics.__version__)
  logger.info('Extracting features with kwarg settings: %s', str(kwargs))

  extractor = featureextractor.RadiomicsFeaturesExtractor(**kwargs)
  extractor.enableInputImages(Original={})
  # extractor.enableInputImages(wavelet= {'level': 2})

  headers = None

  for idx, entry in enumerate(flists, start=1):

    print("(%d/%d) Processing Patient: %s, Study: %s, Reader: %s" % (idx, len(flists), entry[0], entry[1], entry[2]))
    logger.info("(%d/%d) Processing Patient: %s, Study: %s, Reader: %s", idx, len(flists), entry[0], entry[1],
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
        logger.error('FEATURE EXTRACTION FAILED:\n%s', traceback.format_exc())

if __name__ == '__main__':
  main()
