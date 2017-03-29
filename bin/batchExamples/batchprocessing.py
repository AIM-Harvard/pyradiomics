
from __future__ import print_function

import logging
import os

import pandas
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

  try:
    flists = pandas.read_csv(inputCSV).T  # Transpose data so that each column represents one test case
  except Exception:
    logging.error('CSV READ FAILED', exc_info=True)
    exit(-1)

  logging.info('Loading Done')
  logging.info('Patients: %d', len(flists))

  kwargs = {}
  kwargs['binWidth'] = 25
  kwargs['resampledPixelSpacing'] = None  # [3,3,3]
  kwargs['interpolator'] = sitk.sitkBSpline
  kwargs['enableCExtensions'] = True

  logger.info('pyradiomics version: %s', radiomics.__version__)
  logger.info('Extracting features with kwarg settings: %s', str(kwargs))

  extractor = featureextractor.RadiomicsFeaturesExtractor(**kwargs)
  # extractor.enableInputImages(Original={}) # Original enabled by default
  # extractor.enableInputImages(wavelet= {'level': 2})

  results = pandas.DataFrame()

  for entry in flists:  # Loop over all columns (i.e. the test cases)
    logger.info("(%d/%d) Processing Patient (Image: %s, Mask: %s)",
                entry + 1,
                len(flists),
                flists[entry]['Image'],
                flists[entry]['Mask'])

    imageFilepath = flists[entry]['Image']
    maskFilepath = flists[entry]['Mask']

    if (imageFilepath is not None) and (maskFilepath is not None):
      featureVector = flists[entry]
      featureVector['Image'] = os.path.basename(imageFilepath)
      featureVector['Mask'] = os.path.basename(maskFilepath)

      try:
        featureVector = featureVector.append(extractor.execute(imageFilepath, maskFilepath))
      except Exception:
        logger.error('FEATURE EXTRACTION FAILED:', exc_info=True)

      featureVector.name = entry
      results = results.join(featureVector, how='outer')  # If feature extraction failed, results will be all NaN

  logger.info('Extraction complete, writing CSV')
  results.T.to_csv(outputFilepath, index=False, na_rep='NaN')
  logger.info('CSV writing complete')

if __name__ == '__main__':
  main()
