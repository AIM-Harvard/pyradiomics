#!/usr/bin/env python

from __future__ import print_function

import logging
import os

import pandas
import SimpleITK as sitk

import radiomics
from radiomics import featureextractor


def main():
  outPath = r''

  inputCSV = os.path.join(outPath, 'testCases.csv')
  outputFilepath = os.path.join(outPath, 'radiomics_features.csv')
  progress_filename = os.path.join(outPath, 'pyrad_log.txt')
  params = os.path.join(outPath, 'exampleSettings', 'Params.yaml')

  # Configure logging
  rLogger = logging.getLogger('radiomics')

  # Set logging level
  # rLogger.setLevel(logging.INFO)  # Not needed, default log level of logger is INFO

  # Create handler for writing to log file
  handler = logging.FileHandler(filename=progress_filename, mode='w')
  handler.setFormatter(logging.Formatter('%(levelname)s:%(name)s: %(message)s'))
  rLogger.addHandler(handler)

  # Initialize logging for batch log messages
  logger = rLogger.getChild('batch')

  # Set verbosity level for output to stderr (default level = WARNING)
  radiomics.setVerbosity(logging.INFO)

  logger.info('pyradiomics version: %s', radiomics.__version__)
  logger.info('Loading CSV')

  # ####### Up to this point, this script is equal to the 'regular' batchprocessing script ########

  try:
    # Use pandas to read and transpose ('.T') the input data
    # The transposition is needed so that each column represents one test case. This is easier for iteration over
    # the input cases
    flists = pandas.read_csv(inputCSV).T
  except Exception:
    logger.error('CSV READ FAILED', exc_info=True)
    exit(-1)

  logger.info('Loading Done')
  logger.info('Patients: %d', len(flists.columns))

  if os.path.isfile(params):
    extractor = featureextractor.RadiomicsFeatureExtractor(params)
  else:  # Parameter file not found, use hardcoded settings instead
    settings = {}
    settings['binWidth'] = 25
    settings['resampledPixelSpacing'] = None  # [3,3,3]
    settings['interpolator'] = sitk.sitkBSpline
    settings['enableCExtensions'] = True

    extractor = featureextractor.RadiomicsFeatureExtractor(**settings)
    # extractor.enableInputImages(wavelet= {'level': 2})

  logger.info('Enabled input images types: %s', extractor.enabledImagetypes)
  logger.info('Enabled features: %s', extractor.enabledFeatures)
  logger.info('Current settings: %s', extractor.settings)

  # Instantiate a pandas data frame to hold the results of all patients
  results = pandas.DataFrame()

  for entry in flists:  # Loop over all columns (i.e. the test cases)
    logger.info("(%d/%d) Processing Patient (Image: %s, Mask: %s)",
                entry + 1,
                len(flists),
                flists[entry]['Image'],
                flists[entry]['Mask'])

    imageFilepath = flists[entry]['Image']
    maskFilepath = flists[entry]['Mask']
    label = flists[entry].get('Label', None)

    if str(label).isdigit():
      label = int(label)
    else:
      label = None

    if (imageFilepath is not None) and (maskFilepath is not None):
      featureVector = flists[entry]  # This is a pandas Series
      featureVector['Image'] = os.path.basename(imageFilepath)
      featureVector['Mask'] = os.path.basename(maskFilepath)

      try:
        # PyRadiomics returns the result as an ordered dictionary, which can be easily converted to a pandas Series
        # The keys in the dictionary will be used as the index (labels for the rows), with the values of the features
        # as the values in the rows.
        result = pandas.Series(extractor.execute(imageFilepath, maskFilepath, label))
        featureVector = featureVector.append(result)
      except Exception:
        logger.error('FEATURE EXTRACTION FAILED:', exc_info=True)

      # To add the calculated features for this case to our data frame, the series must have a name (which will be the
      # name of the column.
      featureVector.name = entry
      # By specifying an 'outer' join, all calculated features are added to the data frame, including those not
      # calculated for previous cases. This also ensures we don't end up with an empty frame, as for the first patient
      # it is 'joined' with the empty data frame.
      results = results.join(featureVector, how='outer')  # If feature extraction failed, results will be all NaN

  logger.info('Extraction complete, writing CSV')
  # .T transposes the data frame, so that each line will represent one patient, with the extracted features as columns
  results.T.to_csv(outputFilepath, index=False, na_rep='NaN')
  logger.info('CSV writing complete')


if __name__ == '__main__':
  main()
