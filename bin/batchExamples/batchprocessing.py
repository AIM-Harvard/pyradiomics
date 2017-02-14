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
  logging.getLogger().setLevel(logging.INFO)  # Prints out log messages in this script

  # Write out all log entries to log file
  handler = logging.FileHandler(filename=progress_filename, mode='w')
  formatter = logging.Formatter("%(levelname)s:%(name)s: %(message)s")
  handler.setFormatter(formatter)
  logging.getLogger().addHandler(handler)

  # Package logging
  # radiomics.debug()  # Switch on radiomics logging from level=DEBUG (default level=WARNING)
  # Alternative: specify level
  radiomics.logger.setLevel(logging.INFO)

  # Prevent radiomics logger from printing out log entries with level < WARNING to the console
  radiomics.logger.handlers[0].setLevel(logging.WARNING)

  logging.info('Loading CSV')
  print "Loading CSV"

  flists = []
  try:
    with open(inputCSV, 'rb') as inFile:
      cr = csv.reader(inFile, lineterminator='\n')
      flists = [row for row in cr]
  except Exception:
    logging.error('CSV READ FAILED:\n%s', traceback.format_exc())

  print "Loading Done"
  print ("Patients: " + str(len(flists)))

  kwargs = {}
  kwargs['binWidth'] = 25
  kwargs['resampledPixelSpacing'] = None  # [3,3,3]
  kwargs['interpolator'] = sitk.sitkBSpline
  kwargs['verbose'] = True

  logging.info('pyradiomics version: %s', radiomics.__version__)
  logging.info('Extracting features with kwarg settings: %s', str(kwargs))

  extractor = featureextractor.RadiomicsFeaturesExtractor(**kwargs)
  extractor.enableInputImages(original={})
  # extractor.enableInputImages(wavelet= {'level': 2})
  for idx, entry in enumerate(flists, start=1):

    print "(%d/%d) Processing Patient: %s, Study: %s, Reader: %s" % (idx, len(flists), entry[0], entry[1], entry[2])
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
        featureVector.update(extractor.execute(imageFilepath, maskFilepath))

        with open(outputFilepath, 'ab') as outputFile:
          writer = csv.writer(outputFile, lineterminator='\n')
          if idx == 1: writer.writerow(featureVector.keys())
          writer.writerow(featureVector.values())
      except Exception:
        logging.error('FEATURE EXTRACTION FAILED:\n%s', traceback.format_exc())


main()
