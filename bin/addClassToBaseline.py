
from __future__ import print_function

import collections
import csv
import os

import six

from radiomics import featureextractor, imageoperations


def main():
  dataDir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data")
  baselineDir = os.path.join(dataDir, "baseline")

  testCases = []

  kwargs = {'binWidth': 25,
            'interpolator': None,
            'resampledPixelSpacing': None,
            'padDistance': 5,
            'voxelArrayShift': 2000,
            'symmetricalGLCM': False,
            'weightingNorm': None,
            'gldm_a': 0}

  extractor = featureextractor.RadiomicsFeaturesExtractor(**kwargs)

  for cls in extractor.getFeatureClassNames():
    if os.path.exists(os.path.join(baselineDir, 'baseline_%s.csv' % (cls))):
      with open(os.path.join(baselineDir, 'baseline_%s.csv' % (cls)), 'rb') as baselineFile:
        csvReader = csv.reader(baselineFile)
        six.next(csvReader)  # Skip header row
        for testRow in csvReader:
          testCases += [testRow[0]]
      if len(testCases) > 0: break

  if len(testCases) == 0:
    print("No baselinefiles containing testcases found, exiting...")
    exit(-1)

  newClasses = [cls for cls in extractor.getFeatureClassNames() if
                not os.path.exists(os.path.join(baselineDir, 'baseline_%s.csv' % (cls)))]

  if len(newClasses) == 0:
    print("No new classes to add, exiting...")
    exit(0)

  print("Adding new classes: ", newClasses)

  newBaseline = {}

  # When C matrices is merged, force baseline to be build using the Python implementation of the functions.
  # radiomics.pythonMatrixCalculation(True)

  extractor.disableAllFeatures()
  extractor.disableAllInputImages()

  # Calculate new baseline on original image for all new classes
  extractor.enableInputImageByName('original')
  for cls in newClasses:
    newBaseline[cls] = {}

  print("Computing new baseline")
  for testCase in testCases:
    print("\tCalculating test case", testCase)
    imagePath = os.path.join(dataDir, testCase + '_image.nrrd')
    maskPath = os.path.join(dataDir, testCase + '_label.nrrd')
    image, mask = extractor.loadImage(imagePath, maskPath)
    if image is None or mask is None:
      print("Error during loading of image/mask, testcase:", testCase)
      continue  # testImage or mask not found / error during loading

    provenance = extractor.getProvenance(imagePath, maskPath, mask)

    bb, correctedMask = imageoperations.checkMask(image, mask)
    if correctedMask is not None:
      mask = correctedMask
    image, mask = imageoperations.cropToTumorMask(image, mask, bb)
    for cls in newClasses:
      print("\t\tCalculating class", cls)
      newBaseline[cls][testCase] = collections.OrderedDict()
      newBaseline[cls][testCase]["Patient ID"] = testCase
      newBaseline[cls][testCase].update(provenance)
      featureClass = extractor.featureClasses[cls](image, mask, **extractor.settings)
      featureClass.enableAllFeatures()
      featureClass.calculateFeatures()
      newBaseline[cls][testCase].update(featureClass.featureValues)

  print("Writing new baseline")
  for cls in newClasses:
    baselineFile = os.path.join(baselineDir, 'baseline_%s.csv' % (cls))
    with open(baselineFile, 'wb') as baseline:
      csvWriter = csv.writer(baseline)
      header = newBaseline[cls][testCases[0]].keys()
      csvWriter.writerow(header)
      for testCase in testCases:
        row = []
        for h in header:
          row += [newBaseline[cls][testCase][h]]
        csvWriter.writerow(row)
    baseline.close()


if __name__ == '__main__':
  main()
