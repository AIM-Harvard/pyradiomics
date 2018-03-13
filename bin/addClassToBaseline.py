
from __future__ import print_function

import collections
import csv
import os

from radiomics import featureextractor, getTestCase

TEST_CASES = ('brain1', 'brain2', 'breast1', 'lung1', 'lung2')


def main():
  global TEST_CASES
  dataDir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data")
  baselineDir = os.path.join(dataDir, "baseline")

  kwargs = {'binWidth': 25,
            'interpolator': None,
            'resampledPixelSpacing': None,
            'padDistance': 5,
            'voxelArrayShift': 2000,
            'weightingNorm': None,
            'gldm_a': 0}

  extractor = featureextractor.RadiomicsFeaturesExtractor(**kwargs)

  newClasses = [cls for cls in extractor.getFeatureClassNames() if
                not os.path.exists(os.path.join(baselineDir, 'baseline_%s.csv' % (cls)))]

  if len(newClasses) == 0:
    print("No new classes to add, exiting...")
    exit(0)

  print("Adding new classes: ", newClasses)

  newBaseline = {}

  extractor.disableAllFeatures()
  extractor.disableAllImageTypes()

  # Calculate new baseline on original image for all new classes
  extractor.enableImageTypeByName('Original')
  for cls in newClasses:
    newBaseline[cls] = {}

  print("Computing new baseline")
  for testCase in TEST_CASES:
    print("\tCalculating test case", testCase)
    imagePath, maskPath = getTestCase(testCase)
    image, mask = extractor.loadImage(imagePath, maskPath)
    if image is None or mask is None:
      print("Error during loading of image/mask, testcase:", testCase)
      continue  # testImage or mask not found / error during loading

    for cls in newClasses:
      print("\t\tCalculating class", cls)
      extractor.disableAllFeatures()
      extractor.enableFeatureClassByName(cls)
      newBaseline[cls][testCase] = collections.OrderedDict()
      newBaseline[cls][testCase]['general_info_TestCase'] = testCase
      newBaseline[cls][testCase].update(extractor.execute(image, mask))

  print("Writing new baseline")
  for cls in newClasses:
    baselineFile = os.path.join(baselineDir, 'baseline_%s.csv' % (cls))
    with open(baselineFile, 'wb') as baseline:
      csvWriter = csv.writer(baseline)
      header = ['featureName'] + list(TEST_CASES)
      csvWriter.writerow(header)

      features = newBaseline[cls][TEST_CASES[0]].keys()
      for f in features:
        row = [f]
        for testCase in TEST_CASES:
          row.append(newBaseline[cls][testCase].get(f, ''))
        csvWriter.writerow(row)
    baseline.close()


if __name__ == '__main__':
  main()
