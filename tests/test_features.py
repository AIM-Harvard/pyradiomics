# to run this test, from directory above:
# setenv PYTHONPATH /path/to/pyradiomics/radiomics
# nosetests --nocapture -v tests/test_features.py

import logging
import os

from nose_parameterized import parameterized
import six

from radiomics import getFeatureClasses
from testUtils import custom_name_func, RadiomicsTestUtils

testUtils = RadiomicsTestUtils()
tests = sorted(testUtils.getTests())

featureClass = None
featureClasses = getFeatureClasses()


class TestFeatures:

  def generate_scenarios():
    global tests, featureClasses

    for test in tests:
      for featureClassName in featureClasses:
        # Get all feature names for which there is a baseline with current test case
        # Raises an assertion error when the class is not yet present in the baseline
        # Returns None if no baseline is present for this specific test case
        # Returns a list of feature names for which baseline values are present for this test
        baselineFeatureNames = testUtils.getFeatureNames(featureClassName, test)

        if baselineFeatureNames is None:
          continue
        assert (len(baselineFeatureNames) > 0)

        uniqueFeatures = set([f.split('_')[-1] for f in baselineFeatureNames])

        # Get a list of all features for current class
        featureNames = featureClasses[featureClassName].getFeatureNames()
        # Get a list of all non-deprecated features
        activeFeatures = set([f for (f, deprecated) in six.iteritems(featureNames) if not deprecated])
        # Check if all active features have a baseline (exclude deprecated features from this set)
        if len(activeFeatures - uniqueFeatures) > 0:
          raise AssertionError('Missing baseline for active features %s', activeFeatures - uniqueFeatures)
        if len(uniqueFeatures - activeFeatures) > 0:
          raise AssertionError('Missing function(s) for baseline feature(s) %s', uniqueFeatures - activeFeatures)

        logging.debug('generate_scenarios: featureNames = %s', baselineFeatureNames)
        for featureName in baselineFeatureNames:
          yield test, featureName

  @parameterized.expand(generate_scenarios(), testcase_func_name=custom_name_func)
  def test_scenario(self, test, featureName):
    print("")
    global testUtils, featureClass, featureClasses

    featureName = featureName.split('_')

    logging.debug('test_scenario: test = %s, featureClassName = %s, featureName = %s', test, featureName[1],
                  featureName[-1])

    testOrClassChanged = testUtils.setFeatureClassAndTestCase(featureName[1], test)

    testImage = testUtils.getImage(featureName[0])
    testMask = testUtils.getMask(featureName[0])

    if featureClass is None or testOrClassChanged:
      logging.debug('Init %s' % featureName[1])
      featureClass = featureClasses[featureName[1]](testImage, testMask, **testUtils.getSettings())

    assert (featureClass is not None)

    featureClass.disableAllFeatures()
    featureClass.enableFeatureByName(featureName[-1])
    featureClass.execute()
    # get the result and test it
    val = featureClass.featureValues[featureName[-1]]
    testUtils.checkResult(featureName, val)


def teardown_module():
  print("")
  res = testUtils.getResults()
  print('Results:')
  print(res)
  resultsFile = os.path.join(testUtils.getDataDir(), 'PyradiomicsFeatures.csv')
  testUtils.writeCSV(res, resultsFile)
  diff = testUtils.getDiffs()
  print('Differences from baseline:')
  print(diff)
  diffFile = os.path.join(testUtils.getDataDir(), 'Baseline2PyradiomicsFeaturesDiff.csv')
  testUtils.writeCSV(diff, diffFile)
  logging.info(
    "Wrote calculated features to %s, and the differences between the baseline features and the calculated ones to %s.",
    resultsFile, diffFile)
