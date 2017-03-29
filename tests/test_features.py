# to run this test, from directory above:
# setenv PYTHONPATH /path/to/pyradiomics/radiomics
# nosetests --nocapture -v tests/test_features.py

import logging
import os

from nose_parameterized import parameterized
import numpy

from radiomics.featureextractor import RadiomicsFeaturesExtractor
from testUtils import custom_name_func, RadiomicsTestUtils

testUtils = RadiomicsTestUtils()
testCases = testUtils.getTestCases()

extractor = RadiomicsFeaturesExtractor()

featureClass = None


class TestFeatures:
  def generate_scenarios():
    global testCases
    global extractor

    for testCase in testCases:
      for featureClassName in extractor.getFeatureClassNames():
        featureNames = extractor.featureClasses[featureClassName].getFeatureNames()
        assert (featureNames is not None)
        assert (len(featureNames) > 0)
        logging.debug('generate_scenarios: featureNames = %s', featureNames)
        for featureName in featureNames:
          yield testCase, featureClassName, featureName

  global testUtils

  @parameterized.expand(generate_scenarios(), testcase_func_name=custom_name_func)
  def test_scenario(self, testCase, featureClassName, featureName):
    print("")
    global testUtils
    global extractor

    logging.debug('test_scenario: testCase = %s, featureClassName = %s, featureName = %s', testCase, featureClassName,
                  featureName)

    testCaseOrClassChanged = testUtils.setFeatureClassAndTestCase(featureClassName, testCase)

    global featureClass
    testImage = testUtils.getImage()
    testMask = testUtils.getMask()
    if featureClass is None or testCaseOrClassChanged:
      logging.debug('Init %s' % (featureClassName))
      featureClass = extractor.featureClasses[featureClassName](testImage, testMask, **testUtils.getKwargs())

    assert (featureClass is not None)

    try:
      val = getattr(featureClass, 'get%sFeatureValue' % featureName)()
    except Exception:
      logging.error('Feature calculation Failed', exc_info=True)
      val = numpy.nan

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
