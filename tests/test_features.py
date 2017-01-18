# to run this test, from directory above:
# setenv PYTHONPATH /path/to/pyradiomics/radiomics
# nosetests --nocapture -v tests/test_features.py

from radiomics.featureextractor import RadiomicsFeaturesExtractor
from testUtils import RadiomicsTestUtils, custom_name_func
import os
import logging
from nose_parameterized import parameterized

testUtils = RadiomicsTestUtils()
testCases = testUtils.getTestCases()

extractor = RadiomicsFeaturesExtractor()
featureClassNames = extractor.getFeatureClassNames()

featureClass = None


class TestFeatures:
  def generate_scenarios():
    global testCases
    global extractor
    global featureClassNames

    for testCase in testCases:
      for featureClassName in featureClassNames:
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

    featureClass.disableAllFeatures()
    featureClass.enableFeatureByName(featureName)
    featureClass.calculateFeatures()
    # get the result and test it
    val = featureClass.featureValues[featureName]
    testUtils.checkResult(featureName, val)


def teardown_module():
  print("")
  res = testUtils.getResults()
  print 'Results:'
  print res
  resultsFile = os.path.join(testUtils.getDataDir(), 'PyradiomicsFeatures.csv')
  testUtils.writeCSV(res, resultsFile)
  diff = testUtils.getDiffs()
  print 'Differences from baseline:'
  print diff
  diffFile = os.path.join(testUtils.getDataDir(), 'Baseline2PyradiomicsFeaturesDiff.csv')
  testUtils.writeCSV(diff, diffFile)
  logging.info(
    "Wrote calculated features to %s, and the differences between the matlab features and the pyradiomics ones to %s.",
    resultsFile, diffFile)
