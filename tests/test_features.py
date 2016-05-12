# to run this test, from directory above:
# setenv PYTHONPATH /path/to/pyradiomics/radiomics
# nosetests --nocapture -v tests/test_features.py

from radiomics import imageoperations, firstorder, glcm, rlgl, shape, glszm
from testUtils import RadiomicsTestUtils
import sys, os
import logging
from nose_parameterized import parameterized
from itertools import product

testUtils = RadiomicsTestUtils()
defaultTestCases = testUtils.getTestCases()
defaultFeatures = ["firstorder", "glcm", "rlgl", "shape", "glszm"]

testCases = defaultTestCases
features = ["firstorder"]

featureClass = None

class TestFeatures:

  def generate_scenarios():
    global testCases
    global features

    for testCase in testCases:
      for feature in features:
        featureNames = None
        if feature == 'firstorder':
          featureNames = firstorder.RadiomicsFirstOrder.getFeatureNames()
        assert(featureNames != None)
        logging.debug('generate_scenarios: featureNames = %s', featureNames)
        for featureName in featureNames:
          yield testCase, feature, featureName

  global testUtils
  @parameterized.expand(generate_scenarios(), testcase_func_name=testUtils.custom_name_func)
  def test_scenario(self, testCase, featureClassName, featureName):
    print("")
    global testUtils

    testCaseChanged = testUtils.setTestCase(testCase)
    featureChanged = testUtils.setFeatureClassName(featureClassName)

    global featureClass
    if featureClass is None or testCaseChanged:
      if featureClassName == 'firstorder':
        logging.info('making a first order')
        print 'making a first order'
        featureClass = firstorder.RadiomicsFirstOrder(testUtils.getImage(), testUtils.getMask())
    assert (featureClass != None)

    featureClass.disableAllFeatures()
    featureClass.enableFeatureByName(featureName)
    featureClass.calculateFeatures()
    # get the result and test it
    val = featureClass.featureValues[featureName]
    testUtils.checkResult(featureName, val)

def teardown_module():
  print("")
  print("teardown")
  res = testUtils.getResults()
  print res
