# to run this test, from directory above:
# setenv PYTHONPATH /path/to/pyradiomics/radiomics
# nosetests --nocapture -v tests/test_matching.py

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
# testCases = ["breast1"]
features = defaultFeatures
# features = ["firstorder"]

featureClass = None

class TestMapping:

  def generate_scenarios():
    global testCases
    global features

    for testCase in testCases:
      for feature in features:
        featureNames = None
        if feature == 'firstorder':
          featureNames = firstorder.RadiomicsFirstOrder.getFeatureNames()
        elif feature == 'glcm':
          featureNames = glcm.RadiomicsGLCM.getFeatureNames()
        elif feature == 'rlgl':
          featureNames = rlgl.RadiomicsRLGL.getFeatureNames()
        elif feature == 'shape':
          featureNames = shape.RadiomicsShape.getFeatureNames()
        elif feature == 'glszm':
          featureNames = glszm.RadiomicsGLSZM.getFeatureNames()
        assert(featureNames != None)
        logging.debug('generate_scenarios: featureNames = %s', featureNames)
        for featureName in featureNames:
          yield testCase, feature, featureName

  global testUtils
  @parameterized.expand(generate_scenarios(), testcase_func_name=testUtils.custom_name_func)
  def test_scenario(self, testCase, featureClassName, featureName):
    print("")
    global testUtils

    logging.debug('test_scenario: testCase = %s, featureClassName = %s, featureName = %s', testCase, featureClassName, featureName)

    testCaseChanged = testUtils.setTestCase(testCase)
    featureChanged = testUtils.setFeatureClassName(featureClassName)

    matlabString = testUtils.getBaselineFeatureClassAndName(featureName)

    logging.debug('Baseline = %s, for testCase = %s, featureClassName = %s, featureName = %s', matlabString, testCase, featureClassName, featureName)

    assert (matlabString != None)

def teardown_module():
  print("")

