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
# testCases = ["breast1"]
# features = defaultFeatures
features = ["glcm"]

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

    global featureClass
    testImage = testUtils.getImage()
    testMask = testUtils.getMask()
    if featureClass is None or testCaseChanged or featureChanged:
      if featureClassName == 'firstorder':
        logging.debug('Init First Order')
        featureClass = firstorder.RadiomicsFirstOrder(testImage, testMask)
      elif featureClassName == 'glcm':
        logging.debug('Init GLCM')
        featureClass = glcm.RadiomicsGLCM(testImage, testMask)
      elif featureClassName == 'rlgl':
        logging.debug('Init RLGL')
        featureClass = rlgl.RadiomicsRLGL(testImage, testMask)
      elif featureClassName == 'shape':
        logging.debug('Init Shape')
        featureClass = shape.RadiomicsShape(testImage, testMask)
      elif featureClassName == 'glszm':
        logging.debug('Init GLSZM')
        featureClass = glszm.RadiomicsGLSZM(testImage, testMask)
    assert (featureClass != None)

    featureClass.disableAllFeatures()
    featureClass.enableFeatureByName(featureName)
    featureClass.calculateFeatures()
    # get the result and test it
    val = featureClass.featureValues[featureName]
    testUtils.checkResult(featureName, val)

def teardown_module():
  print("")
  res = testUtils.getResults()
  print res
  resultsFile = os.path.join(testUtils.dataDir, 'PyradiomicsFeatures.csv')
  testUtils.writeCSV(res, resultsFile)
  diff = testUtils.getDiffs()
  print 'Differences from baseline:'
  print diff
  diffFile = os.path.join(testUtils.dataDir, 'Matlab2PyradiomicsFeaturesDiff.csv')
  testUtils.writeCSV(diff, diffFile)
  logging.info("Wrote calculated features to %s, and the differences between the matlab features and the pyradiomics ones to %s.", resultsFile, diffFile)