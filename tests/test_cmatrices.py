# to run this test, from directory above:
# setenv PYTHONPATH /path/to/pyradiomics/radiomics
# nosetests --nocapture -v tests/test_features.py

import logging

from nose_parameterized import parameterized
import numpy

from radiomics import glcm, gldm, gldzm, glrlm, glszm, ngtdm, shape
from testUtils import custom_name_func, RadiomicsTestUtils


testUtils = RadiomicsTestUtils()
defaultTestCases = testUtils.getTestCases()
defaultFeatures = ["glcm", "gldm", "ngtdm", "glrlm", "shape", "glszm", "gldzm"]

testCases = defaultTestCases
features = ["glcm", "gldm", "ngtdm", "glrlm", "glszm", "gldzm"]  # defaultFeatures


class TestFeatures:

  def generate_scenarios():
    global testCases
    global features

    for testCase in testCases:
      for feature in features:
        assert(feature is not None)
        logging.debug('generate_scenarios: featureClass = %s', feature)
        yield testCase, feature

  global testUtils

  @parameterized.expand(generate_scenarios(), testcase_func_name=custom_name_func)
  def test_scenario(self, testCase, featureClassName):
    print("")
    global testUtils

    logging.debug('test_scenario: testCase = %s, featureClassName = %s', testCase, featureClassName)

    testUtils.setFeatureClassAndTestCase(featureClassName, testCase)

    testImage = testUtils.getImage()
    testMask = testUtils.getMask()

    pyMat = None
    cMat = None
    featureClass = None

    if featureClassName == 'firstorder':
      logging.debug('No C implementation in firstorder, skipping test')
    elif featureClassName == 'glcm':
      logging.debug('Init GLCM')
      featureClass = glcm.RadiomicsGLCM(testImage, testMask, **testUtils.getKwargs())
      cMat = featureClass.P_glcm
      pyMat = featureClass._calculateGLCM()
    elif featureClassName == 'gldm':
      logging.debug('Init GLDM')
      featureClass = gldm.RadiomicsGLDM(testImage, testMask, **testUtils.getKwargs())
      cMat = featureClass.P_gldm
      pyMat = featureClass._calculateGLDM()
    elif featureClassName == 'ngtdm':
      logging.debug('Init NGTDM')
      featureClass = ngtdm.RadiomicsNGTDM(testImage, testMask, **testUtils.getKwargs())
      cMat = featureClass.P_ngtdm
      pyMat = featureClass._calculateNGTDM()
    elif featureClassName == 'glrlm':
      logging.debug('Init GLRLM')
      featureClass = glrlm.RadiomicsGLRLM(testImage, testMask, **testUtils.getKwargs())
      cMat = featureClass.P_glrlm
      pyMat = featureClass._calculateGLRLM()
    elif featureClassName == 'shape':
      logging.debug('No C implementation in glrlm, skipping test')
      # No C implementation yet, will follow
      # logging.debug('Init Shape')
      # featureClass = shape.RadiomicsShape(testImage, testMask, **testUtils.getKwargs())
    elif featureClassName == 'glszm':
      logging.debug('Init GLSZM')
      featureClass = glszm.RadiomicsGLSZM(testImage, testMask, **testUtils.getKwargs())
      cMat = featureClass.P_glszm
      pyMat = featureClass._calculateGLSZM()
    elif featureClassName == 'gldzm':
      logging.debug('Init GLDZM')
      featureClass = gldzm.RadiomicsGLDZM(testImage, testMask, **testUtils.getKwargs())
      cMat = featureClass.P_gldzm
      pyMat = featureClass._calculateGLDZM()

    assert (featureClass is not None)
    # Check if the calculated arrays match
    assert numpy.max(numpy.abs(pyMat - cMat)) < 1e-3
