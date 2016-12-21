# to run this test, from directory above:
# setenv PYTHONPATH /path/to/pyradiomics/radiomics
# nosetests --nocapture -v tests/test_features.py

import logging

import numpy
import SimpleITK as sitk

from nose_parameterized import parameterized
from radiomics import glcm, glrlm, glszm

from .testUtils import RadiomicsTestUtils


testUtils = RadiomicsTestUtils()
defaultTestCases = testUtils.getTestCases()
defaultFeatures = ["glcm", "glrlm", "shape", "glszm"]

testCases = defaultTestCases
features = ["glcm", "glrlm", "glszm"]  # defaultFeatures

# set testing arguments
kwargs = {}
kwargs['binWidth'] = 25
kwargs['symmetricalGLCM'] = False  # Current baseline is based upon assymetrical GLCM
kwargs['verbose'] = False

testUtils.setResampling(resampledPixelSpacing= None,interpolator= sitk.sitkBSpline)  # resampledPixelSpacing= [3, 3, 3]

class TestFeatures:

  def generate_scenarios():
    global testCases
    global features

    for testCase in testCases:
      for feature in features:
        assert(feature != None)
        logging.debug('generate_scenarios: featureClass = %s', feature)
        yield testCase, feature

  global testUtils
  @parameterized.expand(generate_scenarios(), testcase_func_name=testUtils.custom_name_func)
  def test_scenario(self, testCase, featureClassName):
    print("")
    global testUtils

    logging.debug('test_scenario: testCase = %s, featureClassName = %s', testCase, featureClassName)

    testUtils.setTestCase(testCase)

    testImage = testUtils.getImage()
    testMask = testUtils.getMask()

    pyMat = None
    cMat = None
    featureClass = None

    if featureClassName == 'firstorder':
      logging.debug('No C implementation in firstorder, skipping test')
    elif featureClassName == 'glcm':
      logging.debug('Init GLCM')
      featureClass = glcm.RadiomicsGLCM(testImage, testMask, **kwargs)
      cMat = featureClass.P_glcm
      pyMat = featureClass._calculateGLCM()
    elif featureClassName == 'glrlm':
      logging.debug('Init GLRLM')
      featureClass = glrlm.RadiomicsGLRLM(testImage, testMask, **kwargs)
      cMat = featureClass.P_glrlm
      pyMat = featureClass._calculateGLRLM()
    elif featureClassName == 'shape':
      logging.debug('No C implementation in glrlm, skipping test')
      # No C implementation yet, will follow
      # logging.debug('Init Shape')
      # featureClass = shape.RadiomicsShape(testImage, testMask, **kwargs)
    elif featureClassName == 'glszm':
      logging.debug('Init GLSZM')
      featureClass = glszm.RadiomicsGLSZM(testImage, testMask, **kwargs)
      cMat = featureClass.P_glszm
      pyMat = featureClass._calculateGLSZM()

    assert (featureClass != None)
    # Check if the calculated arrays match
    assert numpy.max(numpy.abs(pyMat - cMat)) < 1e-3
