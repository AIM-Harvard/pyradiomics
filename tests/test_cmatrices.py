# to run this test, from directory above:
# setenv PYTHONPATH /path/to/pyradiomics/radiomics
# nosetests --nocapture -v tests/test_features.py

import radiomics
radiomics.debug(True)  # Ensure debugging is switched off (first calculate C matrices)

from radiomics import glcm, glrlm, shape, glszm,ngtdm, gldm
from testUtils import RadiomicsTestUtils
import logging
from nose_parameterized import parameterized
import SimpleITK as sitk
import numpy


testUtils = RadiomicsTestUtils()
defaultTestCases = testUtils.getTestCases()
defaultFeatures = ["glcm", "gldm", "ngtdm", "glrlm", "shape", "glszm"]

testCases = defaultTestCases
features = ["glcm", "gldm", "ngtdm", "glszm"]  # defaultFeatures

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
      pyMat = featureClass.P_glcm
      cMat = featureClass._calculateCGLCM()
    elif featureClassName == 'gldm':
      logging.debug('Init GLDM')
      featureClass = gldm.RadiomicsGLDM(testImage, testMask, **kwargs)
      pyMat = featureClass.P_gldm
      cMat = featureClass._calculateCGLDM()
    elif featureClassName == 'ngtdm':
      logging.debug('Init NGTDM')
      featureClass = ngtdm.RadiomicsNGTDM(testImage, testMask, **kwargs)
      pyMat = featureClass.P_ngtdm
      cMat = featureClass._calculateCNGTDM()
    elif featureClassName == 'glrlm':
      logging.debug('No C implementation in glrlm, skipping test')
      # No C implementation yet, will follow
      # logging.debug('Init GLRLM')
      # featureClass = glrlm.RadiomicsGLRLM(testImage, testMask, **kwargs)
      # cMat = featureClass.P_glrlm
      # pyMat = featureClass._calculateCGLRLM()
    elif featureClassName == 'shape':
      logging.debug('No C implementation in glrlm, skipping test')
      # No C implementation yet, will follow
      # logging.debug('Init Shape')
      # featureClass = shape.RadiomicsShape(testImage, testMask, **kwargs)
    elif featureClassName == 'glszm':
      logging.debug('Init GLSZM')
      featureClass = glszm.RadiomicsGLSZM(testImage, testMask, **kwargs)
      pyMat = featureClass.P_glszm
      cMat = featureClass._calculateCGLSZM()
    else:
      assert (featureClass != None)

    # Check if the calculated arrays match
    assert numpy.max(numpy.abs(pyMat - cMat)) < 1e-3
