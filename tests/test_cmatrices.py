# to run this test, from directory above:
# setenv PYTHONPATH /path/to/pyradiomics/radiomics
# nosetests --nocapture -v tests/test_features.py

import logging

from nose_parameterized import parameterized
import numpy
import six

from radiomics import cMatsEnabled, getFeatureClasses
from testUtils import custom_name_func, RadiomicsTestUtils


testUtils = RadiomicsTestUtils()

testCases = ('brain1', 'brain2', 'breast1', 'lung1', 'lung2')

featureClasses = getFeatureClasses()

class TestFeatures:

  def generate_scenarios():
    global testCases, featureClasses

    for testCase in testCases:
      for className, featureClass in six.iteritems(featureClasses):
        assert(featureClass is not None)
        if "_calculateCMatrix" in dir(featureClass) or className == "shape":
          logging.debug('generate_scenarios: featureClass = %s', className)
          yield testCase, className

  global testUtils

  @parameterized.expand(generate_scenarios(), testcase_func_name=custom_name_func)
  def test_scenario(self, testCase, featureClassName):
    print("")
    global testUtils, featureClasses

    logging.debug('test_scenario: testCase = %s, featureClassName = %s', testCase, featureClassName)

    assert cMatsEnabled()

    testUtils.setFeatureClassAndTestCase(featureClassName, testCase)

    testImage = testUtils.getImage()
    testMask = testUtils.getMask()

    featureClass = featureClasses[featureClassName](testImage, testMask, **testUtils.getSettings())

    if featureClassName == 'shape':
      cSA = getattr(featureClass, 'SurfaceArea')  # pre-calculated value by C extension
      assert (cSA is not None)

      pySA = getattr(featureClass, '_calculateSurfaceArea')()  # Function, call to calculate SA in full-python mode
      assert (pySA is not None)

      # Check if the calculated values match
      assert (numpy.abs(pySA - cSA)) < 1e-3

    else:
      assert "_calculateMatrix" in dir(featureClass)

      cMat = featureClass._calculateCMatrix()
      assert cMat is not None

      pyMat = featureClass._calculateMatrix()
      assert pyMat is not None

      # Check if the calculated arrays match
      assert numpy.max(numpy.abs(pyMat - cMat)) < 1e-3
