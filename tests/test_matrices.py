# to run this test, from directory above:
# setenv PYTHONPATH /path/to/pyradiomics/radiomics
# nosetests --nocapture -v tests/test_features.py

import logging
import os

from nose_parameterized import parameterized
import numpy
import six

from radiomics import getFeatureClasses, testCases
from testUtils import custom_name_func, RadiomicsTestUtils


testUtils = RadiomicsTestUtils()

featureClasses = getFeatureClasses()


class TestMatrices:

  def generate_scenarios():
    global featureClasses

    for testCase in testCases:
      if testCase.startswith('test'):
        continue
      for className, featureClass in six.iteritems(featureClasses):
        assert(featureClass is not None)
        if "_calculateMatrix" in dir(featureClass):
          logging.debug('generate_scenarios: featureClass = %s', className)
          yield testCase, className

  @parameterized.expand(generate_scenarios(), testcase_func_name=custom_name_func)
  def test_scenario(self, test, featureClassName):
    global testUtils, featureClasses

    logging.debug('test_scenario: testCase = %s, featureClassName = %s', test, featureClassName)

    baselineFile = os.path.join(testUtils.getDataDir(), 'baseline', '%s_%s.npy' % (test, featureClassName))
    assert os.path.isfile(baselineFile)

    baselineMatrix = numpy.load(baselineFile)

    testUtils.setFeatureClassAndTestCase(featureClassName, test)

    testImage = testUtils.getImage('original')
    testMask = testUtils.getMask('original')

    featureClass = featureClasses[featureClassName](testImage, testMask, **testUtils.getSettings())
    featureClass._initCalculation()

    cMat = getattr(featureClass, 'P_%s' % featureClassName)
    assert cMat is not None

    # Check if the calculated arrays match
    assert numpy.max(numpy.abs(baselineMatrix - cMat)) < 1e-3
