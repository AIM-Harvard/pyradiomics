import logging
import os

import numpy
import six

from radiomics import getFeatureClasses, testCases
from testUtils import RadiomicsTestUtils


testUtils = RadiomicsTestUtils()

featureClasses = getFeatureClasses()


def pytest_generate_tests(metafunc):
  metafunc.parametrize(["testCase", "featureClassName"], metafunc.cls.generate_scenarios())


class TestMatrices:

  @staticmethod
  def generate_scenarios():
    global featureClasses

    for testCase in testCases:
      if testCase.startswith('test'):
        continue
      for className, featureClass in six.iteritems(featureClasses):
        assert featureClass is not None
        if "_calculateMatrix" in dir(featureClass):
          logging.debug('generate_scenarios: featureClass = %s', className)
          yield testCase, className

  def test_scenario(self, testCase, featureClassName):
    global testUtils, featureClasses

    logging.debug('test_scenario: testCase = %s, featureClassName = %s', testCase, featureClassName)

    baselineFile = os.path.join(testUtils.getDataDir(), 'baseline', '%s_%s.npy' % (testCase, featureClassName))
    assert os.path.isfile(baselineFile)

    baselineMatrix = numpy.load(baselineFile)

    testUtils.setFeatureClassAndTestCase(featureClassName, testCase)

    testImage = testUtils.getImage('original')
    testMask = testUtils.getMask('original')

    featureClass = featureClasses[featureClassName](testImage, testMask, **testUtils.getSettings())
    featureClass._initCalculation()

    cMat = getattr(featureClass, 'P_%s' % featureClassName)
    assert cMat is not None

    # Check if the calculated arrays match
    assert numpy.max(numpy.abs(baselineMatrix - cMat)) < 1e-3
