# to run this test, from directory above:
# setenv PYTHONPATH /path/to/pyradiomics/radiomics
# nosetests --nocapture -v tests/test_features.py

from radiomics.featureextractor import RadiomicsFeaturesExtractor
from testUtils import RadiomicsTestUtils
import os
import logging
from nose_parameterized import parameterized
import SimpleITK as sitk

testUtils = RadiomicsTestUtils()
testUtils.setResampling(resampledPixelSpacing=None,  # resampledPixelSpacing= [3, 3, 3]
                        interpolator=sitk.sitkBSpline)
testCases = testUtils.getTestCases()

kwargs = {}
kwargs['binWidth'] = 25
kwargs['symmetricalGLCM'] = False  # Current baseline is based upon assymetrical GLCM
kwargs['verbose'] = False

extractor = RadiomicsFeaturesExtractor(**kwargs)

featureClass = None


class TestFeatures:
  def generate_scenarios():
    global testCases

    for testCase in testCases:
      for featureClassName, featureClass in extractor.featureClasses.iteritems():
        featureNames = featureClass.getFeatureNames()
        assert (featureNames != None)
        assert (len(featureNames) > 0)
        logging.debug('generate_scenarios: featureNames = %s', featureNames)
        for featureName in featureNames:
          yield testCase, featureClassName, featureName

  global testUtils

  @parameterized.expand(generate_scenarios(), testcase_func_name=testUtils.custom_name_func)
  def test_scenario(self, testCase, featureClassName, featureName):
    print("")
    global testUtils
    global extractor

    logging.debug('test_scenario: testCase = %s, featureClassName = %s, featureName = %s', testCase, featureClassName,
                  featureName)

    testCaseChanged = testUtils.setTestCase(testCase)
    featureChanged = testUtils.setFeatureClassName(featureClassName)

    global featureClass
    testImage = testUtils.getImage()
    testMask = testUtils.getMask()
    if featureClass is None or testCaseChanged or featureChanged:
      logging.debug('Init %s' % (featureClassName))
      featureClass = extractor.featureClasses[featureClassName](testImage, testMask, **kwargs)

    assert (featureClass != None)

    featureClass.disableAllFeatures()
    featureClass.enableFeatureByName(featureName)
    featureClass.calculateFeatures()
    # get the result and test it
    val = featureClass.featureValues[featureName]
    testUtils.checkResult(featureName, val)


def teardown_module():
  global testUtils
  global extractor
  global featureClass
  print("")
  res = testUtils.getResults()
  print 'Results:'
  print res
  resultsFile = os.path.join(testUtils.dataDir, 'PyradiomicsFeatures.csv')
  testUtils.writeCSV(res, resultsFile)
  diff = testUtils.getDiffs()
  print 'Differences from baseline:'
  print diff
  diffFile = os.path.join(testUtils.dataDir, 'Baseline2PyradiomicsFeaturesDiff.csv')
  testUtils.writeCSV(diff, diffFile)
  logging.info(
    "Wrote calculated features to %s, and the differences between the matlab features and the pyradiomics ones to %s.",
    resultsFile, diffFile)
  # Clean up
  testUtils = None
  extractor = None
  featureClass = None
