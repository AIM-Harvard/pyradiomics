# to run this test, from directory above:
# setenv PYTHONPATH /path/to/pyradiomics/radiomics
# nosetests --nocapture -v tests/test_firstorder.py

from radiomics import firstorder
from testUtils import RadiomicsTestUtils
import sys, os
import logging
from nose_parameterized import parameterized

testUtils = RadiomicsTestUtils('firstorder')
firstOrderFeatures = None

def setup_module(module):
    # runs before anything in this file
    print ("") # this is to get a newline after the dots
    return

class TestFirstOrder:

    #
    # Generate the test cases that test each feature for each test case / data
    # set.
    # The custom function name generation utility ensures that all the features
    # are tested for one data set before moving on to calculating features for
    # another data set.
    #
    def generate_scenarios():
      global testUtils
      # get the list of test cases for which we have baseline information
      testCases = testUtils.getTestCases()
      logging.info('generate_scenarios: testCases = %s', testCases)
      # get the list of features we'll be testing so we can generate the
      # tuples of testcase:feature to iterate over
      featureNames = firstorder.RadiomicsFirstOrder.getFeatureNames()
      logging.info('generate_scenarios: featureNames = %s', featureNames)
      test_tuples = []
      for t in range(len(testCases)):
        for f in range(len(featureNames)):
          test_tuples.append((testCases[t], featureNames[f]))
      logging.info('generate_scenarios: test_tuples = %s', test_tuples)

      for p in range(len(test_tuples)):
        yield (test_tuples[p][0], test_tuples[p][1])

    global testUtils
    @parameterized.expand(generate_scenarios(), testcase_func_name=testUtils.custom_name_func)
    def test_scenario(self, testCase, featureName):
      print("")
      global testUtils
      logging.info('test_scenario: testCase = %s, featureName = %s', testCase, featureName)
      # set the test case and only recalculate features if changed
      testCaseChanged = testUtils.setTestCase(testCase)
      global firstOrderFeatures
      if firstOrderFeatures == None or testCaseChanged:
        logging.info("Instantiating FirstOrder for testCase %s.", testCase)
        firstOrderFeatures = firstorder.RadiomicsFirstOrder(testUtils.getImage(), testUtils.getMask())

      firstOrderFeatures.disableAllFeatures()
      firstOrderFeatures.enableFeatureByName(featureName)
      firstOrderFeatures.calculateFeatures()
      # get the result and test it
      val = firstOrderFeatures.featureValues[featureName]
      testUtils.checkResult(featureName, val)


