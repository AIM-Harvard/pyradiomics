# to run this test, from directory above:
# setenv PYTHONPATH /path/to/pyradiomics/radiomics
# nosetests --nocapture -v tests/test_firstorder.py

from radiomics import firstorder, imageoperations
from testUtils import RadiomicsTestUtils
import SimpleITK as sitk
import sys, os
import logging
from nose_parameterized import parameterized

def setup_module(module):
    # runs before anything in this file
    print ("") # this is to get a newline after the dots
    return

class TestFirstOrder:

    def setup(self):
        # setup before each test method
        print ("") # this is to get a newline after the dots
        self.firstOrderFeatures.disableAllFeatures()

    @classmethod
    def setup_class(self):
        # called before any methods in this class
        print ("") # this is to get a newline after the dots

        # read in the baseline and mapping to matlab features
        self.testUtils = RadiomicsTestUtils('firstorder')
        # set the test case for these files to match the directory and
        # the id in the baseline file
        self.testUtils.setTestCase('brain1')

        dataDir = os.path.dirname(os.path.abspath(__file__)) + os.path.sep + ".." + os.path.sep + "data" + os.path.sep

        imageName = dataDir + self.testUtils.getTestCase() + '_image.nrrd'
        maskName = dataDir + self.testUtils.getTestCase() + '_label.nrrd'

        logging.info("Reading the image and mask.")
        self.image = sitk.ReadImage(imageName)
        self.mask = sitk.ReadImage(maskName)

        logging.info("Instantiating first order features.")
        self.firstOrderFeatures = firstorder.RadiomicsFirstOrder(self.image,self.mask)

    @classmethod
    def teardown_class(self):
        # run after any methods in this class
        print ("") # this is to get a newline after the dots

    def generate_scenarios():
      # get the feature names
      featureNames = firstorder.RadiomicsFirstOrder.getFeatureNames()
      logging.info('generate_scenarios: featureNames = %s', featureNames)
      for f in featureNames:
        yield (f)

    @parameterized.expand(generate_scenarios())
    def test_scenario(self, featureName):
      logging.info('test_scenario: featureName = %s', featureName)
      self.firstOrderFeatures.enableFeatureByName(featureName)
      self.firstOrderFeatures.calculateFeatures()
      val = self.firstOrderFeatures.featureValues[featureName]
      self.testUtils.checkResult(featureName, val)
