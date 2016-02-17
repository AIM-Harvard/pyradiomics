# to run this test, from directory above:
# setenv PYTHONPATH /path/to/pyradiomics/radiomics
# nosetests --nocapture -v tests/test_shape.py

import SimpleITK as sitk
from radiomics import firstorder, shape
from testUtils import RadiomicsTestUtils
import sys, os
import logging
from nose_parameterized import parameterized

def setup_module(module):
    # run before anything in this file
    print ("") # this is to get a newline after the dots
    return

class TestShape:

    def setup(self):
        # setup before each test method
        print ("") # this is to get a newline after the dots
        self.shapeFeatures.disableAllFeatures()

    @classmethod
    def setup_class(self):
        # run before any methods in this class
        print ("") # this is to get a newline after the dots

        # read in the baseline and mapping to matlab features
        self.testUtils = RadiomicsTestUtils('shape')
        # set the test case for these files to match the directory and
        # the id in the baseline file
        self.testUtils.setTestCase('breast1')

        dataDir = os.path.dirname(os.path.abspath(__file__)) + os.path.sep + ".." + os.path.sep + "data" + os.path.sep

        imageName = dataDir + self.testUtils.getTestCase() + '_image.nrrd'
        maskName = dataDir + self.testUtils.getTestCase() + '_label.nrrd'

        logging.info("Reading the image and mask.")
        self.image = sitk.ReadImage(imageName)
        self.mask = sitk.ReadImage(maskName)

        logging.info("Instantiating Shape.")
        self.shapeFeatures = shape.RadiomicsShape(self.image, self.mask)


    @classmethod
    def teardown_class(self):
        # run after any methods in this class
        print ("") # this is to get a newline after the dots

    def generate_scenarios():
      # get the feature names
      featureNames = shape.RadiomicsShape.getFeatureNames()
      logging.info('generate_scenarios: featureNames = %s', featureNames)
      for f in featureNames:
        yield (f)

    @parameterized.expand(generate_scenarios())
    def test_scenario_BreastMRI(self, featureName):
      logging.info('test_scenario: featureName = %s', featureName)
      self.shapeFeatures.enableFeatureByName(featureName)
      self.shapeFeatures.calculateFeatures()
      val = self.shapeFeatures.featureValues[featureName]
      self.testUtils.checkResult(featureName, val)
