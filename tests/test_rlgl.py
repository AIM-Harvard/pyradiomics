# to run this test, from directory above:
# setenv PYTHONPATH /path/to/pyradiomics/radiomics
# nosetests --nocapture -v tests/test_rlgl.py

import SimpleITK as sitk
from radiomics import firstorder, rlgl
from testUtils import RadiomicsTestUtils
import sys, os
import logging

def setup_module(module):
    # run before anything in this file"
    print ("") # this is to get a newline after the dots
    return

class TestRLGL:

    def setup(self):
        # setup before each test method
        print ("") # this is to get a newline after the dots
        # disabling all features
        self.rlglFeatures.disableAllFeatures()

    @classmethod
    def setup_class(self):
        # before any methods in this class"
        print ("") # this is to get a newline after the dots

        # read in the baseline and mapping to matlab features
        self.testUtils = RadiomicsTestUtils('rlgl')
        # set the test case for these files to match the directory and
        # the id in the baseline file
        self.testUtils.setTestCase('brain1')

        dataDir = os.path.dirname(os.path.abspath(__file__)) + os.path.sep + ".." + os.path.sep + "data" + os.path.sep

        imageName = dataDir + self.testUtils.getTestCase() + '_image.nrrd'
        maskName = dataDir + self.testUtils.getTestCase() + '_label.nrrd'

        logging.info("Reading the image and mask.")
        self.image = sitk.ReadImage(imageName)
        self.mask = sitk.ReadImage(maskName)

        logging.info("Instantiating RLGL.")
        self.rlglFeatures = rlgl.RadiomicsRLGL(self.image, self.mask)

    @classmethod
    def teardown_class(self):
        # after any methods in this class
        print ("") # this is to get a newline after the dots

    def test_shortRunEmphasis(self):
        testString = 'ShortRunEmphasis'
        logging.info('Test %s', testString)
        self.rlglFeatures.enableFeatureByName(testString)
        self.rlglFeatures.calculateFeatures()
        val = self.rlglFeatures.featureValues[testString]
        self.testUtils.checkResult(testString, val)

    def test_longRunEmphasis(self):
        testString = 'LongRunEmphasis'
        logging.info('Test %s', testString)
        self.rlglFeatures.enableFeatureByName(testString)
        self.rlglFeatures.calculateFeatures()
        val = self.rlglFeatures.featureValues[testString]
        self.testUtils.checkResult(testString, val)

    def test_grayLevelNonUniformity(self):
        testString = 'GrayLevelNonUniformity'
        logging.info('Test %s', testString)
        self.rlglFeatures.enableFeatureByName(testString)
        self.rlglFeatures.calculateFeatures()
        val = self.rlglFeatures.featureValues[testString]
        self.testUtils.checkResult(testString, val)

    def test_runLengthNonUniformity(self):
        testString = 'RunLengthNonUniformity'
        logging.info('Test %s', testString)
        self.rlglFeatures.enableFeatureByName(testString)
        self.rlglFeatures.calculateFeatures()
        val = self.rlglFeatures.featureValues[testString]
        self.testUtils.checkResult(testString, val)

    def test_runPercentage(self):
        testString = 'RunPercentage'
        logging.info('Test %s', testString)
        self.rlglFeatures.enableFeatureByName(testString)
        self.rlglFeatures.calculateFeatures()
        val = self.rlglFeatures.featureValues[testString]
        self.testUtils.checkResult(testString, val)

    def test_lowGrayLevelRunEmphasis(self):
        testString = 'LowGrayLevelRunEmphasis'
        logging.info('Test %s', testString)
        self.rlglFeatures.enableFeatureByName(testString)
        self.rlglFeatures.calculateFeatures()
        val = self.rlglFeatures.featureValues[testString]
        self.testUtils.checkResult(testString, val)

    def test_highGrayLevelRunEmphasis(self):
        testString = 'HighGrayLevelRunEmphasis'
        logging.info('Test %s', testString)
        self.rlglFeatures.enableFeatureByName(testString)
        self.rlglFeatures.calculateFeatures()
        val = self.rlglFeatures.featureValues[testString]
        self.testUtils.checkResult(testString, val)

    def test_shortRunHighGrayLevelEmphasis(self):
        testString = 'ShortRunHighGrayLevelEmphasis'
        logging.info('Test %s', testString)
        self.rlglFeatures.enableFeatureByName(testString)
        self.rlglFeatures.calculateFeatures()
        val = self.rlglFeatures.featureValues[testString]
        self.testUtils.checkResult(testString, val)

    def test_longRunLowGrayLevelEmphasis(self):
        testString = 'LongRunLowGrayLevelEmphasis'
        logging.info('Test %s', testString)
        self.rlglFeatures.enableFeatureByName(testString)
        self.rlglFeatures.calculateFeatures()
        val = self.rlglFeatures.featureValues[testString]
        self.testUtils.checkResult(testString, val)

    def test_longRunHighGrayLevelEmphasis(self):
        testString = 'LongRunHighGrayLevelEmphasis'
        logging.info('Test %s', testString)
        self.rlglFeatures.enableFeatureByName(testString)
        self.rlglFeatures.calculateFeatures()
        val = self.rlglFeatures.featureValues[testString]
        self.testUtils.checkResult(testString, val)
