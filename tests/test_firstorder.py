# to run this test, from directory above:
# setenv PYTHONPATH /path/to/pyradiomics/radiomics
# nosetests --nocapture -v tests/test_firstorder.py

from radiomics import firstorder, imageoperations
from testUtils import RadiomicsTestUtils
import SimpleITK as sitk
import sys, os
import logging

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

    def test_energy(self):
        testString = 'Energy'
        logging.info('Test %s', testString)
        self.firstOrderFeatures.enableFeatureByName(testString)
        self.firstOrderFeatures.calculateFeatures()
        val = self.firstOrderFeatures.featureValues[testString]
        self.testUtils.checkResult(testString, val)

    def test_totalEnergy(self):
        testString = 'TotalEnergy'
        logging.info('Test %s', testString)
        self.firstOrderFeatures.enableFeatureByName(testString)
        self.firstOrderFeatures.calculateFeatures()
        val = self.firstOrderFeatures.featureValues[testString]
        self.testUtils.checkResult(testString, val)

    def test_entropy(self):
        testString = 'Entropy'
        logging.info('Test %s', testString)
        self.firstOrderFeatures.enableFeatureByName(testString)
        self.firstOrderFeatures.calculateFeatures()
        val = self.firstOrderFeatures.featureValues[testString]
        self.testUtils.checkResult(testString, val)

    def test_minimum(self):
        testString = 'Minimum'
        logging.info('Test %s', testString)
        self.firstOrderFeatures.enableFeatureByName(testString)
        self.firstOrderFeatures.calculateFeatures()
        val = self.firstOrderFeatures.featureValues[testString]
        self.testUtils.checkResult(testString, val)

    def test_maximum(self):
        testString = 'Maximum'
        logging.info('Test %s', testString)
        self.firstOrderFeatures.enableFeatureByName(testString)
        self.firstOrderFeatures.calculateFeatures()
        val = self.firstOrderFeatures.featureValues[testString]
        self.testUtils.checkResult(testString, val)

    def test_mean(self):
        testString = 'Mean'
        logging.info('Test %s', testString)
        self.firstOrderFeatures.enableFeatureByName(testString)
        self.firstOrderFeatures.calculateFeatures()
        val = self.firstOrderFeatures.featureValues[testString]
        self.testUtils.checkResult(testString, val)

    def test_median(self):
        testString = 'Median'
        logging.info('Test %s', testString)
        self.firstOrderFeatures.enableFeatureByName(testString)
        self.firstOrderFeatures.calculateFeatures()
        val = self.firstOrderFeatures.featureValues[testString]
        self.testUtils.checkResult(testString, val)

    def test_range(self):
        testString = 'Range'
        logging.info('Test %s', testString)
        self.firstOrderFeatures.enableFeatureByName(testString)
        self.firstOrderFeatures.calculateFeatures()
        val = self.firstOrderFeatures.featureValues[testString]
        self.testUtils.checkResult(testString, val)

    def test_meanDeviation(self):
        testString = 'MeanDeviation'
        logging.info('Test %s', testString)
        self.firstOrderFeatures.enableFeatureByName(testString)
        self.firstOrderFeatures.calculateFeatures()
        val = self.firstOrderFeatures.featureValues[testString]
        self.testUtils.checkResult(testString, val)

    def test_rootMeanSquared(self):
        testString = 'RootMeanSquared'
        logging.info('Test %s', testString)
        self.firstOrderFeatures.enableFeatureByName(testString)
        self.firstOrderFeatures.calculateFeatures()
        val = self.firstOrderFeatures.featureValues[testString]
        self.testUtils.checkResult(testString, val)

    def test_standardDeviation(self):
        testString = 'StandardDeviation'
        logging.info('Test %s', testString)
        self.firstOrderFeatures.enableFeatureByName(testString)
        self.firstOrderFeatures.calculateFeatures()
        val = self.firstOrderFeatures.featureValues[testString]
        self.testUtils.checkResult(testString, val)

    def test_skewness(self):
        testString = 'Skewness'
        logging.info('Test %s', testString)
        self.firstOrderFeatures.enableFeatureByName(testString)
        self.firstOrderFeatures.calculateFeatures()
        val = self.firstOrderFeatures.featureValues[testString]
        self.testUtils.checkResult(testString, val)

    def test_kurtosis(self):
        testString = 'Kurtosis'
        logging.info('Test %s', testString)
        self.firstOrderFeatures.enableFeatureByName(testString)
        self.firstOrderFeatures.calculateFeatures()
        val = self.firstOrderFeatures.featureValues[testString]
        self.testUtils.checkResult(testString, val)

    def test_variance(self):
        testString = 'Variance'
        logging.info('Test %s', testString)
        self.firstOrderFeatures.enableFeatureByName(testString)
        self.firstOrderFeatures.calculateFeatures()
        val = self.firstOrderFeatures.featureValues[testString]
        self.testUtils.checkResult(testString, val)

    def test_uniformity(self):
        testString = 'Uniformity'
        logging.info('Test %s', testString)
        self.firstOrderFeatures.enableFeatureByName(testString)
        self.firstOrderFeatures.calculateFeatures()
        val = self.firstOrderFeatures.featureValues[testString]
        self.testUtils.checkResult(testString, val)
