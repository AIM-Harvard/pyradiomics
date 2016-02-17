# to run this test, from directory above:
# setenv PYTHONPATH /path/to/pyradiomics/radiomics
# nosetests --nocapture -v tests/test_glcm.py

from radiomics import firstorder, glcm, imageoperations
from testUtils import RadiomicsTestUtils
import SimpleITK as sitk
import sys, os
import logging

def setup_module(module):
    # run before anything in this file
    print ("") # this is to get a newline after the dots
    return

class TestGLCM:

    def setup(self):
        # setup before each test method
        print ("") # this is to get a newline after the dots
        self.glcmFeatures.disableAllFeatures()

    @classmethod
    def setup_class(self):
        # run before any methods in this class
        print ("") # this is to get a newline after the dots

        # read in the baseline and mapping to matlab features
        self.testUtils = RadiomicsTestUtils('glcm')
        # set the test case for these files to match the directory and
        # the id in the baseline file
        self.testUtils.setTestCase('brain1')

        dataDir = os.path.dirname(os.path.abspath(__file__)) + os.path.sep + ".." + os.path.sep + "data" + os.path.sep

        imageName = dataDir + self.testUtils.getTestCase() + '_image.nrrd'
        maskName = dataDir + self.testUtils.getTestCase() + '_label.nrrd'

        logging.info("Reading the image and mask.")
        self.image = sitk.ReadImage(imageName)
        self.mask = sitk.ReadImage(maskName)

        logging.info("Instantiating GLCM.")
        self.glcmFeatures = glcm.RadiomicsGLCM(self.image, self.mask)

    @classmethod
    def teardown_class(self):
        # run after any methods in this class
        print ("") # this is to get a newline after the dots

    def test_autocorrelation(self):
        testString = 'Autocorrelation'
        logging.info('Test %s', testString)
        self.glcmFeatures.enableFeatureByName(testString)
        self.glcmFeatures.calculateFeatures()
        val = self.glcmFeatures.featureValues[testString]
        self.testUtils.checkResult(testString, val)

    def test_clusterProminence(self):
        testString = 'ClusterProminence'
        logging.info('Test %s', testString)
        self.glcmFeatures.enableFeatureByName(testString)
        self.glcmFeatures.calculateFeatures()
        val = self.glcmFeatures.featureValues[testString]
        self.testUtils.checkResult(testString, val)

    def test_clusterShade(self):
        testString = 'ClusterShade'
        logging.info('Test %s', testString)
        self.glcmFeatures.enableFeatureByName(testString)
        self.glcmFeatures.calculateFeatures()
        val = self.glcmFeatures.featureValues[testString]
        self.testUtils.checkResult(testString, val)

    def test_clusterTendency(self):
        testString = 'ClusterTendency'
        logging.info('Test %s', testString)
        self.glcmFeatures.enableFeatureByName(testString)
        self.glcmFeatures.calculateFeatures()
        val = self.glcmFeatures.featureValues[testString]
        self.testUtils.checkResult(testString, val)

    def test_contrast(self):
        testString = 'Contrast'
        logging.info('Test %s', testString)
        self.glcmFeatures.enableFeatureByName(testString)
        self.glcmFeatures.calculateFeatures()
        val = self.glcmFeatures.featureValues[testString]
        self.testUtils.checkResult(testString, val)

    def test_correlation(self):
        testString = 'Correlation'
        logging.info('Test %s', testString)
        self.glcmFeatures.enableFeatureByName(testString)
        self.glcmFeatures.calculateFeatures()
        val = self.glcmFeatures.featureValues[testString]
        self.testUtils.checkResult(testString, val)

    def test_differenceEntropy(self):
        testString = 'DifferenceEntropy'
        logging.info('Test %s', testString)
        self.glcmFeatures.enableFeatureByName(testString)
        self.glcmFeatures.calculateFeatures()
        val = self.glcmFeatures.featureValues[testString]
        self.testUtils.checkResult(testString, val)

    def test_dissimilarity(self):
        testString = 'Dissimilarity'
        logging.info('Test %s', testString)
        self.glcmFeatures.enableFeatureByName(testString)
        self.glcmFeatures.calculateFeatures()
        val = self.glcmFeatures.featureValues[testString]
        self.testUtils.checkResult(testString, val)

    def test_energy(self):
        testString = 'Energy'
        logging.info('Test %s', testString)
        self.glcmFeatures.enableFeatureByName(testString)
        self.glcmFeatures.calculateFeatures()
        val = self.glcmFeatures.featureValues[testString]
        self.testUtils.checkResult(testString, val)

    def test_entropy(self):
        testString = 'Entropy'
        logging.info('Test %s', testString)
        self.glcmFeatures.enableFeatureByName(testString)
        self.glcmFeatures.calculateFeatures()
        val = self.glcmFeatures.featureValues[testString]
        self.testUtils.checkResult(testString, val)

    def test_homogeneity1(self):
        testString = 'Homogeneity1'
        logging.info('Test %s', testString)
        self.glcmFeatures.enableFeatureByName(testString)
        self.glcmFeatures.calculateFeatures()
        val = self.glcmFeatures.featureValues[testString]
        self.testUtils.checkResult(testString, val)

    def test_homogeneity2(self):
        testString = 'Homogeneity2'
        logging.info('Test %s', testString)
        self.glcmFeatures.enableFeatureByName(testString)
        self.glcmFeatures.calculateFeatures()
        val = self.glcmFeatures.featureValues[testString]
        self.testUtils.checkResult(testString, val)

    def test_imc1(self):
        testString = 'Imc1'
        logging.info('Test %s', testString)
        self.glcmFeatures.enableFeatureByName(testString)
        self.glcmFeatures.calculateFeatures()
        val = self.glcmFeatures.featureValues[testString]
        self.testUtils.checkResult(testString, val)

    def test_imc2(self):
        testString = 'Imc2'
        logging.info('Test %s', testString)
        self.glcmFeatures.enableFeatureByName(testString)
        self.glcmFeatures.calculateFeatures()
        val = self.glcmFeatures.featureValues[testString]
        self.testUtils.checkResult(testString, val)

    def test_idmn(self):
        testString = 'Idmn'
        logging.info('Test %s', testString)
        self.glcmFeatures.enableFeatureByName(testString)
        self.glcmFeatures.calculateFeatures()
        val = self.glcmFeatures.featureValues[testString]
        self.testUtils.checkResult(testString, val)

    def test_idn(self):
        testString = 'Idn'
        logging.info('Test %s', testString)
        self.glcmFeatures.enableFeatureByName(testString)
        self.glcmFeatures.calculateFeatures()
        val = self.glcmFeatures.featureValues[testString]
        self.testUtils.checkResult(testString, val)

    def test_inverseVariance(self):
        testString = 'InverseVariance'
        logging.info('Test %s', testString)
        self.glcmFeatures.enableFeatureByName(testString)
        self.glcmFeatures.calculateFeatures()
        val = self.glcmFeatures.featureValues[testString]
        self.testUtils.checkResult(testString, val)

    def test_maximumProbability(self):
        testString = 'MaximumProbability'
        logging.info('Test %s', testString)
        self.glcmFeatures.enableFeatureByName(testString)
        self.glcmFeatures.calculateFeatures()
        val = self.glcmFeatures.featureValues[testString]
        self.testUtils.checkResult(testString, val)

    def test_sumAverage(self):
        testString = 'SumAverage'
        logging.info('Test %s', testString)
        self.glcmFeatures.enableFeatureByName(testString)
        self.glcmFeatures.calculateFeatures()
        val = self.glcmFeatures.featureValues[testString]
        self.testUtils.checkResult(testString, val)

    def test_sumVariance(self):
        testString = 'SumVariance'
        logging.info('Test %s', testString)
        self.glcmFeatures.enableFeatureByName(testString)
        self.glcmFeatures.calculateFeatures()
        val = self.glcmFeatures.featureValues[testString]
        self.testUtils.checkResult(testString, val)

    def test_sumSquares(self):
        testString = 'SumSquares'
        logging.info('Test %s', testString)
        self.glcmFeatures.enableFeatureByName(testString)
        self.glcmFeatures.calculateFeatures()
        val = self.glcmFeatures.featureValues[testString]
        self.testUtils.checkResult(testString, val)
