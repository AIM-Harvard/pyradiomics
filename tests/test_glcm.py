# to run this test, from directory above:
# setenv PYTHONPATH /path/to/pyradiomics/radiomics
# nosetests --nocapture -v tests/test_glcm.py

from radiomics import firstorder, glcm, imageoperations
from testUtils import RadiomicsTestUtils
import SimpleITK as sitk
import sys, os

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

        # set the patient ID for these files to match the directory and
        # the patient id in the baseline file
        self.patientID = 'brain1'

        # read in the baseline and mapping to matlab features
        self.testUtils = RadiomicsTestUtils('glcm')
        self.testUtils.setPatientID(self.patientID)

        dataDir = os.path.dirname(os.path.abspath(__file__)) + os.path.sep + ".." + os.path.sep + "data" + os.path.sep

        imageName = dataDir + self.patientID + '_image.nrrd'
        maskName = dataDir + self.patientID + '_label.nrrd'

        print("Reading the image and mask.")
        self.image = sitk.ReadImage(imageName)
        self.mask = sitk.ReadImage(maskName)

        print("Instantiating GLCM.")
        self.glcmFeatures = glcm.RadiomicsGLCM(self.image, self.mask)

    @classmethod
    def teardown_class(self):
        # run after any methods in this class
        print ("") # this is to get a newline after the dots

    def test_autocorrelation_10(self):
        testString = 'Autocorrelation'
        if self.testUtils.getVerbose():
          print 'Test', testString
        self.glcmFeatures.enableFeatureByName(testString)
        self.glcmFeatures.calculateFeatures()
        val = self.glcmFeatures.featureValues[testString]
        self.testUtils.checkResult(testString, val)

    def test_clusterProminence_10(self):
        testString = 'ClusterProminence'
        if self.testUtils.getVerbose():
          print 'Test', testString
        self.glcmFeatures.enableFeatureByName(testString)
        self.glcmFeatures.calculateFeatures()
        val = self.glcmFeatures.featureValues[testString]
        self.testUtils.checkResult(testString, val)

    def test_clusterShade_10(self):
        testString = 'ClusterShade'
        if self.testUtils.getVerbose():
          print 'Test', testString
        self.glcmFeatures.enableFeatureByName(testString)
        self.glcmFeatures.calculateFeatures()
        val = self.glcmFeatures.featureValues[testString]
        self.testUtils.checkResult(testString, val)

    def test_clusterTendency_10(self):
        testString = 'ClusterTendency'
        if self.testUtils.getVerbose():
          print 'Test', testString
        self.glcmFeatures.enableFeatureByName(testString)
        self.glcmFeatures.calculateFeatures()
        val = self.glcmFeatures.featureValues[testString]
        self.testUtils.checkResult(testString, val)

    def test_contrast_10(self):
        testString = 'Contrast'
        if self.testUtils.getVerbose():
          print 'Test', testString
        self.glcmFeatures.enableFeatureByName(testString)
        self.glcmFeatures.calculateFeatures()
        val = self.glcmFeatures.featureValues[testString]
        self.testUtils.checkResult(testString, val)

    def test_correlation_10(self):
        testString = 'Correlation'
        if self.testUtils.getVerbose():
          print 'Test', testString
        self.glcmFeatures.enableFeatureByName(testString)
        self.glcmFeatures.calculateFeatures()
        val = self.glcmFeatures.featureValues[testString]
        self.testUtils.checkResult(testString, val)

    def test_differenceEntropy_10(self):
        testString = 'DifferenceEntropy'
        if self.testUtils.getVerbose():
          print 'Test', testString
        self.glcmFeatures.enableFeatureByName(testString)
        self.glcmFeatures.calculateFeatures()
        val = self.glcmFeatures.featureValues[testString]
        self.testUtils.checkResult(testString, val)

    def test_dissimilarity_10(self):
        testString = 'Dissimilarity'
        if self.testUtils.getVerbose():
          print 'Test', testString
        self.glcmFeatures.enableFeatureByName(testString)
        self.glcmFeatures.calculateFeatures()
        val = self.glcmFeatures.featureValues[testString]
        self.testUtils.checkResult(testString, val)

    def test_energy_10(self):
        testString = 'Energy'
        if self.testUtils.getVerbose():
          print 'Test', testString
        self.glcmFeatures.enableFeatureByName(testString)
        self.glcmFeatures.calculateFeatures()
        val = self.glcmFeatures.featureValues[testString]
        self.testUtils.checkResult(testString, val)

    def test_entropy_10(self):
        testString = 'Entropy'
        if self.testUtils.getVerbose():
          print 'Test', testString
        self.glcmFeatures.enableFeatureByName(testString)
        self.glcmFeatures.calculateFeatures()
        val = self.glcmFeatures.featureValues[testString]
        self.testUtils.checkResult(testString, val)

    def test_homogeneity1_10(self):
        testString = 'Homogeneity1'
        if self.testUtils.getVerbose():
          print 'Test', testString
        self.glcmFeatures.enableFeatureByName(testString)
        self.glcmFeatures.calculateFeatures()
        val = self.glcmFeatures.featureValues[testString]
        self.testUtils.checkResult(testString, val)

    def test_homogeneity2_10(self):
        testString = 'Homogeneity2'
        if self.testUtils.getVerbose():
          print 'Test', testString
        self.glcmFeatures.enableFeatureByName(testString)
        self.glcmFeatures.calculateFeatures()
        val = self.glcmFeatures.featureValues[testString]
        self.testUtils.checkResult(testString, val)

    def test_imc1_10(self):
        testString = 'Imc1'
        if self.testUtils.getVerbose():
          print 'Test', testString
        self.glcmFeatures.enableFeatureByName(testString)
        self.glcmFeatures.calculateFeatures()
        val = self.glcmFeatures.featureValues[testString]
        self.testUtils.checkResult(testString, val)

    def test_imc2_10(self):
        testString = 'Imc2'
        if self.testUtils.getVerbose():
          print 'Test', testString
        self.glcmFeatures.enableFeatureByName(testString)
        self.glcmFeatures.calculateFeatures()
        val = self.glcmFeatures.featureValues[testString]
        self.testUtils.checkResult(testString, val)

    def test_idmn_10(self):
        testString = 'Idmn'
        if self.testUtils.getVerbose():
          print 'Test', testString
        self.glcmFeatures.enableFeatureByName(testString)
        self.glcmFeatures.calculateFeatures()
        val = self.glcmFeatures.featureValues[testString]
        self.testUtils.checkResult(testString, val)

    def test_idn_10(self):
        testString = 'Idn'
        if self.testUtils.getVerbose():
          print 'Test', testString
        self.glcmFeatures.enableFeatureByName(testString)
        self.glcmFeatures.calculateFeatures()
        val = self.glcmFeatures.featureValues[testString]
        self.testUtils.checkResult(testString, val)

    def test_inverseVariance_10(self):
        testString = 'InverseVariance'
        if self.testUtils.getVerbose():
          print 'Test', testString
        self.glcmFeatures.enableFeatureByName(testString)
        self.glcmFeatures.calculateFeatures()
        val = self.glcmFeatures.featureValues[testString]
        self.testUtils.checkResult(testString, val)

    def test_maximumProbability_10(self):
        testString = 'MaximumProbability'
        if self.testUtils.getVerbose():
          print 'Test', testString
        self.glcmFeatures.enableFeatureByName(testString)
        self.glcmFeatures.calculateFeatures()
        val = self.glcmFeatures.featureValues[testString]
        self.testUtils.checkResult(testString, val)

    def test_sumAverage_10(self):
        testString = 'SumAverage'
        if self.testUtils.getVerbose():
          print 'Test', testString
        self.glcmFeatures.enableFeatureByName(testString)
        self.glcmFeatures.calculateFeatures()
        val = self.glcmFeatures.featureValues[testString]
        self.testUtils.checkResult(testString, val)

    def test_sumVariance_10(self):
        testString = 'SumVariance'
        if self.testUtils.getVerbose():
          print 'Test', testString
        self.glcmFeatures.enableFeatureByName(testString)
        self.glcmFeatures.calculateFeatures()
        val = self.glcmFeatures.featureValues[testString]
        self.testUtils.checkResult(testString, val)

    def test_sumSquares_10(self):
        testString = 'SumSquares'
        if self.testUtils.getVerbose():
          print 'Test', testString
        self.glcmFeatures.enableFeatureByName(testString)
        self.glcmFeatures.calculateFeatures()
        val = self.glcmFeatures.featureValues[testString]
        self.testUtils.checkResult(testString, val)
