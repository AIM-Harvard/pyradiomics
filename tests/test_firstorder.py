# to run this test, from directory above:
# setenv PYTHONPATH /path/to/pyradiomics/radiomics
# nosetests --nocapture -v tests/test_firstorder.py

from radiomics import firstorder, imageoperations
from testUtils import RadiomicsTestUtils
import SimpleITK as sitk
import sys, os

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

        # set the patient ID for these files to match the directory and
        # the patient id in the baseline file
        self.patientID = 'brain1'

        # read in the baseline and mapping to matlab features
        self.testUtils = RadiomicsTestUtils('firstorder')
        self.testUtils.setPatientID(self.patientID)

        dataDir = os.path.dirname(os.path.abspath(__file__)) + os.path.sep + ".." + os.path.sep + "data" + os.path.sep

        imageName = dataDir + self.patientID + '_image.nrrd'
        maskName = dataDir + self.patientID + '_label.nrrd'

        print("Reading the image and mask.")
        self.image = sitk.ReadImage(imageName)
        self.mask = sitk.ReadImage(maskName)

        print("Instantiating first order features.")
        self.firstOrderFeatures = firstorder.RadiomicsFirstOrder(self.image,self.mask)

    @classmethod
    def teardown_class(self):
        # run after any methods in this class
        print ("") # this is to get a newline after the dots

    def test_energy_10(self):
        testString = 'Energy'
        if self.testUtils.getVerbose():
          print 'Test', testString
        self.firstOrderFeatures.enableFeatureByName(testString)
        self.firstOrderFeatures.calculateFeatures()
        val = self.firstOrderFeatures.featureValues[testString]
        self.testUtils.checkResult(testString, val)

    def test_totalEnergy_10(self):
        testString = 'TotalEnergy'
        if self.testUtils.getVerbose():
          print 'Test', testString
        self.firstOrderFeatures.enableFeatureByName(testString)
        self.firstOrderFeatures.calculateFeatures()
        val = self.firstOrderFeatures.featureValues[testString]
        self.testUtils.checkResult(testString, val)

    def test_entropy_10(self):
        testString = 'Entropy'
        if self.testUtils.getVerbose():
          print 'Test', testString
        self.firstOrderFeatures.enableFeatureByName(testString)
        self.firstOrderFeatures.calculateFeatures()
        val = self.firstOrderFeatures.featureValues[testString]
        self.testUtils.checkResult(testString, val)

    def test_minimum_10(self):
        testString = 'Minimum'
        if self.testUtils.getVerbose():
          print 'Test', testString
        self.firstOrderFeatures.enableFeatureByName(testString)
        self.firstOrderFeatures.calculateFeatures()
        val = self.firstOrderFeatures.featureValues[testString]
        self.testUtils.checkResult(testString, val)

    def test_maximum_10(self):
        testString = 'Maximum'
        if self.testUtils.getVerbose():
          print 'Test', testString
        self.firstOrderFeatures.enableFeatureByName(testString)
        self.firstOrderFeatures.calculateFeatures()
        val = self.firstOrderFeatures.featureValues[testString]
        self.testUtils.checkResult(testString, val)

    def test_mean_10(self):
        testString = 'Mean'
        if self.testUtils.getVerbose():
          print 'Test', testString
        self.firstOrderFeatures.enableFeatureByName(testString)
        self.firstOrderFeatures.calculateFeatures()
        val = self.firstOrderFeatures.featureValues[testString]
        self.testUtils.checkResult(testString, val)

    def test_median_10(self):
        testString = 'Median'
        if self.testUtils.getVerbose():
          print 'Test', testString
        self.firstOrderFeatures.enableFeatureByName(testString)
        self.firstOrderFeatures.calculateFeatures()
        val = self.firstOrderFeatures.featureValues[testString]
        self.testUtils.checkResult(testString, val)

    def test_range_10(self):
        testString = 'Range'
        if self.testUtils.getVerbose():
          print 'Test', testString
        self.firstOrderFeatures.enableFeatureByName(testString)
        self.firstOrderFeatures.calculateFeatures()
        val = self.firstOrderFeatures.featureValues[testString]
        self.testUtils.checkResult(testString, val)

    def test_meanDeviation_10(self):
        testString = 'MeanDeviation'
        if self.testUtils.getVerbose():
          print 'Test', testString
        self.firstOrderFeatures.enableFeatureByName(testString)
        self.firstOrderFeatures.calculateFeatures()
        val = self.firstOrderFeatures.featureValues[testString]
        self.testUtils.checkResult(testString, val)

    def test_rootMeanSquared_10(self):
        testString = 'RootMeanSquared'
        if self.testUtils.getVerbose():
          print 'Test', testString
        self.firstOrderFeatures.enableFeatureByName(testString)
        self.firstOrderFeatures.calculateFeatures()
        val = self.firstOrderFeatures.featureValues[testString]
        self.testUtils.checkResult(testString, val)

    def test_standardDeviation_10(self):
        testString = 'StandardDeviation'
        if self.testUtils.getVerbose():
          print 'Test', testString
        self.firstOrderFeatures.enableFeatureByName(testString)
        self.firstOrderFeatures.calculateFeatures()
        val = self.firstOrderFeatures.featureValues[testString]
        self.testUtils.checkResult(testString, val)

    def test_skewness_10(self):
        testString = 'Skewness'
        if self.testUtils.getVerbose():
          print 'Test', testString
        self.firstOrderFeatures.enableFeatureByName(testString)
        self.firstOrderFeatures.calculateFeatures()
        val = self.firstOrderFeatures.featureValues[testString]
        self.testUtils.checkResult(testString, val)

    def test_kurtosis_10(self):
        testString = 'Kurtosis'
        if self.testUtils.getVerbose():
          print 'Test', testString
        self.firstOrderFeatures.enableFeatureByName(testString)
        self.firstOrderFeatures.calculateFeatures()
        val = self.firstOrderFeatures.featureValues[testString]
        self.testUtils.checkResult(testString, val)

    def test_variance_10(self):
        testString = 'Variance'
        if self.testUtils.getVerbose():
          print 'Test', testString
        self.firstOrderFeatures.enableFeatureByName(testString)
        self.firstOrderFeatures.calculateFeatures()
        val = self.firstOrderFeatures.featureValues[testString]
        self.testUtils.checkResult(testString, val)

    def test_uniformity_10(self):
        testString = 'Uniformity'
        if self.testUtils.getVerbose():
          print 'Test', testString
        self.firstOrderFeatures.enableFeatureByName(testString)
        self.firstOrderFeatures.calculateFeatures()
        val = self.firstOrderFeatures.featureValues[testString]
        self.testUtils.checkResult(testString, val)
