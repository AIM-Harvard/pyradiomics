# to run this test, from directory above:
# setenv PYTHONPATH /path/to/pyradiomics/radiomics
# nosetests --nocapture -v tests/test_rlgl.py

import SimpleITK as sitk
from radiomics import firstorder, rlgl
from testUtils import RadiomicsTestUtils
import sys, os

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

        # set the patient ID for these files to match the directory and
        # the patient id in the baseline file
        self.patientID = 'brain1'

        # read in the baseline and mapping to matlab features
        self.testUtils = RadiomicsTestUtils('rlgl')
        self.testUtils.setPatientID(self.patientID)

        dataDir = os.path.dirname(os.path.abspath(__file__)) + os.path.sep + ".." + os.path.sep + "data" + os.path.sep

        imageName = dataDir + self.patientID + '_image.nrrd'
        maskName = dataDir + self.patientID + '_label.nrrd'

        print("Reading the image and mask.")
        self.image = sitk.ReadImage(imageName)
        self.mask = sitk.ReadImage(maskName)

        print("Instantiating RLGL.")
        self.rlglFeatures = rlgl.RadiomicsRLGL(self.image, self.mask)

    @classmethod
    def teardown_class(self):
        # after any methods in this class
        print ("") # this is to get a newline after the dots

    def test_shortRunEmphasis_10(self):
        testString = 'ShortRunEmphasis'
        if self.testUtils.getVerbose():
          print 'Test', testString
        self.rlglFeatures.enableFeatureByName(testString)
        self.rlglFeatures.calculateFeatures()
        val = self.rlglFeatures.featureValues[testString]
        self.testUtils.checkResult(testString, val)

    def test_longRunEmphasis_10(self):
        testString = 'LongRunEmphasis'
        if self.testUtils.getVerbose():
          print 'Test', testString
        self.rlglFeatures.enableFeatureByName(testString)
        self.rlglFeatures.calculateFeatures()
        val = self.rlglFeatures.featureValues[testString]
        self.testUtils.checkResult(testString, val)

    def test_grayLevelNonUniformity_10(self):
        testString = 'GrayLevelNonUniformity'
        if self.testUtils.getVerbose():
          print 'Test', testString
        self.rlglFeatures.enableFeatureByName(testString)
        self.rlglFeatures.calculateFeatures()
        val = self.rlglFeatures.featureValues[testString]
        self.testUtils.checkResult(testString, val)

    def test_runLengthNonUniformity_10(self):
        testString = 'RunLengthNonUniformity'
        if self.testUtils.getVerbose():
          print 'Test', testString
        self.rlglFeatures.enableFeatureByName(testString)
        self.rlglFeatures.calculateFeatures()
        val = self.rlglFeatures.featureValues[testString]
        self.testUtils.checkResult(testString, val)

    def test_runPercentage_10(self):
        testString = 'RunPercentage'
        if self.testUtils.getVerbose():
          print 'Test', testString
        self.rlglFeatures.enableFeatureByName(testString)
        self.rlglFeatures.calculateFeatures()
        val = self.rlglFeatures.featureValues[testString]
        self.testUtils.checkResult(testString, val)

    def test_lowGrayLevelRunEmphasis_10(self):
        testString = 'LowGrayLevelRunEmphasis'
        if self.testUtils.getVerbose():
          print 'Test', testString
        self.rlglFeatures.enableFeatureByName(testString)
        self.rlglFeatures.calculateFeatures()
        val = self.rlglFeatures.featureValues[testString]
        self.testUtils.checkResult(testString, val)

    def test_highGrayLevelRunEmphasis_10(self):
        testString = 'HighGrayLevelRunEmphasis'
        if self.testUtils.getVerbose():
          print 'Test', testString
        self.rlglFeatures.enableFeatureByName(testString)
        self.rlglFeatures.calculateFeatures()
        val = self.rlglFeatures.featureValues[testString]
        self.testUtils.checkResult(testString, val)

    def test_shortRunHighGrayLevelEmphasis_10(self):
        testString = 'ShortRunHighGrayLevelEmphasis'
        if self.testUtils.getVerbose():
          print 'Test', testString
        self.rlglFeatures.enableFeatureByName(testString)
        self.rlglFeatures.calculateFeatures()
        val = self.rlglFeatures.featureValues[testString]
        self.testUtils.checkResult(testString, val)

    def test_longRunLowGrayLevelEmphasis_10(self):
        testString = 'LongRunLowGrayLevelEmphasis'
        if self.testUtils.getVerbose():
          print 'Test', testString
        self.rlglFeatures.enableFeatureByName(testString)
        self.rlglFeatures.calculateFeatures()
        val = self.rlglFeatures.featureValues[testString]
        self.testUtils.checkResult(testString, val)

    def test_longRunHighGrayLevelEmphasis_10(self):
        testString = 'LongRunHighGrayLevelEmphasis'
        if self.testUtils.getVerbose():
          print 'Test', testString
        self.rlglFeatures.enableFeatureByName(testString)
        self.rlglFeatures.calculateFeatures()
        val = self.rlglFeatures.featureValues[testString]
        self.testUtils.checkResult(testString, val)
