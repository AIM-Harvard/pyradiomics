# to run this test, from directory above:
# setenv PYTHONPATH /path/to/pyradiomics/radiomics
# nosetests --nocapture -v tests/rlglTest.py

import SimpleITK as sitk
from radiomics import firstorder, rlgl
from testUtils import RadiomicsTestUtils
import sys, os
import csv

def setup_module(module):
    print ("") # this is to get a newline after the dots
    print ("setup_module before anything in this file")
    return

class TestRLGL:

    def setup(self):
        print ("") # this is to get a newline after the dots
        # print ("setup before each test method, disabling all features")
        self.rlglFeatures.disableAllFeatures()

    @classmethod
    def setup_class(self):
        print ("") # this is to get a newline after the dots
        print ("setup_class() before any methods in this class")

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

        print("Instantiating RLGL.")
        self.rlglFeatures = rlgl.RadiomicsRLGL(self.image, self.mask)

    @classmethod
    def teardown_class(self):
        print ("") # this is to get a newline after the dots
        print ("teardown_class() after any methods in this class")

    def checkResult(self, key, value):
      # use the mapping from the utils
      baseline = self.testUtils.getMatlabValue(key)

      percentDiff = abs(1.0 - (value / float(baseline)))
      print('baseline value = %f, calculated = %f, diff = %f%%' % (float(baseline), value, percentDiff * 100))
      # check for a less than one percent difference
      assert(percentDiff < 0.01)

    def test_shortRunEmphasis_10(self):
        testString = 'ShortRunEmphasis'
        print 'Test', testString
        self.rlglFeatures.enableFeatureByName(testString)
        self.rlglFeatures.calculateFeatures()
        val = self.rlglFeatures.featureValues[testString]
        self.checkResult(testString, val)

    def test_longRunEmphasis_10(self):
        testString = 'LongRunEmphasis'
        print 'Test', testString
        self.rlglFeatures.enableFeatureByName(testString)
        self.rlglFeatures.calculateFeatures()
        val = self.rlglFeatures.featureValues[testString]
        self.checkResult(testString, val)

    def test_grayLevelNonUniformity_10(self):
        testString = 'GrayLevelNonUniformity'
        print 'Test', testString
        self.rlglFeatures.enableFeatureByName(testString)
        self.rlglFeatures.calculateFeatures()
        val = self.rlglFeatures.featureValues[testString]
        self.checkResult(testString, val)

    def test_runLengthNonUniformity_10(self):
        testString = 'RunLengthNonUniformity'
        print 'Test', testString
        self.rlglFeatures.enableFeatureByName(testString)
        self.rlglFeatures.calculateFeatures()
        val = self.rlglFeatures.featureValues[testString]
        self.checkResult(testString, val)

    def test_runPercentage_10(self):
        testString = 'RunPercentage'
        print 'Test', testString
        self.rlglFeatures.enableFeatureByName(testString)
        self.rlglFeatures.calculateFeatures()
        val = self.rlglFeatures.featureValues[testString]
        self.checkResult(testString, val)

    def test_lowGrayLevelRunEmphasis_10(self):
        testString = 'LowGrayLevelRunEmphasis'
        print 'Test', testString
        self.rlglFeatures.enableFeatureByName(testString)
        self.rlglFeatures.calculateFeatures()
        val = self.rlglFeatures.featureValues[testString]
        self.checkResult(testString, val)

    def test_highGrayLevelRunEmphasis_10(self):
        testString = 'HighGrayLevelRunEmphasis'
        print 'Test', testString
        self.rlglFeatures.enableFeatureByName(testString)
        self.rlglFeatures.calculateFeatures()
        val = self.rlglFeatures.featureValues[testString]
        self.checkResult(testString, val)

    def test_shortRunHighGrayLevelEmphasis_10(self):
        testString = 'ShortRunHighGrayLevelEmphasis'
        print 'Test', testString
        self.rlglFeatures.enableFeatureByName(testString)
        self.rlglFeatures.calculateFeatures()
        val = self.rlglFeatures.featureValues[testString]
        self.checkResult(testString, val)

    def test_longRunLowGrayLevelEmphasis_10(self):
        testString = 'LongRunLowGrayLevelEmphasis'
        print 'Test', testString
        self.rlglFeatures.enableFeatureByName(testString)
        self.rlglFeatures.calculateFeatures()
        val = self.rlglFeatures.featureValues[testString]
        self.checkResult(testString, val)

    def test_longRunHighGrayLevelEmphasis_10(self):
        testString = 'LongRunHighGrayLevelEmphasis'
        print 'Test', testString
        self.rlglFeatures.enableFeatureByName(testString)
        self.rlglFeatures.calculateFeatures()
        val = self.rlglFeatures.featureValues[testString]
        self.checkResult(testString, val)
