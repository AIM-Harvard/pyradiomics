# to run this test, from directory above:
# setenv PYTHONPATH /path/to/pyradiomics/radiomics
# nosetests --nocapture -v tests/rlglTest.py

import SimpleITK as sitk
from radiomics import firstorder, rlgl
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
        self.patientID = 'TCGA-02-0003_BrainMRI'

        dataDir = os.path.dirname(os.path.abspath(__file__)) + os.path.sep + ".." + os.path.sep + "data" + os.path.sep

        imageName = dataDir + self.patientID + os.path.sep + 'AXIAL-T1-POST-GD_TCGA-02-0003_TCGA-02-0003_PREOP_SCAN13_Background.nrrd'
        maskName = dataDir + self.patientID + os.path.sep + 'AXIAL-T1-POST-GD_TCGA-02-0003_TCGA-02-0003_PREOP_SCAN13_FSL-Clustered-labelMap_binarized_labelMap.nrrd'

        print("Reading the image and mask.")
        self.image = sitk.ReadImage(imageName)
        self.mask = sitk.ReadImage(maskName)

        print("Instantiating RLGL.")
        self.rlglFeatures = rlgl.RadiomicsRLGL(self.image, self.mask, 10)

        print("Reading expected results")
        baselineFileName = dataDir + 'MatlabFeatures.csv'
        self.foundPatientBaseline = False
        if (os.path.exists(baselineFileName)):
          self.baselineFeatures = {}
          csvFile = open(baselineFileName, 'rb')
          csvFileReader = csv.reader(csvFile)
          # get the column headers
          self.headerRow = csvFileReader.next()
          # search for the patient in the file
          for patientRow in csvFileReader:
            if patientRow[0] == self.patientID:
              print 'Found row for patient ',self.patientID
              # print ', '.join(patientRow)
              self.foundPatientBaseline = True
              columnIndex = 0
              for val in patientRow:
                self.baselineFeatures[self.headerRow[columnIndex]] = val
                columnIndex += 1
              break
            else:
              print 'Not the right patient ID:', patientRow[0]
          # print 'baselineFeatures for ', self.patientID, ' = ',self.baselineFeatures
          if self.foundPatientBaseline == False:
            print 'Unable to find the baseline for patient', patientID, ', conducting the test without evaluating the results.'
        else:
           print 'Unable to find baseline features file ',baselineFileName

    @classmethod
    def teardown_class(self):
        print ("") # this is to get a newline after the dots
        print ("teardown_class() after any methods in this class")

    def checkResult(self, key, value):
      if self.foundPatientBaseline == False:
        print 'Unable to evaluate calculated feature ', key, ' of value ', value
        return
      index = -1
      if key == 'ShortRunEmphasis':
        index = 'LoG_sigma_0_5_mm_2D_rlgl_shortRunEmphasis'
      elif key == 'LongRunEmphasis':
        index = 'LoG_sigma_0_5_mm_2D_rlgl_longRunEmphasis'
      elif key == 'GrayLevelNonUniformity':
          index = 'LoG_sigma_0_5_mm_2D_rlgl_grayLevelNonuniformity'
      elif key == 'RunLengthNonUniformity':
        index = 'LoG_sigma_0_5_mm_2D_rlgl_runLengthNonuniformity'
      elif key == 'RunPercentage':
        index = 'LoG_sigma_0_5_mm_2D_rlgl_runPercentage'
      # elif key == 'LowGrayLevelRunEmphasis':
        # missing from file
      #  index = 'LoG_sigma_0_5_mm_2D_rlgl_lowGrayLevelRunEmphasis'
      elif key == 'HighGrayLevelRunEmphasis':
        index = 'LoG_sigma_0_5_mm_2D_rlgl_highGrayLevelRunEmphasis'
      elif key == 'ShortRunLowGrayLevelEmphasis':
        index = 'LoG_sigma_0_5_mm_2D_rlgl_shortRunLowGrayLevEmpha'
      elif key == 'ShortRunHighGrayLevelEmphasis':
        index = 'LoG_sigma_0_5_mm_2D_rlgl_shortRunHighGrayLevEmpha'
      elif key == 'LongRunLowGrayLevelEmphasis':
        index = 'LoG_sigma_0_5_mm_2D_rlgl_longRunLowGrayLevEmpha'
      elif key == 'LongRunHighGrayLevelEmphasis':
        index = 'LoG_sigma_0_5_mm_2D_rlgl_longRunHighGrayLevEmpha'
      if index == -1:
        print 'Unable to find index for key ',key
        return
      baseline = self.baselineFeatures[index]
      percentDiff = abs(1.0 - (value / float(baseline)))
      print('index = %s, baseline value = %f, calculated = %f, diff = %f%%' % (index, float(baseline), value, percentDiff * 100))
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
