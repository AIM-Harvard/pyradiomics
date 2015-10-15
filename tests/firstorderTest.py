# to run this test, from directory above:
# setenv PYTHONPATH /Users/nicole/Radiomics/github/pyradiomics/source:/Users/nicole/Slicer4-svn/Slicer-build-debug/SimpleITK-build/Wrapping
# setenv DYLD_LIBRARY_PATH /Users/nicole/Slicer4-svn/Slicer-build-debug/python-install/lib:/Users/nicole/Slicer4-svn/Slicer-build-debug/SimpleITK-build/lib
# nosetests --nocapture -v tests/firstorderTest.py

from radiomics import firstorder, preprocessing
import SimpleITK as sitk
import sys, os
import csv

def setup_module(module):
    print ("") # this is to get a newline after the dots
    print ("setup_module before anything in this file")
    return

class TestFirstOrder:

    def setup(self):
        print ("") # this is to get a newline after the dots
        # print ("setup before each test method, disabling all features")
        self.firstOrderFeatures.disableAllFeatures()

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

        print("Instantiating first order features.")
        self.firstOrderFeatures = firstorder.RadiomicsFirstOrder(self.image,self.mask)
        self.firstOrderFeatures.setBinWidth(10)

        print("Reading expected results")
        baselineFileName = dataDir + 'MatlabFeatures.csv'
        self.foundPatientBaseline = False
        if (os.path.exists(baselineFileName)):
          self.baselineFeatures = {}
          csvFile = open(baselineFileName, 'rb')
          csvFileReader = csv.reader(csvFile)
          # get the column headers
          self.headerRow = csvFileReader.next()
          # print self.headerRow
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
      if key == 'Energy':
        index = ''
      elif key == 'TotalEnergy':
        index = ''
      elif key == 'Entropy':
          index = ''
      elif key == 'MinIntensity':
        index = ''
      elif key == 'MaxIntensity':
        index = ''
      elif key == 'MeanIntensity':
        index = ''
      elif key == 'MedianIntensity':
        index = ''
      elif key == 'RangeIntensity':
        index = ''
      elif key == 'MeanDeviation':
        index = ''
      elif key == 'RootMeanSquared':
        index = ''
      elif key == 'StandardDeviation':
        index = ''
      elif key == 'SkewnessValue':
        index = ''
      elif key == 'Kurtosis':
        index = ''
      elif key == 'Variance':
        index = ''
      elif key == 'Uniformity':
        index = ''
      if index == -1:
        print 'Unable to find index for key',key
        return
      if index == '':
        print 'Undefined key',key,' Value =',value
        return
      baseline = self.baselineFeatures[index]
      diff = abs(float(baseline) - value)
      print 'index = ', index, ', baseline value = ', baseline, ', calculated = ', value, ', diff = ', diff
      assert(diff < 0.1)

    def test_energy_10(self):
        self.firstOrderFeatures.enableFeatureByName('Energy')
        print 'Will calculate the following first order features: '
        for f in self.firstOrderFeatures.enabledFeatures.keys():
            print '  ',f
            print eval('self.firstOrderFeatures.get'+f+'FeatureValue.__doc__')

        self.firstOrderFeatures.calculateFeatures()

        print 'Calculated first order features: '
        for (key,val) in self.firstOrderFeatures.featureValues.iteritems():
            print '  ',key,':',val
            if key == 'Energy':
                energy = val

        self.checkResult('Energy', energy)

    def test_totalEnergy_10(self):
        self.firstOrderFeatures.enableFeatureByName('TotalEnergy')
        self.firstOrderFeatures.calculateFeatures()
        print 'Calculated first order features: '
        for (key,val) in self.firstOrderFeatures.featureValues.iteritems():
            print '  ',key,':',val
            if key == 'TotalEnergy':
                totalEnergy = val

        self.checkResult('TotalEnergy', totalEnergy)

    def test_entropy_10(self):
        testString = 'Entropy'
        print 'Test', testString
        self.firstOrderFeatures.enableFeatureByName(testString)
        self.firstOrderFeatures.calculateFeatures()
        val = self.firstOrderFeatures.featureValues[testString]
        self.checkResult(testString, val)

    def test_minIntensity_10(self):
        testString = 'MinIntensity'
        print 'Test', testString
        self.firstOrderFeatures.enableFeatureByName(testString)
        self.firstOrderFeatures.calculateFeatures()
        val = self.firstOrderFeatures.featureValues[testString]
        self.checkResult(testString, val)

    def test_maxIntensity_10(self):
        testString = 'MaxIntensity'
        print 'Test', testString
        self.firstOrderFeatures.enableFeatureByName(testString)
        self.firstOrderFeatures.calculateFeatures()
        val = self.firstOrderFeatures.featureValues[testString]
        self.checkResult(testString, val)

    def test_meanIntensity_10(self):
        testString = 'MeanIntensity'
        print 'Test', testString
        self.firstOrderFeatures.enableFeatureByName(testString)
        self.firstOrderFeatures.calculateFeatures()
        val = self.firstOrderFeatures.featureValues[testString]
        self.checkResult(testString, val)

    def test_medianIntensity_10(self):
        testString = 'MedianIntensity'
        print 'Test', testString
        self.firstOrderFeatures.enableFeatureByName(testString)
        self.firstOrderFeatures.calculateFeatures()
        val = self.firstOrderFeatures.featureValues[testString]
        self.checkResult(testString, val)

    def test_rangeIntensity_10(self):
        testString = 'RangeIntensity'
        print 'Test', testString
        self.firstOrderFeatures.enableFeatureByName(testString)
        self.firstOrderFeatures.calculateFeatures()
        val = self.firstOrderFeatures.featureValues[testString]
        self.checkResult(testString, val)

    def test_meanDeviation_10(self):
        testString = 'MeanDeviation'
        print 'Test', testString
        self.firstOrderFeatures.enableFeatureByName(testString)
        self.firstOrderFeatures.calculateFeatures()
        val = self.firstOrderFeatures.featureValues[testString]
        self.checkResult(testString, val)

    def test_rootMeanSquared_10(self):
        testString = 'RootMeanSquared'
        print 'Test', testString
        self.firstOrderFeatures.enableFeatureByName(testString)
        self.firstOrderFeatures.calculateFeatures()
        val = self.firstOrderFeatures.featureValues[testString]
        self.checkResult(testString, val)

    def test_standardDeviation_10(self):
        testString = 'StandardDeviation'
        print 'Test', testString
        self.firstOrderFeatures.enableFeatureByName(testString)
        self.firstOrderFeatures.calculateFeatures()
        val = self.firstOrderFeatures.featureValues[testString]
        self.checkResult(testString, val)

    def test_skewnessValue_10(self):
        testString = 'SkewnessValue'
        print 'Test', testString
        self.firstOrderFeatures.enableFeatureByName(testString)
        self.firstOrderFeatures.calculateFeatures()
        val = self.firstOrderFeatures.featureValues[testString]
        self.checkResult(testString, val)

    def test_kurtosis_10(self):
        testString = 'Kurtosis'
        print 'Test', testString
        self.firstOrderFeatures.enableFeatureByName(testString)
        self.firstOrderFeatures.calculateFeatures()
        val = self.firstOrderFeatures.featureValues[testString]
        self.checkResult(testString, val)

    def test_variance_10(self):
        testString = 'Variance'
        print 'Test', testString
        self.firstOrderFeatures.enableFeatureByName(testString)
        self.firstOrderFeatures.calculateFeatures()
        val = self.firstOrderFeatures.featureValues[testString]
        self.checkResult(testString, val)

    def test_uniformity_10(self):
        testString = 'Uniformity'
        print 'Test', testString
        self.firstOrderFeatures.enableFeatureByName(testString)
        self.firstOrderFeatures.calculateFeatures()
        val = self.firstOrderFeatures.featureValues[testString]
        self.checkResult(testString, val)
