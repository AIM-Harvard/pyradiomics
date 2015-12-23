# to run this test, from directory above:
# setenv PYTHONPATH /path/to/pyradiomics/radiomics
# nosetests --nocapture -v tests/shapeTest.py

import SimpleITK as sitk
from radiomics import firstorder, shape
import sys, os
import csv

def setup_module(module):
    print ("") # this is to get a newline after the dots
    print ("setup_module before anything in this file")
    return

class TestShape:

    def setup(self):
        print ("") # this is to get a newline after the dots
        # print ("setup before each test method, disabling all features")
        self.shapeFeatures.disableAllFeatures()

    @classmethod
    def setup_class(self):
        print ("") # this is to get a newline after the dots
        print ("setup_class() before any methods in this class")

        # set the patient ID for these files to match the directory and
        # the patient id in the baseline file
        self.patientID = 'TestPatient1_BreastMRI'

        dataDir = os.path.dirname(os.path.abspath(__file__)) + os.path.sep + ".." + os.path.sep + "data" + os.path.sep

        imageName = dataDir + self.patientID + os.path.sep + 'TestPatient1_BreastMRI-Pre-subvolume.nrrd'
        maskName = dataDir + self.patientID + os.path.sep + 'TestPatient1_BreastMRI-Pre_labelMap-subvolume.nrrd'

        print("Reading the image and mask.")
        self.image = sitk.ReadImage(imageName)
        self.mask = sitk.ReadImage(maskName)

        print("Instantiating Shape.")
        self.shapeFeatures = shape.RadiomicsShape(self.image, self.mask)

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
      if key == 'Volume':
        index = 'Shape_volume'
      elif key == 'SurfaceArea':
        index = 'Shape_surface'
      elif key == 'SurfaceVolumeRatio':
          index = 'Shape_surfVolRatio'
      elif key == 'Compactness1':
        index = 'Shape_compactness'
      elif key == 'Compactness2':
        index = 'Shape_compactness2'
      elif key == 'Maximum3DDiameter':
        index = 'Shape_maxDiameter3D'
      elif key == 'SphericalDisproportion':
        index = 'Shape_spherDisprop'
      elif key == 'Sphericity':
        index = 'Shape_sphericity'
      if index == -1:
        print 'Unable to find index for key ',key
        return
      baseline = self.baselineFeatures[index]
      percentDiff = abs(1.0 - (value / float(baseline)))
      print('index = %s, baseline value = %f, calculated = %f, diff = %f%%' % (index, float(baseline), value, percentDiff * 100))
      # check for a less than one percent difference
      assert(percentDiff < 0.01)

    def test_volume_BreastMRI(self):
        testString = 'Volume'
        print 'Test', testString
        self.shapeFeatures.enableFeatureByName(testString)
        self.shapeFeatures.calculateFeatures()
        val = self.shapeFeatures.featureValues[testString]
        self.checkResult(testString, val)

    def test_surfaceArea_BreastMRI(self):
        testString = 'SurfaceArea'
        print 'Test', testString
        self.shapeFeatures.enableFeatureByName(testString)
        self.shapeFeatures.calculateFeatures()
        val = self.shapeFeatures.featureValues[testString]
        self.checkResult(testString, val)

    def test_surfaceVolumeRatio_BreastMRI(self):
        testString = 'SurfaceVolumeRatio'
        print 'Test', testString
        self.shapeFeatures.enableFeatureByName(testString)
        self.shapeFeatures.calculateFeatures()
        val = self.shapeFeatures.featureValues[testString]
        self.checkResult(testString, val)

    def test_compactness1_BreastMRI(self):
        testString = 'Compactness1'
        print 'Test', testString
        self.shapeFeatures.enableFeatureByName(testString)
        self.shapeFeatures.calculateFeatures()
        val = self.shapeFeatures.featureValues[testString]
        self.checkResult(testString, val)

    def test_compactness2_BreastMRI(self):
        testString = 'Compactness2'
        print 'Test', testString
        self.shapeFeatures.enableFeatureByName(testString)
        self.shapeFeatures.calculateFeatures()
        val = self.shapeFeatures.featureValues[testString]
        self.checkResult(testString, val)

    def test_maximum3DDiameter_BreastMRI(self):
        testString = 'Maximum3DDiameter'
        print 'Test', testString
        self.shapeFeatures.enableFeatureByName(testString)
        self.shapeFeatures.calculateFeatures()
        val = self.shapeFeatures.featureValues[testString]
        self.checkResult(testString, val)

    def test_sphericalDisproportion_BreastMRI(self):
        testString = 'SphericalDisproportion'
        print 'Test', testString
        self.shapeFeatures.enableFeatureByName(testString)
        self.shapeFeatures.calculateFeatures()
        val = self.shapeFeatures.featureValues[testString]
        self.checkResult(testString, val)

    def test_sphericity_BreastMRI(self):
        testString = 'Sphericity'
        print 'Test', testString
        self.shapeFeatures.enableFeatureByName(testString)
        self.shapeFeatures.calculateFeatures()
        val = self.shapeFeatures.featureValues[testString]
        self.checkResult(testString, val)
