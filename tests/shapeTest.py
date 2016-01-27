# to run this test, from directory above:
# setenv PYTHONPATH /path/to/pyradiomics/radiomics
# nosetests --nocapture -v tests/shapeTest.py

import SimpleITK as sitk
from radiomics import firstorder, shape
from testUtils import RadiomicsTestUtils
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
        self.patientID = 'brain1'

        # read in the baseline and mapping to matlab features
        self.testUtils = RadiomicsTestUtils('shape')
        self.testUtils.setPatientID(self.patientID)

        dataDir = os.path.dirname(os.path.abspath(__file__)) + os.path.sep + ".." + os.path.sep + "data" + os.path.sep

        imageName = dataDir + self.patientID + '_image.nrrd'
        maskName = dataDir + self.patientID + '_label.nrrd'

        print("Reading the image and mask.")
        self.image = sitk.ReadImage(imageName)
        self.mask = sitk.ReadImage(maskName)

        print("Instantiating Shape.")
        self.shapeFeatures = shape.RadiomicsShape(self.image, self.mask)


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
