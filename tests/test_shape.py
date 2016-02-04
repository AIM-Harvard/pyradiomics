# to run this test, from directory above:
# setenv PYTHONPATH /path/to/pyradiomics/radiomics
# nosetests --nocapture -v tests/test_shape.py

import SimpleITK as sitk
from radiomics import firstorder, shape
from testUtils import RadiomicsTestUtils
import sys, os

def setup_module(module):
    # run before anything in this file
    print ("") # this is to get a newline after the dots
    return

class TestShape:

    def setup(self):
        # setup before each test method
        print ("") # this is to get a newline after the dots
        self.shapeFeatures.disableAllFeatures()

    @classmethod
    def setup_class(self):
        # run before any methods in this class
        print ("") # this is to get a newline after the dots

        # set the patient ID for these files to match the directory and
        # the patient id in the baseline file
        self.patientID = 'breast1'

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
        # run after any methods in this class
        print ("") # this is to get a newline after the dots

    def test_volume_BreastMRI(self):
        testString = 'Volume'
        if self.testUtils.getVerbose():
          print 'Test', testString
        self.shapeFeatures.enableFeatureByName(testString)
        self.shapeFeatures.calculateFeatures()
        val = self.shapeFeatures.featureValues[testString]
        self.testUtils.checkResult(testString, val)

    def test_surfaceArea_BreastMRI(self):
        testString = 'SurfaceArea'
        if self.testUtils.getVerbose():
          print 'Test', testString
        self.shapeFeatures.enableFeatureByName(testString)
        self.shapeFeatures.calculateFeatures()
        val = self.shapeFeatures.featureValues[testString]
        self.testUtils.checkResult(testString, val)

    def test_surfaceVolumeRatio_BreastMRI(self):
        testString = 'SurfaceVolumeRatio'
        if self.testUtils.getVerbose():
          print 'Test', testString
        self.shapeFeatures.enableFeatureByName(testString)
        self.shapeFeatures.calculateFeatures()
        val = self.shapeFeatures.featureValues[testString]
        self.testUtils.checkResult(testString, val)

    def test_compactness1_BreastMRI(self):
        testString = 'Compactness1'
        if self.testUtils.getVerbose():
          print 'Test', testString
        self.shapeFeatures.enableFeatureByName(testString)
        self.shapeFeatures.calculateFeatures()
        val = self.shapeFeatures.featureValues[testString]
        self.testUtils.checkResult(testString, val)

    def test_compactness2_BreastMRI(self):
        testString = 'Compactness2'
        if self.testUtils.getVerbose():
          print 'Test', testString
        self.shapeFeatures.enableFeatureByName(testString)
        self.shapeFeatures.calculateFeatures()
        val = self.shapeFeatures.featureValues[testString]
        self.testUtils.checkResult(testString, val)

    def test_maximum3DDiameter_BreastMRI(self):
        testString = 'Maximum3DDiameter'
        if self.testUtils.getVerbose():
          print 'Test', testString
        self.shapeFeatures.enableFeatureByName(testString)
        self.shapeFeatures.calculateFeatures()
        val = self.shapeFeatures.featureValues[testString]
        self.testUtils.checkResult(testString, val)

    def test_sphericalDisproportion_BreastMRI(self):
        testString = 'SphericalDisproportion'
        if self.testUtils.getVerbose():
          print 'Test', testString
        self.shapeFeatures.enableFeatureByName(testString)
        self.shapeFeatures.calculateFeatures()
        val = self.shapeFeatures.featureValues[testString]
        self.testUtils.checkResult(testString, val)

    def test_sphericity_BreastMRI(self):
        testString = 'Sphericity'
        if self.testUtils.getVerbose():
          print 'Test', testString
        self.shapeFeatures.enableFeatureByName(testString)
        self.shapeFeatures.calculateFeatures()
        val = self.shapeFeatures.featureValues[testString]
        self.testUtils.checkResult(testString, val)
