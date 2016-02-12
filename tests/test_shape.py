# to run this test, from directory above:
# setenv PYTHONPATH /path/to/pyradiomics/radiomics
# nosetests --nocapture -v tests/test_shape.py

import SimpleITK as sitk
from radiomics import firstorder, shape
from testUtils import RadiomicsTestUtils
import sys, os
import logging

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

        # read in the baseline and mapping to matlab features
        self.testUtils = RadiomicsTestUtils('shape')
        # set the test case for these files to match the directory and
        # the id in the baseline file
        self.testUtils.setTestCase('breast1')

        dataDir = os.path.dirname(os.path.abspath(__file__)) + os.path.sep + ".." + os.path.sep + "data" + os.path.sep

        imageName = dataDir + self.testUtils.getTestCase() + '_image.nrrd'
        maskName = dataDir + self.testUtils.getTestCase() + '_label.nrrd'

        logging.info("Reading the image and mask.")
        self.image = sitk.ReadImage(imageName)
        self.mask = sitk.ReadImage(maskName)

        logging.info("Instantiating Shape.")
        self.shapeFeatures = shape.RadiomicsShape(self.image, self.mask)


    @classmethod
    def teardown_class(self):
        # run after any methods in this class
        print ("") # this is to get a newline after the dots

    def test_volume_BreastMRI(self):
        testString = 'Volume'
        logging.info('Test %s', testString)
        self.shapeFeatures.enableFeatureByName(testString)
        self.shapeFeatures.calculateFeatures()
        val = self.shapeFeatures.featureValues[testString]
        self.testUtils.checkResult(testString, val)

    def test_surfaceArea_BreastMRI(self):
        testString = 'SurfaceArea'
        logging.info('Test %s', testString)
        self.shapeFeatures.enableFeatureByName(testString)
        self.shapeFeatures.calculateFeatures()
        val = self.shapeFeatures.featureValues[testString]
        self.testUtils.checkResult(testString, val)

    def test_surfaceVolumeRatio_BreastMRI(self):
        testString = 'SurfaceVolumeRatio'
        logging.info('Test %s', testString)
        self.shapeFeatures.enableFeatureByName(testString)
        self.shapeFeatures.calculateFeatures()
        val = self.shapeFeatures.featureValues[testString]
        self.testUtils.checkResult(testString, val)

    def test_compactness1_BreastMRI(self):
        testString = 'Compactness1'
        logging.info('Test %s', testString)
        self.shapeFeatures.enableFeatureByName(testString)
        self.shapeFeatures.calculateFeatures()
        val = self.shapeFeatures.featureValues[testString]
        self.testUtils.checkResult(testString, val)

    def test_compactness2_BreastMRI(self):
        testString = 'Compactness2'
        logging.info('Test %s', testString)
        self.shapeFeatures.enableFeatureByName(testString)
        self.shapeFeatures.calculateFeatures()
        val = self.shapeFeatures.featureValues[testString]
        self.testUtils.checkResult(testString, val)

    def test_maximum3DDiameter_BreastMRI(self):
        testString = 'Maximum3DDiameter'
        logging.info('Test %s', testString)
        self.shapeFeatures.enableFeatureByName(testString)
        self.shapeFeatures.calculateFeatures()
        val = self.shapeFeatures.featureValues[testString]
        self.testUtils.checkResult(testString, val)

    def test_sphericalDisproportion_BreastMRI(self):
        testString = 'SphericalDisproportion'
        logging.info('Test %s', testString)
        self.shapeFeatures.enableFeatureByName(testString)
        self.shapeFeatures.calculateFeatures()
        val = self.shapeFeatures.featureValues[testString]
        self.testUtils.checkResult(testString, val)

    def test_sphericity_BreastMRI(self):
        testString = 'Sphericity'
        logging.info('Test %s', testString)
        self.shapeFeatures.enableFeatureByName(testString)
        self.shapeFeatures.calculateFeatures()
        val = self.shapeFeatures.featureValues[testString]
        self.testUtils.checkResult(testString, val)
