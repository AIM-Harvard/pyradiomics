# to run this test, from directory above:
# setenv PYTHONPATH /path/to/pyradiomics/radiomics
# nosetests --nocapture -v tests/glcmTest.py

from radiomics import firstorder, glcm, imageoperations
from testUtils import RadiomicsTestUtils
import SimpleITK as sitk
import sys, os
import csv

def setup_module(module):
    print ("") # this is to get a newline after the dots
    print ("setup_module before anything in this file")
    return

class TestGLCM:

    def setup(self):
        print ("") # this is to get a newline after the dots
        # print ("setup before each test method, disabling all features")
        self.glcmFeatures.disableAllFeatures()

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

        print("Instantiating GLCM.")
        self.glcmFeatures = glcm.RadiomicsGLCM(self.image, self.mask)

    @classmethod
    def teardown_class(self):
        print ("") # this is to get a newline after the dots
        print ("teardown_class() after any methods in this class")

    def checkResult(self, key, value):
      # use the mapping from the utils
      baseline = self.testUtils.getMatlabValue(key)
      print('checkResults: for key %s, got baseline = %f' % (key, baseline))
      if baseline == 0.0:
        # avoid divide by zero, the difference is either 0% if the value is also zero, or 100%
        if value - baseline == 0.0:
          percentDiff = 0.0
        else:
          percentDiff = 1.0
      else:
        percentDiff = abs(1.0 - (value / baseline))
      print('baseline value = %f, calculated = %f, diff = %f%%' % (baseline, value, percentDiff * 100))
      # check for a less than one percent difference
      assert(percentDiff < 0.01)

    def test_autocorrelation_10(self):
        self.glcmFeatures.enableFeatureByName('Autocorrelation')
        print 'Will calculate the following GLCM features: '
        for f in self.glcmFeatures.enabledFeatures.keys():
            print '  ',f
            print eval('self.glcmFeatures.get'+f+'FeatureValue.__doc__')

        self.glcmFeatures.calculateFeatures()

        print 'Calculated GLCM features: '
        for (key,val) in self.glcmFeatures.featureValues.iteritems():
            print '  ',key,':',val
            if key == 'Autocorrelation':
                autoCorr = val

        self.checkResult('Autocorrelation', autoCorr)

    def test_clusterProminence_10(self):
        self.glcmFeatures.enableFeatureByName('ClusterProminence')
        self.glcmFeatures.calculateFeatures()
        print 'Calculated GLCM features: '
        for (key,val) in self.glcmFeatures.featureValues.iteritems():
            print '  ',key,':',val
            if key == 'ClusterProminence':
                clusterProminence = val

        self.checkResult('ClusterProminence', clusterProminence)

    def test_clusterShade_10(self):
        testString = 'ClusterShade'
        print 'Test', testString
        self.glcmFeatures.enableFeatureByName(testString)
        self.glcmFeatures.calculateFeatures()
        val = self.glcmFeatures.featureValues[testString]
        self.checkResult(testString, val)

    def test_clusterTendency_10(self):
        testString = 'ClusterTendency'
        print 'Test', testString
        self.glcmFeatures.enableFeatureByName(testString)
        self.glcmFeatures.calculateFeatures()
        val = self.glcmFeatures.featureValues[testString]
        self.checkResult(testString, val)

    def test_contrast_10(self):
        testString = 'Contrast'
        print 'Test', testString
        self.glcmFeatures.enableFeatureByName(testString)
        self.glcmFeatures.calculateFeatures()
        val = self.glcmFeatures.featureValues[testString]
        self.checkResult(testString, val)

    def test_correlation_10(self):
        testString = 'Correlation'
        print 'Test', testString
        self.glcmFeatures.enableFeatureByName(testString)
        self.glcmFeatures.calculateFeatures()
        val = self.glcmFeatures.featureValues[testString]
        self.checkResult(testString, val)

    def test_differenceEntropy_10(self):
        testString = 'DifferenceEntropy'
        print 'Test', testString
        self.glcmFeatures.enableFeatureByName(testString)
        self.glcmFeatures.calculateFeatures()
        val = self.glcmFeatures.featureValues[testString]
        self.checkResult(testString, val)

    def test_dissimilarity_10(self):
        testString = 'Dissimilarity'
        print 'Test', testString
        self.glcmFeatures.enableFeatureByName(testString)
        self.glcmFeatures.calculateFeatures()
        val = self.glcmFeatures.featureValues[testString]
        self.checkResult(testString, val)

    def test_energy_10(self):
        testString = 'Energy'
        print 'Test', testString
        self.glcmFeatures.enableFeatureByName(testString)
        self.glcmFeatures.calculateFeatures()
        val = self.glcmFeatures.featureValues[testString]
        self.checkResult(testString, val)

    def test_entropy_10(self):
        testString = 'Entropy'
        print 'Test', testString
        self.glcmFeatures.enableFeatureByName(testString)
        self.glcmFeatures.calculateFeatures()
        val = self.glcmFeatures.featureValues[testString]
        self.checkResult(testString, val)

    def test_homogeneity1_10(self):
        testString = 'Homogeneity1'
        print 'Test', testString
        self.glcmFeatures.enableFeatureByName(testString)
        self.glcmFeatures.calculateFeatures()
        val = self.glcmFeatures.featureValues[testString]
        self.checkResult(testString, val)

    def test_homogeneity2_10(self):
        testString = 'Homogeneity2'
        print 'Test', testString
        self.glcmFeatures.enableFeatureByName(testString)
        self.glcmFeatures.calculateFeatures()
        val = self.glcmFeatures.featureValues[testString]
        self.checkResult(testString, val)

    def test_imc1_10(self):
        testString = 'Imc1'
        print 'Test', testString
        self.glcmFeatures.enableFeatureByName(testString)
        self.glcmFeatures.calculateFeatures()
        val = self.glcmFeatures.featureValues[testString]
        self.checkResult(testString, val)

    def test_idmn_10(self):
        testString = 'Idmn'
        print 'Test', testString
        self.glcmFeatures.enableFeatureByName(testString)
        self.glcmFeatures.calculateFeatures()
        val = self.glcmFeatures.featureValues[testString]
        self.checkResult(testString, val)

    def test_idn_10(self):
        testString = 'Idn'
        print 'Test', testString
        self.glcmFeatures.enableFeatureByName(testString)
        self.glcmFeatures.calculateFeatures()
        val = self.glcmFeatures.featureValues[testString]
        self.checkResult(testString, val)

    def test_inverseVariance_10(self):
        testString = 'InverseVariance'
        print 'Test', testString
        self.glcmFeatures.enableFeatureByName(testString)
        self.glcmFeatures.calculateFeatures()
        val = self.glcmFeatures.featureValues[testString]
        self.checkResult(testString, val)

    def test_maximumProbability_10(self):
        testString = 'MaximumProbability'
        print 'Test', testString
        self.glcmFeatures.enableFeatureByName(testString)
        self.glcmFeatures.calculateFeatures()
        val = self.glcmFeatures.featureValues[testString]
        self.checkResult(testString, val)

    def test_sumAverage_10(self):
        testString = 'SumAverage'
        print 'Test', testString
        self.glcmFeatures.enableFeatureByName(testString)
        self.glcmFeatures.calculateFeatures()
        val = self.glcmFeatures.featureValues[testString]
        self.checkResult(testString, val)

    def test_sumVariance_10(self):
        testString = 'SumVariance'
        print 'Test', testString
        self.glcmFeatures.enableFeatureByName(testString)
        self.glcmFeatures.calculateFeatures()
        val = self.glcmFeatures.featureValues[testString]
        self.checkResult(testString, val)

    def test_sumSquares_10(self):
        testString = 'SumSquares'
        print 'Test', testString
        self.glcmFeatures.enableFeatureByName(testString)
        self.glcmFeatures.calculateFeatures()
        val = self.glcmFeatures.featureValues[testString]
        self.checkResult(testString, val)
