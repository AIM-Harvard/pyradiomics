# to run this test, from directory above:
# setenv PYTHONPATH /path/to/pyradiomics/radiomics
# nosetests --nocapture -v tests/firstorderTest.py

from radiomics import firstorder, imageoperations
from testUtils import RadiomicsTestUtils
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
        print ("") # this is to get a newline after the dots
        print ("teardown_class() after any methods in this class")

    def checkResult(self, key, value):
      # use the mapping from the utils
      baseline = self.testUtils.getMatlabValue(key)

      percentDiff = abs(1.0 - (value / float(baseline)))
      print('baseline value = %f, calculated = %f, diff = %f%%' % (float(baseline), value, percentDiff * 100))
      # check for a less than one percent difference
      assert(percentDiff < 0.01)

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
        testString = 'Minimum'
        print 'Test', testString
        self.firstOrderFeatures.enableFeatureByName(testString)
        self.firstOrderFeatures.calculateFeatures()
        val = self.firstOrderFeatures.featureValues[testString]
        self.checkResult(testString, val)

    def test_maxIntensity_10(self):
        testString = 'Maximum'
        print 'Test', testString
        self.firstOrderFeatures.enableFeatureByName(testString)
        self.firstOrderFeatures.calculateFeatures()
        val = self.firstOrderFeatures.featureValues[testString]
        self.checkResult(testString, val)

    def test_meanIntensity_10(self):
        testString = 'Mean'
        print 'Test', testString
        self.firstOrderFeatures.enableFeatureByName(testString)
        self.firstOrderFeatures.calculateFeatures()
        val = self.firstOrderFeatures.featureValues[testString]
        self.checkResult(testString, val)

    def test_medianIntensity_10(self):
        testString = 'Median'
        print 'Test', testString
        self.firstOrderFeatures.enableFeatureByName(testString)
        self.firstOrderFeatures.calculateFeatures()
        val = self.firstOrderFeatures.featureValues[testString]
        self.checkResult(testString, val)

    def test_rangeIntensity_10(self):
        testString = 'Range'
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
        testString = 'Skewness'
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
