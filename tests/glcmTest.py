# to run this test, from directory above:
# setenv PYTHONPATH /Users/nicole/Radiomics/github/pyradiomics/radiomics:/Users/nicole/Slicer4-svn/Slicer-build-debug/SimpleITK-build/Wrapping
# setenv DYLD_LIBRARY_PATH /Users/nicole/Slicer4-svn/Slicer-build-debug/python-install/lib:/Users/nicole/Slicer4-svn/Slicer-build-debug/SimpleITK-build/lib
# nosetests --nocapture -v tests/glcmTest.py

from radiomics import firstorder, glcm, preprocessing
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
        self.patientID = 'TCGA-02-0003_BrainMRI'

        dataDir = os.path.dirname(os.path.abspath(__file__)) + os.path.sep + ".." + os.path.sep + "data" + os.path.sep

        imageName = dataDir + self.patientID + os.path.sep + 'AXIAL-T1-POST-GD_TCGA-02-0003_TCGA-02-0003_PREOP_SCAN13_Background.nrrd'
        maskName = dataDir + self.patientID + os.path.sep + 'AXIAL-T1-POST-GD_TCGA-02-0003_TCGA-02-0003_PREOP_SCAN13_FSL-Clustered-labelMap_binarized_labelMap.nrrd'

        print("Reading the image and mask.")
        self.image = sitk.ReadImage(imageName)
        self.mask = sitk.ReadImage(maskName)

        print("Instantiating GLCM.")
        self.glcmFeatures = glcm.RadiomicsGLCM(self.image, self.mask, 10)

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
      if key == 'Autocorrelation':
        index = 'GLCM_autocorr'
      elif key == 'ClusterProminence':
        index = 'GLCM_clusProm'
      elif key == 'ClusterShade':
          index = 'GLCM_clusShade'
      elif key == 'ClusterTendency':
        index = 'GLCM_clusTend'
      elif key == 'Contrast':
        index = 'GLCM_contrast'
      elif key == 'Correlation':
        index = 'GLCM_correl1'
      elif key == 'DifferenceEntropy':
        index = 'GLCM_diffEntro'
      elif key == 'Dissimilarity':
        index = 'GLCM_dissimilar'
      elif key == 'Energy':
        index = 'GLCM_energy'
      elif key == 'Entropy':
        index = 'GLCM_entrop2'
      elif key == 'Homogeneity1':
        index = 'GLCM_homogeneity1'
      elif key == 'Homogeneity2':
        index = 'GLCM_homogeneity2'
      elif key == 'Imc1':
        index = 'GLCM_infoCorr1'
      elif key == 'Idmn':
        index = 'GLCM_invDiffmomnor'
      elif key == 'Idn':
        index = 'GLCM_invDiffnorm'
      elif key == 'InverseVariance':
        index = 'GLCM_inverseVar'
      elif key == 'MaximumProbability':
        index = 'GLCM_maxProb'
      elif key == 'SumAverage':
        index = 'GLCM_sumAvg'
      elif key == 'SumEntropy':
        index = 'GLCM_sumEntro'
      elif key == 'SumVariance':
        index = 'GLCM_sumVar'
      elif key == 'SumSquares':
        index = 'GLCM_sumSquares'
      if index == -1:
        print 'Unable to find index for key ',key
        return
      baseline = self.baselineFeatures[index]
      diff = abs(float(baseline) - value)
      print 'index = ', index, ', baseline value = ', baseline, ', calculated = ', value, ', diff = ', diff
      assert(diff < 0.1)

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
