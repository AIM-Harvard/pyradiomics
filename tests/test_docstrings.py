# to run this test, from directory above:
# setenv PYTHONPATH /path/to/pyradiomics/radiomics
# nosetests --nocapture -v tests/test_docstrings.py

from radiomics import firstorder, glcm, rlgl, shape
from testUtils import RadiomicsTestUtils
import SimpleITK as sitk
import sys, os

def setup_module(module):
    # runs before anything in this file
    print ("") # this is to get a newline after the dots
    return

class TestDocStrings:

    def setup(self):
        # setup before each test method
        print ("") # this is to get a newline after the dots

    @classmethod
    def setup_class(self):
        # called before any methods in this class
        print ("") # this is to get a newline after the dots

        # set the patient ID for these files to match the directory and
        # the patient id in the baseline file
        self.patientID = 'brain1'

        # read in the baseline and mapping to matlab features
        self.testUtils = RadiomicsTestUtils('testing')
        self.testUtils.setPatientID(self.patientID)

        dataDir = os.path.dirname(os.path.abspath(__file__)) + os.path.sep + ".." + os.path.sep + "data" + os.path.sep

        imageName = dataDir + self.patientID + '_image.nrrd'
        maskName = dataDir + self.patientID + '_label.nrrd'

        print("Reading the image and mask.")
        self.image = sitk.ReadImage(imageName)
        self.mask = sitk.ReadImage(maskName)



    @classmethod
    def teardown_class(self):
        # run after any methods in this class
        print ("") # this is to get a newline after the dots


    def test_firstOrder(self):
       print("Instantiating first order features.")
       self.features = firstorder.RadiomicsFirstOrder(self.image, self.mask)
       self.featureNames = firstorder.RadiomicsFirstOrder.getFeatureNames()
       for f in self.featureNames:
           print '  ',f
           doc = eval('self.features.get'+f+'FeatureValue.__doc__')
           print doc
           assert(doc != None)


    def test_glcm(self):
       print("Instantiating GLCM features.")
       self.features = glcm.RadiomicsGLCM(self.image, self.mask)
       self.featureNames = glcm.RadiomicsGLCM.getFeatureNames()
       for f in self.featureNames:
           print '  ',f
           doc = eval('self.features.get'+f+'FeatureValue.__doc__')
           print doc
           assert(doc != None)

    def test_rlgl(self):
       print("Instantiating RLGL features.")
       self.features = rlgl.RadiomicsRLGL(self.image, self.mask)
       self.featureNames = rlgl.RadiomicsRLGL.getFeatureNames()
       for f in self.featureNames:
           print '  ',f
           doc = eval('self.features.get'+f+'FeatureValue.__doc__')
           print doc
           assert(doc != None)

    def test_shape(self):
       print("Instantiating Shape features.")
       self.features = shape.RadiomicsShape(self.image, self.mask)
       self.featureNames = shape.RadiomicsShape.getFeatureNames()
       for f in self.featureNames:
           print '  ',f
           doc = eval('self.features.get'+f+'FeatureValue.__doc__')
           print doc
           assert(doc != None)
