# to run this test, from directory above:
# setenv PYTHONPATH /path/to/pyradiomics/radiomics
# nosetests --nocapture -v tests/test_docstrings.py

from radiomics import firstorder, glcm, glrlm, shape, glszm
from testUtils import RadiomicsTestUtils
import SimpleITK as sitk
import sys, os
import logging
from nose_parameterized import parameterized

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


    @classmethod
    def teardown_class(self):
        # run after any methods in this class
        print ("") # this is to get a newline after the dots

    def generate_scenarios(featureClass):
      logging.info('generate_scenarios %s', featureClass)
      featureNames = featureClass.getFeatureNames()
      for f in featureNames:
        yield (f)

    @parameterized.expand(generate_scenarios(firstorder.RadiomicsFirstOrder))
    def test_firstOrder(self, featureName):
      logging.info('%s', featureName)
      features = firstorder.RadiomicsFirstOrder
      doc = eval('features.get'+featureName+'FeatureValue.__doc__')
      logging.info('%s', doc)
      assert(doc != None)


    @parameterized.expand(generate_scenarios(glcm.RadiomicsGLCM))
    def test_glcm(self, featureName):
      logging.info('%s', featureName)
      features = glcm.RadiomicsGLCM
      doc = eval('features.get'+featureName+'FeatureValue.__doc__')
      logging.info('%s', doc)
      assert(doc != None)


    @parameterized.expand(generate_scenarios(glrlm.RadiomicsGLRLM))
    def test_glrlm(self, featureName):
      logging.info('%s', featureName)
      features = glrlm.RadiomicsGLRLM
      doc = eval('features.get'+featureName+'FeatureValue.__doc__')
      logging.info('%s', doc)
      assert(doc != None)

    @parameterized.expand(generate_scenarios(shape.RadiomicsShape))
    def test_shape(self, featureName):
       logging.info('%s', featureName)
       features = shape.RadiomicsShape
       doc = eval('features.get'+featureName+'FeatureValue.__doc__')
       logging.info('%s', doc)
       assert(doc != None)

    @parameterized.expand(generate_scenarios(glszm.RadiomicsGLSZM))
    def test_glszm(self, featureName):
       logging.info('%s', featureName)
       features = glszm.RadiomicsGLSZM
       doc = eval('features.get'+featureName+'FeatureValue.__doc__')
       logging.info('%s', doc)
       assert(doc != None)
