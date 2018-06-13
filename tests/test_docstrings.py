# to run this test, from directory above:
# setenv PYTHONPATH /path/to/pyradiomics/radiomics
# nosetests --nocapture -v tests/test_docstrings.py

import logging

from nose_parameterized import parameterized
import six

from radiomics import getFeatureClasses
from testUtils import custom_name_func

featureClasses = getFeatureClasses()


def setup_module(module):
  # runs before anything in this file
  print("")  # this is to get a newline after the dots
  return


class TestDocStrings:
  def setup(self):
    # setup before each test method
    print("")  # this is to get a newline after the dots

  @classmethod
  def setup_class(cls):
    # called before any methods in this class
    print("")  # this is to get a newline after the dots

  @classmethod
  def teardown_class(cls):
    # run after any methods in this class
    print("")  # this is to get a newline after the dots

  def generate_scenarios():
    global featureClasses
    for featureClassName, featureClass in six.iteritems(featureClasses):
      logging.info('generate_scenarios %s', featureClassName)
      doc = featureClass.__doc__
      assert (doc is not None)

      featureNames = featureClass.getFeatureNames()
      for f in featureNames:
        yield (featureClassName, f)

  @parameterized.expand(generate_scenarios(), testcase_func_name=custom_name_func)
  def test_class(self, featureClassName, featureName):
    global featureClasses
    logging.info('%s', featureName)
    features = featureClasses[featureClassName]
    doc = getattr(features, "get%sFeatureValue" % featureName).__doc__
    logging.info('%s', doc)
    assert (doc is not None)
