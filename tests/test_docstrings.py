import logging

import six

from radiomics import getFeatureClasses

featureClasses = getFeatureClasses()


def pytest_generate_tests(metafunc):
  metafunc.parametrize(["featureClassName", "featureName"], metafunc.cls.generate_scenarios())


class TestDocStrings:

  @staticmethod
  def generate_scenarios():
    global featureClasses
    for featureClassName, featureClass in six.iteritems(featureClasses):
      logging.info('generate_scenarios %s', featureClassName)
      doc = featureClass.__doc__
      assert (doc is not None)

      featureNames = featureClass.getFeatureNames()
      for f in featureNames:
        yield (featureClassName, f)

  def test_class(self, featureClassName, featureName):
    global featureClasses
    logging.info('%s', featureName)
    features = featureClasses[featureClassName]
    doc = getattr(features, "get%sFeatureValue" % featureName).__doc__
    logging.info('%s', doc)
    assert (doc is not None)
