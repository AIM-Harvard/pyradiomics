import pywt
from radiomics.featureextractor import RadiomicsFeaturesExtractor
from radiomics import c_str_type

featureClasses = RadiomicsFeaturesExtractor.getFeatureClasses()


def checkWavelet(value, rule_obj, path):
  if not isinstance(value, c_str_type):
    raise TypeError('Wavelet not expected type (str)')
  wavelist = pywt.wavelist()
  if value not in wavelist:
    raise ValueError('Wavelet "%s" not available in pyWavelets %s' % (value, wavelist))
  return True


def checkInterpolator(value, rule_obj, path):
  if value is None:
    return True
  if isinstance(value, c_str_type):
    enum = {'sitkNearestNeighbor',
            'sitkLinear',
            'sitkBSpline',
            'sitkGaussian',
            'sitkLabelGaussian',
            'sitkHammingWindowedSinc',
            'sitkCosineWindowedSinc',
            'sitkWelchWindowedSinc',
            'sitkLanczosWindowedSinc',
            'sitkBlackmanWindowedSinc'}
    if value not in enum:
      raise ValueError('Interpolator value "%s" not valid, possible values: %s' % (value, enum))
  elif isinstance(value, int):
    if value < 1 or value > 10:
      raise ValueError('Intepolator value %i, must be in range of [1-10]' % (value))
  else:
    raise TypeError('Interpolator not expected type (str or int)')
  return True


def checkWeighting(value, rule_obj, path):
  if value is None:
    return True
  elif isinstance(value, c_str_type):
    enum = ['euclidean', 'manhattan', 'infinity', 'no_weighting']
    if value not in enum:
      raise ValueError('WeightingNorm value "%s" not valid, possible values: %s' % (value, enum))
  else:
    raise TypeError('WeightingNorm not expected type (str or None)')
  return True


def checkFeatureClass(value, rule_obj, path):
  global featureClasses
  if value is None:
    raise TypeError('featureClass dictionary cannot be None value')
  for className, features in value.items():
    if className not in featureClasses.keys():
      raise ValueError(
        'Feature Class %s is not recognized. Available feature classes are %s' % (className, featureClasses.keys()))
    if features is not None:
      if not isinstance(features, list):
        raise TypeError('Value of feature class %s not expected type (list)' % (className))
      unrecognizedFeatures = set(features) - set(featureClasses[className].getFeatureNames())
      if len(unrecognizedFeatures) > 0:
        raise ValueError('Feature Class %s contains unrecognized features: %s' % (className, str(unrecognizedFeatures)))

  return True
