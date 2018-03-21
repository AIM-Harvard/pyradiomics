import pywt
import six

from radiomics import getFeatureClasses, getImageTypes

featureClasses = getFeatureClasses()
imageTypes = getImageTypes()

def checkWavelet(value, rule_obj, path):
  if not isinstance(value, six.string_types):
    raise TypeError('Wavelet not expected type (str)')
  wavelist = pywt.wavelist()
  if value not in wavelist:
    raise ValueError('Wavelet "%s" not available in pyWavelets %s' % (value, wavelist))
  return True


def checkInterpolator(value, rule_obj, path):
  if value is None:
    return True
  if isinstance(value, six.string_types):
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
  elif isinstance(value, six.string_types):
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
  for className, features in six.iteritems(value):
    if className not in featureClasses.keys():
      raise ValueError(
        'Feature Class %s is not recognized. Available feature classes are %s' % (className, list(featureClasses.keys())))
    if features is not None:
      if not isinstance(features, list):
        raise TypeError('Value of feature class %s not expected type (list)' % (className))
      unrecognizedFeatures = set(features) - set(featureClasses[className].getFeatureNames())
      if len(unrecognizedFeatures) > 0:
        raise ValueError('Feature Class %s contains unrecognized features: %s' % (className, str(unrecognizedFeatures)))

  return True


def checkImageType(value, rule_obj, path):
  global imageTypes
  if value is None:
    raise TypeError('imageType dictionary cannot be None value')

  for im_type in value:
    if im_type not in imageTypes:
      raise ValueError('Image Type %s is not recognized. Available image types are %s' %
                       (im_type, imageTypes))

  return True
