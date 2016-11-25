import pywt


def checkWavelet(value, rule_obj, path):
  if not isinstance(value, basestring):
    raise TypeError('Wavelet not expected type (str)')
  wavelist = pywt.wavelist()
  if value not in wavelist:
    raise ValueError('Wavelet "%s" not available in pyWavelets %s' % (value, wavelist))
  return True


def checkInterpolator(value, rule_obj, path):
  if value is None:
    return True
  if isinstance(value, basestring):
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
  elif isinstance(value, basestring):
    enum = ['euclidean', 'manhattan', 'infinity', 'no_weighting']
    if value not in enum:
      raise ValueError('WeightingNorm value "%s" not valid, possible values: %s' % (value, enum))
  else:
    raise TypeError('WeightingNorm not expected type (str or None)')
  return True
