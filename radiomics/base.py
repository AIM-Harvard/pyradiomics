import SimpleITK as sitk
import inspect
from radiomics import imageoperations

class RadiomicsFeaturesBase(object):
  def __init__(self, inputImage, inputMask, **kwargs):
    '''
    Initialization

    inputImage and inputMask are SimpleITK images. The motivation for using
    SimpleITK images as input is to keep the possibility of reusing the
    optimized feature calculators implemented in SimpleITK in the future.
    '''

    self.binWidth = 25
    self.resampledPixelSpacing = None # no resampling
    self.interpolator = sitk.sitkBSpline
    self.padDistance = 5 # no padding
    self.padFillValue = 0
    self.verbose = True

    for key,value in kwargs.iteritems():
      if key == 'binWidth':
        self.binWidth = value
      elif key == 'resampledPixelSpacing':
        self.resampledPixelSpacing = value
      elif key == 'interpolator':
        self.interpolator = value
      elif key == 'padDistance':
        self.padDistance = value
      elif key == 'padFillValue':
        self.padFillValue = value
      elif key == 'verbose':
        self.verbose = value
      elif self.verbose:
        print 'Warning: unknown parameter:',key

    # all features are disabled by default
    self.disableAllFeatures()

    self.featureNames = self.getFeatureNames()

    if self.interpolator != None and self.resampledPixelSpacing != None:
      self.inputImage, self.inputMask = imageoperations.interpolateImage(inputImage, inputMask, self.resampledPixelSpacing, self.interpolator)
    else:
      self.inputImage = inputImage
      self.inputMask = inputMask

  def enableFeatureByName(self, featureName, enable=True):
    if not featureName in self.featureNames:
      raise LookupError('Feature not found: '+featureName)
    self.enabledFeatures[featureName] = enable

  def enableAllFeatures(self):
    for featureName in self.featureNames:
      self.enableFeatureByName(featureName, True)

  def disableAllFeatures(self):
    self.enabledFeatures = {}
    self.featureValues = {}

  
  def getFeatureNames(self):
    allMembers = dir(self)
    allFeatureNames = [f[3:-12] for f in allMembers if f.endswith('FeatureValue') and f.startswith('get')]
    return allFeatureNames
  
  @classmethod
  def getFeatureNames(c):
    attributes = inspect.getmembers(c)
    features = [a[0][3:-12] for a in attributes if a[0].startswith('get') and a[0].endswith('FeatureValue')]
    return features
  
  def calculateFeatures(self):
    for feature in self.enabledFeatures.keys():
      call = 'self.get'+feature+'FeatureValue()'
      self.featureValues[feature] = eval(call)
