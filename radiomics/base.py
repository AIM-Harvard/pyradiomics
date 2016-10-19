import SimpleITK as sitk
import numpy
import inspect
from radiomics import imageoperations

class RadiomicsFeaturesBase(object):
  def __init__(self, inputImage, inputMask, **kwargs):
    """
    Initialization

    inputImage and inputMask are SimpleITK images. The motivation for using
    SimpleITK images as input is to keep the possibility of reusing the
    optimized feature calculators implemented in SimpleITK in the future.
    """

    self.binWidth = kwargs.get('binWidth', 25)
    self.verbose = kwargs.get('verbose', True)

    # all features are disabled by default
    self.disableAllFeatures()

    self.featureNames = self.getFeatureNames()

    self.inputImage = inputImage
    self.inputMask = inputMask

    if inputImage == None or inputMask == None:
      if self.verbose: print('ERROR: missing input image or mask')
      return

    self.imageArray = sitk.GetArrayFromImage(self.inputImage)
    self.maskArray = sitk.GetArrayFromImage(self.inputMask)

    self.matrix = self.imageArray.astype('int64')
    self.matrixCoordinates = numpy.where(self.maskArray != 0)

    self.targetVoxelArray = self.matrix[self.matrixCoordinates]

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
      try:
        self.featureValues[feature] = eval(call)
      except Exception:
        self.featureValues[feature] = numpy.nan
