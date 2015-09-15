class RadiomicsFeaturesBase:
  def __init__(self):
    self.enabledFeatures = {}

  def setInputImageArray(image):
    self.image = image

  def setInputMaskArray(mask):
    self.mask = mask

  def setInputCoordinatesArray(coord):
    self.coord = coord

  def setPadDistance(padDistance):
    self.padDistance = padDistance

  def setPadFillValue(padFillValue):
    self.padFillValue = padFillValue

  def enableFeatureByName(featureName, enable):
    if not featureName in self.featureNames:
      raise LookupError('Feature not found: '+featureName)
    self.enabledFeatures[featureName] = True

  def getFeatureNames(self):
    allMembers = dir(self)
    allFeatureNames = [f[3:-12] for f in allMembers if f.endswith('FeatureValue') and f.startswith('get')]
    return allFeatureNames
