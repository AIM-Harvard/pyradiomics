import radiomics
import pandas as pd
import itertools

class TestWaveletsOnROIs:

  def setup(self):
    self.testData = []
    self.testData.append({"image": "../data/test-image_64x64x64.nrrd", 
                          "mask": "../data/test-mask_64x64x64.nrrd"})
    self.testData.append({"image": "../data/test-image_37x37x37.nrrd", 
                          "mask": "../data/test-mask_37x37x37.nrrd"})


  def test_features_equal_on_various_rois_of_same_image(self):
    extractors = []
    settings = {}
    settings['binWidth'] = 10
    extractors.append(self.getTestExtrator(settings))
    settings['preCrop'] = True
    settings['padDistance'] = 6
    extractors.append(self.getTestExtrator(settings))

    features = []
    for extractor in extractors:
      for data in self.testData:
        features.append(extractor.execute(data["image"], data["mask"], label = 1))
    
    dfFeatures = pd.DataFrame(features)
    for column in dfFeatures:
      if not column.startswith("general"):
        featList = dfFeatures[column].tolist()
        maxDiff = max(abs(a) - abs(b) for a, b in itertools.combinations(featList, 2))
        assert maxDiff == 0


  def getTestExtrator(self, settings):
    extractor = radiomics.featureextractor.RadiomicsFeaturesExtractor(**settings)
    extractor.disableAllImageTypes()
    extractor.disableAllFeatures()
    extractor.enableImageTypes(Wavelet={})
    extractor.enableFeaturesByName(firstorder=['Median','Mean','Minimum','Maximum'])
    extractor.enableFeaturesByName(glcm=['Contrast','Correlation','JointEnergy','JointEntropy','Idm'])
    return extractor

  


