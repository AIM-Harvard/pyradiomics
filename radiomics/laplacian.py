import numpy
import inspect, importlib
from radiomics import base, preprocessing, firstorder, glcm, rlgl
import SimpleITK as sitk
import pdb

#TODO:
#normalize sigmaX, sigmaY, sigmaZ by *(1/pixelspacing). #is this done in ITK?
#crop image using margin from maximum filter size
#examine how filter rescaling is done in ITK:
"""
Filter rescaling
  The sum of the filter values should always be zero to surpress
  homogeneous regions in the original image. However the response,
  both positive and negative, can be scaled to give filter responses
  in the same range as the original image. This can be achieved by
  multiplying the filter with a scale factor.
"""
#implement Positive-voxel only functions for firstorder (Pos suffix)
#BUG: Entropy calculations failing on LoG images.     
#also implement 2D LoG Filter?
  
class RadiomicsLaplacian(base.RadiomicsFeaturesBase):

  def __init__(self, inputImage, inputMask, **kwargs):
    super(RadiomicsLaplacian,self).__init__(inputImage,inputMask,**kwargs)
    
    self.sigmaValues = numpy.arange(5.0, 0.0, -0.5)[::-1]
    
    #Exclude "Air" (only for CT):
    #Ganeshan paper(s):  exclude -50 HU or lower
    #Kemerink 2008 paper: exclude -400 HU or lower
    #Conservative Threshold: exclude -1000 HU or lower
    #applyHUThreshold set to -1000
    self.applyHUThreshold = -1000
    #set temp to self.padValue after achieving self.padValue=NaN
    self.thresholdedImage, self.thresholdedMask = self.thresholdImageAndMask(temp=-2000) 
    
    self.featureClasses = ['firstorder', 'glcm', 'rlgl']
    
    lrgif = sitk.LaplacianRecursiveGaussianImageFilter()
    lrgif.SetNormalizeAcrossScale(True)
    for sigma in self.sigmaValues:
      lrgif.SetSigma(sigma)
      inputImage_laplacian = lrgif.Execute(self.thresholdedImage)
      
      for featureClass in self.featureClasses:
        featureModule = importlib.import_module("radiomics.%s" % featureClass)
        featureClass_RadiomicsObj_ = inspect.getmembers(featureModule, \
                                      lambda member: (inspect.isclass(member)) \
                                      and ('Radiomics' in member.__name__))[0][1]  
        featureClass_instance = featureClass_RadiomicsObj_(inputImage_laplacian, self.thresholdedMask, **kwargs)
        featureClass_instance.enableAllFeatures()
        
        methods = [method for method in inspect.getmembers(featureClass_instance, inspect.ismethod) if method[0].endswith('FeatureValue')]
        for method in methods:
          featureName = method[0][3:-12]
          featureFunction = method[1]
          self.add_dynamo(sigma, featureClass, featureName, featureFunction)
          
    
    super(RadiomicsLaplacian,self).__init__(inputImage,inputMask,**kwargs)
    #self.InitializeFeatureVector()
    #for f in self.getFeatureNames():
    #  self.enabledFeatures[f] = True

    # TODO: add an option to instantiate the class that reuses initialization

  def thresholdImageAndMask(self, temp=0):
    mmif = sitk.MinimumMaximumImageFilter()
    mmif.Execute(self.inputImage)
    maximumIntensity = mmif.GetMaximum()
    
    tif = sitk.ThresholdImageFilter()
    tif.SetLower(self.applyHUThreshold)
    tif.SetUpper(maximumIntensity)
    tif.SetOutsideValue(temp)
    thresholdedImage = tif.Execute(self.inputImage)
    
    # set label value to 0 where image pixels were thresholded
    thresholdedImageArray = sitk.GetArrayFromImage(thresholdedImage)
    thresholdedVoxels = zip(*numpy.where(thresholdedImageArray==temp))
    
    thresholdedMask = self.inputMask
    for iz, ix, iy in thresholdedVoxels:
      thresholdedMask.SetPixel(ix,iy,iz,0)
    
    return thresholdedImage, thresholdedMask
    
  def add_dynamo(self, sigma, featureClass, featureName, featureFunction):
    """Dynamically Generates Function Definition during run-time"""
    
    laplacianFeatureName = "laplacian_sigma_%s_mm_3D_%s_%s" %(str(sigma).replace('.','_'),featureClass,featureName)   
    def getFeatureValue_dynamo():   
      try:
        featureValue = featureFunction()
      except:
        featureValue = "Failed"
      finally:  
        return featureValue 
       
    getFeatureValue_dynamo.__doc__ = \
      """Applies a Recursive Laplacian of a Gaussian Filter with sigma=%s on the """ \
      """input image and computes %s.%s\n%s""" %(sigma,featureClass,featureName,featureFunction.__doc__)  
    getFeatureValue_dynamo.__name__ = "get%sFeatureValue" % laplacianFeatureName
    setattr(self, getFeatureValue_dynamo.__name__, getFeatureValue_dynamo)