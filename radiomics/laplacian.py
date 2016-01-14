import numpy
import inspect, importlib
from radiomics import base, preprocessing, firstorder, glcm, rlgl
import SimpleITK as sitk
import pdb

def add_dynamo(parent, sigma, featureClass, featureName, featureFunction):
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
  setattr(parent, getFeatureValue_dynamo.__name__, getFeatureValue_dynamo)
      
class RadiomicsLaplacian(base.RadiomicsFeaturesBase):

  def __init__(self, inputImage, inputMask, **kwargs):
    super(RadiomicsLaplacian,self).__init__(inputImage,inputMask,**kwargs)
    
    self.sigmaValues = numpy.arange(5.0, 0.0, -0.5)[::-1]
    self.featureClasses = ['firstorder', 'glcm', 'rlgl']
    
    lrgif = sitk.LaplacianRecursiveGaussianImageFilter()
    for sigma in self.sigmaValues:
      lrgif.SetSigma(sigma)
      inputImage_laplacian = lrgif.Execute(self.inputImage)
      
      for featureClass in self.featureClasses:
        featureModule = importlib.import_module("radiomics.%s" % featureClass)
        featureClass_RadiomicsObj_ = inspect.getmembers(featureModule, \
                                      lambda member: (inspect.isclass(member)) \
                                      and ('Radiomics' in member.__name__))[0][1]  
        featureClass_instance = featureClass_RadiomicsObj_(inputImage_laplacian, self.inputMask, **kwargs)
        featureClass_instance.enableAllFeatures()
        
        methods = [method for method in inspect.getmembers(featureClass_instance, inspect.ismethod) if method[0].endswith('FeatureValue')]
        for method in methods:
          featureName = method[0][3:-12]
          featureFunction = method[1]
          add_dynamo(self, sigma, featureClass, featureName, featureFunction)
          
    
    super(RadiomicsLaplacian,self).__init__(inputImage,inputMask,**kwargs)
    #self.InitializeFeatureVector()
    #for f in self.getFeatureNames():
    #  self.enabledFeatures[f] = True

    # TODO: add an option to instantiate the class that reuses initialization