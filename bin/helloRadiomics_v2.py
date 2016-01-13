from radiomics import firstorder, glcm, preprocessing, shape, rlgl, wavelet

import SimpleITK as sitk
import sys, os
import inspect
import pdb

def getImageAndMask():
  dataDir = os.path.dirname(os.path.abspath(__file__)) + os.path.sep + ".." + os.path.sep + "data"

  #imageName = str(os.path.join(dataDir, 'prostate_phantom_subvolume.nrrd'))
  #maskName = str(os.path.join(dataDir, 'prostate_phantom_subvolume-label.nrrd'))
  imageName = str(os.path.join(dataDir, 'breast1_image.nrrd'))
  maskName = str(os.path.join(dataDir, 'breast1_label.nrrd'))

  if not os.path.exists(imageName):
    print 'Error: problem finding input image',imageName
    os.exit()
  if not os.path.exists(maskName):
    print 'Error: problem finding input image',maskName
    os.exit()

  image = sitk.ReadImage(imageName)
  mask = sitk.ReadImage(maskName)
  return image, mask
  
def getRadiomicsParameters():
    parameters = {}
    parameters['binWidth'] = 25
    parameters['resampledPixelSpacing'] = None
    parameters['interpolator'] = sitk.sitkBSpline
    parameters['padDistance'] = 5
    parameters['padFillValue'] = 0 # TODO: update features so default can be NaN
    #parameters['laplacian_sigmaStart'] = 0
    #parameters['laplacian_sigmaEnd'] = 5
    #parameters['laplacian_sigmaStep'] = 0.5
    #parameters['laplacian_applyHUThreshold'] = -1000
    #parameters['wavelet_waveletType'] = 'coif1'   
    return parameters
  
def printDocstrings_inspect(featureClassName, featureClass_instance):
  # Not Used
  print 'Will calculate the following %s features: ' % featureClassName
  featureMethods = [method for method in inspect.getmembers(featureClass_instance, inspect.ismethod) \
              if (method[0].startswith('get')) \
              and (method[0].endswith('FeatureValue'))]
  
  enabledFeatures = featureClass_instance.enabledFeatures.keys()
  for featureName, featureFunction in featureMethods:
    featureName = featureName[3:-12]
    if featureName in enabledFeatures:
      print "\t%s" % featureName
      print featureFunction.__doc__

def printDocstrings(featureClassName, featureClass_instance):
  print 'Will calculate the following %s features: ' % featureClassName
  featureMethods = [method for method in dir(featureClass_instance) \
            if (method.startswith('get')) \
            and (method.endswith('FeatureValue'))]
            
  for feature in featureClass_instance.enabledFeatures.keys():
    try:
      featureFunctionName = [f for f in featureMethods if feature in f][0]
    except:
      pdb.set_trace()
    print '\t', feature
    print eval('featureClass_instance.%s.__doc__' % featureFunctionName)
      
def calculateFeatures(featureClassName, featureClass_instance):
  print 'Calculating %s features...' % featureClassName
  featureClass_instance.calculateFeatures()
  print 'done'

def printFeatureValues(featureClassName, featureClass_instance):
  print 'Calculated %s features: ' % featureClassName
  for (key,val) in featureClass_instance.featureValues.iteritems():
    print '\t',key,':',val  

    
def main():
  image, mask = getImageAndMask()
  kwargs = getRadiomicsParameters()
  """
  #
  # Show firstorder features
  #
  firstOrderFeatures = firstorder.RadiomicsFirstOrder(image,mask,**kwargs)
  firstOrderFeatures.enableAllFeatures()
  printDocstrings('firstorder', firstOrderFeatures)
  calculateFeatures('firstorder', firstOrderFeatures)
  printFeatureValues('firstorder', firstOrderFeatures)

  #
  # Show shape features
  #
  shapeFeatures = shape.RadiomicsShape(image, mask, **kwargs)
  shapeFeatures.enableAllFeatures()
  printDocstrings('shape', shapeFeatures)
  calculateFeatures('shape', shapeFeatures)
  printFeatureValues('shape', shapeFeatures)

  #
  # Show glcm features
  #
  glcmFeatures = glcm.RadiomicsGLCM(image, mask, **kwargs)
  glcmFeatures.enableAllFeatures()
  printDocstrings('glcm', glcmFeatures)
  calculateFeatures('glcm', glcmFeatures)
  printFeatureValues('glcm', glcmFeatures)

  #
  # Show rlgl features
  #
  rlglFeatures = rlgl.RadiomicsRLGL(image, mask, **kwargs)
  rlglFeatures.enableAllFeatures()
  printDocstrings('rlgl', rlglFeatures)
  calculateFeatures('rlgl', rlglFeatures)
  printFeatureValues('rlgl', rlglFeatures)
  """
  
  #
  # Show Wavelet features
  #
  waveletFeatures = wavelet.RadiomicsWavelet(image, mask, **kwargs)
  waveletFeatures.enableAllFeatures()
  printDocstrings('wavelet', waveletFeatures)
  calculateFeatures('wavelet', waveletFeatures)
  printFeatureValues('wavelet', waveletFeatures)

if __name__=="__main__":
  main()