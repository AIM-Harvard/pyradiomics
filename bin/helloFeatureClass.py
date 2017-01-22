import os
import SimpleITK as sitk
import numpy
from radiomics import firstorder, glcm, imageoperations, shape, glrlm, glszm

# testBinWidth = 25 this is the default bin size
# testResampledPixelSpacing = [3,3,3] no resampling for now.

dataDir = os.path.abspath("") + os.path.sep + ".." + os.path.sep + "data"  # assumes this file in in pyradiomics/bin
# imageName = str(dataDir + os.path.sep + 'prostate_phantom_subvolume.nrrd')
# maskName = str(dataDir + os.path.sep + 'prostate_phantom_subvolume-label.nrrd')
imageName = str(dataDir + os.path.sep + 'breast1_image.nrrd')
maskName = str(dataDir + os.path.sep + 'breast1_label.nrrd')

if not os.path.exists(imageName):
  print 'Error: problem finding input image', imageName
  exit()
if not os.path.exists(maskName):
  print 'Error: problem finding input image', maskName
  exit()

image = sitk.ReadImage(imageName)
mask = sitk.ReadImage(maskName)

applyLog = False
applyWavelet = False

# Setting for the feature calculation.
# Currently, resampling is disabled.
# Can be enabled by setting 'resampledPixelSpacing' to a list of 3 floats (new voxel size in mm for x, y and z)
kwargs = {'binWidth': 25,
          'interpolator': sitk.sitkBSpline,
          'resampledPixelSpacing': None,
          'verbose': True}

#
# If enabled, resample image (resampled image is automatically cropped.
# If resampling is not enabled, crop image instead
#
if kwargs['interpolator'] != None and kwargs['resampledPixelSpacing'] != None:
  image, mask = imageoperations.resampleImage(image, mask, kwargs['resampledPixelSpacing'], kwargs['interpolator'])
else:
  image, mask = imageoperations.cropToTumorMask(image, mask)

#
# Show the first order feature calculations
#
firstOrderFeatures = firstorder.RadiomicsFirstOrder(image, mask, **kwargs)

firstOrderFeatures.enableFeatureByName('Mean', True)
# firstOrderFeatures.enableAllFeatures()

print 'Will calculate the following first order features: '
for f in firstOrderFeatures.enabledFeatures.keys():
  print '  ', f
  print eval('firstOrderFeatures.get' + f + 'FeatureValue.__doc__')

print 'Calculating first order features...',
firstOrderFeatures.calculateFeatures()
print 'done'

print 'Calculated first order features: '
for (key, val) in firstOrderFeatures.featureValues.iteritems():
  print '  ', key, ':', val

#
# Show Shape features
#
shapeFeatures = shape.RadiomicsShape(image, mask, **kwargs)
shapeFeatures.enableAllFeatures()

print 'Will calculate the following Shape features: '
for f in shapeFeatures.enabledFeatures.keys():
  print '  ', f
  print eval('shapeFeatures.get' + f + 'FeatureValue.__doc__')

print 'Calculating Shape features...',
shapeFeatures.calculateFeatures()
print 'done'

print 'Calculated Shape features: '
for (key, val) in shapeFeatures.featureValues.iteritems():
  print '  ', key, ':', val

#
# Show GLCM features
#
glcmFeatures = glcm.RadiomicsGLCM(image, mask, **kwargs)
glcmFeatures.enableAllFeatures()

print 'Will calculate the following GLCM features: '
for f in glcmFeatures.enabledFeatures.keys():
  print '  ', f
  print eval('glcmFeatures.get' + f + 'FeatureValue.__doc__')

print 'Calculating GLCM features...',
glcmFeatures.calculateFeatures()
print 'done'

print 'Calculated GLCM features: '
for (key, val) in glcmFeatures.featureValues.iteritems():
  print '  ', key, ':', val

#
# Show GLRLM features
#
glrlmFeatures = glrlm.RadiomicsGLRLM(image, mask, **kwargs)
glrlmFeatures.enableAllFeatures()

print 'Will calculate the following GLRLM features: '
for f in glrlmFeatures.enabledFeatures.keys():
  print '  ', f
  print eval('glrlmFeatures.get' + f + 'FeatureValue.__doc__')

print 'Calculating GLRLM features...',
glrlmFeatures.calculateFeatures()
print 'done'

print 'Calculated GLRLM features: '
for (key, val) in glrlmFeatures.featureValues.iteritems():
  print '  ', key, ':', val

#
# Show GLSZM features
#
glszmFeatures = glszm.RadiomicsGLSZM(image, mask, **kwargs)
glszmFeatures.enableAllFeatures()

print 'Will calculate the following GLSZM features: '
for f in glszmFeatures.enabledFeatures.keys():
  print '  ', f
  print eval('glszmFeatures.get' + f + 'FeatureValue.__doc__')

print 'Calculating GLSZM features...',
glszmFeatures.calculateFeatures()
print 'done'

print 'Calculated GLSZM features: '
for (key, val) in glszmFeatures.featureValues.iteritems():
  print '  ', key, ':', val

#
# Show FirstOrder features, calculated on a LoG filtered image
#
if applyLog:
  sigmaValues = numpy.arange(5., 0., -.5)[::1]
  for sigma in sigmaValues:
    logImage = imageoperations.applyLoG(image, sigmaValue=sigma)
    logFirstorderFeatures = firstorder.RadiomicsFirstOrder(logImage, mask, **kwargs)
    logFirstorderFeatures.enableAllFeatures()
    logFirstorderFeatures.calculateFeatures()
    print 'Calculated firstorder features with LoG sigma ', sigma
    for (key, val) in logFirstorderFeatures.featureValues.iteritems():
      laplacianFeatureName = 'LoG-sigma-%s_%s' % (str(sigma), key)
      print '  ', laplacianFeatureName, ':', val
#
# Show FirstOrder features, calculated on a wavelet filtered image
#
if applyWavelet:
  ret, approx = imageoperations.swt3(image)

  for idx, wl in enumerate(ret, start=1):
    for decompositionName, decompositionImage in wl.items():
      waveletFirstOrderFeaturs = firstorder.RadiomicsFirstOrder(decompositionImage, mask, **kwargs)
      waveletFirstOrderFeaturs.enableAllFeatures()
      waveletFirstOrderFeaturs.calculateFeatures()
      print 'Calculated firstorder features with wavelet ', decompositionName
      for (key, val) in waveletFirstOrderFeaturs.featureValues.iteritems():
        waveletFeatureName = 'wavelet-%s_%s' % (str(decompositionName), key)
        print '  ', waveletFeatureName, ':', val

  waveletFirstOrderFeaturs = firstorder.RadiomicsFirstOrder(approx, mask, **kwargs)
  waveletFirstOrderFeaturs.enableAllFeatures()
  waveletFirstOrderFeaturs.calculateFeatures()
  print 'Calculated firstorder features with approximation of wavelt (= LLL decomposition)'
  for (key, val) in waveletFirstOrderFeaturs.featureValues.iteritems():
    waveletFeatureName = 'wavelet-LLL_%s' % (key)
    print '  ', waveletFeatureName, ':', val
