from __future__ import print_function, unicode_literals, division, absolute_import
import os
import SimpleITK as sitk
import numpy
from radiomics import firstorder, glcm, imageoperations, shape, glrlm, glszm

# testBinWidth = 25 this is the default bin size
# testResampledPixelSpacing = [3,3,3] no resampling for now.

testCase = 'brain1'
dataDir = os.path.join(os.path.abspath(""), "..", "data")
imageName = os.path.join(dataDir, testCase + '_image.nrrd')
maskName = os.path.join(dataDir, testCase + '_label.nrrd')

if not os.path.exists(imageName):
  print('Error: problem finding input image', imageName)
  exit()
if not os.path.exists(maskName):
  print('Error: problem finding input image', maskName)
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
  image, mask, bb = imageoperations.cropToTumorMask(image, mask)

#
# Show the first order feature calculations
#
firstOrderFeatures = firstorder.RadiomicsFirstOrder(image, mask, **kwargs)

firstOrderFeatures.enableFeatureByName('Mean', True)
# firstOrderFeatures.enableAllFeatures()

print('Will calculate the following first order features: ')
for f in firstOrderFeatures.enabledFeatures.keys():
  print('  ', f)
  print(eval('firstOrderFeatures.get' + f + 'FeatureValue.__doc__'))

print('Calculating first order features...')
firstOrderFeatures.calculateFeatures()
print('done')

print('Calculated first order features: ')
for (key, val) in firstOrderFeatures.featureValues.items():
  print('  ', key, ':', val)

#
# Show Shape features
#
shapeFeatures = shape.RadiomicsShape(image, mask, **kwargs)
shapeFeatures.enableAllFeatures()

print('Will calculate the following Shape features: ')
for f in shapeFeatures.enabledFeatures.keys():
  print('  ', f)
  print(eval('shapeFeatures.get' + f + 'FeatureValue.__doc__'))

print('Calculating Shape features...')
shapeFeatures.calculateFeatures()
print('done')

print('Calculated Shape features: ')
for (key, val) in shapeFeatures.featureValues.items():
  print('  ', key, ':', val)

#
# Show GLCM features
#
glcmFeatures = glcm.RadiomicsGLCM(image, mask, **kwargs)
glcmFeatures.enableAllFeatures()

print('Will calculate the following GLCM features: ')
for f in glcmFeatures.enabledFeatures.keys():
  print('  ', f)
  print(eval('glcmFeatures.get' + f + 'FeatureValue.__doc__'))

print('Calculating GLCM features...')
glcmFeatures.calculateFeatures()
print('done')

print('Calculated GLCM features: ')
for (key, val) in glcmFeatures.featureValues.items():
  print('  ', key, ':', val)

#
# Show GLRLM features
#
glrlmFeatures = glrlm.RadiomicsGLRLM(image, mask, **kwargs)
glrlmFeatures.enableAllFeatures()

print('Will calculate the following GLRLM features: ')
for f in glrlmFeatures.enabledFeatures.keys():
  print('  ', f)
  print(eval('glrlmFeatures.get' + f + 'FeatureValue.__doc__'))

print('Calculating GLRLM features...')
glrlmFeatures.calculateFeatures()
print('done')

print('Calculated GLRLM features: ')
for (key, val) in glrlmFeatures.featureValues.items():
  print('  ', key, ':', val)

#
# Show GLSZM features
#
glszmFeatures = glszm.RadiomicsGLSZM(image, mask, **kwargs)
glszmFeatures.enableAllFeatures()

print('Will calculate the following GLSZM features: ')
for f in glszmFeatures.enabledFeatures.keys():
  print('  ', f)
  print(eval('glszmFeatures.get' + f + 'FeatureValue.__doc__'))

print('Calculating GLSZM features...')
glszmFeatures.calculateFeatures()
print('done')

print('Calculated GLSZM features: ')
for (key, val) in glszmFeatures.featureValues.items():
  print('  ', key, ':', val)

#
# Show FirstOrder features, calculated on a LoG filtered image
#
if applyLog:
  sigmaValues = numpy.arange(5., 0., -.5)[::1]
  for logImage, inputImageName, inputKwargs in imageoperations.applyFilterLoG(image, sigma=sigmaValues, verbose=True):
    logFirstorderFeatures = firstorder.RadiomicsFirstOrder(logImage, mask, **inputKwargs)
    logFirstorderFeatures.enableAllFeatures()
    logFirstorderFeatures.calculateFeatures()
    for (key, val) in logFirstorderFeatures.featureValues.items():
      laplacianFeatureName = '%s_%s' % (inputImageName, key)
      print('  ', laplacianFeatureName, ':', val)
#
# Show FirstOrder features, calculated on a wavelet filtered image
#
if applyWavelet:
  for decompositionImage, decompositionName, inputKwargs in imageoperations.applyFilterWavelet(image):
    waveletFirstOrderFeaturs = firstorder.RadiomicsFirstOrder(decompositionImage, mask, **inputKwargs)
    waveletFirstOrderFeaturs.enableAllFeatures()
    waveletFirstOrderFeaturs.calculateFeatures()
    print('Calculated firstorder features with wavelet ', decompositionName)
    for (key, val) in waveletFirstOrderFeaturs.featureValues.items():
      waveletFeatureName = 'wavelet-%s_%s' % (str(decompositionName), key)
      print('  ', waveletFeatureName, ':', val)
