
from __future__ import print_function

import os

import numpy
import SimpleITK as sitk
import six

from radiomics import firstorder, glcm, glrlm, glszm, imageoperations, shape

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
if kwargs['interpolator'] is not None and kwargs['resampledPixelSpacing'] is not None:
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
for (key, val) in six.iteritems(firstOrderFeatures.featureValues):
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
for (key, val) in six.iteritems(shapeFeatures.featureValues):
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
for (key, val) in six.iteritems(glcmFeatures.featureValues):
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
for (key, val) in six.iteritems(glrlmFeatures.featureValues):
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
for (key, val) in six.iteritems(glszmFeatures.featureValues):
  print('  ', key, ':', val)

#
# Show FirstOrder features, calculated on a LoG filtered image
#
if applyLog:
  sigmaValues = numpy.arange(5., 0., -.5)[::1]
  for logImage, inputImageName, inputKwargs in imageoperations.getLoGImage(image, sigma=sigmaValues, verbose=True):
    logFirstorderFeatures = firstorder.RadiomicsFirstOrder(logImage, mask, **inputKwargs)
    logFirstorderFeatures.enableAllFeatures()
    logFirstorderFeatures.calculateFeatures()
    for (key, val) in six.iteritems(logFirstorderFeatures.featureValues):
      laplacianFeatureName = '%s_%s' % (inputImageName, key)
      print('  ', laplacianFeatureName, ':', val)
#
# Show FirstOrder features, calculated on a wavelet filtered image
#
if applyWavelet:
  for decompositionImage, decompositionName, inputKwargs in imageoperations.getWaveletImage(image):
    waveletFirstOrderFeaturs = firstorder.RadiomicsFirstOrder(decompositionImage, mask, **inputKwargs)
    waveletFirstOrderFeaturs.enableAllFeatures()
    waveletFirstOrderFeaturs.calculateFeatures()
    print('Calculated firstorder features with wavelet ', decompositionName)
    for (key, val) in six.iteritems(waveletFirstOrderFeaturs.featureValues):
      waveletFeatureName = 'wavelet-%s_%s' % (str(decompositionName), key)
      print('  ', waveletFeatureName, ':', val)
