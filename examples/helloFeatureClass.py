#!/usr/bin/env python

from __future__ import print_function

import os

import numpy
import SimpleITK as sitk
import six

from radiomics import firstorder, getTestCase, glcm, glrlm, glszm, imageoperations, shape

# testBinWidth = 25 this is the default bin size
# testResampledPixelSpacing = [3,3,3] no resampling for now.

# Get some test data

# repositoryRoot points to the root of the repository. The following line gets that location if this script is run
# from it's default location in /pyradiomics/examples. Otherwise, it will point to some (invalid) folder, causing the
# getTestCase function to fail to find the test case in the repository. In that case, a test case will be downloaded to
# temporary files and it's location is returned.
repositoryRoot = os.path.abspath(os.path.join(os.getcwd(), ".."))
imageName, maskName = getTestCase('brain1', repositoryRoot)

if imageName is None or maskName is None:  # Something went wrong, in this case PyRadiomics will also log an error
  print('Error getting testcase!')
  exit()

image = sitk.ReadImage(imageName)
mask = sitk.ReadImage(maskName)

applyLog = False
applyWavelet = False

# Setting for the feature calculation.
# Currently, resampling is disabled.
# Can be enabled by setting 'resampledPixelSpacing' to a list of 3 floats (new voxel size in mm for x, y and z)
settings = {'binWidth': 25,
            'interpolator': sitk.sitkBSpline,
            'resampledPixelSpacing': None}

#
# If enabled, resample image (resampled image is automatically cropped.
# If resampling is not enabled, crop image instead
#
if settings['interpolator'] is not None and settings['resampledPixelSpacing'] is not None:
  image, mask = imageoperations.resampleImage(image, mask, settings['resampledPixelSpacing'], settings['interpolator'])
else:
  bb, correctedMask = imageoperations.checkMask(image, mask)
  if correctedMask is not None:
    mask = correctedMask
  image, mask = imageoperations.cropToTumorMask(image, mask, bb)

#
# Show the first order feature calculations
#
firstOrderFeatures = firstorder.RadiomicsFirstOrder(image, mask, **settings)

firstOrderFeatures.enableFeatureByName('Mean', True)
# firstOrderFeatures.enableAllFeatures()

print('Will calculate the following first order features: ')
for f in firstOrderFeatures.enabledFeatures.keys():
  print('  ', f)
  print(getattr(firstOrderFeatures, 'get%sFeatureValue' % f).__doc__)

print('Calculating first order features...')
firstOrderFeatures.calculateFeatures()
print('done')

print('Calculated first order features: ')
for (key, val) in six.iteritems(firstOrderFeatures.featureValues):
  print('  ', key, ':', val)

#
# Show Shape features
#
shapeFeatures = shape.RadiomicsShape(image, mask, **settings)
shapeFeatures.enableAllFeatures()

print('Will calculate the following Shape features: ')
for f in shapeFeatures.enabledFeatures.keys():
  print('  ', f)
  print(getattr(shapeFeatures, 'get%sFeatureValue' % f).__doc__)

print('Calculating Shape features...')
shapeFeatures.calculateFeatures()
print('done')

print('Calculated Shape features: ')
for (key, val) in six.iteritems(shapeFeatures.featureValues):
  print('  ', key, ':', val)

#
# Show GLCM features
#
glcmFeatures = glcm.RadiomicsGLCM(image, mask, **settings)
glcmFeatures.enableAllFeatures()

print('Will calculate the following GLCM features: ')
for f in glcmFeatures.enabledFeatures.keys():
  print('  ', f)
  print(getattr(glcmFeatures, 'get%sFeatureValue' % f).__doc__)

print('Calculating GLCM features...')
glcmFeatures.calculateFeatures()
print('done')

print('Calculated GLCM features: ')
for (key, val) in six.iteritems(glcmFeatures.featureValues):
  print('  ', key, ':', val)

#
# Show GLRLM features
#
glrlmFeatures = glrlm.RadiomicsGLRLM(image, mask, **settings)
glrlmFeatures.enableAllFeatures()

print('Will calculate the following GLRLM features: ')
for f in glrlmFeatures.enabledFeatures.keys():
  print('  ', f)
  print(getattr(glrlmFeatures, 'get%sFeatureValue' % f).__doc__)

print('Calculating GLRLM features...')
glrlmFeatures.calculateFeatures()
print('done')

print('Calculated GLRLM features: ')
for (key, val) in six.iteritems(glrlmFeatures.featureValues):
  print('  ', key, ':', val)

#
# Show GLSZM features
#
glszmFeatures = glszm.RadiomicsGLSZM(image, mask, **settings)
glszmFeatures.enableAllFeatures()

print('Will calculate the following GLSZM features: ')
for f in glszmFeatures.enabledFeatures.keys():
  print('  ', f)
  print(getattr(glszmFeatures, 'get%sFeatureValue' % f).__doc__)

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
  for logImage, imageTypeName, inputKwargs in imageoperations.getLoGImage(image, sigma=sigmaValues):
    logFirstorderFeatures = firstorder.RadiomicsFirstOrder(logImage, mask, **inputKwargs)
    logFirstorderFeatures.enableAllFeatures()
    logFirstorderFeatures.calculateFeatures()
    for (key, val) in six.iteritems(logFirstorderFeatures.featureValues):
      laplacianFeatureName = '%s_%s' % (imageTypeName, key)
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
      waveletFeatureName = '%s_%s' % (str(decompositionName), key)
      print('  ', waveletFeatureName, ':', val)
