from radiomics import firstorder, glcm, imageoperations, shape, rlgl, glszm, wavelet
import SimpleITK as sitk
import sys, os

#testBinWidth = 25 this is the default bin size
#testResampledPixelSpacing = (3,3,3) no resampling for now.

dataDir = os.path.dirname(os.path.abspath(__file__)) + os.path.sep + ".." + os.path.sep + "data"
#imageName = str(dataDir + os.path.sep + 'prostate_phantom_subvolume.nrrd')
#maskName = str(dataDir + os.path.sep + 'prostate_phantom_subvolume-label.nrrd')
imageName = str(dataDir + os.path.sep + 'breast1_image.nrrd')
maskName = str(dataDir + os.path.sep + 'breast1_label.nrrd')

if not os.path.exists(imageName):
  print 'Error: problem finding input image',imageName
  os.exit()
if not os.path.exists(maskName):
  print 'Error: problem finding input image',maskName
  os.exit()

image = sitk.ReadImage(imageName)
mask = sitk.ReadImage(maskName)
"""
#
# Show the first order feature calculations
#
firstOrderFeatures = firstorder.RadiomicsFirstOrder(image,mask)

firstOrderFeatures.enableFeatureByName('Mean', True)
# firstOrderFeatures.enableAllFeatures()

print 'Will calculate the following first order features: '
for f in firstOrderFeatures.enabledFeatures.keys():
  print '  ',f
  print eval('firstOrderFeatures.get'+f+'FeatureValue.__doc__')

print 'Calculating first order features...',
firstOrderFeatures.calculateFeatures()
print 'done'

print 'Calculated first order features: '
for (key,val) in firstOrderFeatures.featureValues.iteritems():
  print '  ',key,':',val

  # now calculate the LoG features, following
  # steps and parameters as specified by @vnarayan13 in #22
  # 1. Get the LoG filtered image
  mmif = sitk.MinimumMaximumImageFilter()
  mmif.Execute(image)
  lowerThreshold = 0
  upperThreshold = mmif.GetMaximum()

  threshImage = imageoperations.applyThreshold(image,lowerThreshold=lowerThreshold, upperThreshold=upperThreshold,outsideValue=0)
  # get the mask of the thresholded pixels
  threshImageMask = imageoperations.applyThreshold(image,lowerThreshold=lowerThreshold, upperThreshold=upperThreshold,outsideValue=0,insideValue=1)
  # only include the voxels that are within the threshold
  threshMask = sitk.Cast(mask,1) & sitk.Cast(threshImageMask,1)
  import numpy
  sigmaValues = numpy.arange(5.,0.,-.5)[::1]
  for sigma in sigmaValues:
    logImage = imageoperations.applyLoG(image,sigmaValue=sigma)
    logFirstorderFeatures = firstorder.RadiomicsFirstOrder(logImage,threshMask)
    logFirstorderFeatures.enableAllFeatures()
    logFirstorderFeatures.calculateFeatures()
    print 'Calculated firstorder features with LoG sigma ',sigma
    for (key,val) in logFirstorderFeatures.featureValues.iteritems():
      laplacianFeatureName = 'LoG_sigma_%s_%s' %(str(sigma),key)
      print '  ',laplacianFeatureName,':',val

#
# Show Shape features
#
shapeFeatures = shape.RadiomicsShape(image, mask)
shapeFeatures.enableAllFeatures()

print 'Will calculate the following Shape features: '
for f in shapeFeatures.enabledFeatures.keys():
  print '  ',f
  print eval('shapeFeatures.get'+f+'FeatureValue.__doc__')

print 'Calculating Shape features...',
shapeFeatures.calculateFeatures()
print 'done'

print 'Calculated Shape features: '
for (key,val) in shapeFeatures.featureValues.iteritems():
  print '  ',key,':',val


#
# Show GLCM features
#
glcmFeatures = glcm.RadiomicsGLCM(image, mask, binWidth=25)
glcmFeatures.enableAllFeatures()

print 'Will calculate the following GLCM features: '
for f in glcmFeatures.enabledFeatures.keys():
  print '  ',f
  print eval('glcmFeatures.get'+f+'FeatureValue.__doc__')

print 'Calculating GLCM features...',
glcmFeatures.calculateFeatures()
print 'done'

print 'Calculated GLCM features: '
for (key,val) in glcmFeatures.featureValues.iteritems():
  print '  ',key,':',val

#
# Show RLGL features
#
rlglFeatures = rlgl.RadiomicsRLGL(image, mask, binWidth=25)
rlglFeatures.enableAllFeatures()

print 'Will calculate the following RLGL features: '
for f in rlglFeatures.enabledFeatures.keys():
  print '  ',f
  print eval('rlglFeatures.get'+f+'FeatureValue.__doc__')

print 'Calculating RLGL features...',
rlglFeatures.calculateFeatures()
print 'done'

print 'Calculated RLGL features: '
for (key,val) in rlglFeatures.featureValues.iteritems():
  print '  ',key,':',val

#
# Show GLSZM features
#
glszmFeatures = glszm.RadiomicsGLSZM(image, mask, binWidth=25)
glszmFeatures.enableAllFeatures()

print 'Will calculate the following GLSZM features: '
for f in glszmFeatures.enabledFeatures.keys():
  print '  ',f
  print eval('glszmFeatures.get'+f+'FeatureValue.__doc__')

print 'Calculating GLSZM features...',
glszmFeatures.calculateFeatures()
print 'done'

print 'Calculated GLSZM features: '
for (key,val) in glszmFeatures.featureValues.iteritems():
  print '  ',key,':',val
"""

#
# Show Wavelet features
#
parameters = {}
parameters['binWidth'] = 25
parameters['resampledPixelSpacing'] = None
parameters['interpolator'] = sitk.sitkBSpline
parameters['padDistance'] = 5
parameters['padFillValue'] = 0
#parameters['wavelet_waveletType'] = 'coif1'

waveletFeatures = wavelet.RadiomicsWavelet(image, mask, **parameters)
waveletFeatures.enableAllFeatures()

print 'Will calculate the following Wavelet features: '
for f in waveletFeatures.enabledFeatures.keys():
  print '  ', f
  print eval('waveletFeatures.get'+f+'FeatureValue.__doc__')
  
print 'Calculating Wavelet features...',
waveletFeatures.calculateFeatures()
print 'done'

print 'Calculated Wavelet features: '
for (key,val) in waveletFeatures.featureValues.iteritems():
  print '  ',key,':',val