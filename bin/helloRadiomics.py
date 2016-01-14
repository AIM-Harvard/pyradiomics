from radiomics import preprocessing, firstorder, shape, glcm, rlgl, laplacian
import SimpleITK as sitk
import sys, os
import pdb

#testBinWidth = 25 this is the default bin size
#testResampledPixelSpacing = (3,3,3) no resampling for now.

dataDir = os.path.dirname(os.path.abspath(__file__)) + os.path.sep + ".." + os.path.sep + "data"
#imageName = str(dataDir + os.path.sep + 'prostate_phantom_subvolume.nrrd')
#maskName = str(dataDir + os.path.sep + 'prostate_phantom_subvolume-label.nrrd')
imageName = os.path.join(dataDir, 'breast1_image.nrrd')
maskName = os.path.join(dataDir, 'breast1_label.nrrd')

if not os.path.exists(imageName):
  print 'Error: problem finding input image',imageName
  os.exit()
if not os.path.exists(maskName):
  print 'Error: problem finding input image',maskName
  os.exit()

image = sitk.ReadImage(imageName)
mask = sitk.ReadImage(maskName)

#
# Show the first order feature calculations
#
firstOrderFeatures = firstorder.RadiomicsFirstOrder(image,mask)

#firstOrderFeatures.enableFeatureByName('MeanIntensity', True)
firstOrderFeatures.enableAllFeatures()

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
glcmFeatures = glcm.RadiomicsGLCM(image, mask)
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
rlglFeatures = rlgl.RadiomicsRLGL(image, mask)
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
# Show Laplacian Of Gaussian features 
#

laplacianFeatures = laplacian.RadiomicsLaplacian(image, mask)
laplacianFeatures.enableAllFeatures()

print 'Will calculate the following Laplacian features: '
for f in laplacianFeatures.enabledFeatures.keys():
  print '  ',f
  print eval('laplacianFeatures.get'+f+'FeatureValue.__doc__')

print 'Calculating Laplacian features...',
laplacianFeatures.calculateFeatures()
print 'done'

print 'Calculated Laplacian features: '
for (key,val) in laplacianFeatures.featureValues.iteritems():
  print '  ',key,':',val  