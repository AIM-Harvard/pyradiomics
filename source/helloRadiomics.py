from radiomics import firstorder, glcm, preprocessing
import SimpleITK as sitk
import sys, os

#imageName = sys.argv[1]
#maskName = sys.argv[2]

dataDir = os.path.dirname(os.path.abspath(__file__)) + os.path.sep + ".." + os.path.sep + "data"
imageName = dataDir + os.path.sep + 'prostate_phantom-subvolume.nrrd'
maskName = dataDir + os.path.sep + 'prostate_phantom_label-subvolume.nrrd'

image = sitk.ReadImage(imageName)
mask = sitk.ReadImage(maskName)

#
# Show the first order feature calculations
#
firstOrderFeatures = firstorder.RadiomicsFirstOrder(image,mask)
firstOrderFeatures.setBinWidth(10)

firstOrderFeatures.enableFeatureByName('MeanIntensity', True)
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


#
# Show GLCM features
#
glcmFeatures = glcm.RadiomicsGLCM(image, mask, 10)
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
