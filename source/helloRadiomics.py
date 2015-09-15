from radiomics import firstorder, preprocessing
import SimpleITK as sitk
import sys

imageName = '/Users/fedorov/github/radiomics-platform/Testing/Data/prostate_phantom.nrrd'
maskName = '/Users/fedorov/github/radiomics-platform/Testing/Data/prostate_phantom_label.nrrd'

#imageName = sys.argv[1]
#maskName = sys.argv[2]

image = sitk.ReadImage(imageName)
mask = sitk.ReadImage(maskName)

firstOrderFeatures = firstorder.RadiomicsFirstOrder(image,mask)
firstOrderFeatures.setBinWidth(10)

firstOrderFeatures.enableFeatureByName('MeanIntensity', True)
firstOrderFeatures.enableAllFeatures()

print 'Will calculate the following features: '
for f in firstOrderFeatures.enabledFeatures.keys():
  print '  ',f
  print eval('firstOrderFeatures.get'+f+'FeatureValue.__doc__')

print 'Calculating...',
firstOrderFeatures.calculateFeatures()
print 'done'

print 'Calculated features: '
for (key,val) in firstOrderFeatures.featureValues.iteritems():
  print '  ',key,':',val
