from radiomics import signatures
import SimpleITK as sitk
import sys, os

dataDir = os.path.dirname(os.path.abspath(__file__)) + os.path.sep + ".." + os.path.sep + "data"
#imageName = str(dataDir + os.path.sep + 'prostate_phantom_subvolume.nrrd')
#maskName = str(dataDir + os.path.sep + 'prostate_phantom_subvolume-label.nrrd')
imageName = str(dataDir + os.path.sep + 'brain1_image.nrrd')
maskName = str(dataDir + os.path.sep + 'brain1_label.nrrd')

if not os.path.exists(imageName):
  print 'Error: problem finding input image',imageName
  exit()
if not os.path.exists(maskName):
  print 'Error: problem finding input labelmap',maskName
  exit()

# Define settings for signature calculation
# These are currently set equal to the respective default values
kwargs = {}
kwargs['binWidth'] = 25
kwargs['resampledPixelSpacing'] = None  # [3,3,3] is an example for defining resampling (voxels with size 3x3x3mm)
kwargs['interpolator'] = sitk.sitkBSpline
kwargs['verbose'] = True

# Initialize wrapperClass to generate signature
Sig = signatures.RadiomicsSignature(**kwargs)

# Only enable the original image as input (disables LoG, Wavelet)
Sig.enableInputImages(original= {})

# Disable all classes except firstorder
Sig.disableAllFeatures()

# Enable all features in firstorder
# Sig.enableFeatureClassByName('firstorder')

# Only enable mean and skewness in firstorder
Sig.enableFeaturesByName(firstorder= ['Mean', 'Skewness'])

for f in Sig.getFeaturesFromClass('firstorder'):
  print "Feature: %s" %(f)
  print eval('Sig.featureClasses["firstorder"].get'+f+'FeatureValue.__doc__')

print "Calculating features"
featureVector = Sig.computeSignature(imageName, maskName)

for featureName in featureVector.keys():
  print "Computed %s: %s" %(featureName, featureVector[featureName])
