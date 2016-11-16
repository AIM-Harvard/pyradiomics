import os
import logging
import SimpleITK as sitk
from radiomics import featureextractor
import radiomics

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
extractor = featureextractor.RadiomicsFeaturesExtractor(**kwargs)

# By default, only original is enabled. Optionally enable some filters:
# extractor.enableInputImages(original={}, log={}, wavelet={})

# Disable all classes except firstorder
extractor.disableAllFeatures()

# Enable all features in firstorder
# extractor.enableFeatureClassByName('firstorder')

# Only enable mean and skewness in firstorder
extractor.enableFeaturesByName(firstorder= ['Mean', 'Skewness'])

# Enable writing out the log using radiomics logger
radiomics.debug()  # Switch on radiomics logging from level=DEBUG (default level=WARNING)

# Prevent radiomics logger from printing out log entries with level < WARNING to the console
logger = logging.getLogger('radiomics')
logger.handlers[0].setLevel(logging.WARNING)
logger.propagate = False  # Do not pass log messages on to root logger

# Write out all log entries to a file
handler = logging.FileHandler(filename='testLog.txt', mode='w')
formatter = logging.Formatter("%(levelname)s:%(name)s: %(message)s")
handler.setFormatter(formatter)
logger.addHandler(handler)

for f in extractor.getFeaturesNames('firstorder'):
  print "Feature: %s" %(f)
  print eval('extractor.featureClasses["firstorder"].get'+f+'FeatureValue.__doc__')

print "Calculating features"
featureVector = extractor.execute(imageName, maskName)

for featureName in featureVector.keys():
  print "Computed %s: %s" %(featureName, featureVector[featureName])
