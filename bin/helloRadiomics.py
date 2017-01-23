import os
import logging
import SimpleITK as sitk
from radiomics import featureextractor
import radiomics

testCase = 'brain1'
dataDir = os.path.join(os.path.abspath(""), "..", "data")
imageName = os.path.join(dataDir, testCase + '_image.nrrd')
maskName = os.path.join(dataDir, testCase + '_label.nrrd')

if not os.path.exists(imageName):
  print 'Error: problem finding input image', imageName
  exit()
if not os.path.exists(maskName):
  print 'Error: problem finding input labelmap', maskName
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
# extractor.enableInputImages(Original={}, LoG={}, Wavelet={})

# Disable all classes except firstorder
extractor.disableAllFeatures()

# Enable all features in firstorder
# extractor.enableFeatureClassByName('firstorder')

# Only enable mean and skewness in firstorder
extractor.enableFeaturesByName(firstorder=['Mean', 'Skewness'])

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

print "Active features:"
for cls, features in extractor.enabledFeatures.iteritems():
  if len(features) == 0:
    features = extractor.getFeatureNames(cls)
  for f in features:
    print f
    print eval('extractor.featureClasses["%s"].get%sFeatureValue.__doc__' % (cls, f))

print "Calculating features"
featureVector = extractor.execute(imageName, maskName)

for featureName in featureVector.keys():
  print "Computed %s: %s" % (featureName, featureVector[featureName])
