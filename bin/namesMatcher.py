'''
This is a helper function to extract classes/features names from both
Matlab saved results and those supported by pyradiomics, and save them
into indexed text files of name format <baseline|pyradiomics>_<feature_class>.txt.
The extraction process aims to separate the prefixes related to the
pre-processing done on the images from the feature calculation.

The files have format <feature_number>:<feature_name>.

The idea is next to create manually files of name format
baseline2pyradiomics_<feature_class>.txt with the content structured as

<baseline_feature_number>:<pyradiomics_feature_number>
'''


from radiomics import base, firstorder, glcm, imageoperations, shape, glrlm, glszm
import sys, os

dataDir = os.path.join(os.path.dirname(os.path.abspath(__file__)),"..","data")
matlabFeaturesFile = os.path.join(dataDir,"MatlabBaselineFeatures.csv")
addBaselineFeaturesFile = os.path.join(dataDir, "AdditionalBaselineFeatures.csv")
outputDir = os.path.join(dataDir,'mapping')

# mapping from the feature class names used in Matlab into pyradiomics class
# names. Note that Matlab implementation does not name feature classes
# consistently, thus multiple Matlab feature classes may map into the same
# pyradiomics feature class (at least this was the assumption of @fedorov)
classMap = {"firstorder":"firstorder","Shape":"shape","RLGL":"glrlm","GLCM":"glcm","GLSZM":"glszm","rlgl":"glrlm","glcm":"glcm","glszm":"glszm","Stats":"firstorder","shape":"shape"}

# Read header from Matlab baseline features file
# If present, update with header from Additional baseline features file
baselineFeaturesNamesList = set(open(matlabFeaturesFile,'r').readline()[:-1].split(','))
if os.path.exists(addBaselineFeaturesFile):
  baselineFeaturesNamesList.update(set(open(addBaselineFeaturesFile, 'r').readline()[:-1].split(',')))

# parse names from baseline csv file(s)
baselineNames = {}
baseline2pyradiomics_features = {}
for n in baselineFeaturesNamesList:
  # try to figure out the class name first
  classNameFound = False
  for baselineClassName in classMap.keys():
    if n.find(baselineClassName)>= 0:
      n = n.split(baselineClassName+'_')[1]
      classNameFound = True
      break

  # log a watning if none of the classes that are known were recognized
  if not classNameFound:
    print "Failed to recognize class from the feature code",n
    continue

  # save this particular set of feature name mappings
  baseline2pyradiomics_features[baselineClassName] = classMap[baselineClassName]

  category = baselineClassName
  name = n

  try:
    baselineNames[classMap[category]].add(name)
  except:
    baselineNames[classMap[category]] = set()
    baselineNames[classMap[category]].add(name)

# save the feature class mapping
baseline2pyradiomics_featuresFile = os.path.join(outputDir, 'baseline2pyradiomics_featureClasses.txt')
f = open(baseline2pyradiomics_featuresFile,'w')
for m,p in baseline2pyradiomics_features.iteritems():
  f.write(m+":"+p+"\n")
f.close()

# populate pyradiomics names
pyradiomicsNames = {}
pyradiomicsNames['firstorder'] = firstorder.RadiomicsFirstOrder.getFeatureNames()
pyradiomicsNames['glcm'] = glcm.RadiomicsGLCM.getFeatureNames()
pyradiomicsNames['glszm'] = glszm.RadiomicsGLSZM.getFeatureNames()
pyradiomicsNames['shape'] = shape.RadiomicsShape.getFeatureNames()
pyradiomicsNames['glrlm'] = glrlm.RadiomicsGLRLM.getFeatureNames()

# save feature names for individual classes, indexed
#  with the keys to be used for mapping
for k,v in baselineNames.iteritems():
  i = 0
  featureListFile = os.path.join(outputDir,'baseline_'+k+'.txt')
  f = open(featureListFile,'w')
  mappingFile = os.path.join(outputDir,'baseline2pyradiomics_'+k+'.txt')
  m = open(mappingFile,'w')
  for fn in v:
    f.write(str(i)+":"+fn+"\n")
    i = i+1
    # if there is an exact match, save the mapping
    if fn in pyradiomicsNames[k]:
      m.write(str(list(baselineNames[k]).index(fn)))
      m.write(':')
      m.write(str(list(pyradiomicsNames[k]).index(fn)))
      m.write('\n')

  m.close()
  f.close()

for k,v in pyradiomicsNames.iteritems():
  i = 0
  featureListFile = os.path.join(outputDir,'pyradiomics_'+k+'.txt')
  f = open(featureListFile,'w')
  for fn in v:
    f.write(str(i)+":"+fn+"\n")
    i = i+1
f.close()
