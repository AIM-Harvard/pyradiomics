'''
This is a helper function to extract classes/features names from both
Matlab saved results and those supported by pyradiomics, and save them
into indexed text files of name format <matlab|pyradiomics>_<feature_class>.txt.
The extraction process aims to separate the prefixes related to the
pre-processing done on the images from the feature calculation.

The files have format <feature_number>:<feature_name>.

The idea is next to create manually files of name format
matlab2pyradiomics_<feature_class>.txt with the content structured as

<matlab_feature_number>:<pyradiomics_feature_number>
'''


from radiomics import base, firstorder, glcm, preprocessing, shape, rlgl
import sys, os

dataDir = os.path.join(os.path.dirname(os.path.abspath(__file__)),"..","data")
matlabFeaturesFile = os.path.join(dataDir,"MatlabFeatures.csv")
outputDir = os.path.join(dataDir,'mapping')

# mapping from the feature class names used in Matlab into pyradiomics class
# names. Note that Matlab implementation does not name feature classes
# consistently, thus multiple Matlab feature classes may map into the same
# pyradiomics feature class (at least this was the assumption of @fedorov)
classMap = {"stats":"firstorder","Shape":"shape","RLGL":"rlgl","GLCM":"glcm","rlgl":"rlgl","glcm":"glcm","Stats":"firstorder"}

matlabFeaturesNamesList = open(matlabFeaturesFile,'r').readline()[:-1].split(',')

# parse names from Matlab csv file
matlabNames = {}
for n in matlabFeaturesNamesList:
  # try to figure out the class name first
  classNameFound = False
  for matlabClassName in classMap.keys():
    if n.find(matlabClassName)>= 0:
      n = n.split(matlabClassName+'_')[1]
      classNameFound = True
      break

  # log a watning if none of the classes that are known were recognized
  if not classNameFound:
    print "Failed to recognize class from the feature code",n
    continue

  category = matlabClassName
  name = n

  try:
    matlabNames[classMap[category]].add(name)
  except:
    matlabNames[classMap[category]] = set()
    matlabNames[classMap[category]].add(name)


# populate pyradiomics names
pyradiomicsNames = {}
pyradiomicsNames['firstorder'] = firstorder.RadiomicsFirstOrder.getFeatureNames()
pyradiomicsNames['glcm'] = glcm.RadiomicsGLCM.getFeatureNames()
pyradiomicsNames['shape'] = shape.RadiomicsShape.getFeatureNames()
pyradiomicsNames['rlgl'] = rlgl.RadiomicsRLGL.getFeatureNames()

# save feature names for individual classes, indexed
#  with the keys to be used for mapping
for k,v in matlabNames.iteritems():
  i = 0
  featureListFile = os.path.join(outputDir,'matlab_'+k+'.txt')
  f = open(featureListFile,'w')
  for fn in v:
    f.write(str(i)+":"+fn+"\n")
    i = i+1
f.close()

for k,v in pyradiomicsNames.iteritems():
  i = 0
  featureListFile = os.path.join(outputDir,'pyradiomics_'+k+'.txt')
  f = open(featureListFile,'w')
  for fn in v:
    f.write(str(i)+":"+fn+"\n")
    i = i+1
f.close()
