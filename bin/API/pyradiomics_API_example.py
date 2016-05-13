# -*- coding: utf-8 -*-
import sys
import os
import glob
import collections
import csv

import traceback
import pdb

import numpy as np
#import pandas as pd
import SimpleITK as sitk

from batchprocessing import DatasetHierarchyReader
from radiomics import firstorder, glcm, imageoperations, shape, rlgl, glszm


def main():
    dataDirectory = r"/Users/nicole/Downloads/TCGA-PyRadiomics-Data/"
    inputDirectory = dataDirectory
    # Put the output in a direcotry one up from the data directory (a recursive
    # file grep is used to find data to process), creating the output path if
    # necessary.
    outputDirectory = os.path.join(dataDirectory, "..", "pyradiomics-processing")
    if not os.path.exists(outputDirectory):
      os.mkdir(outputDirectory)
      print "Created output directory ",outputDirectory
    outputFilepath = os.path.join(outputDirectory, "radiomics_features.csv")
    progress_filename = os.path.join(outputDirectory, "pyrad_log.txt")
    print "Output file = ", outputFilepath
    print "Progress file = ", progress_filename
    filetype = ".nrrd"
    
    keywordSettings = {} 
    keywordSettings['image'] = ""
    keywordSettings['imageExclusion'] = "label"
    keywordSettings['mask'] = "label"
    keywordSettings['maskExclusion'] = ""
    
    kwargs = {}
    kwargs['binWidth'] = 25
    kwargs['resampledPixelSpacing'] = None
    kwargs['interpolator'] = sitk.sitkBSpline
    kwargs['padDistance'] = 5
    kwargs['padFillValue'] = 0
    
    datasetReader = DatasetHierarchyReader(inputDirectory, filetype=filetype)   
    datasetHierarchyDict = datasetReader.ReadDatasetHierarchy()
    
    for patientIndex, patientDirectory in enumerate(datasetHierarchyDict):
        patientID = os.path.basename(patientDirectory)
        
        for studyDirectory in datasetHierarchyDict[patientDirectory]:
            studyDate = os.path.basename(studyDirectory)
            with open(progress_filename,mode='a') as printfile:
                print "(%s/%s) Processing Patient: %s, Study: %s" %(str(patientIndex+1), str(len(datasetHierarchyDict.keys())), patientID, studyDate)
                printfile.write("(%s/%s) Processing Patient: %s, Study: %s" %(str(patientIndex+1), str(len(datasetHierarchyDict.keys())), patientID, studyDate) + "\n")
                
            imageFilepaths = datasetHierarchyDict[patientDirectory][studyDirectory]["reconstructions"]
            maskFilepaths = datasetHierarchyDict[patientDirectory][studyDirectory]["segmentations"]

            imageFilepath, maskFilepath = datasetReader.findImageAndLabelPair(imageFilepaths, maskFilepaths, keywordSettings)            
            if (imageFilepath is not None) and (maskFilepath is not None):
                featureVector = collections.OrderedDict()
                featureVector['PatientID'] = patientID
                featureVector['Study'] = studyDate
                featureVector['image'] = os.path.basename(imageFilepath)
                featureVector['mask'] = os.path.basename(maskFilepath)
                
                try:
                    featureVector.update( computeFeatures(imageFilepath, maskFilepath, **kwargs) )                  
                    with open(outputFilepath, 'ab') as outputFile:
                        writer = csv.writer(outputFile, lineterminator = '\n')
                        if patientIndex==0: writer.writerow(featureVector.keys())
                        writer.writerow(featureVector.values())
                except Exception, e:
                    with open(progress_filename,mode='a') as printfile:    
                        printfile.write('\tFAILED%s\n\n' %(traceback.format_exc()))
                
def computeFeatures(ImageFilePath, MaskFilePath, **kwargs):

    if isinstance(ImageFilePath, basestring) and os.path.exists(ImageFilePath):   
        image = sitk.ReadImage(ImageFilePath)
    elif isinstance(ImageFilePath, sitk.SimpleITK.Image):
        image = ImageFilePath
    else:
        print "Error reading image Filepath or SimpleITK object"
        
    if isinstance(MaskFilePath, basestring) and os.path.exists(MaskFilePath):    
        mask = sitk.ReadImage(MaskFilePath)
    elif isinstance(ImageFilePath, sitk.SimpleITK.Image):
        mask = MaskFilePath
    else:
        print "Error reading mask Filepath or SimpleITK object"
        
      
    """
    imageArray = sitk.GetArrayFromImage(image)
    maskArray = sitk.GetArrayFromImage(mask)
    tumorVoxels = imageArray[np.where(maskArray==1)]
    mean = np.mean(tumorVoxels)
    std = np.std(tumorVoxels)
    minBound = mean - 3*std
    maxBound = mean + 3*std
    maskArray[np.where(imageArray < minBound)] = 0  
    maskArray[np.where(imageArray > maxBound)] = 0
    mask_normalized = sitk.GetImageFromArray(maskArray)
    mask_normalized.CopyInformation(mask)
    mask = mask_normalized
    """
    
    FeatureVector = collections.OrderedDict()
    
    print "\tComputing First Order"
    FeatureVector.update( computeFirstOrder(image, mask, **kwargs) )
    
    print "\tComputing Shape"
    FeatureVector.update( computeShape(image, mask, **kwargs) )
    
    print "\tComputing GLCM"
    FeatureVector.update( computeGLCM(image, mask, **kwargs) )
    
    print "\tComputing RLGL"
    FeatureVector.update( computeRLGL(image, mask, **kwargs) )
    
    print "\tComputing GLSZM"
    FeatureVector.update( computeGLSZM(image, mask, **kwargs) )
    
    print "\tComputing LoG"
    FeatureVector.update( computeLoG(image, mask, **kwargs) )
    
    print "\tComputing Wavelet"
    FeatureVector.update( computeWavelet(image, mask, **kwargs) )
    
    return FeatureVector
    
def computeFirstOrder(image, mask, **kwargs):
    firstOrderFeatures = firstorder.RadiomicsFirstOrder(image, mask, **kwargs)
    firstOrderFeatures.enableAllFeatures()
    firstOrderFeatures.calculateFeatures()

    firstOrderFeatureVector = collections.OrderedDict()
    for (featureName, featureValue) in firstOrderFeatures.featureValues.iteritems():
        firstOrderFeatureName = "firstorder_%s" %(featureName)
        firstOrderFeatureVector[firstOrderFeatureName] = featureValue
        
    return firstOrderFeatureVector 

def computeShape(image, mask, **kwargs):    
    shapeFeatures = shape.RadiomicsShape(image, mask, **kwargs)
    shapeFeatures.enableAllFeatures()
    shapeFeatures.calculateFeatures()

    shapeFeatureVector = collections.OrderedDict()
    for (featureName, featureValue) in shapeFeatures.featureValues.iteritems():
        shapeFeatureName = "shape_%s" %(featureName)
        shapeFeatureVector[shapeFeatureName] = featureValue
        
    return shapeFeatureVector 

def computeGLCM(image, mask, **kwargs):
    glcmFeatures = glcm.RadiomicsGLCM(image, mask, **kwargs)
    glcmFeatures.enableAllFeatures()     
    glcmFeatures.calculateFeatures()
    
    glcmFeatureVector = collections.OrderedDict()
    for (featureName, featureValue) in glcmFeatures.featureValues.iteritems():
        glcmFeatureName = "glcm_%s" %(featureName)
        glcmFeatureVector[glcmFeatureName] = featureValue
        
    return glcmFeatureVector 
    
def computeRLGL(image, mask, **kwargs):    
    rlglFeatures = rlgl.RadiomicsRLGL(image, mask, **kwargs)
    rlglFeatures.enableAllFeatures()
    rlglFeatures.calculateFeatures()
    
    rlglFeatureVector = collections.OrderedDict()
    for (featureName, featureValue) in rlglFeatures.featureValues.iteritems():
        rlglFeatureName = "rlgl_%s" %(featureName)
        rlglFeatureVector[rlglFeatureName] = featureValue
        
    return rlglFeatureVector 

def computeGLSZM(image, mask, **kwargs):    
    glszmFeatures = glszm.RadiomicsGLSZM(image, mask, **kwargs)
    glszmFeatures.enableAllFeatures()
    glszmFeatures.calculateFeatures()
    
    glszmFeatureVector = collections.OrderedDict()
    for (featureName, featureValue) in glszmFeatures.featureValues.iteritems():
        glszmFeatureName = "glszm_%s" %(featureName)
        glszmFeatureVector[glszmFeatureName] = featureValue
        
    return glszmFeatureVector   

def computeLoG(image, mask, **kwargs):
    logFeatureVector = collections.OrderedDict()
    
    mmif = sitk.MinimumMaximumImageFilter()
    mmif.Execute(image)
    lowerThreshold = 0
    upperThreshold = mmif.GetMaximum()

    threshImage = imageoperations.applyThreshold(image,lowerThreshold=lowerThreshold, upperThreshold=upperThreshold,outsideValue=0)
    threshImageMask = imageoperations.applyThreshold(image,lowerThreshold=lowerThreshold, upperThreshold=upperThreshold,outsideValue=0,insideValue=1)
    threshMask = sitk.Cast(mask,1) & sitk.Cast(threshImageMask,1)
    
    sigmaValues = np.arange(5.,0.,-.5)[::1]
    for sigma in sigmaValues:
        logImage = imageoperations.applyLoG(image, sigmaValue=sigma)
        logFirstorderFeatures = firstorder.RadiomicsFirstOrder(logImage, threshMask)

        logSigmaFeatureVector = collections.OrderedDict()
        logSigmaFeatureVector.update( computeFirstOrder(logImage, threshMask, **kwargs) )
        logSigmaFeatureVector.update( computeGLCM(logImage, threshMask, **kwargs) )
        logSigmaFeatureVector.update( computeRLGL(logImage, threshMask, **kwargs) )
        logSigmaFeatureVector.update( computeGLSZM(logImage, threshMask, **kwargs) )
        
        for (featureName, featureValue) in logSigmaFeatureVector.iteritems():
            laplacianFeatureName = "log_sigma_%s_mm_3D_%s" %(str(sigma).replace('.','_'), featureName)
            logFeatureVector[laplacianFeatureName] = featureValue
            
    return logFeatureVector

def computeWavelet(image, mask, **kwargs):
    waveletFeatureVector = collections.OrderedDict()
    
    apporx, ret = imageoperations.swt3(image)
    for decompositionName, decompositionImage in ret[0].items():
        #sitk.WriteImage(decompositionImage, "%s.nrrd" % decompositionName)
        #continue
        
        waveletDecompositionFeatureVector = collections.OrderedDict()
        waveletDecompositionFeatureVector.update( computeFirstOrder(decompositionImage, mask, **kwargs) )
        waveletDecompositionFeatureVector.update( computeGLCM(decompositionImage, mask, **kwargs) )
        waveletDecompositionFeatureVector.update( computeRLGL(decompositionImage, mask, **kwargs) )
        waveletDecompositionFeatureVector.update( computeGLSZM(decompositionImage, mask, **kwargs) )
        
        for (featureName, featureValue) in waveletDecompositionFeatureVector.iteritems():
            waveletFeatureName = "wavelet_%s_%s" %(decompositionName, featureName)
            waveletFeatureVector[waveletFeatureName] = featureValue
    
    return waveletFeatureVector
    
main()
