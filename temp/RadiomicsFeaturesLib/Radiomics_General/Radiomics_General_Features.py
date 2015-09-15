import os
import numpy
import collections

def imageNodeFileName(Settings):
    return(os.path.basename(Settings["imagefilepath"]))
  
def labelNodeFileName(Settings):
    return(os.path.basename(Settings["labelfilepath"]))
  
def basePixelSpacing(Settings):
    return(Settings["basepixelspacing"])
  
def resampledPixelSpacing(Settings):
    return(Settings["resampledpixelspacing"])
    
def baseDimensions(Settings):
    return(Settings["dimensions"])
    
def binWidth(Settings):
    return(Settings["binwidth"])
  
def labelLevels(Settings):
    return(Settings["levels"])
  
def imageModality(Settings):
    return(Settings["modality"])
    
def patientID(Settings):
    return(Settings["patientid"])
    
def studyDate(Settings):
    return(Settings["studydate"])
    
def seriesDescription(Settings):
    return(Settings["seriesdescription"])