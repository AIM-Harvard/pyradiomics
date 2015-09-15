import os
import numpy
import collections

import Radiomics_General_Features

class Radiomics_General:

    def __init__(self, Settings):    
        self.Settings = Settings
        
        self.InitializeFeatureVector()
      
    def InitializeFeatureVector(self):
        self.prefix = "General_"
        self.radiomics_general_FeatureVector = collections.OrderedDict()
        self.radiomics_general_FeatureVector[self.prefix+"patientID"] = "Radiomics_General_Features.patientID(self.Settings)"
        self.radiomics_general_FeatureVector[self.prefix+"studyDate"] = "Radiomics_General_Features.studyDate(self.Settings)"
        self.radiomics_general_FeatureVector[self.prefix+"seriesDescription"] = "Radiomics_General_Features.seriesDescription(self.Settings)"
        self.radiomics_general_FeatureVector[self.prefix+"imageFile"] = "Radiomics_General_Features.imageNodeFileName(self.Settings)"
        self.radiomics_general_FeatureVector[self.prefix+"labelFile"] = "Radiomics_General_Features.labelNodeFileName(self.Settings)"
        self.radiomics_general_FeatureVector[self.prefix+"imageModality"] = "Radiomics_General_Features.imageModality(self.Settings)"
        self.radiomics_general_FeatureVector[self.prefix+"basePixelSpacing"] = "Radiomics_General_Features.basePixelSpacing(self.Settings)"
        self.radiomics_general_FeatureVector[self.prefix+"resampledPixelSpacing"] = "Radiomics_General_Features.resampledPixelSpacing(self.Settings)"
        self.radiomics_general_FeatureVector[self.prefix+"baseDimensions"] = "Radiomics_General_Features.baseDimensions(self.Settings)"
        self.radiomics_general_FeatureVector[self.prefix+"labelLevels"] = "Radiomics_General_Features.labelLevels(self.Settings)"
        self.radiomics_general_FeatureVector[self.prefix+"binWidth"] = "Radiomics_General_Features.binWidth(self.Settings)"
      
    def EvaluateFeatures(self):
        for feature in self.radiomics_general_FeatureVector:
            try:
                self.radiomics_general_FeatureVector[feature] = eval(self.radiomics_general_FeatureVector[feature])
            except AttributeError:
                self.radiomics_general_FeatureVector[feature] = "Function Does Not Exist"
                
        return(self.radiomics_general_FeatureVector)
