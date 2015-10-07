import numpy
import collections
import pdb

import RadiomicsPlatform.RadiomicsImageArrayLib
import Radiomics_RLGL_Matrix
import Radiomics_RLGL_Features
import Radiomics_RLGL_Coefficients

class Radiomics_RLGL:

    def __init__(self,  matrix, matrixCoordinates, binwidth):
        self.matrix = matrix
        self.matrixCoordinates = matrixCoordinates
        self.binwidth = binwidth
        
        self.targetVoxelArray = self.matrix[self.matrixCoordinates]
        self.Coefficients = {}
        
        # binning
        self.matrix, self.histogram = RadiomicsPlatform.RadiomicsImageArrayLib.BinImage(self.binwidth, self.targetVoxelArray, self.matrix, self.matrixCoordinates)
        self.Coefficients['Ng'] = len(self.histogram[0])
        self.Coefficients['grayLevels'] = numpy.linspace(1,self.Coefficients['Ng'],num=self.Coefficients['Ng'])
        self.Coefficients['Nr'] = numpy.max(self.matrix.shape)
        self.Coefficients['Np'] = self.targetVoxelArray.size
        
        self.P_rlgl = Radiomics_RLGL_Matrix.CreateRLGL(self.Coefficients['Ng'], self.Coefficients['Nr'], self.Coefficients['grayLevels'], self.matrix, self.matrixCoordinates, angles=13)
        
        self.Coefficients.update(Radiomics_RLGL_Coefficients.CalculateCoefficients(self.P_rlgl))
        self.InitializeFeatureVector()

    def InitializeFeatureVector(self):
        self.prefix = "RLGL_"
        self.radiomics_rlgl_FeatureVector = collections.OrderedDict()
        self.radiomics_rlgl_FeatureVector[self.prefix+"shortRunEmphasis"] = "Radiomics_RLGL_Features.shortRunEmphasis(self.Coefficients['pr'], self.Coefficients['jvector'], self.Coefficients['sumP_rlgl'])"
        self.radiomics_rlgl_FeatureVector[self.prefix+"longRunEmphasis"] = "Radiomics_RLGL_Features.longRunEmphasis(self.Coefficients['pr'], self.Coefficients['jvector'], self.Coefficients['sumP_rlgl'])"
        self.radiomics_rlgl_FeatureVector[self.prefix+"grayLevelNonuniformity"] = "Radiomics_RLGL_Features.grayLevelNonUniformity(self.Coefficients['pg'], self.Coefficients['sumP_rlgl'])"
        self.radiomics_rlgl_FeatureVector[self.prefix+"runLengthNonuniformity"] = "Radiomics_RLGL_Features.runLengthNonUniformity(self.Coefficients['pr'], self.Coefficients['sumP_rlgl'])"
        self.radiomics_rlgl_FeatureVector[self.prefix+"runPercentage"] = "Radiomics_RLGL_Features.runPercentage(self.P_rlgl, self.Coefficients['Np'])"
        
        #These do not match Matlab values
        self.radiomics_rlgl_FeatureVector[self.prefix+"lowGrayLevelRunEmphasis"] = "Radiomics_RLGL_Features.lowGrayLevelRunEmphasis(self.Coefficients['pg'], self.Coefficients['ivector'], self.Coefficients['sumP_rlgl'])"
        self.radiomics_rlgl_FeatureVector[self.prefix+"highGrayLevelRunEmphasis"] = "Radiomics_RLGL_Features.highGrayLevelRunEmphasis(self.Coefficients['pg'], self.Coefficients['ivector'], self.Coefficients['sumP_rlgl'])"
        self.radiomics_rlgl_FeatureVector[self.prefix+"shortRunLowGrayLevEmpha"] = "Radiomics_RLGL_Features.shortRunLowGrayLevelEmphasis(self.P_rlgl, self.Coefficients['ivector'], self.Coefficients['jvector'], self.Coefficients['sumP_rlgl'])"
        self.radiomics_rlgl_FeatureVector[self.prefix+"shortRunHighGrayLevEmpha"] = "Radiomics_RLGL_Features.shortRunHighGrayLevelEmphasis(self.P_rlgl, self.Coefficients['ivector'], self.Coefficients['jvector'], self.Coefficients['sumP_rlgl'])"
        self.radiomics_rlgl_FeatureVector[self.prefix+"longRunLowGrayLevEmpha"] = "Radiomics_RLGL_Features.longRunLowGrayLevelEmphasis(self.P_rlgl, self.Coefficients['ivector'], self.Coefficients['jvector'], self.Coefficients['sumP_rlgl'])"
        self.radiomics_rlgl_FeatureVector[self.prefix+"longRunHighGrayLevEmpha"] = "Radiomics_RLGL_Features.longRunHighGrayLevelEmphasis(self.P_rlgl, self.Coefficients['ivector'], self.Coefficients['jvector'], self.Coefficients['sumP_rlgl'])"
           
    def EvaluateFeatures(self):  
        for feature in self.radiomics_rlgl_FeatureVector:
            try:
                self.radiomics_rlgl_FeatureVector[feature] = eval(self.radiomics_rlgl_FeatureVector[feature])
            except AttributeError:
                self.radiomics_rlgl_FeatureVector[feature] = "Function Does Not Exist"
              
        return(self.radiomics_rlgl_FeatureVector)

