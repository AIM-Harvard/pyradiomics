import numpy
import collections
import pdb

import RadiomicsPlatform.RadiomicsImageArrayLib
import Radiomics_GLSZM_Coefficients
import Radiomics_GLSZM_Features
import Radiomics_GLSZM_Matrix

class TextureGLSZM:

    def __init__(self, matrix, matrixCoordinates, binwidth):
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
        
        #self.P_glszm = Radiomics_GLSZM_Matrix.CreateGLSZM(self.Coefficients['Ng'], self.Coefficients['Nr'], self.Coefficients['grayLevels'], self.matrix, self.matrixCoordinates, angles=13)
        
        #self.Coefficients.update(Radiomics_RLGL_Coefficients.CalculateCoefficients(self.P_rlgl))
        self.InitializeFeatureVector()
        
        self.CalculateCoefficients()
      
    def InitializeFeatureVector(self):
        self.prefix = "GLSZM_"
        self.radiomics_glszm_FeatureVector = collections.OrderedDict()
        self.radiomics_glszm_FeatureVector[self.prefix+"smallAreaEmphasis"] = "self.smallAreaEmphasis(self.P_glszm, self.jvector, self.sumP_glszm)"
        self.radiomics_glszm_FeatureVector[self.prefix+"largeAreaEmphasis"] = "self.largeAreaEmphasis(self.P_glszm, self.jvector, self.sumP_glszm)"
        self.radiomics_glszm_FeatureVector[self.prefix+"intensityVariability"] = "self.intensityVariability(self.P_glszm, self.sumP_glszm)"
        self.radiomics_glszm_FeatureVector[self.prefix+"sizeZoneVariability"] = "self.sizeZoneVariability(self.P_glszm, self.sumP_glszm)"
        self.radiomics_glszm_FeatureVector[self.prefix+"zonePercentage"] = "self.zonePercentage(self.P_glszm, self.Np)"
        
        #These do not match Matlab values
        self.radiomics_glszm_FeatureVector[self.prefix+"lowIntensityEmphasis"] = "self.lowIntensityEmphasis(self.P_glszm, self.ivector, self.sumP_glszm)"
        self.radiomics_glszm_FeatureVector[self.prefix+"highIntensityEmphasis"] = "self.highIntensityEmphasis(self.P_glszm, self.ivector, self.sumP_glszm)"
        self.radiomics_glszm_FeatureVector[self.prefix+"lowIntensitySmallAreaEmp"] = "self.lowIntensitySmallAreaEmphasis(self.P_glszm, self.ivector, self.jvector, self.sumP_glszm)"
        self.radiomics_glszm_FeatureVector[self.prefix+"highIntensitySmallAreaEmp"] = "self.highIntensitySmallAreaEmphasis(self.P_glszm, self.ivector, self.jvector, self.sumP_glszm)"
        self.radiomics_glszm_FeatureVector[self.prefix+"lowIntensityLargeAreaEmp"] = "self.lowIntensityLargeAreaEmphasis(self.P_glszm, self.ivector, self.jvector, self.sumP_glszm)"
        self.radiomics_glszm_FeatureVector[self.prefix+"highIntensityLarteAreaEmp"] = "self.highIntensityLargeAreaEmphasis(self.P_glszm, self.ivector, self.jvector, self.sumP_glszm)"
        
    
    def CalculateCoefficients(self):
        # binning
        self.binedges = numpy.union1d(numpy.arange(0,min(self.parameterValues)-self.binwidth,-self.binwidth), numpy.arange(0,max(self.parameterValues)+self.binwidth,self.binwidth))
        self.hist = numpy.histogram(self.parameterValues, bins=self.binedges)
        self.Ng = len(self.hist[0])
        self.grayLevels = numpy.linspace(1,self.Ng,num=self.Ng)
        self.parameterMatrix[self.parameterMatrixCoordinates] = numpy.digitize(self.parameterValues,self.binedges)
      
        #setup
        self.angles = 13
        self.Nr = numpy.max(self.parameterMatrix.shape)
        self.Np = self.parameterValues.size
            
        self.P_glszm = numpy.zeros((self.Ng, self.Nr, self.angles))
        self.P_glszm = self.calculate_glszm(self.grayLevels, self.parameterMatrix, self.parameterMatrixCoordinates, self.angles, self.P_glszm)
         
        self.sumP_glszm = numpy.sum( numpy.sum(self.P_glszm, 0), 0 )
        self.sumP_glszm[self.sumP_glszm==0] = 1
        
        self.pr = numpy.sum(self.P_glszm, 0)
        self.pg = numpy.sum(self.P_glszm, 1)
        
        self.ivector = numpy.arange(1, self.P_glszm.shape[0] + 1)
        self.jvector = numpy.arange(1, self.P_glszm.shape[1] + 1)
         
    def calculate_glszm(self, grayLevels, matrix, matrixCoordinates, angles, P_out):
        return P_out  
      
    def EvaluateFeatures(self):

        #Evaluate dictionary elements corresponding to user selected keys
        for key in self.keys:
            self.radiomics_glszm_FeatureVector[key] = eval(self.radiomics_glszm_FeatureVector[key])
            
        return(self.radiomics_glszm_FeatureVector)  
     
