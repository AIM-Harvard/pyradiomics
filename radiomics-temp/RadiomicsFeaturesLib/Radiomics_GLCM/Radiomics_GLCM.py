import numpy
import collections
import pdb

import RadiomicsPlatform.RadiomicsImageArrayLib
import Radiomics_GLCM_Matrix
import Radiomics_GLCM_Features
import Radiomics_GLCM_Coefficients

class Radiomics_GLCM:

    def __init__(self, matrix, matrixCoordinates, binwidth):
        self.matrix = matrix
        self.matrixCoordinates = matrixCoordinates
        self.binwidth = binwidth
        
        self.targetVoxelArray = self.matrix[self.matrixCoordinates]
        self.Coefficients = {}
        
        # binning
        self.matrix, self.histogram = RadiomicsPlatform.RadiomicsImageArrayLib.BinImage(self.binwidth, self.targetVoxelArray, self.matrix, self.matrixCoordinates)
        self.Coefficients['Ng'] = len(self.histogram[0])
        
        self.P_glcm = Radiomics_GLCM_Matrix.CreateGLCM(self.Coefficients['Ng'], self.matrix, self.matrixCoordinates)
        
        self.Coefficients.update( Radiomics_GLCM_Coefficients.CalculateCoefficients(self.Coefficients['Ng'], self.P_glcm) )
        self.InitializeFeatureVector()    
      
    def InitializeFeatureVector(self):
        self.prefix = "GLCM_"
        self.radiomics_glcm_FeatureVector = collections.OrderedDict()       
        self.radiomics_glcm_FeatureVector[self.prefix+"autocorr"] = "Radiomics_GLCM_Features.autocorrelation(self.P_glcm, self.Coefficients['i'], self.Coefficients['j'])"
        self.radiomics_glcm_FeatureVector[self.prefix+"clusProm"] = "Radiomics_GLCM_Features.clusterProminence(self.P_glcm, self.Coefficients['i'], self.Coefficients['j'], self.Coefficients['ux'], self.Coefficients['uy'])"
        self.radiomics_glcm_FeatureVector[self.prefix+"clusShade"] = "Radiomics_GLCM_Features.clusterShade(self.P_glcm, self.Coefficients['i'], self.Coefficients['j'], self.Coefficients['ux'], self.Coefficients['uy'])"
        self.radiomics_glcm_FeatureVector[self.prefix+"clusTend"] = "Radiomics_GLCM_Features.clusterTendency(self.P_glcm, self.Coefficients['i'], self.Coefficients['j'], self.Coefficients['ux'], self.Coefficients['uy'])"
        self.radiomics_glcm_FeatureVector[self.prefix+"contrast"] = "Radiomics_GLCM_Features.contrast(self.P_glcm, self.Coefficients['i'], self.Coefficients['j'])"
        self.radiomics_glcm_FeatureVector[self.prefix+"correl1"] = "Radiomics_GLCM_Features.correlation(self.P_glcm, self.Coefficients['i'], self.Coefficients['j'], self.Coefficients['ux'], self.Coefficients['uy'], self.Coefficients['sigx'], self.Coefficients['sigy'])"
        self.radiomics_glcm_FeatureVector[self.prefix+"diffEntro"] = "Radiomics_GLCM_Features.differenceEntropy(self.Coefficients['pxSuby'], self.Coefficients['eps'])"
        self.radiomics_glcm_FeatureVector[self.prefix+"dissimilar"] = "Radiomics_GLCM_Features.dissimilarity(self.P_glcm, self.Coefficients['i'], self.Coefficients['j'])"
        self.radiomics_glcm_FeatureVector[self.prefix+"energy"] = "Radiomics_GLCM_Features.energy(self.P_glcm)"
        self.radiomics_glcm_FeatureVector[self.prefix+"entrop2"] = "Radiomics_GLCM_Features.entropy(self.P_glcm, self.Coefficients['eps'])" 
        self.radiomics_glcm_FeatureVector[self.prefix+"homogeneity1"] = "Radiomics_GLCM_Features.homogeneity1(self.P_glcm, self.Coefficients['i'], self.Coefficients['j'])"
        self.radiomics_glcm_FeatureVector[self.prefix+"homogeneity2"] = "Radiomics_GLCM_Features.homogeneity2(self.P_glcm, self.Coefficients['i'], self.Coefficients['j'])"
        self.radiomics_glcm_FeatureVector[self.prefix+"infoCorr1"] = "Radiomics_GLCM_Features.imc1(self.Coefficients['HXY'], self.Coefficients['HXY1'], self.Coefficients['HX'], self.Coefficients['HY'])"
        self.radiomics_glcm_FeatureVector[self.prefix+"invDiffmomnor"] = "Radiomics_GLCM_Features.idmn(self.P_glcm, self.Coefficients['i'], self.Coefficients['j'], self.Coefficients['Ng'])"
        self.radiomics_glcm_FeatureVector[self.prefix+"invDiffnor"] = "Radiomics_GLCM_Features.idn(self.P_glcm, self.Coefficients['i'], self.Coefficients['j'], self.Coefficients['Ng'])"
        self.radiomics_glcm_FeatureVector[self.prefix+"inverseVar"] = "Radiomics_GLCM_Features.inverseVariance(self.P_glcm, self.Coefficients['i'], self.Coefficients['j'], self.Coefficients['Ng'])"
        self.radiomics_glcm_FeatureVector[self.prefix+"maxProb"] = "Radiomics_GLCM_Features.maximumProbability(self.P_glcm)"
        self.radiomics_glcm_FeatureVector[self.prefix+"sumAvg"] = "Radiomics_GLCM_Features.sumAverage(self.Coefficients['pxAddy'], self.Coefficients['kValuesSum'])"
        self.radiomics_glcm_FeatureVector[self.prefix+"sumEntro"] = "Radiomics_GLCM_Features.sumEntropy(self.Coefficients['pxAddy'], self.Coefficients['eps'])" 
        self.radiomics_glcm_FeatureVector[self.prefix+"sumVar"] = "Radiomics_GLCM_Features.sumVariance(self.Coefficients['pxAddy'], self.Coefficients['kValuesSum'])"   
        self.radiomics_glcm_FeatureVector[self.prefix+"sumSquares"] = "Radiomics_GLCM_Features.sumSquares(self.P_glcm, self.Coefficients['i'], self.Coefficients['u'])"
        
        #self.radiomics_glcm_FeatureVector[self.prefix+"infoCorr2"] = "sum(imc2)/len(imc2)" #"self.imc2(self,)"  # produces a calculation error
      
    def EvaluateFeatures(self):
        for feature in self.radiomics_glcm_FeatureVector:
            try:
                self.radiomics_glcm_FeatureVector[feature] = eval(self.radiomics_glcm_FeatureVector[feature])
            except AttributeError:
                self.radiomics_glcm_FeatureVector[feature] = "Function Does Not Exist"
                
        return(self.radiomics_glcm_FeatureVector)  
   
