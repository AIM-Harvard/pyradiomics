import numpy
import operator
import collections

import RadiomicsPlatform.RadiomicsImageArrayLib
import Radiomics_Shape_Features

class Radiomics_Shape:

    def __init__(self, matrix, matrixCoordinates, pixelSpacing):
        self.matrix = matrix
        self.matrixCoordinates = matrixCoordinates
        self.pixelSpacing = pixelSpacing
      
        self.targetVoxelArray = self.matrix[self.matrixCoordinates]
        self.cubicMMPerVoxel = reduce(lambda x,y: x*y , self.pixelSpacing)
        
        # Pad tumor matrix by 10 voxels in each of the three dimensions
        self.maxDimsSA = tuple(map(operator.add, self.matrix.shape, ([10,10,10])))
        self.matrixSA, self.matrixSACoordinates = RadiomicsPlatform.RadiomicsImageArrayLib.PadCubicMatrix(self.matrix, self.matrixCoordinates, self.maxDimsSA)
        
        # Volume and Surface Area are pre-calculated
        self.Volume = Radiomics_Shape_Features.volumeMM3(self.targetVoxelArray, self.cubicMMPerVoxel)
        self.SurfaceArea = Radiomics_Shape_Features.surfaceArea(self.matrixSA, self.matrixSACoordinates, self.targetVoxelArray, self.pixelSpacing)
        
        self.InitializeFeatureVector()
      
    def InitializeFeatureVector(self):
        self.prefix = "Shape_"
        self.radiomics_shape_FeatureVector = collections.OrderedDict()
        self.radiomics_shape_FeatureVector[self.prefix+"voxelNum"] = 'Radiomics_Shape_Features.voxelNumber(self.targetVoxelArray)'
        self.radiomics_shape_FeatureVector[self.prefix+"volume"] = 'self.Volume'
        self.radiomics_shape_FeatureVector[self.prefix+"surface"] = 'self.SurfaceArea'
        self.radiomics_shape_FeatureVector[self.prefix+"surfVolRatio"] = 'Radiomics_Shape_Features.surfaceVolumeRatio(self.SurfaceArea, self.Volume)'
        self.radiomics_shape_FeatureVector[self.prefix+"compactness"] = 'Radiomics_Shape_Features.compactness1(self.SurfaceArea, self.Volume)'
        self.radiomics_shape_FeatureVector[self.prefix+"compactness2"] = 'Radiomics_Shape_Features.compactness2(self.SurfaceArea, self.Volume)'
        self.radiomics_shape_FeatureVector[self.prefix+"maxDiameter3D"] = 'Radiomics_Shape_Features.maximum3DDiameter(self.matrixSA, self.matrixSACoordinates, self.pixelSpacing)'
        self.radiomics_shape_FeatureVector[self.prefix+"spherDisprop"] = 'Radiomics_Shape_Features.sphericalDisproportion(self.SurfaceArea, self.Volume)'
        self.radiomics_shape_FeatureVector[self.prefix+"sphericity"] = 'Radiomics_Shape_Features.sphericityValue(self.SurfaceArea, self.Volume)'
        
        #self.radiomics_shape_FeatureVector[self.prefix+"maxDiameter3D_pdist"] = 
        #self.radiomics_shape_FeatureVector[self.prefix+"maxDiameter2Dx"] = 
        #self.radiomics_shape_FeatureVector[self.prefix+"maxDiameter2Dy"] = 
        #self.radiomics_shape_FeatureVector[self.prefix+"maxDiameter2Dz"] = 
      
    def EvaluateFeatures(self):    
        for feature in self.radiomics_shape_FeatureVector:
            try:
                self.radiomics_shape_FeatureVector[feature] = eval(self.radiomics_shape_FeatureVector[feature])
            except AttributeError:
                self.radiomics_shape_FeatureVector[feature] = "Function Does Not Exist"
                
        return(self.radiomics_shape_FeatureVector)   
