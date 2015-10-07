from __future__ import print_function

import os
import numpy
import collections
import itertools
import operator
import SimpleITK as sitk

import RadiomicsPlatform.RadiomicsImageArrayLib
import RadiomicsPlatform.RadiomicsFeaturesLib
import pdb

"""
'Settings' is dict object containing the following data:

Settings["levels"]: 
  A List of the target label values as strings in the 
  label-map image i.e. [1] or [1,8]. 
  TODO: Radiomics features will be computed separately for each label value.
  
Settings["resampledpixelspacing"]:
  A 3-tuple representing the final pixel-spacing to resample the image
  prior to computing Radiomics features i.e. (1,1,1) or (3,3,3).
  If this field is set to the bool value, False, there will be no resampling
  and the original (base) pixel-spacing will be used.

Settings["binwidth"]:
  An integer representing the bin width (a range of pixel intensity values). Default is 25.
  Prior to computing features, The voxel intensity values within the tumor region 
  will be assigned to a bin index, with bin edges centered at the 0 value.
    
Settings["modality"]:
  A string representing the modality of the image i.e. "CT", "MRI", "PET". 

Settings["imagefilepath"]:
  A string representing the absolute filepath of the image file
  i.e. "C://Dataset1//234056//07112015_Study1//Reconstructions//CT_Lung.nrrd"

Settings["labelfilepath"]:
  A string representing the absolute filepath of the label-map file
  i.e. "C://Dataset1//234056//07112015_Study1//Segmentations//CT_Lung_IndexedLabelmap.nrrd"
  
Settings["basepixelspacing"]:
  A 3-tuple representing the original pixel spacing of the image i.e. (1.1, 1.1, 3.87).

Settings["dimensions"]
  A 3-tuple representing the original dimensions of the image in (i,j,k) i.e. (512,512,34).

Settings["seriesdescription"]:
  A string representing the series description for the image. 
  Ideally, the series description is the basename of the image filepath.
  i.e. "C://Dataset1//234056//07112015_Study1//Reconstructions//CT_Lung.nrrd"
    Series Description is "CT_Lung"
  
Settings["studydate"]:
  A string representing the study date for the image. 
  Ideally, the study date is the basename of the directory two directories up from the image file.
  i.e. "C://Dataset1//234056//07112015_Study1//Reconstructions//CT_Lung.nrrd"
    The study date is "07112015_Study1"
    
Settings["patientid"]:
  A string representing the PatientID for the image. 
  Ideally, the PatientID is the basename of the directory three directories up from the image file.
  i.e. "C://Dataset1//234056//07112015_Study1//Reconstructions//CT_Lung.nrrd"
    The PatientID is "234056"
"""

class RadiomicsPreProcessing:

    def __init__(self, Settings):
        self.Settings = {}
        self.Settings.update(Settings)
        self.RadiomicsFeatureVector = collections.OrderedDict()
        
        if self.Settings["resampledpixelspacing"]:
            pass
            #TODO: SimpleITK interpolation
            #self.imageNode = Interpolate3D.interpolateScalarVolumeNode3D(selImageNode, outputPixelSpacing=Settings["resampledpixelspacing"])
            #self.labelNode = Interpolate3D.interpolateScalarVolumeNode3D(selLabelNode, outputPixelSpacing=Settings["resampledpixelspacing"])
        else:
            self.Settings["resampledpixelspacing"] = self.Settings["basepixelspacing"]
         
        # create Numpy Arrays 
        self.imageNodeArray = RadiomicsPlatform.RadiomicsImageArrayLib.ImageToNumpyArray(self.Settings["imagefilepath"])
        self.labelNodeArray = RadiomicsPlatform.RadiomicsImageArrayLib.ImageToNumpyArray(self.Settings["labelfilepath"])
        
        # create a rectangular matrix with the dimensions of the tumor
        # center the tumor voxels in the matrix and set the non-tumor voxels to a pad value 
        self.matrix, self.matrixCoordinates = RadiomicsPlatform.RadiomicsImageArrayLib.PadTumorMaskToCube(self.imageNodeArray, self.labelNodeArray) 
        
        """
        from __main__ import vtk, qt, ctk, slicer
        self.InitializeProgressBar(self.Settings["seriesdescription"])
        """
        
    def ComputeRadiomicsFeatures(self):
        # Node Information     
        #self.UpdateProgressBar(self.Settings["seriesdescription"], "General Features")         
        self.generalFeatures = RadiomicsPlatform.RadiomicsFeaturesLib.Radiomics_General(self.Settings)
        self.RadiomicsFeatureVector.update( self.generalFeatures.EvaluateFeatures() )         
        
        # First Order Statistics    
        #self.UpdateProgressBar(self.Settings["seriesdescription"], "First Order Statistics")
        self.firstOrderStatistics = RadiomicsPlatform.RadiomicsFeaturesLib.Radiomics_First_Order(self.matrix, self.matrixCoordinates, self.Settings["binwidth"], self.Settings["resampledpixelspacing"])
        self.RadiomicsFeatureVector.update( self.firstOrderStatistics.EvaluateFeatures() )
        # Variance and Standard Deviation in numpy are different from Matlab values
        
        # Shape Features
        #self.UpdateProgressBar(self.Settings["seriesdescription"], "Shape Features")
        self.shapeFeatures = RadiomicsPlatform.RadiomicsFeaturesLib.Radiomics_Shape(self.matrix, self.matrixCoordinates, self.Settings["resampledpixelspacing"])
        self.RadiomicsFeatureVector.update( self.shapeFeatures.EvaluateFeatures() )
        
        # Texture Features(GLCM)
        #self.UpdateProgressBar(self.Settings["seriesdescription"], "GLCM Texture Features")      
        self.textureFeaturesGLCM = RadiomicsPlatform.RadiomicsFeaturesLib.Radiomics_GLCM(self.matrix, self.matrixCoordinates, self.Settings["binwidth"])   
        self.RadiomicsFeatureVector.update( self.textureFeaturesGLCM.EvaluateFeatures() )
          
        # Texture Features(GLRL)  
        #self.UpdateProgressBar(self.Settings["seriesdescription"], "GLRL Texture Features")
        self.textureFeaturesGLRL = RadiomicsPlatform.RadiomicsFeaturesLib.Radiomics_RLGL(self.matrix, self.matrixCoordinates, self.Settings["binwidth"])
        self.RadiomicsFeatureVector.update( self.textureFeaturesGLRL.EvaluateFeatures() )
        # Recheck feature computations
        
        """
        # Texture Features(GLSZM) 
        #self.UpdateProgressBar(self.Settings["seriesdescription"], "GLSZM Texture Features")
        self.textureFeaturesGLSZM = RadiomicsPlatform.RadiomicsFeaturesLib.TextureGLSZM(self.matrix, self.matrixCoordinates, self.targetVoxelArray)
        self.RadiomicsFeatureVector.update( self.textureFeaturesGLSZM.EvaluateFeatures() )
        """
        
        
        # Laplacian of a Gaussian Features(LoG)  
        #self.UpdateProgressBar(self.Settings["seriesdescription"], "LoG Features")
        self.laplacianOfGaussianFeatures = RadiomicsPlatform.RadiomicsFeaturesLib.LoGFeatures(self.Settings["imagefilepath"], self.Settings["labelfilepath"], self.Settings["binwidth"], self.Settings["resampledpixelspacing"])
        self.RadiomicsFeatureVector.update( self.laplacianOfGaussianFeatures.EvaluateFeatures() )
        
      
        """
        # Wavelet Features  
        #self.UpdateProgressBar(self.Settings["seriesdescription"], "Wavelet Features")
        self.waveletFeatures = RadiomicsPlatform.RadiomicsFeaturesLib.WaveletFeatures(self.matrix, self.matrixCoordinates, self.targetVoxelArray)
        self.RadiomicsFeatureVector.update( self.waveletFeatures.EvaluateFeatures() )
        """
      
        # close progress bar
        # self.UpdateProgressBar(self.Settings["seriesdescription"], "Populating Summary Table")
        # self.progressBar.close()
        # self.progressBar = None
    
    """
    def InitializeProgressBar(self, imageSeriesDescription):
        # initialize Progress Bar
        self.progressBar = qt.QProgressDialog(slicer.util.mainWindow())
        self.progressBar.minimumDuration = 0
        self.progressBar.show()
        self.progressBar.setValue(0)
        self.progressBar.setMaximum(5)
        self.progressBar.labelText = 'Calculating for %s: ' % imageSeriesDescription
      
    def UpdateProgressBar(self, nodeName, nextFeatureString):
        self.progressBar.labelText = 'Calculating %s: %s' % (nodeName, nextFeatureString)
        self.progressBar.setValue(self.progressBar.value + 1)
        slicer.app.processEvents()
    """
    
    def GetFeatureVector(self):
        return (self.RadiomicsFeatureVector)