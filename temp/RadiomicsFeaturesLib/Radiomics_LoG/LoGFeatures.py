import numpy
import collections
import SimpleITK as sitk

import RadiomicsPlatform.RadiomicsImageArrayLib
import pdb

class LoGFeatures:

    def __init__(self, imageFilePath, labelFilePath, binwidth, pixelSpacing):
        self.sigmaValues = numpy.arange(5.0, 0.0, -0.5)[::-1]

        self.imageFilePath = imageFilePath
        self.labelFilePath = labelFilePath
        self.binwidth = binwidth
        self.pixelSpacing = pixelSpacing
        
        self.sitkImageNode = sitk.ReadImage(imageFilePath)
        self.sitkLabelNode = sitk.ReadImage(labelFilePath)
        self.labelNodeArray = sitk.GetArrayFromImage(self.sitkLabelNode)  
        
        #self.bincount = numpy.ceil((numpy.max(self.parameterValues) - numpy.min(self.parameterValues))/float(self.binwidth))
        #self.cubicMMPerVoxel = reduce(lambda x,y: x*y , self.pixelSpacing)
        
        self.InitializeFeatureVector()
    
    def InitializeFeatureVector(self):
        self.prefix = "LoG_"
        self.laplacian_gaussian_FeatureVector = collections.OrderedDict()
        
        for sigma in self.sigmaValues:
            matrix_LoGFiltered, matrixCoordinates_LoGFiltered = self.ApplyLoGFilter(self.sitkImageNode, self.labelNodeArray, sigma)
            
            qwert = matrix_LoGFiltered.copy()
            
            pdb.set_trace()
            
            try:
                LoGFeatureVector = collections.OrderedDict()
                
                # entropyPos, meanPos, uniformityPos use the positive values in the filtered image array only
                # filteredImageValuesPos = filteredImageValues[filteredImageValues>=0]
                # later, when computing
            
                LoGFirstOrderStatistics = RadiomicsPlatform.RadiomicsFeaturesLib.Radiomics_First_Order(matrix_LoGFiltered, matrixCoordinates_LoGFiltered, self.binwidth, self.pixelSpacing)
                LoGFeatureVector.update( LoGFirstOrderStatistics.EvaluateFeatures() )
                
                LoGTextureFeaturesGLCM = RadiomicsPlatform.RadiomicsFeaturesLib.Radiomics_GLCM(matrix_LoGFiltered, matrixCoordinates_LoGFiltered, self.binwidth)   
                LoGFeatureVector.update( LoGTextureFeaturesGLCM.EvaluateFeatures() )
            
                LoGTextureFeaturesGLRL = RadiomicsPlatform.RadiomicsFeaturesLib.Radiomics_RLGL(matrix_LoGFiltered, matrixCoordinates_LoGFiltered, self.binwidth)
                LoGFeatureVector.update( LoGTextureFeaturesGLRL.EvaluateFeatures() )
            
                #LoGTextureFeaturesGLSZM = RadiomicsPlatform.RadiomicsFeaturesLib.TextureGLSZM(matrix_LoGFiltered, matrixCoordinates_LoGFiltered, self.binwidth)
                #LoGFeatureVector.update( LoGTextureFeaturesGLSZM.EvaluateFeatures() )
            except IndexError:
                continue
            
            for radiomicsLoGFeature in LoGFeatureVector:
                self.laplacian_gaussian_FeatureVector[self.prefix + str(sigma).replace(".","_") + "_mm_3D_" + radiomicsLoGFeature] = LoGFeatureVector[radiomicsLoGFeature]
    
    def ApplyLoGFilter(self, sitkImageNode, labelNodeArray, sigma):
        LoGFilter = sitk.LaplacianRecursiveGaussianImageFilter()
        LoGFilter.SetSigma(sigma)
        
        sitkImageNode_LoGFiltered = LoGFilter.Execute(sitkImageNode)        
        imageNodeArray_LoGFiltered = sitk.GetArrayFromImage(sitkImageNode_LoGFiltered)    
        matrix_LoGFiltered, matrixCoordinates_LoGFiltered = RadiomicsPlatform.RadiomicsImageArrayLib.PadTumorMaskToCube(imageNodeArray_LoGFiltered, labelNodeArray)

        return matrix_LoGFiltered, matrixCoordinates_LoGFiltered
    
    def entropyValue(self, parameterArray, bincount):
        bins = numpy.histogram(parameterArray, bins=bincount)[0]
        bins = bins/float(bins.sum())
        return (-1.0 * numpy.sum(bins*numpy.where(bins!=0,numpy.log2(bins),0)))
      
        
    def meanIntensity(self, parameterArray):
        return (numpy.mean(parameterArray))
        
    
    def standardDeviation(self, parameterArray):
        return (numpy.std(parameterArray))
    
    def _moment(self, a, moment=1, axis=0):
        if moment == 1:
            return numpy.float64(0.0)
        else:
            mn = numpy.expand_dims(numpy.mean(a,axis), axis)
            s = numpy.power((a-mn), moment)
            return numpy.mean(s, axis)

    def skewnessValue(self, a, axis=0):
        m2 = self._moment(a, 2, axis)
        m3 = self._moment(a, 3, axis)
        
        # Control Flow: if m2==0 then vals = 0; else vals = m3/m2**1.5
        zero = (m2 == 0)
        vals = numpy.where(zero, 0, m3 / m2**1.5)
        
        if vals.ndim == 0:
            return vals.item()
            
        return vals

    def kurtosisValue(self, a, axis=0):
        m2 = self._moment(a,2,axis)
        m4 = self._moment(a,4,axis)
        zero = (m2 == 0)
        
        # Set Floating-Point Error Handling
        olderr = numpy.seterr(all='ignore')
        try:
            vals = numpy.where(zero, 0, m4 / m2**2.0)
        finally:
            numpy.seterr(**olderr)
        if vals.ndim == 0:
            vals = vals.item() # array scalar

        return vals
        
    
    def uniformityValue(self, parameterArray):
        bins = numpy.histogram(parameterArray, bins=self.bincount)[0]
        bins = bins/float(bins.sum())
        return (numpy.sum(bins**2))
      
    def EvaluateFeatures(self):
        for feature in self.laplacian_gaussian_FeatureVector:
            try:
                self.laplacian_gaussian_FeatureVector[feature] = str(self.laplacian_gaussian_FeatureVector[feature])
            except AttributeError:
                self.laplacian_gaussian_FeatureVector[feature] = "Function Does Not Exist"
                
        return(self.laplacian_gaussian_FeatureVector)  
              
