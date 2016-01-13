from itertools import chain
import inspect, importlib

import numpy
import SimpleITK as sitk
import pywt

from radiomics import base, imageoperations, firstorder, glcm, rlgl, glszm

class RadiomicsWavelet(base.RadiomicsFeaturesBase):

  def __init__(self, inputImage, inputMask, **kwargs):
    super(RadiomicsWavelet,self).__init__(inputImage,inputMask,**kwargs)
    
    self.featureClasses = ['firstorder','glszm'] #['glcm','rlgl']
    self.imageArray = sitk.GetArrayFromImage(inputImage)
    self.maskArray = sitk.GetArrayFromImage(inputMask)   
    
    self.decompositionLevel = 1
    self.waveletFilter = 'coif1' #'haar'      
    self.approximation, self.decompositions = self.swt3(self.imageArray, self.waveletFilter, level=self.decompositionLevel)
    
    for decompositionName, decompositionArray in self.decompositions[0].iteritems():
      decompositionImage = sitk.GetImageFromArray(decompositionArray)
      decompositionImage.CopyInformation(inputImage)
      for featureClass in self.featureClasses:
        featureModule = importlib.import_module("radiomics.%s" % featureClass)
        featureClass_RadiomicsObj_ = inspect.getmembers(featureModule, \
                                      lambda member: (inspect.isclass(member)) \
                                      and ('Radiomics' in member.__name__))[0][1]  
        
        # new binwidth will be used for each decomposition to maintain consistent bin count
        (matrix, matrixCoordinates) = \
          imageoperations.padTumorMaskToCube(decompositionArray,self.maskArray)
        targetVoxelArray = matrix[matrixCoordinates]
        matrix_binned, histogram = imageoperations.binImage(self.binWidth, targetVoxelArray, matrix, matrixCoordinates)
        binCount = len(histogram[0])
    
        dec_binWidth = numpy.ptp(matrix[matrixCoordinates])/(binCount-1)
        dec_kwargs = {}
        dec_kwargs.update(kwargs)
        dec_kwargs['binWidth'] = dec_binWidth
        
        featureClass_instance = featureClass_RadiomicsObj_(decompositionImage, self.inputMask, **dec_kwargs)
        featureClass_instance.enableAllFeatures()
        
        methods = [method for method in inspect.getmembers(featureClass_instance, inspect.ismethod) \
                    if (method[0].startswith('get')) \
                    and (method[0].endswith('FeatureValue'))]
                    
        for method in methods:
          featureName = method[0].lstrip('get').rstrip('FeatureValue')
          featureFunction = method[1]
          self.add_dynamo(decompositionName, featureClass, featureName, featureFunction)
    
    super(RadiomicsWavelet,self).__init__(inputImage,inputMask,**kwargs)
  
  def swt3(self, matrix, wavelet, level=1, start_level=0):
    matrix = numpy.asarray(matrix)
    data = matrix.copy()
    if data.ndim != 3:
      raise ValueError("Expected 3D data array")
    
    original_shape = matrix.shape
    adjusted_shape = tuple([dim+1 if dim % 2 != 0 else dim for dim in original_shape])
    data.resize(adjusted_shape)    

    if not isinstance(wavelet, pywt.Wavelet):
      wavelet = pywt.Wavelet(wavelet)

    ret = []
    for i in range(start_level, start_level + level):
      H, L = self.decompose_i(data, wavelet)

      HH, HL = self.decompose_j(H, wavelet)
      LH, LL = self.decompose_j(L, wavelet)
      
      HHH, HHL = self.decompose_k(HH, wavelet)
      HLH, HLL = self.decompose_k(HL, wavelet)
      LHH, LHL = self.decompose_k(LH, wavelet)
      LLH, LLL = self.decompose_k(LL, wavelet)

      approximation = LLL.copy()
      approximation.resize(original_shape)
      dec = {'HHH': HHH,
             'HHL': HHL,
             'HLH': HLH,
             'HLL': HLL,
             'LHH': LHH,
             'LHL': LHL,
             'LLH': LLH}
      for decName, decImage in dec.iteritems(): 
        decTemp = decImage.copy()
        decTemp.resize(original_shape)
        dec[decName] = decTemp  
      ret.append(dec)
    return approximation, ret    
  
  def decompose_i(self, data, wavelet):
    #process in i:
    H, L = [], []
    i_arrays = chain.from_iterable(numpy.transpose(data,(0,1,2)))
    for i_array in i_arrays:
      cA, cD = pywt.swt(i_array, wavelet, level=1, start_level=0)[0]
      H.append(cD)
      L.append(cA)
    H = numpy.hstack(H).reshape(data.shape)
    L = numpy.hstack(L).reshape(data.shape)
    return H, L  

  def decompose_j(self, data, wavelet):
    #process in j:
    H, L = [], []
    j_arrays = chain.from_iterable(numpy.transpose(data,(0,1,2)))
    for j_array in j_arrays: 
      cA, cD = pywt.swt(j_array, wavelet, level=1, start_level=0)[0]
      H.append(cD)
      L.append(cA)
    H = numpy.asarray( [slice.T for slice in numpy.split(numpy.vstack(H), data.shape[0])] ) 
    L = numpy.asarray( [slice.T for slice in numpy.split(numpy.vstack(L), data.shape[0])] ) 
    return H, L
    
  def decompose_k(self, data, wavelet):
    #process in k:
    H, L = [], []
    k_arrays = chain.from_iterable(numpy.transpose(data,(1,2,0)))
    for k_array in k_arrays: 
      cA, cD = pywt.swt(k_array, wavelet, level=1, start_level=0)[0]
      H.append(cD)
      L.append(cA)
    H = numpy.dstack(H).reshape(data.shape)
    L = numpy.dstack(L).reshape(data.shape)
    return H, L
   
  def add_dynamo(self, decName, featureClass, featureName, featureFunction):
    # can receive a parameter dict as argument to use in docstring?
    waveletFeatureName = "wavelet_%s_%s_%s" %(decName, featureClass, featureName)
    
    def getFeatureValue_dynamo():   
      try:
        featureValue = featureFunction()
      except:
        featureValue = "Failed"
      finally:  
        return featureValue 
       
    getFeatureValue_dynamo.__doc__ = \
      """After a single-level, non-decimated, 3D discrete wavelet transform is applied on the """ \
      """input image, the coefficients for the %s decomposition are reconstructed""" \
      """into a new image with a inverse discrete wavelet transform and %s.%s is computed""" \
      """on the tumor region:\n%s""" %(decName, featureClass, featureName, featureFunction.__doc__)  
    getFeatureValue_dynamo.__name__ = "get%sFeatureValue" % waveletFeatureName
    setattr(self, getFeatureValue_dynamo.__name__, getFeatureValue_dynamo)  