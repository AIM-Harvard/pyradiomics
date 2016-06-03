import SimpleITK as sitk
import numpy, operator, pywt, logging
from itertools import chain

def binImage(binwidth, parameterValues, parameterMatrix, parameterMatrixCoordinates):
  lowBound = min(parameterValues)
  highBound = max(parameterValues) + binwidth

  binedges = numpy.arange(lowBound, highBound, binwidth)
  binedges[-1] += 1 # ensures that max(parametervalues) is binned to upper bin by numpy.digitize

  histogram = numpy.histogram(parameterValues, bins=binedges)
  parameterMatrix[parameterMatrixCoordinates] = numpy.digitize(parameterValues,binedges)

  return parameterMatrix, histogram

def imageToNumpyArray(imageFilePath):
  sitkImage = sitk.ReadImage(imageFilePath)
  sitkImageArray = sitk.GetArrayFromImage(sitkImage)
  return sitkImageArray

def padCubicMatrix(a, matrixCoordinates, dims):
  # pads matrix 'a' with zeros and resizes 'a' to a cube with dimensions increased to the next greatest power of 2
  # numpy version 1.7 has numpy.pad function
  # center coordinates onto padded matrix
  # consider padding with NaN or eps = numpy.spacing(1)
  pad = tuple(map(operator.div, tuple(map(operator.sub, dims, a.shape)), ([2,2,2])))
  matrixCoordinatesPadded = tuple(map(operator.add, matrixCoordinates, pad))
  matrix2 = numpy.zeros(dims)
  matrix2[matrixCoordinatesPadded] = a[matrixCoordinates]
  return (matrix2, matrixCoordinatesPadded)

def padTumorMaskToCube(imageNodeArray, labelNodeArray, padDistance=0):
  """
  Create a 3D matrix of the segmented region of the image based on the input label.

  Create a 3D matrix of the labelled region of the image, padded to have a
  cuboid shape equal to the ijk boundaries of the label.
  """
  targetVoxelsCoordinates = numpy.where(labelNodeArray != 0)
  ijkMinBounds = numpy.min(targetVoxelsCoordinates, 1)
  ijkMaxBounds = numpy.max(targetVoxelsCoordinates, 1)
  matrix = numpy.zeros(ijkMaxBounds - ijkMinBounds + 1)
  matrixCoordinates = tuple(map(operator.sub, targetVoxelsCoordinates, tuple(ijkMinBounds)))
  matrix[matrixCoordinates] = imageNodeArray[targetVoxelsCoordinates].astype('int64')
  if padDistance > 0:
    matrix, matrixCoordinates = padCubicMatrix(matrix, matrixCoordinates, padDistance)
  return(matrix, matrixCoordinates)

def padCubicMatrix(a, matrixCoordinates, padDistance):
  """
  Pad a 3D matrix in all dimensions with unit length equal to padDistance.

  Pads matrix 'a' with zeros and resizes 'a' to a cube with dimensions increased to
  the next greatest power of 2 and centers coordinates onto the padded matrix.
  """
  dims = tuple(map(operator.add, a.shape, tuple(numpy.tile(padDistance,3))))
  pad = tuple(map(operator.div, tuple(map(operator.sub, dims, a.shape)), ([2,2,2])))

  matrixCoordinatesPadded = tuple(map(operator.add, matrixCoordinates, pad))
  matrix2 = numpy.zeros(dims)
  matrix2[matrixCoordinatesPadded] = a[matrixCoordinates]
  return (matrix2, matrixCoordinatesPadded)

def interpolateImage(imageNode, resampledPixelSpacing, interpolator=sitk.sitkBSpline):
  """Resamples image or label to the specified pixel spacing (The default interpolator is Bspline)

  'imageNode' is a SimpleITK Object, and 'resampledPixelSpacing' is the output pixel spacing.
  Enumerator references for interpolator:
  0 - sitkNearestNeighbor
  1 - sitkLinear
  2 - sitkBSpline
  3 - sitkGaussian
  """
  rif = sitk.ResampleImageFilter()
  rif.SetOutputSpacing(resampledPixelSpacing)
  rif.SetInterpolator(interpolator)
  resampledImageNode = rif.Execute(imageNode)
  return resampledImageNode

#
# Use the SimpleITK LaplacianRecursiveGaussianImageFilter
# on the input image with the given sigmaValue and return
# the filtered image.
# If sigmaValue is not greater than zero, return the input image.
#
def applyLoG(inputImage, sigmaValue=0.5):
  if sigmaValue > 0.0:
    lrgif = sitk.LaplacianRecursiveGaussianImageFilter()
    lrgif.SetNormalizeAcrossScale(True)
    lrgif.SetSigma(sigmaValue)
    return lrgif.Execute(inputImage)
  else:
    logging.info('applyLoG: sigma must be greater than 0.0: %g', sigmaValue)
    return inputImage

def applyThreshold(inputImage, lowerThreshold, upperThreshold, insideValue=None, outsideValue=0):
  # this mode is useful to generate the mask of thresholded voxels
  if insideValue:
    tif = sitk.BinaryThresholdImageFilter()
    tif.SetInsideValue(insideValue)
    tif.SetLowerThreshold(lowerThreshold)
    tif.SetUpperThreshold(upperThreshold)
  else:
    tif = sitk.ThresholdImageFilter()
    tif.SetLower(lowerThreshold)
    tif.SetUpper(upperThreshold)
  tif.SetOutsideValue(outsideValue)
  return tif.Execute(inputImage)

def swt3(inputImage, wavelet="coif1", level=1, start_level=0):
  matrix = sitk.GetArrayFromImage(inputImage)
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
    H, L = decompose_i(data, wavelet)

    HH, HL = decompose_j(H, wavelet)
    LH, LL = decompose_j(L, wavelet)

    HHH, HHL = decompose_k(HH, wavelet)
    HLH, HLL = decompose_k(HL, wavelet)
    LHH, LHL = decompose_k(LH, wavelet)
    LLH, LLL = decompose_k(LL, wavelet)

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
      sitkImage = sitk.GetImageFromArray(decTemp)
      sitkImage.CopyInformation(inputImage)
      dec[decName] = sitkImage

    ret.append(dec)

  return approximation, ret

def decompose_i(data, wavelet):
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

def decompose_j(data, wavelet):
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

def decompose_k(data, wavelet):
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
