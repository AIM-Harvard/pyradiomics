import SimpleITK as sitk
import numpy, operator
import logging

def binImage(binwidth, parameterValues, parameterMatrix, parameterMatrixCoordinates):
  lowBound = min(parameterValues) - binwidth
  highBound = max(parameterValues) + binwidth

  binedges = numpy.union1d(numpy.arange(0,lowBound,-binwidth), numpy.arange(0,highBound,binwidth))
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
