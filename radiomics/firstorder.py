import numpy
import collections
from radiomics import base, imageoperations
import SimpleITK as sitk

class RadiomicsFirstOrder(base.RadiomicsFeaturesBase):

  def __init__(self, inputImage, inputMask, **kwargs):
    super(RadiomicsFirstOrder,self).__init__(inputImage,inputMask,**kwargs)

    if inputImage == None or inputMask == None:
      print ('ERROR FirstOrder: missing input image or mask')
      return

    self.pixelSpacing = inputImage.GetSpacing()

    #self.featureNames = self.getFeatureNames()

    self.imageArray = sitk.GetArrayFromImage(inputImage)
    self.maskArray = sitk.GetArrayFromImage(inputMask)

    (self.matrix, self.matrixCoordinates) = \
      imageoperations.padTumorMaskToCube(self.imageArray,self.maskArray)

    self.targetVoxelArray = self.matrix[self.matrixCoordinates]

    #self.InitializeFeatureVector()
    #for f in self.getFeatureNames():
    #  self.enabledFeatures[f] = True

    # TODO: add an option to instantiate the class that reuses initialization

  def _moment(self, a, moment=1, axis=0):
    """Calculate n-order moment of an array for a given axis"""
    if moment == 1:
      return numpy.float64(0.0)
    else:
      mn = numpy.expand_dims(numpy.mean(a,axis), axis)
      s = numpy.power((a-mn), moment)
      return numpy.mean(s, axis)

  def getEnergyFeatureValue(self):
    """
    Calculate the Energy of the image array.

    Energy is a measure of the magnitude of voxel values in
    an image. A larger values implies a greater sum of the
    squares of these values.
    """
    shiftedParameterArray = self.targetVoxelArray + 2000
    return (numpy.sum(shiftedParameterArray**2))

  def getTotalEnergyFeatureValue(self):
    """
    Calculate the Total Energy of the image array.

    Total Energy is a measure of the magnitude of voxel values
    and voxel volumes in an image. A larger values implies
    a greater sum of the squares of these values.
    """
    shiftedParameterArray = self.targetVoxelArray + 2000
    cubicMMPerVoxel = reduce(lambda x,y: x*y , self.pixelSpacing)
    return(cubicMMPerVoxel*numpy.sum(shiftedParameterArray**2))

  def getEntropyFeatureValue(self):
    """
    Calculate the Entropy of the image array.

    Entropy Specifies the uncertainty/randomness in the
    image values. It measures the average amount of
    information required to encode the image values
    """
    ##check for binning centered at 0
    bincount = numpy.ceil((numpy.max(self.targetVoxelArray) - numpy.min(self.targetVoxelArray))/float(self.binWidth))
    bins = numpy.histogram(self.targetVoxelArray, bins=bincount)[0]
    bins = bins/float(bins.sum())
    return (-1.0 * numpy.sum(bins*numpy.where(bins!=0,numpy.log2(bins),0)))

  def getMinimumFeatureValue(self):
    """Calculate the Minimum Value in the image array."""
    return (numpy.min(self.targetVoxelArray))

  def getMaximumFeatureValue(self):
    """Calculate the Maximum Value in the image array."""
    return (numpy.max(self.targetVoxelArray))

  def getMeanFeatureValue(self):
    """Calculate the Mean Value for the image array."""
    return (numpy.mean(self.targetVoxelArray))

  def getMedianFeatureValue (self):
    """Calculate the Median Value for the image array."""
    return (numpy.median(self.targetVoxelArray))

  def getRangeFeatureValue (self):
    """Calculate the Range of Values in the image array."""
    return (numpy.max(self.targetVoxelArray) - numpy.min(self.targetVoxelArray))

  def getMeanDeviationFeatureValue(self):
    """
    Calculate the Mean Deviation for the image array.

    Mean Deviation is the mean distance of all intensity values
    from the Mean Value of the image array.
    """
    return ( numpy.mean(numpy.absolute( (numpy.mean(self.targetVoxelArray) - self.targetVoxelArray) )) )

  def getRootMeanSquaredFeatureValue(self):
    """
    Calculate the Root Mean Squared of the image array.

    RMS is the square-root of the mean of all the squared
    intensity values. It is another measure of the magnitude
    of the image values.
    """
    shiftedParameterArray = self.targetVoxelArray + 2000
    return ( numpy.sqrt((numpy.sum(shiftedParameterArray**2))/float(shiftedParameterArray.size)) )

  def getStandardDeviationFeatureValue(self):
    """
    Calculate the Standard Deviation of the image array.

    Standard Deviation measures the amount of variation
    or dispersion from the Mean Value.
    """
    return (numpy.std(self.targetVoxelArray))

  def getSkewnessFeatureValue(self, axis=0):
    """
    Calculate the Skewness of the image array.

    Skewness measures the asymmetry of the distribution of
    values about the Mean value. Depending
    on where the tail is elongated and the mass of the distribution
    is concentrated, this value can be positive or negative.
    """
    m2 = self._moment(self.targetVoxelArray, 2, axis)
    m3 = self._moment(self.targetVoxelArray, 3, axis)
    zero = (m2 == 0)
    vals = numpy.where(zero, 0, m3 / m2**1.5)

    if vals.ndim == 0:
      return vals.item()

    return vals

  def getKurtosisFeatureValue(self, axis=0):
    """
    Calculate the Kurtosis of the image array.

    Kurtosis is a measure of the 'peakedness' of the distribution
    of values in the image ROI. A higher kurtosis implies that the
    mass of the distribution is concentrated towards the tail(s)
    rather than towards the mean. A lower kurtosis implies the reverse:
    that the mass of the distribution is concentrated towards a
    spike near the Mean value.
    """
    m2 = self._moment(self.targetVoxelArray,2,axis)
    m4 = self._moment(self.targetVoxelArray,4,axis)
    zero = (m2 == 0)
    olderr = numpy.seterr(all='ignore')

    try:
      vals = numpy.where(zero, 0, m4 / m2**2.0)
    finally:
      numpy.seterr(**olderr)
    if vals.ndim == 0:
      vals = vals.item()

    return vals

  def getVarianceFeatureValue(self):
    """
    Calculate the Variance in the image array.

    Variance is the the mean of the squared distances of each intensity
    value from the Mean value. This is a measure of the spread
    of the distribution about the mean..
    """
    return (numpy.std(self.targetVoxelArray)**2)

  def getUniformityFeatureValue(self):
    """
    Calculate the Uniformity of the image array.

    Uniformity is a measure of the sum of the squares of each intensity
    value. This is a measure of the heterogeneity of the image array,
    where a greater uniformity implies a greater heterogeneity or a
    greater range of discrete intensity values.
    """
    bincount = numpy.ceil((numpy.max(self.targetVoxelArray) - numpy.min(self.targetVoxelArray))/float(self.binWidth))
    bins = numpy.histogram(self.targetVoxelArray, bins=bincount)[0]
    bins = bins/float(bins.sum())
    return (numpy.sum(bins**2))
