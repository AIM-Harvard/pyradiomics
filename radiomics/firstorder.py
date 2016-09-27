import numpy
import collections
from radiomics import base, imageoperations
import SimpleITK as sitk

class RadiomicsFirstOrder(base.RadiomicsFeaturesBase):
  r"""
  First-order statistics describe the distribution of voxel intensities within the image region defined by the mask through commonly used and basic metrics.
  Let :math:`\textbf{X}` denote the three dimensional image matrix with :math:`N` voxels and :math:`\textbf{P}` the first order histogram with :math:`N_l` discrete intensity levels, where :math:`l` is defined by the number of levels is calculated based on the binWidth parameter of the constructor.

  Based on the definitions above, the following first order statistics can be extracted.
  """

  def __init__(self, inputImage, inputMask, **kwargs):
    super(RadiomicsFirstOrder,self).__init__(inputImage,inputMask,**kwargs)

    self.pixelSpacing = inputImage.GetSpacing()

  def _moment(self, a, moment=1, axis=0):
    r"""
    Calculate n-order moment of an array for a given axis
    """

    if moment == 1:
      return numpy.float64(0.0)
    else:
      mn = numpy.mean(a,axis, keepdims= True)
      s = numpy.power((a-mn), moment)
      return numpy.mean(s, axis)

  def getEnergyFeatureValue(self):
    r"""
    Calculate the Energy of the image array.

    :math:`energy = \displaystyle\sum^{N}_{i=1}{\textbf{X}(i)^2}`

    Energy is a measure of the magnitude of voxel values in
    an image. A larger values implies a greater sum of the
    squares of these values.
    """

    shiftedParameterArray = self.targetVoxelArray + 2000
    return (numpy.sum(shiftedParameterArray**2))

  def getTotalEnergyFeatureValue(self):
    r"""
    Calculate the Total Energy of the image array.

    :math:`total\ energy = V_{voxel}\displaystyle\sum^{N}_{i=1}{\textbf{X}(i)^2}`

    Total Energy is the value of Energy feature scaled by the volume of the voxel in cubic mm.
    """

    shiftedParameterArray = self.targetVoxelArray + 2000
    cubicMMPerVoxel = reduce(lambda x,y: x*y , self.pixelSpacing)
    return(cubicMMPerVoxel*numpy.sum(shiftedParameterArray**2))

  def getEntropyFeatureValue(self):
    r"""
    Calculate the Entropy of the image array.

    :math:`entropy = -\displaystyle\sum^{N_l}_{i=1}{\textbf{P}(i)\log_2\textbf{P}(i)}`

    Entropy specifies the uncertainty/randomness in the
    image values. It measures the average amount of
    information required to encode the image values.
    """

    eps = numpy.spacing(1)

    bins = imageoperations.getHistogram(self.binWidth, self.targetVoxelArray)[0]
    bins = bins + eps
    bins = bins/float(bins.sum())
    return (-1.0 * numpy.sum(bins*numpy.log2(bins)))

  def getMinimumFeatureValue(self):
    r"""
    Calculate the Minimum Value in the image array.
    """

    return (numpy.min(self.targetVoxelArray))

  def getMaximumFeatureValue(self):
    r"""
    Calculate the Maximum Value in the image array.
    """

    return (numpy.max(self.targetVoxelArray))

  def getMeanFeatureValue(self):
    r"""
    Calculate the Mean Value for the image array.

    :math:`mean = \frac{1}{N}\displaystyle\sum^{N}_{i=1}{\textbf{X}(i)}`
    """

    return (numpy.mean(self.targetVoxelArray))

  def getMedianFeatureValue (self):
    r"""
    Calculate the Median Value for the image array.
    """

    return (numpy.median(self.targetVoxelArray))

  def getRangeFeatureValue (self):
    r"""
    Calculate the Range of Values in the image array.

    :math:`range = \max(X) - \min(X)`
    """

    return (numpy.max(self.targetVoxelArray) - numpy.min(self.targetVoxelArray))

  def getMeanDeviationFeatureValue(self):
    r"""
    Calculate the Mean Deviation for the image array.

    :math:`mean\ deviation = \frac{1}{N}\displaystyle\sum^{N}_{i=1}{|\textbf{X}(i)-\bar{X}|}`

    Mean Deviation is the mean distance of all intensity values
    from the Mean Value of the image array.
    """

    return ( numpy.mean(numpy.absolute( (numpy.mean(self.targetVoxelArray) - self.targetVoxelArray) )) )

  def getRootMeanSquaredFeatureValue(self):
    r"""
    Calculate the Root Mean Squared of the image array.

    :math:`RMS = \sqrt{\frac{1}{N}\sum^{N}_{i=1}{\textbf{X}(i)^2}}`

    RMS is the square-root of the mean of all the squared
    intensity values. It is another measure of the magnitude
    of the image values.
    """

    shiftedParameterArray = self.targetVoxelArray + 2000
    return ( numpy.sqrt((numpy.sum(shiftedParameterArray**2))/float(shiftedParameterArray.size)) )

  def getStandardDeviationFeatureValue(self):
    r"""
    Calculate the Standard Deviation of the image array.

    :math:`standard\ deviation = \sqrt{\frac{1}{N-1}\sum^{N}_{i=1}{(\textbf{X}(i)-\bar{X})^2}}`

    Standard Deviation measures the amount of variation
    or dispersion from the Mean Value.
    """

    return (numpy.std(self.targetVoxelArray, ddof= 1))

  def getSkewnessFeatureValue(self, axis=0):
    r"""
    Calculate the Skewness of the image array.

    :math:`\gamma_1 = \frac{\frac{1}{N}\sum^{N}_{i=1}{(\textbf{X}(i)-\bar{X}})^3}{(\sum^{N}_{i=1}{(\textbf{X}(i)-\bar{X}})^2)^3}`

    Skewness measures the asymmetry of the distribution of
    values about the Mean value. Depending
    on where the tail is elongated and the mass of the distribution
    is concentrated, this value can be positive or negative.

    Related links:

    https://en.wikipedia.org/wiki/Skewness
    """

    m2 = self._moment(self.targetVoxelArray, 2, axis)
    m3 = self._moment(self.targetVoxelArray, 3, axis)

    if (m2 == 0): return numpy.core.nan

    return m3 / m2**1.5

  def getKurtosisFeatureValue(self, axis=0):
    r"""
    Calculate the Kurtosis of the image array.

    :math:`kurtosis = \frac{\frac{1}{N}\sum^{N}_{i=1}{(\textbf{X}(i)-\bar{X}})^4}{(\sum^{N}_{i=1}{(\textbf{X}(i)-\bar{X}})^2)^2}`

    Kurtosis is a measure of the 'peakedness' of the distribution
    of values in the image ROI. A higher kurtosis implies that the
    mass of the distribution is concentrated towards the tail(s)
    rather than towards the mean. A lower kurtosis implies the reverse:
    that the mass of the distribution is concentrated towards a
    spike near the Mean value.

    Related links:

    https://en.wikipedia.org/wiki/Kurtosis
    """

    m2 = self._moment(self.targetVoxelArray,2,axis)
    m4 = self._moment(self.targetVoxelArray,4,axis)

    if (m2==0): return numpy.core.nan

    return m4 / m2**2.0

  def getVarianceFeatureValue(self):
    r"""
    Calculate the Variance in the image array.

    :math:`variance = \frac{1}{N-1}\displaystyle\sum^{N}_{i=1}{(\textbf{X}(i)-\bar{X})^2}`

    Variance is the the mean of the squared distances of each intensity
    value from the Mean value. This is a measure of the spread
    of the distribution about the mean..
    """

    return (numpy.std(self.targetVoxelArray, ddof= 1)**2)

  def getUniformityFeatureValue(self):
    r"""
    Calculate the Uniformity of the image array.

    :math:`uniformity = \displaystyle\sum^{N_l}_{i=1}{\textbf{P}(i)^2}`

    Uniformity is a measure of the sum of the squares of each intensity
    value. This is a measure of the heterogeneity of the image array,
    where a greater uniformity implies a greater heterogeneity or a
    greater range of discrete intensity values.
    """

    bins = imageoperations.getHistogram(self.binWidth, self.targetVoxelArray)[0]
    try:
      bins = bins/(float(bins.sum()))
      return (numpy.sum(bins**2))
    except ZeroDivisionError:
      return numpy.core.nan
