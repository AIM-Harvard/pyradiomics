import numpy
from radiomics import base, imageoperations


class RadiomicsFirstOrder(base.RadiomicsFeaturesBase):
  r"""
  First-order statistics describe the distribution of voxel intensities within the image region defined by the mask through commonly used and basic metrics.

  Let:

  :math:`\textbf{X}` denote the three dimensional image matrix with :math:`N` voxels

  :math:`\textbf{P}(i)` the first order histogram with :math:`N_l` discrete intensity levels,
  where :math:`l` is defined by the number of levels is calculated based on the binWidth parameter of the constructor.

  :math:`p(i)` be the normalized first order histogram and equal to :math:`\frac{\textbf{P}(i)}{\sum{\textbf{P}(i)}}`

  Following addiotional settings are possible:

  - voxelArrayShift [2000]: This amount is added to the gray level intensity in Energy, Total Energy and RMS, this is to prevent negative values from occuring when using CT data.
  """

  def __init__(self, inputImage, inputMask, **kwargs):
    super(RadiomicsFirstOrder, self).__init__(inputImage, inputMask, **kwargs)

    self.pixelSpacing = inputImage.GetSpacing()
    self.voxelArrayShift = kwargs.get('voxelArrayShift', 2000)

  def _moment(self, a, moment=1, axis=0):
    r"""
    Calculate n-order moment of an array for a given axis
    """

    if moment == 1:
      return numpy.float64(0.0)
    else:
      mn = numpy.mean(a, axis, keepdims=True)
      s = numpy.power((a - mn), moment)
      return numpy.mean(s, axis)

  def getEnergyFeatureValue(self):
    r"""
    Calculate the Energy of the image array.

    :math:`energy = \displaystyle\sum^{N}_{i=1}{\textbf{X}(i)^2}`

    Energy is a measure of the magnitude of voxel values in
    an image. A larger values implies a greater sum of the
    squares of these values.
    """

    shiftedParameterArray = self.targetVoxelArray + self.voxelArrayShift
    return (numpy.sum(shiftedParameterArray ** 2))

  def getTotalEnergyFeatureValue(self):
    r"""
    Calculate the Total Energy of the image array.

    :math:`total\ energy = V_{voxel}\displaystyle\sum^{N}_{i=1}{\textbf{X}(i)^2}`

    Total Energy is the value of Energy feature scaled by the volume of the voxel in cubic mm.
    """
    x, y, z = self.pixelSpacing
    cubicMMPerVoxel = x * y * z
    return (cubicMMPerVoxel * self.getEnergyFeatureValue())

  def getEntropyFeatureValue(self):
    r"""
    Calculate the Entropy of the image array.

    :math:`entropy = -\displaystyle\sum^{N_l}_{i=1}{p(i)\log_2\big(p(i)+\epsilon\big)}`

    Entropy specifies the uncertainty/randomness in the
    image values. It measures the average amount of
    information required to encode the image values.
    """

    eps = numpy.spacing(1)

    bins = imageoperations.getHistogram(self.binWidth, self.targetVoxelArray)[0]

    try:
      bins = bins + eps
      bins = bins / float(bins.sum())
      return (-1.0 * numpy.sum(bins * numpy.log2(bins)))
    except ZeroDivisionError:
      return numpy.core.nan

  def getMinimumFeatureValue(self):
    r"""
    Calculate the Minimum Value in the image array.
    """

    return (numpy.min(self.targetVoxelArray))

  def get10PercentileFeatureValue(self):
    r"""
    Calculate the 10\ :sup:`th` percentile in the image array.
    """

    return (numpy.percentile(self.targetVoxelArray, 10))

  def get90PercentileFeatureValue(self):
    r"""
    Calculate the 90\ :sup:`th` percentile in the image array.
    """

    return (numpy.percentile(self.targetVoxelArray, 90))

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

  def getMedianFeatureValue(self):
    r"""
    Calculate the Median Value for the image array.
    """

    return (numpy.median(self.targetVoxelArray))

  def getInterquartileRangeFeatureValue(self):
    r"""
    Calculate the interquartile range of the image array.

    :math:`interquartile\ range = \textbf{P}_{75} - \textbf{P}_{25}`, where :math:`\textbf{P}_{25}` and
    :math:`\textbf{P}_{75}` are the 25\ :sup:`th` and 75\ :sup:`th` percentile of the image array, respectively.
    """

    return numpy.percentile(self.targetVoxelArray, 75) - numpy.percentile(self.targetVoxelArray, 25)

  def getRangeFeatureValue(self):
    r"""
    Calculate the Range of Values in the image array.

    :math:`range = \max(\textbf{X}) - \min(\textbf{X})`
    """

    return (numpy.max(self.targetVoxelArray) - numpy.min(self.targetVoxelArray))

  def getMeanAbsoluteDeviationFeatureValue(self):
    r"""
    Calculate the Mean Absolute Deviation for the image array.

    :math:`mean\ absolute\ deviation = \frac{1}{N}\displaystyle\sum^{N}_{i=1}{|\textbf{X}(i)-\bar{X}|}`

    Mean Absolute Deviation is the mean distance of all intensity values
    from the Mean Value of the image array.
    """

    return (numpy.mean(numpy.absolute((numpy.mean(self.targetVoxelArray) - self.targetVoxelArray))))

  def getRobustMeanAbsoluteDeviationFeatureValue(self):
    r"""
    Calculate the Robust Mean Absolute Deviation for the image array.

    :math:`robust\ mean\ absolute\ deviation = \frac{1}{N_{10-90}}\displaystyle\sum^{N_{10-90}}_{i=1}{|\textbf{X}_{10-90}(i)-\bar{X}_{10-90}|}`

    Robust Mean Absolute Deviation is the mean distance of all intensity values
    from the Mean Value calculated on the subset of image array with gray levels in between, or equal
    to the 10\ :sub:`th` and 90\ :sub:`th` percentile.
    """
    prcnt10 = self.get10PercentileFeatureValue()
    prcnt90 = self.get90PercentileFeatureValue()
    percentileArray = self.targetVoxelArray[(self.targetVoxelArray >= prcnt10) * (self.targetVoxelArray <= prcnt90)]
    return numpy.mean(numpy.absolute(percentileArray - numpy.mean(percentileArray)))

  def getRootMeanSquaredFeatureValue(self):
    r"""
    Calculate the Root Mean Squared of the image array.

    :math:`RMS = \sqrt{\frac{1}{N}\sum^{N}_{i=1}{\textbf{X}(i)^2}}`

    RMS is the square-root of the mean of all the squared
    intensity values. It is another measure of the magnitude
    of the image values.
    """

    shiftedParameterArray = self.targetVoxelArray + self.voxelArrayShift
    return (numpy.sqrt((numpy.sum(shiftedParameterArray ** 2)) / float(shiftedParameterArray.size)))

  def getStandardDeviationFeatureValue(self):
    r"""
    Calculate the Standard Deviation of the image array.

    :math:`standard\ deviation = \sqrt{\frac{1}{N}\sum^{N}_{i=1}{(\textbf{X}(i)-\bar{X})^2}}`

    Standard Deviation measures the amount of variation
    or dispersion from the Mean Value.
    """

    return (numpy.std(self.targetVoxelArray))

  def getSkewnessFeatureValue(self, axis=0):
    r"""
    Calculate the Skewness of the image array.

    :math:`skewness = \displaystyle\frac{\mu_3}{\sigma^3}
    = \frac{\frac{1}{N}\sum^{N}_{i=1}{(\textbf{X}(i)-\bar{X})^3}}
    {\left(\sqrt{\frac{1}{N}\sum^{N}_{i=1}{(\textbf{X}(i)-\bar{X})^2}}\right)^3}`

    Where :math:`\mu_3` is the 3\ :sup:`rd` central moment.

    Skewness measures the asymmetry of the distribution of values about the Mean value. Depending
    on where the tail is elongated and the mass of the distribution
    is concentrated, this value can be positive or negative.

    Related links:

    https://en.wikipedia.org/wiki/Skewness
    """

    m2 = self._moment(self.targetVoxelArray, 2, axis)
    m3 = self._moment(self.targetVoxelArray, 3, axis)

    if (m2 == 0): return numpy.core.nan

    return m3 / m2 ** 1.5

  def getKurtosisFeatureValue(self, axis=0):
    r"""
    Calculate the Kurtosis of the image array.

    :math:`kurtosis = \displaystyle\frac{\mu_4}{\sigma^4}
    = \frac{\frac{1}{N}\sum^{N}_{i=1}{(\textbf{X}(i)-\bar{X})^4}}
    {\left(\frac{1}{N}\sum^{N}_{i=1}{(\textbf{X}(i)-\bar{X}})^2\right)^2}`

    Where :math:`\mu_4` is the 4\ :sup:`th` central moment.

    Kurtosis is a measure of the 'peakedness' of the distribution
    of values in the image ROI. A higher kurtosis implies that the
    mass of the distribution is concentrated towards the tail(s)
    rather than towards the mean. A lower kurtosis implies the reverse:
    that the mass of the distribution is concentrated towards a
    spike near the Mean value.

    Related links:

    https://en.wikipedia.org/wiki/Kurtosis
    """

    m2 = self._moment(self.targetVoxelArray, 2, axis)
    m4 = self._moment(self.targetVoxelArray, 4, axis)

    if (m2 == 0): return numpy.core.nan

    return m4 / m2 ** 2.0

  def getVarianceFeatureValue(self):
    r"""
    Calculate the Variance in the image array.

    :math:`variance = \sigma^2 = \frac{1}{N}\displaystyle\sum^{N}_{i=1}{(\textbf{X}(i)-\bar{X})^2}`

    Variance is the the mean of the squared distances of each intensity
    value from the Mean value. This is a measure of the spread
    of the distribution about the mean.
    """

    return (numpy.std(self.targetVoxelArray) ** 2)

  def getUniformityFeatureValue(self):
    r"""
    Calculate the Uniformity of the image array.

    :math:`uniformity = \displaystyle\sum^{N_l}_{i=1}{p(i)^2}`

    Uniformity is a measure of the sum of the squares of each intensity
    value. This is a measure of the heterogeneity of the image array,
    where a greater uniformity implies a greater heterogeneity or a
    greater range of discrete intensity values.
    """

    bins = imageoperations.getHistogram(self.binWidth, self.targetVoxelArray)[0]
    try:
      bins = bins / (float(bins.sum()))
      return (numpy.sum(bins ** 2))
    except ZeroDivisionError:
      return numpy.core.nan
