import numpy
from six.moves import range

from radiomics import base, cMatrices, deprecated


class RadiomicsFirstOrder(base.RadiomicsFeaturesBase):
  r"""
  First-order statistics describe the distribution of voxel intensities within the image region defined by the mask
  through commonly used and basic metrics.

  Let:

  - :math:`\textbf{X}` be a set of :math:`N_p` voxels included in the ROI
  - :math:`\textbf{P}(i)` be the first order histogram with :math:`N_g` discrete intensity levels,
    where :math:`N_g` is the number of non-zero bins, equally spaced from 0 with a width defined in the ``binWidth``
    parameter.
  - :math:`p(i)` be the normalized first order histogram and equal to :math:`\frac{\textbf{P}(i)}{N_p}`

  Following additional settings are possible:

  - voxelArrayShift [0]: Integer, This amount is added to the gray level intensity in features Energy, Total Energy and
    RMS, this is to prevent negative values. *If using CT data, or data normalized with mean 0, consider setting this
    parameter to a fixed value (e.g. 2000) that ensures non-negative numbers in the image. Bear in mind however, that
    the larger the value, the larger the volume confounding effect will be.*

  .. note::
    In the IBSI feature definitions, no correction for negative gray values is implemented. To achieve similar behaviour
    in PyRadiomics, set ``voxelArrayShift`` to 0.
  """

  def __init__(self, inputImage, inputMask, **kwargs):
    super(RadiomicsFirstOrder, self).__init__(inputImage, inputMask, **kwargs)

    self.pixelSpacing = inputImage.GetSpacing()
    self.voxelArrayShift = kwargs.get('voxelArrayShift', 0)
    self.discretizedImageArray = self._applyBinning(self.imageArray.copy())

  def _initVoxelBasedCalculation(self):
    super(RadiomicsFirstOrder, self)._initVoxelBasedCalculation()

    kernelRadius = self.settings.get('kernelRadius', 1)

    # Get the size of the input, which depends on whether it is in masked mode or not
    if self.masked:
      size = numpy.max(self.labelledVoxelCoordinates, 1) - numpy.min(self.labelledVoxelCoordinates, 1) + 1
    else:
      size = numpy.array(self.imageArray.shape)

    # Take the minimum size along each dimension from either the size of the ROI or the kernel
    boundingBoxSize = numpy.minimum(size, kernelRadius * 2 + 1)

    # Calculate the offsets, which can be used to generate a list of kernel Coordinates. Shape (Nd, Nk)
    self.kernelOffsets = cMatrices.generate_angles(boundingBoxSize,
                                                   numpy.array(range(1, kernelRadius + 1)),
                                                   True,  # Bi-directional
                                                   self.settings.get('force2D', False),
                                                   self.settings.get('force2Ddimension', 0))
    self.kernelOffsets = numpy.append(self.kernelOffsets, [[0, 0, 0]], axis=0)  # add center voxel
    self.kernelOffsets = self.kernelOffsets.transpose((1, 0))

    self.imageArray = self.imageArray.astype('float')
    self.imageArray[~self.maskArray] = numpy.nan
    self.imageArray = numpy.pad(self.imageArray,
                                pad_width=self.settings.get('kernelRadius', 1),
                                mode='constant', constant_values=numpy.nan)
    self.maskArray = numpy.pad(self.maskArray,
                               pad_width=self.settings.get('kernelRadius', 1),
                               mode='constant', constant_values=False)

  def _initCalculation(self, voxelCoordinates=None):

    if voxelCoordinates is None:
      self.targetVoxelArray = self.imageArray[self.maskArray].astype('float').reshape((1, -1))
      _, p_i = numpy.unique(self.discretizedImageArray[self.maskArray], return_counts=True)
      p_i = p_i.reshape((1, -1))
    else:
      # voxelCoordinates shape (Nd, Nvox)
      voxelCoordinates = voxelCoordinates.copy() + self.settings.get('kernelRadius', 1)  # adjust for padding
      kernelCoords = self.kernelOffsets[:, None, :] + voxelCoordinates[:, :, None]  # Shape (Nd, Nvox, Nk)
      kernelCoords = tuple(kernelCoords)  # shape (Nd, (Nvox, Nk))

      self.targetVoxelArray = self.imageArray[kernelCoords]  # shape (Nvox, Nk)

      p_i = numpy.empty((voxelCoordinates.shape[1], len(self.coefficients['grayLevels'])))  # shape (Nvox, Ng)
      for gl_idx, gl in enumerate(self.coefficients['grayLevels']):
        p_i[:, gl_idx] = numpy.nansum(self.discretizedImageArray[kernelCoords] == gl, 1)

    sumBins = numpy.sum(p_i, 1, keepdims=True).astype('float')
    sumBins[sumBins == 0] = 1  # Prevent division by 0 errors
    p_i = p_i.astype('float') / sumBins
    self.coefficients['p_i'] = p_i

    self.logger.debug('First order feature class initialized')

  @staticmethod
  def _moment(a, moment=1):
    r"""
    Calculate n-order moment of an array for a given axis
    """

    if moment == 1:
      return numpy.float(0.0)
    else:
      mn = numpy.nanmean(a, 1, keepdims=True)
      s = numpy.power((a - mn), moment)
      return numpy.nanmean(s, 1)

  def getEnergyFeatureValue(self):
    r"""
    **1. Energy**

    .. math::
      \textit{energy} = \displaystyle\sum^{N_p}_{i=1}{(\textbf{X}(i) + c)^2}

    Here, :math:`c` is optional value, defined by ``voxelArrayShift``, which shifts the intensities to prevent negative
    values in :math:`\textbf{X}`. This ensures that voxels with the lowest gray values contribute the least to Energy,
    instead of voxels with gray level intensity closest to 0.

    Energy is a measure of the magnitude of voxel values in an image. A larger values implies a greater sum of the
    squares of these values.

    .. note::
      This feature is volume-confounded, a larger value of :math:`c` increases the effect of volume-confounding.
    """

    shiftedParameterArray = self.targetVoxelArray + self.voxelArrayShift

    return numpy.nansum(shiftedParameterArray ** 2, 1)

  def getTotalEnergyFeatureValue(self):
    r"""
    **2. Total Energy**

    .. math::
      \textit{total energy} = V_{voxel}\displaystyle\sum^{N_p}_{i=1}{(\textbf{X}(i) + c)^2}

    Here, :math:`c` is optional value, defined by ``voxelArrayShift``, which shifts the intensities to prevent negative
    values in :math:`\textbf{X}`. This ensures that voxels with the lowest gray values contribute the least to Energy,
    instead of voxels with gray level intensity closest to 0.

    Total Energy is the value of Energy feature scaled by the volume of the voxel in cubic mm.

    .. note::
      This feature is volume-confounded, a larger value of :math:`c` increases the effect of volume-confounding.

    .. note::
      Not present in IBSI feature definitions
    """

    cubicMMPerVoxel = numpy.multiply.reduce(self.pixelSpacing)

    return self.getEnergyFeatureValue() * cubicMMPerVoxel

  def getEntropyFeatureValue(self):
    r"""
    **3. Entropy**

    .. math::
      \textit{entropy} = -\displaystyle\sum^{N_g}_{i=1}{p(i)\log_2\big(p(i)+\epsilon\big)}

    Here, :math:`\epsilon` is an arbitrarily small positive number (:math:`\approx 2.2\times10^{-16}`).

    Entropy specifies the uncertainty/randomness in the image values. It measures the average amount of information
    required to encode the image values.

    .. note::
      Defined by IBSI as Intensity Histogram Entropy.
    """
    p_i = self.coefficients['p_i']

    eps = numpy.spacing(1)
    return -1.0 * numpy.sum(p_i * numpy.log2(p_i + eps), 1)

  def getMinimumFeatureValue(self):
    r"""
    **4. Minimum**

    .. math::
      \textit{minimum} = \min(\textbf{X})
    """

    return numpy.nanmin(self.targetVoxelArray, 1)

  def get10PercentileFeatureValue(self):
    r"""
    **5. 10th percentile**

    The 10\ :sup:`th` percentile of :math:`\textbf{X}`
    """
    return numpy.nanpercentile(self.targetVoxelArray, 10, axis=1)

  def get90PercentileFeatureValue(self):
    r"""
    **6. 90th percentile**

    The 90\ :sup:`th` percentile of :math:`\textbf{X}`
    """

    return numpy.nanpercentile(self.targetVoxelArray, 90, axis=1)

  def getMaximumFeatureValue(self):
    r"""
    **7. Maximum**

    .. math::
      \textit{maximum} = \max(\textbf{X})

    The maximum gray level intensity within the ROI.
    """

    return numpy.nanmax(self.targetVoxelArray, 1)

  def getMeanFeatureValue(self):
    r"""
    **8. Mean**

    .. math::
      \textit{mean} = \frac{1}{N_p}\displaystyle\sum^{N_p}_{i=1}{\textbf{X}(i)}

    The average gray level intensity within the ROI.
    """

    return numpy.nanmean(self.targetVoxelArray, 1)

  def getMedianFeatureValue(self):
    r"""
    **9. Median**

    The median gray level intensity within the ROI.
    """

    return numpy.nanmedian(self.targetVoxelArray, 1)

  def getInterquartileRangeFeatureValue(self):
    r"""
    **10. Interquartile Range**

    .. math::
      \textit{interquartile range} = \textbf{P}_{75} - \textbf{P}_{25}

    Here :math:`\textbf{P}_{25}` and :math:`\textbf{P}_{75}` are the 25\ :sup:`th` and 75\ :sup:`th` percentile of the
    image array, respectively.
    """

    return numpy.nanpercentile(self.targetVoxelArray, 75, 1) - numpy.nanpercentile(self.targetVoxelArray, 25, 1)

  def getRangeFeatureValue(self):
    r"""
    **11. Range**

    .. math::
      \textit{range} = \max(\textbf{X}) - \min(\textbf{X})

    The range of gray values in the ROI.
    """

    return numpy.nanmax(self.targetVoxelArray, 1) - numpy.nanmin(self.targetVoxelArray, 1)

  def getMeanAbsoluteDeviationFeatureValue(self):
    r"""
    **12. Mean Absolute Deviation (MAD)**

    .. math::
      \textit{MAD} = \frac{1}{N_p}\displaystyle\sum^{N_p}_{i=1}{|\textbf{X}(i)-\bar{X}|}

    Mean Absolute Deviation is the mean distance of all intensity values from the Mean Value of the image array.
    """

    u_x = numpy.nanmean(self.targetVoxelArray, 1, keepdims=True)
    return numpy.nanmean(numpy.absolute(self.targetVoxelArray - u_x), 1)

  def getRobustMeanAbsoluteDeviationFeatureValue(self):
    r"""
    **13. Robust Mean Absolute Deviation (rMAD)**

    .. math::
      \textit{rMAD} = \frac{1}{N_{10-90}}\displaystyle\sum^{N_{10-90}}_{i=1}
      {|\textbf{X}_{10-90}(i)-\bar{X}_{10-90}|}

    Robust Mean Absolute Deviation is the mean distance of all intensity values
    from the Mean Value calculated on the subset of image array with gray levels in between, or equal
    to the 10\ :sup:`th` and 90\ :sup:`th` percentile.
    """

    prcnt10 = self.get10PercentileFeatureValue()
    prcnt90 = self.get90PercentileFeatureValue()
    percentileArray = self.targetVoxelArray.copy()

    # First get a mask for all valid voxels
    msk = ~numpy.isnan(percentileArray)
    # Then, update the mask to reflect all valid voxels that are outside the the closed 10-90th percentile range
    msk[msk] = ((percentileArray - prcnt10[:, None])[msk] < 0) | ((percentileArray - prcnt90[:, None])[msk] > 0)
    # Finally, exclude the invalid voxels by setting them to numpy.nan.
    percentileArray[msk] = numpy.nan

    return numpy.nanmean(numpy.absolute(percentileArray - numpy.nanmean(percentileArray, 1, keepdims=True)), 1)

  def getRootMeanSquaredFeatureValue(self):
    r"""
    **14. Root Mean Squared (RMS)**

    .. math::
      \textit{RMS} = \sqrt{\frac{1}{N_p}\sum^{N_p}_{i=1}{(\textbf{X}(i) + c)^2}}

    Here, :math:`c` is optional value, defined by ``voxelArrayShift``, which shifts the intensities to prevent negative
    values in :math:`\textbf{X}`. This ensures that voxels with the lowest gray values contribute the least to RMS,
    instead of voxels with gray level intensity closest to 0.

    RMS is the square-root of the mean of all the squared intensity values. It is another measure of the magnitude of
    the image values. This feature is volume-confounded, a larger value of :math:`c` increases the effect of
    volume-confounding.
    """

    # If no voxels are segmented, prevent division by 0 and return 0
    if self.targetVoxelArray.size == 0:
      return 0

    shiftedParameterArray = self.targetVoxelArray + self.voxelArrayShift
    Nvox = numpy.sum(~numpy.isnan(self.targetVoxelArray), 1).astype('float')
    return numpy.sqrt(numpy.nansum(shiftedParameterArray ** 2, 1) / Nvox)

  @deprecated
  def getStandardDeviationFeatureValue(self):
    r"""
    **15. Standard Deviation**

    .. math::
      \textit{standard deviation} = \sqrt{\frac{1}{N_p}\sum^{N_p}_{i=1}{(\textbf{X}(i)-\bar{X})^2}}

    Standard Deviation measures the amount of variation or dispersion from the Mean Value. By definition,
    :math:`\textit{standard deviation} = \sqrt{\textit{variance}}`

    .. note::
      As this feature is correlated with variance, it is marked so it is not enabled by default.
      To include this feature in the extraction, specify it by name in the enabled features
      (i.e. this feature will not be enabled if no individual features are specified (enabling 'all' features),
      but will be enabled when individual features are specified, including this feature).
      Not present in IBSI feature definitions (correlated with variance)
    """

    return numpy.nanstd(self.targetVoxelArray, axis=1)

  def getSkewnessFeatureValue(self):
    r"""
    **16. Skewness**

    .. math::
      \textit{skewness} = \displaystyle\frac{\mu_3}{\sigma^3} =
      \frac{\frac{1}{N_p}\sum^{N_p}_{i=1}{(\textbf{X}(i)-\bar{X})^3}}
      {\left(\sqrt{\frac{1}{N_p}\sum^{N_p}_{i=1}{(\textbf{X}(i)-\bar{X})^2}}\right)^3}

    Where :math:`\mu_3` is the 3\ :sup:`rd` central moment.

    Skewness measures the asymmetry of the distribution of values about the Mean value. Depending on where the tail is
    elongated and the mass of the distribution is concentrated, this value can be positive or negative.

    Related links:

    https://en.wikipedia.org/wiki/Skewness

    .. note::
      In case of a flat region, the standard deviation and 4\ :sup:`rd` central moment will be both 0. In this case, a
      value of 0 is returned.
    """

    m2 = self._moment(self.targetVoxelArray, 2)
    m3 = self._moment(self.targetVoxelArray, 3)

    m2[m2 == 0] = 1  # Flat Region, prevent division by 0 errors
    m3[m2 == 0] = 0  # ensure Flat Regions are returned as 0

    return m3 / m2 ** 1.5

  def getKurtosisFeatureValue(self):
    r"""
    **17. Kurtosis**

    .. math::
      \textit{kurtosis} = \displaystyle\frac{\mu_4}{\sigma^4} =
      \frac{\frac{1}{N_p}\sum^{N_p}_{i=1}{(\textbf{X}(i)-\bar{X})^4}}
      {\left(\frac{1}{N_p}\sum^{N_p}_{i=1}{(\textbf{X}(i)-\bar{X}})^2\right)^2}

    Where :math:`\mu_4` is the 4\ :sup:`th` central moment.

    Kurtosis is a measure of the 'peakedness' of the distribution of values in the image ROI. A higher kurtosis implies
    that the mass of the distribution is concentrated towards the tail(s) rather than towards the mean. A lower kurtosis
    implies the reverse: that the mass of the distribution is concentrated towards a spike near the Mean value.

    Related links:

    https://en.wikipedia.org/wiki/Kurtosis

    .. note::
      In case of a flat region, the standard deviation and 4\ :sup:`rd` central moment will be both 0. In this case, a
      value of 0 is returned.

    .. note::
      The IBSI feature definition implements excess kurtosis, where kurtosis is corrected by -3, yielding 0 for normal
      distributions. The PyRadiomics kurtosis is not corrected, yielding a value 3 higher than the IBSI kurtosis.
    """

    m2 = self._moment(self.targetVoxelArray, 2)
    m4 = self._moment(self.targetVoxelArray, 4)

    m2[m2 == 0] = 1  # Flat Region, prevent division by 0 errors
    m4[m2 == 0] = 0  # ensure Flat Regions are returned as 0

    return m4 / m2 ** 2.0

  def getVarianceFeatureValue(self):
    r"""
    **18. Variance**

    .. math::
      \textit{variance} = \frac{1}{N_p}\displaystyle\sum^{N_p}_{i=1}{(\textbf{X}(i)-\bar{X})^2}

    Variance is the the mean of the squared distances of each intensity value from the Mean value. This is a measure of
    the spread of the distribution about the mean. By definition, :math:`\textit{variance} = \sigma^2`
    """

    return numpy.nanstd(self.targetVoxelArray, 1) ** 2

  def getUniformityFeatureValue(self):
    r"""
    **19. Uniformity**

    .. math::
      \textit{uniformity} = \displaystyle\sum^{N_g}_{i=1}{p(i)^2}

    Uniformity is a measure of the sum of the squares of each intensity value. This is a measure of the homogeneity of
    the image array, where a greater uniformity implies a greater homogeneity or a smaller range of discrete intensity
    values.

    .. note::
      Defined by IBSI as Intensity Histogram Uniformity.
    """
    p_i = self.coefficients['p_i']
    return numpy.nansum(p_i ** 2, 1)
