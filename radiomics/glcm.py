import numpy
from six.moves import range

from radiomics import base, cMatrices, deprecated


class RadiomicsGLCM(base.RadiomicsFeaturesBase):
  r"""
  A Gray Level Co-occurrence Matrix (GLCM) of size :math:`N_g \times N_g` describes the second-order joint probability
  function of an image region constrained by the mask and is defined as :math:`\textbf{P}(i,j|\delta,\theta)`.
  The :math:`(i,j)^{\text{th}}` element of this matrix represents the number of times the combination of
  levels :math:`i` and :math:`j` occur in two pixels in the image, that are separated by a distance of :math:`\delta`
  pixels along angle :math:`\theta`.
  The distance :math:`\delta` from the center voxel is defined as the distance according to the infinity norm.
  For :math:`\delta=1`, this results in 2 neighbors for each of 13 angles in 3D (26-connectivity) and for
  :math:`\delta=2` a 98-connectivity (49 unique angles).

  Note that pyradiomics by default computes symmetrical GLCM!

  As a two dimensional example, let the following matrix :math:`\textbf{I}` represent a 5x5 image, having 5 discrete
  grey levels:

  .. math::
    \textbf{I} = \begin{bmatrix}
    1 & 2 & 5 & 2 & 3\\
    3 & 2 & 1 & 3 & 1\\
    1 & 3 & 5 & 5 & 2\\
    1 & 1 & 1 & 1 & 2\\
    1 & 2 & 4 & 3 & 5 \end{bmatrix}

  For distance :math:`\delta = 1` (considering pixels with a distance of 1 pixel from each other)
  and angle :math:`\theta=0^\circ` (horizontal plane, i.e. voxels to the left and right of the center voxel),
  the following symmetrical GLCM is obtained:

  .. math::
    \textbf{P} = \begin{bmatrix}
    6 & 4 & 3 & 0 & 0\\
    4 & 0 & 2 & 1 & 3\\
    3 & 2 & 0 & 1 & 2\\
    0 & 1 & 1 & 0 & 0\\
    0 & 3 & 2 & 0 & 2 \end{bmatrix}

  Let:

  - :math:`\epsilon` be an arbitrarily small positive number (:math:`\approx 2.2\times10^{-16}`)
  - :math:`\textbf{P}(i,j)` be the co-occurence matrix for an arbitrary :math:`\delta` and :math:`\theta`
  - :math:`p(i,j)` be the normalized co-occurence matrix and equal to
    :math:`\frac{\textbf{P}(i,j)}{\sum{\textbf{P}(i,j)}}`
  - :math:`N_g` be the number of discrete intensity levels in the image
  - :math:`p_x(i) = \sum^{N_g}_{j=1}{P(i,j)}` be the marginal row probabilities
  - :math:`p_y(j) = \sum^{N_g}_{i=1}{P(i,j)}` be the marginal column probabilities
  - :math:`\mu_x` be the mean gray level intensity of :math:`p_x` and defined as
    :math:`\mu_x = \displaystyle\sum^{N_g}_{i=1}{p_x(i)i}`
  - :math:`\mu_y` be the mean gray level intensity of :math:`p_y` and defined as
    :math:`\mu_y = \displaystyle\sum^{N_g}_{j=1}{p_y(j)j}`
  - :math:`\sigma_x` be the standard deviation of :math:`p_x`
  - :math:`\sigma_y` be the standard deviation of :math:`p_y`
  - :math:`p_{x+y}(k) = \sum^{N_g}_{i=1}\sum^{N_g}_{j=1}{p(i,j)},\text{ where }i+j=k,\text{ and }k=2,3,\dots,2N_g`
  - :math:`p_{x-y}(k) = \sum^{N_g}_{i=1}\sum^{N_g}_{j=1}{p(i,j)},\text{ where }|i-j|=k,\text{ and }k=0,1,\dots,N_g-1`
  - :math:`HX =  -\sum^{N_g}_{i=1}{p_x(i)\log_2\big(p_x(i)+\epsilon\big)}` be the entropy of :math:`p_x`
  - :math:`HY =  -\sum^{N_g}_{j=1}{p_y(j)\log_2\big(p_y(j)+\epsilon\big)}` be the entropy of :math:`p_y`
  - :math:`HXY =  -\sum^{N_g}_{i=1}\sum^{N_g}_{j=1}{p(i,j)\log_2\big(p(i,j)+\epsilon\big)}` be the entropy of
    :math:`p(i,j)`
  - :math:`HXY1 =  -\sum^{N_g}_{i=1}\sum^{N_g}_{j=1}{p(i,j)\log_2\big(p_x(i)p_y(j)+\epsilon\big)}`
  - :math:`HXY2 =  -\sum^{N_g}_{i=1}\sum^{N_g}_{j=1}{p_x(i)p_y(j)\log_2\big(p_x(i)p_y(j)+\epsilon\big)}`

  By default, the value of a feature is calculated on the GLCM for each angle separately, after which the mean of these
  values is returned. If distance weighting is enabled, GLCM matrices are weighted by weighting factor W and
  then summed and normalised. Features are then calculated on the resultant matrix.
  Weighting factor W is calculated for the distance between neighbouring voxels by:

  :math:`W = e^{-\|d\|^2}`, where d is the distance for the associated angle according
  to the norm specified in setting 'weightingNorm'.

  The following class specific settings are possible:

  - distances [[1]]: List of integers. This specifies the distances between the center voxel and the neighbor, for which
    angles should be generated.
  - symmetricalGLCM [True]: boolean, indicates whether co-occurrences should be assessed in two directions per angle,
    which results in a symmetrical matrix, with equal distributions for :math:`i` and :math:`j`. A symmetrical matrix
    corresponds to the GLCM as defined by Haralick et al.
  - weightingNorm [None]: string, indicates which norm should be used when applying distance weighting.
    Enumerated setting, possible values:

    - 'manhattan': first order norm
    - 'euclidean': second order norm
    - 'infinity': infinity norm.
    - 'no_weighting': GLCMs are weighted by factor 1 and summed
    - None: Applies no weighting, mean of values calculated on separate matrices is returned.

    In case of other values, an warning is logged and option 'no_weighting' is used.

  References

  - Haralick, R., Shanmugan, K., Dinstein, I; Textural features for image classification;
    IEEE Transactions on Systems, Man and Cybernetics; 1973(3), p610-621
  - `<https://en.wikipedia.org/wiki/Co-occurrence_matrix>`_
  - `<http://www.fp.ucalgary.ca/mhallbey/the_glcm.htm>`_
  """

  def __init__(self, inputImage, inputMask, **kwargs):
    super(RadiomicsGLCM, self).__init__(inputImage, inputMask, **kwargs)

    self.symmetricalGLCM = kwargs.get('symmetricalGLCM', True)
    self.weightingNorm = kwargs.get('weightingNorm', None)  # manhattan, euclidean, infinity

    self.P_glcm = None
    self.imageArray = self._applyBinning(self.imageArray)

  def _initCalculation(self, voxelCoordinates=None):
    self.P_glcm = self._calculateMatrix(voxelCoordinates)

    self._calculateCoefficients()

    self.logger.debug('GLCM feature class initialized, calculated GLCM with shape %s', self.P_glcm.shape)

  def _calculateMatrix(self, voxelCoordinates=None):
    r"""
    Compute GLCMs for the input image for every direction in 3D.
    Calculated GLCMs are placed in array P_glcm with shape (i/j, a)
    i/j = total gray-level bins for image array,
    a = directions in 3D (generated by imageoperations.generateAngles)
    """
    self.logger.debug('Calculating GLCM matrix in C')

    Ng = self.coefficients['Ng']

    matrix_args = [
      self.imageArray,
      self.maskArray,
      numpy.array(self.settings.get('distances', [1])),
      Ng,
      self.settings.get('force2D', False),
      self.settings.get('force2Ddimension', 0)
    ]
    if self.voxelBased:
      matrix_args += [self.settings.get('kernelRadius', 1), voxelCoordinates]

    P_glcm, angles = cMatrices.calculate_glcm(*matrix_args)

    self.logger.debug('Process calculated matrix')

    # Delete rows and columns that specify gray levels not present in the ROI
    NgVector = range(1, Ng + 1)  # All possible gray values
    GrayLevels = self.coefficients['grayLevels']  # Gray values present in ROI
    emptyGrayLevels = numpy.array(list(set(NgVector) - set(GrayLevels)), dtype=int)  # Gray values NOT present in ROI

    P_glcm = numpy.delete(P_glcm, emptyGrayLevels - 1, 1)
    P_glcm = numpy.delete(P_glcm, emptyGrayLevels - 1, 2)

    # Optionally make GLCMs symmetrical for each angle
    if self.symmetricalGLCM:
      self.logger.debug('Create symmetrical matrix')
      # Transpose and copy GLCM and add it to P_glcm. Numpy.transpose returns a view if possible, use .copy() to ensure
      # a copy of the array is used and not just a view (otherwise erroneous additions can occur)
      P_glcm += numpy.transpose(P_glcm, (0, 2, 1, 3)).copy()

    # Optionally apply a weighting factor
    if self.weightingNorm is not None:
      self.logger.debug('Applying weighting (%s)', self.weightingNorm)
      pixelSpacing = self.inputImage.GetSpacing()[::-1]
      weights = numpy.empty(len(angles))
      for a_idx, a in enumerate(angles):
        if self.weightingNorm == 'infinity':
          weights[a_idx] = numpy.exp(-max(numpy.abs(a) * pixelSpacing) ** 2)
        elif self.weightingNorm == 'euclidean':
          weights[a_idx] = numpy.exp(-numpy.sum((numpy.abs(a) * pixelSpacing) ** 2))  # sqrt ^ 2 = 1
        elif self.weightingNorm == 'manhattan':
          weights[a_idx] = numpy.exp(-numpy.sum(numpy.abs(a) * pixelSpacing) ** 2)
        elif self.weightingNorm == 'no_weighting':
          weights[a_idx] = 1
        else:
          self.logger.warning('weigthing norm "%s" is unknown, W is set to 1', self.weightingNorm)
          weights[a_idx] = 1

      P_glcm = numpy.sum(P_glcm * weights[None, None, None, :], 3, keepdims=True)

    sumP_glcm = numpy.sum(P_glcm, (1, 2))

    # Delete empty angles if no weighting is applied
    if P_glcm.shape[3] > 1:
      emptyAngles = numpy.where(numpy.sum(sumP_glcm, 0) == 0)
      if len(emptyAngles[0]) > 0:  # One or more angles are 'empty'
        self.logger.debug('Deleting %d empty angles:\n%s', len(emptyAngles[0]), angles[emptyAngles])
        P_glcm = numpy.delete(P_glcm, emptyAngles, 3)
        sumP_glcm = numpy.delete(sumP_glcm, emptyAngles, 1)
      else:
        self.logger.debug('No empty angles')

    # Mark empty angles with NaN, allowing them to be ignored in feature calculation
    sumP_glcm[sumP_glcm == 0] = numpy.nan
    # Normalize each glcm
    P_glcm /= sumP_glcm[:, None, None, :]

    return P_glcm

  # check if ivector and jvector can be replaced
  def _calculateCoefficients(self):
    r"""
    Calculate and fill in the coefficients dict.
    """
    self.logger.debug('Calculating GLCM coefficients')

    Ng = self.coefficients['Ng']
    eps = numpy.spacing(1)

    NgVector = self.coefficients['grayLevels'].astype('float')
    # shape = (Ng, Ng)
    i, j = numpy.meshgrid(NgVector, NgVector, indexing='ij', sparse=True)

    # shape = (2*Ng-1)
    kValuesSum = numpy.arange(2, (Ng * 2) + 1, dtype='float')
    # shape = (Ng-1)
    kValuesDiff = numpy.arange(0, Ng, dtype='float')

    # marginal row probabilities #shape = (Nv, Ng, 1, angles)
    px = self.P_glcm.sum(2, keepdims=True)
    # marginal column probabilities #shape = (Nv, 1, Ng, angles)
    py = self.P_glcm.sum(1, keepdims=True)

    # shape = (Nv, 1, 1, angles)
    ux = numpy.sum(i[None, :, :, None] * self.P_glcm, (1, 2), keepdims=True)
    uy = numpy.sum(j[None, :, :, None] * self.P_glcm, (1, 2), keepdims=True)

    # shape = (Nv, 2*Ng-1, angles)
    pxAddy = numpy.array([numpy.sum(self.P_glcm[:, i + j == k, :], 1) for k in kValuesSum]).transpose((1, 0, 2))
    # shape = (Nv, Ng, angles)
    pxSuby = numpy.array([numpy.sum(self.P_glcm[:, numpy.abs(i - j) == k, :], 1) for k in kValuesDiff]).transpose((1, 0, 2))

    # shape = (Nv, angles)
    HXY = (-1) * numpy.sum((self.P_glcm * numpy.log2(self.P_glcm + eps)), (1, 2))

    self.coefficients['eps'] = eps
    self.coefficients['i'] = i
    self.coefficients['j'] = j
    self.coefficients['kValuesSum'] = kValuesSum
    self.coefficients['kValuesDiff'] = kValuesDiff
    self.coefficients['px'] = px
    self.coefficients['py'] = py
    self.coefficients['ux'] = ux
    self.coefficients['uy'] = uy
    self.coefficients['pxAddy'] = pxAddy
    self.coefficients['pxSuby'] = pxSuby
    self.coefficients['HXY'] = HXY

  def getAutocorrelationFeatureValue(self):
    r"""
    **1. Autocorrelation**

    .. math::
      \textit{autocorrelation} = \displaystyle\sum^{N_g}_{i=1}\displaystyle\sum^{N_g}_{j=1}{p(i,j)ij}

    Autocorrelation is a measure of the magnitude of the fineness and coarseness of texture.
    """
    i = self.coefficients['i']
    j = self.coefficients['j']
    ac = numpy.sum(self.P_glcm * (i * j)[None, :, :, None], (1, 2))
    return numpy.nanmean(ac, 1)

  def getJointAverageFeatureValue(self):
    r"""
    **2. Joint Average**

    .. math::
      \textit{joint average} = \mu_x = \displaystyle\sum^{N_g}_{i=1}\displaystyle\sum^{N_g}_{j=1}{p(i,j)i}

    Returns the mean gray level intensity of the :math:`i` distribution.

    .. warning::
      As this formula represents the average of the distribution of :math:`i`, it is independent from the
      distribution of :math:`j`. Therefore, only use this formula if the GLCM is symmetrical, where
      :math:`p_x(i) = p_y(j) \text{, where } i = j`.
    """
    if not self.symmetricalGLCM:
      self.logger.warning('The formula for GLCM - Joint Average assumes that the GLCM is symmetrical, but this is not the case.')
    return self.coefficients['ux'].mean((1, 2, 3))

  def getClusterProminenceFeatureValue(self):
    r"""
    **3. Cluster Prominence**

    .. math::
      \textit{cluster prominence} = \displaystyle\sum^{N_g}_{i=1}\displaystyle\sum^{N_g}_{j=1}
      {\big( i+j-\mu_x-\mu_y\big)^4p(i,j)}

    Cluster Prominence is a measure of the skewness and asymmetry of the GLCM. A higher values implies more asymmetry
    about the mean while a lower value indicates a peak near the mean value and less variation about the mean.
    """
    i = self.coefficients['i']
    j = self.coefficients['j']
    ux = self.coefficients['ux']
    uy = self.coefficients['uy']
    cp = numpy.sum((self.P_glcm * (((i + j)[None, :, :, None] - ux - uy) ** 4)), (1, 2))
    return numpy.nanmean(cp, 1)

  def getClusterShadeFeatureValue(self):
    r"""
    **4. Cluster Shade**

    .. math::
      \textit{cluster shade} = \displaystyle\sum^{N_g}_{i=1}\displaystyle\sum^{N_g}_{j=1}
      {\big(i+j-\mu_x-\mu_y\big)^3p(i,j)}

    Cluster Shade is a measure of the skewness and uniformity of the GLCM.
    A higher cluster shade implies greater asymmetry about the mean.
    """
    i = self.coefficients['i']
    j = self.coefficients['j']
    ux = self.coefficients['ux']
    uy = self.coefficients['uy']
    cs = numpy.sum((self.P_glcm * (((i + j)[None, :, :, None] - ux - uy) ** 3)), (1, 2))
    return numpy.nanmean(cs, 1)

  def getClusterTendencyFeatureValue(self):
    r"""
    **5. Cluster Tendency**

    .. math::
      \textit{cluster tendency} = \displaystyle\sum^{N_g}_{i=1}\displaystyle\sum^{N_g}_{j=1}
      {\big(i+j-\mu_x-\mu_y\big)^2p(i,j)}

    Cluster Tendency is a measure of groupings of voxels with similar gray-level values.
    """
    i = self.coefficients['i']
    j = self.coefficients['j']
    ux = self.coefficients['ux']
    uy = self.coefficients['uy']
    ct = numpy.sum((self.P_glcm * (((i + j)[None, :, :, None] - ux - uy) ** 2)), (1, 2))
    return numpy.nanmean(ct, 1)

  def getContrastFeatureValue(self):
    r"""
    **6. Contrast**

    .. math::
      \textit{contrast} = \displaystyle\sum^{N_g}_{i=1}\displaystyle\sum^{N_g}_{j=1}{(i-j)^2p(i,j)}

    Contrast is a measure of the local intensity variation, favoring values away from the diagonal :math:`(i = j)`. A
    larger value correlates with a greater disparity in intensity values among neighboring voxels.
    """
    i = self.coefficients['i']
    j = self.coefficients['j']
    cont = numpy.sum((self.P_glcm * ((numpy.abs(i - j))[None, :, :, None] ** 2)), (1, 2))
    return numpy.nanmean(cont, 1)

  def getCorrelationFeatureValue(self):
    r"""
    **7. Correlation**

    .. math::
      \textit{correlation} = \frac{\sum^{N_g}_{i=1}\sum^{N_g}_{j=1}{p(i,j)ij-\mu_x\mu_y}}{\sigma_x(i)\sigma_y(j)}

    Correlation is a value between 0 (uncorrelated) and 1 (perfectly correlated) showing the
    linear dependency of gray level values to their respective voxels in the GLCM.

    .. note::
      When there is only 1 discreet gray value in the ROI (flat region), :math:`\sigma_x` and :math:`\sigma_y` will be
      0. In this case, an arbitrary value of 1 is returned instead. This is assessed on a per-angle basis.
    """
    eps = self.coefficients['eps']
    i = self.coefficients['i']
    j = self.coefficients['j']
    ux = self.coefficients['ux']
    uy = self.coefficients['uy']

    # shape = (Nv, 1, 1, angles)
    sigx = numpy.sum(self.P_glcm * ((i[None, :, :, None] - ux) ** 2), (1, 2), keepdims=True) ** 0.5
    # shape = (Nv, 1, 1, angles)
    sigy = numpy.sum(self.P_glcm * ((j[None, :, :, None] - uy) ** 2), (1, 2), keepdims=True) ** 0.5

    corm = numpy.sum(self.P_glcm * (i[None, :, :, None] - ux) * (j[None, :, :, None] - uy), (1, 2), keepdims=True)
    corr = corm / (sigx * sigy + eps)
    corr[sigx * sigy == 0] = 1  # Set elements that would be divided by 0 to 1.
    return numpy.nanmean(corr, (1, 2, 3))

  def getDifferenceAverageFeatureValue(self):
    r"""
    **8. Difference Average**

    .. math::
      \textit{difference average} = \displaystyle\sum^{N_g-1}_{k=0}{kp_{x-y}(k)}

    Difference Average measures the relationship between occurrences of pairs
    with similar intensity values and occurrences of pairs with differing intensity
    values.
    """
    pxSuby = self.coefficients['pxSuby']
    kValuesDiff = self.coefficients['kValuesDiff']
    diffavg = numpy.sum((kValuesDiff[None, :, None] * pxSuby), 1)
    return numpy.nanmean(diffavg, 1)

  def getDifferenceEntropyFeatureValue(self):
    r"""
    **9. Difference Entropy**

    .. math::
      \textit{difference entropy} = \displaystyle\sum^{N_g-1}_{k=0}{p_{x-y}(k)\log_2\big(p_{x-y}(k)+\epsilon\big)}

    Difference Entropy is a measure of the randomness/variability
    in neighborhood intensity value differences.
    """
    pxSuby = self.coefficients['pxSuby']
    eps = self.coefficients['eps']
    difent = (-1) * numpy.sum((pxSuby * numpy.log2(pxSuby + eps)), 1)
    return numpy.nanmean(difent, 1)

  def getDifferenceVarianceFeatureValue(self):
    r"""
    **10. Difference Variance**

    .. math::
      \textit{difference variance} = \displaystyle\sum^{N_g-1}_{k=0}{(k-DA)^2p_{x-y}(k)}

    Difference Variance is a measure of heterogeneity that places higher weights on
    differing intensity level pairs that deviate more from the mean.
    """
    pxSuby = self.coefficients['pxSuby']
    kValuesDiff = self.coefficients['kValuesDiff']
    diffavg = numpy.sum((kValuesDiff[None, :, None] * pxSuby), 1, keepdims=True)
    diffvar = numpy.sum((pxSuby * ((kValuesDiff[None, :, None] - diffavg) ** 2)), 1)
    return numpy.nanmean(diffvar, 1)

  @deprecated
  def getDissimilarityFeatureValue(self):
    r"""
    **DEPRECATED. Dissimilarity**

    .. math::

      \textit{dissimilarity} = \displaystyle\sum^{N_g}_{i=1}\displaystyle\sum^{N_g}_{j=1}{|i-j|p(i,j)}

    .. warning::
      This feature has been deprecated, as it is mathematically equal to Difference Average
      :py:func:`~radiomics.glcm.RadiomicsGLCM.getDifferenceAverageFeatureValue()`.
      See :ref:`here <radiomics-excluded-dissimilarity-label>` for the proof. **Enabling this feature will result in the
      logging of a DeprecationWarning (does not interrupt extraction of other features), no value is calculated for this features**
    """
    raise DeprecationWarning('GLCM - Dissimilarity is mathematically equal to GLCM - Difference Average, '
                             'see http://pyradiomics.readthedocs.io/en/latest/removedfeatures.html for more details')

  def getJointEnergyFeatureValue(self):
    r"""
    **11. Joint Energy**

    .. math::
      \textit{joint energy} = \displaystyle\sum^{N_g}_{i=1}\displaystyle\sum^{N_g}_{j=1}{\big(p(i,j)\big)^2}

    Energy is a measure of homogeneous patterns
    in the image. A greater Energy implies that there are more instances
    of intensity value pairs in the image that neighbor each other at
    higher frequencies.

    .. note::
      Defined by IBSI as Angular Second Moment.
    """
    ene = numpy.sum((self.P_glcm ** 2), (1, 2))
    return numpy.nanmean(ene, 1)

  def getJointEntropyFeatureValue(self):
    r"""
    **12. Joint Entropy**

    .. math::
      \textit{joint entropy} = -\displaystyle\sum^{N_g}_{i=1}\displaystyle\sum^{N_g}_{j=1}
      {p(i,j)\log_2\big(p(i,j)+\epsilon\big)}


    Joint entropy is a measure of the randomness/variability in neighborhood intensity values.

    .. note::
      Defined by IBSI as Joint entropy
    """
    ent = self.coefficients['HXY']
    return numpy.nanmean(ent, 1)

  @deprecated
  def getHomogeneity1FeatureValue(self):
    r"""
    **DEPRECATED. Homogeneity 1**

    .. math::

      \textit{homogeneity 1} = \displaystyle\sum^{N_g}_{i=1}\displaystyle\sum^{N_g}_{j=1}{\frac{p(i,j)}{1+|i-j|}}

    .. warning::
      This feature has been deprecated, as it is mathematically equal to Inverse Difference
      :py:func:`~radiomics.glcm.RadiomicsGLCM.getIdFeatureValue()`. **Enabling this feature will result in the
      logging of a DeprecationWarning (does not interrupt extraction of other features), no value is calculated for this features**
    """
    raise DeprecationWarning('GLCM - Homogeneity 1 is mathematically equal to GLCM - Inverse Difference, '
                             'see documentation of the GLCM feature class (section "Radiomic Features") for more details')

  @deprecated
  def getHomogeneity2FeatureValue(self):
    r"""
    **DEPRECATED. Homogeneity 2**

    .. math::

      \textit{homogeneity 2} = \displaystyle\sum^{N_g}_{i=1}\displaystyle\sum^{N_g}_{j=1}{\frac{p(i,j)}{1+|i-j|^2}}

    .. warning::
      This feature has been deprecated, as it is mathematically equal to Inverse Difference Moment
      :py:func:`~radiomics.glcm.RadiomicsGLCM.getIdmFeatureValue()`. **Enabling this feature will result in the
      logging of a DeprecationWarning (does not interrupt extraction of other features), no value is calculated for this features**
    """
    raise DeprecationWarning('GLCM - Homogeneity 2 is mathematically equal to GLCM - Inverse Difference Moment, '
                             'see documentation of the GLCM feature class (section "Radiomic Features") for more details')

  def getImc1FeatureValue(self):
    r"""
    **13. Informational Measure of Correlation (IMC) 1**

    .. math::

      \textit{IMC 1} = \displaystyle\frac{HXY-HXY1}{\max\{HX,HY\}}

    IMC1 assesses the correlation between the probability distributions of :math:`i` and :math:`j` (quantifying the
    complexity of the texture), using mutual information I(x, y):

    .. math::

      I(i, j) = \sum^{N_g}_{i=1}\sum^{N_g}_{j=1}{p(i,j)\log_2\big(\frac{p(i,j)}{p_x(i)p_y(j)}\big)}

              = \sum^{N_g}_{i=1}\sum^{N_g}_{j=1}{p(i,j)\big(\log_2 (p(i,j)) - \log_2 (p_x(i)p_y(j))\big)}

              = \sum^{N_g}_{i=1}\sum^{N_g}_{j=1}{p(i,j)\log_2 \big(p(i,j)\big)} -
                \sum^{N_g}_{i=1}\sum^{N_g}_{j=1}{p(i,j)\log_2 \big(p_x(i)p_y(j)\big)}

              = -HXY + HXY1

    However, in this formula, the numerator is defined as HXY - HXY1 (i.e. :math:`-I(x, y)`), and is
    therefore :math:`\leq 0`. This reflects how this feature is defined in the original Haralick paper.

    In the case where the distributions are independent, there is no mutual information and the result will therefore be
    0. In the case of uniform distribution with complete dependence, mutual information will be equal to
    :math:`\log_2(N_g)`.

    Finally, :math:`HXY - HXY1` is divided by the maximum of the 2 marginal entropies, where in the latter case of
    complete dependence (not necessarily uniform; low complexity) it will result in :math:`IMC1 = -1`, as
    :math:`HX = HY = I(i, j)`.

    .. note::

      In the case where both HX and HY are 0 (as is the case in a flat region), an arbitrary value of 0 is returned to
      prevent a division by 0. This is done on a per-angle basis (i.e. prior to any averaging).
    """
    HXY = self.coefficients['HXY']

    eps = self.coefficients['eps']
    px = self.coefficients['px']
    py = self.coefficients['py']

    # entropy of px # shape = (Nv, angles)
    HX = (-1) * numpy.sum((px * numpy.log2(px + eps)), (1, 2))
    # entropy of py # shape = (Nv, angles)
    HY = (-1) * numpy.sum((py * numpy.log2(py + eps)), (1, 2))
    # shape = (Nv, angles)
    HXY1 = (-1) * numpy.sum((self.P_glcm * numpy.log2(px * py + eps)), (1, 2))

    div = numpy.fmax(HX, HY)

    imc1 = HXY - HXY1
    imc1[div != 0] /= div[div != 0]
    imc1[div == 0] = 0  # Set elements that would be divided by 0 to 0

    return numpy.nanmean(imc1, 1)

  def getImc2FeatureValue(self):
    r"""
    **14. Informational Measure of Correlation (IMC) 2**

    .. math::

      \textit{IMC 2} = \displaystyle\sqrt{1-e^{-2(HXY2-HXY)}}

    IMC2 also assesses the correlation between the probability distributions of :math:`i` and :math:`j` (quantifying the
    complexity of the texture). Of interest is to note that :math:`HXY1 = HXY2` and that :math:`HXY2 - HXY \geq 0`
    represents the mutual information of the 2 distributions. Therefore, the range of IMC2 = [0, 1), with 0 representing
    the case of 2 independent distributions (no mutual information) and the maximum value representing the case of 2
    fully dependent and uniform distributions (maximal mutual information, equal to :math:`\log_2(N_g)`). In this latter
    case, the maximum value is then equal to :math:`\displaystyle\sqrt{1-e^{-2\log_2(N_g)}}`, approaching 1.

    .. note::

      Due to machine precision errors, it is possble that HXY > HXY2, which would result in returning complex numbers.
      In these cases, a value of 0 is returned for IMC2. This is done on a per-angle basis (i.e. prior to any
      averaging).
    """
    HXY = self.coefficients['HXY']

    eps = self.coefficients['eps']
    px = self.coefficients['px']
    py = self.coefficients['py']

    # shape = (Nv, angles)
    HXY2 = (-1) * numpy.sum(((px * py) * numpy.log2(px * py + eps)), (1, 2))

    imc2 = (1 - numpy.e ** (-2 * (HXY2 - HXY))) ** 0.5
    imc2[HXY2 == HXY] = 0

    return numpy.nanmean(imc2, 1)

  def getIdmFeatureValue(self):
    r"""
    **15. Inverse Difference Moment (IDM)**

    .. math::
      \textit{IDM} = \displaystyle\sum^{N_g-1}_{k=0}{\frac{p_{x-y}(k)}{1+k^2}}

    IDM (a.k.a Homogeneity 2) is a measure of the local
    homogeneity of an image. IDM weights are the inverse of the Contrast
    weights (decreasing exponentially from the diagonal i=j in the GLCM).
    """
    pxSuby = self.coefficients['pxSuby']
    kValuesDiff = self.coefficients['kValuesDiff']
    idm = numpy.sum(pxSuby / (1 + (kValuesDiff[None, :, None] ** 2)), 1)
    return numpy.nanmean(idm, 1)

  def getMCCFeatureValue(self):
    r"""
    **16. Maximal Correlation Coefficient (MCC)**

    .. math::
      \textit{MCC} = \sqrt{\text{second largest eigenvalue of Q}}

      Q(i, j) = \displaystyle\sum^{N_g}_{k=0}{\frac{p(i,k)p(j, k)}{p_x(i)p_y(k)}}

    The Maximal Correlation Coefficient is a measure of complexity of the texture and :math:`0 \leq MCC \leq 1`.

    In case of a flat region, each GLCM matrix has shape (1, 1), resulting in just 1 eigenvalue. In this case, an
    arbitrary value of 1 is returned.
    """
    px = self.coefficients['px']
    py = self.coefficients['py']
    eps = self.coefficients['eps']

    # Calculate Q (shape (i, i, d)). To prevent division by 0, add epsilon (such a division can occur when in a ROI
    # along a certain angle, voxels with gray level i do not have neighbors
    Q = ((self.P_glcm[:, :, None, 0, :] * self.P_glcm[:, None, :, 0, :]) /  # slice: v, i, j, k, d
         (px[:, :, None, 0, :] * py[:, None, :, 0, :] + eps))  # sum over k (4th axis --> index 3)

    for gl in range(1, self.P_glcm.shape[1]):
      Q += ((self.P_glcm[:, :, None, gl, :] * self.P_glcm[:, None, :, gl, :]) /  # slice: v, i, j, k, d
            (px[:, :, None, 0, :] * py[:, None, :, gl, :] + eps))  # sum over k (4th axis --> index 3)

    # calculation of eigenvalues if performed on last 2 dimensions, therefore, move the angles dimension (d) forward
    Q_eigenValue = numpy.linalg.eigvals(Q.transpose((0, 3, 1, 2)))
    Q_eigenValue.sort()  # sorts along last axis --> eigenvalues, low to high

    if Q_eigenValue.shape[2] < 2:
      return 1  # flat region

    MCC = numpy.sqrt(Q_eigenValue[:, :, -2])  # 2nd highest eigenvalue

    return numpy.nanmean(MCC, 1).real

  def getIdmnFeatureValue(self):
    r"""
    **17. Inverse Difference Moment Normalized (IDMN)**

    .. math::
      \textit{IDMN} = \displaystyle\sum^{N_g-1}_{k=0}{ \frac{p_{x-y}(k)}{1+\left(\frac{k^2}{N_g^2}\right)} }

    IDMN (inverse difference moment normalized)  is a measure of the local
    homogeneity of an image. IDMN weights are the inverse of the Contrast
    weights (decreasing exponentially from the diagonal :math:`i=j` in the GLCM).
    Unlike Homogeneity2, IDMN normalizes the square of the difference between
    neighboring intensity values by dividing over the square of the total
    number of discrete intensity values.
    """
    pxSuby = self.coefficients['pxSuby']
    kValuesDiff = self.coefficients['kValuesDiff']
    Ng = self.coefficients['Ng']
    idmn = numpy.sum(pxSuby / (1 + ((kValuesDiff[None, :, None] ** 2) / (Ng ** 2))), 1)
    return numpy.nanmean(idmn, 1)

  def getIdFeatureValue(self):
    r"""
    **18. Inverse Difference (ID)**

    .. math::
      \textit{ID} = \displaystyle\sum^{N_g-1}_{k=0}{\frac{p_{x-y}(k)}{1+k}}

    ID (a.k.a. Homogeneity 1) is another measure of the local homogeneity of an image.
    With more uniform gray levels, the denominator will remain low, resulting in a higher overall value.
    """
    pxSuby = self.coefficients['pxSuby']
    kValuesDiff = self.coefficients['kValuesDiff']
    invDiff = numpy.sum(pxSuby / (1 + kValuesDiff[None, :, None]), 1)
    return numpy.nanmean(invDiff, 1)

  def getIdnFeatureValue(self):
    r"""
    **19. Inverse Difference Normalized (IDN)**

    .. math::
      \textit{IDN} = \displaystyle\sum^{N_g-1}_{k=0}{ \frac{p_{x-y}(k)}{1+\left(\frac{k}{N_g}\right)} }

    IDN (inverse difference normalized) is another measure of the local
    homogeneity of an image. Unlike Homogeneity1, IDN normalizes the difference
    between the neighboring intensity values by dividing over the total number
    of discrete intensity values.
    """
    pxSuby = self.coefficients['pxSuby']
    kValuesDiff = self.coefficients['kValuesDiff']
    Ng = self.coefficients['Ng']
    idn = numpy.sum(pxSuby / (1 + (kValuesDiff[None, :, None] / Ng)), 1)
    return numpy.nanmean(idn, 1)

  def getInverseVarianceFeatureValue(self):
    r"""
    **20. Inverse Variance**

    .. math::
      \textit{inverse variance} = \displaystyle\sum^{N_g-1}_{k=1}{\frac{p_{x-y}(k)}{k^2}}

    Note that :math:`k=0` is skipped, as this would result in a division by 0.
    """
    pxSuby = self.coefficients['pxSuby']
    kValuesDiff = self.coefficients['kValuesDiff']
    inv = numpy.sum(pxSuby[:, 1:, :] / kValuesDiff[None, 1:, None] ** 2, 1)  # Skip k = 0 (division by 0)
    return numpy.nanmean(inv, 1)

  def getMaximumProbabilityFeatureValue(self):
    r"""
    **21. Maximum Probability**

    .. math::

      \textit{maximum probability} = \max\big(p(i,j)\big)

    Maximum Probability is occurrences of the most predominant pair of
    neighboring intensity values.

    .. note::
      Defined by IBSI as Joint maximum
    """
    maxprob = numpy.amax(self.P_glcm, (1, 2))
    return numpy.nanmean(maxprob, 1)

  def getSumAverageFeatureValue(self):
    r"""
    **22. Sum Average**

    .. math::

      \textit{sum average} = \displaystyle\sum^{2N_g}_{k=2}{p_{x+y}(k)k}

    Sum Average measures the relationship between occurrences of pairs
    with lower intensity values and occurrences of pairs with higher intensity
    values.

    .. warning::
      When GLCM is symmetrical, :math:`\mu_x = \mu_y`, and therefore :math:`\text{Sum Average} = \mu_x + \mu_y =
      2 \mu_x = 2 * Joint Average`. See formulas (4.), (5.) and (6.) defined
      :ref:`here <radiomics-excluded-sumvariance-label>` for the proof that :math:`\text{Sum Average} = \mu_x + \mu_y`.
      In the default parameter files provided in the ``examples/exampleSettings``, this feature has been disabled.
    """
    # warn the user if the GLCM is symmetrical and this feature is calculated (as it is then linearly correlated to Joint Average)
    if self.symmetricalGLCM:
      self.logger.warning('GLCM is symmetrical, therefore Sum Average = 2 * Joint Average, only 1 needs to be calculated')

    pxAddy = self.coefficients['pxAddy']
    kValuesSum = self.coefficients['kValuesSum']
    sumavg = numpy.sum((kValuesSum[None, :, None] * pxAddy), 1)
    return numpy.nanmean(sumavg, 1)

  @deprecated
  def getSumVarianceFeatureValue(self):
    r"""
    **DEPRECATED. Sum Variance**

    .. math::
      \textit{sum variance} = \displaystyle\sum^{2N_g}_{k=2}{(k-SA)^2p_{x+y}(k)}

    .. warning::
      This feature has been deprecated, as it is mathematically equal to Cluster Tendency
      :py:func:`~radiomics.glcm.RadiomicsGLCM.getClusterTendencyFeatureValue()`.
      See :ref:`here <radiomics-excluded-sumvariance-label>` for the proof. **Enabling this feature will result in the
      logging of a DeprecationWarning (does not interrupt extraction of other features), no value is calculated for this features**
    """
    raise DeprecationWarning('GLCM - Sum Variance is mathematically equal to GLCM - Cluster Tendency, '
                             'see http://pyradiomics.readthedocs.io/en/latest/removedfeatures.html for more details')

  def getSumEntropyFeatureValue(self):
    r"""
    **23. Sum Entropy**

    .. math::

      \textit{sum entropy} = \displaystyle\sum^{2N_g}_{k=2}{p_{x+y}(k)\log_2\big(p_{x+y}(k)+\epsilon\big)}

    Sum Entropy is a sum of neighborhood intensity value differences.
    """
    pxAddy = self.coefficients['pxAddy']
    eps = self.coefficients['eps']
    sumentr = (-1) * numpy.sum((pxAddy * numpy.log2(pxAddy + eps)), 1)
    return numpy.nanmean(sumentr, 1)

  def getSumSquaresFeatureValue(self):
    r"""
    **24. Sum of Squares**

    .. math::

      \textit{sum squares} = \displaystyle\sum^{N_g}_{i=1}\displaystyle\sum^{N_g}_{j=1}{(i-\mu_x)^2p(i,j)}

    Sum of Squares or Variance is a measure in the distribution of neigboring intensity level pairs
    about the mean intensity level in the GLCM.

    .. warning::

      This formula represents the variance of the distribution of :math:`i` and is independent from the distribution
      of :math:`j`. Therefore, only use this formula if the GLCM is symmetrical, where
      :math:`p_x(i) = p_y(j) \text{, where } i = j`

    .. note::
      Defined by IBSI as Joint Variance
    """
    if not self.symmetricalGLCM:
      self.logger.warning('The formula for GLCM - Sum of Squares assumes that the GLCM is symmetrical, but this is not the case.')
    i = self.coefficients['i']
    ux = self.coefficients['ux']
    # Also known as Variance
    ss = numpy.sum((self.P_glcm * ((i[None, :, :, None] - ux) ** 2)), (1, 2))
    return numpy.nanmean(ss, 1)
