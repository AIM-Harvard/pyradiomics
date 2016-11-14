import numpy
from six.moves import range
from tqdm import trange

import radiomics
from radiomics import base, imageoperations

class RadiomicsGLCM(base.RadiomicsFeaturesBase):
  r"""
  A Gray Level Co-occurrence Matrix (GLCM) of size :math:`N_g \times N_g` describes the second-order joint probability function of an image region
  constrained by the mask and is defined as :math:`\textbf{P}(i,j|\delta,\alpha)`.
  The :math:`(i,j)`\ :sup:`th` element of this matrix represents the number of times the combination of
  levels :math:`i` and :math:`j` occur in two pixels in the image,
  that are separated by a distance of :math:`\delta` pixels in direction :math:`\alpha`, and :math:`N_g`
  is the number of discrete gray level intensities.
  The distance :math:`\delta` from the center voxel is defined as the distance according to the infinity norm.
  For :math:`\delta=1`, this assumes 26-connectivity in 3D and for :math:`\delta=2` a 98-connectivity.

  Note that pyradiomics by default computes symmetrical GLCM!

  As a two dimensional example, let the following matrix :math:`\textbf{I}` represent a 5x5 image, having 5 discrete grey levels:

  .. math::
    \textbf{I} = \begin{bmatrix}
    1 & 2 & 5 & 2 & 3\\
    3 & 2 & 1 & 3 & 1\\
    1 & 3 & 5 & 5 & 2\\
    1 & 1 & 1 & 1 & 2\\
    1 & 2 & 4 & 3 & 5 \end{bmatrix}

  For distance :math:`\delta = 1` (considering pixels with a distance of 1 pixel from each other)
  in directions :math:`\alpha=0^\circ` and opposite :math:`\alpha=180^\circ`
  (i.e., to the left and right from the pixel with the given value), the following symmetrical GLCM is obtained:

  .. math::
    \textbf{P} = \begin{bmatrix}
    6 & 4 & 3 & 0 & 0\\
    4 & 0 & 2 & 1 & 3\\
    3 & 2 & 0 & 1 & 2\\
    0 & 1 & 1 & 0 & 0\\
    0 & 3 & 2 & 0 & 2 \end{bmatrix}

  Let:

  :math:`\epsilon` be an arbitrarily small positive number (:math:`\approx 2.2\times10^{-16}`)

  :math:`\textbf{P}(i,j)` be the co-occurence matrix for an arbitrary :math:`\delta` and :math:`\alpha`

  :math:`p(i,j)` be the normalized co-occurence matrix and equal to :math:`\frac{\textbf{P}(i,j)}{\sum{\textbf{P}(i,j)}}`

  :math:`N_g` be the number of discrete intensity levels in the image

  :math:`p_x(i) = \sum^{N_g}_{j=1}{P(i,j)}` be the marginal row probabilities

  :math:`p_y(j) = \sum^{N_g}_{i=1}{P(i,j)}` be the marginal column probabilities

  :math:`\mu_x` be the mean gray level intensity of :math:`p_x` and defined as
  :math:`\mu_x = \displaystyle\sum^{N_g}_{i=1}\displaystyle\sum^{N_g}_{j=1}{p(i,j)i}`

  :math:`\mu_y` be the mean gray level intensity of :math:`p_y` and defined as
  :math:`\mu_y = \displaystyle\sum^{N_g}_{i=1}\displaystyle\sum^{N_g}_{j=1}{p(i,j)j}`

  :math:`\sigma_x` be the standard deviation of :math:`p_x`

  :math:`\sigma_y` be the standard deviation of :math:`p_y`

  :math:`p_{x+y}(k) = \sum^{N_g}_{i=1}\sum^{N_g}_{j=1}{p(i,j)},\text{ where }i+j=k,\text{ and }k=2,3,\dots,2N_g`

  :math:`p_{x-y}(k) = \sum^{N_g}_{i=1}\sum^{N_g}_{j=1}{p(i,j)},\text{ where }|i-j|=k,\text{ and }k=0,1,\dots,N_g-1`

  :math:`HX =  -\sum^{N_g}_{i=1}{p_x(i)\log_2\big(p_x(i)+\epsilon\big)}` be the entropy of :math:`p_x`

  :math:`HY =  -\sum^{N_g}_{j=1}{p_y(j)\log_2\big(p_y(j)+\epsilon\big)}` be the entropy of :math:`p_y`

  :math:`HXY =  -\sum^{N_g}_{i=1}\sum^{N_g}_{j=1}{p(i,j)\log_2\big(p(i,j)+\epsilon\big)}` be the entropy of :math:`p(i,j)`

  :math:`HXY1 =  -\sum^{N_g}_{i=1}\sum^{N_g}_{j=1}{p(i,j)\log_2\big(p_x(i)p_y(j)+\epsilon\big)}`

  :math:`HXY2 =  -\sum^{N_g}_{i=1}\sum^{N_g}_{j=1}{p_x(i)p_y(j)\log_2\big(p_x(i)p_y(j)+\epsilon\big)}`

  By default, the value of a feature is calculated on the GLCM for each angle separately, after which the mean of these
  values is returned. If distance weighting is enabled, GLCM matrices are weighted by weighting factor W and
  then summed and normalised. Features are then calculated on the resultant matrix.
  Weighting factor W is calculated for the distance between neighbouring voxels by:

  :math:`W = e^{-\|d\|^2}`, where d is the distance for the associated angle according
  to the norm specified in setting 'weightingNorm'.

  The following class specific settings are possible:

  - symmetricalGLCM [True]: boolean, indicates whether co-occurrences should be assessed in two directions per angle,
    which results in a symmetrical matrix, with equal distributions for :math:`i` and :math:`j`.
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

  - https://en.wikipedia.org/wiki/Co-occurrence_matrix

  - http://www.fp.ucalgary.ca/mhallbey/the_glcm.htm
  """

  def __init__(self, inputImage, inputMask, **kwargs):
    super(RadiomicsGLCM, self).__init__(inputImage, inputMask, **kwargs)

    self.symmetricalGLCM = kwargs.get('symmetricalGLCM', True)
    self.weightingNorm = kwargs.get('weightingNorm', None)  # manhattan, euclidean, infinity

    self.coefficients = {}
    self.P_glcm = None

    # binning
    self.matrix, self.histogram = imageoperations.binImage(self.binWidth, self.matrix, self.matrixCoordinates)
    self.coefficients['Ng'] = self.histogram[1].shape[0] - 1

    if radiomics.cMatsEnabled:
      self.P_glcm = self._calculateCGLCM()
    else:
      self.P_glcm = self._calculateGLCM()
    self._calculateCoefficients()

  def _calculateGLCM(self):
    r"""
    Compute GLCMs for the input image for every direction in 3D.
    Calculated GLCMs are placed in array P_glcm with shape (i/j, a)
    i/j = total gray-level bins for image array,
    a = directions in 3D (generated by imageoperations.generateAngles)
    """
    Ng = self.coefficients['Ng']

    # Exclude voxels outside segmentation, due to binning, no negative values will be encountered inside the mask
    self.matrix[self.maskArray == 0] = -1

    size = numpy.max(self.matrixCoordinates, 1) - numpy.min(self.matrixCoordinates, 1) + 1
    angles = imageoperations.generateAngles(size)

    P_glcm = numpy.zeros((Ng, Ng, int(angles.shape[0])), dtype='float64')

    if self.verbose: bar = trange(Ng, desc='calculate GLCM')

    # iterate over gray levels for center voxel
    for i in xrange(1, Ng + 1):
      # give some progress
      if self.verbose: bar.update()

      # get the indices to all voxels which have the current gray level i
      i_indices = numpy.where(self.matrix == i)

      # iterate over gray levels for neighbouring voxel
      for j in range(1, Ng + 1):
        # get the indices to all voxels which have the current gray level j
        j_indices = set(zip(*numpy.where(self.matrix == j)))

        for a_idx, a in enumerate(angles):
          # get the corresponding indices of the neighbours for angle a
          neighbour_indices = set(zip(*(i_indices + a[:, None])))

          # The following intersection yields the indices to voxels with gray level j
          # that are also a neighbour of a voxel with gray level i for angle a.
          # The number of indices is then equal to the total number of pairs with gray level i and j for angle a
          count = len(neighbour_indices.intersection(j_indices))
          P_glcm[i - 1, j - 1, a_idx] = count
    if self.verbose: bar.close()

    # Optionally make GLCMs symmetrical for each angle
    if self.symmetricalGLCM:
      P_glcm += numpy.transpose(P_glcm, (1, 0, 2))

    # Optionally apply a weighting factor
    if self.weightingNorm is not None:
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

      P_glcm = numpy.sum(P_glcm * weights[None, None, :], 2, keepdims=True)

    sumP_glcm = numpy.sum(P_glcm, (0, 1), keepdims=True)  # , keepdims=True)

    # Delete empty angles if no weighting is applied
    if P_glcm.shape[2] > 1:
      P_glcm = numpy.delete(P_glcm, numpy.where(sumP_glcm == 0), 2)
      sumP_glcm = numpy.delete(sumP_glcm, numpy.where(sumP_glcm == 0), 2)

    # Normalize each glcm
    return P_glcm / sumP_glcm

  def _calculateCGLCM(self):
    size = numpy.max(self.matrixCoordinates, 1) - numpy.min(self.matrixCoordinates, 1) + 1
    angles = imageoperations.generateAngles(size)
    Ng = self.coefficients['Ng']

    P_glcm = radiomics.cMatrices.calculate_glcm(self.matrix, self.maskArray, angles, Ng)

    # Optionally make GLCMs symmetrical for each angle
    if self.symmetricalGLCM:
      P_glcm += numpy.transpose(P_glcm, (1, 0, 2))

    # Optionally apply a weighting factor
    if self.weightingNorm is not None:
      pixelSpacing = self.inputImage.GetSpacing()[::-1]
      weights = numpy.empty(len(angles))
      for a_idx, a in enumerate(angles):
        if self.weightingNorm == 'infinity':
          weights[a_idx] = numpy.exp(-max(numpy.abs(a) * pixelSpacing) ** 2)
        elif self.weightingNorm == 'euclidean':
          weights[a_idx] = numpy.exp(-numpy.sum((numpy.abs(a) * pixelSpacing) ** 2))  # sqrt ^ 2 = 1
        elif self.weightingNorm == 'manhattan':
          weights[a_idx] = numpy.exp(-numpy.sum(numpy.abs(a) * pixelSpacing) ** 2)
        else:
          self.logger.warning('weigthing norm "%s" is unknown, W is set to 1', self.weightingNorm)
          weights[a_idx] = 1

      P_glcm = numpy.sum(P_glcm * weights[None, None, :], 2, keepdims=True)

    sumP_glcm = numpy.sum(P_glcm, (0, 1), keepdims=True)  # , keepdims=True)

    # Delete empty angles if no weighting is applied
    if P_glcm.shape[2] > 1:
      P_glcm = numpy.delete(P_glcm, numpy.where(sumP_glcm == 0), 2)
      sumP_glcm = numpy.delete(sumP_glcm, numpy.where(sumP_glcm == 0), 2)

    # Normalize each glcm
    return P_glcm / sumP_glcm

  # check if ivector and jvector can be replaced
  def _calculateCoefficients(self):
    r"""
    Calculate and fill in the coefficients dict.
    """
    Ng = self.coefficients['Ng']
    eps = numpy.spacing(1)

    NgVector = numpy.arange(1, self.P_glcm.shape[0] + 1, dtype='float64')
    # shape = (Ng, Ng)
    i, j = numpy.meshgrid(NgVector, NgVector, indexing='ij')

    # shape = (2*Ng-1)
    kValuesSum = numpy.arange(2, (Ng * 2) + 1)
    # shape = (Ng-1)
    kValuesDiff = numpy.arange(0, Ng)

    # marginal row probabilities #shape = (Ng, 1, angles)
    px = self.P_glcm.sum(1, keepdims=True)
    # marginal column probabilities #shape = (1, Ng, angles)
    py = self.P_glcm.sum(0, keepdims=True)

    # shape = (1, 1, angles)
    ux = numpy.sum(i[:, :, None] * self.P_glcm, (0, 1), keepdims=True)
    uy = numpy.sum(j[:, :, None] * self.P_glcm, (0, 1), keepdims=True)

    # shape = (1, 1, angles)
    sigx = numpy.sum(self.P_glcm * ((i[:, :, None] - ux) ** 2), (0, 1), keepdims=True) ** 0.5
    # shape = (1, 1, angles)
    sigy = numpy.sum(self.P_glcm * ((j[:, :, None] - uy) ** 2), (0, 1), keepdims=True) ** 0.5

    # shape = (2*Ng-1, angles)
    pxAddy = numpy.array([numpy.sum(self.P_glcm[i + j == k], 0) for k in kValuesSum])
    # shape = (Ng, angles)
    pxSuby = numpy.array([numpy.sum(self.P_glcm[numpy.abs(i - j) == k], 0) for k in kValuesDiff])

    # entropy of px # shape = (angles)
    HX = (-1) * numpy.sum((px * numpy.log2(px + eps)), (0, 1))
    # entropy of py # shape = (angles)
    HY = (-1) * numpy.sum((py * numpy.log2(py + eps)), (0, 1))
    # shape = (angles)
    HXY = (-1) * numpy.sum((self.P_glcm * numpy.log2(self.P_glcm + eps)), (0, 1))

    # shape = (angles)
    HXY1 = (-1) * numpy.sum((self.P_glcm * numpy.log2(px * py + eps)), (0, 1))
    # shape = (angles)
    HXY2 = (-1) * numpy.sum(((px * py) * numpy.log2(px * py + eps)), (0, 1))

    self.coefficients['eps'] = eps
    self.coefficients['i'] = i
    self.coefficients['j'] = j
    self.coefficients['kValuesSum'] = kValuesSum
    self.coefficients['kValuesDiff'] = kValuesDiff
    self.coefficients['px'] = px
    self.coefficients['py'] = py
    self.coefficients['ux'] = ux
    self.coefficients['uy'] = uy
    self.coefficients['sigx'] = sigx
    self.coefficients['sigy'] = sigy
    self.coefficients['pxAddy'] = pxAddy
    self.coefficients['pxSuby'] = pxSuby
    self.coefficients['HX'] = HX
    self.coefficients['HY'] = HY
    self.coefficients['HXY'] = HXY
    self.coefficients['HXY1'] = HXY1
    self.coefficients['HXY2'] = HXY2

  def getAutocorrelationFeatureValue(self):
    r"""
    Calculate and return the mean Autocorrelation.

    :math:`autocorrelation = \displaystyle\sum^{N_g}_{i=1}\displaystyle\sum^{N_g}_{j=1}{p(i,j)ij}`

    Autocorrelation is a measure of the magnitude of the
    fineness and coarseness of texture.
    """
    i = self.coefficients['i']
    j = self.coefficients['j']
    ac = numpy.sum(self.P_glcm * (i * j)[:, :, None], (0, 1))
    return (ac.mean())

  def getAverageIntensityFeatureValue(self):
    r"""
    Return the mean gray level intensity of the :math:`i` distribution.

    :math:`\mu_x = \displaystyle\sum^{N_g}_{i=1}\displaystyle\sum^{N_g}_{j=1}{p(i,j)i}`

    N.B. As this formula represents the average of the distribution of :math:`i`, it is independent from the
    distribution of :math:`j`. Therefore, only use this formula if the GLCM is symmetrical, where both distrubutions
    are equal.
    """

    return self.coefficients['ux'].mean()

  def getClusterProminenceFeatureValue(self):
    r"""
    Using coefficients :math:`\mu_x` and :math:`\mu_y`, calculate and return the mean Cluster Prominence.

    :math:`cluster\ prominence = \displaystyle\sum^{N_g}_{i=1}\displaystyle\sum^{N_g}_{j=1}{\big( i+j-\mu_x(i)-\mu_y(j)\big)^4p(i,j)}`

    Cluster Prominence is a measure of the skewness and asymmetry of the GLCM.
    A higher values implies more asymmetry about the mean while a lower value
    indicates a peak near the mean value and less variation about the mean.
    """
    i = self.coefficients['i']
    j = self.coefficients['j']
    ux = self.coefficients['ux']
    uy = self.coefficients['uy']
    cp = numpy.sum((self.P_glcm * (((i + j)[:, :, None] - ux - uy) ** 4)), (0, 1))
    return (cp.mean())

  def getClusterShadeFeatureValue(self):
    r"""
    Using coefficients :math:`\mu_x` and :math:`\mu_y`, calculate and return the mean Cluster Shade.

    :math:`cluster\ shade = \displaystyle\sum^{N_g}_{i=1}\displaystyle\sum^{N_g}_{j=1}{\big(i+j-\mu_x(i)-\mu_y(j)\big)^3p(i,j)}`

    Cluster Shade is a measure of the skewness and uniformity of the GLCM.
    A higher cluster shade implies greater asymmetry about the mean.
    """
    i = self.coefficients['i']
    j = self.coefficients['j']
    ux = self.coefficients['ux']
    uy = self.coefficients['uy']
    cs = numpy.sum((self.P_glcm * (((i + j)[:, :, None] - ux - uy) ** 3)), (0, 1))
    return (cs.mean())

  def getClusterTendencyFeatureValue(self):
    r"""
    Using coefficients :math:`\mu_x` and :math:`\mu_y`, calculate and return the mean Cluster Tendency.

    :math:`cluster\ prominence = \displaystyle\sum^{N_g}_{i=1}\displaystyle\sum^{N_g}_{j=1}{\big(i+j-\mu_x(i)-\mu_y(j)\big)^2p(i,j)}`

    Cluster Tendency is a measure of groupings of voxels with similar gray-level values.
    """
    i = self.coefficients['i']
    j = self.coefficients['j']
    ux = self.coefficients['ux']
    uy = self.coefficients['uy']
    ct = numpy.sum((self.P_glcm * (((i + j)[:, :, None] - ux - uy) ** 2)), (0, 1))
    return (ct.mean())

  def getContrastFeatureValue(self):
    r"""
    Using the squared difference between gray values of neighbouring paris, calculate and return the mean Contrast.

    :math:`contrast = \displaystyle\sum^{N_g}_{i=1}\displaystyle\sum^{N_g}_{j=1}{(i-j)^2p(i,j)}`

    Contrast is a measure of the local intensity variation, favoring :math:`P(i,j)`
    values away from the diagonal :math:`(i = j)`. A larger value correlates with
    a greater disparity in intensity values among neighboring voxels.
    """
    i = self.coefficients['i']
    j = self.coefficients['j']
    cont = numpy.sum((self.P_glcm * ((numpy.abs(i - j))[:, :, None] ** 2)), (0, 1))
    return (cont.mean())

  def getCorrelationFeatureValue(self):
    r"""
    Using coefficients :math:`\mu_x`, :math:`\mu_y`, :math:`\sigma_x` and :math:`\sigma_y`, calculate and return the
    mean Correlation.

    :math:`correlation = \frac{\sum^{N_g}_{i=1}\sum^{N_g}_{j=1}{p(i,j)ij-\mu_x(i)\mu_y(j)}}{\sigma_x(i)\sigma_y(j)}`

    Correlation is a value between 0 (uncorrelated) and 1 (perfectly correlated) showing the
    linear dependency of gray level values to their respective voxels in the GLCM.
    """
    i = self.coefficients['i']
    j = self.coefficients['j']
    ux = self.coefficients['ux']
    uy = self.coefficients['uy']
    sigx = self.coefficients['sigx']
    sigy = self.coefficients['sigy']

    try:
      corm = numpy.sum(self.P_glcm * (i[:, :, None] - ux) * (j[:, :, None] - uy), (0, 1), keepdims=True)
      corr = corm / (sigx * sigy)
      return (corr.mean())
    except ZeroDivisionError:
      return numpy.core.nan

  def getDifferenceAverageFeatureValue(self):
    r"""
    Using coefficient :math:`p_{x-y}`, calculate and return the mean Difference Average.

    :math:`Difference\ average = \displaystyle\sum^{N_g-1}_{k=0}{k\textbf{P}_{x-y}(k)}`

    Difference Average measures the relationship between occurrences of pairs
    with similar intensity values and occurrences of pairs with differing intensity
    values.
    """
    pxSuby = self.coefficients['pxSuby']
    kValuesDiff = self.coefficients['kValuesDiff']
    diffavg = numpy.sum((kValuesDiff[:, None] * pxSuby), 0)
    return (diffavg.mean())

  def getDifferenceEntropyFeatureValue(self):
    r"""
    Using coefficient :math:`p_{x-y}`, calculate and return the mean Difference Entropy.

    :math:`difference\ entropy = \displaystyle\sum^{N_g-1}_{k=0}{p_{x-y}(k)\log_2\big(p_{x-y}(k)+\epsilon\big)}`

    Difference Entropy is a measure of the randomness/variability
    in neighborhood intensity value differences.
    """
    pxSuby = self.coefficients['pxSuby']
    eps = self.coefficients['eps']
    difent = (-1) * numpy.sum((pxSuby * numpy.log2(pxSuby + eps)), 0)
    return (difent.mean())

  def getDifferenceVarianceFeatureValue(self):
    r"""
    Using coefficients :math:`p_{x-y}` and DifferenceAverage (DA) calculate and return the mean Difference Variance.

    :math:`Difference\ variance = \displaystyle\sum^{N_g-1}_{k=0}{(1-DA)^2\textbf{P}_{x-y}(k)}`

    Difference Variance is a measure of heterogeneity that places higher weights on
    differing intensity level pairs that deviate more from the mean.
    """
    pxSuby = self.coefficients['pxSuby']
    kValuesDiff = self.coefficients['kValuesDiff']
    diffavg = numpy.sum((kValuesDiff[:, None] * pxSuby), 0, keepdims=True)
    diffvar = numpy.sum((pxSuby * ((kValuesDiff[:, None] - diffavg) ** 2)), 0)
    return (diffvar.mean())

  def getDissimilarityFeatureValue(self):
    r"""
    Calculate and return the mean Dissimilarity.

    :math:`dissimilarity = \displaystyle\sum^{N_g}_{i=1}\displaystyle\sum^{N_g}_{j=1}{|i-j|p(i,j)}`

    Dissimilarity is a measure of local intensity variation defined as the mean absolute difference between the
    neighbouring pairs. A larger value correlates with a greater disparity in intensity values
    among neighboring voxels.
    """
    i = self.coefficients['i']
    j = self.coefficients['j']
    dis = numpy.sum((self.P_glcm * (numpy.abs(i - j))[:, :, None]), (0, 1))
    return (dis.mean())

  def getEnergyFeatureValue(self):
    r"""
    Calculate and return the mean Energy.

    :math:`energy = \displaystyle\sum^{N_g}_{i=1}\displaystyle\sum^{N_g}_{j=1}{\big(p(i,j)\big)^2}`

    Energy (or Angular Second Moment)is a measure of homogeneous patterns
    in the image. A greater Energy implies that there are more instances
    of intensity value pairs in the image that neighbor each other at
    higher frequencies.
    """
    ene = numpy.sum((self.P_glcm ** 2), (0, 1))
    return (ene.mean())

  def getEntropyFeatureValue(self):
    r"""
    Calculate and return the mean Entropy.

    :math:`entropy = -\displaystyle\sum^{N_g}_{i=1}\displaystyle\sum^{N_g}_{j=1}{p(i,j)\log_2\big(p(i,j)+\epsilon\big)}`

    Entropy is a measure of the randomness/variability in neighborhood intensity values.
    """
    ent = self.coefficients['HXY']
    return (ent.mean())

  def getHomogeneity1FeatureValue(self):
    r"""
    Calculate and return the mean Homogeneity 1.

    :math:`homogeneity\ 1 = \displaystyle\sum^{N_g}_{i=1}\displaystyle\sum^{N_g}_{j=1}{\frac{p(i,j)}{1+|i-j|}}`

    Homogeneity 1 is a measure of the similarity in intensity values for
    neighboring voxels. It is a measure of local homogeneity that increases
    with less contrast in the window.
    """
    i = self.coefficients['i']
    j = self.coefficients['j']
    homo1 = numpy.sum((self.P_glcm / (1 + (numpy.abs(i - j))[:, :, None])), (0, 1))
    return (homo1.mean())

  def getHomogeneity2FeatureValue(self):
    r"""
    Calculate and return the mean Homogeneity 2.

    :math:`homogeneity\ 2 = \displaystyle\sum^{N_g}_{i=1}\displaystyle\sum^{N_g}_{j=1}{\frac{p(i,j)}{1+|i-j|^2}}`

    Homogeneity 2 is a measure of the similarity in intensity values
    for neighboring voxels.
    """
    i = self.coefficients['i']
    j = self.coefficients['j']
    homo2 = numpy.sum((self.P_glcm / (1 + (numpy.abs(i - j))[:, :, None] ** 2)), (0, 1))
    return (homo2.mean())

  def getImc1FeatureValue(self):
    r"""
    Using coefficients :math:`HX`, :math:`HY`, :math:`HXY` and :math:`HXY1`, calculate and return the mean Informal
    Measure of Correlation 1.

    :math:`IMC\ 1 = \frac{HXY-HXY1}{\max\{HX,HY\}}`
    """
    HX = self.coefficients['HX']
    HY = self.coefficients['HY']
    HXY = self.coefficients['HXY']
    HXY1 = self.coefficients['HXY1']
    imc1 = (HXY - HXY1) / numpy.max(([HX, HY]), 0)
    return (imc1.mean())

  def getImc2FeatureValue(self):
    r"""
    Using coefficients :math:`HXY` and :math:`HXY2`, calculate and return the mean Informal Measure of Correlation 2.

    :math:`IMC\ 2 = \sqrt{1-e^{-2(HXY2-HXY)}}`
    """
    HXY = self.coefficients['HXY']
    HXY2 = self.coefficients['HXY2']

    imc2 = (1 - numpy.e ** (-2 * (HXY2 - HXY))) ** (0.5)  # matlab:(1-exp(-2*(hxy2-hxy)))^0.5;

    return (imc2.mean())

  def getIdmFeatureValue(self):
    r"""
    Calculate and return the mean Inverse Difference Moment.

    :math:`IDM = \displaystyle\sum^{N_g}_{i=1}\displaystyle\sum^{N_g}_{j=1}{ \frac{\textbf{P}(i,j)}{1+|i-j|^2} }`

    IDM (inverse difference moment)  is a measure of the local
    homogeneity of an image. IDM weights are the inverse of the Contrast
    weights (decreasing exponentially from the diagonal i=j in the GLCM).
    """
    i = self.coefficients['i']
    j = self.coefficients['j']
    idm = numpy.sum((self.P_glcm / (1 + (((numpy.abs(i - j))[:, :, None] ** 2)))), (0, 1))
    return (idm.mean())

  def getIdmnFeatureValue(self):
    r"""
    Calculate and return the mean Inverse Difference Moment Normalized.

    :math:`IDMN = \displaystyle\sum^{N_g}_{i=1}\displaystyle\sum^{N_g}_{j=1}{ \frac{p(i,j)}{1+\left(\frac{|i-j|^2}{N_g^2}\right)} }`

    IDMN (inverse difference moment normalized)  is a measure of the local
    homogeneity of an image. IDMN weights are the inverse of the Contrast
    weights (decreasing exponentially from the diagonal :math:`i=j` in the GLCM).
    Unlike Homogeneity2, IDMN normalizes the square of the difference between
    neighboring intensity values by dividing over the square of the total
    number of discrete intensity values.
    """
    i = self.coefficients['i']
    j = self.coefficients['j']
    Ng = self.coefficients['Ng']
    idmn = numpy.sum((self.P_glcm / (1 + (((numpy.abs(i - j))[:, :, None] ** 2) / (Ng ** 2)))), (0, 1))
    return (idmn.mean())

  def getIdFeatureValue(self):
    r"""
    Calculate and return the mean Inverse Difference.

    :math:`ID = \displaystyle\sum^{N_g}_{i=1}\displaystyle\sum^{N_g}_{j=1}{ \frac{\textbf{P}(i,j)}{1+|i-j|} }`

    ID (inverse difference) is another measure of the local homogeneity of an image.
    With more uniform gray levels, the denominator will remain low, resulting in a higher overall value.
    """
    i = self.coefficients['i']
    j = self.coefficients['j']
    id = numpy.sum((self.P_glcm / (1 + ((numpy.abs(i - j))[:, :, None]))), (0, 1))
    return (id.mean())

  def getIdnFeatureValue(self):
    r"""
    Calculate and return the mean Inverse Difference Normalized.

    :math:`IDN = \displaystyle\sum^{N_g}_{i=1}\displaystyle\sum^{N_g}_{j=1}{ \frac{p(i,j)}{1+\left(\frac{|i-j|}{N_g}\right)} }`

    IDN (inverse difference normalized) is another measure of the local
    homogeneity of an image. Unlike Homogeneity1, IDN normalizes the difference
    between the neighboring intensity values by dividing over the total number
    of discrete intensity values.
    """
    i = self.coefficients['i']
    j = self.coefficients['j']
    Ng = self.coefficients['Ng']
    idn = numpy.sum((self.P_glcm / (1 + ((numpy.abs(i - j))[:, :, None] / Ng))), (0, 1))
    return (idn.mean())

  def getInverseVarianceFeatureValue(self):
    r"""
    Calculate and return the mean Inverse Variance.

    :math:`inverse\ variance = \displaystyle\sum^{N_g}_{i=1}\displaystyle\sum^{N_g}_{j=1}{\frac{p(i,j)}{|i-j|^2}}, i \neq j`
    """
    i = self.coefficients['i']
    j = self.coefficients['j']
    maskDiags = numpy.abs(i - j) > 0
    inv = numpy.sum((self.P_glcm[maskDiags] / ((numpy.abs(i - j))[:, :, None] ** 2)[maskDiags]), 0)
    return (inv.mean())

  def getMaximumProbabilityFeatureValue(self):
    r"""
    Calculate and return the mean Maximum Probability.

    :math:`maximum\ probability = \max\big(p(i,j)\big)`

    Maximum Probability is occurrences of the most predominant pair of
    neighboring intensity values.
    """
    maxprob = self.P_glcm.max((0, 1))
    return (maxprob.mean())

  def getSumAverageFeatureValue(self):
    r"""
    Coefficient :math:`p_{x+y}`, calculate and return the mean Sum Average.

    :math:`sum\ average = \displaystyle\sum^{2N_g}_{k=2}{p_{x+y}(k)k}`

    Sum Average measures the relationship between occurrences of pairs
    with lower intensity values and occurrences of pairs with higher intensity
    values.
    """
    pxAddy = self.coefficients['pxAddy']
    kValuesSum = self.coefficients['kValuesSum']
    sumavg = numpy.sum((kValuesSum[:, None] * pxAddy), 0)
    return (sumavg.mean())

  def getSumEntropyFeatureValue(self):
    r"""
    Using coefficient :math:`p_{x+y}`, calculate and return the mean Sum Entropy.

    :math:`sum\ entropy = \displaystyle\sum^{2N_g}_{k=2}{p_{x+y}(k)\log_2\big(p_{x+y}(k)+\epsilon\big)}`

    Sum Entropy is a sum of neighborhood intensity value differences.
    """
    pxAddy = self.coefficients['pxAddy']
    eps = self.coefficients['eps']
    sumentr = (-1) * numpy.sum((pxAddy * numpy.log2(pxAddy + eps)), 0)
    return (sumentr.mean())

  def getSumVarianceFeatureValue(self):
    r"""
    Using coefficients :math:`p_{x+y}` and SumEntropy (SE) calculate and return the mean Sum Variance.

    :math:`sum\ variance = \displaystyle\sum^{2N_g}_{k=2}{(k-SE)^2p_{x+y}(k)}`

    Sum Variance is a measure of heterogeneity that places higher weights on
    neighboring intensity level pairs that deviate more from the mean.
    """
    eps = self.coefficients['eps']
    pxAddy = self.coefficients['pxAddy']
    kValuesSum = self.coefficients['kValuesSum']
    sumentr = (-1) * numpy.sum((pxAddy * numpy.log2(pxAddy + eps)), 0, keepdims=True)
    sumvar = numpy.sum((pxAddy * ((kValuesSum[:, None] - sumentr) ** 2)), 0)
    return (sumvar.mean())

  def getSumVariance2FeatureValue(self):
    r"""
    Using coefficients :math:`p_{x+y}` and SumAvarage (SA) calculate and return the mean Sum Variance 2.
    :math:`sum\ variance\ 2 = \displaystyle\sum^{2N_g}_{k=2}{(k-SA)^2p_{x+y}(k)}`

    Sum Variance 2 is a measure of heterogeneity that places higher weights on
    neighboring intensity level pairs that deviate more from the mean.

    This formula differs from SumVariance in that instead of subtracting the SumEntropy from the intensity,
    it subtracts the SumAvarage, which is the mean of intensities and not its entropy
    """
    pxAddy = self.coefficients['pxAddy']
    kValuesSum = self.coefficients['kValuesSum']
    sumavg = numpy.sum((kValuesSum[:, None] * pxAddy), 0, keepdims=True)
    sumvar = numpy.sum((pxAddy * ((kValuesSum[:, None] - sumavg) ** 2)), 0)
    return (sumvar.mean())

  def getSumSquaresFeatureValue(self):
    r"""
    Using coefficients :math:`i` and math:`\mu_x`, calculate and return the mean Sum of Squares (also known as
    Variance) of the :math:`i` distribution.

    :math:`sum\ squares = \displaystyle\sum^{N_g}_{i=1}\displaystyle\sum^{N_g}_{j=1}{(i-\mu_x)^2p(i,j)}`

    Sum of Squares or Variance is a measure in the distribution of neigboring intensity level pairs
    about the mean intensity level in the GLCM.

    N.B. This formula represents the variance of the distribution of :math:`i` and is independent from the distribution
    of :math:`j`. Therefore, only use this formula if the GLCM is symmetrical, where VAR(i) is equal to VAR(j).
    """
    i = self.coefficients['i']
    ux = self.coefficients['ux']
    # Also known as Variance
    ss = numpy.sum((self.P_glcm * ((i[:, :, None] - ux) ** 2)), (0, 1))
    return (ss.mean())
