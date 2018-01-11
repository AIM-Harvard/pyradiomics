import numpy

from radiomics import base, cMatrices


class RadiomicsGLRLM(base.RadiomicsFeaturesBase):
  r"""
  A Gray Level Run Length Matrix (GLRLM) quantifies gray level runs, which are defined as the length in number of
  pixels, of consecutive pixels that have the same gray level value. In a gray level run length matrix
  :math:`\textbf{P}(i,j|\theta)`, the :math:`(i,j)^{\text{th}}` element describes the number of runs with gray level
  :math:`i` and length :math:`j` occur in the image (ROI) along angle :math:`\theta`.

  As a two dimensional example, consider the following 5x5 image, with 5 discrete gray levels:

  .. math::
    \textbf{I} = \begin{bmatrix}
    5 & 2 & 5 & 4 & 4\\
    3 & 3 & 3 & 1 & 3\\
    2 & 1 & 1 & 1 & 3\\
    4 & 2 & 2 & 2 & 3\\
    3 & 5 & 3 & 3 & 2 \end{bmatrix}

  The GLRLM for :math:`\theta = 0`, where 0 degrees is the horizontal direction, then becomes:

  .. math::
    \textbf{P} = \begin{bmatrix}
    1 & 0 & 1 & 0 & 0\\
    3 & 0 & 1 & 0 & 0\\
    4 & 1 & 1 & 0 & 0\\
    1 & 1 & 0 & 0 & 0\\
    3 & 0 & 0 & 0 & 0 \end{bmatrix}

  Let:

  - :math:`N_g` be the number of discreet intensity values in the image
  - :math:`N_r` be the number of discreet run lengths in the image
  - :math:`N_p` be the number of voxels in the image
  - :math:`N_z(\theta)` be the number of runs in the image along angle :math:`\theta`, which is equal to
    :math:`\sum^{N_g}_{i=1}\sum^{N_r}_{j=1}{\textbf{P}(i,j|\theta)}` and :math:`1 \leq N_z(\theta) \leq N_p`
  - :math:`\textbf{P}(i,j|\theta)` be the run length matrix for an arbitrary direction :math:`\theta`
  - :math:`p(i,j|\theta)` be the normalized run length matrix, defined as :math:`p(i,j|\theta) =
    \frac{\textbf{P}(i,j|\theta)}{N_z(\theta)}`

  By default, the value of a feature is calculated on the GLRLM for each angle separately, after which the mean of these
  values is returned. If distance weighting is enabled, GLRLMs are weighted by the distance between neighbouring voxels
  and then summed and normalised. Features are then calculated on the resultant matrix. The distance between
  neighbouring voxels is calculated for each angle using the norm specified in 'weightingNorm'.

  The following class specific settings are possible:

  - weightingNorm [None]: string, indicates which norm should be used when applying distance weighting.
    Enumerated setting, possible values:

    - 'manhattan': first order norm
    - 'euclidean': second order norm
    - 'infinity': infinity norm.
    - 'no_weighting': GLCMs are weighted by factor 1 and summed
    - None: Applies no weighting, mean of values calculated on separate matrices is returned.

    In case of other values, an warning is logged and option 'no_weighting' is used.

  References

  - Galloway MM. 1975. Texture analysis using gray level run lengths. Computer Graphics and Image Processing,
    4(2):172-179.
  - Chu A., Sehgal C.M., Greenleaf J. F. 1990. Use of gray value distribution of run length for texture analysis.
    Pattern Recognition Letters, 11(6):415-419
  - Xu D., Kurani A., Furst J., Raicu D. 2004. Run-Length Encoding For Volumetric Texture. International Conference on
    Visualization, Imaging and Image Processing (VIIP), p. 452-458
  - Tang X. 1998. Texture information in run-length matrices. IEEE Transactions on Image Processing 7(11):1602-1609.
  - `Tustison N., Gee J. Run-Length Matrices For Texture Analysis. Insight Journal 2008 January - June.
    <http://www.insight-journal.org/browse/publication/231>`_
  """

  def __init__(self, inputImage, inputMask, **kwargs):
    super(RadiomicsGLRLM, self).__init__(inputImage, inputMask, **kwargs)

    self.weightingNorm = kwargs.get('weightingNorm', None)  # manhattan, euclidean, infinity

    self.P_glrlm = None
    self._applyBinning()

  def _initCalculation(self):
    self.coefficients['Nr'] = numpy.max(self.matrix.shape)
    self.coefficients['Np'] = len(self.labelledVoxelCoordinates[0])

    self.P_glrlm = self._calculateMatrix()

    self._calculateCoefficients()

    self.logger.debug('GLRLM feature class initialized, calculated GLRLM with shape %s', self.P_glrlm.shape)

  def _calculateMatrix(self):
    self.logger.debug('Calculating GLRLM matrix in C')
    P_glrlm, angles = cMatrices.calculate_glrlm(self.matrix,
                                                self.maskArray,
                                                self.coefficients['Ng'],
                                                self.coefficients['Nr'],
                                                self.settings.get('force2D', False),
                                                self.settings.get('force2Ddimension', 0))

    self.logger.debug('Process calculated matrix')

    # Optionally apply a weighting factor
    if self.weightingNorm is not None:
      self.logger.debug('Applying weighting (%s)', self.weightingNorm)
      # Correct the number of voxels for the number of times it is used (once per angle), affects run percentage
      self.coefficients['Np'] *= len(angles)

      pixelSpacing = self.inputImage.GetSpacing()[::-1]
      weights = numpy.empty(len(angles))
      for a_idx, a in enumerate(angles):
        if self.weightingNorm == 'infinity':
          weights[a_idx] = max(numpy.abs(a) * pixelSpacing)
        elif self.weightingNorm == 'euclidean':
          weights[a_idx] = numpy.sqrt(numpy.sum((numpy.abs(a) * pixelSpacing) ** 2))
        elif self.weightingNorm == 'manhattan':
          weights[a_idx] = numpy.sum(numpy.abs(a) * pixelSpacing)
        elif self.weightingNorm == 'no_weighting':
          weights[a_idx] = 1
        else:
          self.logger.warning('weigthing norm "%s" is unknown, weighting factor is set to 1', self.weightingNorm)
          weights[a_idx] = 1

      P_glrlm = numpy.sum(P_glrlm * weights[None, None, :], 2, keepdims=True)

    Nz = numpy.sum(P_glrlm, (0, 1))

    # Delete empty angles if no weighting is applied
    if P_glrlm.shape[2] > 1:
      emptyAngles = numpy.where(Nz == 0)
      if len(emptyAngles[0]) > 0:  # One or more angles are 'empty'
        self.logger.debug('Deleting %d empty angles:\n%s', len(emptyAngles[0]), angles[emptyAngles])
        P_glrlm = numpy.delete(P_glrlm, emptyAngles, 2)
        Nz = numpy.delete(Nz, emptyAngles, 0)
      else:
        self.logger.debug('No empty angles')

    self.coefficients['Nz'] = Nz

    return P_glrlm

  def _calculateCoefficients(self):
    self.logger.debug('Calculating GLRLM coefficients')

    pr = numpy.sum(self.P_glrlm, 0)
    pg = numpy.sum(self.P_glrlm, 1)

    ivector = self.coefficients['grayLevels']
    jvector = numpy.arange(1, self.P_glrlm.shape[1] + 1, dtype=numpy.float64)

    emptyGrayLevels = numpy.where(numpy.sum(pg, 1) == 0)
    emptyRunLenghts = numpy.where(numpy.sum(pr, 1) == 0)

    self.P_glrlm = numpy.delete(self.P_glrlm, emptyGrayLevels, 0)
    self.P_glrlm = numpy.delete(self.P_glrlm, emptyRunLenghts, 1)
    jvector = numpy.delete(jvector, emptyRunLenghts)

    pg = numpy.delete(pg, emptyGrayLevels, 0)
    pr = numpy.delete(pr, emptyRunLenghts, 0)

    self.coefficients['pr'] = pr
    self.coefficients['pg'] = pg
    self.coefficients['ivector'] = ivector
    self.coefficients['jvector'] = jvector

  def getShortRunEmphasisFeatureValue(self):
    r"""
    **1. Short Run Emphasis (SRE)**

    .. math::
      \textit{SRE} = \frac{\sum^{N_g}_{i=1}\sum^{N_r}_{j=1}{\frac{\textbf{P}(i,j|\theta)}{j^2}}}{N_z(\theta)}

    SRE is a measure of the distribution of short run lengths, with a greater value indicative of shorter run lengths
    and more fine textural textures.
    """
    pr = self.coefficients['pr']
    jvector = self.coefficients['jvector']
    Nz = self.coefficients['Nz']

    sre = numpy.sum((pr / (jvector[:, None] ** 2)), 0) / Nz
    return sre.mean()

  def getLongRunEmphasisFeatureValue(self):
    r"""
    **2. Long Run Emphasis (LRE)**

    .. math::
      \textit{LRE} = \frac{\sum^{N_g}_{i=1}\sum^{N_r}_{j=1}{\textbf{P}(i,j|\theta)j^2}}{N_z(\theta)}

    LRE is a measure of the distribution of long run lengths, with a greater value indicative of longer run lengths and
    more coarse structural textures.
    """
    pr = self.coefficients['pr']
    jvector = self.coefficients['jvector']
    Nz = self.coefficients['Nz']

    lre = numpy.sum((pr * (jvector[:, None] ** 2)), 0) / Nz
    return lre.mean()

  def getGrayLevelNonUniformityFeatureValue(self):
    r"""
    **3. Gray Level Non-Uniformity (GLN)**

    .. math::
      \textit{GLN} = \frac{\sum^{N_g}_{i=1}\left(\sum^{N_r}_{j=1}{\textbf{P}(i,j|\theta)}\right)^2}{N_z(\theta)}

    GLN measures the similarity of gray-level intensity values in the image, where a lower GLN value correlates with a
    greater similarity in intensity values.
    """
    pg = self.coefficients['pg']
    Nz = self.coefficients['Nz']

    gln = numpy.sum((pg ** 2), 0) / Nz
    return gln.mean()

  def getGrayLevelNonUniformityNormalizedFeatureValue(self):
    r"""
    **4. Gray Level Non-Uniformity Normalized (GLNN)**

    .. math::
      \textit{GLNN} = \frac{\sum^{N_g}_{i=1}\left(\sum^{N_r}_{j=1}{\textbf{P}(i,j|\theta)}\right)^2}{N_z(\theta)^2}

    GLNN measures the similarity of gray-level intensity values in the image, where a lower GLNN value correlates with a
    greater similarity in intensity values. This is the normalized version of the GLN formula.
    """
    pg = self.coefficients['pg']
    Nz = self.coefficients['Nz']

    glnn = numpy.sum(pg ** 2, 0) / (Nz ** 2)
    return glnn.mean()

  def getRunLengthNonUniformityFeatureValue(self):
    r"""
    **5. Run Length Non-Uniformity (RLN)**

    .. math::
      \textit{RLN} = \frac{\sum^{N_r}_{j=1}\left(\sum^{N_g}_{i=1}{\textbf{P}(i,j|\theta)}\right)^2}{N_z(\theta)}

    RLN measures the similarity of run lengths throughout the image, with a lower value indicating more homogeneity
    among run lengths in the image.
    """
    pr = self.coefficients['pr']
    Nz = self.coefficients['Nz']

    rln = numpy.sum((pr ** 2), 0) / Nz
    return rln.mean()

  def getRunLengthNonUniformityNormalizedFeatureValue(self):
    r"""
    **6. Run Length Non-Uniformity Normalized (RLNN)**

    .. math::
      \textit{RLNN} = \frac{\sum^{N_r}_{j=1}\left(\sum^{N_g}_{i=1}{\textbf{P}(i,j|\theta)}\right)^2}{N_z(\theta)^2}

    RLNN measures the similarity of run lengths throughout the image, with a lower value indicating more homogeneity
    among run lengths in the image. This is the normalized version of the RLN formula.
    """
    pr = self.coefficients['pr']
    Nz = self.coefficients['Nz']

    rlnn = numpy.sum((pr ** 2), 0) / Nz ** 2
    return rlnn.mean()

  def getRunPercentageFeatureValue(self):
    r"""
    **7. Run Percentage (RP)**

    .. math::
      \textit{RP} = {\frac{N_z(\theta)}{N_p}}

    RP measures the coarseness of the texture by taking the ratio of number of runs and number of voxels in the ROI.

    Values are in range :math:`\frac{1}{N_p} \leq RP \leq 1`, with higher values indicating a larger portion of the ROI
    consists of short runs (indicates a more fine texture).

    .. note::
      Note that when weighting is applied and matrices are merged before calculation, :math:`N_p` is multiplied by
      :math:`n` number of matrices merged to ensure correct normalization (as each voxel is considered :math:`n` times)
    """
    Np = self.coefficients['Np']
    Nz = self.coefficients['Nz']

    rp = Nz / Np
    return rp.mean()

  def getGrayLevelVarianceFeatureValue(self):
    r"""
    **8. Gray Level Variance (GLV)**

    .. math::
      \textit{GLV} = \displaystyle\sum^{N_g}_{i=1}\displaystyle\sum^{N_r}_{j=1}{p(i,j|\theta)(i - \mu)^2}

    Here, :math:`\mu = \displaystyle\sum^{N_g}_{i=1}\displaystyle\sum^{N_r}_{j=1}{p(i,j|\theta)i}`

    GLV measures the variance in gray level intensity for the runs.
    """
    ivector = self.coefficients['ivector']
    Nz = self.coefficients['Nz']
    pg = self.coefficients['pg'] / Nz  # divide by Nz to get the normalized matrix

    u_i = numpy.sum(pg * ivector[:, None], 0)
    glv = numpy.sum(pg * (ivector[:, None] - u_i[None, :]) ** 2, 0)
    return glv.mean()

  def getRunVarianceFeatureValue(self):
    r"""
    **9. Run Variance (RV)**

    .. math::
      \textit{RV} = \displaystyle\sum^{N_g}_{i=1}\displaystyle\sum^{N_r}_{j=1}{p(i,j|\theta)(j - \mu)^2}

    Here, :math:`\mu = \displaystyle\sum^{N_g}_{i=1}\displaystyle\sum^{N_r}_{j=1}{p(i,j|\theta)j}`

    RV is a measure of the variance in runs for the run lengths.
    """
    jvector = self.coefficients['jvector']
    Nz = self.coefficients['Nz']
    pr = self.coefficients['pr'] / Nz   # divide by Nz to get the normalized matrix

    u_j = numpy.sum(pr * jvector[:, None], 0)
    rv = numpy.sum(pr * (jvector[:, None] - u_j[None, :]) ** 2, 0)
    return rv.mean()

  def getRunEntropyFeatureValue(self):
    r"""
    **10. Run Entropy (RE)**

    .. math::
      \textit{RE} = -\displaystyle\sum^{N_g}_{i=1}\displaystyle\sum^{N_r}_{j=1}
      {p(i,j|\theta)\log_{2}(p(i,j|\theta)+\epsilon)}

    Here, :math:`\epsilon` is an arbitrarily small positive number (:math:`\approx 2.2\times10^{-16}`).

    RE measures the uncertainty/randomness in the distribution of run lengths and gray levels. A higher value indicates
    more heterogeneity in the texture patterns.
    """
    eps = numpy.spacing(1)
    Nz = self.coefficients['Nz']
    p_glrlm = self.P_glrlm / Nz  # divide by Nz to get the normalized matrix

    re = -numpy.sum(p_glrlm * numpy.log2(p_glrlm + eps), (0, 1))
    return re.mean()

  def getLowGrayLevelRunEmphasisFeatureValue(self):
    r"""
    **11. Low Gray Level Run Emphasis (LGLRE)**

    .. math::
      \textit{LGLRE} = \frac{\sum^{N_g}_{i=1}\sum^{N_r}_{j=1}{\frac{\textbf{P}(i,j|\theta)}{i^2}}}{N_z(\theta)}

    LGLRE measures the distribution of low gray-level values, with a higher value indicating a greater concentration of
    low gray-level values in the image.
    """
    pg = self.coefficients['pg']
    ivector = self.coefficients['ivector']
    Nz = self.coefficients['Nz']

    lglre = numpy.sum((pg / (ivector[:, None] ** 2)), 0) / Nz
    return lglre.mean()

  def getHighGrayLevelRunEmphasisFeatureValue(self):
    r"""
    **12. High Gray Level Run Emphasis (HGLRE)**

    .. math::
      \textit{HGLRE} = \frac{\sum^{N_g}_{i=1}\sum^{N_r}_{j=1}{\textbf{P}(i,j|\theta)i^2}}{N_z(\theta)}

    HGLRE measures the distribution of the higher gray-level values, with a higher value indicating a greater
    concentration of high gray-level values in the image.
    """
    pg = self.coefficients['pg']
    ivector = self.coefficients['ivector']
    Nz = self.coefficients['Nz']

    hglre = numpy.sum((pg * (ivector[:, None] ** 2)), 0) / Nz
    return hglre.mean()

  def getShortRunLowGrayLevelEmphasisFeatureValue(self):
    r"""
    **13. Short Run Low Gray Level Emphasis (SRLGLE)**

    .. math::
      \textit{SRLGLE} = \frac{\sum^{N_g}_{i=1}\sum^{N_r}_{j=1}{\frac{\textbf{P}(i,j|\theta)}{i^2j^2}}}{N_z(\theta)}

    SRLGLE measures the joint distribution of shorter run lengths with lower gray-level values.
    """
    ivector = self.coefficients['ivector']
    jvector = self.coefficients['jvector']
    Nz = self.coefficients['Nz']

    srlgle = numpy.sum((self.P_glrlm / ((ivector[:, None, None] ** 2) * (jvector[None, :, None] ** 2))),
                       (0, 1)) / Nz
    return srlgle.mean()

  def getShortRunHighGrayLevelEmphasisFeatureValue(self):
    r"""
    **14. Short Run High Gray Level Emphasis (SRHGLE)**

    .. math::
      \textit{SRHGLE} = \frac{\sum^{N_g}_{i=1}\sum^{N_r}_{j=1}{\frac{\textbf{P}(i,j|\theta)i^2}{j^2}}}{N_z(\theta)}

    SRHGLE measures the joint distribution of shorter run lengths with higher gray-level values.
    """
    ivector = self.coefficients['ivector']
    jvector = self.coefficients['jvector']
    Nz = self.coefficients['Nz']

    srhgle = numpy.sum((self.P_glrlm * (ivector[:, None, None] ** 2) / (jvector[None, :, None] ** 2)),
                       (0, 1)) / Nz
    return srhgle.mean()

  def getLongRunLowGrayLevelEmphasisFeatureValue(self):
    r"""
    **15. Long Run Low Gray Level Emphasis (LRLGLE)**

    .. math::
      \textit{LRLGLRE} = \frac{\sum^{N_g}_{i=1}\sum^{N_r}_{j=1}{\frac{\textbf{P}(i,j|\theta)j^2}{i^2}}}{N_z(\theta)}

    LRLGLRE measures the joint distribution of long run lengths with lower gray-level values.
    """
    ivector = self.coefficients['ivector']
    jvector = self.coefficients['jvector']
    Nz = self.coefficients['Nz']

    lrlgle = numpy.sum((self.P_glrlm * (jvector[None, :, None] ** 2) / (ivector[:, None, None] ** 2)),
                       (0, 1)) / Nz
    return lrlgle.mean()

  def getLongRunHighGrayLevelEmphasisFeatureValue(self):
    r"""
    **16. Long Run High Gray Level Emphasis (LRHGLE)**

    .. math::
      \textit{LRHGLRE} = \frac{\sum^{N_g}_{i=1}\sum^{N_r}_{j=1}{\textbf{P}(i,j|\theta)i^2j^2}}{N_z(\theta)}

    LRHGLRE measures the joint distribution of long run lengths with higher gray-level values.
    """
    ivector = self.coefficients['ivector']
    jvector = self.coefficients['jvector']
    Nz = self.coefficients['Nz']

    lrhgle = numpy.sum((self.P_glrlm * ((jvector[None, :, None] ** 2) * (ivector[:, None, None] ** 2))),
                       (0, 1)) / Nz
    return lrhgle.mean()
