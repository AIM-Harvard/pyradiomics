import numpy

from radiomics import base, cMatrices, deprecated


class RadiomicsGLDM(base.RadiomicsFeaturesBase):
  r"""
  A Gray Level Dependence Matrix (GLDM) quantifies gray level dependencies in an image.
  A gray level dependency is defined as a the number of connected voxels within distance :math:`\delta` that are
  dependent on the center voxel.
  A neighbouring voxel with gray level :math:`j` is considered dependent on center voxel with gray level :math:`i`
  if :math:`|i-j|\le\alpha`. In a gray level dependence matrix :math:`\textbf{P}(i,j)` the :math:`(i,j)`\ :sup:`th`
  element describes the number of times a voxel with gray level :math:`i` with :math:`j` dependent voxels
  in its neighbourhood appears in image.

  As a two dimensional example, consider the following 5x5 image, with 5 discrete gray levels:

  .. math::
    \textbf{I} = \begin{bmatrix}
    5 & 2 & 5 & 4 & 4\\
    3 & 3 & 3 & 1 & 3\\
    2 & 1 & 1 & 1 & 3\\
    4 & 2 & 2 & 2 & 3\\
    3 & 5 & 3 & 3 & 2 \end{bmatrix}

  For :math:`\alpha=0` and :math:`\delta = 1`, the GLDM then becomes:

  .. math::
    \textbf{P} = \begin{bmatrix}
    0 & 1 & 2 & 1 \\
    1 & 2 & 3 & 0 \\
    1 & 4 & 4 & 0 \\
    1 & 2 & 0 & 0 \\
    3 & 0 & 0 & 0 \end{bmatrix}

  Let:

  - :math:`N_g` be the number of discreet intensity values in the image
  - :math:`N_d` be the number of discreet dependency sizes in the image
  - :math:`N_z` be the number of dependency zones in the image, which is equal to
    :math:`\sum^{N_g}_{i=1}\sum^{N_d}_{j=1}{\textbf{P}(i,j)}`
  - :math:`\textbf{P}(i,j)` be the dependence matrix
  - :math:`p(i,j)` be the normalized dependence matrix, defined as :math:`p(i,j) = \frac{\textbf{P}(i,j)}{N_z}`

  .. note::
    Because incomplete zones are allowed, every voxel in the ROI has a dependency zone. Therefore, :math:`N_z = N_p`,
    where :math:`N_p` is the number of voxels in the image.
    Due to the fact that :math:`Nz = N_p`, the Dependence Percentage and Gray Level Non-Uniformity Normalized (GLNN)
    have been removed. The first because it would always compute to 1, the latter because it is mathematically equal to
    first order - Uniformity (see :py:func:`~radiomics.firstorder.RadiomicsFirstOrder.getUniformityFeatureValue()`). For
    mathematical proofs, see :ref:`here <radiomics-excluded-gldm-label>`.

  The following class specific settings are possible:

  - distances [[1]]: List of integers. This specifies the distances between the center voxel and the neighbor, for which
    angles should be generated.
  - gldm_a [0]: float, :math:`\alpha` cutoff value for dependence. A neighbouring voxel with gray level :math:`j` is considered
    dependent on center voxel with gray level :math:`i` if :math:`|i-j|\le\alpha`

  References:

  - Sun C, Wee WG. Neighboring Gray Level Dependence Matrix for Texture Classification. Comput Vision,
    Graph Image Process. 1983;23:341-352
  """

  def __init__(self, inputImage, inputMask, **kwargs):
    super(RadiomicsGLDM, self).__init__(inputImage, inputMask, **kwargs)

    self.gldm_a = kwargs.get('gldm_a', 0)

    self.P_gldm = None
    self._applyBinning()

  def _initCalculation(self):
    self.coefficients['Np'] = len(self.labelledVoxelCoordinates[0])

    self.P_gldm = self._calculateMatrix()

    self.logger.debug('Feature class initialized, calculated GLDM with shape %s', self.P_gldm.shape)

  def _calculateMatrix(self):
    P_gldm = cMatrices.calculate_gldm(self.matrix,
                                      self.maskArray,
                                      numpy.array(self.settings.get('distances', [1])),
                                      self.coefficients['Ng'],
                                      self.gldm_a,
                                      self.settings.get('force2D', False),
                                      self.settings.get('force2Ddimension', 0))

    jvector = numpy.arange(1, P_gldm.shape[1] + 1, dtype='float64')

    # Delete rows and columns that specify gray levels not present in the ROI
    pd = numpy.sum(P_gldm, 0)
    pg = numpy.sum(P_gldm, 1)

    P_gldm = numpy.delete(P_gldm, numpy.where(pg == 0), 0)
    P_gldm = numpy.delete(P_gldm, numpy.where(pd == 0), 1)

    jvector = numpy.delete(jvector, numpy.where(pd == 0))

    pg = numpy.delete(pg, numpy.where(pg == 0))
    pd = numpy.delete(pd, numpy.where(pd == 0))

    self.coefficients['pd'] = pd
    self.coefficients['pg'] = pg

    self.coefficients['ivector'] = self.coefficients['grayLevels']
    self.coefficients['jvector'] = jvector

    return P_gldm

  def getSmallDependenceEmphasisFeatureValue(self):
    r"""
    **1. Small Dependence Emphasis (SDE)**

    .. math::
      SDE = \frac{\sum^{N_g}_{i=1}\sum^{N_d}_{j=1}{\frac{\textbf{P}(i,j)}{i^2}}}{N_z}

    A measure of the distribution of small dependencies, with a greater value indicative
    of smaller dependence and less homogeneous textures.
    """
    pd = self.coefficients['pd']
    jvector = self.coefficients['jvector']
    Nz = self.coefficients['Np']  # Nz = Np, see class docstring

    try:
      sde = numpy.sum(pd / (jvector ** 2)) / Nz
      return sde
    except ZeroDivisionError:
      return numpy.core.nan

  def getLargeDependenceEmphasisFeatureValue(self):
    r"""
    **2. Large Dependence Emphasis (LDE)**

    .. math::
      LDE = \frac{\sum^{N_g}_{i=1}\sum^{N_d}_{j=1}{\textbf{P}(i,j)j^2}}{N_z}

    A measure of the distribution of large dependencies, with a greater value indicative
    of larger dependence and more homogeneous textures.
    """
    pd = self.coefficients['pd']
    jvector = self.coefficients['jvector']
    Nz = self.coefficients['Np']  # Nz = Np, see class docstring

    try:
      lre = numpy.sum(pd * (jvector ** 2)) / Nz
      return lre
    except ZeroDivisionError:
      return numpy.core.nan

  def getGrayLevelNonUniformityFeatureValue(self):
    r"""
    **3. Gray Level Non-Uniformity (GLN)**

    .. math::
      GLN = \frac{\sum^{N_g}_{i=1}\left(\sum^{N_d}_{j=1}{\textbf{P}(i,j)}\right)^2}{N_z}

    Measures the similarity of gray-level intensity values in the image, where a lower GLN value
    correlates with a greater similarity in intensity values.
    """
    pg = self.coefficients['pg']
    Nz = self.coefficients['Np']  # Nz = Np, see class docstring

    try:
      gln = numpy.sum(pg ** 2) / Nz
      return gln
    except ZeroDivisionError:
      return numpy.core.nan

  @deprecated
  def getGrayLevelNonUniformityNormalizedFeatureValue(self):
    r"""
    **DEPRECATED. Gray Level Non-Uniformity Normalized (GLNN)**

    :math:`GLNN = \frac{\sum^{N_g}_{i=1}\left(\sum^{N_d}_{j=1}{\textbf{P}(i,j)}\right)^2}{\sum^{N_g}_{i=1}\sum^{N_d}_{j=1}{\textbf{P}(i,j)}^2}`

    .. warning::
      This feature has been deprecated, as it is mathematically equal to First Order - Uniformity
      :py:func:`~radiomics.firstorder.RadiomicsFirstOrder.getUniformityFeatureValue()`.
      See :ref:`here <radiomics-excluded-gldm-glnn-label>` for the proof. **Enabling this feature will result in the
      logging of a DeprecationWarning (does not interrupt extraction of other features), no value is calculated for this features**
    """
    raise DeprecationWarning('GLDM - Gray Level Non-Uniformity Normalized is mathematically equal to First Order - Uniformity, '
                             'see http://pyradiomics.readthedocs.io/en/latest/removedfeatures.html for more details')

  def getDependenceNonUniformityFeatureValue(self):
    r"""
    **4. Dependence Non-Uniformity (DN)**

    .. math::
      DN = \frac{\sum^{N_d}_{j=1}\left(\sum^{N_g}_{i=1}{\textbf{P}(i,j)}\right)^2}{N_z}

    Measures the similarity of dependence throughout the image, with a lower value indicating
    more homogeneity among dependencies in the image.
    """
    pd = self.coefficients['pd']
    Nz = self.coefficients['Np']  # Nz = Np, see class docstring

    try:
      dn = numpy.sum(pd ** 2) / Nz
      return dn
    except ZeroDivisionError:
      return numpy.core.nan

  def getDependenceNonUniformityNormalizedFeatureValue(self):
    r"""
    **5. Dependence Non-Uniformity Normalized (DNN)**

    .. math::
      DNN = \frac{\sum^{N_d}_{j=1}\left(\sum^{N_g}_{i=1}{\textbf{P}(i,j)}\right)^2}{N_z^2}

    Measures the similarity of dependence throughout the image, with a lower value indicating
    more homogeneity among dependencies in the image. This is the normalized version of the DLN formula.
    """
    pd = self.coefficients['pd']
    Nz = self.coefficients['Np']  # Nz = Np, see class docstring

    try:
      dnn = numpy.sum(pd ** 2) / (Nz ** 2)
      return dnn
    except ZeroDivisionError:
      return numpy.core.nan

  def getGrayLevelVarianceFeatureValue(self):
    r"""
    **6. Gray Level Variance (GLV)**

    .. math::
      GLV = \displaystyle\sum^{N_g}_{i=1}\displaystyle\sum^{N_d}_{j=1}{p(i,j)(i - \mu)^2} \text{, where}
      \mu = \displaystyle\sum^{N_g}_{i=1}\displaystyle\sum^{N_d}_{j=1}{ip(i,j)}

    Measures the variance in grey level in the image.
    """
    ivector = self.coefficients['ivector']
    Nz = self.coefficients['Np']  # Nz = Np, see class docstring
    pg = self.coefficients['pg'] / Nz  # divide by Nz to get the normalized matrix

    u_i = numpy.sum(pg * ivector)
    glv = numpy.sum(pg * (ivector - u_i) ** 2)
    return glv

  def getDependenceVarianceFeatureValue(self):
    r"""
    **7. Dependence Variance (DV)**

    .. math::
      DV = \displaystyle\sum^{N_g}_{i=1}\displaystyle\sum^{N_d}_{j=1}{p(i,j)(j - \mu)^2} \text{, where}
      \mu = \displaystyle\sum^{N_g}_{i=1}\displaystyle\sum^{N_d}_{j=1}{jp(i,j)}

    Measures the variance in dependence size in the image.
    """
    jvector = self.coefficients['jvector']
    Nz = self.coefficients['Np']  # Nz = Np, see class docstring
    pd = self.coefficients['pd'] / Nz  # divide by Nz to get the normalized matrix

    u_j = numpy.sum(pd * jvector)
    dv = numpy.sum(pd * (jvector - u_j) ** 2)
    return dv

  def getDependenceEntropyFeatureValue(self):
    r"""
    **8. Dependence Entropy (DE)**

    .. math::
      Dependence Entropy = -\displaystyle\sum^{N_g}_{i=1}\displaystyle\sum^{N_d}_{j=1}{p(i,j)\log_{2}(p(i,j)+\epsilon)}
    """
    eps = numpy.spacing(1)
    Nz = self.coefficients['Np']  # Nz = Np, see class docstring
    p_gldm = self.P_gldm / Nz  # divide by Nz to get the normalized matrix

    return -numpy.sum(p_gldm * numpy.log2(p_gldm + eps))

  @deprecated
  def getDependencePercentageFeatureValue(self):
    r"""
    **DEPRECATED. Dependence Percentage**

    .. math::
      \textit{dependence percentage} = \frac{N_z}{N_p}

    .. warning::
      This feature has been deprecated, as it would always compute 1. See
      :ref:`here <radiomics-excluded-gldm-dependence-percentage-label>` for more details. **Enabling this feature will result in the
      logging of a DeprecationWarning (does not interrupt extraction of other features), no value is calculated for this features**
    """
    raise DeprecationWarning('GLDM - Dependence Percentage always computes 1, '
                             'see http://pyradiomics.readthedocs.io/en/latest/removedfeatures.html for more details')

  def getLowGrayLevelEmphasisFeatureValue(self):
    r"""
    **9. Low Gray Level Emphasis (LGLE)**

    .. math::
      LGLE = \frac{\sum^{N_g}_{i=1}\sum^{N_d}_{j=1}{\frac{\textbf{P}(i,j)}{i^2}}}{N_z}

    Measures the distribution of low gray-level values, with a higher value indicating a greater
    concentration of low gray-level values in the image.
    """
    pg = self.coefficients['pg']
    ivector = self.coefficients['ivector']
    Nz = self.coefficients['Np']

    try:
      lgle = numpy.sum(pg / (ivector ** 2)) / Nz
      return lgle
    except ZeroDivisionError:
      return numpy.core.nan

  def getHighGrayLevelEmphasisFeatureValue(self):
    r"""
    **10. High Gray Level Emphasis (HGLE)**

    .. math::
      HGLE = \frac{\sum^{N_g}_{i=1}\sum^{N_d}_{j=1}{\textbf{P}(i,j)i^2}}{N_z}

    Measures the distribution of the higher gray-level values, with a higher value indicating
    a greater concentration of high gray-level values in the image.
    """
    pg = self.coefficients['pg']
    ivector = self.coefficients['ivector']
    Nz = self.coefficients['Np']

    try:
      hgle = numpy.sum(pg * (ivector ** 2)) / Nz
      return hgle
    except ZeroDivisionError:
      return numpy.core.nan

  def getSmallDependenceLowGrayLevelEmphasisFeatureValue(self):
    r"""
    **11. Small Dependence Low Gray Level Emphasis (SDLGLE)**

    .. math::
      SDLGLE = \frac{\sum^{N_g}_{i=1}\sum^{N_d}_{j=1}{\frac{\textbf{P}(i,j)}{i^2j^2}}}{N_z}

    Measures the joint distribution of small dependence with lower gray-level values.
    """
    ivector = self.coefficients['ivector']
    jvector = self.coefficients['jvector']
    Nz = self.coefficients['Np']

    try:
      sdlgle = numpy.sum(self.P_gldm / ((ivector[:, None] ** 2) * (jvector[None, :] ** 2))) / Nz
      return sdlgle
    except ZeroDivisionError:
      return numpy.core.nan

  def getSmallDependenceHighGrayLevelEmphasisFeatureValue(self):
    r"""
    **12. Small Dependence High Gray Level Emphasis (SDHGLE)**

    .. math:
      SDHGLE = \frac{\sum^{N_g}_{i=1}\sum^{N_d}_{j=1}{\frac{\textbf{P}(i,j)i^2}{j^2}}}{N_z}

    Measures the joint distribution of small dependence with higher gray-level values.
    """
    ivector = self.coefficients['ivector']
    jvector = self.coefficients['jvector']
    Nz = self.coefficients['Np']

    try:
      sdhgle = numpy.sum(self.P_gldm * (ivector[:, None] ** 2) / (jvector[None, :] ** 2)) / Nz
      return sdhgle
    except ZeroDivisionError:
      return numpy.core.nan

  def getLargeDependenceLowGrayLevelEmphasisFeatureValue(self):
    r"""
    **13. Large Dependence Low Gray Level Emphasis (LDLGLE)**

    .. math::
      LDLGLE = \frac{\sum^{N_g}_{i=1}\sum^{N_d}_{j=1}{\frac{\textbf{P}(i,j)j^2}{i^2}}}{N_z}

    Measures the joint distribution of large dependence with lower gray-level values.
    """
    ivector = self.coefficients['ivector']
    jvector = self.coefficients['jvector']
    Nz = self.coefficients['Np']

    try:
      ldlgle = numpy.sum(self.P_gldm * (jvector[None, :] ** 2) / (ivector[:, None] ** 2)) / Nz
      return ldlgle
    except ZeroDivisionError:
      return numpy.core.nan

  def getLargeDependenceHighGrayLevelEmphasisFeatureValue(self):
    r"""
    **14. Large Dependence High Gray Level Emphasis (LDHGLE)**

    .. math::
      LDHGLE = \frac{\sum^{N_g}_{i=1}\sum^{N_d}_{j=1}{\textbf{P}(i,j)i^2j^2}}{N_z}

    Measures the joint distribution of large dependence with higher gray-level values.
    """
    ivector = self.coefficients['ivector']
    jvector = self.coefficients['jvector']
    Nz = self.coefficients['Np']

    try:
      ldhgle = numpy.sum(self.P_gldm * ((jvector[None, :] ** 2) * (ivector[:, None] ** 2))) / Nz
      return ldhgle
    except ZeroDivisionError:
      return numpy.core.nan
