import numpy
from six.moves import range

from radiomics import base, cMatrices


class RadiomicsGLSZM(base.RadiomicsFeaturesBase):
  r"""
  A Gray Level Size Zone (GLSZM) quantifies gray level zones in an image. A gray level zone is defined as a the number
  of connected voxels that share the same gray level intensity. A voxel is considered connected if the distance is 1
  according to the infinity norm (26-connected region in a 3D, 8-connected region in 2D).
  In a gray level size zone matrix :math:`P(i,j)` the :math:`(i,j)^{\text{th}}` element equals the number of zones
  with gray level :math:`i` and size :math:`j` appear in image. Contrary to GLCM and GLRLM, the GLSZM is rotation
  independent, with only one matrix calculated for all directions in the ROI.

  As a two dimensional example, consider the following 5x5 image, with 5 discrete gray levels:

  .. math::
    \textbf{I} = \begin{bmatrix}
    5 & 2 & 5 & 4 & 4\\
    3 & 3 & 3 & 1 & 3\\
    2 & 1 & 1 & 1 & 3\\
    4 & 2 & 2 & 2 & 3\\
    3 & 5 & 3 & 3 & 2 \end{bmatrix}

  The GLSZM then becomes:

  .. math::
    \textbf{P} = \begin{bmatrix}
    0 & 0 & 0 & 1 & 0\\
    1 & 0 & 0 & 0 & 1\\
    1 & 0 & 1 & 0 & 1\\
    1 & 1 & 0 & 0 & 0\\
    3 & 0 & 0 & 0 & 0 \end{bmatrix}

  Let:

  - :math:`N_g` be the number of discreet intensity values in the image
  - :math:`N_s` be the number of discreet zone sizes in the image
  - :math:`N_p` be the number of voxels in the image
  - :math:`N_z` be the number of zones in the ROI, which is equal to :math:`\sum^{N_g}_{i=1}\sum^{N_s}_{j=1}
    {\textbf{P}(i,j)}` and :math:`1 \leq N_z \leq N_p`
  - :math:`\textbf{P}(i,j)` be the size zone matrix
  - :math:`p(i,j)` be the normalized size zone matrix, defined as :math:`p(i,j) = \frac{\textbf{P}(i,j)}{N_z}`

  .. note::
    The mathematical formulas that define the GLSZM features correspond to the definitions of features extracted from
    the GLRLM.

  References

  - Guillaume Thibault; Bernard Fertil; Claire Navarro; Sandrine Pereira; Pierre Cau; Nicolas Levy; Jean Sequeira;
    Jean-Luc Mari (2009). "Texture Indexes and Gray Level Size Zone Matrix. Application to Cell Nuclei Classification".
    Pattern Recognition and Information Processing (PRIP): 140-145.
  - `<https://en.wikipedia.org/wiki/Gray_level_size_zone_matrix>`_
  """

  def __init__(self, inputImage, inputMask, **kwargs):
    super(RadiomicsGLSZM, self).__init__(inputImage, inputMask, **kwargs)

    self.P_glszm = None
    self._applyBinning()

  def _initCalculation(self):
    self.coefficients['Np'] = len(self.labelledVoxelCoordinates[0])

    self.P_glszm = self._calculateMatrix()

    self._calculateCoefficients()

    self.logger.debug('GLSZM feature class initialized, calculated GLSZM with shape %s', self.P_glszm.shape)

  def _calculateMatrix(self):
    """
    Number of times a region with a
    gray level and voxel count occurs in an image. P_glszm[level, voxel_count] = # occurrences

    For 3D-images this concerns a 26-connected region, for 2D an 8-connected region
    """
    self.logger.debug('Calculating GLSZM matrix in C')
    Ng = self.coefficients['Ng']
    Ns = self.coefficients['Np']

    P_glszm = cMatrices.calculate_glszm(self.matrix,
                                        self.maskArray,
                                        Ng,
                                        Ns,
                                        self.settings.get('force2D', False),
                                        self.settings.get('force2Ddimension', 0))

    # Delete rows that specify gray levels not present in the ROI
    NgVector = range(1, Ng + 1)  # All possible gray values
    GrayLevels = self.coefficients['grayLevels']  # Gray values present in ROI
    emptyGrayLevels = numpy.array(list(set(NgVector) - set(GrayLevels)))  # Gray values NOT present in ROI

    P_glszm = numpy.delete(P_glszm, emptyGrayLevels - 1, 0)

    return P_glszm

  def _calculateCoefficients(self):
    self.logger.debug('Calculating GLSZM coefficients')

    Nz = numpy.sum(self.P_glszm, (0, 1))

    # set sum to numpy.spacing(1) if sum is 0?
    if Nz == 0:
      Nz = 1

    ps = numpy.sum(self.P_glszm, 0)
    pg = numpy.sum(self.P_glszm, 1)

    ivector = self.coefficients['grayLevels']
    jvector = numpy.arange(1, self.P_glszm.shape[1] + 1, dtype=numpy.float64)

    # Delete columns that specify zone sizes not present in the ROI
    emptyZoneSizes = numpy.where(ps == 0)
    self.P_glszm = numpy.delete(self.P_glszm, emptyZoneSizes, 1)
    jvector = numpy.delete(jvector, emptyZoneSizes)
    ps = numpy.delete(ps, emptyZoneSizes)

    self.coefficients['Nz'] = Nz
    self.coefficients['ps'] = ps
    self.coefficients['pg'] = pg
    self.coefficients['ivector'] = ivector
    self.coefficients['jvector'] = jvector

  def getSmallAreaEmphasisFeatureValue(self):
    r"""
    **1. Small Area Emphasis (SAE)**

    .. math::
      \textit{SAE} = \frac{\sum^{N_g}_{i=1}\sum^{N_s}_{j=1}{\frac{\textbf{P}(i,j)}{j^2}}}{N_z}

    SAE is a measure of the distribution of small size zones, with a greater value indicative of more smaller size zones
    and more fine textures.
    """
    ps = self.coefficients['ps']
    jvector = self.coefficients['jvector']
    Nz = self.coefficients['Nz']

    sae = numpy.sum(ps / (jvector ** 2)) / Nz
    return sae

  def getLargeAreaEmphasisFeatureValue(self):
    r"""
    **2. Large Area Emphasis (LAE)**

    .. math::
      \textit{LAE} = \frac{\sum^{N_g}_{i=1}\sum^{N_s}_{j=1}{\textbf{P}(i,j)j^2}}{N_z}

    LAE is a measure of the distribution of large area size zones, with a greater value indicative of more larger size
    zones and more coarse textures.
    """
    ps = self.coefficients['ps']
    jvector = self.coefficients['jvector']
    Nz = self.coefficients['Nz']

    lae = numpy.sum(ps * (jvector ** 2)) / Nz
    return lae

  def getGrayLevelNonUniformityFeatureValue(self):
    r"""
    **3. Gray Level Non-Uniformity (GLN)**

    .. math::
      \textit{GLN} = \frac{\sum^{N_g}_{i=1}\left(\sum^{N_s}_{j=1}{\textbf{P}(i,j)}\right)^2}{N_z}

    GLN measures the variability of gray-level intensity values in the image, with a lower value indicating more
    homogeneity in intensity values.
    """
    pg = self.coefficients['pg']
    Nz = self.coefficients['Nz']

    iv = numpy.sum(pg ** 2) / Nz
    return iv

  def getGrayLevelNonUniformityNormalizedFeatureValue(self):
    r"""
    **4. Gray Level Non-Uniformity Normalized (GLNN)**

    .. math::
      \textit{GLNN} = \frac{\sum^{N_g}_{i=1}\left(\sum^{N_s}_{j=1}{\textbf{P}(i,j)}\right)^2}{N_z^2}

    GLNN measures the variability of gray-level intensity values in the image, with a lower value indicating a greater
    similarity in intensity values. This is the normalized version of the GLN formula.
    """
    pg = self.coefficients['pg']
    Nz = self.coefficients['Nz']

    ivn = numpy.sum(pg ** 2) / Nz ** 2
    return ivn

  def getSizeZoneNonUniformityFeatureValue(self):
    r"""
    **5. Size-Zone Non-Uniformity (SZN)**

    .. math::
      \textit{SZN} = \frac{\sum^{N_s}_{j=1}\left(\sum^{N_g}_{i=1}{\textbf{P}(i,j)}\right)^2}{N_z}

    SZN measures the variability of size zone volumes in the image, with a lower value indicating more homogeneity in
    size zone volumes.
    """
    ps = self.coefficients['ps']
    Nz = self.coefficients['Nz']

    szv = numpy.sum(ps ** 2) / Nz
    return szv

  def getSizeZoneNonUniformityNormalizedFeatureValue(self):
    r"""
    **6. Size-Zone Non-Uniformity Normalized (SZNN)**

    .. math::
      \textit{SZNN} = \frac{\sum^{N_s}_{j=1}\left(\sum^{N_g}_{i=1}{\textbf{P}(i,j)}\right)^2}{N_z^2}

    SZNN measures the variability of size zone volumes throughout the image, with a lower value indicating more
    homogeneity among zone size volumes in the image. This is the normalized version of the SZN formula.
    """
    ps = self.coefficients['ps']
    Nz = self.coefficients['Nz']

    szvn = numpy.sum(ps ** 2) / Nz ** 2
    return szvn

  def getZonePercentageFeatureValue(self):
    r"""
    **7. Zone Percentage (ZP)**

    .. math::
      \textit{ZP} = \frac{N_z}{N_p}

    ZP measures the coarseness of the texture by taking the ratio of number of zones and number of voxels in the ROI.

    Values are in range :math:`\frac{1}{N_p} \leq ZP \leq 1`, with higher values indicating a larger portion of the ROI
    consists of small zones (indicates a more fine texture).
    """
    Nz = self.coefficients['Nz']
    Np = self.coefficients['Np']

    zp = Nz / Np
    return zp

  def getGrayLevelVarianceFeatureValue(self):
    r"""
    **8. Gray Level Variance (GLV)**

    .. math::
      \textit{GLV} = \displaystyle\sum^{N_g}_{i=1}\displaystyle\sum^{N_s}_{j=1}{p(i,j)(i - \mu)^2}

    Here, :math:`\mu = \displaystyle\sum^{N_g}_{i=1}\displaystyle\sum^{N_s}_{j=1}{p(i,j)i}`

    GLV measures the variance in gray level intensities for the zones.
    """
    ivector = self.coefficients['ivector']
    Nz = self.coefficients['Nz']
    pg = self.coefficients['pg'] / Nz  # divide by Nz to get the normalized matrix

    u_i = numpy.sum(pg * ivector)
    glv = numpy.sum(pg * (ivector - u_i) ** 2)
    return glv

  def getZoneVarianceFeatureValue(self):
    r"""
    **9. Zone Variance (ZV)**

    .. math::
      \textit{ZV} = \displaystyle\sum^{N_g}_{i=1}\displaystyle\sum^{N_s}_{j=1}{p(i,j)(j - \mu)^2}

    Here, :math:`\mu = \displaystyle\sum^{N_g}_{i=1}\displaystyle\sum^{N_s}_{j=1}{p(i,j)j}`

    ZV measures the variance in zone size volumes for the zones.
    """
    jvector = self.coefficients['jvector']
    Nz = self.coefficients['Nz']
    ps = self.coefficients['ps'] / Nz  # divide by Nz to get the normalized matrix

    u_j = numpy.sum(ps * jvector)
    zv = numpy.sum(ps * (jvector - u_j) ** 2)
    return zv

  def getZoneEntropyFeatureValue(self):
    r"""
    **10. Zone Entropy (ZE)**

    .. math::
      \textit{ZE} = -\displaystyle\sum^{N_g}_{i=1}\displaystyle\sum^{N_s}_{j=1}{p(i,j)\log_{2}(p(i,j)+\epsilon)}

    Here, :math:`\epsilon` is an arbitrarily small positive number (:math:`\approx 2.2\times10^{-16}`).

    ZE measures the uncertainty/randomness in the distribution of zone sizes and gray levels. A higher value indicates
    more heterogeneneity in the texture patterns.
    """
    eps = numpy.spacing(1)
    Nz = self.coefficients['Nz']
    p_glszm = self.P_glszm / Nz  # divide by Nz to get the normalized matrix

    ze = -numpy.sum(p_glszm * numpy.log2(p_glszm + eps))
    return ze

  def getLowGrayLevelZoneEmphasisFeatureValue(self):
    r"""
    **11. Low Gray Level Zone Emphasis (LGLZE)**

    .. math::
      \textit{LGLZE} = \frac{\sum^{N_g}_{i=1}\sum^{N_s}_{j=1}{\frac{\textbf{P}(i,j)}{i^2}}}{N_z}

    LGLZE measures the distribution of lower gray-level size zones, with a higher value indicating a greater proportion
    of lower gray-level values and size zones in the image.
    """
    pg = self.coefficients['pg']
    ivector = self.coefficients['ivector']
    Nz = self.coefficients['Nz']

    lie = numpy.sum(pg / (ivector ** 2)) / Nz
    return lie

  def getHighGrayLevelZoneEmphasisFeatureValue(self):
    r"""
    **12. High Gray Level Zone Emphasis (HGLZE)**

    .. math::
      \textit{HGLZE} = \frac{\sum^{N_g}_{i=1}\sum^{N_s}_{j=1}{\textbf{P}(i,j)i^2}}{N_z}

    HGLZE measures the distribution of the higher gray-level values, with a higher value indicating a greater proportion
    of higher gray-level values and size zones in the image.
    """
    pg = self.coefficients['pg']
    ivector = self.coefficients['ivector']
    Nz = self.coefficients['Nz']

    hie = numpy.sum(pg * (ivector ** 2)) / Nz
    return hie

  def getSmallAreaLowGrayLevelEmphasisFeatureValue(self):
    r"""
    **13. Small Area Low Gray Level Emphasis (SALGLE)**

    .. math::
      \textit{SALGLE} = \frac{\sum^{N_g}_{i=1}\sum^{N_s}_{j=1}{\frac{\textbf{P}(i,j)}{i^2j^2}}}{N_z}

    SALGLE measures the proportion in the image of the joint distribution of smaller size zones with lower gray-level
    values.
    """
    ivector = self.coefficients['ivector']
    jvector = self.coefficients['jvector']
    Nz = self.coefficients['Nz']

    lisae = numpy.sum(self.P_glszm / ((ivector[:, None] ** 2) * (jvector[None, :] ** 2))) / Nz
    return lisae

  def getSmallAreaHighGrayLevelEmphasisFeatureValue(self):
    r"""
    **14. Small Area High Gray Level Emphasis (SAHGLE)**

    .. math::
      \textit{SAHGLE} = \frac{\sum^{N_g}_{i=1}\sum^{N_s}_{j=1}{\frac{\textbf{P}(i,j)i^2}{j^2}}}{N_z}

    SAHGLE measures the proportion in the image of the joint distribution of smaller size zones with higher gray-level
    values.
    """
    ivector = self.coefficients['ivector']
    jvector = self.coefficients['jvector']
    Nz = self.coefficients['Nz']

    hisae = numpy.sum(self.P_glszm * (ivector[:, None] ** 2) / (jvector[None, :] ** 2)) / Nz
    return hisae

  def getLargeAreaLowGrayLevelEmphasisFeatureValue(self):
    r"""
    **15. Large Area Low Gray Level Emphasis (LALGLE)**

    .. math::
      \textit{LALGLE} = \frac{\sum^{N_g}_{i=1}\sum^{N_s}_{j=1}{\frac{\textbf{P}(i,j)j^2}{i^2}}}{N_z}

    LALGLE measures the proportion in the image of the joint distribution of larger size zones with lower gray-level
    values.
    """
    ivector = self.coefficients['ivector']
    jvector = self.coefficients['jvector']
    Nz = self.coefficients['Nz']

    lilae = numpy.sum(self.P_glszm * (jvector[None, :] ** 2) / (ivector[:, None] ** 2)) / Nz
    return lilae

  def getLargeAreaHighGrayLevelEmphasisFeatureValue(self):
    r"""
    **16. Large Area High Gray Level Emphasis (LAHGLE)**

    .. math::
      \textit{LAHGLE} = \frac{\sum^{N_g}_{i=1}\sum^{N_s}_{j=1}{\textbf{P}(i,j)i^2j^2}}{N_z}

    LAHGLE measures the proportion in the image of the joint distribution of larger size zones with higher gray-level
    values.
    """
    ivector = self.coefficients['ivector']
    jvector = self.coefficients['jvector']
    Nz = self.coefficients['Nz']

    hilae = numpy.sum(self.P_glszm * (ivector[:, None] ** 2) * (jvector[None, :] ** 2)) / Nz
    return hilae
