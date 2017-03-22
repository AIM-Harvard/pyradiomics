import numpy
from six.moves import range
from tqdm import trange

from radiomics import base, cMatrices, cMatsEnabled, imageoperations


class RadiomicsGLSZM(base.RadiomicsFeaturesBase):
  r"""
  A Gray Level Size Zone (GLSZM) quantifies gray level zones in an image.
  A gray level zone is defined as a the number of connected voxels that share the same
  gray level intensity. A voxel is considered connected if the distance is 1 according to the infinity norm. This
  yields a 26-connected region in a 3D image, and an 8-connected region in a 2D image.
  In a gray level size zone matrix :math:`P(i,j)` the :math:`(i,j)^{\text{th}}` element describes the number of times
  a gray level zone with gray level :math:`i` and size :math:`j` appears in image.

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

  :math:`\textbf{P}(i,j)` be the size zone matrix

  :math:`p(i,j)` be the normalized size zone matrix, defined as :math:`p(i,j) = \frac{\textbf{P}(i,j)}{\sum{\textbf{P}(i,j)}}`

  :math:`N_g` be the number of discreet intensity values in the image

  :math:`N_s` be the number of discreet zone sizes in the image

  :math:`N_p` be the number of voxels in the image

  References

  - Guillaume Thibault; Bernard Fertil; Claire Navarro; Sandrine Pereira; Pierre Cau; Nicolas Levy; Jean Sequeira;
    Jean-Luc Mari (2009). "Texture Indexes and Gray Level Size Zone Matrix. Application to Cell Nuclei Classification".
    Pattern Recognition and Information Processing (PRIP): 140-145.
  - `<https://en.wikipedia.org/wiki/Gray_level_size_zone_matrix>`_
  """

  def __init__(self, inputImage, inputMask, **kwargs):
    super(RadiomicsGLSZM, self).__init__(inputImage, inputMask, **kwargs)

    self.coefficients = {}
    self.P_glszm = {}

    # binning
    self.matrix, self.binEdges = imageoperations.binImage(self.binWidth, self.matrix, self.matrixCoordinates)
    self.coefficients['Ng'] = int(numpy.max(self.matrix[self.matrixCoordinates]))  # max gray level in the ROI
    self.coefficients['Np'] = self.targetVoxelArray.size

    if cMatsEnabled():
      self.P_glszm = self._calculateCMatrix()
    else:
      self.P_glszm = self._calculateMatrix()

    self._calculateCoefficients()

    self.logger.debug('Feature class initialized, calculated GLSZM with shape %s', self.P_glszm.shape)

  def _calculateMatrix(self):
    """
    Number of times a region with a
    gray level and voxel count occurs in an image. P_glszm[level, voxel_count] = # occurrences

    For 3D-images this concerns a 26-connected region, for 2D an 8-connected region
    """
    self.logger.debug("Calculating GLSZM matrix in Python")

    size = numpy.max(self.matrixCoordinates, 1) - numpy.min(self.matrixCoordinates, 1) + 1
    angles = imageoperations.generateAngles(size, **self.kwargs)

    # Empty GLSZ matrix
    P_glszm = numpy.zeros((self.coefficients['Ng'], self.coefficients['Np']))

    # Iterate over all gray levels in the image
    numGrayLevels = self.coefficients['Ng'] + 1

    if self.verbose: bar = trange(numGrayLevels - 1, desc='calculate GLSZM')

    for i in range(1, numGrayLevels):
      # give some progress
      if self.verbose: bar.update()

      ind = zip(*numpy.where(self.matrix == i))
      ind = list(set(ind).intersection(set(zip(*self.matrixCoordinates))))

      while ind:  # check if ind is not empty: unprocessed regions for current gray level
        # Pop first coordinate of an unprocessed zone, start new stack
        ind_region = [ind.pop()]

        # Define regionSize
        regionSize = 0

        # Grow zone for item popped from stack of region indices, loop until stack of region indices is exhausted
        # Each loop represents one voxel belonging to current zone. Therefore, count number of loops as regionSize
        while ind_region:
          regionSize += 1

          # Use pop to remove next node for set of unprocessed region indices
          ind_node = ind_region.pop()

          # get all coordinates in the 26-connected region, 2 voxels per angle
          region_full = [tuple(sum(a) for a in zip(ind_node, angle_i)) for angle_i in angles]
          region_full += [tuple(sum(a) for a in zip(ind_node, angle_i)) for angle_i in angles * -1]

          # get all unprocessed coordinates in the 26-connected region with same gray level
          region_level = list(set(ind).intersection(set(region_full)))

          # Remove already processed indices to prevent reprocessing
          ind = list(set(ind) - set(region_level))

          # Add all found neighbours to the total stack of unprocessed neighbours
          ind_region.extend(region_level)

        # Update the gray level size zone matrix
        P_glszm[i - 1, regionSize - 1] += 1

    if self.verbose: bar.close()

    # Crop gray-level axis of GLSZM matrix to between minimum and maximum observed gray-levels
    # Crop size-zone area axis of GLSZM matrix up to maximum observed size-zone area
    self.logger.debug('Cropping calculated matrix to observed gray levels and maximum observed zone size')
    P_glszm_bounds = numpy.argwhere(P_glszm)
    (xstart, ystart), (xstop, ystop) = P_glszm_bounds.min(0), P_glszm_bounds.max(0) + 1  # noqa: F841
    return P_glszm[xstart:xstop, :ystop]

  def _calculateCMatrix(self):
    self.logger.debug("Calculating GLSZM matrix in C")

    size = numpy.max(self.matrixCoordinates, 1) - numpy.min(self.matrixCoordinates, 1) + 1
    angles = imageoperations.generateAngles(size, **self.kwargs)
    Ng = self.coefficients['Ng']
    Ns = self.coefficients['Np']

    return cMatrices.calculate_glszm(self.matrix, self.maskArray, angles, Ng, Ns)

  def _calculateCoefficients(self):
    self.logger.debug("Calculating GLSZM coefficients")

    sumP_glszm = numpy.sum(self.P_glszm, (0, 1))

    # set sum to numpy.spacing(1) if sum is 0?
    if sumP_glszm == 0:
      sumP_glszm = 1

    pr = numpy.sum(self.P_glszm, 0)
    pg = numpy.sum(self.P_glszm, 1)

    ivector = numpy.arange(1, self.P_glszm.shape[0] + 1, dtype=numpy.float64)
    jvector = numpy.arange(1, self.P_glszm.shape[1] + 1, dtype=numpy.float64)

    self.coefficients['sumP_glszm'] = sumP_glszm
    self.coefficients['pr'] = pr
    self.coefficients['pg'] = pg
    self.coefficients['ivector'] = ivector
    self.coefficients['jvector'] = jvector

  def getSmallAreaEmphasisFeatureValue(self):
    r"""
    Calculate and return the Small Area Emphasis (SAE) value.

    :math:`SAE = \frac{\sum^{N_g}_{i=1}\sum^{N_s}_{j=1}{\frac{\textbf{P}(i,j)}{j^2}}}{\sum^{N_g}_{i=1}\sum^{N_s}_{j=1}{\textbf{P}(i,j)}}`

    A measure of the distribution of small size zones, with a greater value indicative
    of more smaller size zones and more fine textures.
    """
    try:
      sae = numpy.sum(self.coefficients['pr'] / (self.coefficients['jvector'] ** 2)) / self.coefficients['sumP_glszm']
    except ZeroDivisionError:
      sae = numpy.core.numeric.NaN
    return (sae)

  def getLargeAreaEmphasisFeatureValue(self):
    r"""
    Calculate and return the Large Area Emphasis (LAE) value.

    :math:`LAE = \frac{\sum^{N_g}_{i=1}\sum^{N_s}_{j=1}{\textbf{P}(i,j)j^2}}{\sum^{N_g}_{i=1}\sum^{N_s}_{j=1}{\textbf{P}(i,j)}}`

    A measure of the distribution of large area size zones, with a greater value indicative
    of more larger size zones and more coarse textures.
    """
    try:
      lae = numpy.sum(self.coefficients['pr'] * (self.coefficients['jvector'] ** 2)) / self.coefficients['sumP_glszm']
    except ZeroDivisionError:
      lae = numpy.core.numeric.NaN
    return (lae)

  def getGrayLevelNonUniformityFeatureValue(self):
    r"""
    Calculate and return the Gray Level Non-Uniformity (GLN) value.

    :math:`GLN = \frac{\sum^{N_g}_{i=1}\left(\sum^{N_s}_{j=1}{\textbf{P}(i,j)}\right)^2}{\sum^{N_g}_{i=1}\sum^{N_s}_{j=1}{\textbf{P}(i,j)}}`

    Measures the variability of gray-level intensity values in the image, with a lower value indicating
    more homogeneity in intensity values.
    """
    try:
      iv = numpy.sum(self.coefficients['pg'] ** 2) / self.coefficients['sumP_glszm']
    except ZeroDivisionError:
      iv = numpy.core.numeric.NaN
    return (iv)

  def getGrayLevelNonUniformityNormalizedFeatureValue(self):
    r"""
    Calculate and return the Gray Level Non-Uniformity Normalized (GLNN) value.

    :math:`GLNN = \frac{\sum^{N_g}_{i=1}\left(\sum^{N_s}_{j=1}{\textbf{P}(i,j)}\right)^2}{\sum^{N_g}_{i=1}\sum^{N_d}_{j=1}{\textbf{P}(i,j)}^2}`

    Measures the variability of gray-level intensity values in the image, with a lower value indicating
    a greater similarity in intensity values. This is the normalized version of the GLN formula.
    """
    try:
      ivn = numpy.sum(self.coefficients['pg'] ** 2) / self.coefficients['sumP_glszm'] ** 2
    except ZeroDivisionError:
      ivn = numpy.core.numeric.NaN
    return (ivn)

  def getSizeZoneNonUniformityFeatureValue(self):
    r"""
    Calculate and return the Size-Zone Non-Uniformity (SZN) value.

    :math:`SZN = \frac{\sum^{N_s}_{j=1}\left(\sum^{N_g}_{i=1}{\textbf{P}(i,j)}\right)^2}{\sum^{N_g}_{i=1}\sum^{N_s}_{j=1}{\textbf{P}(i,j)}}`

    Measures the variability of size zone volumes in the image, with a lower value indicating
    more homogeneity in size zone volumes.
    """
    try:
      szv = numpy.sum(self.coefficients['pr'] ** 2) / self.coefficients['sumP_glszm']
    except ZeroDivisionError:
      szv = numpy.core.numeric.NaN
    return (szv)

  def getSizeZoneNonUniformityNormalizedFeatureValue(self):
    r"""
    Calculate and return the Size-Zone Non-Uniformity Normalized (SZNN) value.

    :math:`SZNN = \frac{\sum^{N_s}_{j=1}\left(\sum^{N_g}_{i=1}{\textbf{P}(i,j)}\right)^2}{\sum^{N_g}_{i=1}\sum^{N_d}_{j=1}{\textbf{P}(i,j)}^2}`

    Measures the variability of size zone volumes throughout the image, with a lower value indicating
    more homogeneity among zone size volumes in the image. This is the normalized version of the SZN formula.
    """
    try:
      szvn = numpy.sum(self.coefficients['pr'] ** 2) / self.coefficients['sumP_glszm'] ** 2
    except ZeroDivisionError:
      szvn = numpy.core.numeric.NaN
    return (szvn)

  def getZonePercentageFeatureValue(self):
    r"""
    Calculate and return the Zone Percentage (ZP) value.

    :math:`ZP = \sum^{N_g}_{i=1}\sum^{N_s}_{j=1}{\frac{\textbf{P}(i,j)}{N_p}}`

    Measures the coarseness of the texture by taking the ratio of number of zones and number of voxels in the ROI.
    Values are in range :math:`\frac{1}{N_p} \leq ZP \leq 1`, with higher values indicating a larger portion of the ROI
    consists of small zones (indicates a more fine texture).
    """
    try:
      zp = self.coefficients['sumP_glszm'] / self.coefficients['Np']
    except ZeroDivisionError:
      zp = numpy.core.numeric.NaN
    return (zp)

  def getGrayLevelVarianceFeatureValue(self):
    r"""
    Calculate and return the Gray Level Variance (GLV) value.

    :math:`GLV = \displaystyle\sum^{N_g}_{i=1}\displaystyle\sum^{N_s}_{j=1}{p(i,j)(i - \mu)^2}`, where

    :math:`\mu = \displaystyle\sum^{N_g}_{i=1}\displaystyle\sum^{N_s}_{j=1}{p(i,j)i}`

    Measures the variance in gray level intensities for the zones.
    """
    ivector = self.coefficients['ivector']
    sumP_glszm = self.coefficients['sumP_glszm']
    u_i = numpy.sum(self.coefficients['pg'] * ivector) / sumP_glszm
    glv = numpy.sum(self.coefficients['pg'] * (ivector - u_i) ** 2) / sumP_glszm
    return glv

  def getZoneVarianceFeatureValue(self):
    r"""
    Calculate and return the Zone Variance (ZV) value.

    :math:`ZV = \displaystyle\sum^{N_g}_{i=1}\displaystyle\sum^{N_s}_{j=1}{p(i,j)(j - \mu)^2}`, where

    :math:`\mu = \displaystyle\sum^{N_g}_{i=1}\displaystyle\sum^{N_s}_{j=1}{p(i,j)j}`

    Measures the variance in zone size volumes for the zones.
    """
    jvector = self.coefficients['jvector']
    sumP_glszm = self.coefficients['sumP_glszm']
    u_j = numpy.sum(self.coefficients['pr'] * jvector) / sumP_glszm
    zv = numpy.sum(self.coefficients['pr'] * (jvector - u_j) ** 2) / sumP_glszm
    return zv

  def getZoneEntropyFeatureValue(self):
    r"""
    Calculate and return the Zone Entropy (ZE) value.

    :math:`ZE = -\displaystyle\sum^{N_g}_{i=1}\displaystyle\sum^{N_s}_{j=1}{p(i,j)\log_{2}(p(i,j)+\epsilon)}`

    Here, :math:`\epsilon` is an arbitrarily small positive number (:math:`\approx 2.2\times10^{-16}`).

    Measures the uncertainty/randomness in the distribution of zone sizes and gray levels. A higher value indicates more
    heterogeneneity in the texture patterns.
    """
    eps = numpy.spacing(1)
    sumP_glszm = self.coefficients['sumP_glszm']
    p_glszm = self.P_glszm / sumP_glszm
    return -numpy.sum(p_glszm * numpy.log2(p_glszm + eps))

  def getLowGrayLevelZoneEmphasisFeatureValue(self):
    r"""
    Calculate and return the Low Gray Level Zone Emphasis (LGLZE) value.

    :math:`LGLZE = \frac{\sum^{N_g}_{i=1}\sum^{N_s}_{j=1}{\frac{\textbf{P}(i,j)}{i^2}}}{\sum^{N_g}_{i=1}\sum^{N_s}_{j=1}{\textbf{P}(i,j)}}`

    Measures the distribution of lower gray-level size zones, with a higher value indicating a greater proportion of
    lower gray-level values and size zones in the image.
    """
    try:
      lie = numpy.sum((self.coefficients['pg'] / (self.coefficients['ivector'] ** 2))) / self.coefficients['sumP_glszm']
    except ZeroDivisionError:
      lie = numpy.core.numeric.NaN
    return (lie)

  def getHighGrayLevelZoneEmphasisFeatureValue(self):
    r"""
    Calculate and return the High Gray Level Zone Emphasis (HGLZE) value.

    :math:`HGLZE = \frac{\sum^{N_g}_{i=1}\sum^{N_s}_{j=1}{\textbf{P}(i,j)i^2}}{\sum^{N_g}_{i=1}\sum^{N_s}_{j=1}{\textbf{P}(i,j)}}`

    Measures the distribution of the higher gray-level values, with a higher value indicating a greater proportion of
    higher gray-level values and size zones in the image.
    """
    try:
      hie = numpy.sum((self.coefficients['pg'] * (self.coefficients['ivector'] ** 2))) / self.coefficients['sumP_glszm']
    except ZeroDivisionError:
      hie = numpy.core.numeric.NaN
    return (hie)

  def getSmallAreaLowGrayLevelEmphasisFeatureValue(self):
    r"""
    Calculate and return the Small Area Low Gray Level Emphasis (SALGLE) value.

    :math:`SALGLE = \frac{\sum^{N_g}_{i=1}\sum^{N_s}_{j=1}{\frac{\textbf{P}(i,j)}{i^2j^2}}}{\sum^{N_g}_{i=1}\sum^{N_s}_{j=1}{\textbf{P}(i,j)}}`

    Measures the proportion in the image of the joint distribution of smaller size zones with lower gray-level values.
    """
    try:
      lisae = numpy.sum(
        (self.P_glszm / ((self.coefficients['ivector'][:, None] ** 2) * (self.coefficients['jvector'][None, :] ** 2))),
        (0, 1)) / self.coefficients['sumP_glszm']
    except ZeroDivisionError:
      lisae = numpy.core.numeric.NaN
    return (lisae)

  def getSmallAreaHighGrayLevelEmphasisFeatureValue(self):
    r"""
    Calculate and return the Small Area High Gray Level Emphasis (SAHGLE) value.

    :math:`SAHGLE = \frac{\sum^{N_g}_{i=1}\sum^{N_s}_{j=1}{\frac{\textbf{P}(i,j)i^2}{j^2}}}{\sum^{N_g}_{i=1}\sum^{N_s}_{j=1}{\textbf{P}(i,j)}}`

    Measures the proportion in the image of the joint distribution of smaller size zones with higher gray-level values.
    """
    try:
      hisae = numpy.sum(
        (self.P_glszm * (self.coefficients['ivector'][:, None] ** 2) / (self.coefficients['jvector'][None, :] ** 2)),
        (0, 1)) / self.coefficients['sumP_glszm']
    except ZeroDivisionError:
      hisae = numpy.core.numeric.NaN
    return (hisae)

  def getLargeAreaLowGrayLevelEmphasisFeatureValue(self):
    r"""
    Calculate and return the Large Area Low Gray Level Emphasis (LALGLE) value.

    :math:`LALGLE = \frac{\sum^{N_g}_{i=1}\sum^{N_s}_{j=1}{\frac{\textbf{P}(i,j)j^2}{i^2}}}{\sum^{N_g}_{i=1}\sum^{N_s}_{j=1}{\textbf{P}(i,j)}}`

    Measures the proportion in the image of the joint distribution of larger size zones with lower gray-level values.
    """
    try:
      lilae = numpy.sum(
        (self.P_glszm * (self.coefficients['jvector'][None, :] ** 2) / (self.coefficients['ivector'][:, None] ** 2)),
        (0, 1)) / self.coefficients['sumP_glszm']
    except ZeroDivisionError:
      lilae = numpy.core.numeric.NaN
    return (lilae)

  def getLargeAreaHighGrayLevelEmphasisFeatureValue(self):
    r"""
    Calculate and return the Large Area High Gray Level Emphasis (LAHGLE) value.

    :math:`LAHGLE = \frac{\sum^{N_g}_{i=1}\sum^{N_s}_{j=1}{\textbf{P}(i,j)i^2j^2}}{\sum^{N_g}_{i=1}\sum^{N_s}_{j=1}{\textbf{P}(i,j)}}`

    Measures the proportion in the image of the joint distribution of larger size zones with higher gray-level values.
    """
    try:
      hilae = numpy.sum(
        (self.P_glszm * ((self.coefficients['jvector'][None, :] ** 2) * (self.coefficients['ivector'][:, None] ** 2))),
        (0, 1)) / self.coefficients['sumP_glszm']
    except ZeroDivisionError:
      hilae = numpy.core.numeric.NaN
    return (hilae)
