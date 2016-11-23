from itertools import chain
import numpy
import SimpleITK as sitk
from radiomics import base, imageoperations
from tqdm import trange

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

  :math:`\textbf{P}(i,j)` be the dependence matrix

  :math:`p(i,j)` be the normalized dependence matrix, defined as :math:`p(i,j) = \frac{\textbf{P}(i,j)}{\sum{\textbf{P}(i,j)}}`

  :math:`N_g` be the number of discreet intensity values in the image

  :math:`N_d` be the number of discreet dependency sizes in the image

  :math:`N_p` be the number of voxels in the image

  The following class specific settings are possible:

  - gldm_a [0]: float, :math:`\alpha` cutoff value for dependence. A neighbouring voxel with gray level :math:`j` is considered
    dependent on center voxel with gray level :math:`i` if :math:`|i-j|\le\alpha`
  """

  def __init__(self, inputImage, inputMask, **kwargs):
    super(RadiomicsGLDM,self).__init__(inputImage, inputMask, **kwargs)

    self.gldm_a = kwargs.get('gldm_a', 0)

    self.coefficients = {}
    self.P_gldm = {}

    # binning
    self.matrix, self.histogram = imageoperations.binImage(self.binWidth, self.targetVoxelArray, self.matrix, self.matrixCoordinates)
    self.coefficients['Ng'] = self.histogram[1].shape[0] - 1

    self._calculateGLDM()
    self._calculateCoefficients()

  def _calculateGLDM(self):
    Ng = self.coefficients['Ng']

    self.matrix = self.matrix.astype('float')

    # Set voxels outside delineation to padding value
    padVal = numpy.nan
    self.matrix[(self.maskArray != self.label)] = padVal

    size = numpy.max(self.matrixCoordinates, 1) - numpy.min(self.matrixCoordinates, 1) + 1
    angles = imageoperations.generateAngles(size)
    angles = numpy.concatenate((angles, angles * -1))

    depMat = numpy.zeros(self.matrix.shape, dtype='int')

    if self.verbose: bar = trange(len(angles), desc='Calculate shifted matrices (GLDM)')
    for a_idx, a in enumerate(angles):
      if self.verbose: bar.update()
      # create shifted array (by angle), so that for an index idx, angMat[idx] is the neigbour of self.matrix[idx]
      # for the current angle.
      angMat = numpy.roll(numpy.roll(numpy.roll(self.matrix, -a[0], 0), -a[1], 1), -a[2], 2) - self.matrix
      if a[0] > 0:
        angMat[-a[0]:, :, :] = padVal
      elif a[0] < 0:
        angMat[:-a[0], :, :] = padVal

      if a[1] > 0:
        angMat[:, -a[1]:, :] = padVal
      elif a[1] < 0:
        angMat[:, :-a[1], :] = padVal

      if a[2] > 0:
        angMat[:, :, -a[2]:] = padVal
      elif a[2] < 0:
        angMat[:, :, :-a[2]] = padVal

      nanMask = numpy.isnan(angMat)
      depMat[~nanMask] += (numpy.abs(angMat[~nanMask]) <= self.gldm_a)

    if self.verbose: bar.close()

    Nd = numpy.max(depMat)
    P_gldm = numpy.zeros((Ng, Nd + 1))
    grayLevels = numpy.unique(self.matrix[self.matrixCoordinates])

    if self.verbose: bar = trange(len(grayLevels), desc= 'calculate GLDM')
    for i in grayLevels:
        if self.verbose: bar.update()
        i_mat = (self.matrix == i)
        for d in numpy.unique(depMat[i_mat]):
            # By multiplying i_mat and depMat == d, a boolean area is obtained,
            # where the number of elements that are true (1) is equal to the number of voxels
            # with gray level i and dependence d.
            P_gldm[i-1, d] = numpy.sum(i_mat * (depMat == d))
    if self.verbose: bar.close()

    sumP_gldm = numpy.sum(self.P_gldm, (0, 1))
    self.coefficients['sumP_gldm'] = sumP_gldm

  def _calculateCoefficients(self):

    pd = numpy.sum(self.P_gldm, 0)
    pg = numpy.sum(self.P_gldm, 1)

    ivector = numpy.arange(1, self.P_gldm.shape[0] + 1, dtype=numpy.float64)
    jvector = numpy.arange(1, self.P_gldm.shape[1] + 1, dtype=numpy.float64)

    self.coefficients['pd'] = pd
    self.coefficients['pg'] = pg
    self.coefficients['ivector'] = ivector
    self.coefficients['jvector'] = jvector

  def getSmallDependenceEmphasisFeatureValue(self):
    r"""
    Calculate and return the Small Dependence Emphasis (SDE) value.

    :math:`SDE = \frac{\sum^{N_g}_{i=1}\sum^{N_d}_{j=1}{\frac{\textbf{P}(i,j)}{i^2}}}{\sum^{N_g}_{i=1}\sum^{N_d}_{j=1}{\textbf{P}(i,j)}}`

    A measure of the distribution of small dependencies, with a greater value indicative
    of smaller dependence and less homogeneous textures.
    """
    pd = self.coefficients['pd']
    jvector = self.coefficients['jvector']
    sumP_gldm = self.coefficients['sumP_gldm']

    try:
      sde = numpy.sum(pd/(jvector**2)) / sumP_gldm
      return sde
    except ZeroDivisionError:
      return numpy.core.nan

  def getLargeDependenceEmphasisFeatureValue(self):
    r"""
    Calculate and return the Large Dependence Emphasis (LDE) value.

    :math:`LDE = \frac{\sum^{N_g}_{i=1}\sum^{N_d}_{j=1}{\textbf{P}(i,j)j^2}}{\sum^{N_g}_{i=1}\sum^{N_d}_{j=1}{\textbf{P}(i,j)}}`

    A measure of the distribution of large dependencies, with a greater value indicative
    of larger dependence and more homogeneous textures.
    """
    pd =  self.coefficients['pd']
    jvector = self.coefficients['jvector']
    sumP_gldm = self.coefficients['sumP_gldm']

    try:
      lre = numpy.sum(pd*(jvector**2)) / sumP_gldm
      return lre
    except ZeroDivisionError:
      return numpy.core.nan

  def getGrayLevelNonUniformityFeatureValue(self):
    r"""
    Calculate and return the Gray Level Non-Uniformity (GLN) value.

    :math:`GLN = \frac{\sum^{N_g}_{i=1}\left(\sum^{N_d}_{j=1}{\textbf{P}(i,j)}\right)^2}{\sum^{N_g}_{i=1}\sum^{N_d}_{j=1}{\textbf{P}(i,j)}}`

    Measures the similarity of gray-level intensity values in the image, where a lower GLN value
    correlates with a greater similarity in intensity values.
    """
    pg = self.coefficients['pg']
    sumP_gldm = self.coefficients['sumP_gldm']

    try:
      gln = numpy.sum(pg**2) / sumP_gldm
      return gln
    except ZeroDivisionError:
      return numpy.core.nan

  def getGrayLevelNonUniformityNormalizedFeatureValue(self):
    r"""
    Calculate and return the Gray Level Non-Uniformity Normalized (GLNN) value.

    :math:`GLNN = \frac{\sum^{N_g}_{i=1}\left(\sum^{N_d}_{j=1}{\textbf{P}(i,j)}\right)^2}{\sum^{N_g}_{i=1}\sum^{N_d}_{j=1}{\textbf{P}(i,j)}^2}`

    Measures the similarity of gray-level intensity values in the image, where a lower GLNN value
    correlates with a greater similarity in intensity values. This is the normalized version of the GLN formula.
    """
    pg = self.coefficients['pg']
    sumP_gldm = self.coefficients['sumP_gldm']

    try:
      glnn = numpy.sum(pg**2) / (sumP_gldm**2)
      return glnn
    except ZeroDivisionError:
      return numpy.core.nan

  def getDependenceNonUniformityFeatureValue(self):
    r"""
    Calculate and return the Dependence Non-Uniformity (DN) value.

    :math:`DN = \frac{\sum^{N_d}_{j=1}\left(\sum^{N_g}_{i=1}{\textbf{P}(i,j)}\right)^2}{\sum^{N_g}_{i=1}\sum^{N_d}_{j=1}{\textbf{P}(i,j)}}`

    Measures the similarity of dependence throughout the image, with a lower value indicating
    more homogeneity among dependencies in the image.
    """
    pd = self.coefficients['pd']
    sumP_gldm = self.coefficients['sumP_gldm']

    try:
      dn = numpy.sum(pd**2) / sumP_gldm
      return dn
    except ZeroDivisionError:
      return numpy.core.nan

  def getDependenceNonUniformityNormalizedFeatureValue(self):
    r"""
    Calculate and return the Dependence Non-Uniformity Normalized (DNN) value.

    :math:`DNN = \frac{\sum^{N_d}_{j=1}\left(\sum^{N_g}_{i=1}{\textbf{P}(i,j)}\right)^2}{\sum^{N_g}_{i=1}\sum^{N_d}_{j=1}{\textbf{P}(i,j)}^2}`

    Measures the similarity of dependence throughout the image, with a lower value indicating
    more homogeneity among dependencies in the image. This is the normalized version of the DLN formula.
    """
    pd = self.coefficients['pd']
    sumP_gldm = self.coefficients['sumP_gldm']

    try:
      dnn = numpy.sum(pd**2) / (sumP_gldm**2)
      return dnn
    except ZeroDivisionError:
      return numpy.core.nan

  def getGrayLevelVarianceFeatureValue(self):
    r"""
    Calculate and return the Gray Level Variance (GLV) value.

    :math:`GLV = \displaystyle\sum^{N_g}_{i=1}\displaystyle\sum^{N_d}_{j=1}{p(i,j)(i - \mu)^2}`, where

    :math:`\mu = \displaystyle\sum^{N_g}_{i=1}\displaystyle\sum^{N_d}_{j=1}{ip(i,j)}`

    Measures the variance in grey level in the image.
    """
    ivector = self.coefficients['ivector']
    sumP_gldm = self.coefficients['sumP_gldm']

    u_i = numpy.sum(self.coefficients['pg'] * ivector) / sumP_gldm
    glv = numpy.sum(self.coefficients['pg'] * (ivector - u_i)**2) / sumP_gldm
    return glv

  def getDependenceVarianceFeatureValue(self):
    r"""
    Calculate and return the Dependence Variance (DV) value.

    :math:`DV = \displaystyle\sum^{N_g}_{i=1}\displaystyle\sum^{N_d}_{j=1}{p(i,j)(j - \mu)^2}`, where

    :math:`\mu = \displaystyle\sum^{N_g}_{i=1}\displaystyle\sum^{N_d}_{j=1}{jp(i,j)}`

    Measures the variance in dependence size in the image.
    """
    jvector = self.coefficients['jvector']
    sumP_gldm = self.coefficients['sumP_gldm']
    u_j = numpy.sum(self.coefficients['pd'] * jvector) / sumP_gldm
    dv = numpy.sum(self.coefficients['pd'] * (jvector - u_j)**2) / sumP_gldm
    return dv

  def getDependenceEntropyFeatureValue(self):
    r"""
    Calculate and return the Dependence Entropy (DE) value.

    :math:`Dependence Entropy = -\displaystyle\sum^{N_g}_{i=1}\displaystyle\sum^{N_d}_{j=1}{p(i,j)\log_{2}(p(i,j)+\epsilon)}`
    """
    eps = numpy.spacing(1)
    p_gldm = self.P_gldm / self.coefficients['sumP_gldm']
    return -numpy.sum(p_gldm * numpy.log2(p_gldm + eps))

  def getLowGrayLevelEmphasisFeatureValue(self):
    r"""
    Calculate and return the Low Gray Level Emphasis (LGLE) value.

    :math:`LGLE = \frac{\sum^{N_g}_{i=1}\sum^{N_d}_{j=1}{\frac{\textbf{P}(i,j)}{i^2}}}{\sum^{N_g}_{i=1}\sum^{N_d}_{j=1}{\textbf{P}(i,j)}}`

    Measures the distribution of low gray-level values, with a higher value indicating a greater
    concentration of low gray-level values in the image.
    """
    pg = self.coefficients['pg']
    ivector = self.coefficients['ivector']
    sumP_gldm = self.coefficients['sumP_gldm']

    try:
      lgle = numpy.sum(pg/(ivector**2)) / sumP_gldm
      return lgle
    except ZeroDivisionError:
      return numpy.core.nan

  def getHighGrayLevelEmphasisFeatureValue(self):
    r"""
    Calculate and return the High Gray Level Emphasis (HGLE) value.

    :math:`HGLE = \frac{\sum^{N_g}_{i=1}\sum^{N_d}_{j=1}{\textbf{P}(i,j)i^2}}{\sum^{N_g}_{i=1}\sum^{N_d}_{j=1}{\textbf{P}(i,j)}}`

    Measures the distribution of the higher gray-level values, with a higher value indicating
    a greater concentration of high gray-level values in the image.
    """
    pg = self.coefficients['pg']
    ivector = self.coefficients['ivector']
    sumP_gldm = self.coefficients['sumP_gldm']

    try:
      hgle = numpy.sum(pg*(ivector**2)) / sumP_gldm
      return hgle
    except ZeroDivisionError:
      return numpy.core.nan

  def getSmallDependenceLowGrayLevelEmphasisFeatureValue(self):
    r"""
    Calculate and return the Small Dependence Low Gray Level Emphasis (SDLGLE) value.

    :math:`SDLGLE = \frac{\sum^{N_g}_{i=1}\sum^{N_d}_{j=1}{\frac{\textbf{P}(i,j)}{i^2j^2}}}{\sum^{N_g}_{i=1}\sum^{N_d}_{j=1}{\textbf{P}(i,j)}}`

    Measures the joint distribution of small dependence with lower gray-level values.
    """
    ivector = self.coefficients['ivector']
    jvector = self.coefficients['jvector']
    sumP_gldm = self.coefficients['sumP_gldm']

    try:
      sdlgle = numpy.sum(self.P_gldm/((ivector[:,None]**2)*(jvector[None,:]**2))) / sumP_gldm
      return sdlgle
    except ZeroDivisionError:
      return numpy.core.nan

  def getSmallDependenceHighGrayLevelEmphasisFeatureValue(self):
    r"""
    Calculate and return the Small Dependence High Gray Level Emphasis (SDHGLE) value.

    :math:`SDHGLE = \frac{\sum^{N_g}_{i=1}\sum^{N_d}_{j=1}{\frac{\textbf{P}(i,j)i^2}{j^2}}}{\sum^{N_g}_{i=1}\sum^{N_d}_{j=1}{\textbf{P}(i,j)}}`

    Measures the joint distribution of small dependence with higher gray-level values.
    """
    ivector = self.coefficients['ivector']
    jvector = self.coefficients['jvector']
    sumP_gldm = self.coefficients['sumP_gldm']

    try:
      sdhgle = numpy.sum(self.P_gldm*(ivector[:,None]**2)/(jvector[None,:]**2)) / sumP_gldm
      return sdhgle
    except ZeroDivisionError:
      return numpy.core.nan

  def getLargeDependenceLowGrayLevelEmphasisFeatureValue(self):
    r"""
    Calculate and return the Large Dependence Low Gray Level Emphasis (LDLGLE) value.

    :math:`LDLGLE = \frac{\sum^{N_g}_{i=1}\sum^{N_d}_{j=1}{\frac{\textbf{P}(i,j)j^2}{i^2}}}{\sum^{N_g}_{i=1}\sum^{N_d}_{j=1}{\textbf{P}(i,j)}}`

    Measures the joint distribution of large dependence with lower gray-level values.
    """
    ivector = self.coefficients['ivector']
    jvector = self.coefficients['jvector']
    sumP_gldm = self.coefficients['sumP_gldm']

    try:
      ldlgle = numpy.sum(self.P_gldm*(jvector[None,:]**2)/(ivector[:,None]**2)) / sumP_gldm
      return ldlgle
    except ZeroDivisionError:
      return numpy.core.nan

  def getLargeDependenceHighGrayLevelEmphasisFeatureValue(self):
    r"""
    Calculate and return the Large Dependence High Gray Level Emphasis (LDHGLE) value.

    :math:`LDHGLE = \frac{\sum^{N_g}_{i=1}\sum^{N_d}_{j=1}{\textbf{P}(i,j)i^2j^2}}{\sum^{N_g}_{i=1}\sum^{N_d}_{j=1}{\textbf{P}(i,j)}}`

    Measures the joint distribution of large dependence with higher gray-level values.
    """
    ivector = self.coefficients['ivector']
    jvector = self.coefficients['jvector']
    sumP_gldm = self.coefficients['sumP_gldm']

    try:
      ldhgle = numpy.sum(self.P_gldm*((jvector[None,:]**2)*(ivector[:,None]**2))) / sumP_gldm
      return ldhgle
    except ZeroDivisionError:
      return numpy.core.nan
