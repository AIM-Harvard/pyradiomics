from itertools import chain
import numpy
import SimpleITK as sitk
from radiomics import base, imageoperations


class RadiomicsGLRLM(base.RadiomicsFeaturesBase):
  r"""
  A Gray Level Run Length Matrix (GLRLM) quantifies gray level runs in an image.
  A gray level run is defined as the length in number of pixels,
  of consecutive pixels that have the same gray level value. In a gray level run length matrix
  :math:`\textbf{P}(i,j|\theta)`, the :math:`(i,j)`\ :sup:`th` element describes the number of times
  a gray level :math:`i` appears consecutively :math:`j` times in the direction specified by :math:`\theta`.

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

  :math:`\textbf{P}(i,j|\theta)` be the run length matrix for an arbitrary direction :math:`\theta`

  :math:`p(i,j|\theta)` be the normalized run length matrix, defined as :math:`p(i,j|\theta) =
  \frac{\textbf{P}(i,j|\theta)}{\sum{\textbf{P}(i,j|\theta)}}`

  :math:`N_g` be the number of discreet intensity values in the image

  :math:`N_r` be the number of discreet run lengths in the image

  :math:`N_p` be the number of voxels in the image

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

    In case of other values, an warning is logged and GLCMs are all weighted by factor 1 and summed.

  References

  - Galloway MM. 1975. Texture analysis using gray level run lengths. Computer Graphics and Image Processing 4:172-179.

  - Tang X. 1998. Texture information in run-length matrices. IEEE Transactions on Image Processing 7(11):1602-1609.
  """

  def __init__(self, inputImage, inputMask, **kwargs):
    super(RadiomicsGLRLM, self).__init__(inputImage, inputMask, **kwargs)

    self.weightingNorm = kwargs.get('weightingNorm', None)  # manhattan, euclidean, infinity

    self.coefficients = {}
    self.P_glrlm = {}

    # binning
    self.matrix, self.histogram = imageoperations.binImage(self.binWidth, self.targetVoxelArray, self.matrix,
                                                           self.matrixCoordinates)
    self.coefficients['Ng'] = self.histogram[1].shape[0] - 1
    self.coefficients['Nr'] = numpy.max(self.matrix.shape)
    self.coefficients['Np'] = self.targetVoxelArray.size

    self._calculateGLRLM()
    self._calculateCoefficients()

  def _calculateGLRLM(self):
    Ng = self.coefficients['Ng']
    Nr = self.coefficients['Nr']

    padVal = -2000  # use eps or NaN to pad matrix
    self.matrix[(self.maskArray == 0)] = padVal

    matrixDiagonals = []

    size = numpy.max(self.matrixCoordinates, 1) - numpy.min(self.matrixCoordinates, 1) + 1
    angles = imageoperations.generateAngles(size)

    for angle in angles:
      staticDims, = numpy.where(angle == 0)  # indices for static dimensions for current angle (z, y, x)
      movingDims, = numpy.where(angle != 0)  # indices for moving dimensions for current angle (z, y, x)

      if len(movingDims) == 1:  # movement in one dimension, e.g. angle (0, 0, 1)
        T = tuple(numpy.append(staticDims, movingDims))
        diags = chain.from_iterable(numpy.transpose(self.matrix, T))

      elif len(movingDims) == 2:  # movement in two dimension, e.g. angle (0, 1, 1)
        d1 = movingDims[0]
        d2 = movingDims[1]
        direction = numpy.where(angle < 0, -1, 1)
        diags = chain.from_iterable([self.matrix[::direction[0], ::direction[1], ::direction[2]].diagonal(a, d1, d2)
                                     for a in xrange(-self.matrix.shape[d1] + 1, self.matrix.shape[d2])])

      else:  # movement in 3 dimensions, e.g. angle (1, 1, 1)
        diags = []
        direction = numpy.where(angle < 0, -1, 1)
        for h in [self.matrix[::direction[0], ::direction[1], ::direction[2]].diagonal(a, 0, 1)
                  for a in xrange(-self.matrix.shape[0] + 1, self.matrix.shape[1])]:
          diags.extend([h.diagonal(b, 0, 1) for b in xrange(-h.shape[0] + 1, h.shape[1])])

      matrixDiagonals.append(filter(lambda diag: numpy.nonzero(diag != padVal)[0].size > 0, diags))

    P_glrlm = numpy.zeros((Ng, Nr, int(len(matrixDiagonals))))

    # Run-Length Encoding (rle) for the list of diagonals
    # (1 list per direction/angle)
    for angle_idx, angle in enumerate(matrixDiagonals):
      P = P_glrlm[:, :, angle_idx]
      # Check whether delineation is 2D for current angle (all diagonals contain 0 or 1 non-pad value)
      isMultiElement = False
      for d in angle:
        if numpy.where(d != padVal)[0].shape[0] > 1:
          isMultiElement = True
          break
      if isMultiElement:
        for diagonal in angle:
          pos, = numpy.where(numpy.diff(diagonal) != 0)
          pos = numpy.concatenate(([0], pos + 1, [len(diagonal)]))
          rle = zip([int(n) for n in diagonal[pos[:-1]]], pos[1:] - pos[:-1])
          for level, run_length in rle:
            if level != padVal:
              P[level - 1, run_length - 1] += 1

    # Crop gray-level axis of GLRLMs to between minimum and maximum observed gray-levels
    # Crop run-length axis of GLRLMs up to maximum observed run-length
    P_glrlm_bounds = numpy.argwhere(P_glrlm)
    (xstart, ystart, zstart), (xstop, ystop, zstop) = P_glrlm_bounds.min(0), P_glrlm_bounds.max(0) + 1
    self.P_glrlm = P_glrlm[xstart:xstop, :ystop, :]

    # Optionally apply a weighting factor
    if not self.weightingNorm is None:
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

      self.P_glrlm = numpy.sum(self.P_glrlm * weights[None, None, :], 2, keepdims=True)

    sumP_glrlm = numpy.sum(self.P_glrlm, (0, 1))

    # Delete empty angles if no weighting is applied
    if self.P_glrlm.shape[2] > 1:
      self.P_glrlm = numpy.delete(self.P_glrlm, numpy.where(sumP_glrlm == 0), 2)
      sumP_glrlm = numpy.delete(sumP_glrlm, numpy.where(sumP_glrlm == 0), 0)

    self.coefficients['sumP_glrlm'] = sumP_glrlm

  def _calculateCoefficients(self):

    pr = numpy.sum(self.P_glrlm, 0)
    pg = numpy.sum(self.P_glrlm, 1)

    ivector = numpy.arange(1, self.P_glrlm.shape[0] + 1, dtype=numpy.float64)
    jvector = numpy.arange(1, self.P_glrlm.shape[1] + 1, dtype=numpy.float64)

    self.coefficients['pr'] = pr
    self.coefficients['pg'] = pg
    self.coefficients['ivector'] = ivector
    self.coefficients['jvector'] = jvector

  def getShortRunEmphasisFeatureValue(self):
    r"""
    Calculate and return the mean Short Run Emphasis (SRE) value for all GLRLMs.

    :math:`SRE = \frac{\sum^{N_g}_{i=1}\sum^{N_r}_{j=1}{\frac{\textbf{P}(i,j|\theta)}{i^2}}}{\sum^{N_g}_{i=1}\sum^{N_r}_{j=1}{\textbf{P}(i,j|\theta)}}`

    A measure of the distribution of short run lengths, with a greater value indicative
    of shorter run lengths and more fine textural textures.
    """
    pr = self.coefficients['pr']
    jvector = self.coefficients['jvector']
    sumP_glrlm = self.coefficients['sumP_glrlm']

    try:
      sre = numpy.sum((pr / (jvector[:, None] ** 2)), 0) / sumP_glrlm
      return (sre.mean())
    except ZeroDivisionError:
      return numpy.core.nan

  def getLongRunEmphasisFeatureValue(self):
    r"""
    Calculate and return the mean Long Run Emphasis (LRE) value for all GLRLMs.

    :math:`LRE = \frac{\sum^{N_g}_{i=1}\sum^{N_r}_{j=1}{\textbf{P}(i,j|\theta)j^2}}{\sum^{N_g}_{i=1}\sum^{N_r}_{j=1}{\textbf{P}(i,j|\theta)}}`

    A measure of the distribution of long run lengths, with a greater value indicative
    of longer run lengths and more coarse structural textures.
    """
    pr = self.coefficients['pr']
    jvector = self.coefficients['jvector']
    sumP_glrlm = self.coefficients['sumP_glrlm']

    try:
      lre = numpy.sum((pr * (jvector[:, None] ** 2)), 0) / sumP_glrlm
      return (lre.mean())
    except ZeroDivisionError:
      return numpy.core.nan

  def getGrayLevelNonUniformityFeatureValue(self):
    r"""
    Calculate and return the mean Gray Level Non-Uniformity (GLN) value for all GLRLMs.

    :math:`GLN = \frac{\sum^{N_g}_{i=1}\left(\sum^{N_r}_{j=1}{\textbf{P}(i,j|\theta)}\right)^2}{\sum^{N_g}_{i=1}\sum^{N_r}_{j=1}{\textbf{P}(i,j|\theta)}}`

    Measures the similarity of gray-level intensity values in the image, where a lower GLN value
    correlates with a greater similarity in intensity values.
    """
    pg = self.coefficients['pg']
    sumP_glrlm = self.coefficients['sumP_glrlm']

    try:
      gln = numpy.sum((pg ** 2), 0) / sumP_glrlm
      return (gln.mean())
    except ZeroDivisionError:
      return numpy.core.nan

  def getGrayLevelNonUniformityNormalizedFeatureValue(self):
    r"""
    Calculate and return the Gray Level Non-Uniformity Normalized (GLNN) value.

    :math:`GLNN = \frac{\sum^{N_g}_{i=1}\left(\sum^{N_r}_{j=1}{\textbf{P}(i,j|\theta)}\right)^2}{\sum^{N_g}_{i=1}\sum^{N_r}_{j=1}{\textbf{P}(i,j|\theta)}^2}`

    Measures the similarity of gray-level intensity values in the image, where a lower GLNN value
    correlates with a greater similarity in intensity values. This is the normalized version of the GLN formula.
    """
    pg = self.coefficients['pg']
    sumP_gldm = self.coefficients['sumP_glrlm']

    try:
      glnn = numpy.sum(pg ** 2, 0) / (sumP_gldm ** 2)
      return glnn.mean()
    except ZeroDivisionError:
      return numpy.core.nan

  def getRunLengthNonUniformityFeatureValue(self):
    r"""
    Calculate and return the mean Run Length Non-Uniformity (RLN) value for all GLRLMs.

    :math:`RLN = \frac{\sum^{N_r}_{j=1}\left(\sum^{N_g}_{i=1}{\textbf{P}(i,j|\theta)}\right)^2}{\sum^{N_g}_{i=1}\sum^{N_r}_{j=1}{\textbf{P}(i,j|\theta)}}`

    Measures the similarity of run lengths throughout the image, with a lower value indicating
    more homogeneity among run lengths in the image.
    """
    pr = self.coefficients['pr']
    sumP_glrlm = self.coefficients['sumP_glrlm']

    try:
      rln = numpy.sum((pr ** 2), 0) / sumP_glrlm
      return (rln.mean())
    except ZeroDivisionError:
      return numpy.core.nan

  def getRunLengthNonUniformityNormalizedFeatureValue(self):
    r"""
    Calculate and return the mean Run Length Non-Uniformity Normalized (RLNN) value for all GLRLMs.

    :math:`RLNN = \frac{\sum^{N_r}_{j=1}\left(\sum^{N_g}_{i=1}{\textbf{P}(i,j|\theta)}\right)^2}{\sum^{N_g}_{i=1}\sum^{N_r}_{j=1}{\textbf{P}(i,j|\theta)}}`

    Measures the similarity of run lengths throughout the image, with a lower value indicating
    more homogeneity among run lengths in the image. This is the normalized version of the RLN formula.
    """
    pr = self.coefficients['pr']
    sumP_glrlm = self.coefficients['sumP_glrlm']

    try:
      rlnn = numpy.sum((pr ** 2), 0) / sumP_glrlm ** 2
      return (rlnn.mean())
    except ZeroDivisionError:
      return numpy.core.nan

  def getRunPercentageFeatureValue(self):
    r"""
    Calculate and return the mean Run Percentage (RP) value for all GLRLMs.

    :math:`RP = \displaystyle\sum^{N_g}_{i=1}\displaystyle\sum^{N_r}_{j=1}{\frac{\textbf{P}(i,j|\theta)}{N_p}}`

    Measures the homogeneity and distribution of runs of an image.
    """
    Np = self.coefficients['Np']

    try:
      rp = numpy.sum((self.P_glrlm / (Np)), (0, 1))
      return (rp.mean())
    except ZeroDivisionError:
      return numpy.core.nan

  def getGrayLevelVarianceFeatureValue(self):
    r"""
    Calculate and return the Gray Level Variance (GLV) value.

    :math:`GLV = \displaystyle\sum^{N_g}_{i=1}\displaystyle\sum^{N_r}_{j=1}{p(i,j|\theta)(i - \mu)^2}`, where

    :math:`\mu = \displaystyle\sum^{N_g}_{i=1}\displaystyle\sum^{N_r}_{j=1}{p(i,j|\theta)i}`

    Measures the variance in gray level intensity for the runs.
    """
    ivector = self.coefficients['ivector']
    sumP_glrlm = self.coefficients['sumP_glrlm']
    u_i = numpy.sum(self.coefficients['pg'] * ivector[:, None], 0) / sumP_glrlm
    glv = numpy.sum(self.coefficients['pg'] * (ivector[:, None] - u_i[None, :]) ** 2, 0) / sumP_glrlm
    return glv.mean()

  def getRunVarianceFeatureValue(self):
    r"""
    Calculate and return the Run Variance (RV) value.

    :math:`RV = \displaystyle\sum^{N_g}_{i=1}\displaystyle\sum^{N_r}_{j=1}{p(i,j|\theta)(j - \mu)^2}`, where

    :math:`\mu = \displaystyle\sum^{N_g}_{i=1}\displaystyle\sum^{N_r}_{j=1}{p(i,j|\theta)j}`

    Measures the variance in runs for the run lengths.
    """
    jvector = self.coefficients['jvector']
    sumP_glrlm = self.coefficients['sumP_glrlm']
    u_j = numpy.sum(self.coefficients['pr'] * jvector[:, None], 0) / sumP_glrlm
    rv = numpy.sum(self.coefficients['pr'] * (jvector[:, None] - u_j[None, :]) ** 2, 0) / sumP_glrlm
    return rv.mean()

  def getRunEntropyFeatureValue(self):
    r"""1
    Calculate and return the Run Entropy (RE) value.

    :math:`RE = -\displaystyle\sum^{N_g}_{i=1}\displaystyle\sum^{N_r}_{j=1}{p(i,j|\theta)\log_{2}(p(i,j|\theta)+\epsilon)}`
    """
    eps = numpy.spacing(1)
    p_glrlm = self.P_glrlm / self.coefficients['sumP_glrlm']
    re = -numpy.sum(p_glrlm * numpy.log2(p_glrlm + eps), (0, 1))
    return re.mean()

  def getLowGrayLevelRunEmphasisFeatureValue(self):
    r"""
    Calculate and return the mean Low Gray Level Run Emphasis (LGLRE) value for all GLRLMs.

    :math:`LGLRE = \frac{\sum^{N_g}_{i=1}\sum^{N_r}_{j=1}{\frac{\textbf{P}(i,j|\theta)}{i^2}}}{\sum^{N_g}_{i=1}\sum^{N_r}_{j=1}{\textbf{P}(i,j|\theta)}}`

    Measures the distribution of low gray-level values, with a higher value indicating a greater
    concentration of low gray-level values in the image.
    """
    pg = self.coefficients['pg']
    ivector = self.coefficients['ivector']
    sumP_glrlm = self.coefficients['sumP_glrlm']

    try:
      lglre = numpy.sum((pg / (ivector[:, None] ** 2)), 0) / sumP_glrlm
      return (lglre.mean())
    except ZeroDivisionError:
      return numpy.core.nan

  def getHighGrayLevelRunEmphasisFeatureValue(self):
    r"""
    Calculate and return the mean High Gray Level Run Emphasis (HGLRE) value for all GLRLMs.

    :math:`HGLRE = \frac{\sum^{N_g}_{i=1}\sum^{N_r}_{j=1}{\textbf{P}(i,j|\theta)i^2}}{\sum^{N_g}_{i=1}\sum^{N_r}_{j=1}{\textbf{P}(i,j|\theta)}}`

    Measures the distribution of the higher gray-level values, with a higher value indicating
    a greater concentration of high gray-level values in the image.
    """
    pg = self.coefficients['pg']
    ivector = self.coefficients['ivector']
    sumP_glrlm = self.coefficients['sumP_glrlm']

    try:
      hglre = numpy.sum((pg * (ivector[:, None] ** 2)), 0) / sumP_glrlm
      return (hglre.mean())
    except ZeroDivisionError:
      return numpy.core.nan

  def getShortRunLowGrayLevelEmphasisFeatureValue(self):
    r"""
    Calculate and return the mean Short Run Low Gray Level Emphasis (SRLGLE) value for all GLRLMs.

    :math:`SRLGLE = \frac{\sum^{N_g}_{i=1}\sum^{N_r}_{j=1}{\frac{\textbf{P}(i,j|\theta)}{i^2j^2}}}{\sum^{N_g}_{i=1}\sum^{N_r}_{j=1}{\textbf{P}(i,j|\theta)}}`

    Measures the joint distribution of shorter run lengths with lower gray-level values.
    """
    ivector = self.coefficients['ivector']
    jvector = self.coefficients['jvector']
    sumP_glrlm = self.coefficients['sumP_glrlm']

    try:
      srlgle = numpy.sum((self.P_glrlm / ((ivector[:, None, None] ** 2) * (jvector[None, :, None] ** 2))),
                         (0, 1)) / sumP_glrlm
      return (srlgle.mean())
    except ZeroDivisionError:
      return numpy.core.nan

  def getShortRunHighGrayLevelEmphasisFeatureValue(self):
    r"""
    Calculate and return the mean Short Run High Gray Level Emphasis (SRHGLE) value for all GLRLMs.

    :math:`SRHGLE = \frac{\sum^{N_g}_{i=1}\sum^{N_r}_{j=1}{\frac{\textbf{P}(i,j|\theta)i^2}{j^2}}}{\sum^{N_g}_{i=1}\sum^{N_r}_{j=1}{\textbf{P}(i,j|\theta)}}`

    Measures the joint distribution of shorter run lengths with higher gray-level values.
    """
    ivector = self.coefficients['ivector']
    jvector = self.coefficients['jvector']
    sumP_glrlm = self.coefficients['sumP_glrlm']

    try:
      srhgle = numpy.sum((self.P_glrlm * (ivector[:, None, None] ** 2) / (jvector[None, :, None] ** 2)),
                         (0, 1)) / sumP_glrlm
      return (srhgle.mean())
    except ZeroDivisionError:
      return numpy.core.nan

  def getLongRunLowGrayLevelEmphasisFeatureValue(self):
    r"""
    Calculate and return the mean Long Run Low Gray Level Emphasis (LRLGLE) value for all GLRLMs.

    :math:`LRLGLRE = \frac{\sum^{N_g}_{i=1}\sum^{N_r}_{j=1}{\frac{\textbf{P}(i,j|\theta)j^2}{i^2}}}{\sum^{N_g}_{i=1}\sum^{N_r}_{j=1}{\textbf{P}(i,j|\theta)}}`

    Measures the joint distribution of long run lengths with lower gray-level values.
    """
    ivector = self.coefficients['ivector']
    jvector = self.coefficients['jvector']
    sumP_glrlm = self.coefficients['sumP_glrlm']

    try:
      lrlgle = numpy.sum((self.P_glrlm * (jvector[None, :, None] ** 2) / (ivector[:, None, None] ** 2)),
                         (0, 1)) / sumP_glrlm
      return (lrlgle.mean())
    except ZeroDivisionError:
      return numpy.core.nan

  def getLongRunHighGrayLevelEmphasisFeatureValue(self):
    r"""
    Calculate and return the mean Long Run High Gray Level Emphasis (LRHGLE) value for all GLRLMs.

    :math:`LRHGLRE = \frac{\sum^{N_g}_{i=1}\sum^{N_r}_{j=1}{\textbf{P}(i,j|\theta)i^2j^2}}{\sum^{N_g}_{i=1}\sum^{N_r}_{j=1}{\textbf{P}(i,j|\theta)}}`

    Measures the joint distribution of long run lengths with higher gray-level values.
    """
    ivector = self.coefficients['ivector']
    jvector = self.coefficients['jvector']
    sumP_glrlm = self.coefficients['sumP_glrlm']

    try:
      lrhgle = numpy.sum((self.P_glrlm * ((jvector[None, :, None] ** 2) * (ivector[:, None, None] ** 2))),
                         (0, 1)) / sumP_glrlm
      return (lrhgle.mean())
    except ZeroDivisionError:
      return numpy.core.nan
