import numpy
from tqdm import trange
import SimpleITK as sitk
from radiomics import base, imageoperations


class RadiomicsGLDZM(base.RadiomicsFeaturesBase):
  r"""
  The Gray Level Distance Zone Matrix (GLDZM) quantifies the relationship between the location and gray level of the
  connected zones in an image. A zone consists of connected neighbours which have the same gray level intensity.
  Furthermore, a connected neighbour is classified using 26-connectedness for 3D and 8-connectedness for 2D images
  (fully-connected).

  Using SimpleITK functionality, a signed Maurer distance map is generated, which yields the euclidean distance of a
  voxel to the edge of the Region of Interest (ROI). This distance is then rounded to an integer, and 1 is added, so
  that voxels on the edge of the ROI are given a distance of 1. This ensures a discreet number of distances in the
  matrix, while zones on the edge are not ignored (where a multiplication by the distance is needed).

  The location of a zone is the defined minimum distance value of the voxels contained within the zone.

  In a GLDZM :math:`\textbf{P}(i,j)`, the :math:`(i,j)`\ :sup:`th` element describes the number of zones in an image
  with gray level :math:`i` and located at distance :math:`j` from the edge of the ROI.

  As a 2 dimensional example, consider the following 7x7 image, with 5 discreet gray values and where a value of
  :math:`N` indicates a pixel not belonging to the ROI:

  .. math::
    \textbf{I} = \begin{bmatrix}
    N & N & N & 4 & 4 & 4 & N\\
    N & N & 3 & 1 & 3 & 4 & N\\
    2 & 1 & 1 & 1 & 3 & 2 & N\\
    4 & 4 & 2 & 2 & 3 & 2 & 1\\
    3 & 5 & 3 & 3 & 2 & 1 & 1\\
    3 & 5 & 3 & 3 & 2 & 4 & N\\
    3 & 1 & N & N & N & 4 & N\end{bmatrix}

  Then the corresponding distance map is:

  .. math::
    \textbf{D} = \begin{bmatrix}
    N & N & N & 1 & 1 & 1 & N\\
    N & N & 1 & 1 & 2 & 1 & N\\
    1 & 1 & 2 & 2 & 2 & 1 & N\\
    1 & 2 & 2 & 2 & 2 & 1 & 1\\
    1 & 2 & 2 & 2 & 2 & 1 & 1\\
    1 & 1 & 1 & 1 & 1 & 1 & N\\
    1 & 1 & N & N & N & 1 & N\end{bmatrix}

  And the GLDZM then becomes:

  .. math::
    \textbf{P}=\begin{bmatrix}
    3 & 0\\
    3 & 1\\
    3 & 1\\
    2 & 0\\
    1 & 1\end{bmatrix}

  Let

  :math:`\textbf{P}(i,j)` be the distance zone matrix

  :math:`p(i,j)` be the normalized distance zone matrix, defined as :math:`p(i,j) =
  \frac{\textbf{P}(i,j)}{\sum{\textbf{P}(i,j)}}`

  :math:`N_g` be the number of discreet intensity values in the image

  :math:`N_d` be the number of discreet distances in the image

  :math:`N_p` be the number of voxels in the image

  N.B. For the distance map to calculate correctly in 2D images, image and mask must be cropped on the bounding box of
  the image and the mask! This can be done by a call to :py:func:`~imageoperations.cropToTumorMask`.

  References

    - Thibault, G., Angulo, J., and Meyer, F. (2014); Advanced statistical matrices for texture characterization:
      application to cell classification; IEEE transactions on bio-medical engineering, 61(3):630-7.
  """

  def __init__(self, inputImage, inputMask, **kwargs):
    super(RadiomicsGLDZM, self).__init__(inputImage, inputMask, **kwargs)

    self.coefficients = {}
    self.P_gldzm = {}

    # Pad inputMask to prevent errors in the distancemap (voxels on the border are not considered to be on the edge of
    # the ROI. Do not pad in 2D directions, as this would result in a distance map where all distances are 0 (i.e.
    # bordering on the next slice.
    cpif = sitk.ConstantPadImageFilter()

    # Only apply padding where the size is more than 1. Flip the direction of the matrix, as padding will be
    # applied on the image (x, y, z) and the matrix is (z, y, x)
    padding = numpy.where(numpy.array(self.matrix.shape)[::-1] == 1, 0, 1)

    cpif.SetPadLowerBound(padding)
    cpif.SetPadUpperBound(padding)

    self.inputImage = cpif.Execute(self.inputImage)
    self.inputMask = cpif.Execute(self.inputMask)

    # Reassign matrix
    self.matrix = sitk.GetArrayFromImage(self.inputImage).astype('float')
    # Reassign self.maskArray using the now-padded self.inputMask and make it binary
    self.maskArray = (sitk.GetArrayFromImage(self.inputMask) == self.label).astype('int')
    self.matrixCoordinates = numpy.where(self.maskArray != 0)

    # binning
    self.matrix, self.histogram = imageoperations.binImage(self.binWidth, self.matrix, self.matrixCoordinates)
    self.coefficients['Ng'] = self.histogram[1].shape[0] - 1
    self.coefficients['Np'] = self.targetVoxelArray.size

    self._calculateGLDZM()
    self._calculateCoefficients()

  def _calculateGLDZM(self):
    """
    Number of times a region with a gray level :math:`i` and occurs with a minimum distance :math:`j` in an image.
    P_gldzm[level, distance] = # occurrences.
    """
    size = numpy.max(self.matrixCoordinates, 1) - numpy.min(self.matrixCoordinates, 1) + 1
    angles = imageoperations.generateAngles(size)

    smdmif = sitk.SignedMaurerDistanceMapImageFilter()
    smdmif.SquaredDistanceOff()
    smdmif.InsideIsPositiveOn()
    distImage = smdmif.Execute(self.inputMask)
    distMap = numpy.round(sitk.GetArrayFromImage(distImage), 0)  # Round distances to make them usable as indices

    # Empty GLDZ matrix
    P_gldzm = numpy.zeros((self.coefficients['Ng'], int(numpy.max(distMap) + 1)))

    # Iterate over all gray levels in the image
    grayLevels = numpy.unique(self.matrix[self.matrixCoordinates])

    if self.verbose: bar = trange(len(grayLevels), desc='calculate GLDZM')

    for i in grayLevels:
      # give some progress
      if self.verbose: bar.update()

      ind = zip(*numpy.where(self.matrix == i))
      ind = list(set(ind).intersection(set(zip(*self.matrixCoordinates))))

      while ind:  # check if ind is not empty: unprocessed regions for current gray level
        # Pop first coordinate of an unprocessed zone, start new stack
        ind_region = [ind.pop()]

        # Define minDistance
        minDistance = -1

        # Grow zone for item popped from stack of region indices, loop until stack of region indices is exhausted
        # store minimum distance
        while ind_region:
          # Use pop to remove next node for set of unprocessed region indices
          ind_node = ind_region.pop()

          # get all coordinates in the 6-connected region, 2 voxels per angle
          region_full = [tuple(sum(a) for a in zip(ind_node, angle_i)) for angle_i in angles]
          region_full += [tuple(sum(a) for a in zip(ind_node, angle_i)) for angle_i in angles * -1]

          # get all unprocessed coordinates in the 26-connected region with same gray level
          region_level = list(set(ind).intersection(set(region_full)))

          # Remove already processed indices to prevent reprocessing
          ind = list(set(ind) - set(region_level))

          # Add all found neighbours to the total stack of unprocessed neighbours
          ind_region.extend(region_level)

          if minDistance < 0 or distMap[ind_node] < minDistance:
            # minDistance is not set or new minimum is found
            minDistance = distMap[ind_node]

        # Update the gray level distance zone matrix, minDistance starts at 0 (voxels on the edge of the ROI.
        P_gldzm[int(i - 1), int(minDistance - 1)] += 1

    if self.verbose: bar.close()

    # Crop gray-level axis of GLDZM matrix to between minimum and maximum observed gray-levels
    # Crop distance-zone axis of GLDZM matrix up to maximum observed distance
    P_gldzm_bounds = numpy.argwhere(P_gldzm)
    (xstart, ystart), (xstop, ystop) = P_gldzm_bounds.min(0), P_gldzm_bounds.max(0) + 1
    self.P_gldzm = P_gldzm[xstart:xstop, :ystop]

  def _calculateCoefficients(self):
    sumP_gldzm = numpy.sum(self.P_gldzm, (0, 1))

    # set sum to numpy.spacing(1) if sum is 0?
    if sumP_gldzm == 0:
      sumP_gldzm = 1

    pd = numpy.sum(self.P_gldzm, 0)
    pg = numpy.sum(self.P_gldzm, 1)

    ivector = numpy.arange(1, self.P_gldzm.shape[0] + 1, dtype=numpy.float64)
    jvector = numpy.arange(1, self.P_gldzm.shape[1] + 1, dtype=numpy.float64)  # ensure distances start at 1

    self.coefficients['sumP_gldzm'] = sumP_gldzm
    self.coefficients['pd'] = pd
    self.coefficients['pg'] = pg
    self.coefficients['ivector'] = ivector
    self.coefficients['jvector'] = jvector

  def getSmallDistanceEmphasisFeatureValue(self):
    r"""
    Calculate and return the Small Distance Emphasis (SDE) value.

    :math:`SDE = \frac{\sum^{N_g}_{i=1}\sum^{N_d}_{j=1}{\frac{\textbf{P}(i,j)}{j^2}}}{\sum^{N_g}_{i=1}\sum^{N_d}_{j=1}{\textbf{P}(i,j)}}`

    A measure of the distribution of small distance zones, with a greater value indicative
    of more smaller distances.
    """
    try:
      sde = numpy.sum(self.coefficients['pd'] / (self.coefficients['jvector'] ** 2)) / self.coefficients['sumP_gldzm']
    except ZeroDivisionError:
      sde = numpy.core.numeric.NaN
    return (sde)

  def getLargeDistanceEmphasisFeatureValue(self):
    r"""
    Calculate and return the Large Distance Emphasis (LDE) value.

    :math:`LDE = \frac{\sum^{N_g}_{i=1}\sum^{N_d}_{j=1}{\textbf{P}(i,j)j^2}}{\sum^{N_g}_{i=1}\sum^{N_d}_{j=1}{\textbf{P}(i,j)}}`

    A measure of the distribution of large distance zones, with a greater value indicative
    of larger distances.
    """
    try:
      lde = numpy.sum(self.coefficients['pd'] * (self.coefficients['jvector'] ** 2)) / self.coefficients['sumP_gldzm']
    except ZeroDivisionError:
      lde = numpy.core.numeric.NaN
    return (lde)

  def getIntensityVariabilityFeatureValue(self):
    r"""
    Calculate and return the Intensity Variability (IV) value.

    :math:`IV = \frac{\sum^{N_g}_{i=1}\left(\sum^{N_d}_{j=1}{\textbf{P}(i,j)}\right)^2}{\sum^{N_g}_{i=1}\sum^{N_d}_{j=1}{\textbf{P}(i,j)}}`

    Measures the variability of gray-level intensity values in the image, where a lower IV value
    correlates with more homogeneity in intensity values.
    """
    try:
      iv = numpy.sum(self.coefficients['pg'] ** 2) / self.coefficients['sumP_gldzm']
    except ZeroDivisionError:
      iv = numpy.core.numeric.NaN
    return (iv)

  def getIntensityVariabilityNormalizedFeatureValue(self):
    r"""
    Calculate and return the Intensity Variability Normalized (IVN) value.

    :math:`IVN = \frac{\sum^{N_g}_{i=1}\left(\sum^{N_d}_{j=1}{\textbf{P}(i,j)}\right)^2}{\sum^{N_g}_{i=1}\sum^{N_d}_{j=1}{\textbf{P}(i,j)}^2}`

    Measures the variability of gray-level intensity values in the image, where a lower IVN value
    correlates with a greater similarity in intensity values.
    """
    try:
      iv = numpy.sum(self.coefficients['pg'] ** 2) / self.coefficients['sumP_gldzm'] ** 2
    except ZeroDivisionError:
      iv = numpy.core.numeric.NaN
    return (iv)

  def getDistanceZoneVariabilityFeatureValue(self):
    r"""
    Calculate and return the Distance-Zone Variability (DZV) value.

    :math:`DZV = \frac{\sum^{N_d}_{j=1}\left(\sum^{N_g}_{i=1}{\textbf{P}(i,j)}\right)^2}{\sum^{N_g}_{i=1}\sum^{N_d}_{j=1}{\textbf{P}(i,j)}}`

    Measures the variability of distances in the image.
    """
    try:
      szv = numpy.sum(self.coefficients['pd'] ** 2) / self.coefficients['sumP_gldzm']
    except ZeroDivisionError:
      szv = numpy.core.numeric.NaN
    return (szv)

  def getDistanceZoneVariabilityNormalizedFeatureValue(self):
    r"""
    Calculate and return the Distance-Zone Variability Normalized (DZVN) value.

    :math:`DZVN = \frac{\sum^{N_d}_{j=1}\left(\sum^{N_g}_{i=1}{\textbf{P}(i,j)}\right)^2}{\sum^{N_g}_{i=1}\sum^{N_d}_{j=1}{\textbf{P}(i,j)}^2}`

    Measures the variability of distances throughout the image, with a lower value indicating
    more homogeneity among distances in the image. This is the normalized version of the DZVN formula.
    """
    try:
      szv = numpy.sum(self.coefficients['pd'] ** 2) / self.coefficients['sumP_gldzm'] ** 2
    except ZeroDivisionError:
      szv = numpy.core.numeric.NaN
    return (szv)

  def getZonePercentageFeatureValue(self):
    r"""
    Calculate and return the Zone Percentage (ZP) value.

    :math:`ZP = \sum^{N_g}_{i=1}\sum^{N_d}_{j=1}{\frac{\textbf{P}(i,j)}{N_p}}`

    Measures the homogeneity of the distribution of distances in an image among the observed gray-levels.
    """
    try:
      zp = self.coefficients['sumP_gldzm'] / self.coefficients['Np']
    except ZeroDivisionError:
      zp = numpy.core.numeric.NaN
    return (zp)

  def getGrayLevelVariance(self):
    r"""
    Calculate and return the Gray Level Variance (GLV) value.

    :math:`GLV = \displaystyle\sum^{N_g}_{i=1}\displaystyle\sum^{N_d}_{j=1}{p(i,j)(i - \mu)^2}`, where

    :math:`\mu = \displaystyle\sum^{N_g}_{i=1}\displaystyle\sum^{N_d}_{j=1}{ip(i,j)}`

    Measures the variance in distance counts for the grey levels.
    """
    ivector = self.coefficients['ivector']
    u_i = numpy.sum(self.coefficients['pg'] * ivector[:]) / self.coefficients['sumP_gldzm']
    glv = numpy.sum(self.coefficients['pg'] * (ivector - u_i) ** 2) / self.coefficients['sumP_gldzm']
    return glv

  def getDistanceVariance(self):
    r"""
    Calculate and return the Distance Variance (DV) value.

    :math:`DV = \displaystyle\sum^{N_g}_{i=1}\displaystyle\sum^{N_d}_{j=1}{p(i,j)(j - \mu)^2}`, where

    :math:`\mu = \displaystyle\sum^{N_g}_{i=1}\displaystyle\sum^{N_d}_{j=1}{jp(i,j)}`

    Measures the variance in distance counts for the dependence sizes.
    """
    jvector = self.coefficients['jvector']
    u_j = numpy.sum(self.coefficients['pd'] * jvector) / self.coefficients['sumP_gldzm']
    dv = numpy.sum(self.coefficients['pd'] * (jvector - u_j) ** 2) / self.coefficients['sumP_gldzm']
    return dv

  def getDistanceEntropy(self):
    r"""
    Calculate and return the Distance Entropy (DE) value.

    :math:`DE = -\displaystyle\sum^{N_g}_{i=1}\displaystyle\sum^{N_d}_{j=1}{p(i,j)\log_{2}(p(i,j)+\epsilon)}`
    """
    eps = numpy.spacing(1)
    p_gldzm = self.P_gldzm / self.coefficients['sumP_gldzm']
    return -numpy.sum(p_gldzm * numpy.log2(p_gldzm + eps))

  def getLowIntensityEmphasisFeatureValue(self):
    r"""
    Calculate and return the Low Intensity Emphasis (LIE) value.

    :math:`LIE = \frac{\sum^{N_g}_{i=1}\sum^{N_d}_{j=1}{\frac{\textbf{P}(i,j)}{i^2}}}{\sum^{N_g}_{i=1}\sum^{N_d}_{j=1}{\textbf{P}(i,j)}}`

    Measures the distribution of lower gray-level distance zones, with a higher value indicating a greater
    proportion of lower gray-level values and distance zones in the image.
    """
    try:
      lie = numpy.sum((self.coefficients['pg'] / (self.coefficients['ivector'] ** 2))) / self.coefficients['sumP_gldzm']
    except ZeroDivisionError:
      lie = numpy.core.numeric.NaN
    return (lie)

  def getHighIntensityEmphasisFeatureValue(self):
    r"""
    Calculate and return the High Intensity Emphasis (HIE) value.

    :math:`HIE = \frac{\sum^{N_g}_{i=1}\sum^{N_d}_{j=1}{\textbf{P}(i,j)i^2}}{\sum^{N_g}_{i=1}\sum^{N_d}_{j=1}{\textbf{P}(i,j)}}`

    Measures the distribution of the higher gray-level values, with a higher value indicating
    a greater proportion of higher gray-level values and distance zones in the image.
    """
    try:
      hie = numpy.sum((self.coefficients['pg'] * (self.coefficients['ivector'] ** 2))) / self.coefficients['sumP_gldzm']
    except ZeroDivisionError:
      hie = numpy.core.numeric.NaN
    return (hie)

  def getLowIntensitySmallDistanceEmphasisFeatureValue(self):
    r"""
    Calculate and return the Low Intensity Distance Area Emphases (LISDE) value.

    :math:`LISDE = \frac{\sum^{N_g}_{i=1}\sum^{N_d}_{j=1}{\frac{\textbf{P}(i,j)}{i^2j^2}}}{\sum^{N_g}_{i=1}\sum^{N_d}_{j=1}{\textbf{P}(i,j)}}`

    Measures the proportion in the image of the joint distribution of smaller distance zones with lower gray-level values.
    """
    try:
      lisde = numpy.sum(
        (self.P_gldzm / ((self.coefficients['ivector'][:, None] ** 2) * (self.coefficients['jvector'][None, :] ** 2))),
        (0, 1)) / self.coefficients['sumP_gldzm']
    except ZeroDivisionError:
      lisde = numpy.core.numeric.NaN
    return (lisde)

  def getHighIntensitySmallDistanceEmphasisFeatureValue(self):
    r"""
    Calculate and return the High Intensity Distance Area Emphases (HISDE) value.

    :math:`HISDE = \frac{\sum^{N_g}_{i=1}\sum^{N_d}_{j=1}{\frac{\textbf{P}(i,j)i^2}{j^2}}}{\sum^{N_g}_{i=1}\sum^{N_d}_{j=1}{\textbf{P}(i,j)}}`

    Measures the proportion in the image of the joint distribution of smaller distance zones with higher gray-level values.
    """
    try:
      hisde = numpy.sum(
        (self.P_gldzm * (self.coefficients['ivector'][:, None] ** 2) / (self.coefficients['jvector'][None, :] ** 2)),
        (0, 1)) / self.coefficients['sumP_gldzm']
    except ZeroDivisionError:
      hisde = numpy.core.numeric.NaN
    return (hisde)

  def getLowIntensityLargeDistanceEmphasisFeatureValue(self):
    r"""
    Calculate and return the Low Intensity Large Distance Emphases (LILDE) value.

    :math:`LILDE = \frac{\sum^{N_g}_{i=1}\sum^{N_d}_{j=1}{\frac{\textbf{P}(i,j)j^2}{i^2}}}{\sum^{N_g}_{i=1}\sum^{N_d}_{j=1}{\textbf{P}(i,j)}}`

    Measures the proportion in the image of the joint distribution of larger distance zones with lower gray-level values.
    """
    try:
      lilde = numpy.sum(
        (self.P_gldzm * (self.coefficients['jvector'][None, :] ** 2) / (self.coefficients['ivector'][:, None] ** 2)),
        (0, 1)) / self.coefficients['sumP_gldzm']
    except ZeroDivisionError:
      lilde = numpy.core.numeric.NaN
    return (lilde)

  def getHighIntensityLargeDistanceEmphasisFeatureValue(self):
    r"""
    Calculate and return the High Intensity Large Area Emphases (HILDE) value.

    :math:`HILDE = \frac{\sum^{N_g}_{i=1}\sum^{N_d}_{j=1}{\textbf{P}(i,j)i^2j^2}}{\sum^{N_g}_{i=1}\sum^{N_d}_{j=1}{\textbf{P}(i,j)}}`

    Measures the proportion in the image of the joint distribution of larger distance zones with higher gray-level values.
    """
    try:
      hilde = numpy.sum(
        (self.P_gldzm * ((self.coefficients['jvector'][None, :] ** 2) * (self.coefficients['ivector'][:, None] ** 2))),
        (0, 1)) / self.coefficients['sumP_gldzm']
    except ZeroDivisionError:
      hilde = numpy.core.numeric.NaN
    return (hilde)
