from itertools import chain
import numpy
import SimpleITK as sitk
from radiomics import base, imageoperations


class RadiomicsNGTDM(base.RadiomicsFeaturesBase):
  r"""
  A Neighbouring Gray Tone Difference Matrix quantifies the difference between a gray value and the average gray value
  of its neighbours. The sum of absolute differences for gray level i is stored in the matrix
  Amadasun and King (1989)
  """

  def __init__(self, inputImage, inputMask, **kwargs):
    super(RadiomicsNGTDM, self).__init__(inputImage, inputMask, **kwargs)

    self.gldm_a = kwargs.get('gldm_a', 0)

    self.coefficients = {}

    # binning
    self.matrix, self.histogram = imageoperations.binImage(self.binWidth, self.targetVoxelArray, self.matrix,
                                                           self.matrixCoordinates)
    self.coefficients['Ng'] = self.histogram[1].shape[0] - 1
    self.coefficients['Np'] = self.targetVoxelArray.size

    self._calculateGLDM()

  def _calculateGLDM(self):
    Ng = self.coefficients['Ng']

    self.matrix = self.matrix.astype('float')

    # Set voxels outside delineation to padding value
    padVal = numpy.nan
    self.matrix[(self.maskArray != self.label)] = padVal

    size = numpy.max(self.matrixCoordinates, 1) - numpy.min(self.matrixCoordinates, 1) + 1
    angles = imageoperations.generateAngles(size)
    angles = numpy.concatenate((angles, angles * -1))
    Nd = len(angles)

    angMats = numpy.empty((self.matrix.shape + (Nd,)), dtype='int')

    for a_idx, a in enumerate(angles):
      # create shifted array (by angle), so that for an index idx, angMat[idx] is the neigbour of self.matrix[idx]
      # for the current angle.
      angMat = angMats[:, :, :, a_idx]
      angMat[:, :, :] = numpy.roll(numpy.roll(numpy.roll(self.matrix, -a[0], 0), -a[1], 1), -a[2], 2)
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

    # Create a difference matrix by taking the absolut of self.matrix - the neighbourhood mean, excluding NaNs
    # in the calculation of the mean
    diffMat = numpy.abs(self.matrix - numpy.nanmean(angMats, axis=3))

    # For each gray level present in self.matrix:
    # element 0 = probability of gray level (p_i),
    # element 1 = sum of the absolute differences (s_i),
    # element 2 = gray level (i)
    self.P_ngtdm = numpy.zeros((Ng, 3), dtype='float')
    for i in numpy.unique(self.matrix):
        if not numpy.isnan(i):
          i_ind = numpy.where(self.matrix == i)
          self.P_ngtdm[int(i-1), 0] = len(i_ind[0])
          self.P_ngtdm[int(i-1), 1] = numpy.sum(diffMat[i_ind])

    # Fill in gray levels (needed as empty gray level slices will be deleted)
    self.P_ngtdm[:, 2] = numpy.arange(1, Ng+1)

    # Delete empty grey levels
    self.P_ngtdm = numpy.delete(self.P_ngtdm, numpy.where(self.P_ngtdm[:, 0] == 0), 0)

    # Normalize P_ngtdm[:, 0] (= p_i)
    self.P_ngtdm[:, 0] /= self.coefficients['Np']

    # Ngp = number of graylevels, for which p_i > 0
    self.coefficients['Ngp'] = self.P_ngtdm.shape[0]

  def getCoarsenessFeatureValue(self):
    r"""
    Calculate and return the coarseness.


    :math:`Coarseness = \frac{1}{/sum^{N_g}_{i=1}{p_{i}s_{i}}}`

    N.B. :math:`/sum^{N_g}_{i=1}{p_{i}s_{i}}` potentially evaluates to 0 (in case of a completely homogeneous image).
    If this is the case, an arbitrary value of :math:`10^6` is returned.
    """
    p_i = self.P_ngtdm[:, 0]
    s_i = self.P_ngtdm[:, 1]
    sum_coarse = numpy.sum(p_i * s_i)
    if sum_coarse == 0:
      return 1000000
    else:
      return 1 / sum_coarse

  def getContrastFeatureValue(self):
    r"""
    Calculate and return the contrast.


    :math:`Contrast = \big(\frac{1}{N_{g,p}(N_{g,p}-1)\displaystyle\sum^{N_g}_{i=1}\displaystyle\sum^{N_g}_{j=1}{p_{i}p_{j}(i-j)^2}\big)\big(\frac{1}{N_p^2}\displaystyle\sum^{N_g}_{i=1}{s_i}\big)\text{ where }p_i \neq 0, p_j \neq 0`
    """
    Ngp = self.coefficients['Ngp']
    Np = self.coefficients['Np']
    p_i = self.P_ngtdm[:, 0]
    s_i = self.P_ngtdm[:, 1]
    i = self.P_ngtdm[:, 2]
    contrast = (numpy.sum(p_i[:, None] * p_i[None, :] * (i[:, None] - i[None, :])**2)/(Ngp*(Ngp-1))) * \
               ((numpy.sum(s_i))/(Np**2))

    return contrast

  def getBusynessFeatureValue(self):
    r"""
    Calculate and return the busyness.


    :math:`Busyness = \frac{\sum^{N_g}_{i = 1}{p_{i}s_{i}}{\sum^{N_g}_{i = 1}\sum^{N_g}_{j = 1}{|ip_i - jp_j|}}\text{ where }p_i \neq 0, p_j \neq 0`
    """
    p_i = self.P_ngtdm[:, 0]
    s_i = self.P_ngtdm[:, 1]
    i = self.P_ngtdm[:, 2]
    busyness = numpy.sum(p_i * s_i) / numpy.sum(numpy.abs((i*p_i)[:, None] * (i*p_i)[None, :]))
    return busyness

  def getComplexityFeatureValue(self):
    r"""
    Calculate and return the complexity.


    :math:`Complexity = \frac{1}{N_p^2}\displaystyle\sum^{N_g}_{i = 1}\displaystyle\sum^{N_g}_{j = 1}{|i - j|\frac{p_{i}s_{i} + p_{j}s_{j}}{p_i + p_j}}\text{ where }p_i \neq 0, p_j \neq 0`
    """
    Np = self.coefficients['Np']
    p_i = self.P_ngtdm[:, 0]
    s_i = self.P_ngtdm[:, 1]
    i = self.P_ngtdm[:, 2]
    complexity = numpy.sum(numpy.abs(i[:, None] - i[None, :]) * (((p_i*s_i)[:, None] + (p_i*s_i)[None, :])/(p_i[:, None] + p_i[None, :]))) / Np**2
    return complexity

  def getStrengthFeatureValue(self):
    r"""
    Calculate and return the strength.


    :math:`Strength = \frac{\sum^{N_g}_{i = 1}\sum^{N_g}_{j = 1}{(p_i + p_j)(i-j)^2}}{\sum^{N_g}_{i = 1}{s_i}}\text{ where }p_i \neq 0, p_j \neq 0`

    N.B. :math:`/sum^{N_g}_{i=1}{s_i}` potentially evaluates to 0 (in case of a completely homogeneous image).
    If this is the case, 0 is returned.
    """
    p_i = self.P_ngtdm[:, 0]
    s_i = self.P_ngtdm[:, 1]
    i = self.P_ngtdm[:, 2]
    sum_s_i = numpy.sum(s_i)
    if sum_s_i == 0:
      return 0
    else:
      strength = numpy.sum((p_i[:, None] + p_i[None, :]) * (i[:, None] - i[None, :])**2)/sum_s_i
      return strength
