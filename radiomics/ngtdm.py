import numpy

from radiomics import base, cMatrices


class RadiomicsNGTDM(base.RadiomicsFeaturesBase):
  r"""
  A Neighbouring Gray Tone Difference Matrix quantifies the difference between a gray value and the average gray value
  of its neighbours within distance :math:`\delta`. The sum of absolute differences for gray level :math:`i` is stored in the matrix.
  Let :math:`\textbf{X}_{gl}` be a set of segmented voxels and :math:`x_{gl}(j_x,j_y,j_z) \in \textbf{X}_{gl}` be the gray level of a voxel at postion
  :math:`(j_x,j_y,j_z)`, then the average gray level of the neigbourhood is:

  .. math::

    \bar{A}_i &= \bar{A}(j_x, j_y, j_z) \\
    &= \displaystyle\frac{1}{W} \displaystyle\sum_{k_x=-\delta}^{\delta}\displaystyle\sum_{k_y=-\delta}^{\delta}
    \displaystyle\sum_{k_z=-\delta}^{\delta}{x_{gl}(j_x+k_x, j_y+k_y, j_z+k_z)}, \\
    &\mbox{where }(k_x,k_y,k_z)\neq(0,0,0)\mbox{ and } x_{gl}(j_x+k_x, j_y+k_y, j_z+k_z) \in \textbf{X}_{gl}

  Here, :math:`W` is the number of voxels in the neighbourhood that are also in :math:`\textbf{X}_{gl}`.

  As a two dimensional example, let the following matrix :math:`\textbf{I}` represent a 4x4 image,
  having 5 discrete grey levels, but no voxels with gray level 4:

  .. math::
    \textbf{I} = \begin{bmatrix}
    1 & 2 & 5 & 2\\
    3 & 5 & 1 & 3\\
    1 & 3 & 5 & 5\\
    3 & 1 & 1 & 1\end{bmatrix}

  The following NGTDM can be obtained:

  .. math::
    \begin{array}{cccc}
    i & n_i & p_i & s_i\\
    \hline
    1 & 6 & 0.375 & 13.35\\
    2 & 2 & 0.125 & 2.00\\
    3 & 4 & 0.25  & 2.63\\
    4 & 0 & 0.00  & 0.00\\
    5 & 4 & 0.25  & 10.075\end{array}

  6 pixels have gray level 1, therefore:

  :math:`s_1 = |1-10/3| + |1-30/8| + |1-15/5| + |1-13/5| + |1-15/5| + |1-11/3| = 13.35`

  For gray level 2, there are 2 pixels, therefore:

  :math:`s_2 = |2-15/5| + |2-15/5| = 2`

  Similar for gray values 3 and 5:

  :math:`s_3 = |3-12/5| + |3-18/5| + |3-20/8| + |3-5/3| = 3.03 \\
  s_5 = |5-14/5| + |5-18/5| + |5-20/8| + |5-11/5| = 10.075`

  Let:

  :math:`n_i` be the number of voxels in :math:`X_{gl}` with gray level :math:`i`

  :math:`N_{v,p}` be the total number of voxels in :math:`X_{gl}` and equal to :math:`\sum{n_i}` (i.e. the number of voxels
  with a valid region; at least 1 neighbor). :math:`N_{v,p} \leq N_p`, where :math:`N_p` is the total number of voxels in the ROI.

  :math:`p_i` be the gray level probability and equal to :math:`n_i/N_v`

  :math:`s_i = \left\{ {\begin{array} {rcl}
  \sum^{n_i}{|i-\bar{A}_i|} & \mbox{for} & n_i \neq 0 \\
  0 & \mbox{for} & n_i = 0 \end{array}}\right.`
  be the sum of absolute differences for gray level :math:`i`

  :math:`N_g` be the number of discreet gray levels

  :math:`N_{g,p}` be the number of gray levels where :math:`p_i \neq 0`

  The following class specific settings are possible:

  - distances [[1]]: List of integers. This specifies the distances between the center voxel and the neighbor, for which
    angles should be generated.

  References

  - Amadasun M, King R; Textural features corresponding to textural properties;
    Systems, Man and Cybernetics, IEEE Transactions on 19:1264-1274 (1989). doi: 10.1109/21.44046
  """

  def __init__(self, inputImage, inputMask, **kwargs):
    super(RadiomicsNGTDM, self).__init__(inputImage, inputMask, **kwargs)

    self.P_ngtdm = None
    self._applyBinning()

  def _initCalculation(self):
    self.coefficients['Np'] = len(self.labelledVoxelCoordinates[0])
    self.P_ngtdm = self._calculateMatrix()
    self._calculateCoefficients()

  def _calculateMatrix(self):
    P_ngtdm = cMatrices.calculate_ngtdm(self.matrix,
                                        self.maskArray,
                                        numpy.array(self.settings.get('distances', [1])),
                                        self.coefficients['Ng'],
                                        self.settings.get('force2D', False),
                                        self.settings.get('force2Ddimension', 0))

    # Delete empty grey levels
    P_ngtdm = numpy.delete(P_ngtdm, numpy.where(P_ngtdm[:, 0] == 0), 0)

    return P_ngtdm

  def _calculateCoefficients(self):
    Nvp = numpy.sum(self.P_ngtdm[:, 0])  # = No of voxels that have a valid region, lesser equal to Np
    self.coefficients['Nvp'] = Nvp
    if Nvp < self.coefficients['Np']:
      self.logger.debug('Detected %d voxels without valid neighbors', self.coefficients['Np'] - Nvp)

    # Normalize P_ngtdm[:, 0] (= n_i) to obtain p_i
    self.coefficients['p_i'] = self.P_ngtdm[:, 0] / Nvp

    self.coefficients['s_i'] = self.P_ngtdm[:, 1]
    self.coefficients['ivector'] = self.P_ngtdm[:, 2]

    # Ngp = number of graylevels, for which p_i > 0
    self.coefficients['Ngp'] = self.P_ngtdm.shape[0]

  def getCoarsenessFeatureValue(self):
    r"""
    Calculate and return the coarseness.

    :math:`Coarseness = \frac{1}{\sum^{N_g}_{i=1}{p_{i}s_{i}}}`

    Coarseness is a measure of average difference between the center voxel and its neighbourhood and is an indication
    of the spatial rate of change. A higher value indicates a lower spatial change rate and a locally more uniform texture.

    N.B. :math:`\sum^{N_g}_{i=1}{p_{i}s_{i}}` potentially evaluates to 0 (in case of a completely homogeneous image).
    If this is the case, an arbitrary value of :math:`10^6` is returned.
    """
    p_i = self.coefficients['p_i']
    s_i = self.coefficients['s_i']
    sum_coarse = numpy.sum(p_i * s_i)
    if sum_coarse == 0:
      return 1000000
    else:
      return 1 / sum_coarse

  def getContrastFeatureValue(self):
    r"""
    Calculate and return the contrast.

    :math:`Contrast = \left(\frac{1}{N_{g,p}(N_{g,p}-1)}\displaystyle\sum^{N_g}_{i=1}\displaystyle\sum^{N_g}_{j=1}{p_{i}p_{j}(i-j)^2}\right)
    \left(\frac{1}{N_{v,p}}\displaystyle\sum^{N_g}_{i=1}{s_i}\right)\text{, where }p_i \neq 0, p_j \neq 0`

    Contrast is a measure of the spatial intensity change, but is also dependent on the overall gray level dynamic range.
    Contrast is high when both the dynamic range and the spatial change rate are high, i.e. an image with a large range
    of gray levels, with large changes between voxels and their neighbourhood.

    N.B. In case of a completely homogeneous image, :math:`N_{g,p} = 1`, which would result in a division by 0. In this
    case, an arbitray value of 0 is returned.
    """
    Ngp = self.coefficients['Ngp']
    Nvp = self.coefficients['Nvp']
    p_i = self.coefficients['p_i']
    s_i = self.coefficients['s_i']
    i = self.coefficients['ivector']

    if Ngp <= 1:
      return 0

    contrast = (numpy.sum(p_i[:, None] * p_i[None, :] * (i[:, None] - i[None, :]) ** 2) / (Ngp * (Ngp - 1))) * \
               ((numpy.sum(s_i)) / Nvp)

    return contrast

  def getBusynessFeatureValue(self):
    r"""
    Calculate and return the busyness.

    :math:`Busyness = \frac{\sum^{N_g}_{i = 1}{p_{i}s_{i}}}{\sum^{N_g}_{i = 1}\sum^{N_g}_{j = 1}{|ip_i - jp_j|}}\text{, where }p_i \neq 0, p_j \neq 0`

    A measure of the change from a pixel to its neighbour. A high value for busyness indicates a 'busy' image, with rapid
    changes of intensity between pixels and its neighbourhood.

    N.B. if :math:`N_{g,p} = 1`, then :math:`busyness = \frac{0}{0}`. If this is the case, 0 is returned, as it concerns
    a fully homogeneous region.
    """
    p_i = self.coefficients['p_i']
    s_i = self.coefficients['s_i']
    i = self.coefficients['ivector']
    absdiff = numpy.sum(numpy.abs((i * p_i)[:, None] - (i * p_i)[None, :]))
    if absdiff == 0:
      return 0
    else:
      busyness = numpy.sum(p_i * s_i) / absdiff
    return busyness

  def getComplexityFeatureValue(self):
    r"""
    Calculate and return the complexity.

    :math:`Complexity = \frac{1}{N_{v,p}}\displaystyle\sum^{N_g}_{i = 1}\displaystyle\sum^{N_g}_{j = 1}{|i - j|
    \frac{p_{i}s_{i} + p_{j}s_{j}}{p_i + p_j}}\text{, where }p_i \neq 0, p_j \neq 0`

    An image is considered complex when there are many primitive components in the image, i.e. the image is non-uniform
    and there are many rapid changes in gray level intensity.
    """
    Nvp = self.coefficients['Nvp']
    p_i = self.coefficients['p_i']
    s_i = self.coefficients['s_i']
    i = self.coefficients['ivector']
    complexity = numpy.sum(numpy.abs(i[:, None] - i[None, :]) * (
      ((p_i * s_i)[:, None] + (p_i * s_i)[None, :]) / (p_i[:, None] + p_i[None, :]))) / Nvp
    return complexity

  def getStrengthFeatureValue(self):
    r"""
    Calculate and return the strength.

    :math:`Strength = \frac{\sum^{N_g}_{i = 1}\sum^{N_g}_{j = 1}{(p_i + p_j)(i-j)^2}}{\sum^{N_g}_{i = 1}{s_i}}\text{, where }p_i \neq 0, p_j \neq 0`

    Strenght is a measure of the primitives in an image. Its value is high when the primitives are easily defined and
    visible, i.e. an image with slow change in intensity but more large coarse differences in gray level intensities.

    N.B. :math:`\sum^{N_g}_{i=1}{s_i}` potentially evaluates to 0 (in case of a completely homogeneous image).
    If this is the case, 0 is returned.
    """
    p_i = self.coefficients['p_i']
    s_i = self.coefficients['s_i']
    i = self.coefficients['ivector']
    sum_s_i = numpy.sum(s_i)
    if sum_s_i == 0:
      return 0
    else:
      strength = numpy.sum((p_i[:, None] + p_i[None, :]) * (i[:, None] - i[None, :]) ** 2) / sum_s_i
      return strength
