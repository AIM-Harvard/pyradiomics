import numpy
from tqdm import trange
from radiomics import base, imageoperations

class RadiomicsDWGLCM(base.RadiomicsFeaturesBase):
  """
  Distance weighted Gray Level Co-occurrence Matrix (DWGLCM) feature calculation.

  Compute 1 GLCM matrix for the input image all directions in 3D.
  Contributions of voxels can be weighted by weigthing factor W, which is calculated from distance d by:

  :math:`W = e^{-\|d\|^2}`

  The weighting function to use for d can be specified by setting weightingNorm.
  The following settings are allowed:
  'infinity': Infinity Norm
  'euclidean': Euclidean Norm
  'manhattan': Manhattan (or TaxiCab) Norm
  None / Other value: No weighting (W = 1).

  default is "infinity".
  """
  def __init__(self, inputImage, inputMask, **kwargs):
    super(RadiomicsDWGLCM,self).__init__(inputImage, inputMask, **kwargs)

    self.pixelSpacing = self.inputImage.GetSpacing()
    self.symmetricalGLCM = kwargs.get('symmetricalGLCM', True)
    self.weightingNorm = kwargs.get('weightingNorm', 'infinity')  # manhattan, euclidean, infinity

    self.coefficients = {}
    self.P_glcm = {}

    # binning
    self.matrix, self.histogram = imageoperations.binImage(self.binWidth, self.targetVoxelArray, self.matrix, self.matrixCoordinates)
    self.coefficients['Ng'] = self.histogram[1].shape[0] - 1

    self._calculateDWGLCM()
    self.calculateCoefficients()

  def _calculateDWGLCM(self):
    Ng = self.coefficients['Ng']

    # Exclude voxels outside segmentation, due to binning, no negative values will be encountered inside the mask
    self.matrix[self.maskArray != self.label] = -1

    size = numpy.max(self.matrixCoordinates, 1) - numpy.min(self.matrixCoordinates, 1) + 1
    angles = imageoperations.generateAngles(size)

    self.P_glcm = numpy.zeros( (Ng, Ng), dtype='float64' )

    if self.verbose: bar = trange(len(angles), desc= 'calculate GLCM')

    weights = numpy.empty(len(angles))
    for a_idx, a in enumerate(angles):
      if self.verbose: bar.update()

      if self.weightingNorm == 'infinity':
        weights[a_idx] = numpy.exp(-max(numpy.abs(a))**2)
      elif self.weightingNorm == 'euclidean':
        weights[a_idx] = numpy.exp(-numpy.sum((numpy.abs(a)*self.pixelSpacing[::-1])**2))
      elif self.weightingNorm == 'manhattan':
        weights[a_idx] = numpy.exp(-numpy.sum(numpy.abs(a)*self.pixelSpacing[::-1])**2)
      else:
        weights[a_idx] = 1

    if self.verbose: bar = trange(Ng, desc='calculate GLCM')

    # iterate over gray levels for center voxel
    for i in xrange(1, Ng + 1):
      # give some progress
      if self.verbose: bar.update()

      # get the indices to all voxels which have the current gray level i
      i_indices = numpy.where(self.matrix == i)

      # iterate over gray levels for neighbouring voxel
      for j in xrange(1, Ng + 1):
        # get the indices to all voxels which have the current gray level j
        j_indices = set(zip(*numpy.where(self.matrix == j)))

        for a_idx, a in enumerate(angles):
          # get the corresponding indices of the neighbours for angle a
          neighbour_indices = set(zip(*(i_indices + a[:, None])))

          # The following intersection yields the indices to voxels with gray level j
          # that are also a neighbour of a voxel with gray level i for angle a.
          # The number of indices is then equal to the total number of pairs with gray level i and j for angle a
          count = len(neighbour_indices.intersection(j_indices))
          self.P_glcm[i - 1, j - 1] += count * weights[a_idx]

    if self.verbose: bar.close()

    # Optionally make DWGLCM symmetrical for each angle
    if self.symmetricalGLCM:
      self.P_glcm += numpy.transpose(self.P_glcm, (1, 0))

    # Normalize each glcm
    sumGlcm = numpy.sum(self.P_glcm, (0,1))
    self.P_glcm = self.P_glcm/sumGlcm

  # check if ivector and jvector can be replaced
  def calculateCoefficients(self):
    r"""
    Calculate and fill in the coefficients dict.
    :math:`p_{x+y} = \displaystyle\sum^{N_g}_{i=1}\displaystyle\sum^{N_g}_{j=1}{\textbf{P}(i,j)},i+j=k,k=2,3,\dots,2N_g`

    :math:`p_{x-y} = \displaystyle\sum^{N_g}_{i=1}\displaystyle\sum^{N_g}_{j=1}{\textbf{P}(i,j)},|i-j|=k,k=0,1,\dots,N_g-1`

    :math:`HX =  -\displaystyle\sum^{N_g}_{i=1}{p_x(i)\log_2\big(p_x(i)\big)}`

    :math:`HY =  -\displaystyle\sum^{N_g}_{j=1}{p_y(j)\log_2\big(p_y(j)\big)}`

    :math:`HXY =  -\displaystyle\sum^{N_g}_{i=1}\displaystyle\sum^{N_g}_{j=1}{\textbf{P}(i,j)\log_2\big(\textbf{P}(i,j)\big)}`

    :math:`HXY1 =  -\displaystyle\sum^{N_g}_{i=1}\displaystyle\sum^{N_g}_{j=1}{\textbf{P}(i,j)\log_2\big(p_x(i)p_y(j)\big)}`

    :math:`HXY2 =  -\displaystyle\sum^{N_g}_{i=1}\displaystyle\sum^{N_g}_{j=1}{p_x(i)p_y(j)\log_2\big(p_x(i)p_y(j)\big)}`
    """
    Ng = self.coefficients['Ng']
    eps = numpy.spacing(1)


    NgVector = numpy.arange(1,self.P_glcm.shape[0]+1, dtype= 'float64')
    # shape = (Ng, Ng)
    i, j = numpy.meshgrid(NgVector, NgVector, indexing= 'ij')

    # shape = (2*Ng-1)
    kValuesSum = numpy.arange(2, (Ng*2)+1)
    # shape = (Ng-1)
    kValuesDiff = numpy.arange(0,Ng)

    # shape = (1, 1)
    u = numpy.sum( i*self.P_glcm, (0, 1), keepdims= True )
    # marginal row probabilities #shape = (Ng, 1)
    px = self.P_glcm.sum(1, keepdims= True)
    # marginal column probabilities #shape = (1, Ng)
    py = self.P_glcm.sum(0, keepdims= True)

    # shape = (1, 1)
    ux = numpy.sum( i*self.P_glcm, (0, 1), keepdims= True )
    uy = numpy.sum( j*self.P_glcm, (0, 1), keepdims= True )

    # shape = (1, 1)
    sigx = numpy.sum(self.P_glcm*((i-ux)**2), (0, 1), keepdims= True )**0.5
    # shape = (1, 1)
    sigy = numpy.sum(self.P_glcm*((j-uy)**2), (0, 1), keepdims= True )**0.5

    # shape = (2*Ng-1)
    pxAddy = numpy.array([ numpy.sum(self.P_glcm[i+j == k], 0) for k in kValuesSum ])
    # shape = (Ng)
    pxSuby = numpy.array([ numpy.sum(self.P_glcm[numpy.abs(i-j) == k], 0) for k in kValuesDiff ])

    # entropy of px # shape = (1)
    HX = (-1) * numpy.sum( (px * numpy.log2(px+eps)), (0, 1))
    # entropy of py # shape = (1)
    HY = (-1) * numpy.sum( (py * numpy.log2(py+eps)), (0, 1))
    # shape = (1)
    HXY = (-1) * numpy.sum( (self.P_glcm * numpy.log2(self.P_glcm+eps)), (0, 1) )

    # shape = (1)
    HXY1 = (-1) * numpy.sum( (self.P_glcm * numpy.log2(px*py+eps)), (0, 1) )
    # shape = (1)
    HXY2 = (-1) * numpy.sum( ((px*py) * numpy.log2(px*py+eps)), (0, 1) )

    self.coefficients['Ng'] = Ng
    self.coefficients['eps'] = eps
    self.coefficients['i'] = i
    self.coefficients['j'] = j
    self.coefficients['kValuesSum'] = kValuesSum
    self.coefficients['kValuesDiff'] = kValuesDiff
    self.coefficients['u'] = u
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
    Using the i and j arrays, calculate and return the Autocorrelation.

    :math:`autocorrelation = \displaystyle\sum^{N_g}_{i=1}\displaystyle\sum^{N_g}_{j=1}{ij\textbf{P}(i,j)}`

    Autocorrelation is a measure of the magnitude of the
    fineness and coarseness of texture.
    """
    i = self.coefficients['i']
    j = self.coefficients['j']
    ac = numpy.sum( self.P_glcm * (i*j), (0, 1) )
    return ac

  def getClusterProminenceFeatureValue(self):
    r"""
    Using coefficients i, j, ux, uy, calculate and return the Cluster Prominence.

    :math:`cluster\ prominence = \displaystyle\sum^{N_g}_{i=1}\displaystyle\sum^{N_g}_{j=1}{\big( i+j-\mu_x(i)-\mu_y(j)\big)^4\textbf{P}(i,j)}`

    Cluster Prominence is a measure of the skewness and asymmetry of the DWGLCM.
    A higher values implies more asymmetry about the mean while a lower value
    indicates a peak near the mean value and less variation about the mean.
    """
    i = self.coefficients['i']
    j = self.coefficients['j']
    ux = self.coefficients['ux']
    uy = self.coefficients['uy']
    cp = numpy.sum( (self.P_glcm * (((i+j) - ux - uy)**4)), (0, 1) )
    return cp

  def getClusterShadeFeatureValue(self):
    r"""
    Using coefficients i, j, ux, uy, calculate and return the Cluster Shade.

    :math:`cluster\ shade = \displaystyle\sum^{N_g}_{i=1}\displaystyle\sum^{N_g}_{j=1}{\big(i+j-\mu_x(i)-\mu_y(j)\big)^3\textbf{P}(i,j)}`

    Cluster Shade is a measure of the skewness and uniformity of the DWGLCM.
    A higher cluster shade implies greater asymmetry about the mean.
    """
    i = self.coefficients['i']
    j = self.coefficients['j']
    ux = self.coefficients['ux']
    uy = self.coefficients['uy']
    cs = numpy.sum( (self.P_glcm * (((i+j) - ux - uy)**3)), (0, 1) )
    return cs

  def getClusterTendencyFeatureValue(self):
    r"""
    Using coefficients i, j, ux, uy, calculate and return the Cluster Tendency.

    :math:`cluster\ prominence = \displaystyle\sum^{N_g}_{i=1}\displaystyle\sum^{N_g}_{j=1}{\big(i+j-\mu_x(i)-\mu_y(j)\big)^2\textbf{P}(i,j)}`

    Cluster Tendency is a measure of groupings of voxels with similar gray-level values.
    """
    i = self.coefficients['i']
    j = self.coefficients['j']
    ux = self.coefficients['ux']
    uy = self.coefficients['uy']
    ct = numpy.sum( (self.P_glcm * (((i+j) - ux - uy)**2)), (0, 1) )
    return ct

  def getContrastFeatureValue(self):
    r"""
    Using coefficients i, j, calculate and return the Contrast.

    :math:`contrast = \displaystyle\sum^{N_g}_{i=1}\displaystyle\sum^{N_g}_{j=1}{(i-j)^2\textbf{P}(i,j)}`

    Contrast is a measure of the local intensity variation, favoring P(i,j)
    values away from the diagonal (i != j). A larger value correlates with
    a greater disparity in intensity values among neighboring voxels.
    """
    i = self.coefficients['i']
    j = self.coefficients['j']
    cont = numpy.sum( (self.P_glcm * (numpy.abs(i-j)**2)), (0, 1) )
    return cont

  def getCorrelationFeatureValue(self):
    r"""
    Using coefficients i, j, ux, uy, sigx, sigy, calculate and return the Correlation value.

    :math:`correlation = \frac{\sum^{N_g}_{i=1}\sum^{N_g}_{j=1}{ij\textbf{P}(i,j)-\mu_x(i)\mu_y(j)}}{\sigma_x(i)\sigma_y(j)}`

    Correlation is a value between 0 (uncorrelated) and 1 (perfectly correlated) showing the
    linear dependency of gray level values to their respective voxels in the DWGLCM.
    """
    i = self.coefficients['i']
    j = self.coefficients['j']
    ux = self.coefficients['ux']
    uy = self.coefficients['uy']
    sigx = self.coefficients['sigx']
    sigy = self.coefficients['sigy']

    try:
      corm = numpy.sum( self.P_glcm*(i-ux)*(j-uy), (0, 1), keepdims= True )
      corr = corm/(sigx*sigy)
      return corr.mean()
    except ZeroDivisionError:
      return numpy.core.nan

  def getDifferenceEntropyFeatureValue(self):
    r"""
    Using coefficients pxSuby, eps, calculate and return the Difference Entropy.

    :math:`difference\ entropy = \displaystyle\sum^{N_g-1}_{i=0}{\textbf{P}_{x-y}(i)\log_2\big(\textbf{P}_{x-y}(i)\big)}`

    Difference Entropy is a measure of the randomness/variability
    in neighborhood intensity value differences.
    """
    pxSuby = self.coefficients['pxSuby']
    eps = self.coefficients['eps']
    difent = (-1) * numpy.sum( (pxSuby*numpy.log2(pxSuby+eps)), 0 )
    return difent

  def getDissimilarityFeatureValue(self):
    r"""
    Using coefficients i, j, calculate and return the Dissimilarity.

    :math:`dissimilarity = \displaystyle\sum^{N_g}_{i=1}\displaystyle\sum^{N_g}_{j=1}{|i-j|\textbf{P}(i,j)}`

    Dissimilarity is a measure of local intensity variation. A larger
    value correlates with a greater disparity in intensity values
    among neighboring voxels.
    """
    i = self.coefficients['i']
    j = self.coefficients['j']
    dis = numpy.sum( (self.P_glcm * numpy.abs(i-j)), (0, 1) )
    return dis

  def getEnergyFeatureValue(self):
    r"""
    Using P_glcm, calculate and return the Energy.

    :math:`energy = \displaystyle\sum^{N_g}_{i=1}\displaystyle\sum^{N_g}_{j=1}{\big(\textbf{P}(i,j)\big)^2}`

    Energy (or Angular Second Moment)is a measure of homogeneous patterns
    in the image. A greater Energy implies that there are more instances
    of intensity value pairs in the image that neighbor each other at
    higher frequencies.
    """
    ene = numpy.sum( (self.P_glcm**2), (0, 1))
    return ene

  def getEntropyFeatureValue(self):
    r"""
    Using coefficients eps, calculate and return the Entropy.

    :math:`entropy = -\displaystyle\sum^{N_g}_{i=1}\displaystyle\sum^{N_g}_{j=1}{\textbf{P}(i,j)\log_2\big(\textbf{P}(i,j)\big)}`

    Entropy is a measure of the randomness/variability in neighborhood intensity values.
    """
    ent = self.coefficients['HXY']
    return ent

  def getHomogeneity1FeatureValue(self):
    r"""
    Using coefficients i, j, calculate and return the Homogeneity 1.

    :math:`homogeneity\ 1 = \displaystyle\sum^{N_g}_{i=1}\displaystyle\sum^{N_g}_{j=1}{\frac{\textbf{P}(i,j)}{1+|i-j|}}`

    Homogeneity 1 is a measure of the similarity in intensity values for
    neighboring voxels. It is a measure of local homogeneity that increases
    with less contrast in the window.
    """
    i = self.coefficients['i']
    j = self.coefficients['j']
    homo1 = numpy.sum( (self.P_glcm / (1 + numpy.abs(i-j))), (0, 1) )
    return homo1

  def getHomogeneity2FeatureValue(self):
    r"""
    Using coefficients i, j, calculate and return the Homogeneity 2.

    :math:`homogeneity\ 2 = \displaystyle\sum^{N_g}_{i=1}\displaystyle\sum^{N_g}_{j=1}{\frac{\textbf{P}(i,j)}{1+|i-j|^2}}`

    Homogeneity 2 is a measure of the similarity in intensity values
    for neighboring voxels.
    """
    i = self.coefficients['i']
    j = self.coefficients['j']
    homo2 = numpy.sum( (self.P_glcm / (1 + numpy.abs(i-j)**2)), (0, 1) )
    return homo2

  def getImc1FeatureValue(self):
    r"""
    Using coefficients HX, HY, HXY, HXY1, calculate and return the
    Informal Measure of Correlation 1.

    :math:`IMC\ 1 = \frac{HXY-HXY1}{\max\{HX,HY\}}`
    """
    HX = self.coefficients['HX']
    HY = self.coefficients['HY']
    HXY = self.coefficients['HXY']
    HXY1 = self.coefficients['HXY1']
    imc1 = (HXY - HXY1)/numpy.max(([HX,HY]),0)
    return imc1

  def getImc2FeatureValue(self):
    r"""
    Using coefficients HXY, HXY2, calculate and return the Informal Measure
    of Correlation 2.

    :math:`IMC\ 2 = \sqrt{1-e^{-2(HXY2-HXY)}}`
    """
    HXY = self.coefficients['HXY']
    HXY2 = self.coefficients['HXY2']

    imc2 = ( 1-numpy.e**(-2*(HXY2-HXY)) )**(0.5)    # matlab:(1-exp(-2*(hxy2-hxy)))^0.5;

    return imc2

  def getIdmnFeatureValue(self):
    r"""
    Using coefficients i, j, Ng, calculate and return the IDMN.

    :math:`IDMN = \displaystyle\sum^{N_g}_{i=1}\displaystyle\sum^{N_g}_{j=1}{ \frac{\textbf{P}(i,j)}{1+\left(\frac{|i-j|^2}{N_g^2}\right)} }`

    IDMN (inverse difference moment normalized)  is a measure of the local
    homogeneity of an image. IDMN weights are the inverse of the Contrast
    weights (decreasing exponentially from the diagonal i=j in the DWGLCM).
    Unlike Homogeneity2, IDMN normalizes the square of the difference between
    neighboring intensity values by dividing over the square of the total
    number of discrete intensity values.
    """
    i = self.coefficients['i']
    j = self.coefficients['j']
    Ng = self.coefficients['Ng']
    idmn = numpy.sum( (self.P_glcm / (1 + (((numpy.abs(i-j))**2)/(Ng**2)))), (0, 1) )
    return idmn

  def getIdnFeatureValue(self):
    r"""
    Using coefficients i, j, Ng, calculate and return the IDN.

    :math:`IDN = \displaystyle\sum^{N_g}_{i=1}\displaystyle\sum^{N_g}_{j=1}{ \frac{\textbf{P}(i,j)}{1+\left(\frac{|i-j|}{N_g}\right)} }`

    IDN (inverse difference normalized) is another measure of the local
    homogeneity of an image. Unlike Homogeneity1, IDN normalizes the difference
    between the neighboring intensity values by dividing over the total number
    of discrete intensity values.
    """
    i = self.coefficients['i']
    j = self.coefficients['j']
    Ng = self.coefficients['Ng']
    idn  = numpy.sum( (self.P_glcm / (1 + ((numpy.abs(i-j))/Ng))), (0, 1) )
    return idn

  def getInverseVarianceFeatureValue(self):
    r"""Using the i, j coeffients, calculate and return the Inverse Variance.

    :math:`inverse\ variance = \displaystyle\sum^{N_g}_{i=1}\displaystyle\sum^{N_g}_{j=1}{\frac{\textbf{P}(i,j)}{|i-j|^2}}, i \neq j`
    """
    i = self.coefficients['i']
    j = self.coefficients['j']
    maskDiags = numpy.abs(i-j) > 0
    inv = numpy.sum( (self.P_glcm[maskDiags] / ((numpy.abs(i-j))**2)[maskDiags]), 0 )
    return inv

  def getMaximumProbabilityFeatureValue(self):
    r"""
    Using P_glcm, calculate and return the Maximum Probability.

    :math:`maximum\ probability = \max\big(\textbf{P}(i,j)\big)`

    Maximum Probability is occurrences of the most predominant pair of
    neighboring intensity values.
    """
    maxprob = self.P_glcm.max((0, 1))
    return maxprob

  def getSumAverageFeatureValue(self):
    r"""
    Using coefficients pxAddy, kValuesSum, calculate and return the Sum Average.

    :math:`sum\ average = \displaystyle\sum^{2N_g}_{i=2}{i\textbf{P}_{x+y}(i)}`

    Sum Average measures the relationship between occurrences of pairs
    with lower intensity values and occurrences of pairs with higher intensity
    values.
    """
    pxAddy = self.coefficients['pxAddy']
    kValuesSum = self.coefficients['kValuesSum']
    sumavg =  numpy.sum( (kValuesSum*pxAddy), 0 )
    return sumavg

  def getSumEntropyFeatureValue(self):
    r"""
    Using coefficients pxAddy, eps, calculate and return the Sum Entropy.

    :math:`sum\ entropy = \displaystyle\sum^{2N_g}_{i=2}{\textbf{P}_{x+y}(i)\log_2\big(\textbf{P}_{x+y}(i)\big)}`

    Sum Entropy is a sum of neighborhood intensity value differences.
    """
    pxAddy = self.coefficients['pxAddy']
    eps = self.coefficients['eps']
    sumentr = (-1) * numpy.sum( (pxAddy * numpy.log2(pxAddy+eps)), 0 )
    return sumentr

  def getSumVarianceFeatureValue(self):
    r"""
    Using coefficients pxAddy, kValuesSum, SumAvarage calculate and return the Sum Variance.

    :math:`sum\ variance\ 2 = \displaystyle\sum^{2N_g}_{i=2}{(1-SA)^2\textbf{P}_{x+y}(i)}`

    Sum Variance is a measure of heterogeneity that places higher weights on
    neighboring intensity level pairs that deviate more from the mean.

    This formula differs from SumVariance in that instead of subtracting the SumEntropy from the intensity,
    it subtracts the SumAvarage, which is the mean of intensities and not its entropy
    """
    pxAddy = self.coefficients['pxAddy']
    kValuesSum = self.coefficients['kValuesSum']
    sumavg = self.getSumAverageFeatureValue()
    sumvar = numpy.sum( (pxAddy*((kValuesSum - sumavg)**2)), 0 )
    return sumvar

  def getSumSquaresFeatureValue(self):
    r"""
    Using coefficients i, u, calculate and return the sum of Squares (Variance).

    :math:`sum\ squares = \displaystyle\sum^{N_g}_{i=1}\displaystyle\sum^{N_g}_{j=1}{(i-\mu)^2\textbf{P}(i,j)}`

    Sum of Squares or Variance is a measure in the distribution of neigboring
    intensity level pairs about the mean intensity level in the DWGLCM.
    """
    i = self.coefficients['j']
    u = self.coefficients['u']
    # Also known as Variance
    ss = numpy.sum( (self.P_glcm * ((i-u)**2)), (0, 1) )
    return ss