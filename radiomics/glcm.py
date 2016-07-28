import numpy
import collections
from radiomics import base, imageoperations
import SimpleITK as sitk
from tqdm import  trange

class RadiomicsGLCM(base.RadiomicsFeaturesBase):
  """GLCM feature calculation."""
  def __init__(self, inputImage, inputMask, **kwargs):
    super(RadiomicsGLCM,self).__init__(inputImage, inputMask, **kwargs)

    if inputImage == None or inputMask == None:
      if self.verbose: print('ERROR GLCM: missing input image or mask')
      return

    self.imageArray = sitk.GetArrayFromImage(self.inputImage)
    self.maskArray = sitk.GetArrayFromImage(self.inputMask)

    (self.matrix, self.matrixCoordinates) = \
      imageoperations.padTumorMaskToCube(self.imageArray,self.maskArray)

    self.targetVoxelArray = self.matrix[self.matrixCoordinates]
    self.coefficients = {}
    self.P_glcm = {}

    # binning
    self.matrix, self.histogram = imageoperations.binImage(self.binWidth, self.targetVoxelArray, self.matrix, self.matrixCoordinates)
    self.coefficients['Ng'] = self.histogram[1].shape[0] - 1

    self.createGLCM()

    self.calculateCoefficients()

  def calculateGLCM(self, distances, angles):
    """
    Compute 13 GLCM matrices for the input image for every direction in 3D.

    (13 GLCMs for each neighboring voxel from a reference voxel centered in a 3x3 cube)
    """
    maxrows = self.matrix.shape[2]
    maxcols = self.matrix.shape[1]
    maxheight = self.matrix.shape[0]
    indices = zip(*self.matrixCoordinates)

    if self.verbose: bar = trange(len(indices), desc= 'calculate GLCM')

    for h, c, r in indices:
      if self.verbose: bar.update()

      for angles_idx, angle in enumerate(angles):
        for distances_idx, distance in enumerate(distances):
          i = self.matrix[h, c, r]
          i_idx = int(i-1)

          row = r + angle[2]
          col = c + angle[1]
          height = h + angle[0]

          if row >= 0 and row < maxrows and col >= 0 and col < maxcols:
            if tuple((height, col, row)) in indices:
              j = self.matrix[height, col, row]
              j_idx = int(j-1)
              self.P_glcm[i_idx, j_idx, distances_idx, angles_idx] += 1

    if self.verbose: bar.close()

  def createGLCM(self):
    """
    Generate container for GLCM Matrices: P_glcm as an array of shape: (i,j,gamma,a)

    i/j = total gray-level bins for image array,
    gamma = 1 (distance in voxels),
    a = 1...13 (directions in 3D)
    """
    Ng = self.coefficients['Ng']

    distances=numpy.array([1])
    angles = numpy.array([ (0, 1, 0),
                           (-1, 1, 0),
                           (-1, 0, 0),
                           (-1, -1, 0),
                           (0, 1, -1),
                           (0, 0, -1),
                           (0, -1, -1),
                           (-1, 0, -1),
                           (1, 0, -1),
                           (-1, 1, -1),
                           (1, -1, -1),
                           (-1, -1, -1),
                           (1, 1, -1) ])
    self.P_glcm = numpy.zeros( (Ng, Ng, distances.size, int(angles.shape[0])), dtype='float64' )
    self.calculateGLCM(distances, angles)

    # Normalize each glcm
    for q in xrange(int(angles.shape[0])):
      self.P_glcm[:,:,0,q] = self.P_glcm[:,:,0,q]/(self.P_glcm[:,:,0,q].sum())

  # check if ivector and jvector can be replaced
  def calculateCoefficients(self):
    """Calculate and fill in the coefficients dict"""
    Ng = self.coefficients['Ng']
    eps = numpy.spacing(1)


    NgVector = numpy.arange(1,self.P_glcm.shape[0]+1, dtype= 'float64')
    # shape = (Ng, Ng)
    i,j = numpy.meshgrid(NgVector, NgVector, indexing= 'ij')

    # shape = (2*Ng-1)
    kValuesSum = numpy.arange(2, (Ng*2)+1)
    # shape = (Ng-1)
    kValuesDiff = numpy.arange(0,Ng)

    # shape = (1, 1, distances.size, angles)
    u = self.P_glcm.mean((0, 1), keepdims= True)
    # marginal row probabilities #shape = (Ng, 1, distances.size, angles)
    px = self.P_glcm.sum(1, keepdims= True)
    # marginal column probabilities #shape = (1, Ng, distances.size, angles)
    py = self.P_glcm.sum(0, keepdims= True)

    # shape = (1, 1, distances.size, angles)
    ux = numpy.sum( i[:,:,None,None]*self.P_glcm, (0, 1), keepdims= True )
    uy = numpy.sum( j[:,:,None,None]*self.P_glcm, (0, 1), keepdims= True )

    # shape = (1, 1, distances.size, angles)
    sigx = numpy.sum(self.P_glcm*((i[:,:,None,None]-ux)**2), (0, 1), keepdims= True )**0.5
    # shape = (1, 1, distances.size, angles)
    sigy = numpy.sum(self.P_glcm*((j[:,:,None,None]-uy)**2), (0, 1), keepdims= True )**0.5

    # shape = (2*Ng-1, distances.size, angles)
    pxAddy = numpy.array([ numpy.sum(self.P_glcm[i+j == k], 0) for k in kValuesSum ])
    # shape = (Ng, distances.size, angles)
    pxSuby = numpy.array([ numpy.sum(self.P_glcm[numpy.abs(i-j) == k], 0) for k in kValuesDiff ])

    # entropy of px # shape = (distances.size, angles)
    HX = (-1) * numpy.sum( (px * numpy.log2(px+eps)), (0, 1))
    # entropy of py # shape = (distances.size, angles)
    HY = (-1) * numpy.sum( (py * numpy.log2(py+eps)), (0, 1))
    # shape = (distances.size, angles)
    HXY = (-1) * numpy.sum( (self.P_glcm * numpy.log2(self.P_glcm+eps)), (0, 1) )

    # shape = (distances.size, angles)
    HXY1 = (-1) * numpy.sum( (self.P_glcm * numpy.log2(px*py+eps)), (0, 1) )
    # shape = (distances.size, angles)
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
    """
    Using the i and j arrays, calculate and return the mean Autocorrelation for all 13 GLCMs.

    Autocorrelation is a measure of the magnitude of the
    fineness and coarseness of texture.
    """
    i = self.coefficients['i']
    j = self.coefficients['j']
    ac = numpy.sum( self.P_glcm * (i*j)[:,:,None,None], (0, 1) )
    return (ac.mean())

  def getClusterProminenceFeatureValue(self):
    """
    Using coefficients i, j, ux, uy, calculate and return the mean Cluster Prominence for all 13 GLCMs.

    Cluster Prominence is a measure of the skewness and asymmetry of the GLCM.
    A higher values implies more asymmetry about the mean while a lower value
    indicates a peak near the mean value and less variation about the mean.
    """
    i = self.coefficients['i']
    j = self.coefficients['j']
    ux = self.coefficients['ux']
    uy = self.coefficients['uy']
    cp = numpy.sum( (self.P_glcm * (((i+j)[:,:,None,None] - ux - uy)**4)), (0, 1) )
    return (cp.mean())

  def getClusterShadeFeatureValue(self):
    """
    Using coefficients i, j, ux, uy, calculate and return the mean Cluster Shade for all 13 GLCMs.

    Cluster Shade is a measure of the skewness and uniformity of the GLCM.
    A higher cluster shade implies greater asymmetry about the mean.
    """
    i = self.coefficients['i']
    j = self.coefficients['j']
    ux = self.coefficients['ux']
    uy = self.coefficients['uy']
    cs = numpy.sum( (self.P_glcm * (((i+j)[:,:,None,None] - ux - uy)**3)), (0, 1) )
    return (cs.mean())

  def getClusterTendencyFeatureValue(self):
    """
    Using coefficients i, j, ux, uy, calculate and return the mean Cluster Tendency for all 13 GLCMs.

    Cluster Tendency is a measure of groupings of voxels with similar gray-level values.
    """
    i = self.coefficients['i']
    j = self.coefficients['j']
    ux = self.coefficients['ux']
    uy = self.coefficients['uy']
    ct = numpy.sum( (self.P_glcm * (((i+j)[:,:,None,None] - ux - uy)**2)), (0, 1) )
    return (ct.mean())

  def getContrastFeatureValue(self):
    """
    Using coefficients i, j, calculate and return the mean Contrast for all 13 GLCMs.

    Contrast is a measure of the local intensity variation, favoring P(i,j)
    values away from the diagonal (i != j). A larger value correlates with
    a greater disparity in intensity values among neighboring voxels.
    """
    i = self.coefficients['i']
    j = self.coefficients['j']
    cont = numpy.sum( (self.P_glcm * ((numpy.abs(i-j))[:,:,None,None]**2)), (0, 1) )
    return (cont.mean())

  def getCorrelationFeatureValue(self):
    """
    Using coefficients i, j, ux, uy, sigx, sigy, calculate and return the mean Correlation value for all 13 GLCMs.

    Correlation is a value between 0 (uncorrelated) and 1 (perfectly correlated) showing the
    linear dependency of gray level values to their respective voxels in the GLCM.
    """
    i = self.coefficients['i']
    j = self.coefficients['j']
    ux = self.coefficients['ux']
    uy = self.coefficients['uy']
    sigx = self.coefficients['sigx']
    sigy = self.coefficients['sigy']

    try:
      corm = numpy.sum( self.P_glcm*(i[:,:,None,None]-ux)*(j[:,:,None,None]-uy), (0, 1), keepdims= True )
      corr = corm/(sigx*sigy)
      return (corr.mean())
    except ZeroDivisionError:
      return numpy.core.nan

  def getDifferenceEntropyFeatureValue(self):
    """
    Using coefficients pxSuby, eps, calculate and return the mean Difference Entropy for all 13 GLCMs.

    Difference Entropy is a measure of the randomness/variability
    in neighborhood intensity value differences.
    """
    pxSuby = self.coefficients['pxSuby']
    eps = self.coefficients['eps']
    difent = (-1) * numpy.sum( (pxSuby*numpy.log2(pxSuby+eps)), 0 )
    return (difent.mean())

  def getDissimilarityFeatureValue(self):
    """
    Using coefficients i, j, calculate and return the mean Dissimilarity for all 13 GLCMs.

    Dissimilarity is a measure of local intensity variation. A larger
    value correlates with a greater disparity in intensity values
    among neighboring voxels.
    """
    i = self.coefficients['i']
    j = self.coefficients['j']
    dis = numpy.sum( (self.P_glcm * (numpy.abs(i-j))[:,:,None,None]), (0, 1) )
    return (dis.mean())

  def getEnergyFeatureValue(self):
    """
    Using P_glcm, calculate and return the mean Energy for all 13 GLCMs.

    Energy (or Angular Second Moment)is a measure of homogeneous patterns
    in the image. A greater Energy implies that there are more instances
    of intensity value pairs in the image that neighbor each other at
    higher frequencies.
    """
    ene = numpy.sum( (self.P_glcm**2), (0, 1))
    return (ene.mean())

  def getEntropyFeatureValue(self):
    """
    Using coefficients eps, calculate and return the mean Entropy for all 13 GLCMs.

    Entropy is a measure of the randomness/variability in neighborhood intensity values.
    """
    ent = self.coefficients['HXY']
    return (ent.mean())

  def getHomogeneity1FeatureValue(self):
    """
    Using coefficients i, j, calculate and return the mean Homogeneity 1 for all 13 GLCMs.

    Homogeneity 1 is a measure of the similarity in intensity values for
    neighboring voxels. It is a measure of local homogeneity that increases
    with less contrast in the window.
    """
    i = self.coefficients['i']
    j = self.coefficients['j']
    homo1 = numpy.sum( (self.P_glcm / (1 + (numpy.abs(i-j))[:,:,None,None])), (0, 1) )
    return (homo1.mean())

  def getHomogeneity2FeatureValue(self):
    """
    Using coefficients i, j, calculate and return the mean Homogeneity 2 for all 13 GLCMs.

    Homogeneity 2 is a measure of the similarity in intensity values
    for neighboring voxels.
    """
    i = self.coefficients['i']
    j = self.coefficients['j']
    homo2 = numpy.sum( (self.P_glcm / (1 + (numpy.abs(i-j))[:,:,None,None]**2)), (0, 1) )
    return (homo2.mean())

  def getImc1FeatureValue(self):
    """
    Using coefficients HX, HY, HXY, HXY1, calculate and return the mean
    Informal Measure of Correlation 1 for all 13 GLCMs.
    """
    HX = self.coefficients['HX']
    HY = self.coefficients['HY']
    HXY = self.coefficients['HXY']
    HXY1 = self.coefficients['HXY1']
    imc1 = (HXY - HXY1)/numpy.max(([HX,HY]),0)
    return (imc1.mean())

  def getImc2FeatureValue(self):
    """
    Using coefficients HXY, HXY2, calculate and return the mean Informal Measure
    of Correlation 2 for all 13 GLCMs.
    """
    HXY = self.coefficients['HXY']
    HXY2 = self.coefficients['HXY2']

    imc2 = ( 1-numpy.e**(-2*(HXY2-HXY)) )**(0.5)    # matlab:(1-exp(-2*(hxy2-hxy)))^0.5;

    return (imc2.mean())

  def getIdmnFeatureValue(self):
    """
    Using coefficients i, j, Ng, calculate and return the mean IDMN for all 13 GLCMs.

    IDMN (inverse difference moment normalized)  is a measure of the local
    homogeneity of an image. IDMN weights are the inverse of the Contrast
    weights (decreasing exponentially from the diagonal i=j in the GLCM).
    Unlike Homogeneity2, IDMN normalizes the square of the difference between
    neighboring intensity values by dividing over the square of the total
    number of discrete intensity values.
    """
    i = self.coefficients['i']
    j = self.coefficients['j']
    Ng = self.coefficients['Ng']
    idmn = numpy.sum( (self.P_glcm / (1 + (((numpy.abs(i-j))[:,:,None,None]**2)/(Ng**2)))), (0, 1) )
    return (idmn.mean())

  def getIdnFeatureValue(self):
    """
    Using coefficients i, j, Ng, calculate and return the mean IDN for all 13 GLCMs.

    IDN (inverse difference normalized) is another measure of the local
    homogeneity of an image. Unlike Homogeneity1, IDN normalizes the difference
    between the neighboring intensity values by dividing over the total number
    of discrete intensity values.
    """
    i = self.coefficients['i']
    j = self.coefficients['j']
    Ng = self.coefficients['Ng']
    idn  = numpy.sum( (self.P_glcm / (1 + ((numpy.abs(i-j))[:,:,None,None]/Ng))), (0, 1) )
    return (idn.mean())

  def getInverseVarianceFeatureValue(self):
    """Using the i, j coeffients, calculate and return the mean Inverse Variance for all 13 GLCMs."""
    i = self.coefficients['i']
    j = self.coefficients['j']
    maskDiags = numpy.abs(i-j) > 0
    inv = numpy.sum( (self.P_glcm[maskDiags] / ((numpy.abs(i-j))[:,:,None,None]**2)[maskDiags]), 0 )
    return (inv.mean())

  def getMaximumProbabilityFeatureValue(self):
    """
    Using P_glcm, calculate and return the mean Maximum Probability for all 13 GLCMs.

    Maximum Probability is occurrences of the most predominant pair of
    neighboring intensity values.
    """
    maxprob = self.P_glcm.max((0, 1))
    return (maxprob.mean())

  def getSumAverageFeatureValue(self):
    """
    Using coefficients pxAddy, kValuesSum, calculate and return the mean Sum Average for all 13 GLCMs.

    Sum Average measures the relationship between occurrences of pairs
    with lower intensity values and occurrences of pairs with higher intensity
    values.
    """
    pxAddy = self.coefficients['pxAddy']
    kValuesSum = self.coefficients['kValuesSum']
    sumavg =  numpy.sum( (kValuesSum[:,None,None]*pxAddy), 0 )
    return (sumavg.mean())

  def getSumEntropyFeatureValue(self):
    """
    Using coefficients pxAddy, eps, calculate and return the mean Sum Entropy for all 13 GLCMs.

    Sum Entropy is a sum of neighborhood intensity value differences.
    """
    pxAddy = self.coefficients['pxAddy']
    eps = self.coefficients['eps']
    sumentr = (-1) * numpy.sum( (pxAddy * numpy.log2(pxAddy+eps)), 0 )
    return (sumentr.mean())

  def getSumVarianceFeatureValue(self):
    """
    Using coefficients pxAddy, kValuesSum, calculate and return the mean Sum Variance for all 13 GLCMs.

    Sum Variance is a measure of heterogeneity that places higher weights on
    neighboring intensity level pairs that deviate more from the mean.
    """
    eps = self.coefficients['eps']
    pxAddy = self.coefficients['pxAddy']
    kValuesSum = self.coefficients['kValuesSum']
    sumentr = (-1) * numpy.sum( (pxAddy * numpy.log2(pxAddy+eps)), 0, keepdims= True )
    sumvar = numpy.sum( (pxAddy*((kValuesSum[:,None,None] - sumentr)**2)), 0 )
    return (sumvar.mean())

  def getSumSquaresFeatureValue(self):
    """
    Using coefficients j, u, calculate and return the mean um of Squares (Variance) for all 13 GLCMs.

    Sum of Squares or Variance is a measure in the distribution of neigboring
    intensity level pairs about the mean intensity level in the GLCM.
    """
    j = self.coefficients['j']
    u = self.coefficients['u']
    # Also known as Variance
    ss = numpy.sum( (self.P_glcm * ((j[:,:,None,None]-u)**2)), (0, 1) )
    return (ss.mean())
