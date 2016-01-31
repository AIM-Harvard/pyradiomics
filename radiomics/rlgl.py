import numpy
import SimpleITK as sitk
from radiomics import base, imageoperations


class RadiomicsRLGL(base.RadiomicsFeaturesBase):
  """RLGL feature calculation."""
  def __init__(self, inputImage, inputMask, **kwargs):
      super(RadiomicsRLGL,self).__init__(inputImage, inputMask, **kwargs)

      self.imageArray = sitk.GetArrayFromImage(inputImage)
      self.maskArray = sitk.GetArrayFromImage(inputMask)

      (self.matrix, self.matrixCoordinates) = \
        imageoperations.padTumorMaskToCube(self.imageArray,self.maskArray)

      self.targetVoxelArray = self.matrix[self.matrixCoordinates]
      self.coefficients = {}
      self.P_rlgl = {}

      # binning
      self.matrix, self.histogram = imageoperations.binImage(self.binWidth, self.targetVoxelArray, self.matrix, self.matrixCoordinates)
      self.coefficients['Ng'] = len(self.histogram[0])
      self.coefficients['grayLevels'] = numpy.linspace(1,self.coefficients['Ng'],num=self.coefficients['Ng'])
      self.coefficients['Nr'] = numpy.max(self.matrix.shape)
      self.coefficients['Np'] = self.targetVoxelArray.size

      self.calculateRLGL(13)

      self.calculateCoefficients()

  def calculateRLGL(self, angles):
      Ng = self.coefficients['Ng']
      Nr = self.coefficients['Nr']
      grayLevels = self.coefficients['grayLevels']

      self.P_rlgl = numpy.zeros((Ng, Nr, angles))

      padVal = -1000   #use eps or NaN to pad matrix
      padMask = numpy.zeros(self.matrix.shape,dtype=bool)
      padMask[self.matrixCoordinates] = True
      self.matrix[~padMask] = padVal
      matrixDiagonals = list()

      #(1,0,0), #(-1,0,0),
      aDiags = reduce(lambda x,y: x+y, [a.tolist() for a in numpy.transpose(self.matrix,(1,2,0))])
      matrixDiagonals.append( filter(lambda x: numpy.nonzero(x)[0].size>1, aDiags) )

      #(0,1,0), #(0,-1,0),
      bDiags = reduce(lambda x,y: x+y, [a.tolist() for a in numpy.transpose(self.matrix,(0,2,1))])
      matrixDiagonals.append( filter(lambda x: numpy.nonzero(x)[0].size>1, bDiags) )

      #(0,0,1), #(0,0,-1),
      cDiags = reduce(lambda x,y: x+y, [a.tolist() for a in numpy.transpose(self.matrix,(0,1,2))])
      matrixDiagonals.append( filter(lambda x: numpy.nonzero(x)[0].size>1, cDiags) )

      #(1,1,0),#(-1,-1,0),
      lowBound = -self.matrix.shape[0]+1
      highBound = self.matrix.shape[1]

      dDiags = reduce(lambda x,y: x+y, [self.matrix.diagonal(a,0,1).tolist() for a in xrange(lowBound, highBound)])
      matrixDiagonals.append( filter(lambda x: numpy.nonzero(x)[0].size>1, dDiags) )

      #(1,0,1), #(-1,0-1),
      lowBound = -self.matrix.shape[0]+1
      highBound = self.matrix.shape[2]

      eDiags = reduce(lambda x,y: x+y, [self.matrix.diagonal(a,0,2).tolist() for a in xrange(lowBound, highBound)])
      matrixDiagonals.append( filter(lambda x: numpy.nonzero(x)[0].size>1, eDiags) )

      #(0,1,1), #(0,-1,-1),
      lowBound = -self.matrix.shape[1]+1
      highBound = self.matrix.shape[2]

      fDiags = reduce(lambda x,y: x+y, [self.matrix.diagonal(a,1,2).tolist() for a in xrange(lowBound, highBound)])
      matrixDiagonals.append( filter(lambda x: numpy.nonzero(x)[0].size>1, fDiags) )

      #(1,-1,0), #(-1,1,0),
      lowBound = -self.matrix.shape[0]+1
      highBound = self.matrix.shape[1]

      gDiags = reduce(lambda x,y: x+y, [self.matrix[:,::-1,:].diagonal(a,0,1).tolist() for a in xrange(lowBound, highBound)])
      matrixDiagonals.append( filter(lambda x: numpy.nonzero(x)[0].size>1, gDiags) )

      #(-1,0,1), #(1,0,-1),
      lowBound = -self.matrix.shape[0]+1
      highBound = self.matrix.shape[2]

      hDiags = reduce(lambda x,y: x+y, [self.matrix[:,:,::-1].diagonal(a,0,2).tolist() for a in xrange(lowBound, highBound)])
      matrixDiagonals.append( filter(lambda x: numpy.nonzero(x)[0].size>1, hDiags) )

      #(0,1,-1), #(0,-1,1),
      lowBound = -self.matrix.shape[1]+1
      highBound = self.matrix.shape[2]

      iDiags = reduce(lambda x,y: x+y, [self.matrix[:,:,::-1].diagonal(a,1,2).tolist() for a in xrange(lowBound, highBound)])
      matrixDiagonals.append( filter(lambda x: numpy.nonzero(x)[0].size>1, iDiags) )

      #(1,1,1), #(-1,-1,-1)
      lowBound = -self.matrix.shape[0]+1
      highBound = self.matrix.shape[1]

      jDiags = [ numpy.diagonal(h,x,0,1).tolist() for h in [self.matrix.diagonal(a,0,1) for a in xrange(lowBound, highBound)] for x in xrange(-h.shape[0]+1, h.shape[1]) ]
      matrixDiagonals.append( filter(lambda x: numpy.nonzero(x)[0].size>1, jDiags) )

      #(-1,1,-1), #(1,-1,1),
      lowBound = -self.matrix.shape[0]+1
      highBound = self.matrix.shape[1]

      kDiags = [ numpy.diagonal(h,x,0,1).tolist() for h in [self.matrix[:,::-1,:].diagonal(a,0,1) for a in xrange(lowBound, highBound)] for x in xrange(-h.shape[0]+1, h.shape[1]) ]
      matrixDiagonals.append( filter(lambda x: numpy.nonzero(x)[0].size>1, kDiags) )

      #(1,1,-1), #(-1,-1,1),
      lowBound = -self.matrix.shape[0]+1
      highBound = self.matrix.shape[1]

      lDiags = [ numpy.diagonal(h,x,0,1).tolist() for h in [self.matrix[:,:,::-1].diagonal(a,0,1) for a in xrange(lowBound, highBound)] for x in xrange(-h.shape[0]+1, h.shape[1]) ]
      matrixDiagonals.append( filter(lambda x: numpy.nonzero(x)[0].size>1, lDiags) )

      #(-1,1,1), #(1,-1,-1),
      lowBound = -self.matrix.shape[0]+1
      highBound = self.matrix.shape[1]

      mDiags = [ numpy.diagonal(h,x,0,1).tolist() for h in [self.matrix[:,::-1,::-1].diagonal(a,0,1) for a in xrange(lowBound, highBound)] for x in xrange(-h.shape[0]+1, h.shape[1]) ]
      matrixDiagonals.append( filter(lambda x: numpy.nonzero(x)[0].size>1, mDiags) )

      # Run-Length Encoding (rle) for the 13 list of diagonals
      # (1 list per 3D direction/angle)
      for angle in xrange (len(matrixDiagonals)):
        P = self.P_rlgl[:,:,angle]
        for diagonal in matrixDiagonals[angle]:
          diagonal = numpy.array(diagonal, dtype='int')
          pos, = numpy.where(numpy.diff(diagonal) != 0)
          pos = numpy.concatenate(([0], pos+1, [len(diagonal)]))

          rle = zip([n for n in diagonal[pos[:-1]]], pos[1:] - pos[:-1])
          for level, run_length in rle:
            if level in grayLevels:
              level_idx = level-1
              run_length_idx = run_length-1
              P[level_idx, run_length_idx] += 1


  def calculateCoefficients(self):
      sumP_rlgl = numpy.sum( numpy.sum(self.P_rlgl, 0), 0 )
      sumP_rlgl[sumP_rlgl==0] = 1

      pr = numpy.sum(self.P_rlgl, 0)
      pg = numpy.sum(self.P_rlgl, 1)

      ivector = numpy.arange(1, self.P_rlgl.shape[0] + 1)
      jvector = numpy.arange(1, self.P_rlgl.shape[1] + 1)

      self.coefficients['sumP_rlgl'] = sumP_rlgl
      self.coefficients['pr'] = pr
      self.coefficients['pg'] = pg
      self.coefficients['ivector'] = ivector
      self.coefficients['jvector'] = jvector

  def getShortRunEmphasisFeatureValue(self):
    pr = self.coefficients['pr']
    jvector = self.coefficients['jvector']
    sumP_rlgl = self.coefficients['sumP_rlgl']

    try:
      sre = numpy.sum( (pr/(jvector[:,None]**2)), 0 ) / sumP_rlgl
    except ZeroDivisionError:
      sre = 0
    return (sre.mean())

  def getLongRunEmphasisFeatureValue(self):
    pr =  self.coefficients['pr']
    jvector = self.coefficients['jvector']
    sumP_rlgl = self.coefficients['sumP_rlgl']

    try:
      lre = numpy.sum( (pr*(jvector[:,None]**2)), 0 ) / sumP_rlgl
    except ZeroDivisionError:
      lre = 0
    return (lre.mean())

  def getGrayLevelNonUniformityFeatureValue(self):
    pg = self.coefficients['pg']
    sumP_rlgl = self.coefficients['sumP_rlgl']

    try:
      gln = numpy.sum( (pg**2) , 0 ) / sumP_rlgl
    except ZeroDivisionError:
      gln = 0
    return (gln.mean())

  def getRunLengthNonUniformityFeatureValue(self):
    pr = self.coefficients['pr']
    sumP_rlgl = self.coefficients['sumP_rlgl']

    try:
      rln = numpy.sum( (pr**2) , 0 ) / sumP_rlgl
    except ZeroDivisionError:
      rln = 0
    return (rln.mean())

  def getRunPercentageFeatureValue(self):
    Np = self.coefficients['Np']

    try:
      rp = numpy.sum( numpy.sum( (self.P_rlgl/(Np)) , 0 ), 0 )
    except ZeroDivisionError:
      rp = 0
    return (rp.mean())

  def getLowGrayLevelRunEmphasisFeatureValue(self):
    pg = self.coefficients['pg']
    ivector = self.coefficients['ivector']
    sumP_rlgl = self.coefficients['sumP_rlgl']

    try:
      lglre = numpy.sum( (pg/(ivector[:,None]**2)) , 0 ) / sumP_rlgl
    except ZeroDivisionError:
      lglre = 0
    return (lglre.mean())

  def getHighGrayLevelRunEmphasisFeatureValue(self):
    pg = self.coefficients['pg']
    ivector = self.coefficients['ivector']
    sumP_rlgl = self.coefficients['sumP_rlgl']

    try:
      hglre = numpy.sum( (pg*(ivector[:,None]**2)) , 0 ) / sumP_rlgl
    except ZeroDivisionError:
      hglre = 0
    return (hglre.mean())

  def getShortRunLowGrayLevelEmphasisFeatureValue(self):
    ivector = self.coefficients['ivector']
    jvector = self.coefficients['jvector']
    sumP_rlgl = self.coefficients['sumP_rlgl']

    try:
      srlgle = numpy.sum( numpy.sum( (self.P_rlgl/((ivector[:,None,None]**2)*(jvector[None,:,None]**2))) , 0 ), 0 ) / sumP_rlgl
    except ZeroDivisionError:
      srlgle = 0
    return (srlgle.mean())

  def getShortRunHighGrayLevelEmphasisFeatureValue(self):
    ivector = self.coefficients['ivector']
    jvector = self.coefficients['jvector']
    sumP_rlgl = self.coefficients['sumP_rlgl']

    try:
      srhgle = numpy.sum( numpy.sum( (self.P_rlgl*(ivector[:,None,None]**2)/(jvector[None,:,None]**2)) , 0 ), 0 ) / sumP_rlgl
    except ZeroDivisionError:
      srhgle = 0
    return (srhgle.mean())

  def getLongRunLowGrayLevelEmphasisFeatureValue(self):
    ivector = self.coefficients['ivector']
    jvector = self.coefficients['jvector']
    sumP_rlgl = self.coefficients['sumP_rlgl']

    try:
      lrlgle = numpy.sum( numpy.sum( (self.P_rlgl*(jvector[None,:,None]**2)/(ivector[:,None,None]**2)) , 0 ), 0 ) / sumP_rlgl
    except ZeroDivisionError:
      lrlgle = 0
    return (lrlgle.mean())

  def getLongRunHighGrayLevelEmphasisFeatureValue(self):
    ivector = self.coefficients['ivector']
    jvector = self.coefficients['jvector']
    sumP_rlgl = self.coefficients['sumP_rlgl']

    try:
      lrhgle = numpy.sum( numpy.sum( (self.P_rlgl*((jvector[None,:,None]**2)*(ivector[:,None,None]**2))) , 0 ), 0 ) / sumP_rlgl
    except ZeroDivisionError:
      lrhgle = 0
    return (lrhgle.mean())
