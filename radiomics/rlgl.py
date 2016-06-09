from itertools import chain
import numpy
import SimpleITK as sitk
from radiomics import base, imageoperations

class RadiomicsRLGL(base.RadiomicsFeaturesBase):
  """RLGL feature calculation."""
  def __init__(self, inputImage, inputMask, **kwargs):
      super(RadiomicsRLGL,self).__init__(inputImage, inputMask, **kwargs)

      if inputImage == None or inputMask == None:
        if self.verbose: print('ERROR RLGL: missing input image or mask')
        return

      self.imageArray = sitk.GetArrayFromImage(self.inputImage)
      self.maskArray = sitk.GetArrayFromImage(self.inputMask)

      (self.matrix, self.matrixCoordinates) = \
        imageoperations.padTumorMaskToCube(self.imageArray,self.maskArray)

      self.targetVoxelArray = self.matrix[self.matrixCoordinates]
      self.coefficients = {}
      self.P_rlgl = {}

      # binning
      self.matrix, self.histogram = imageoperations.binImage(self.binWidth, self.targetVoxelArray, self.matrix, self.matrixCoordinates)
      self.coefficients['Ng'] = self.histogram[1].shape[0] - 1
      self.coefficients['grayLevels'] = numpy.linspace(1,self.coefficients['Ng'],num=self.coefficients['Ng'])
      self.coefficients['Nr'] = numpy.max(self.matrix.shape)
      self.coefficients['Np'] = self.targetVoxelArray.size

      self.calculateRLGL()
      self.calculateCoefficients()

  def calculateRLGL(self):
      Ng = self.coefficients['Ng']
      Nr = self.coefficients['Nr']
      grayLevels = self.coefficients['grayLevels']

      P_rlgl = numpy.zeros((Ng, Nr, 13))

      padVal = -2000   #use eps or NaN to pad matrix
      padMask = numpy.zeros(self.matrix.shape,dtype=bool)
      padMask[self.matrixCoordinates] = True
      self.matrix[~padMask] = padVal
      
      matrixDiagonals = []

      #(1,0,0), (-1,0,0)
      aDiags = chain.from_iterable(numpy.transpose(self.matrix,(1,2,0)))
      matrixDiagonals.append( filter(lambda x: numpy.nonzero(x != padVal)[0].size>0, aDiags) )

      #(0,1,0), (0,-1,0)
      bDiags = chain.from_iterable(numpy.transpose(self.matrix,(0,2,1)))
      matrixDiagonals.append( filter(lambda x: numpy.nonzero(x != padVal)[0].size>0, bDiags) )

      #(0,0,1), (0,0,-1)
      cDiags = chain.from_iterable(numpy.transpose(self.matrix,(0,1,2)))
      matrixDiagonals.append( filter(lambda x: numpy.nonzero(x != padVal)[0].size>0, cDiags) )

      #(1,1,0), (-1,-1,0)
      lowBound = -self.matrix.shape[0]+1
      highBound = self.matrix.shape[1]
      dDiags = chain.from_iterable([self.matrix.diagonal(a,0,1) for a in xrange(lowBound, highBound)])
      matrixDiagonals.append( filter(lambda x: numpy.nonzero(x != padVal)[0].size>0, dDiags) )

      #(1,0,1), (-1,0-1)
      lowBound = -self.matrix.shape[0]+1
      highBound = self.matrix.shape[2]
      eDiags = chain.from_iterable([self.matrix.diagonal(a,0,2) for a in xrange(lowBound, highBound)])
      matrixDiagonals.append( filter(lambda x: numpy.nonzero(x != padVal)[0].size>0, eDiags) )

      #(0,1,1), (0,-1,-1)
      lowBound = -self.matrix.shape[1]+1
      highBound = self.matrix.shape[2]
      fDiags = chain.from_iterable([self.matrix.diagonal(a,1,2) for a in xrange(lowBound, highBound)])
      matrixDiagonals.append( filter(lambda x: numpy.nonzero(x != padVal)[0].size>0, fDiags) )

      #(1,-1,0), (-1,1,0)
      lowBound = -self.matrix.shape[0]+1
      highBound = self.matrix.shape[1]
      gDiags = chain.from_iterable([self.matrix[:,::-1,:].diagonal(a,0,1) for a in xrange(lowBound, highBound)])
      matrixDiagonals.append( filter(lambda x: numpy.nonzero(x != padVal)[0].size>0, gDiags) )

      #(-1,0,1), (1,0,-1)
      lowBound = -self.matrix.shape[0]+1
      highBound = self.matrix.shape[2]
      hDiags = chain.from_iterable([self.matrix[:,:,::-1].diagonal(a,0,2) for a in xrange(lowBound, highBound)])
      matrixDiagonals.append( filter(lambda x: numpy.nonzero(x != padVal)[0].size>0, hDiags) )

      #(0,1,-1), (0,-1,1)
      lowBound = -self.matrix.shape[1]+1
      highBound = self.matrix.shape[2]
      iDiags = chain.from_iterable([self.matrix[:,:,::-1].diagonal(a,1,2) for a in xrange(lowBound, highBound)])
      matrixDiagonals.append( filter(lambda x: numpy.nonzero(x != padVal)[0].size>0, iDiags) )

      #(1,1,1), (-1,-1,-1)
      lowBound = -self.matrix.shape[0]+1
      highBound = self.matrix.shape[1]
      jDiags = []
      for h in [self.matrix.diagonal(a,0,1) for a in xrange(lowBound, highBound)]:
        for x in xrange(-h.shape[0]+1, h.shape[1]):
          jDiags.append(numpy.diagonal(h,x,0,1))    
      matrixDiagonals.append( filter(lambda x: numpy.nonzero(x != padVal)[0].size>0, jDiags) )

      #(-1,1,-1), #(1,-1,1),
      lowBound = -self.matrix.shape[0]+1
      highBound = self.matrix.shape[1]
      kDiags = []
      for h in [self.matrix[:,::-1,:].diagonal(a,0,1) for a in xrange(lowBound, highBound)]:
        for x in xrange(-h.shape[0]+1, h.shape[1]):
          kDiags.append(numpy.diagonal(h,x,0,1))   
      matrixDiagonals.append( filter(lambda x: numpy.nonzero(x != padVal)[0].size>0, kDiags) )

      #(1,1,-1), #(-1,-1,1),
      lowBound = -self.matrix.shape[0]+1
      highBound = self.matrix.shape[1]
      lDiags = []
      for h in [self.matrix[:,:,::-1].diagonal(a,0,1) for a in xrange(lowBound, highBound)]:
        for x in xrange(-h.shape[0]+1, h.shape[1]):
          lDiags.append(numpy.diagonal(h,x,0,1))    
      matrixDiagonals.append( filter(lambda x: numpy.nonzero(x != padVal)[0].size>0, lDiags) )

      #(-1,1,1), #(1,-1,-1),
      lowBound = -self.matrix.shape[0]+1
      highBound = self.matrix.shape[1]
      mDiags = []
      for h in [self.matrix[:,::-1,::-1].diagonal(a,0,1) for a in xrange(lowBound, highBound)]:
        for x in xrange(-h.shape[0]+1, h.shape[1]):
          mDiags.append(numpy.diagonal(h,x,0,1))    
      matrixDiagonals.append( filter(lambda x: numpy.nonzero(x != padVal)[0].size>0, mDiags) )
      
      # Run-Length Encoding (rle) for the 13 list of diagonals
      # (1 list per 3D direction/angle)
      for angle_idx, angle in enumerate(matrixDiagonals):
        P = P_rlgl[:,:,angle_idx]
        for diagonal in angle:
          pos, = numpy.where(numpy.diff(diagonal) != 0)
          pos = numpy.concatenate(([0], pos+1, [len(diagonal)]))
          rle = zip([int(n) for n in diagonal[pos[:-1]]], pos[1:] - pos[:-1])
          for level, run_length in rle:
            if level != padVal:
              P[level-1, run_length-1] += 1
      
      # Crop gray-level axis of RLGL matrix to between minimum and maximum observed gray-levels
      # Crop run-length axis of RLGL matrix up to maximum observed run-length
      P_rlgl_bounds = numpy.argwhere(P_rlgl)
      (xstart, ystart, zstart), (xstop, ystop, zstop) = P_rlgl_bounds.min(0), P_rlgl_bounds.max(0) + 1
      self.P_rlgl = P_rlgl[xstart:xstop,:ystop,:]
      
  def calculateCoefficients(self):
      sumP_rlgl = numpy.sum( numpy.sum(self.P_rlgl, 0), 0 )

      pr = numpy.sum(self.P_rlgl, 0)
      pg = numpy.sum(self.P_rlgl, 1)

      ivector = numpy.arange(1, self.P_rlgl.shape[0] + 1, dtype=numpy.float64)
      jvector = numpy.arange(1, self.P_rlgl.shape[1] + 1, dtype=numpy.float64)

      self.coefficients['sumP_rlgl'] = sumP_rlgl
      self.coefficients['pr'] = pr
      self.coefficients['pg'] = pg
      self.coefficients['ivector'] = ivector
      self.coefficients['jvector'] = jvector

  def getShortRunEmphasisFeatureValue(self):
    """Calculate and return the mean Short Run Emphasis (SRE) value for all 13 RLGL matrices.
    
    A measure of the distribution of short run lengths, with a greater value indicative 
    of shorter run lengths and more fine textural textures.
    """
    pr = self.coefficients['pr']
    jvector = self.coefficients['jvector']
    sumP_rlgl = self.coefficients['sumP_rlgl']

    try:
      sre = numpy.sum( (pr/(jvector[:,None]**2)), 0 ) / sumP_rlgl
    except ZeroDivisionError:
      sre = 0
    return (sre.mean())

  def getLongRunEmphasisFeatureValue(self):
    """Calculate and return the mean Long Run Emphasis (LRE) value for all 13 RLGL matrices.
    
    A measure of the distribution of long run lengths, with a greater value indicative 
    of longer run lengths and more coarse structural textures.
    """  
    pr =  self.coefficients['pr']
    jvector = self.coefficients['jvector']
    sumP_rlgl = self.coefficients['sumP_rlgl']

    try:
      lre = numpy.sum( (pr*(jvector[:,None]**2)), 0 ) / sumP_rlgl
    except ZeroDivisionError:
      lre = 0
    return (lre.mean())

  def getGrayLevelNonUniformityFeatureValue(self):
    """Calculate and return the mean Gray Level Non-Uniformity (GLN) value for all 13 RLGL matrices.
    
    Measures the similarity of gray-level intensity values in the image, where a lower GLN value 
    correlates with a greater similarity in intensity values.
    """
    pg = self.coefficients['pg']
    sumP_rlgl = self.coefficients['sumP_rlgl']

    try:
      gln = numpy.sum( (pg**2) , 0 ) / sumP_rlgl
    except ZeroDivisionError:
      gln = 0
    return (gln.mean())

  def getRunLengthNonUniformityFeatureValue(self):
    """Calculate and return the mean Run Length Non-Uniformity (RLN) value for all 13 RLGL matrices.
    
    Measures the similarity of run lengths throughout the image, with a lower value indicating 
    more homogeneity among run lengths in the image.
    """
    pr = self.coefficients['pr']
    sumP_rlgl = self.coefficients['sumP_rlgl']

    try:
      rln = numpy.sum( (pr**2) , 0 ) / sumP_rlgl
    except ZeroDivisionError:
      rln = 0
    return (rln.mean())

  def getRunPercentageFeatureValue(self):
    """Calculate and return the mean Run Percentage (RP) value for all 13 RLGL matrices.
    
    Measures the homogeneity and distribution of runs of an image for a certain direction.
    """
    Np = self.coefficients['Np']

    try:
      rp = numpy.sum( (self.P_rlgl/(Np)) , (0, 1) )
    except ZeroDivisionError:
      rp = 0
    return (rp.mean())

  def getLowGrayLevelRunEmphasisFeatureValue(self):
    """Calculate and return the mean Low Gray Level Run Emphasis (LGLRE) value for all 13 RLGL matrices.
    
    Measures the distribution of low gray-level values, with a higher value indicating a greater 
    concentration of low gray-level values in the image.
    """    
    pg = self.coefficients['pg']
    ivector = self.coefficients['ivector']
    sumP_rlgl = self.coefficients['sumP_rlgl']

    try:
      lglre = numpy.sum( (pg/(ivector[:,None]**2)) , 0 ) / sumP_rlgl
    except ZeroDivisionError:
      lglre = 0
    return (lglre.mean())

  def getHighGrayLevelRunEmphasisFeatureValue(self):
    """Calculate and return the mean High Gray Level Run Emphasis (HGLRE) value for all 13 RLGL matrices.
    
    Measures the distribution of the higher gray-level values, with a higher value indicating  
    a greater concentration of high gray-level values in the image.
    """    
    pg = self.coefficients['pg']
    ivector = self.coefficients['ivector']
    sumP_rlgl = self.coefficients['sumP_rlgl']

    try:
      hglre = numpy.sum( (pg*(ivector[:,None]**2)) , 0 ) / sumP_rlgl
    except ZeroDivisionError:
      hglre = 0
    return (hglre.mean())

  def getShortRunLowGrayLevelEmphasisFeatureValue(self):
    """Calculate and return the mean Short Run Low Gray Level Emphasis (SRLGLE) value for all 13 RLGL matrices.
    
    Measures the joint distribution of shorter run lengths with lower gray-level values.
    """    
    ivector = self.coefficients['ivector']
    jvector = self.coefficients['jvector']
    sumP_rlgl = self.coefficients['sumP_rlgl']

    try:
      srlgle = numpy.sum( (self.P_rlgl/((ivector[:,None,None]**2)*(jvector[None,:,None]**2))) , (0, 1) ) / sumP_rlgl
    except ZeroDivisionError:
      srlgle = 0
    return (srlgle.mean())

  def getShortRunHighGrayLevelEmphasisFeatureValue(self):
    """Calculate and return the mean Short Run High Gray Level Emphasis (SRHGLE) value for all 13 RLGL matrices.
    
    Measures the joint distribution of shorter run lengths with higher gray-level values.
    """    
    ivector = self.coefficients['ivector']
    jvector = self.coefficients['jvector']
    sumP_rlgl = self.coefficients['sumP_rlgl']

    try:
      srhgle = numpy.sum( (self.P_rlgl*(ivector[:,None,None]**2)/(jvector[None,:,None]**2)) , (0, 1) ) / sumP_rlgl
    except ZeroDivisionError:
      srhgle = 0
    return (srhgle.mean())

  def getLongRunLowGrayLevelEmphasisFeatureValue(self):
    """Calculate and return the mean Long Run Low Gray Level Emphasis (LRLGLE) value for all 13 RLGL matrices.
    
    Measures the joint distribution of long run lengths with lower gray-level values.
    """    
    ivector = self.coefficients['ivector']
    jvector = self.coefficients['jvector']
    sumP_rlgl = self.coefficients['sumP_rlgl']

    try:
      lrlgle = numpy.sum( (self.P_rlgl*(jvector[None,:,None]**2)/(ivector[:,None,None]**2)) , (0, 1) ) / sumP_rlgl
    except ZeroDivisionError:
      lrlgle = 0
    return (lrlgle.mean())

  def getLongRunHighGrayLevelEmphasisFeatureValue(self):
    """Calculate and return the mean Long Run High Gray Level Emphasis (LRHGLE) value for all 13 RLGL matrices.
    
    Measures the joint distribution of long run lengths with higher gray-level values.
    """    
    ivector = self.coefficients['ivector']
    jvector = self.coefficients['jvector']
    sumP_rlgl = self.coefficients['sumP_rlgl']

    try:
      lrhgle = numpy.sum( (self.P_rlgl*((jvector[None,:,None]**2)*(ivector[:,None,None]**2))) , (0, 1) ) / sumP_rlgl
    except ZeroDivisionError:
      lrhgle = 0
    return (lrhgle.mean())
