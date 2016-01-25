import numpy
import SimpleITK as sitk
from scipy import ndimage
from radiomics import base, preprocessing


class RadiomicsGLSZM(base.RadiomicsFeaturesBase):
  """GLSZM feature calculation."""
  def __init__(self, inputImage, inputMask, **kwargs):
      super(RadiomicsGLSZM,self).__init__(inputImage, inputMask, **kwargs)

      self.imageArray = sitk.GetArrayFromImage(inputImage)
      self.maskArray = sitk.GetArrayFromImage(inputMask)

      (self.matrix, self.matrixCoordinates) = \
        preprocessing.RadiomicsHelpers.padTumorMaskToCube(self.imageArray,self.maskArray)

      self.targetVoxelArray = self.matrix[self.matrixCoordinates]
      self.coefficients = {}
      self.P_glszm = {}

      # binning
      self.matrix, self.histogram = preprocessing.RadiomicsHelpers.binImage(self.binWidth, self.targetVoxelArray, self.matrix, self.matrixCoordinates)
      self.coefficients['Ng'] = len(self.histogram[0])
      self.coefficients['grayLevels'] = numpy.linspace(1,self.coefficients['Ng'],num=self.coefficients['Ng'])
      self.coefficients['Nr'] = numpy.max(self.matrix.shape)
      self.coefficients['Np'] = self.targetVoxelArray.size

      self.calculateGLSZM(angles=13)

      self.calculateCoefficients()

      
  def calculateGLSZM(self, angles=13):
    """
    Number of times a 26-connected region with a 
    gray level and voxel count occurs in an image. P_glszm[level, voxel_count] = # occurrences
    """
    
    # Kernel for 26-connected neighborhood
    B = numpy.ones((3,3,3))
    
    # Maximum size of a zone is the total volume of the image
    maxSize = self.matrix[self.matrixCoordinates].size
    
    # Empty GLSZ matrix
    self.P_glszm = numpy.zeros((self.coefficients['grayLevels'], maxSize))
    
    # Iterate over all gray levels in the image
    for i in xrange(1, self.coefficients['grayLevels']+1):
      #dataTemp = self.matrix[numpy.where(self.matrix==i)] #should this be 3D?
      dataTemp = numpy.where(self.matrix==i, 1, 0)
      ind = zip(*numpy.where(self.matrix==i))  
      labels = numpy.zeros(dataTemp.shape)
      n = 0
      while len(ind) > 0: # check if ind is not empty
        # Current label number and first coordinate
        n = n+1
        ind = ind[0]
        X = numpy.zeros(dataTemp.shape)
        X[ind] = 1
        
        # Iterative dilating with kernel       
        Y = dataTemp and (ndimage.convolve(X,B,mode='constant') >= 1) # what is this output?
        while X != Y:
            X = Y;
            Y = dataTemp and (ndimage.convolve(X,B,mode='constant') >= 1)
        # Y appears to be a boolean array. the final "region".
        # intersection of convolution and dataTemp until it equals X?
        
        # Set already processed indices to zero
        dataTemp[Y] = 0 #is Y a set of coordinates now? 
        # or [dataTemp[pos] = 0 for pos in Y]
        
        # Assign the label n
        labels[Y] = n #is Y a label value?
        # or [labels[pos] = n for pos in Y]
        
        # Size of the region (# of voxels in region)
        regionSize = len(zip(*numpy.where(Y!=0)))
        
        # Update the gray level size zone matrix
        self.P_glszm[i,regionSize] = self.P_glszm[i,regionSize] += 1
        
        # Find unprocessed nonzero positions for current gray level
        ind = numpy.where(dataTemp==i)
              
  def calculateCoefficients(self):
      sumP_glszm = numpy.sum( numpy.sum(self.P_rlgl, 0), 0 )
      sumP_glszm[sumP_glszm==0] = 1

      pr = numpy.sum(self.P_glszm, 0)
      pg = numpy.sum(self.P_glszm, 1)

      ivector = numpy.arange(1, self.P_glszm.shape[0] + 1)
      jvector = numpy.arange(1, self.P_glszm.shape[1] + 1)

      self.coefficients['sumP_glszm'] = sumP_glszm
      self.coefficients['pr'] = pr
      self.coefficients['pg'] = pg
      self.coefficients['ivector'] = ivector
      self.coefficients['jvector'] = jvector     
    
  def smallAreaEmphasis(self):
      try:
        #sre = numpy.sum( numpy.sum( (self.P_glszm/(self.j[:,:,None]**2)) , 0 ), 0 ) / (self.coefficients['sumP_glszm'][None,None,:])
        sre = numpy.sum( (self.coefficients['pr']/(self.coefficients['jvector'][:,None]**2)), 0 ) / self.coefficients['sumP_glszm']
      except ZeroDivisionError:
        sre = 0    
      return (sre.mean())
    
  def largeAreaEmphasis(self):
      try:
        #lre = numpy.sum( numpy.sum( (self.P_glszm*(self.j[:,:,None]**2)) , 0 ), 0 ) / (self.coefficients['sumP_glszm'][None,None,:])
        lre = numpy.sum( (self.coefficients['pr']*(self.coefficients['jvector'][:,None]**2)), 0 ) / self.coefficients['sumP_glszm']
      except ZeroDivisionError:
        lre = 0
      return (lre.mean())

  def intensityVariability(self):
      try:
        gln = numpy.sum( (self.coefficients['pg']**2) , 0 ) / self.coefficients['sumP_glszm']   
      except ZeroDivisionError:
        gln = 0  
      return (gln.mean())
    
  def sizeZoneVariability(self):
      try:
        rln = numpy.sum( (self.coefficients['pr']**2) , 0 ) / self.coefficients['sumP_glszm'] 
      except ZeroDivisionError:
        rln = 0
      return (rln.mean())

  def zonePercentage(self):
      try:
        rp = numpy.sum( numpy.sum( (self.P_glszm/(self.coefficients['Np'])) , 0 ), 0 )
      except ZeroDivisionError:
        rp = 0
      return (rp.mean())

  def lowIntensityEmphasis(self):
      try:
        #pdb.set_trace()
        #lglre = numpy.sum( numpy.sum( (self.P_glszm/(self.i[:,:,None]**2)) , 0 ), 0 )/ (self.coefficients['sumP_glszm'][None,None,:])
        lglre = numpy.sum( (self.coefficients['pg']/(self.coefficients['ivector'][:,None]**2)) , 0 ) / self.coefficients['sumP_glszm']
      except ZeroDivisionError:
        lglre = 0 
      return (lglre.mean())
      
  def highIntensityEmphasis(self):
    try:
      #hglre = numpy.sum( numpy.sum( (self.P_glszm*(self.i[:,:,None]**2)) , 0 ), 0 ) / (self.coefficients['sumP_glszm'][None,None,:])
      hglre = numpy.sum( (self.coefficients['pg']*(self.coefficients['ivector'][:,None]**2)) , 0 ) / self.coefficients['sumP_glszm']
    except ZeroDivisionError:
      hglre = 0 
    return (hglre.mean())
    
  def lowIntensitySmallAreaEmphasis(self):
    try:
      srlgle = numpy.sum( numpy.sum( (self.P_glszm/((self.coefficients['ivector'][:,None,None]**2)*(self.coefficients['jvector'][None,:,None]**2))) , 0 ), 0 ) / self.coefficients['sumP_glszm']
    except ZeroDivisionError:
      srlgle = 0
    return (srlgle.mean())
     
  def highIntensitySmallAreaEmphasis(self):
    try:
      srhgle = numpy.sum( numpy.sum( (self.P_glszm*(self.coefficients['ivector'][:,None,None]**2)/(self.coefficients['jvector'][None,:,None]**2)) , 0 ), 0 ) / self.coefficients['sumP_glszm']
    except ZeroDivisionError:
      srhgle = 0   
    return (srhgle.mean())
        
  def lowIntensityLargeAreaEmphasis(self):
    try:
      lrlgle = numpy.sum( numpy.sum( (self.P_glszm*(self.coefficients['jvector'][None,:,None]**2)/(self.coefficients['ivector'][:,None,None]**2)) , 0 ), 0 ) / self.coefficients['sumP_glszm']
    except ZeroDivisionError:
      lrlgle = 0
    return (lrlgle.mean())
            
  def highIntensityLargeAreaEmphasis(self):
    try:
      lrhgle = numpy.sum( numpy.sum( (self.P_glszm*((self.coefficients['jvector'][None,:,None]**2)*(self.coefficients['ivector'][:,None,None]**2))) , 0 ), 0 ) / self.coefficients['sumP_glszm']
    except ZeroDivisionError:
      lrhgle = 0 
    return (lrhgle.mean())    