import numpy
import SimpleITK as sitk
from radiomics import base, preprocessing
import pdb

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

      self.calculateGLSZM()

      self.calculateCoefficients()

      
  def calculateGLSZM(self):
    """
    Number of times a 26-connected region with a 
    gray level and voxel count occurs in an image. P_glszm[level, voxel_count] = # occurrences
    """
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
    
    # Kernel for 26-connected neighborhood
    B = numpy.ones((3,3,3))
    
    # Empty GLSZ matrix
    self.P_glszm = numpy.zeros((self.coefficients['grayLevels'].size, self.coefficients['Np']))
    
    # Iterate over all gray levels in the image
    for i in xrange(1, self.coefficients['grayLevels'].size+1):
      dataTemp = numpy.where(self.matrix==i, 1, 0)
      ind = zip(*numpy.where(dataTemp==1))
      ind = list(set(ind).intersection(set(zip(*self.matrixCoordinates))))      
      labels = numpy.zeros(dataTemp.shape)
      n = 0
      while ind: # check if ind is not empty
        # Current label number and first coordinate
        n = n+1
        ind_node = ind[0]
        
        # get all coordinates in the 26-connected region
        region_full = [ind_node] + [tuple(sum(a) for a in zip(ind_node,angle_i)) for angle_i in angles]
        
        # get coordinates in 26-connected region with same grey level
        region_level = list(set(ind).intersection(set(region_full)))
        
        # Set already processed indices to zero
        for pos in region_level:
          dataTemp[pos] = 0
        
        # Assign the label n
        for pos in region_level:
          labels[pos] = n
        
        # Size of the region (# of voxels in region)
        #regionSize = len(zip(*numpy.where(Y!=0)))
        regionSize = len(region_level)
        
        # Update the gray level size zone matrix
        self.P_glszm[i-1,regionSize-1] += 1
        
        # Find unprocessed nonzero positions for current gray level
        ind = zip(*numpy.where(dataTemp==1))
        ind = list(set(ind).intersection(set(zip(*self.matrixCoordinates))))
              
  def calculateCoefficients(self):
      sumP_glszm = numpy.sum( numpy.sum(self.P_glszm, 0), 0 )
      
      # set sum to numpy.spacing(1) if sum is 0?
      if sumP_glszm == 0:
        sumP_glszm = 1

      pr = numpy.sum(self.P_glszm, 0)
      pg = numpy.sum(self.P_glszm, 1)

      ivector = numpy.arange(1, self.P_glszm.shape[0] + 1)
      jvector = numpy.arange(1, self.P_glszm.shape[1] + 1)

      self.coefficients['sumP_glszm'] = sumP_glszm
      self.coefficients['pr'] = pr
      self.coefficients['pg'] = pg
      self.coefficients['ivector'] = ivector
      self.coefficients['jvector'] = jvector     
    
  def getSmallAreaEmphasisFeatureValue(self):
      try:
        #sum(pr./(j_vector.^2))/nruns
        #sae = numpy.sum( numpy.sum( (self.P_glszm/(self.j[:,:,None]**2)) , 0 ), 0 ) / (self.coefficients['sumP_glszm'][None,None,:])
        sae = numpy.sum(self.coefficients['pr']/(self.coefficients['jvector']**2)) / self.coefficients['sumP_glszm']
      except ZeroDivisionError:
        sae = 0    
      return (sae)
    
  def getLargeAreaEmphasisFeatureValue(self):
      try:
        #lae = numpy.sum( numpy.sum( (self.P_glszm*(self.j[:,:,None]**2)) , 0 ), 0 ) / (self.coefficients['sumP_glszm'][None,None,:])
        lae = numpy.sum(self.coefficients['pr']*(self.coefficients['jvector']**2)) / self.coefficients['sumP_glszm']
      except ZeroDivisionError:
        lae = 0
      return (lae)

  def getIntensityVariabilityFeatureValue(self):
      try:
        iv = numpy.sum(self.coefficients['pg']**2) / self.coefficients['sumP_glszm']   
      except ZeroDivisionError:
        iv = 0  
      return (iv)
    
  def getSizeZoneVariabilityFeatureValue(self):
      try:
        szv = numpy.sum(self.coefficients['pr']**2) / self.coefficients['sumP_glszm'] 
      except ZeroDivisionError:
        szv = 0
      return (szv)

  def getZonePercentageFeatureValue(self):
      try:
        #zp = numpy.sum( (self.P_glszm/(self.coefficients['Np'])) , 0 )
        zp = self.coefficients['sumP_glszm'] / numpy.sum(self.coefficients['pr']*self.coefficients['jvector'])
      except ZeroDivisionError:
        zp = 0
      return (zp)

  def getLowIntensityEmphasisFeatureValue(self):
      try:
        #pdb.set_trace()
        #lie = numpy.sum( numpy.sum( (self.P_glszm/(self.i[:,:,None]**2)) , 0 ), 0 )/ (self.coefficients['sumP_glszm'][None,None,:])
        lie = numpy.sum( (self.coefficients['pg']/(self.coefficients['ivector']**2)) ) / self.coefficients['sumP_glszm']
      except ZeroDivisionError:
        lie = 0 
      return (lie)
      
  def getHighIntensityEmphasisFeatureValue(self):
    try:
      #hie = numpy.sum( numpy.sum( (self.P_glszm*(self.i[:,:,None]**2)) , 0 ), 0 ) / (self.coefficients['sumP_glszm'][None,None,:])
      hie = numpy.sum( (self.coefficients['pg']*(self.coefficients['ivector']**2)) ) / self.coefficients['sumP_glszm']
    except ZeroDivisionError:
      hie = 0 
    return (hie)
    
  def getLowIntensitySmallAreaEmphasisFeatureValue(self):
    try:
      lisae = numpy.sum( numpy.sum( (self.P_glszm/((self.coefficients['ivector'][:,None]**2)*(self.coefficients['jvector'][None,:]**2))) , 0 ), 0 ) / self.coefficients['sumP_glszm']
    except ZeroDivisionError:
      lisae = 0
    return (lisae)
     
  def getHighIntensitySmallAreaEmphasisFeatureValue(self):
    try:
      hisae = numpy.sum( numpy.sum( (self.P_glszm*(self.coefficients['ivector'][:,None]**2)/(self.coefficients['jvector'][None,:]**2)) , 0 ), 0 ) / self.coefficients['sumP_glszm']
    except ZeroDivisionError:
      hisae = 0   
    return (hisae)
        
  def getLowIntensityLargeAreaEmphasisFeatureValue(self):
    try:
      lilae = numpy.sum( numpy.sum( (self.P_glszm*(self.coefficients['jvector'][None,:]**2)/(self.coefficients['ivector'][:,None]**2)) , 0 ), 0 ) / self.coefficients['sumP_glszm']
    except ZeroDivisionError:
      lilae = 0
    return (lilae)
            
  def getHighIntensityLargeAreaEmphasisFeatureValue(self):
    try:
      hilae = numpy.sum( numpy.sum( (self.P_glszm*((self.coefficients['jvector'][None,:]**2)*(self.coefficients['ivector'][:,None]**2))) , 0 ), 0 ) / self.coefficients['sumP_glszm']
    except ZeroDivisionError:
      hilae = 0 
    return (hilae)    