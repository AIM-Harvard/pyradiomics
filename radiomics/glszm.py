import numpy
import SimpleITK as sitk
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

      self.calculateGLSZM(13)

      self.calculateCoefficients()

  def calculateGLSZM(self, angles):
  def calculate_glszm(grayLevels, matrix, matrixCoordinates, angles, P_rlgl):
    # Kernel for 26-connected neighborhood
    B = numpy.ones((3,3,3))
    maxSize = matrix[matrixCoordinates].size
    
    P_glszm = numpy.zeros((grayLevels, maxSize))
    
    for i in xrange(grayLevels):
      level = i + 1
      grayLevelCoordinates = zip(*numpy.where(matrix==level))
      dataTemp = matrix[numpy.where(matrix==level)]
      labels = numpy.zeros(len(grayLevelCoordinates))
      n = 0
      while len(grayLevelCoordinates) > 0:
        n = n+1
        ind = grayLevelCoordinates[0]
        X = numpy.zeros(dataTemp.size)
        """
        #Matlab
        X(ind) = 1;
        %   Iterative dilating with kernel
        Y = dataTemp&(convn(X,B,'same')>=1);
        while ~isequal(X,Y)
            X = Y;
            Y = dataTemp&(convn(X,B,'same')>=1);
        end
        %   Set already processed indices to zero
        dataTemp(Y) = 0;        
        %   Assign the label n
        labels(Y) = n;
        %   Size of the region
        regionSize = length(find(Y));
        %   Update the gray level size zone matrix
        glszMat(i,regionSize) = glszMat(i,regionSize)+1;
        %   Find unprocessed nonzero positions for current gray level
        ind = find(dataTemp==i);
        end
        %   Now we have a labeled image with n 26-neighborhood connected
        %   components, having a gray level of i.
        %   For now "labels" is not used, but maybe something for the future.
        end
      """      
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