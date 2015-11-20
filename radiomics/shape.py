import numpy
import operator
import collections
from radiomics import base, preprocessing
import SimpleITK as sitk

class RadiomicsShape(base.RadiomicsFeaturesBase):

  def __init__(self, inputImage, inputMask):
    super(RadiomicsFirstOrder,self).__init__(inputImage,inputMask)

    self.pixelSpacing = inputImage.GetSpacing()
    self.cubicMMPerVoxel = reduce(lambda x,y: x*y , self.pixelSpacing)
    
    #self.featureNames = self.getFeatureNames()

    self.imageArray = sitk.GetArrayFromImage(inputImage)
    self.maskArray = sitk.GetArrayFromImage(inputMask)
        
    # Generate a cuboid matrix of the tumor region and pad by 10 voxels in three directions 
    # for surface area calculation
    (self.matrix, self.matrixCoordinates) = \
      preprocessing.RadiomicsHelpers.padTumorMaskToCube(self.imageArray, self.maskArray, padDistance=10)
    self.targetVoxelArray = self.matrix[self.matrixCoordinates]

    # Volume and Surface Area are pre-calculated
    self.Volume = self.getVolumeMM3FeatureValue(self.targetVoxelArray, self.cubicMMPerVoxel)
    self.SurfaceArea = self.getSurfaceAreaFeatureValue(self.matrix, self.matrixCoordinates, self.targetVoxelArray, self.pixelSpacing)

    #self.InitializeFeatureVector()
    #for f in self.getFeatureNames():
    #  self.enabledFeatures[f] = True

    # TODO: add an option to instantiate the class that reuses initialization
       
  def getVolumeMM3FeatureValue(self, matrix, cubicMMPerVoxel):
    """Calculate the volume of the tumor region in cubic millimeters."""
    return (matrix.size * cubicMMPerVoxel)

  def getSurfaceAreaFeatureValue(self, matrix, matrixCoordinates, matrixValues, pixelSpacing):
    """Calculate the surface area of the tumor region in square millimeters."""
    x, y, z = pixelSpacing
    xz = x*z
    yz = y*z
    xy = x*y
    voxelTotalSA = (2*xy + 2*xz + 2*yz)
    totalSA = matrixValues.size * voxelTotalSA
    
    # in matrixCoordinates:
    # i corresponds to height (z)
    # j corresponds to vertical (y)
    # k corresponds to horizontal (x)
    
    i, j, k = 0, 0, 0
    surfaceArea = 0   
    for voxel in xrange(0, matrixValues.size):
      i, j, k = matrixCoordinates[0][voxel], matrixCoordinates[1][voxel], matrixCoordinates[2][voxel]      
      fxy = (numpy.array([ matrix[i+1,j,k], matrix[i-1,j,k] ]) == 0) # evaluate to 1 if true, 0 if false
      fyz = (numpy.array([ matrix[i,j+1,k], matrix[i,j-1,k] ]) == 0) # evaluate to 1 if true, 0 if false
      fxz = (numpy.array([ matrix[i,j,k+1], matrix[i,j,k-1] ]) == 0) # evaluate to 1 if true, 0 if false  
      surface = (numpy.sum(fxz) * xz) + (numpy.sum(fyz) * yz) + (numpy.sum(fxy) * xy)     
      surfaceArea += surface
      
    return (surfaceArea)
     
  def getSurfaceVolumeRatioFeatureValue(self):
    """Calculate the surface area to volume ratio of the tumor region"""
    return (self.SurfaceArea/self.Volume)
           
  def getCompactness1FeatureValue(self):
    """
    Calculate the compactness (1) of the tumor region.
    
    Compactness 1 is a measure of how compact the shape of the tumor is relative to a sphere (most compact).
    """
    return ( (self.Volume) / ((self.SurfaceArea)**(2.0/3.0) * math.sqrt(math.pi)) )
     
  def getCompactness2FeatureValue(self):
    """
    Calculate the Compactness (2) of the tumor region.
    
    Compactness 2 is a measure of how compact the shape of the tumor is relative to a sphere (most compact).
    """  
    return ((36.0 * math.pi) * ((self.Volume)**2.0)/((self.SurfaceArea)**3.0)) 

  def getMaximum3DDiameterFeatureValue(self, matrix, matrixCoordinates, pixelSpacing):
    """Calculate the largest pairwise euclidean distance between tumor surface voxels"""   
    x, y, z = pixelSpacing
    
    minBounds = numpy.array([numpy.min(matrixCoordinates[0]), numpy.min(matrixCoordinates[1]), numpy.min(matrixCoordinates[2])])
    maxBounds = numpy.array([numpy.max(matrixCoordinates[0]), numpy.max(matrixCoordinates[1]), numpy.max(matrixCoordinates[2])])
    
    a = numpy.array(zip(*matrixCoordinates))
    edgeVoxelsMinCoords = numpy.vstack([a[a[:,0]==minBounds[0]], a[a[:,1]==minBounds[1]], a[a[:,2]==minBounds[2]]]) * [z,y,x]
    edgeVoxelsMaxCoords = numpy.vstack([(a[a[:,0]==maxBounds[0]]+1), (a[a[:,1]==maxBounds[1]]+1), (a[a[:,2]==maxBounds[2]]+1)]) * [z,y,x]
    
    maxDiameter = 1   
    for voxel1 in edgeVoxelsMaxCoords:
      voxelDistances = numpy.sqrt(numpy.sum((edgeVoxelsMinCoords-voxel1)**2))
      if voxelDistance.max() > maxDiameter: 
        maxDiameter = voxelDistance.max()
      
    return(maxDiameter)     
      
  def getSphericalDisproportionFeatureValue(self):
    """
    Calculate the Spherical Disproportion of the tumor region.
    
    Spherical Disproportion is the ratio of the surface area of the
    tumor region to the surface area of a sphere with the same 
    volume as the tumor region.
    """ 
    R = ( (3.0*self.Volume)/(4.0*math.pi) )**(1.0/3.0)   
    return ( (self.SurfaceArea)/(4.0*math.pi*(R**2.0)) )
      
  def getSphericityFeatureValue(self):
    """
    Calculate the Sphericity of the tumor region.
    
    Sphericity is a measure of the roundness of the shape of the tumor region
    relative to a sphere. This is another measure of the compactness of a tumor.
    """   
    return ( ((math.pi)**(1.0/3.0) * (6.0 * self.Volume)**(2.0/3.0)) / (self.SurfaceArea) ) 