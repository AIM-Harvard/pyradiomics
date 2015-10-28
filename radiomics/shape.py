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

    (self.matrix, self.matrixCoordinates) = \
      preprocessing.RadiomicsHelpers.padTumorMaskToCube(self.imageArray,self.maskArray)

    self.targetVoxelArray = self.matrix[self.matrixCoordinates]

    # Pad and center the inputMask matrix by 10 voxels in three directions.
    self.padding = ([10,10,10])
    self.maxDimsSA = tuple(map(operator.add, self.matrix.shape, self.padding))
    self.matrixSA, self.matrixSACoordinates = \
      preprocessing.RadiomicsHelpers.padCubicMatrix(self.matrix, self.matrixCoordinates, self.maxDimsSA)

    # Volume and Surface Area are pre-calculated
    self.Volume = self.volumeMM3(self.targetVoxelArray, self.cubicMMPerVoxel)
    self.SurfaceArea = self.surfaceArea(self.matrixSA, self.matrixSACoordinates, self.targetVoxelArray, self.pixelSpacing)

    #self.InitializeFeatureVector()
    #for f in self.getFeatureNames():
    #  self.enabledFeatures[f] = True

    # TODO: add an option to instantiate the class that reuses initialization
       
  def getVoxelNumberFeatureValue(targetVoxelArray):
    return (targetVoxelArray.size)
    
  def getVolumeMM3FeatureValue(matrixSA, cubicMMPerVoxel):      
    return (matrixSA.size * cubicMMPerVoxel)
    
  """
  def getSurfaceAreaFeatureValueNew(a, matrixSACoordinates, matrixSAValues, pixelSpacing):
    x, y, z = pixelSpacing 
    p, q, r = numpy.meshgrid(y*numpy.arange(a.shape[2]), x*numpy.arange(a.shape[1]), z*numpy.arange(a.shape[0]))
    faces, vertices = calculateIsosurface(p,q,r,a,0.5)
    
    a = vertices[faces[:, 2], :] - vertices[faces[:, 1], :];
    b = vertices[faces[:, 3], :] - vertices[faces[:, 1], :];
    c = cross(a, b, 2);

    return( 0.5 * numpy.sum(numpy.sqrt(numpy.sum(numpy.float(c)**2, axis=2))))

  def _calculateIsosurface(x,y,z,v,f):
    pass
  """  

  def getSurfaceAreaFeatureValue(matrixSA, matrixSACoordinates, matrixSAValues, pixelSpacing):
    x, y, z = pixelSpacing
    xz = x*z
    yz = y*z
    xy = x*y
    voxelTotalSA = (2*xy + 2*xz + 2*yz)
    totalSA = matrixSAValues.size * voxelTotalSA
    
    # in matrixSACoordinates:
    # i corresponds to height (z)
    # j corresponds to vertical (y)
    # k corresponds to horizontal (x)
    
    i, j, k = 0, 0, 0
    surfaceArea = 0   
    for voxel in xrange(0, matrixSAValues.size):
      i, j, k = matrixSACoordinates[0][voxel], matrixSACoordinates[1][voxel], matrixSACoordinates[2][voxel]      
      fxy = (numpy.array([ matrixSA[i+1,j,k], matrixSA[i-1,j,k] ]) == 0) # evaluate to 1 if true, 0 if false
      fyz = (numpy.array([ matrixSA[i,j+1,k], matrixSA[i,j-1,k] ]) == 0) # evaluate to 1 if true, 0 if false
      fxz = (numpy.array([ matrixSA[i,j,k+1], matrixSA[i,j,k-1] ]) == 0) # evaluate to 1 if true, 0 if false  
      surface = (numpy.sum(fxz) * xz) + (numpy.sum(fyz) * yz) + (numpy.sum(fxy) * xy)     
      surfaceArea += surface
      
    return (surfaceArea)
     
  def getSurfaceVolumeRatioFeatureValue(surfaceArea, volumeMM3):      
    return (surfaceArea/volumeMM3)
           
  def getCompactness1FeatureValue(surfaceArea, volumeMM3):      
    return ( (volumeMM3) / ((surfaceArea)**(2.0/3.0) * math.sqrt(math.pi)) )
     
  def getCompactness2FeatureValue(surfaceArea, volumeMM3):      
    return ((36.0 * math.pi) * ((volumeMM3)**2.0)/((surfaceArea)**3.0)) 

  def getMaximum3DDiameterFeatureValue(matrixSA, matrixSACoordinates, pixelSpacing):
    # largest pairwise euclidean distance between tumor surface voxels   
    x, y, z = pixelSpacing
    
    minBounds = numpy.array([numpy.min(matrixSACoordinates[0]), numpy.min(matrixSACoordinates[1]), numpy.min(matrixSACoordinates[2])])
    maxBounds = numpy.array([numpy.max(matrixSACoordinates[0]), numpy.max(matrixSACoordinates[1]), numpy.max(matrixSACoordinates[2])])
    
    a = numpy.array(zip(*matrixSACoordinates))
    edgeVoxelsMinCoords = numpy.vstack([a[a[:,0]==minBounds[0]], a[a[:,1]==minBounds[1]], a[a[:,2]==minBounds[2]]]) * [z,y,x]
    edgeVoxelsMaxCoords = numpy.vstack([(a[a[:,0]==maxBounds[0]]+1), (a[a[:,1]==maxBounds[1]]+1), (a[a[:,2]==maxBounds[2]]+1)]) * [z,y,x]
    
    maxDiameter = 1   
    for voxel1 in edgeVoxelsMaxCoords:
      voxelDistances = numpy.sqrt(numpy.sum((edgeVoxelsMinCoords-voxel1)**2))
      if voxelDistance.max() > maxDiameter: 
        maxDiameter = voxelDistance.max()
      
    return(maxDiameter)     
      
  def getSphericalDisproportionFeatureValue(surfaceArea, volumeMM3):
    R = ( (3.0*volumeMM3)/(4.0*math.pi) )**(1.0/3.0)   
    return ( (surfaceArea)/(4.0*math.pi*(R**2.0)) )
      
  def getSphericityFeatureValue(surfaceArea, volumeMM3):      
    return ( ((math.pi)**(1.0/3.0) * (6.0 * volumeMM3)**(2.0/3.0)) / (surfaceArea) ) 