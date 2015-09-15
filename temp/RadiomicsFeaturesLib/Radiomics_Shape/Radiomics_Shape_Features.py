import string
import numpy
import math
import operator
import collections
    
def voxelNumber(parameterArray):
    return (parameterArray.size)
  
def volumeMM3 (matrixSA, cubicMMPerVoxel):      
    return (matrixSA.size * cubicMMPerVoxel)
  
"""
def surfaceAreaNew(a, matrixSACoordinates, matrixSAValues, pixelSpacing):
    x, y, z = pixelSpacing 
    p, q, r = numpy.meshgrid(y*numpy.arange(a.shape[2]), x*numpy.arange(a.shape[1]), z*numpy.arange(a.shape[0]))
    faces, vertices = calculateIsosurface(p,q,r,a,0.5)
    
    a = vertices[faces[:, 2], :] - vertices[faces[:, 1], :];
    b = vertices[faces[:, 3], :] - vertices[faces[:, 1], :];
    c = cross(a, b, 2);

    return( 0.5 * numpy.sum(numpy.sqrt(numpy.sum(numpy.float(c)**2, axis=2))))

def calculateIsosurface(x,y,z,v,f):
    pass
"""  

def surfaceArea(a, matrixSACoordinates, matrixSAValues, pixelSpacing):
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
        fxy = (numpy.array([ a[i+1,j,k], a[i-1,j,k] ]) == 0) # evaluate to 1 if true, 0 if false
        fyz = (numpy.array([ a[i,j+1,k], a[i,j-1,k] ]) == 0) # evaluate to 1 if true, 0 if false
        fxz = (numpy.array([ a[i,j,k+1], a[i,j,k-1] ]) == 0) # evaluate to 1 if true, 0 if false  
        surface = (numpy.sum(fxz) * xz) + (numpy.sum(fyz) * yz) + (numpy.sum(fxy) * xy)     
        surfaceArea += surface
        
    return (surfaceArea)
     
def surfaceVolumeRatio (surfaceArea, volumeMM3):      
    return (surfaceArea/volumeMM3)
         
def compactness1 (surfaceArea, volumeMM3):      
    return ( (volumeMM3) / ((surfaceArea)**(2/3.0) * math.sqrt(math.pi)) )
   
def compactness2 (surfaceArea, volumeMM3):      
    return ((36 * math.pi) * ((volumeMM3)**2)/((surfaceArea)**3)) 

def maximum3DDiameter(matrixSA, matrixSACoordinates, pixelSpacing):
    # largest pairwise euclidean distance between tumor surface voxels
     
    x, y, z = pixelSpacing
    
    minBounds = numpy.array([numpy.min(matrixSACoordinates[0]), numpy.min(matrixSACoordinates[1]), numpy.min(matrixSACoordinates[2])])
    maxBounds = numpy.array([numpy.max(matrixSACoordinates[0]), numpy.max(matrixSACoordinates[1]), numpy.max(matrixSACoordinates[2])])
    
    a = numpy.array(zip(*matrixSACoordinates))
    edgeVoxelsMinCoords = numpy.vstack([a[a[:,0]==minBounds[0]], a[a[:,1]==minBounds[1]], a[a[:,2]==minBounds[2]]]) * [z,y,x]
    edgeVoxelsMaxCoords = numpy.vstack([(a[a[:,0]==maxBounds[0]]+1), (a[a[:,1]==maxBounds[1]]+1), (a[a[:,2]==maxBounds[2]]+1)]) * [z,y,x]
    
    maxDiameter = 1
    for voxel1 in edgeVoxelsMaxCoords:
        for voxel2 in edgeVoxelsMinCoords:       
            voxelDistance = numpy.sqrt(numpy.sum((voxel2-voxel1)**2))
            if voxelDistance > maxDiameter:
                maxDiameter= voxelDistance
                
    return(maxDiameter)     
      
def sphericalDisproportion (surfaceArea, volumeMM3):
    R = ( (3*volumeMM3)/(4*math.pi) )**(1/3.0)   
    return ( (surfaceArea)/(4*math.pi*(R**2)) )
      
def sphericityValue(surfaceArea, volumeMM3):      
    return ( ((math.pi)**(1/3.0) * (6 * volumeMM3)**(2/3.0)) / (surfaceArea) ) 