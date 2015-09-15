import numpy
import collections
import pdb

def CalculateGLCM(matrix, matrixCoordinates, distances, angles, P_glcm):
    # 13 GLCM matrices for each image for every direction from the voxel 
    # (13 for each neighboring voxel from a reference voxel centered in a 3x3 cube)
    # for GLCM matrices P(i,j;gamma, a), gamma = 1, a = 1...13
    
    maxrows = matrix.shape[2]
    maxcols = matrix.shape[1]
    maxheight = matrix.shape[0]
    indices = zip(*matrixCoordinates)
            
    for h, c, r in indices:
        for angles_idx, angle in enumerate(angles):
            for distances_idx, distance in enumerate(distances):
                i = matrix[h, c, r]
                i_idx = int(i-1)

                row = r + angle[2]
                col = c + angle[1]
                height = h + angle[0]
                
                if row >= 0 and row < maxrows and col >= 0 and col < maxcols:
                    if tuple((height, col, row)) in indices:
                        j = matrix[height, col, row]
                        j_idx = int(j-1)
                        P_glcm[i_idx, j_idx, distances_idx, angles_idx] += 1          
              
    return (P_glcm)
  
def CreateGLCM(Ng, matrix, matrixCoordinates):
    # Generate container for GLCM Matrices: P_glcm
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
    P_glcm = numpy.zeros( (Ng, Ng, distances.size, int(angles.shape[0])), dtype='float32' )
    P_glcm = CalculateGLCM(matrix, matrixCoordinates, distances, angles, P_glcm) 
    
    #Normalize each glcm
    for q in xrange(int(angles.shape[0])):
        P_glcm[:,:,0,q] = P_glcm[:,:,0,q]/(P_glcm[:,:,0,q].sum())
      
    return P_glcm  

