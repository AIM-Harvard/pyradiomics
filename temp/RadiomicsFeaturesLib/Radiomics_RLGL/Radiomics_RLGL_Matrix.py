import numpy
import collections
import pdb
  
def calculate_rlgl(grayLevels, matrix, matrixCoordinates, angles, P_rlgl):    
    padVal = -1000   #use eps or NaN to pad matrix
    padMask = numpy.zeros(matrix.shape,dtype=bool)
    padMask[matrixCoordinates] = True
    matrix[~padMask] = padVal
    matrixDiagonals = list()
    
    # TODO: try using itertools list merging with lists of rlgl diagonal    
    # i.e.: self.heterogeneityFeatureWidgets = list(itertools.chain.from_iterable(self.featureWidgets.values()))
    
    # For a single direction or diagonal (aDiags, bDiags...lDiags, mDiags):   
    # Generate a 1D array for each valid offset of the diagonal, a, in the range specified by lowBound and highBound  
    # Convert each 1D array to a python list ( matrix.diagonal(a,,).tolist() ) 
    # Join lists using reduce(lamda x,y: x+y, ...) to represent all 1D arrays for the direction/diagonal       
    # Use filter(lambda x: numpy.nonzero(x)[0].size>1, ....) to filter 1D arrays of size < 2 or value == 0 or padValue
    
    # Should change from nonzero() to filter for the padValue specifically (NaN, eps, etc)
    
    #(1,0,0), #(-1,0,0),
    aDiags = reduce(lambda x,y: x+y, [a.tolist() for a in numpy.transpose(matrix,(1,2,0))])  
    matrixDiagonals.append( filter(lambda x: numpy.nonzero(x)[0].size>1, aDiags) )
                 
    #(0,1,0), #(0,-1,0),
    bDiags = reduce(lambda x,y: x+y, [a.tolist() for a in numpy.transpose(matrix,(0,2,1))])
    matrixDiagonals.append( filter(lambda x: numpy.nonzero(x)[0].size>1, bDiags) )
               
    #(0,0,1), #(0,0,-1), 
    cDiags = reduce(lambda x,y: x+y, [a.tolist() for a in numpy.transpose(matrix,(0,1,2))])
    matrixDiagonals.append( filter(lambda x: numpy.nonzero(x)[0].size>1, cDiags) )
                
    #(1,1,0),#(-1,-1,0),
    lowBound = -matrix.shape[0]+1
    highBound = matrix.shape[1]
    
    dDiags = reduce(lambda x,y: x+y, [matrix.diagonal(a,0,1).tolist() for a in xrange(lowBound, highBound)])
    matrixDiagonals.append( filter(lambda x: numpy.nonzero(x)[0].size>1, dDiags) )
                         
    #(1,0,1), #(-1,0-1),
    lowBound = -matrix.shape[0]+1
    highBound = matrix.shape[2]
      
    eDiags = reduce(lambda x,y: x+y, [matrix.diagonal(a,0,2).tolist() for a in xrange(lowBound, highBound)])
    matrixDiagonals.append( filter(lambda x: numpy.nonzero(x)[0].size>1, eDiags) )
        
    #(0,1,1), #(0,-1,-1),
    lowBound = -matrix.shape[1]+1
    highBound = matrix.shape[2]
    
    fDiags = reduce(lambda x,y: x+y, [matrix.diagonal(a,1,2).tolist() for a in xrange(lowBound, highBound)])
    matrixDiagonals.append( filter(lambda x: numpy.nonzero(x)[0].size>1, fDiags) )
                             
    #(1,-1,0), #(-1,1,0),    
    lowBound = -matrix.shape[0]+1
    highBound = matrix.shape[1]
    
    gDiags = reduce(lambda x,y: x+y, [matrix[:,::-1,:].diagonal(a,0,1).tolist() for a in xrange(lowBound, highBound)])
    matrixDiagonals.append( filter(lambda x: numpy.nonzero(x)[0].size>1, gDiags) )
     
    #(-1,0,1), #(1,0,-1),
    lowBound = -matrix.shape[0]+1
    highBound = matrix.shape[2]
    
    hDiags = reduce(lambda x,y: x+y, [matrix[:,:,::-1].diagonal(a,0,2).tolist() for a in xrange(lowBound, highBound)])
    matrixDiagonals.append( filter(lambda x: numpy.nonzero(x)[0].size>1, hDiags) )
                  
    #(0,1,-1), #(0,-1,1),
    lowBound = -matrix.shape[1]+1
    highBound = matrix.shape[2]
            
    iDiags = reduce(lambda x,y: x+y, [matrix[:,:,::-1].diagonal(a,1,2).tolist() for a in xrange(lowBound, highBound)])
    matrixDiagonals.append( filter(lambda x: numpy.nonzero(x)[0].size>1, iDiags) )
               
    #(1,1,1), #(-1,-1,-1)
    lowBound = -matrix.shape[0]+1
    highBound = matrix.shape[1]
    
    jDiags = [ numpy.diagonal(h,x,0,1).tolist() for h in [matrix.diagonal(a,0,1) for a in xrange(lowBound, highBound)] for x in xrange(-h.shape[0]+1, h.shape[1]) ]
    matrixDiagonals.append( filter(lambda x: numpy.nonzero(x)[0].size>1, jDiags) )
                             
    #(-1,1,-1), #(1,-1,1),
    lowBound = -matrix.shape[0]+1
    highBound = matrix.shape[1]
    
    kDiags = [ numpy.diagonal(h,x,0,1).tolist() for h in [matrix[:,::-1,:].diagonal(a,0,1) for a in xrange(lowBound, highBound)] for x in xrange(-h.shape[0]+1, h.shape[1]) ]
    matrixDiagonals.append( filter(lambda x: numpy.nonzero(x)[0].size>1, kDiags) )
                          
    #(1,1,-1), #(-1,-1,1),
    lowBound = -matrix.shape[0]+1
    highBound = matrix.shape[1]
    
    lDiags = [ numpy.diagonal(h,x,0,1).tolist() for h in [matrix[:,:,::-1].diagonal(a,0,1) for a in xrange(lowBound, highBound)] for x in xrange(-h.shape[0]+1, h.shape[1]) ]
    matrixDiagonals.append( filter(lambda x: numpy.nonzero(x)[0].size>1, lDiags) )
                         
    #(-1,1,1), #(1,-1,-1),
    lowBound = -matrix.shape[0]+1
    highBound = matrix.shape[1]
    
    mDiags = [ numpy.diagonal(h,x,0,1).tolist() for h in [matrix[:,::-1,::-1].diagonal(a,0,1) for a in xrange(lowBound, highBound)] for x in xrange(-h.shape[0]+1, h.shape[1]) ]
    matrixDiagonals.append( filter(lambda x: numpy.nonzero(x)[0].size>1, mDiags) )
                
    #[n for n in mDiags if numpy.nonzero(n)[0].size>1] instead of filter(lambda x: numpy.nonzero(x)[0].size>1, mDiags)?
    
    #Run-Length Encoding (rle) for the 13 list of diagonals (1 list per 3D direction/angle)
    for angle in xrange (len(matrixDiagonals)):
        P = P_rlgl[:,:,angle]     
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
         
    return (P_rlgl)

def CreateRLGL(Ng, Nr, grayLevels, parameterMatrix, parameterMatrixCoordinates, angles=13):  
    P_rlgl = numpy.zeros((Ng, Nr, angles))
    P_rlgl = calculate_rlgl(grayLevels, parameterMatrix, parameterMatrixCoordinates, angles, P_rlgl)
    
    #normalize like glcm?
    """
    for q in xrange(angles):
      P_rlgl[:,:,0,q] = P_rlgl[:,:,0,q]/(P_rlgl[:,:,0,q].sum())  
    """
    return P_rlgl    