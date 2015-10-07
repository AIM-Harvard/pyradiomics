import numpy
import collections
import pdb

#check if ivector and jvector can be replaced
def CalculateCoefficients(Ng, P_glcm): 
    eps = numpy.spacing(1)
    
    ##Calculate GLCM Coefficients
    ivector = numpy.arange(1,Ng+1) #shape = (Ng, distances.size, angles)
    jvector = numpy.arange(1,Ng+1) #shape = (Ng, distances.size, angles)
    
    i,j = numpy.meshgrid(numpy.arange(1,P_glcm.shape[0]+1), numpy.arange(1,P_glcm.shape[1]+1))
    
    prodMatrix = numpy.multiply.outer(ivector, jvector) #shape = (Ng, Ng)
    sumMatrix = numpy.add.outer(ivector, jvector) #shape = (Ng, Ng) #sum-1 so 111 in corner
    diffMatrix = numpy.absolute(numpy.subtract.outer(ivector, jvector)) #shape = (Ng, Ng) #absval(subtract)+1 so 1-1 and 56-56 in corners
    
    kValuesSum = numpy.arange(2, (Ng*2)+1) - 1 #shape = (2*Ng-1) #1-111
    kValuesDiff = numpy.arange(0,Ng) + 1 #shape = (Ng-1) #0-55
    
    # shape = (distances.size, angles)
    u = P_glcm.mean(0).mean(0) 
    # marginal row probabilities #shape = (Ng, distances.size, angles)
    px = P_glcm.sum(1)
    # marginal column probabilities #shape = (Ng, distances.size, angles)
    py = P_glcm.sum(0) 
    
    # shape = (distances.size, angles)
    ux = numpy.sum(numpy.sum( ivector[None,:,None,None]*P_glcm, 0 ), 0 )
    uy = numpy.sum(numpy.sum( jvector[:,None,None,None]*P_glcm, 0 ), 0 )
    
    # shape = (distances.size, angles)
    sigx = numpy.sum(numpy.sum(P_glcm*((i[:,:,None,None]-ux[None,None,:,:])**2), 0 ), 0 )**0.5 
    # shape = (distances.size, angles)
    sigy = numpy.sum(numpy.sum(P_glcm*((j[:,:,None,None]-uy[None,None,:,:])**2), 0 ), 0 )**0.5
    
    # shape = (2*Ng-1, distances.size, angles) #111 long
    pxAddy = numpy.array([ numpy.sum(P_glcm[i+j == k], 0) for k in kValuesSum ]) 
    # shape = (Ng, distances.size, angles) #56 long
    pxSuby = numpy.array([ numpy.sum(P_glcm[numpy.abs(i-j) == k], 0) for k in kValuesDiff ])
    
    # entropy of px #shape = (distances.size, angles)
    HX = (-1) * numpy.sum( (px * numpy.log(px+eps)), 0)
    # entropy of py #shape = (distances.size, angles)
    HY = (-1) * numpy.sum( (py * numpy.log(py+eps)), 0)
    # shape = (distances.size, angles)
    HXY = (-1) * numpy.sum( numpy.sum( (P_glcm * numpy.log(P_glcm+eps)), 0 ), 0 )
        
    HXY1 = (-1) * numpy.sum( numpy.sum( (P_glcm * numpy.log(P_glcm+eps)), 0 ), 0) #shape = (distances.size, angles)  
    HXY2 = (-1) * numpy.sum( numpy.sum( ((px*py) * numpy.log(px*py+eps)), 0 ), 0) #shape = (distances.size, angles)
    
    Coefficients = {}
    Coefficients['Ng'] = Ng
    Coefficients['eps'] = eps
    Coefficients['i'] = i
    Coefficients['j'] = j
    Coefficients['kValuesSum'] = kValuesSum
    Coefficients['kValuesDiff'] = kValuesDiff
    Coefficients['u'] = u
    Coefficients['px'] = px
    Coefficients['py'] = py
    Coefficients['ux'] = ux
    Coefficients['uy'] = uy
    Coefficients['sigx'] = sigx
    Coefficients['sigy'] = sigy
    Coefficients['pxAddy'] = pxAddy
    Coefficients['pxSuby'] = pxSuby
    Coefficients['HX'] = HX
    Coefficients['HY'] = HY
    Coefficients['HXY'] = HXY
    Coefficients['HXY1'] = HXY1
    Coefficients['HXY2'] = HXY2

    return Coefficients