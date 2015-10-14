import numpy
import collections
import pdb
               
def autocorrelation(P_glcm, i, j):
    ac = numpy.sum(numpy.sum(P_glcm * (i*j)[:,:,None,None], 0 ), 0 )
    return (ac.mean())
   
def clusterProminence(P_glcm, i, j, ux, uy):
    cp = numpy.sum( numpy.sum( (P_glcm * (((i+j)[:,:,None,None] - ux[None,None,:,:] - uy[None,None,:,:])**4)), 0 ), 0 )
    return (cp.mean())
    
def clusterShade(P_glcm, i, j, ux, uy):
    cs = numpy.sum( numpy.sum( (P_glcm * (((i+j)[:,:,None,None] - ux[None,None,:,:] - uy[None,None,:,:])**3)), 0 ), 0 )
    return (cs.mean())
  
def clusterTendency(P_glcm, i, j, ux, uy):
    ct = numpy.sum( numpy.sum( (P_glcm * (((i+j)[:,:,None,None] - ux[None,None,:,:] - uy[None,None,:,:])**2)), 0 ), 0 )
    return (ct.mean())
  
def contrast(P_glcm, i, j):
    cont = numpy.sum( numpy.sum( (P_glcm * ((numpy.abs(i-j))[:,:,None,None]**2)), 0 ), 0 )
    return (cont.mean())
  
def correlation(P_glcm, i, j, ux, uy, sigx, sigy):
    corm = numpy.sum(numpy.sum(P_glcm*(i[:,:,None,None]-ux[None,None,:,:])*(j[:,:,None,None]-uy[None,None,:,:]), 0 ), 0 )
    corr = corm/(sigx*sigy)
    return (corr.mean())
  
def differenceEntropy(pxSuby, eps):
    difent = (-1) * numpy.sum( (pxSuby*numpy.log(pxSuby+eps)), 0 )
    return (difent.mean())
  
def dissimilarity(P_glcm, i, j):
    dis = numpy.sum( numpy.sum( (P_glcm * (numpy.abs(i-j))[:,:,None,None]), 0 ), 0 )
    return (dis.mean())
    
def energy(P_glcm):
    ene = numpy.sum( numpy.sum( (P_glcm**2), 0), 0 )
    return (ene.mean())
  
def entropy(P_glcm, eps):
    ent = (-1) * numpy.sum( numpy.sum( (P_glcm * numpy.log(P_glcm+eps)), 0), 0) 
    return (ent.mean())
  
def homogeneity1(P_glcm, i, j):
    homo1 = numpy.sum( numpy.sum( (P_glcm / (1 + (numpy.abs(i-j))[:,:,None,None])), 0 ), 0 )
    return (homo1.mean())
    
def homogeneity2(P_glcm, i, j):
    homo2 = numpy.sum( numpy.sum( (P_glcm / (1 + (numpy.abs(i-j))[:,:,None,None]**2)), 0 ), 0 )
    return (homo2.mean())
  
def imc1(HXY, HXY1, HX, HY):
    imc1 = (HXY - HXY1)/numpy.max(([HX,HY]),0)
    return (imc1.mean())  
 
def imc2():
    # Produces Error
    
    imc2[g,a] = ( 1-numpy.e**(-2*(HXY2[g,a]-HXY[g,a])) )**(0.5) #nan value too high
    # matlab:(1-exp(-2*(hxy2-hxy)))^0.5; 
    
    #produces Nan(square root of a negative)
    #exponent = decimal.Decimal( -2*(HXY2[g,a]-self.HXY[g,a]) )      
    #imc2.append( ( decimal.Decimal(1)-decimal.Decimal(numpy.e)**(exponent) )**(decimal.Decimal(0.5)) )
    
    #if meanFlag:
      #return (homo2.mean())
    #else:
      #return homo2
   
def idmn(P_glcm, i, j, Ng):
    idmn = numpy.sum( numpy.sum( (P_glcm / (1 + (((numpy.abs(i-j))[:,:,None,None]**2)/(Ng**2)))), 0 ), 0 )
    return (idmn.mean())
   
def idn(P_glcm, i, j, Ng):
    idn  = numpy.sum( numpy.sum( (P_glcm / (1 + ((numpy.abs(i-j))[:,:,None,None]/Ng))), 0 ), 0 )
    return (idn.mean())
   
def inverseVariance(P_glcm, i, j, Ng):
    maskDiags = numpy.ones((numpy.abs(i-j)).shape, dtype = bool)
    maskDiags[numpy.diag_indices(Ng)] = False         
    inv = numpy.sum( (P_glcm[maskDiags] / ((numpy.abs(i-j))[:,:,None,None]**2)[maskDiags]), 0 )
    return (inv.mean())
     
def maximumProbability(P_glcm):
    maxprob = P_glcm.max(0).max(0)
    return (maxprob.mean())
   
def sumAverage(pxAddy, kValuesSum):
    sumavg =  numpy.sum( (kValuesSum[:,None,None]*pxAddy), 0 )
    return (sumavg.mean())
   
def sumEntropy(pxAddy, eps):
    sumentr = (-1) * numpy.sum( (pxAddy * numpy.log(pxAddy+eps)), 0 )
    return (sumentr.mean())
   
def sumVariance(pxAddy, kValuesSum):
    sumvar = numpy.sum( (pxAddy*((kValuesSum[:,None,None] - kValuesSum[:,None,None]*pxAddy)**2)), 0 )
    return (sumvar.mean())
 
def sumSquares(P_glcm, j, u):
    # Also known as Variance
    ss = numpy.sum( numpy.sum( (P_glcm * ((j[:,:,None,None]-u[None,None,:,:])**2)), 0 ), 0 )
    return (ss.mean())