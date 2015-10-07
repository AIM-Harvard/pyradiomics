import numpy

def shortRunEmphasis(pr, jvector, sumP_rlgl):
    try:
        #sre = numpy.sum( numpy.sum( (P_rlgl/(self.j[:,:,None]**2)) , 0 ), 0 ) / (sumP_rlgl[None,None,:])
        sre = numpy.sum( (pr/(jvector[:,None]**2)), 0 ) / sumP_rlgl
    except ZeroDivisionError:
        sre = 0    
    return (sre.mean())

def longRunEmphasis(pr, jvector, sumP_rlgl):
    try:
        #lre = numpy.sum( numpy.sum( (P_rlgl*(self.j[:,:,None]**2)) , 0 ), 0 ) / (sumP_rlgl[None,None,:])
        lre = numpy.sum( (pr*(jvector[:,None]**2)), 0 ) / sumP_rlgl
    except ZeroDivisionError:
        lre = 0
    return (lre.mean())

def grayLevelNonUniformity(pg, sumP_rlgl):
    try:
        gln = numpy.sum( (pg**2) , 0 ) / sumP_rlgl   
    except ZeroDivisionError:
        gln = 0  
    return (gln.mean())
  
def runLengthNonUniformity(pr, sumP_rlgl):
    try:
        rln = numpy.sum( (pr**2) , 0 ) / sumP_rlgl 
    except ZeroDivisionError:
        rln = 0
    return (rln.mean())

def runPercentage(P_rlgl, Np):
    try:
        rp = numpy.sum( numpy.sum( (P_rlgl/(Np)) , 0 ), 0 )
    except ZeroDivisionError:
        rp = 0
    return (rp.mean())

def lowGrayLevelRunEmphasis(pg, ivector, sumP_rlgl):
    try:
        #lglre = numpy.sum( numpy.sum( (P_rlgl/(self.i[:,:,None]**2)) , 0 ), 0 )/ (sumP_rlgl[None,None,:])
        lglre = numpy.sum( (pg/(ivector[:,None]**2)) , 0 ) / sumP_rlgl
    except ZeroDivisionError:
        lglre = 0 
    return (lglre.mean())
    
def highGrayLevelRunEmphasis(pg, ivector, sumP_rlgl):
    try:
        #hglre = numpy.sum( numpy.sum( (P_rlgl*(self.i[:,:,None]**2)) , 0 ), 0 ) / (sumP_rlgl[None,None,:])
        hglre = numpy.sum( (pg*(ivector[:,None]**2)) , 0 ) / sumP_rlgl
    except ZeroDivisionError:
        hglre = 0 
    return (hglre.mean())
  
def shortRunLowGrayLevelEmphasis(P_rlgl, ivector, jvector, sumP_rlgl):
    try:
        srlgle = numpy.sum( numpy.sum( (P_rlgl/((ivector[:,None,None]**2)*(jvector[None,:,None]**2))) , 0 ), 0 ) / sumP_rlgl
    except ZeroDivisionError:
        srlgle = 0
    return (srlgle.mean())
   
def shortRunHighGrayLevelEmphasis(P_rlgl, ivector, jvector, sumP_rlgl):
    try:
        srhgle = numpy.sum( numpy.sum( (P_rlgl*(ivector[:,None,None]**2)/(jvector[None,:,None]**2)) , 0 ), 0 ) / sumP_rlgl
    except ZeroDivisionError:
        srhgle = 0   
    return (srhgle.mean())
      
def longRunLowGrayLevelEmphasis(P_rlgl, ivector, jvector, sumP_rlgl):
    try:
        lrlgle = numpy.sum( numpy.sum( (P_rlgl*(jvector[None,:,None]**2)/(ivector[:,None,None]**2)) , 0 ), 0 ) / sumP_rlgl
    except ZeroDivisionError:
        lrlgle = 0
    return (lrlgle.mean())
          
def longRunHighGrayLevelEmphasis(P_rlgl, ivector, jvector, sumP_rlgl):
    try:
        lrhgle = numpy.sum( numpy.sum( (P_rlgl*((jvector[None,:,None]**2)*(ivector[:,None,None]**2))) , 0 ), 0 ) / sumP_rlgl
    except ZeroDivisionError:
        lrhgle = 0 
    return (lrhgle.mean()) 