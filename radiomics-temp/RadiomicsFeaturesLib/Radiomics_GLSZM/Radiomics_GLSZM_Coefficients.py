import numpy

def CalculateCoefficients(P_rlgl):
    sumP_rlgl = numpy.sum( numpy.sum(P_rlgl, 0), 0 )
    sumP_rlgl[sumP_rlgl==0] = 1
    
    pr = numpy.sum(P_rlgl, 0)
    pg = numpy.sum(P_rlgl, 1)
    
    ivector = numpy.arange(1, P_rlgl.shape[0] + 1)
    jvector = numpy.arange(1, P_rlgl.shape[1] + 1)
    
    Coefficients = {}
    Coefficients['sumP_rlgl'] = sumP_rlgl
    Coefficients['pr'] = pr
    Coefficients['pg'] = pg
    Coefficients['ivector'] = ivector
    Coefficients['jvector'] = jvector
    
    return Coefficients