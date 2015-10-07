import numpy
import collections

def _moment(a, moment=1, axis=0):
    if moment == 1:
        return numpy.float64(0.0)
    else:
        mn = numpy.expand_dims(numpy.mean(a,axis), axis)
        s = numpy.power((a-mn), moment)
        return numpy.mean(s, axis)
 
def energyValue(targetVoxelArray):
    shiftedParameterArray = targetVoxelArray + 2000
    return (numpy.sum(shiftedParameterArray**2))

def totalEnergyValue(targetVoxelArray, pixelSpacing):
    shiftedParameterArray = targetVoxelArray + 2000
    cubicMMPerVoxel = reduce(lambda x,y: x*y , pixelSpacing)
    return(cubicMMPerVoxel*numpy.sum(shiftedParameterArray**2))
  
def entropyValue(targetVoxelArray, binwidth):
    ##check for binning centered at 0
    bincount = numpy.ceil((numpy.max(targetVoxelArray) - numpy.min(targetVoxelArray))/float(binwidth))
    bins = numpy.histogram(targetVoxelArray, bins=bincount)[0]
    bins = bins/float(bins.sum())
    return (-1.0 * numpy.sum(bins*numpy.where(bins!=0,numpy.log2(bins),0)))
  
def minIntensity(targetVoxelArray):
    return (numpy.min(targetVoxelArray))

def maxIntensity(targetVoxelArray):
    return (numpy.max(targetVoxelArray))
    
def meanIntensity(targetVoxelArray):
    return (numpy.mean(targetVoxelArray))
    
def medianIntensity (targetVoxelArray):
    return (numpy.median(targetVoxelArray))
    
def rangeIntensity (targetVoxelArray):
    return (numpy.max(targetVoxelArray) - numpy.min(targetVoxelArray))

def meanDeviation(targetVoxelArray):
    return ( numpy.mean(numpy.absolute( (numpy.mean(targetVoxelArray) - targetVoxelArray) )) )

def rootMeanSquared(targetVoxelArray):
    shiftedParameterArray = targetVoxelArray + 2000  
    return ( numpy.sqrt((numpy.sum(shiftedParameterArray**2))/float(shiftedParameterArray.size)) )

def standardDeviation(targetVoxelArray):
    return (numpy.std(targetVoxelArray))

def skewnessValue(a, axis=0):
    m2 = _moment(a, 2, axis)
    m3 = _moment(a, 3, axis)
    
    zero = (m2 == 0)
    vals = numpy.where(zero, 0, m3 / m2**1.5)
  
    if vals.ndim == 0:
        return vals.item()
        
    return vals

def kurtosisValue(a, axis=0):
    m2 = _moment(a,2,axis)
    m4 = _moment(a,4,axis)
    zero = (m2 == 0)
    
    olderr = numpy.seterr(all='ignore')
    
    try:
        vals = numpy.where(zero, 0, m4 / m2**2.0)
    finally:
        numpy.seterr(**olderr)
    if vals.ndim == 0:
        vals = vals.item()
    
    return vals
    
def varianceValue(targetVoxelArray):
    return (numpy.std(targetVoxelArray)**2)

def uniformityValue(targetVoxelArray, binwidth):
    bincount = numpy.ceil((numpy.max(targetVoxelArray) - numpy.min(targetVoxelArray))/float(binwidth))
    bins = numpy.histogram(targetVoxelArray, bins=bincount)[0]
    bins = bins/float(bins.sum())
    return (numpy.sum(bins**2))