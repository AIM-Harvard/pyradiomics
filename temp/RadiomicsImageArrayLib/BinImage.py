import numpy

def BinImage(binwidth, parameterValues, parameterMatrix, parameterMatrixCoordinates):
    lowBound = min(parameterValues) - binwidth
    highBound = max(parameterValues)+binwidth
    
    binedges = numpy.union1d(numpy.arange(0,lowBound,-binwidth), numpy.arange(0,highBound,binwidth))
    histogram = numpy.histogram(parameterValues, bins=binedges)
    parameterMatrix[parameterMatrixCoordinates] = numpy.digitize(parameterValues,binedges)
    
    return parameterMatrix, histogram