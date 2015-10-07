import numpy
import operator

def PadTumorMaskToCube(imageNodeArray, labelNodeArray):
    targetVoxelsCoordinates = numpy.where(labelNodeArray != 0)
    ijkMinBounds = numpy.min(targetVoxelsCoordinates, 1)
    ijkMaxBounds = numpy.max(targetVoxelsCoordinates, 1) 
    matrix = numpy.zeros(ijkMaxBounds - ijkMinBounds + 1)
    matrixCoordinates = tuple(map(operator.sub, targetVoxelsCoordinates, tuple(ijkMinBounds)))
    matrix[matrixCoordinates] = imageNodeArray[targetVoxelsCoordinates].astype('int64')
    return(matrix, matrixCoordinates)