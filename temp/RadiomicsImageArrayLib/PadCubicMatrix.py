import numpy
import operator

def PadCubicMatrix(a, matrixCoordinates, dims):
    # pads matrix 'a' with zeros and resizes 'a' to a cube with dimensions increased to the next greatest power of 2
    # numpy version 1.7 has numpy.pad function
         
    # center coordinates onto padded matrix  # consider padding with NaN or eps = numpy.spacing(1)
    pad = tuple(map(operator.div, tuple(map(operator.sub, dims, a.shape)), ([2,2,2])))
    matrixCoordinatesPadded = tuple(map(operator.add, matrixCoordinates, pad))
    matrix2 = numpy.zeros(dims)
    matrix2[matrixCoordinatesPadded] = a[matrixCoordinates]
    return (matrix2, matrixCoordinatesPadded)