import numpy
import collections
from radiomics import base

class RadiomicsFirstOrder(base.RadiomicsFeaturesBase):

    def __init__(self, matrix, matrixCoordinates, binwidth, pixelSpacing):
        self.matrix = matrix
        self.matrixCoordinates = matrixCoordinates
        self.targetVoxelArray = self.matrix[self.matrixCoordinates]

        self.binwidth = binwidth
        self.pixelSpacing = pixelSpacing

        self.featureNames = self.getFeatureNames()

        #self.InitializeFeatureVector()
        #for f in self.getFeatureNames():
        #  self.enabledFeatures[f] = True

    # TODO: add an option to instantiate the class that reuses initialization

    def initializeFeatureVector(self):
        self.prefix = "Stats_"
        self.radiomics_first_order_FeatureVector = collections.OrderedDict()
        self.radiomics_first_order_FeatureVector[self.prefix+"energy"] = "self.energyValue(self.targetVoxelArray)"
        self.radiomics_first_order_FeatureVector[self.prefix+"totalenergy"] = "self.totalEnergyValue(self.targetVoxelArray, self.pixelSpacing)"
        self.radiomics_first_order_FeatureVector[self.prefix+"entropy"] = "self.entropyValue(self.targetVoxelArray, self.binwidth)"
        self.radiomics_first_order_FeatureVector[self.prefix+"min"] = "self.minIntensity(self.targetVoxelArray)"
        self.radiomics_first_order_FeatureVector[self.prefix+"max"] = "self.maxIntensity(self.targetVoxelArray)"
        self.radiomics_first_order_FeatureVector[self.prefix+"mean"] = "self.meanIntensity(self.targetVoxelArray)"
        self.radiomics_first_order_FeatureVector[self.prefix+"median"] = "self.medianIntensity(self.targetVoxelArray)"
        self.radiomics_first_order_FeatureVector[self.prefix+"range"] = "self.rangeIntensity(self.targetVoxelArray)"
        self.radiomics_first_order_FeatureVector[self.prefix+"md"] = "self.meanDeviation(self.targetVoxelArray)"
        self.radiomics_first_order_FeatureVector[self.prefix+"rms"] = "self.rootMeanSquared(self.targetVoxelArray)"
        self.radiomics_first_order_FeatureVector[self.prefix+"std"] = "self.standardDeviation(self.targetVoxelArray)"
        self.radiomics_first_order_FeatureVector[self.prefix+"skewness"] = "self.skewnessValue(self.targetVoxelArray)"
        self.radiomics_first_order_FeatureVector[self.prefix+"kurtosis"] = "self.kurtosisValue(self.targetVoxelArray)"
        self.radiomics_first_order_FeatureVector[self.prefix+"var"] = "self.varianceValue(self.targetVoxelArray)"
        self.radiomics_first_order_FeatureVector[self.prefix+"uniformity"] = "self.uniformityValue(self.targetVoxelArray, self.binwidth)"

    def evaluateFeatures(self):
        for feature in self.radiomics_first_order_FeatureVector:
            try:
                self.radiomics_first_order_FeatureVector[feature] = eval(self.radiomics_first_order_FeatureVector[feature])
            except AttributeError:
                self.radiomics_first_order_FeatureVector[feature] = "Function Does Not Exist"
        return(self.radiomics_first_order_FeatureVector)

    def _getMomentFeatureValue(self, a, moment=1, axis=0):
      if moment == 1:
        return numpy.float64(0.0)
      else:
        mn = numpy.expand_dims(numpy.mean(a,axis), axis)
        s = numpy.power((a-mn), moment)
        return numpy.mean(s, axis)

    def getEnergyFeatureValue(self, targetVoxelArray):
      """
      Calculate Energy of the array

      Energy is defined as ...
      See reference below for details:
      """
      shiftedParameterArray = targetVoxelArray + 2000
      return (numpy.sum(shiftedParameterArray**2))

    def getTotalEnergyFeatureValue(self, targetVoxelArray, pixelSpacing):
        shiftedParameterArray = targetVoxelArray + 2000
        cubicMMPerVoxel = reduce(lambda x,y: x*y , pixelSpacing)
        return(cubicMMPerVoxel*numpy.sum(shiftedParameterArray**2))

    def getEntropyFeatureValue(self, targetVoxelArray, binwidth):
        ##check for binning centered at 0
        bincount = numpy.ceil((numpy.max(targetVoxelArray) - numpy.min(targetVoxelArray))/float(binwidth))
        bins = numpy.histogram(targetVoxelArray, bins=bincount)[0]
        bins = bins/float(bins.sum())
        return (-1.0 * numpy.sum(bins*numpy.where(bins!=0,numpy.log2(bins),0)))

    def getMinIntensityFeatureValue(self, targetVoxelArray):
        return (numpy.min(targetVoxelArray))

    def getMaxIntensityFeatureValue(self, targetVoxelArray):
        return (numpy.max(targetVoxelArray))

    def getMeanIntensityFeatureValue(self, targetVoxelArray):
        return (numpy.mean(targetVoxelArray))

    def getMedianIntensityFeatureValue (self, targetVoxelArray):
        return (numpy.median(targetVoxelArray))

    def getRangeIntensityFeatureValue (self, targetVoxelArray):
        return (numpy.max(targetVoxelArray) - numpy.min(targetVoxelArray))

    def getMeanDeviationFeatureValue(self, targetVoxelArray):
        return ( numpy.mean(numpy.absolute( (numpy.mean(targetVoxelArray) - targetVoxelArray) )) )

    def getRootMeanSquaredFeatureValue(self, targetVoxelArray):
        shiftedParameterArray = targetVoxelArray + 2000
        return ( numpy.sqrt((numpy.sum(shiftedParameterArray**2))/float(shiftedParameterArray.size)) )

    def getStandardDeviationFeatureValue(self, targetVoxelArray):
        return (numpy.std(targetVoxelArray))

    def getSkewnessValueFeatureValue(self, a, axis=0):
        m2 = _moment(a, 2, axis)
        m3 = _moment(a, 3, axis)

        zero = (m2 == 0)
        vals = numpy.where(zero, 0, m3 / m2**1.5)

        if vals.ndim == 0:
            return vals.item()

        return vals

    def getKurtosisFeatureValue(self, a, axis=0):
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

    def getVarianceFeatureValue(self, targetVoxelArray):
        return (numpy.std(targetVoxelArray)**2)

    def getUniformityFeatureValue(self, targetVoxelArray, binwidth):
        bincount = numpy.ceil((numpy.max(targetVoxelArray) - numpy.min(targetVoxelArray))/float(binwidth))
        bins = numpy.histogram(targetVoxelArray, bins=bincount)[0]
        bins = bins/float(bins.sum())
        return (numpy.sum(bins**2))
