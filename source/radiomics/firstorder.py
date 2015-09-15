import numpy
import collections
from radiomics import base, preprocessing
import SimpleITK as sitk

class RadiomicsFirstOrder(base.RadiomicsFeaturesBase):

    def __init__(self, inputImage, inputMask):
        super(RadiomicsFirstOrder,self).__init__(inputImage,inputMask)

        self.pixelSpacing = inputImage.GetSpacing()

        #self.featureNames = self.getFeatureNames()

        self.imageArray = sitk.GetArrayFromImage(inputImage)
        self.maskArray = sitk.GetArrayFromImage(inputMask)

        (self.matrix, self.matrixCoordinates) = \
          preprocessing.RadiomicsHelpers.padTumorMaskToCube(self.imageArray,self.maskArray)

        self.targetVoxelArray = self.matrix[self.matrixCoordinates]

        #self.InitializeFeatureVector()
        #for f in self.getFeatureNames():
        #  self.enabledFeatures[f] = True

        # TODO: add an option to instantiate the class that reuses initialization

    def _moment(self, a, moment=1, axis=0):
      if moment == 1:
        return numpy.float64(0.0)
      else:
        mn = numpy.expand_dims(numpy.mean(a,axis), axis)
        s = numpy.power((a-mn), moment)
        return numpy.mean(s, axis)

    def getEnergyFeatureValue(self):
      """
      Calculate Energy of the array

      Energy is defined as ...
      See reference below for details:
      """
      shiftedParameterArray = self.targetVoxelArray + 2000
      return (numpy.sum(shiftedParameterArray**2))

    def getTotalEnergyFeatureValue(self):
        shiftedParameterArray = self.targetVoxelArray + 2000
        cubicMMPerVoxel = reduce(lambda x,y: x*y , self.pixelSpacing)
        return(cubicMMPerVoxel*numpy.sum(shiftedParameterArray**2))

    def getEntropyFeatureValue(self):
        ##check for binning centered at 0
        bincount = numpy.ceil((numpy.max(self.targetVoxelArray) - numpy.min(self.targetVoxelArray))/float(self.binWidth))
        bins = numpy.histogram(self.targetVoxelArray, bins=bincount)[0]
        bins = bins/float(bins.sum())
        return (-1.0 * numpy.sum(bins*numpy.where(bins!=0,numpy.log2(bins),0)))

    def getMinIntensityFeatureValue(self):
        return (numpy.min(self.targetVoxelArray))

    def getMaxIntensityFeatureValue(self):
        return (numpy.max(self.targetVoxelArray))

    def getMeanIntensityFeatureValue(self):
        return (numpy.mean(self.targetVoxelArray))

    def getMedianIntensityFeatureValue (self):
        return (numpy.median(self.targetVoxelArray))

    def getRangeIntensityFeatureValue (self):
        return (numpy.max(self.targetVoxelArray) - numpy.min(self.targetVoxelArray))

    def getMeanDeviationFeatureValue(self):
        return ( numpy.mean(numpy.absolute( (numpy.mean(self.targetVoxelArray) - self.targetVoxelArray) )) )

    def getRootMeanSquaredFeatureValue(self):
        shiftedParameterArray = self.targetVoxelArray + 2000
        return ( numpy.sqrt((numpy.sum(shiftedParameterArray**2))/float(shiftedParameterArray.size)) )

    def getStandardDeviationFeatureValue(self):
        return (numpy.std(self.targetVoxelArray))

    def getSkewnessValueFeatureValue(self, axis=0):
        m2 = self._moment(self.targetVoxelArray, 2, axis)
        m3 = self._moment(self.targetVoxelArray, 3, axis)

        zero = (m2 == 0)
        vals = numpy.where(zero, 0, m3 / m2**1.5)

        if vals.ndim == 0:
            return vals.item()

        return vals

    def getKurtosisFeatureValue(self, axis=0):
        m2 = self._moment(self.targetVoxelArray,2,axis)
        m4 = self._moment(self.targetVoxelArray,4,axis)
        zero = (m2 == 0)

        olderr = numpy.seterr(all='ignore')

        try:
            vals = numpy.where(zero, 0, m4 / m2**2.0)
        finally:
            numpy.seterr(**olderr)
        if vals.ndim == 0:
            vals = vals.item()

        return vals

    def getVarianceFeatureValue(self):
        return (numpy.std(self.targetVoxelArray)**2)

    def getUniformityFeatureValue(self):
        bincount = numpy.ceil((numpy.max(self.targetVoxelArray) - numpy.min(self.targetVoxelArray))/float(self.binWidth))
        bins = numpy.histogram(self.targetVoxelArray, bins=bincount)[0]
        bins = bins/float(bins.sum())
        return (numpy.sum(bins**2))
