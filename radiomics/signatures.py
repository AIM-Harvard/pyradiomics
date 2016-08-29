# -*- coding: utf-8 -*-
import os
import collections
import numpy as np
import SimpleITK as sitk

import pkgutil
import inspect
import radiomics
from radiomics import base, imageoperations


class RadiomicsSignature():

    def __init__(self, **kwargs):
        self.featureClasses = self.getFeatureClasses()

        self.kwargs = kwargs

        # Try get values for interpolation and verbose. If not present in kwargs, use defaults
        self.resampledPixelSpacing = self.kwargs.get('resampledPixelSpacing', None) #  no resampling by default
        self.interpolator = self.kwargs.get('interpolator', sitk.sitkBSpline)
        self.padDistance = self.kwargs.get('padDistance', 5)
        self.verbose = self.kwargs.get('verbose', True)

        self.inputImages = {}
        self.inputImages['original'] = {}
        self.inputImages['log'] = {'sigma': np.arange(5.,0.,-.5)}
        self.inputImages['wavelet'] = {}

        self.enabledFeatures = {}

    def enableInputImages(self, **inputImages):
        self.inputImages = inputImages

    def enableAllFeatures(self):
        """
        Enable all classes and all features.
        """

        self.enabledFeatures = {}

    def enableFeatureClassByName(self, featureClass, enabled= True):
        """
        Enable or disable all features in given class.
        """

        if enabled:
            if featureClass in self.enabledFeatures: del self.enabledFeatures[featureClass]
        else:
            self.enabledFeatures[featureClass] = []

    def enableFeaturesByName(self, **enabledFeatures):
        """
        Specify which features to enable. Key is feature class name, value is a list of enabled feature names.

        To disable all features for a class, provide the class name with an empty list as value.
        Only settings for feature classes specified in enabledFeatures.keys are updated.
        To enable the entire class, use enableFeatureClassByName instead.
        """

        self.enabledFeatures.update(enabledFeatures)

    def computeSignature(self, imageFilepath, maskFilepath):
        featureVector = collections.OrderedDict()

        image, mask = self.loadImage(imageFilepath, maskFilepath)

        if 'original' in self.inputImages:
            args = self.kwargs.copy()
            args.update(self.inputImages['original'])
            featureVector.update(self.computeFeatures(image, mask, **args))
        if 'log' in self.inputImages:
            if self.verbose: print "\tComputing LoG"
            args = self.kwargs.copy()
            args.update(self.inputImages['log'])
            featureVector.update(self.computeLoG(image, mask, **args))
        if 'wavelet' in self.inputImages:
            if self.verbose: print "\tComputing Wavelet"
            args = self.kwargs.copy()
            args.update(self.inputImages['wavelet'])
            featureVector.update(self.computeWavelet(image, mask, **args))

        return featureVector

    def loadImage(self, ImageFilePath, MaskFilePath):
        if isinstance(ImageFilePath, basestring) and os.path.exists(ImageFilePath):
            image = sitk.ReadImage(ImageFilePath)
        elif isinstance(ImageFilePath, sitk.SimpleITK.Image):
            image = ImageFilePath
        else:
            if self.verbose: print "Error reading image Filepath or SimpleITK object"
            image = None

        if isinstance(MaskFilePath, basestring) and os.path.exists(MaskFilePath):
            mask = sitk.ReadImage(MaskFilePath)
        elif isinstance(ImageFilePath, sitk.SimpleITK.Image):
            mask = MaskFilePath
        else:
            if self.verbose: print "Error reading mask Filepath or SimpleITK object"
            mask = None

        """
        imageArray = sitk.GetArrayFromImage(image)
        maskArray = sitk.GetArrayFromImage(mask)
        tumorVoxels = imageArray[np.where(maskArray==1)]
        mean = np.mean(tumorVoxels)
        std = np.std(tumorVoxels)
        minBound = mean - 3*std
        maxBound = mean + 3*std
        maskArray[np.where(imageArray < minBound)] = 0
        maskArray[np.where(imageArray > maxBound)] = 0
        mask_normalized = sitk.GetImageFromArray(maskArray)
        mask_normalized.CopyInformation(mask)
        mask = mask_normalized
        """

        if self.interpolator != None and self.resampledPixelSpacing != None:
            image, mask = imageoperations.resampleImage(image, mask, self.resampledPixelSpacing, self.interpolator, self.padDistance)
        else:
            image,mask = imageoperations.cropTumorMaskToCube(image, mask, self.padDistance)

        return image, mask

    def computeFeatures(self, image, mask, **kwargs):

        featureVector = collections.OrderedDict()

        for featureClassName in self.getFeatureClassNames():
            if featureClassName not in self.enabledFeatures.keys():
                featureClass = self.featureClasses[featureClassName](image, mask, **kwargs)
                featureClass.enableAllFeatures()

            elif len(self.enabledFeatures[featureClassName]) > 0:
                featureClass = self.featureClasses[featureClassName](image, mask, **kwargs)
                for enabledFeature in self.enabledFeatures[featureClassName]:
                    featureClass.enableFeatureByName(enabledFeature)

            else:
                featureClass = None

            if featureClass != None:
                if self.verbose: print "\tComputing %s" %(featureClassName)
                featureClass.calculateFeatures()
                for (featureName, featureValue) in featureClass.featureValues.iteritems():
                    shapeFeatureName = "%s_%s" %(featureClassName, featureName)
                    featureVector[shapeFeatureName] = featureValue

        return featureVector

    def computeLoG(self, image, mask, **kwargs):

        logFeatureVector = collections.OrderedDict()
        sigmaValues = kwargs.get('sigma', [])

        for sigma in sigmaValues:

            logSigmaFeatureVector = collections.OrderedDict()

            logImage = imageoperations.applyLoG(image, sigmaValue=sigma)

            logSigmaFeatureVector.update( self.computeFeatures(logImage, mask, **kwargs))

            for (featureName, featureValue) in logSigmaFeatureVector.iteritems():

                laplacianFeatureName = "log_sigma_%s_mm_3D_%s" %(str(sigma).replace('.','_'), featureName)
                logFeatureVector[laplacianFeatureName] = featureValue

        return logFeatureVector

    def computeWavelet(self, image, mask, **kwargs):

        waveletArgs = {}
        waveletArgs['wavelet'] = kwargs.get('wavelet', 'coif1')
        waveletArgs['level'] = kwargs.get('level', 1)
        waveletArgs['start_level'] = kwargs.get('start_level', 0)

        waveletFeatureVector = collections.OrderedDict()

        approx, ret = imageoperations.swt3(image, **waveletArgs)

        for idx, wl in enumerate(ret, start= 1):
            for decompositionName, decompositionImage in wl.items():
                waveletDecompositionFeatureVector = collections.OrderedDict()
                waveletDecompositionFeatureVector.update( self.computeFeatures(decompositionImage, mask, **kwargs) )

                for (featureName, featureValue) in waveletDecompositionFeatureVector.iteritems():
                    if idx == 1: waveletFeatureName = "wavelet_%s_%s" %(decompositionName, featureName)
                    else: waveletFeatureName = "wavelet%s_%s_%s" %(idx, decompositionName, featureName)
                    waveletFeatureVector[waveletFeatureName] = featureValue

        waveletDecompositionFeatureVector = collections.OrderedDict()
        waveletDecompositionFeatureVector.update( self.computeFeatures(approx, mask, **kwargs) )

        for (featureName, featureValue) in waveletDecompositionFeatureVector.iteritems():
            if len(ret) == 1: waveletFeatureName = "wavelet_LLL_%s" %(featureName)
            else: waveletFeatureName = "wavelet%s_LLL_%s" %(len(ret), featureName)
            waveletFeatureVector[waveletFeatureName] = featureValue

        return waveletFeatureVector

    def getFeatureClasses(self):
        """
        Iterates over all modules of the radiomics package using pkgutil and subsequently imports those modules.

        Return a dictionary of all modules containing featureClasses, with modulename as key, abstract
        class object of the featureClass as value. Assumes only one featureClass per module

        This is achieved by inspect.getmembers. Modules are added if it contains a memeber that is a class,
        with name starting with 'Radiomics' and is inherited from radiomics.base.RadiomicsFeaturesBase.
        """

        featureClasses = {}
        for _,mod,_ in pkgutil.iter_modules(radiomics.__path__):
            __import__('radiomics.' + mod)
            attributes = inspect.getmembers(eval('radiomics.' + mod), inspect.isclass)
            for a in attributes:
                if a[0].startswith('Radiomics'):
                    if radiomics.base.RadiomicsFeaturesBase in inspect.getmro(a[1])[1:]:
                        featureClasses[mod] = a[1]


        return featureClasses

    def getFeatureClassNames(self):
        return self.featureClasses.keys()
