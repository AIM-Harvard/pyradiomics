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
        self.verbose = self.kwargs.get('verbose', True)

        self.inputImages = {}
        self.inputImages['original'] = {}
        self.inputImages['log'] = {'sigma': np.arange(5.,0.,-.5)}
        self.inputImages['wavelet'] = {}

        self.enabledFeatures = {}
        for featureClassName in self.getFeatureClassNames():
            self.enabledFeatures[featureClassName] = []

    def enableInputImages(self, **inputImages):
        self.inputImages = inputImages

    def enableAllFeatures(self):
        """
        Enable all classes and all features.
        """

        for featureClassName in self.getFeatureClassNames():
            self.enabledFeatures[featureClassName] = []

    def disableAllFeatures(self):
        """
        Disable all classes.
        """

        self.enabledFeatures = {}

    def enableFeatureClassByName(self, featureClass, enabled= True):
        """
        Enable or disable all features in given class.
        """

        if enabled:
            self.enabledFeatures[featureClass] = []
        else:
            if featureClass in self.enabledFeatures: del self.enabledFeatures[featureClass]

    def enableFeaturesByName(self, **enabledFeatures):
        """
        Specify which features to enable. Key is feature class name, value is a list of enabled feature names.

        To enable all features for a class, provide the class name with an empty list as value.
        Settings for feature classes specified in enabledFeatures.keys are updated, settings for feature classes
        not yet present in enabledFeatures.keys are added.
        To disable the entire class, use disableAllFeatures or enableFeatureClassByName instead.
        """

        self.enabledFeatures.update(enabledFeatures)

    def computeSignature(self, imageFilepath, maskFilepath):
        featureVector = collections.OrderedDict()

        image, mask = self.loadImage(imageFilepath, maskFilepath)

        if 'original' in self.inputImages:
            if self.verbose: print "Computing Original"
            args = self.kwargs.copy()
            args.update(self.inputImages['original'])
            featureVector.update(self.computeFeatures(image, mask, 'original', **args))
        if 'log' in self.inputImages:
            if self.verbose: print "Computing LoG"
            args = self.kwargs.copy()
            args.update(self.inputImages['log'])
            featureVector.update(self.computeLoG(image, mask, **args))
        if 'wavelet' in self.inputImages:
            if self.verbose: print "Computing Wavelet"
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
            image, mask = imageoperations.resampleImage(image, mask, self.resampledPixelSpacing, self.interpolator)
        else:
            image,mask = imageoperations.cropToTumorMask(image, mask)

        return image, mask

    def computeFeatures(self, image, mask, inputImageName, **kwargs):

        featureVector = collections.OrderedDict()
        for featureClassName, enabledFeatures in self.enabledFeatures.iteritems():
            if featureClassName in self.getFeatureClassNames():
                featureClass = self.featureClasses[featureClassName](image, mask, **kwargs)

                if len(enabledFeatures) == 0:
                    featureClass.enableAllFeatures()
                else:
                    for feature in enabledFeatures:
                        featureClass.enableFeatureByName(feature)

                if self.verbose: print "\t\tComputing %s" %(featureClassName)
                featureClass.calculateFeatures()
                for (featureName, featureValue) in featureClass.featureValues.iteritems():
                    newFeatureName = "%s_%s_%s" %(inputImageName, featureClassName, featureName)
                    featureVector[newFeatureName] = featureValue

        return featureVector

    def computeLoG(self, image, mask, **kwargs):

        logFeatureVector = collections.OrderedDict()
        sigmaValues = kwargs.get('sigma', [])

        for sigma in sigmaValues:
            if self.verbose: print "\tComputing LoG with sigma %s" %(str(sigma))
            logImage = imageoperations.applyLoG(image, sigmaValue=sigma)

            inputImageName = "log-sigma-%s-mm-3D" %(str(sigma).replace('.','-'))
            logFeatureVector.update( self.computeFeatures(logImage, mask, inputImageName, **kwargs))

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
                if self.verbose: print "\tComputing Wavelet %s" % (decompositionName)

                if idx == 1:
                    inputImageName = 'wavelet-%s' % (decompositionName)
                else:
                    inputImageName = 'wavelet%s-%s' % (idx, decompositionName)
                waveletFeatureVector.update(self.computeFeatures(decompositionImage, mask, inputImageName, **kwargs))

        if len(ret) == 1:
            inputImageName = 'wavelet-LLL'
        else:
            inputImageName = 'wavelet%s-LLL' % (len(ret))
        waveletFeatureVector.update(self.computeFeatures(approx, mask, inputImageName, **kwargs))

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

    def GetFeaturesFromClass(self, featureClassName):
        """
        Returns a list of all possible features in provided featureClass
        """

        return self.featureClasses[featureClassName].getFeatureNames()
