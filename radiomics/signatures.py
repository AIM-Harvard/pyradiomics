# -*- coding: utf-8 -*-
import os
import collections
import numpy as np
import SimpleITK as sitk
from radiomics import firstorder, glcm, imageoperations, shape, rlgl, glszm


class RadiomicsSignature():

    def __init__(self, **kwargs):
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
            image, mask = imageoperations.resampleImage(image, mask, self.resampledPixelSpacing, self.interpolator)

        return image, mask

    def computeFeatures(self, image, mask, **kwargs):

        FeatureVector = collections.OrderedDict()

        if self.verbose: print "\tComputing First Order"
        FeatureVector.update( self.computeFirstOrder(image, mask, **kwargs) )

        if self.verbose: print "\tComputing Shape"
        FeatureVector.update( self.computeShape(image, mask, **kwargs) )

        if self.verbose: print "\tComputing GLCM"
        FeatureVector.update( self.computeGLCM(image, mask, **kwargs) )

        if self.verbose: print "\tComputing RLGL"
        FeatureVector.update( self.computeRLGL(image, mask, **kwargs) )

        if self.verbose: print "\tComputing GLSZM"
        FeatureVector.update( self.computeGLSZM(image, mask, **kwargs) )

        return FeatureVector

    def computeFirstOrder(self, image, mask, **kwargs):
        firstOrderFeatures = firstorder.RadiomicsFirstOrder(image, mask, **kwargs)
        if 'firstorder' in self.enabledFeatures:
            for enabledFeature in self.enabledFeatures['firstorder']:
                firstOrderFeatures.enableFeatureByName(enabledFeature)
        else:
            firstOrderFeatures.enableAllFeatures()
        firstOrderFeatures.calculateFeatures()

        firstOrderFeatureVector = collections.OrderedDict()
        for (featureName, featureValue) in firstOrderFeatures.featureValues.iteritems():
            firstOrderFeatureName = "firstorder_%s" %(featureName)
            firstOrderFeatureVector[firstOrderFeatureName] = featureValue

        return firstOrderFeatureVector

    def computeShape(self, image, mask, **kwargs):
        shapeFeatures = shape.RadiomicsShape(image, mask, **kwargs)
        if 'shape' in self.enabledFeatures:
            for enabledFeature in self.enabledFeatures['shape']:
                shapeFeatures.enableFeatureByName(enabledFeature)
        else:
            shapeFeatures.enableAllFeatures()
        shapeFeatures.calculateFeatures()

        shapeFeatureVector = collections.OrderedDict()
        for (featureName, featureValue) in shapeFeatures.featureValues.iteritems():
            shapeFeatureName = "shape_%s" %(featureName)
            shapeFeatureVector[shapeFeatureName] = featureValue

        return shapeFeatureVector

    def computeGLCM(self, image, mask, **kwargs):
        glcmFeatures = glcm.RadiomicsGLCM(image, mask, **kwargs)
        if 'glcm' in self.enabledFeatures:
            for enabledFeature in self.enabledFeatures['glcm']:
                glcmFeatures.enableFeatureByName(enabledFeature)
        else:
            glcmFeatures.enableAllFeatures()
        glcmFeatures.calculateFeatures()

        glcmFeatureVector = collections.OrderedDict()
        for (featureName, featureValue) in glcmFeatures.featureValues.iteritems():
            glcmFeatureName = "glcm_%s" %(featureName)
            glcmFeatureVector[glcmFeatureName] = featureValue

        return glcmFeatureVector

    def computeRLGL(self, image, mask, **kwargs):
        rlglFeatures = rlgl.RadiomicsRLGL(image, mask, **kwargs)
        if 'rlgl' in self.enabledFeatures:
            for enabledFeature in self.enabledFeatures['rlgl']:
                rlglFeatures.enableFeatureByName(enabledFeature)
        else:
            rlglFeatures.enableAllFeatures()
        rlglFeatures.calculateFeatures()

        rlglFeatureVector = collections.OrderedDict()
        for (featureName, featureValue) in rlglFeatures.featureValues.iteritems():
            rlglFeatureName = "rlgl_%s" %(featureName)
            rlglFeatureVector[rlglFeatureName] = featureValue

        return rlglFeatureVector

    def computeGLSZM(self, image, mask, **kwargs):
        glszmFeatures = glszm.RadiomicsGLSZM(image, mask, **kwargs)
        if 'glszm' in self.enabledFeatures:
            for enabledFeature in self.enabledFeatures['glszm']:
                glszmFeatures.enableFeatureByName(enabledFeature)
        else:
            glszmFeatures.enableAllFeatures()
        glszmFeatures.calculateFeatures()

        glszmFeatureVector = collections.OrderedDict()
        for (featureName, featureValue) in glszmFeatures.featureValues.iteritems():
            glszmFeatureName = "glszm_%s" %(featureName)
            glszmFeatureVector[glszmFeatureName] = featureValue

        return glszmFeatureVector

    def computeLoG(self, image, mask, **kwargs):

        logFeatureVector = collections.OrderedDict()
        sigmaValues = kwargs.get('sigma', [])
        if 'sigma' in kwargs: del kwargs['sigma']

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

        if 'wavelet' in kwargs: del kwargs['wavelet']
        if 'level' in kwargs: del kwargs['level']
        if 'start_level' in kwargs: del kwargs['start_level']

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
