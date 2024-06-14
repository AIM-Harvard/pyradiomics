#!/usr/bin/env python
# %% Import libraries
from __future__ import print_function

import logging
import os

import six

import radiomics
from radiomics import featureextractor, getFeatureClasses, shape2D
import matplotlib.pyplot as plt
import SimpleITK as sitk
import nibabel as nib
import numpy as np

# %% Start using radiomics
# Get some test data
path_control = '/home/hoyeh3/GitHub-wsl2/pyradiomics/XH_Texture_Analysis/Lam_test'
imagefile_control= 'VentilationImage.nii.gz'
maskfile_control = 'MaskImage.nii.gz'
imageName = sitk.ReadImage(os.path.join(path_control, imagefile_control))
maskName = sitk.ReadImage(os.path.join(path_control, maskfile_control))

# Get the location of the example settings file
paramsFile = os.path.abspath(os.path.join('Settings', 'Params.yaml'))

if imageName is None or maskName is None:  # Something went wrong, in this case PyRadiomics will also log an error
  print('Error getting testcase!')
  exit()

# Initialize feature extractor using the settings file
extractor = featureextractor.RadiomicsFeatureExtractor(paramsFile)
featureClasses = getFeatureClasses()

print("Active features:")
for cls, features in six.iteritems(extractor.enabledFeatures):
  if features is None or len(features) == 0:
    features = [f for f, deprecated in six.iteritems(featureClasses[cls].getFeatureNames()) if not deprecated]
  for f in features:
    print(f)
    print(getattr(featureClasses[cls], 'get%sFeatureValue' % f).__doc__)

print("Calculating features")
featureVector = extractor.execute(imageName, maskName)

for featureName in featureVector.keys():
  print("Computed %s: %s" % (featureName, featureVector[featureName]))

imageVent_control = sitk.GetArrayFromImage(imageName)
maskVent_control = sitk.GetArrayFromImage(maskName)
for numImage in range(imageVent_control.shape[0]):
  plt.figure(figsize=(10,10))
  plt.subplot(2,2,1)
  plt.imshow(imageVent_control[numImage,:,:], cmap="gray")
  plt.title(f"Ventilation #{numImage}")
  plt.subplot(2,2,2)
  plt.imshow(maskVent_control[numImage,:,:])        
  plt.title(f"Mask #{numImage}")



# %% Use nibabel instead of SimpleItk

from radiomics import shape2D


path = '/home/hoyeh3/GitHub-wsl2/pyradiomics/XH_Texture_Analysis/Lam_test'
imagefile= 'VentilationImage.nii.gz'
maskfile = 'MaskImage.nii.gz'
imageName = nib.load(os.path.join(path, imagefile)).get_fdata()
maskName = nib.load(os.path.join(path, maskfile)).get_fdata()

# Create a folder to store all the *.nii files
xh_folder = os.path.join(path, r'xh_analysis')
try:
    os.mkdir(xh_folder)
    print("Folder %s created!" % xh_folder)
except FileExistsError:
    print("Folder %s already exists" % xh_folder)

# Separate the slices and the mask into independent *.nii files
for i in range(imageName.shape[2]):
    output = np.rot90(imageName[:,:,i],2)
    ni_img = nib.Nifti1Image(output, affine=np.eye(4))
    nib.save(ni_img, xh_folder+'/'+'VentilationImage_{0}.nii.gz'.format(i))

    # Rotate the images as nibabel pulls the images upside-down and stores them in the folder created
    output = np.rot90(maskName[:,:,i],2)
    ni_img = nib.Nifti1Image(output, affine=np.eye(4))
    nib.save(ni_img, xh_folder+'/'+'MaskImage_{0}.nii.gz'.format(i))

    print('VentilationImage_{0}.nii.gz and MaskImage_{0}.nii.gz'.format(i))

# To obtain the features slice by slice, 2D. To do so force2D=True
for numImage in range(3,16): #range(imageName.shape[2]):
    path_image = os.path.join(xh_folder,'VentilationImage_{0}.nii.gz').format(numImage) 
    path_mask = os.path.join(xh_folder,'MaskImage_{0}.nii.gz').format(numImage)

    imageName = sitk.ReadImage(path_image)
    maskName = sitk.ReadImage(path_mask)

    plt.figure(figsize=(10,10))
    plt.subplot(2,2,1)
    plt.imshow(np.rot90(sitk.GetArrayFromImage(imageName),2), cmap="gray")
    plt.title(f"Ventilation #{numImage}")
    plt.subplot(2,2,2)
    plt.imshow(np.rot90(sitk.GetArrayFromImage(maskName),2))        
    plt.title(f"Mask #{numImage}")

    paramsFile = os.path.abspath(os.path.join('Settings', 'Params.yaml'))
    extractor = featureextractor.RadiomicsFeatureExtractor(paramsFile)
    featureClasses = getFeatureClasses()

    print(f"Calculating features for ventilation slice: {numImage}")
    featureVector = extractor.execute(imageName, maskName)

    for featureName in featureVector.keys():
        print("Computed %s: %s" % (featureName, featureVector[featureName]))
    print(' ')



# %% Read the images one by one and obtain different properties, 
# extractor was 'shape2D.RadiomicsShape2D'

path = '/home/hoyeh3/GitHub-wsl2/pyradiomics/XH_Texture_Analysis/Lam_test'
imagefile= 'VentilationImage.nii.gz'
maskfile = 'MaskImage.nii.gz'
imageName = nib.load(os.path.join(path, imagefile)).get_fdata()
maskName = nib.load(os.path.join(path, maskfile)).get_fdata()

# Create a folder to store all the *.nii files
xh_folder = os.path.join(path, r'xh_analysis')
try:
    os.mkdir(xh_folder)
    print("Folder %s created!" % xh_folder)
except FileExistsError:
    print("Folder %s already exists" % xh_folder)

# Separate the slices and the mask into independent *.nii files
for i in range(imageName.shape[2]):
    output = np.rot90(imageName[:,:,i],2)
    ni_img = nib.Nifti1Image(output, affine=np.eye(4))
    nib.save(ni_img, xh_folder+'/'+'VentilationImage_{0}.nii.gz'.format(i))

    # Rotate the images as nibabel pulls the images upside-down and stores them in the folder created
    output = np.rot90(maskName[:,:,i],2)
    ni_img = nib.Nifti1Image(output, affine=np.eye(4))
    nib.save(ni_img, xh_folder+'/'+'MaskImage_{0}.nii.gz'.format(i))

    print('VentilationImage_{0}.nii.gz and MaskImage_{0}.nii.gz'.format(i))

print(' ')
# To obtain the features slice by slice, 2D.
for numImage in range(imageName.shape[2]):
  path_image = os.path.join(xh_folder,'VentilationImage_{0}.nii.gz').format(numImage) 
  path_mask = os.path.join(xh_folder,'MaskImage_{0}.nii.gz').format(numImage)

  imageName = sitk.ReadImage(path_image)
  maskName = sitk.ReadImage(path_mask)

  extractor = shape2D.RadiomicsShape2D(imageName,maskName)
  featureClasses = getFeatureClasses()

  print(f"Calculating features for ventilation slice: {numImage}")
  featureVector = extractor.execute()

  for featureName in featureVector.keys():
    print("Computed %s: %s" % (featureName, featureVector[featureName]))
  print(' ')

  plt.figure(figsize=(10,10))
  plt.subplot(2,2,1)
  plt.imshow(np.rot90(sitk.GetArrayFromImage(imageName),2), cmap="gray")
  plt.title(f"Ventilation #{numImage}")
  plt.subplot(2,2,2)
  plt.imshow(np.rot90(sitk.GetArrayFromImage(maskName),2))        
  plt.title(f"Mask #{numImage}")

