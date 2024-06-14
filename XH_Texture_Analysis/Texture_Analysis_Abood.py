#!/usr/bin/env python
# %% Import libraries
from __future__ import print_function

import logging
import os

import six

import radiomics
from radiomics import featureextractor, getFeatureClasses

import matplotlib.pyplot as plt
import SimpleITK as sitk
import numpy as np
import nibabel as nib
from scipy.io import loadmat
import pandas as pd
from openpyxl import Workbook

# %% Import ventilation image and mask

# Control patient
print('Control patient')
path = '/home/hoyeh3/GitHub-wsl2/pyradiomics/XH_Texture_Analysis'
path_control = os.path.join(path, 'Healthy_test')

imagefile_control= 'VentilationImage.nii.gz'
maskfile_control = 'MaskImage.nii.gz'

imageName_control = sitk.ReadImage(os.path.join(path_control, imagefile_control))
maskName_control = sitk.ReadImage(os.path.join(path_control, maskfile_control))

# Sanity check
imageVent_control = sitk.GetArrayFromImage(imageName_control)
maskVent_control = sitk.GetArrayFromImage(maskName_control)
for numImage in range(maskVent_control.shape[0]):
  plt.figure(figsize=(10,10))
  plt.subplot(1,3,1)
  plt.imshow(imageVent_control[numImage,:,:], cmap="gray")
  plt.title(f"Ventilation #{numImage}")
  plt.subplot(1,3,2)
  plt.imshow(maskVent_control[numImage,:,:])        
  plt.title(f"Mask #{numImage}")
  plt.subplot(1,3,3)
  plt.imshow(imageVent_control[numImage,:,:]*maskVent_control[numImage,:,:])        
  plt.title(f"Mask and ventilation #{numImage}")

# %% Import the ventilation defects from Matlab
defect_path = os.path.join(path_control, 'VDPThresholdAnalysis.mat')
threshold_mat = loadmat(defect_path)
defect = np.array(threshold_mat['defectArray'])
defect[defect != 0] = 1
mask_defect_control = np.zeros_like(maskVent_control)

# Sanity check to observe the original mask without the ventilation defects
for numImage in range(maskVent_control.shape[0]):
  plt.figure(figsize=(10,10))
  plt.subplot(1,3,1)
  plt.title(f"Defect #{numImage}")
  plt.imshow(defect[:,:,numImage],cmap='gray')
  plt.subplot(1,3,2)
  plt.imshow(maskVent_control[numImage,:,:],cmap='gray')        
  plt.title(f"Mask #{numImage}")
  plt.subplot(1,3,3)
  mask_defect_control[numImage,:,:] = abs(maskVent_control[numImage,:,:]-defect[:,:,numImage])
  plt.imshow(imageVent_control[numImage,:,:]*mask_defect_control[numImage,:,:])        
  plt.title(f"Mask-defect #{numImage}")

# Transform the array from the mask minus the vent. defect into an sitk Image
mask_defect_sitk_control = sitk.GetImageFromArray(mask_defect_control)

# %% Obtain radiomics properties
# Get the location of the settings file
paramsFile = os.path.abspath(os.path.join(path,'Settings', 'Params.yaml'))
if imageName_control is None or maskName_control is None:
  print('Error getting testcase!')
  exit()

# Initialize feature extractor using the settings file 'Params.yaml'
extractor = featureextractor.RadiomicsFeatureExtractor(paramsFile)
featureClasses = getFeatureClasses()

# Print the features description from radiomics
print("Active features:")
for cls, features in six.iteritems(extractor.enabledFeatures):
  if features is None or len(features) == 0:
    features = [f for f, deprecated in six.iteritems(featureClasses[cls].getFeatureNames()) if not deprecated]
  for f in features:
    print(f)
    print(getattr(featureClasses[cls], 'get%sFeatureValue' % f).__doc__)

# %%Print each calculation from the features
print("Calculating features")
featureVector = extractor.execute(imageName_control, mask_defect_sitk_control)

for featureName in featureVector.keys():
  print("Computed feature: %s " % (featureName))

df = pd.DataFrame(featureVector.keys(),columns=['Features'])
df_2 = pd.DataFrame(featureVector.values(),columns=['Measurements'])
df_merged_control = pd.concat([df, df_2], axis=1)
# df_merged.to_csv(path_control+'/healthy_features.csv',index = False)
# df_merged.to_excel(path_control+'/Healthy_features.xlsx', sheet_name='Healthy',index=False)


# %% Import ventilation image and mask
# Lam patient 
print('Lam patient')
path = '/home/hoyeh3/GitHub-wsl2/pyradiomics/XH_Texture_Analysis'
path_Lam = os.path.join(path, 'Lam_test')

imagefile_Lam= 'VentilationImage.nii.gz'
maskfile_Lam = 'MaskImage.nii.gz'

imageName_Lam = sitk.ReadImage(os.path.join(path_Lam, imagefile_Lam))
maskName_Lam = sitk.ReadImage(os.path.join(path_Lam, maskfile_Lam))

# Sanity check
imageVent_Lam = sitk.GetArrayFromImage(imageName_Lam)
maskVent_Lam = sitk.GetArrayFromImage(maskName_Lam)
for numImage in range(maskVent_Lam.shape[0]):
  plt.figure(figsize=(10,10))
  plt.subplot(1,3,1)
  plt.imshow(imageVent_Lam[numImage,:,:], cmap="gray")
  plt.title(f"Ventilation #{numImage}")
  plt.subplot(1,3,2)
  plt.imshow(maskVent_Lam[numImage,:,:])        
  plt.title(f"Mask #{numImage}")
  plt.subplot(1,3,3)
  plt.imshow(imageVent_Lam[numImage,:,:]*maskVent_Lam[numImage,:,:])        
  plt.title(f"Mask and ventilation #{numImage}")

# %% Import the ventilation defects from Matlab
defect_path = os.path.join(path_Lam, 'VDPThresholdAnalysis.mat')
threshold_mat = loadmat(defect_path)
defect = np.array(threshold_mat['defectArray'])
defect[defect != 0] = 1
mask_defect_Lam = np.zeros_like(maskVent_Lam)

# Sanity check to observe the original mask without the ventilation defects
for numImage in range(maskVent_Lam.shape[0]):
  plt.figure(figsize=(10,10))
  plt.subplot(1,3,1)
  plt.title(f"Defect #{numImage}")
  plt.imshow(defect[:,:,numImage],cmap='gray')
  plt.subplot(1,3,2)
  plt.imshow(maskVent_Lam[numImage,:,:],cmap='gray')        
  plt.title(f"Mask #{numImage}")
  plt.subplot(1,3,3)
  mask_defect_Lam[numImage,:,:] = abs(maskVent_Lam[numImage,:,:]-defect[:,:,numImage])
  plt.imshow(imageVent_Lam[numImage,:,:]*mask_defect_Lam[numImage,:,:])        
  plt.title(f"Mask-defect #{numImage}")

# Transform the array from the mask minus the vent. defect into an sitk Image
mask_defect_sitk_Lam = sitk.GetImageFromArray(mask_defect_Lam)

# %% Obtain radiomics properties
# Get the location of the settings file
paramsFile = os.path.abspath(os.path.join(path,'Settings', 'Params.yaml'))
if imageName_Lam is None or maskName_Lam is None:
  print('Error getting testcase!')
  exit()

# Initialize feature extractor using the settings file 'Params.yaml'
extractor = featureextractor.RadiomicsFeatureExtractor(paramsFile)
featureClasses = getFeatureClasses()

# Print each calculation from the features
print("Calculating features")
featureVector = extractor.execute(imageName_Lam, mask_defect_sitk_Lam)

for featureName in featureVector.keys():
  print("Computed feature: %s " % (featureName))

df = pd.DataFrame(featureVector.keys(),columns=['Features'])
df_2 = pd.DataFrame(featureVector.values(),columns=['Measurements'])
df_merged_Lam = pd.concat([df, df_2], axis=1)


# %% Import ventilation image and mask
# Very patchy Lam patient 
print('Lam patient patchy')
path = '/home/hoyeh3/GitHub-wsl2/pyradiomics/XH_Texture_Analysis'
path_Lam_patchy = os.path.join(path, 'Lam_test_patchy')

imagefile_Lam_patchy = 'VentilationImage.nii.gz'
maskfile_Lam_patchy = 'MaskImage.nii.gz'

imageName_Lam_patchy = sitk.ReadImage(os.path.join(path_Lam_patchy, imagefile_Lam_patchy))
maskName_Lam_patchy = sitk.ReadImage(os.path.join(path_Lam_patchy, maskfile_Lam_patchy))

# Sanity check
imageVent_Lam_patchy = sitk.GetArrayFromImage(imageName_Lam_patchy)
maskVent_Lam_patchy = sitk.GetArrayFromImage(maskName_Lam_patchy)
for numImage in range(maskVent_Lam_patchy.shape[0]):
  plt.figure(figsize=(10,10))
  plt.subplot(1,3,1)
  plt.imshow(imageVent_Lam_patchy[numImage,:,:], cmap="gray")
  plt.title(f"Ventilation #{numImage}")
  plt.subplot(1,3,2)
  plt.imshow(maskVent_Lam_patchy[numImage,:,:])        
  plt.title(f"Mask #{numImage}")
  plt.subplot(1,3,3)
  plt.imshow(imageVent_Lam_patchy[numImage,:,:]*maskVent_Lam_patchy[numImage,:,:])        
  plt.title(f"Mask and ventilation #{numImage}")

# %% Import the ventilation defects from Matlab
defect_path = os.path.join(path_Lam_patchy, 'VDPThresholdAnalysis.mat')
threshold_mat = loadmat(defect_path)
defect = np.array(threshold_mat['defectArray'])
defect[defect != 0] = 1
mask_defect_Lam_patchy = np.zeros_like(maskVent_Lam_patchy)

# Sanity check to observe the original mask without the ventilation defects
for numImage in range(maskVent_Lam_patchy.shape[0]):
  plt.figure(figsize=(10,10))
  plt.subplot(1,3,1)
  plt.title(f"Defect #{numImage}")
  plt.imshow(defect[:,:,numImage],cmap='gray')
  plt.subplot(1,3,2)
  plt.imshow(maskVent_Lam_patchy[numImage,:,:],cmap='gray')        
  plt.title(f"Mask #{numImage}")
  plt.subplot(1,3,3)
  mask_defect_Lam_patchy[numImage,:,:] = abs(maskVent_Lam_patchy[numImage,:,:]-defect[:,:,numImage])
  plt.imshow(imageVent_Lam_patchy[numImage,:,:]*mask_defect_Lam_patchy[numImage,:,:])        
  plt.title(f"Mask-defect #{numImage}")

# Transform the array from the mask minus the vent. defect into an sitk Image
mask_defect_sitk_Lam_patchy = sitk.GetImageFromArray(mask_defect_Lam_patchy)

# %% Obtain radiomics properties
# Get the location of the settings file
paramsFile = os.path.abspath(os.path.join(path,'Settings', 'Params.yaml'))
if imageName_Lam_patchy is None or maskName_Lam_patchy is None:
  print('Error getting testcase!')
  exit()

# Initialize feature extractor using the settings file 'Params.yaml'
extractor = featureextractor.RadiomicsFeatureExtractor(paramsFile)
featureClasses = getFeatureClasses()

# Print each calculation from the features
print("Calculating features")
featureVector = extractor.execute(imageName_Lam_patchy, mask_defect_sitk_Lam_patchy)

for featureName in featureVector.keys():
  print("Computed feature: %s " % (featureName))

df = pd.DataFrame(featureVector.keys(),columns=['Features'])
df_2 = pd.DataFrame(featureVector.values(),columns=['Measurements'])
df_merged_Lam_patchy = pd.concat([df, df_2], axis=1)






# %% Store features into an Excel document

#Remove packages versions and original image and mask data, NOT radiomics
df_merged_control_filt = df_merged_control.iloc[22:]
df_merged_Lam_filt = df_merged_Lam.iloc[22:]
df_merged_Lam_filt_patchy = df_merged_Lam_patchy.iloc[22:]

# Save data in excel
with pd.ExcelWriter(path+'/Features.xlsx') as writer:
    df_merged_control.to_excel(writer, sheet_name="Control", index=False)
    df_merged_Lam.to_excel(writer, sheet_name="Lam", index=False)
    df_merged_Lam_patchy.to_excel(writer, sheet_name="Lam", index=False)
    df_merged_control_filt.to_excel(writer, sheet_name="Control_filt", index=False)
    df_merged_Lam_filt.to_excel(writer, sheet_name="Lam_filt", index=False)
    df_merged_Lam_filt_patchy.to_excel(writer, sheet_name="Lam_filt", index=False)

