from radiomics import firstorder, preprocessing
import SimpleITK as sitk
import sys

imageName = '/Users/fedorov/github/radiomics-platform/Testing/Data/prostate_phantom.nrrd'
maskName = '/Users/fedorov/github/radiomics-platform/Testing/Data/prostate_phantom_label.nrrd'

#imageName = sys.argv[1]
#maskName = sys.argv[2]

image = sitk.ReadImage(imageName)
mask = sitk.ReadImage(maskName)

imageArray = sitk.GetArrayFromImage(image)
maskArray = sitk.GetArrayFromImage(mask)

(matrix, matrixCoordinates) = preprocessing.RadiomicsHelpers.padTumorMaskToCube(imageArray,maskArray)

f = firstorder.RadiomicsFirstOrder(matrix, matrixCoordinates, 10, (1,1,1) )
print f.getFeatureNames()

#f.getAllFeatureNames()
