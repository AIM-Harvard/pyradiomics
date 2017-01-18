from radiomics import imageoperations
import SimpleITK as sitk
import sys

image = sitk.ReadImage(sys.argv[1])
mask = sitk.ReadImage(sys.argv[2])

(ii, im) = imageoperations.resampleImage(image, mask, [2, 2, 2])

sitk.WriteImage(ii, sys.argv[3])
sitk.WriteImage(im, sys.argv[4])
