"""
This script is an example of how a mask image can be resampled to the geometry of a reference image. This enables the
a mask that has been generated on one image, to be used for feature extraction of another image. Please note that this
still requires the mask to have a similar physical space as the reference image.

This script can be used directly from the command line by calling ``python resampleMask.py <args>``, where <args> are
the reference image, mask to resample and a filename to store the result (the resampled mask).
"""
import argparse

import SimpleITK as sitk

parser = argparse.ArgumentParser()
parser.add_argument('image', metavar='Image', help='Reference image to resample the mask to')
parser.add_argument('mask', metavar='Mask', help='Input mask to resample')
parser.add_argument('resMask', metavar='Out', help='Filename to store resampled mask')

def main():
  args = parser.parse_args()
  image = sitk.ReadImage(args.image)
  mask = sitk.ReadImage(args.mask)

  rif = sitk.ResampleImageFilter()
  rif.SetReferenceImage(image)
  rif.SetOutputPixelType(mask.GetPixelID())
  rif.SetInterpolator(sitk.sitkNearestNeighbor)
  resMask = rif.Execute(mask)

  sitk.WriteImage(resMask, args.resMask, True)  # True enables compression when saving the resampled mask

if __name__ == '__main__':
    main()
