import numpy
import SimpleITK as sitk

"""
# Legacy Function
def ImageToNumpyArraySlicer (imageNode):
    from __main__ import vtk, qt, ctk, slicer
    # Generate Numpy Array from vtkMRMLScalarVolumeNode
    imageData = vtk.vtkImageData()
    imageData = imageNode.GetImageData()
    shapeData = list(imageData.GetDimensions())
    shapeData.reverse()
    return (vtk.util.numpy_support.vtk_to_numpy(imageData.GetPointData().GetScalars()).reshape(shapeData))
"""

def ImageToNumpyArray(imageFilePath):
    sitkImage = sitk.ReadImage(imageFilePath)
    sitkImageArray = sitk.GetArrayFromImage(sitkImage)
    return sitkImageArray