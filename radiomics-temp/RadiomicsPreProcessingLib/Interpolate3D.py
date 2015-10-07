from __future__ import print_function

def interpolateScalarVolumeNode3D(selNode, outputPixelSpacing=(3,3,3)):
  from __main__ import vtk, qt, ctk, slicer, os
  resampledNode = slicer.vtkMRMLScalarVolumeNode()
  resampledImageNode.SetName(selNode.GetName() + "_resampled")
  slicer.mrmlScene.AddNode(resampledNode)
  resampleScalarVolume(selNode,resampledNode,outputPixelSpacing)
  return resampledNode

def resampleScalarVolume(node,resampledNode,outputPixelSpacing):
  parameters = {}
  parameters["outputPixelSpacing"] = outputPixelSpacing
  parameters["interpolationType"] = "linear"
  parameters["InputVolume"] = node.GetID() 
  parameters["OutputVolume"] = resampledNode.GetID()
  resampleVolume = slicer.modules.resamplescalarvolume
  return (slicer.cli.run(resampleVolume, None, parameters, wait_for_completion = True))