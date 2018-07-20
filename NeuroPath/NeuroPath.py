import os
import unittest
import vtk, qt, ctk, slicer
from slicer.ScriptedLoadableModule import *
import logging
import numpy as np
import numpy.linalg as la
import math
import datetime
import csv
global renderer
renderer = slicer.app.layoutManager().threeDWidget(0).threeDView().renderWindow().GetRenderers().GetFirstRenderer()
renderWindow = slicer.app.layoutManager().threeDWidget(0).threeDView().renderWindow()
renderWindow.AddRenderer(renderer)
renderWindowInteractor = slicer.app.layoutManager().threeDWidget(0).threeDView().interactor()
renderWindowInteractor.SetRenderWindow(renderWindow)
 
#
# NeuroPath
#
 
class NeuroPath(ScriptedLoadableModule):
  """Uses ScriptedLoadableModule base class, available at:
  https://github.com/Slicer/Slicer/blob/master/Base/Python/slicer/ScriptedLoadableModule.py
  """
 
  def __init__(self, parent):
    ScriptedLoadableModule.__init__(self, parent)
    self.parent.title = "NeuroPath" # TODO make this more human readable by adding spaces
    self.parent.categories = ["Examples"]
    self.parent.dependencies = []
    self.parent.contributors = ["John Doe (AnyWare Corp.)"] # replace with "Firstname Lastname (Organization)"
    self.parent.helpText = """
    This is an example of scripted loadable module bundled in an extension.
    It performs a simple thresholding on the input volume and optionally captures a screenshot.
    """
    self.parent.acknowledgementText = """
    This file was originally developed by Jean-Christophe Fillion-Robin, Kitware Inc.
    and Steve Pieper, Isomics, Inc. and was partially funded by NIH grant 3P41RR013218-12S1.
""" # replace with organization, grant and thanks.
 
#
# NeuroPathWidget
#
 
class NeuroPathWidget(ScriptedLoadableModuleWidget):
  """Uses ScriptedLoadableModuleWidget base class, available at:
  https://github.com/Slicer/Slicer/blob/master/Base/Python/slicer/ScriptedLoadableModule.py
  """
 
  def setup(self):
    ScriptedLoadableModuleWidget.setup(self)
 
    # Instantiate and connect widgets ...
 
    #
    # Parameters Area
    #
    parametersCollapsibleButton = ctk.ctkCollapsibleButton()
    parametersCollapsibleButton.text = "Parameters"
    self.layout.addWidget(parametersCollapsibleButton)
 
    # Layout within the dummy collapsible button
    parametersFormLayout = qt.QFormLayout(parametersCollapsibleButton)
    #
    # Fiducial  Selector
    #
    self.ROISelector = slicer.qMRMLNodeComboBox()
    self.ROISelector.nodeTypes = ["vtkMRMLMarkupsFiducialNode"]
    self.ROISelector.selectNodeUponCreation = True
    self.ROISelector.addEnabled = True
    self.ROISelector.removeEnabled = True
    self.ROISelector.noneEnabled = True
    self.ROISelector.showHidden = False
    self.ROISelector.showChildNodeTypes = False
    self.ROISelector.setMRMLScene( slicer.mrmlScene )
    self.ROISelector.setToolTip( "Pick List of Markups Start & Targets" )
    parametersFormLayout.addRow("Fiducial : ", self.ROISelector)
    self.ROISelector.connect("currentNodeChanged(vtkMRMLNode*)", self.OnSelect_ROI)

    #
    # Sphere Radius
    #
    self.SliderSphere = ctk.ctkSliderWidget()
    self.SliderSphere.singleStep = 0.1
    self.SliderSphere.minimum = 0
    self.SliderSphere.maximum = 20
    self.SliderSphere.value = 10
    self.SliderSphere.setToolTip("Sphere")
    self.SliderSphere.valueChanged.connect(self.onSliderChange)
    parametersFormLayout.addRow("Sphere Radius", self.SliderSphere)
    
    #
    # Cylinder Radius
    #
    self.SliderCylinder = ctk.ctkSliderWidget()
    self.SliderCylinder.singleStep = 0.1
    self.SliderCylinder.minimum = 0
    self.SliderCylinder.maximum = 20
    self.SliderCylinder.value = 7.5
    self.SliderCylinder.setToolTip("Set threshold value for computing the output image. Voxels that have intensities lower than this value will set to zero.")
    self.SliderCylinder.valueChanged.connect(self.onSliderChange)
    parametersFormLayout.addRow("Cylinder Radius", self.SliderCylinder)

    #
    # Cylinder Height
    #
    self.SliderCylinder_2 = ctk.ctkSliderWidget()
    self.SliderCylinder_2.singleStep = 0.1
    self.SliderCylinder_2.minimum = 5
    self.SliderCylinder_2.maximum = 120
    self.SliderCylinder_2.value = 30
    self.SliderCylinder_2.setToolTip("Set threshold value for computing the output image. Voxels that have intensities lower than this value will set to zero.")
    self.SliderCylinder_2.valueChanged.connect(self.onSliderChange)
    parametersFormLayout.addRow("Cylinder Height", self.SliderCylinder_2)

    #
    # ROI Selector get the box around the brain
    #
    self.ROISelector_2 = slicer.qMRMLNodeComboBox()
    self.ROISelector_2.nodeTypes = ["vtkMRMLAnnotationROINode"]
    self.ROISelector_2.selectNodeUponCreation = True
    self.ROISelector_2.addEnabled = True
    self.ROISelector_2.removeEnabled = True
    self.ROISelector_2.noneEnabled = True
    self.ROISelector_2.showHidden = False
    self.ROISelector_2.showChildNodeTypes = False
    self.ROISelector_2.setMRMLScene( slicer.mrmlScene )
    self.ROISelector_2.setToolTip( "Pick the ROI Interactive Box to define the region for labeling of Fibers." )
    parametersFormLayout.addRow("ROI for Labeling: ", self.ROISelector_2)
    self.ROISelector_2.connect("currentNodeChanged(vtkMRMLNode*)", self.OnSelect_OptBox)
  

  def onSliderChange(self):
      global R_sphere,R_cylinder,H_cylinder
      R_sphere =  self.SliderSphere.value
      R_cylinder = self.SliderCylinder.value
      H_cylinder = self.SliderCylinder_2.value


  ## Uses Fiducial Selector, Sphere, Cylinder Radius & Cylynder Height
  def OnSelect_ROI(self):

  ## Supposed to get the sphere to be at fiudcial 3 and the cylinder between fiducial 1 and fiducial 2
      global ft12,markupNode,sphere,cylinder,logic, plane,model_C, model_S, transF
      logic = NeuroPathLogic()
      fiducial1 = [0,0,0] 
      fiducial2 = [0,0,0]
      fiducial3 = [0,0,0]
      #gets the position in RAS of 3 fiducials
      markupNode = self.ROISelector.currentNode()
      markupNode.GetNthFiducialPosition(0,fiducial1)
      markupNode.GetNthFiducialPosition(1,fiducial2)
      markupNode.GetNthFiducialPosition(2,fiducial3)
      #print ftone 
      #print fttwo
      #print fss 
     
      ## not sure this is necessary since its already defined at onSliderChange
      R_sphere =  self.SliderSphere.value
      R_cylinder = self.SliderCylinder.value
      H_cylinder = self.SliderCylinder_2.value
      
      #----------------------------
      #creates an array
      f1 = np.array(fiducial1)
      f2 = np.array(fiducial2)
      f3 = np.array(fiducial3)
      deltaf1f2 = f1-f2  #position of fiducial 1 - fiducial 2. Why??
      Y_axis = np.array([0,1,0])
      C_Cent = ((f1+f2)/2) 
     
      #------------------------------
      sphere = vtk.vtkSphereSource()
      sphere.SetCenter(f3) # Set the center of the sphere.
      sphere.SetThetaResolution(30) #set the number of points in the longitude direction
      sphere.SetPhiResolution(30) #set the number of points in the latitude direction
      sphere.SetRadius(R_sphere) #set radius of sphere

      #
      Mapper = vtk.vtkPolyDataMapper()
      Mapper.SetInputConnection(sphere.GetOutputPort())


      #represents an object in a rendering scence -> gets properties from mapper
      Actor = vtk.vtkActor()
      Actor.SetMapper(Mapper)
      Actor.GetProperty().SetColor(0,1,0)

      #clips polygonal data returning everything inside teh specified implicit function
      clip = vtk.vtkClipPolyData()
      clip.SetValue(0);
      clip.GenerateClippedOutputOn();
      clip.SetInputConnection(sphere.GetOutputPort())

      plane = vtk.vtkPlane() #plane computations
      normal = deltaf1f2/np.sqrt(sum(deltaf1f2*deltaf1f2))
      fend = C_Cent+(H_cylinder/2)*normal
      plane.SetNormal(normal)
      plane.SetOrigin(fend)
      clip.SetClipFunction(plane)

      cylinder = vtk.vtkCylinderSource()
      cylinder.SetRadius(R_cylinder)
      cylinder.SetHeight(H_cylinder)
      cylinder.SetCenter(C_Cent)
      cylinder.SetResolution(100)
      cylinder.CappingOff()

      Transform = vtk.vtkTransform()
      Transform.PostMultiply()
      Transform.Translate(-C_Cent[0],-C_Cent[1],-C_Cent[2])
      #Transform.RotateWXYZ(-1*np.degrees(logic.py_ang(deltaf1f2,Y_axis)),np.cross(deltaf1f2,Y_axis))
      Transform.Translate(C_Cent[0],C_Cent[1],C_Cent[2])
      transF = vtk.vtkTransformPolyDataFilter()
      transF.SetInputConnection(cylinder.GetOutputPort())
      transF.SetTransform(Transform)
      transF.Update()

      modelsLogic = slicer.modules.models.logic()
      model_C = modelsLogic.AddModel(transF.GetOutputPort())
      model_C.GetDisplayNode().SetSliceIntersectionVisibility(True)
      model_C.GetDisplayNode().SetSliceIntersectionThickness(2)
      model_C.GetDisplayNode().SetColor(0,0,1)
      model_C.SetName('Cylinder')

      model_S = modelsLogic.AddModel(clip.GetOutputPort())
      model_S.GetDisplayNode().SetSliceIntersectionVisibility(True)
      model_S.GetDisplayNode().SetSliceIntersectionThickness(2)
      model_S.GetDisplayNode().SetColor(0,1,0)
      model_S.SetName('Sphere')

  ### holds the ROI coordinates found in Annotation Properities
  def OnSelect_OptBox(self):
    global bounds,markupNode_2
    bounds = [0,0,0,0,0,0] # [L,R,P,A,I,S]
    markupNode_2 = self.ROISelector_2.currentNode()
    markupNode_2.GetRASBounds(bounds)
    print bounds


# NeuroPathLogic
#
 
class NeuroPathLogic(ScriptedLoadableModuleLogic):
  """This class should implement all the actual
  computation done by your module.  The interface
  should be such that other python code can import
  this class and make use of the functionality without
  requiring an instance of the Widget.
  Uses ScriptedLoadableModuleLogic base class, available at:
  https://github.com/Slicer/Slicer/blob/master/Base/Python/slicer/ScriptedLoadableModule.py
  """
  global Score, newArrays
  def run_3(self):
    X_1 = int(bounds[0]/10)  #takes the first item in the bounds array and divides by 10
    X_2 = int(bounds[1]/10)
    Y_1 = int(bounds[2]/10)
    Y_2 = int(bounds[3]/10)
    Z_1 = int(bounds[4]/10)
    Z_2 = int(bounds[5]/10)
    Num_sample = (X_2-X_1+1)*(Y_2-Y_1+1)+2*((Y_2-Y_1+1)*(Z_2-Z_1+1)+(X_2-X_1+1)*(Z_2-Z_1+1))
    Lg_Table = np.zeros((Num_sample+1,10)) # 10 columns 
    Lg_Table[0,4] = X_1
    Lg_Table[0,5] = X_2
    Lg_Table[0,6] = Y_1
    Lg_Table[0,7] = Y_2
    Lg_Table[0,8] = Z_1
    Lg_Table[0,9] = Z_2
    k = 1








