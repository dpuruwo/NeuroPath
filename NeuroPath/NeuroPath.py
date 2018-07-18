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
    self.ROISelector.setToolTip( "Pick the ROI Interactive Box to define the region for labeling of Fibers." )
    parametersFormLayout.addRow("Fiducial : ", self.ROISelector)
    self.ROISelector.connect("currentNodeChanged(vtkMRMLNode*)", self.OnSelect_ROI)
 
    #
    # Sphere value
    #
    self.SliderSphere = ctk.ctkSliderWidget()
    self.SliderSphere.singleStep = 0.1
    self.SliderSphere.minimum = 0
    self.SliderSphere.maximum = 20
    self.SliderSphere.value = 10
    self.SliderSphere.setToolTip("Set threshold value for computing the output image. Voxels that have intensities lower than this value will set to zero.")
    self.SliderSphere.valueChanged.connect(self.onSliderChange)
    parametersFormLayout.addRow("Sphere Radius", self.SliderSphere)
     
    #
    # Cylinder value
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
    # Cylinder value
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
    # Apply Button
    #
    self.applyButton_7 = qt.QCheckBox("Update the path")
    self.applyButton_7.toolTip = "Run the algorithm!!."
    #self.applyButton_5.enabled = True
    parametersFormLayout.addRow("Path Refresh", self.applyButton_7)
    # connections
    self.applyButton_7.stateChanged.connect(self.onApplyButton_7)
	
    #
    # Apply Button
    #
    self.applyButton_8 = qt.QCheckBox("Update the fibers")
    self.applyButton_8.toolTip = "Run the algorithm!!."
    #self.applyButton_5.enabled = True
    parametersFormLayout.addRow("Fiber and Grey Matter Refresh", self.applyButton_8)
    # connections
    self.applyButton_8.stateChanged.connect(self.onApplyButton_8)
	
	#
    # Apply Button
    #
    self.applyButton_9 = qt.QCheckBox("Record the scores")
    self.applyButton_9.toolTip = "Run the algorithm!!."
    #self.applyButton_5.enabled = True
    parametersFormLayout.addRow("Calculate and save scores", self.applyButton_9)
    # connections
    self.applyButton_9.stateChanged.connect(self.onApplyButton_9)
    #
    # Input the Labeled Volume
    #
    self.Label = slicer.qMRMLNodeComboBox()
    self.Label.nodeTypes = ["vtkMRMLScalarVolumeNode"]
    self.Label.selectNodeUponCreation = True
    self.Label.addEnabled = False
    self.Label.removeEnabled = False
    self.Label.noneEnabled = True
    self.Label.showHidden = False
    self.Label.showChildNodeTypes = False
    self.Label.setMRMLScene( slicer.mrmlScene )
    self.Label.setToolTip( "Pick the Atlas-based Segmented Brain Volume." )
    parametersFormLayout.addRow("Input the Atlas Labeled Volume: ", self.Label)
    self.Label.connect("currentNodeChanged(vtkMRMLNode*)", self.OnSelect_label) # missing self.OnSelect_label function
	#
    # Generate a rectangular Volume
    #
    self.Rec_Buck = slicer.qMRMLNodeComboBox()
    self.Rec_Buck.nodeTypes = ["vtkMRMLScalarVolumeNode"]
    self.Rec_Buck.selectNodeUponCreation = True
    self.Rec_Buck.addEnabled = True
    self.Rec_Buck.removeEnabled = True
    self.Rec_Buck.noneEnabled = True
    self.Rec_Buck.showHidden = False
    self.Rec_Buck.showChildNodeTypes = False
    self.Rec_Buck.setMRMLScene( slicer.mrmlScene )
    self.Rec_Buck.setToolTip( "Pick the Atlas-based Segmented Brain Volume." )
    #parametersFormLayout.addRow("Generate the labeled Volume: ", self.Rec_Buck)
    self.Rec_Buck.connect("currentNodeChanged(vtkMRMLNode*)", self.OnSelect_Rec_Buck) # missing self.OnSelect_label function
	
    #
    # Generate a rectangular Box for Score Guidance
    #
    self.Rec_Box = slicer.qMRMLNodeComboBox()
    self.Rec_Box.nodeTypes = ["vtkMRMLScalarVolumeNode"]
    self.Rec_Box.selectNodeUponCreation = True
    self.Rec_Box.addEnabled = True
    self.Rec_Box.removeEnabled = True
    self.Rec_Box.noneEnabled = True
    self.Rec_Box.showHidden = False
    self.Rec_Box.showChildNodeTypes = False
    self.Rec_Box.setMRMLScene( slicer.mrmlScene )
    self.Rec_Box.setToolTip( "Pick the Atlas-based Segmented Brain Volume." )
    parametersFormLayout.addRow("Generate the Score Box: ", self.Rec_Box)
    self.Rec_Box.connect("currentNodeChanged(vtkMRMLNode*)", self.OnSelect_Rec_Box) # missing self.OnSelect_label function
	
    #
    # Import the rectangular volume
    #
    self.Rec_Buck_2 = slicer.qMRMLNodeComboBox()
    self.Rec_Buck_2.nodeTypes = ["vtkMRMLScalarVolumeNode"]
    self.Rec_Buck_2.selectNodeUponCreation = False
    self.Rec_Buck_2.addEnabled = False
    self.Rec_Buck_2.removeEnabled = False
    self.Rec_Buck_2.noneEnabled = True
    self.Rec_Buck_2.showHidden = False
    self.Rec_Buck_2.showChildNodeTypes = False
    self.Rec_Buck_2.setMRMLScene( slicer.mrmlScene )
    self.Rec_Buck_2.setToolTip( "Pick the Atlas-based Segmented Brain Volume." )
    parametersFormLayout.addRow("Input the Skull Mask Volume: ", self.Rec_Buck_2)
    self.Rec_Buck_2.connect("currentNodeChanged(vtkMRMLNode*)", self.OnSelect_Rec_Buck_2) 
   
    #
	#Input tractography file selector
	#
    self.inputSelector = slicer.qMRMLNodeComboBox()
    self.inputSelector.nodeTypes = ["vtkMRMLFiberBundleNode"]
    self.inputSelector.selectNodeUponCreation = True	
    self.inputSelector.addEnabled = False
    self.inputSelector.removeEnabled = False
    self.inputSelector.noneEnabled = True
    self.inputSelector.showHidden = False
    self.inputSelector.showChildNodeTypes = False
    self.inputSelector.setMRMLScene( slicer.mrmlScene )
    self.inputSelector.setToolTip( "Pick the Full-brain Tractography file." )
    parametersFormLayout.addRow("Input Fiber Bundle: ", self.inputSelector)
    self.inputSelector.connect("currentNodeChanged(vtkMRMLNode*)", self.OnSelect_fb) #fb = self.inputSelector.currentNode() add this to applybutton
	
    #
    # Apply Button
    #
    self.applyButton_4 = qt.QPushButton("Generate and Save Arrays")
    self.applyButton_4.toolTip = "Run the algorithm!!."
    #self.applyButton_4.enabled = True
    #FunctionalFormLayout.addWidget(self.applyButton_4,1,0,1,1)
    parametersFormLayout.addRow(self.applyButton_4)
    # connections
    self.applyButton_4.connect('clicked(bool)', self.onApplyButton_4)
	
	#
    # Apply Button
    #
    self.applyButton_5 = qt.QPushButton("Load Saved Array")
    self.applyButton_5.toolTip = "Run the algorithm!!."
    #self.applyButton_5.enabled = True
    #FunctionalFormLayout.addWidget(self.applyButton_5,1,1,1,1)
    parametersFormLayout.addRow(self.applyButton_5)
    # connections
    self.applyButton_5.connect('clicked(bool)', self.onApplyButton_5)

    
	#
    # output fiberbundle generator
    #
    self.outputSelector = slicer.qMRMLNodeComboBox()
    self.outputSelector.nodeTypes = ["vtkMRMLFiberBundleNode"]
    self.outputSelector.selectNodeUponCreation = True
    self.outputSelector.addEnabled = True
    self.outputSelector.removeEnabled = True
    self.outputSelector.noneEnabled = True
    self.outputSelector.showHidden = False
    self.outputSelector.showChildNodeTypes = False
    self.outputSelector.setMRMLScene( slicer.mrmlScene )
    self.outputSelector.setToolTip( "Create or pick the output FiberBundle file." )
    parametersFormLayout.addRow("Output Fiberbundle: ", self.outputSelector)
    self.outputSelector.connect("currentNodeChanged(vtkMRMLNode*)", self.OnSelect_selected_fibers)
	
    #
    # Input the Models Hierarchy Node
    #
    self.Models = slicer.qMRMLNodeComboBox()
    self.Models.nodeTypes = ["vtkMRMLModelHierarchyNode"]
    self.Models.selectNodeUponCreation = True
    self.Models.addEnabled = False
    self.Models.removeEnabled = False
    self.Models.noneEnabled = True
    self.Models.showHidden = False
    self.Models.showChildNodeTypes = False
    self.Models.setMRMLScene( slicer.mrmlScene )
    self.Models.setToolTip( "Pick the Hierarchy of Brain Surface Models." )
    parametersFormLayout.addRow("Input the Generated Models: ", self.Models)
    self.Models.connect("currentNodeChanged(vtkMRMLNode*)", self.OnSelect_Mod)
    
	#
    # Apply Button
    #
    self.RadioLayout = qt.QHBoxLayout()
    Button_grp = qt.QButtonGroup()
    RadioButtons=[]
    Button_grp = qt.QButtonGroup()
    RadioButtons.append(qt.QRadioButton('FC_GM')) 
    RadioButtons.append(qt.QRadioButton('SC_GM'))
    RadioButtons.append(qt.QRadioButton('Opacity_reset')) 
    Button_grp.addButton(RadioButtons[0],0)
    Button_grp.addButton(RadioButtons[1],1)
    Button_grp.addButton(RadioButtons[2],2)
    self.RadioLayout.addWidget(RadioButtons[0],2)
    self.RadioLayout.addWidget(RadioButtons[1],2)
    self.RadioLayout.addWidget(RadioButtons[2],2)
    parametersFormLayout.addRow("Select the Grey Matter Score to paint: ",self.RadioLayout)
    RadioButtons[0].connect('clicked()',self.radio_button_FC_clicked)
    RadioButtons[1].connect('clicked()',self.radio_button_SC_clicked)
    RadioButtons[2].connect('clicked()',self.radio_button_OpacityRest)

    #
    # ROI Selector
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
	
    #
    # Apply Button
    #
    self.applyButton_3= qt.QCheckBox("Run the Optimization ")
    self.applyButton_3.toolTip = "Run the algorithm!!."
    #self.applyButton_5.enabled = True
    parametersFormLayout.addRow("Search Table", self.applyButton_3)
    # connections
    self.applyButton_3.stateChanged.connect(self.onApplyButton_3)

	
  def OnSelect_ROI(self):
      global ftone,fttwo,ft12,fss,markupNode,sphere,cylinder,logic, plane,model_C, model_S, transF
      global ft12
      logic = NeuroPathLogic()
      ftone = [0,0,0]
      fttwo = [0,0,0]
      fss = [0,0,0]
      markupNode = self.ROISelector.currentNode()
      markupNode.GetNthFiducialPosition(0,ftone)
      markupNode.GetNthFiducialPosition(1,fttwo)
      markupNode.GetNthFiducialPosition(2,fss)
      R_sphere =  self.SliderSphere.value
      R_cylinder = self.SliderCylinder.value
      H_cylinder = self.SliderCylinder_2.value
      ft1 = np.array(ftone)
      ft2 = np.array(fttwo)
      fs = np.array(fss)
      ft12  = ft1-ft2
      Y_axis = np.array([0,1,0])
      C_Cent = ((ft1+ft2)/2)
      sphere = vtk.vtkSphereSource()
      sphere.SetCenter(fs)
      sphere.SetThetaResolution(30)
      sphere.SetPhiResolution(30)
      sphere.SetRadius(R_sphere)
      Mapper = vtk.vtkPolyDataMapper()
      Mapper.SetInputConnection(sphere.GetOutputPort())
      Actor = vtk.vtkActor()
      Actor.SetMapper(Mapper)
      Actor.GetProperty().SetColor(0,1,0)
      clip = vtk.vtkClipPolyData()
      clip.SetValue(0);
      clip.GenerateClippedOutputOn();
      clip.SetInputConnection(sphere.GetOutputPort())
      plane = vtk.vtkPlane()
      normal = ft12/np.sqrt(sum(ft12*ft12))
      fend = C_Cent+(H_cylinder/2)*normal
      plane.SetNormal(normal)
      plane.SetOrigin(fend)
      clip.SetClipFunction(plane)
      #polyDataMapper = vtk.vtkPolyDataMapper()
      #polyDataMapper.SetInputConnection(clip.GetOutputPort()) 
      #actor3 = vtk.vtkActor()
      #actor3.SetMapper(polyDataMapper)
      #actor3.GetProperty().SetColor(0,1,0)
      #renderer.AddActor(actor3)
      cylinder = vtk.vtkCylinderSource()
      cylinder.SetRadius(R_cylinder)
      cylinder.SetHeight(H_cylinder)
      cylinder.SetCenter(C_Cent)
      cylinder.SetResolution(100)
      cylinder.CappingOff()

      Transform = vtk.vtkTransform()
      Transform.PostMultiply()
      Transform.Translate(-C_Cent[0],-C_Cent[1],-C_Cent[2])
      Transform.RotateWXYZ(-1*np.degrees(logic.py_ang(ft12,Y_axis)),np.cross(ft12,Y_axis))
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





  def cleanup(self):
    pass
  def onSliderChange(self):
    global R_sphere,R_cylinder,H_cylinder
    R_sphere =  self.SliderSphere.value
    R_cylinder = self.SliderCylinder.value
    H_cylinder = self.SliderCylinder_2.value
    logic.run_7()
 
  def onApplyButton_7(self):#Compute
     if self.applyButton_7.isChecked():
        logic.run_7()

  def onApplyButton_8(self):#Compute
     if self.applyButton_8.isChecked():
        logic.run_8()
		
  def onApplyButton_9(self):#Compute
     if self.applyButton_9.isChecked():
        logic.run_9()

  def OnSelect_label(self): # For loading input label, two main purposes generating grid and model handling
    self.logic.setAtlasImage(self.InputLabel_A.curretNode())
    
    global InputLabel_A,imageRec_Buck,imageCastLabel_A,imageRec_Box
    InputLabel_A = self.Label.currentNode()
    imageRec_Buck = vtk.vtkImageCast()
    imageRec_Buck.SetInputData(InputLabel_A.GetImageData())
    imageRec_Buck.Update()
    labelDims = imageRec_Buck.GetOutput().GetDimensions()
    imageCastLabel_A = vtk.vtkImageCast()# image dedicated for model handling
    imageCastLabel_A.SetInputData(InputLabel_A.GetImageData())
    imageCastLabel_A.Update()

  def OnSelect_Rec_Buck(self):# For generating grid image
    global a,b,c,space,labelDims,Rec_Buck_node,trans
    Rec_Buck_node = self.Rec_Buck.currentNode()
    Label_A_RASToIJK = vtk.vtkMatrix4x4()
    InputLabel_A.GetRASToIJKMatrix(Label_A_RASToIJK)
    Rec_Buck_node.SetRASToIJKMatrix(Label_A_RASToIJK)
    Origin = InputLabel_A.GetOrigin()
    Rec_Buck_node.SetOrigin(Origin)
    labelDims = imageRec_Buck.GetOutput().GetDimensions()
    trans = vtk.vtkTransform()
    trans.Identity()
    trans.PreMultiply()
    trans.SetMatrix(Label_A_RASToIJK)
    count = 0
    space = 5
    a = int(labelDims[0]/space)
    b = int(labelDims[1]/space)
    c = int(labelDims[2]/space)
    for r in range(0,labelDims[0]):
      for s in range(0,labelDims[1]):
         for t in range(0,labelDims[2]):
              imageRec_Buck.GetOutput().SetScalarComponentFromDouble(r,s,t,0,b*c*int(r/space)+c*int(s/space)+np.ceil(t/space))
              #imageRec_Buck.GetOutput().SetScalarComponentFromDouble(r,s,t,0,0)
    Rec_Buck_node.SetAndObserveImageData(imageRec_Buck.GetOutput())
	
	
  def OnSelect_Rec_Box(self):# For generating grid image
    self.logic.generateDamageScoreFromInput(self.InputLabel_A.currentNode())

    global a,b,c,space,labelDims,Rec_Box_node,trans
    Rec_Box_node = self.Rec_Box.currentNode()
    Label_A_RASToIJK = vtk.vtkMatrix4x4()
    InputLabel_A.GetRASToIJKMatrix(Label_A_RASToIJK)
    Rec_Box_node.SetRASToIJKMatrix(Label_A_RASToIJK)
    Origin = InputLabel_A.GetOrigin()
    Rec_Box_node.SetOrigin(Origin)
    imageRec_Box = vtk.vtkImageCast()	
    imageRec_Box.SetInputData(InputLabel_A.GetImageData())
    imageRec_Box.Update()
    #labelDims = imageRec_Box.GetOutput().GetDimensions()
    file_name_3 = qt.QFileDialog.getOpenFileName()
    f = open(file_name_3,'rb')
    csv_f = csv.reader(f)
    ar = np.genfromtxt (file_name_3, delimiter=",")
    sz = ar.shape
    global r,c,item
    r = sz[0]
    c = sz[1]   
    Lg_Table = np.zeros((r,c))
    i=0
    for row in csv_f:
        j=0
        for field in row:
          Lg_Table[i,j] = float(field)  
          j=j+1		  
        i = i+1
    Lg_Table_Skull = Lg_Table
    NumSpace = r
    trans = vtk.vtkTransform()
    trans.Identity()
    trans.PreMultiply()
    trans.SetMatrix(Label_A_RASToIJK)
    pIJK = 3*[0]
    pIJK_far = 3*[0]
    pt = 3*[int()]
    ftarget = np.array([int(Lg_Table[0,0]),int(Lg_Table[0,1]),int(Lg_Table[0,2])])
    k=0
    for s in range(1,NumSpace): 
         fstart_p = np.array([int(Lg_Table[s,0]),int(Lg_Table[s,1]),int(Lg_Table[s,2])])
         line_length = np.sqrt(sum((fstart_p-ftarget)*(fstart_p-ftarget)))
         normal_2 = (ftarget-fstart_p)/line_length
         for i in range(0,int(line_length)):
             point_p = fstart_p+i*normal_2
             point_p_far = fstart_p+(i-4)*normal_2
             trans.TransformPoint(point_p,pIJK) #RAS to IJK transformation 
             trans.TransformPoint(point_p_far,pIJK_far)
             pt[0]= int(pIJK[0])
             pt[1]= int(pIJK[1])
             pt[2]= int(pIJK[2])
             if imageRec_Main.GetOutput().GetScalarComponentAsDouble(pt[0],pt[1],pt[2],0)>0:
                k=k+1
                print(k)                   
                imageRec_Box.GetOutput().SetScalarComponentFromDouble(int(pIJK_far[0]),int(pIJK_far[1]),int(pIJK_far[2]),0,int(Lg_Table[s,9]+1))
                #imageRec_Box.GetOutput().SetScalarComponentFromDouble(int(pIJK[0]),int(pIJK[1]),int(pIJK[2]),0,int(Lg_Table[s,9]+1))
                #imageRec_Box.GetOutput().SetScalarComponentFromDouble(int(pIJK[0])+1,int(pIJK[1]),int(pIJK[2]),0,int(Lg_Table[s,9]+1))
                #imageRec_Box.GetOutput().SetScalarComponentFromDouble(int(pIJK[0])-1,int(pIJK[1]),int(pIJK[2]),0,int(Lg_Table[s,9]+1))
                #imageRec_Box.GetOutput().SetScalarComponentFromDouble(int(pIJK[0]),int(pIJK[1])+1,int(pIJK[2]),0,int(Lg_Table[s,9]+1))
                #imageRec_Box.GetOutput().SetScalarComponentFromDouble(int(pIJK[0]),int(pIJK[1])-1,int(pIJK[2]),0,int(Lg_Table[s,9]+1))
                #imageRec_Box.GetOutput().SetScalarComponentFromDouble(int(pIJK[0]),int(pIJK[1]),int(pIJK[2])+1,0,int(Lg_Table[s,9]+1))
                #imageRec_Box.GetOutput().SetScalarComponentFromDouble(int(pIJK[0]),int(pIJK[1]),int(pIJK[2])-1,0,int(Lg_Table[s,9]+1))
                Lg_Table_Skull[k,0:3] = point_p_far
                break
    Rec_Box_node.SetAndObserveImageData(imageRec_Box.GetOutput())
    Dir_name_save = qt.QFileDialog.getExistingDirectory()
    os.chdir(Dir_name_save)
    np.savetxt("Lg_Table_Skull.csv",Lg_Table_Skull,delimiter=",")
    
    print("Done")

  def OnSelect_Rec_Buck_2(self):
    global a,b,c,space,labelDims,Rec_Buck_node,imageRec_Main,trans
    Rec_Buck_node = self.Rec_Buck_2.currentNode()
    imageRec_Main = vtk.vtkImageCast()
    imageRec_Main.SetInputData(Rec_Buck_node.GetImageData())
    imageRec_Main.Update()
    labelDims = imageRec_Main.GetOutput().GetDimensions()
    Label_A_RASToIJK = vtk.vtkMatrix4x4()
    Rec_Buck_node.GetRASToIJKMatrix(Label_A_RASToIJK)
    trans = vtk.vtkTransform()
    trans.Identity()
    trans.PreMultiply()
    trans.SetMatrix(Label_A_RASToIJK)
    space = 5
    a = int(labelDims[0]/space)
    b = int(labelDims[1]/space)
    c = int(labelDims[2]/space)
    	
  def OnSelect_fb(self):
    global fb_Node,input,inPts,inLines,logic,oldTensors,oldArrays, Num_arrays, Arrays_list
    fb_Node = self.inputSelector.currentNode()
    input = fb_Node.GetPolyData()
    oldTensors = input.GetPointData().GetTensors()
    inPts =input.GetPoints()
    inLines = input.GetLines()
    ptids = vtk.vtkIdList()
    inLines.InitTraversal()
    logic = NeuroPathLogic()
    Num_arrays = input.GetPointData().GetNumberOfArrays()-1 # demark tensors which count as an array
    Arrays_list =[]
    for i in range(1,Num_arrays+1): #Array  0 is the tensors
        Arrays_list.append(input.GetPointData().GetArray(i))

  def onApplyButton_4(self):
    global Dir_name_save
    Dir_name_save = qt.QFileDialog.getExistingDirectory() 
    logic.run_4()
  def onApplyButton_5(self):
    global Directory_name
    Directory_name = qt.QFileDialog.getExistingDirectory()
    logic.run_5()
  def OnSelect_selected_fibers(self):
    global Selected_fibers_Node
    Selected_fibers_Node = self.outputSelector.currentNode()
    Selected_fibers_Node.SetAndObserveTransformNodeID(fb_Node.GetTransformNodeID())
    Selected_fibers_Node.CreateDefaultDisplayNodes()
    Selected_fibers_Node.CreateDefaultStorageNode()
    Selected_fibers_Node.SetName('Selected_Fibers')
    
    #logic.run_5()#Copy selected fibers
  def OnSelect_Mod(self):
    global Node_Models, node_names, num_models
    node_names =[]
    Node_Models = self.Models.currentNode()
    num_models = Node_Models.GetNumberOfChildrenNodes()
    for t in range(0,num_models):
       child = Node_Models.GetNthChildNode(t)
       node = child.GetModelNode()	   
       name = node.GetName()
       name = name[6:name.__len__()]# Remove first part of the name which is usually Model_	   
       node_names.append(int(name[0:name.index('_')]))
  def radio_button_FC_clicked(self):
       for t in range(0,num_models):
           child = Node_Models.GetNthChildNode(t)
           node = child.GetModelNode()
           name = node.GetName()        
           d_node =node.GetDisplayNode()
           d_node.SetOpacity(FC_GM[t])
           if FC_GM[t]<0.5:
                    d_node.SetColor(2*FC_GM[t],1,0)
           else:
                    d_node.SetColor(1,2-2*FC_GM[t],0)
  def radio_button_SC_clicked(self):
       for t in range(0,num_models):
           child = Node_Models.GetNthChildNode(t)
           node = child.GetModelNode()
           name = node.GetName()        
           d_node =node.GetDisplayNode()
           d_node.SetOpacity(SC_GM[t])
           if SC_GM[t]<0.5:
                    d_node.SetColor(2*SC_GM[t],0,1)
           else:
                    d_node.SetColor(1,0,2-2*SC_GM[t])

  def radio_button_OpacityRest(self):
       for t in range(0,num_models):
           child = Node_Models.GetNthChildNode(t)
           node = child.GetModelNode()
           name = node.GetName()        
           d_node =node.GetDisplayNode()
           d_node.SetOpacity(1)

					
  def OnSelect_OptBox(self):
    global bounds,markupNode_2
    bounds = [0,0,0,0,0,0]
    markupNode_2 = self.ROISelector_2.currentNode()
    markupNode_2.GetRASBounds(bounds)
    
	
  def onApplyButton_3(self):
     if self.applyButton_3.isChecked():
        logic.run_3()
	    

#
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
    X_1 = int(bounds[0]/10)
    X_2 = int(bounds[1]/10)
    Y_1 = int(bounds[2]/10)
    Y_2 = int(bounds[3]/10)
    Z_1 = int(bounds[4]/10)
    Z_2 = int(bounds[5]/10)
    Num_sample = (X_2-X_1+1)*(Y_2-Y_1+1)+2*((Y_2-Y_1+1)*(Z_2-Z_1+1)+(X_2-X_1+1)*(Z_2-Z_1+1))
    Lg_Table = np.zeros((Num_sample+1,10))
    Lg_Table[0,4] = X_1
    Lg_Table[0,5] = X_2
    Lg_Table[0,6] = Y_1
    Lg_Table[0,7] = Y_2
    Lg_Table[0,8] = Z_1
    Lg_Table[0,9] = Z_2
    k = 1
    for i in range(0,(X_2-X_1)+1):
     for j in range(0,(Y_2-Y_1)+1):
        Lg_Table[k,0:3] = [(X_1+i)*10,(Y_1+j)*10,Z_2*10]
        k=k+1
      
    for i in range(0,(Y_2-Y_1)+1):
     for j in range(0,(Z_2-Z_1)+1):
        Lg_Table[k,0:3] = [X_1*10,(Y_1+i)*10,(Z_1+j)*10]
        k=k+1  
    		
    for i in range(0,(Y_2-Y_1)+1):
     for j in range(0,(Z_2-Z_1)+1):
        Lg_Table[k,0:3] = [X_2*10,(Y_1+i)*10,(Z_1+j)*10]
        k=k+1
              
    for i in range(0,(X_2-X_1)+1):
     for j in range(0,(Z_2-Z_1)+1):
        Lg_Table[k,0:3] = [(X_1+i)*10,Y_1*10,(Z_1+j)*10]
        k=k+1
        
    for i in range(0,(X_2-X_1)+1):
     for j in range(0,(Z_2-Z_1)+1):
        Lg_Table[k,0:3] = [(X_1+i)*10,Y_2*10,(Z_1+j)*10]
        k=k+1
      
    markupNode.GetNthFiducialPosition(0,ftone)
    markupNode.GetNthFiducialPosition(2,fss)
    ft1 = np.array(ftone)
    fs = np.array(fss)
    Lg_Table[0,0:3] = ft1  
    for m in range(1,Num_sample+1):     

      Score=[0,0,0,0,0,0]      
      ft2 = np.array(Lg_Table[m,0:3])      
      ft12  = (ft2-ft1)
      H_cylinder_SEARCH = np.sqrt(sum(ft12*ft12))
      normal = ft12/H_cylinder_SEARCH
      C_Cent = ft1 + (H_cylinder_SEARCH/2)*normal
      fend = C_Cent+(H_cylinder_SEARCH/2)*normal
      fstart = C_Cent-(H_cylinder_SEARCH/2)*normal
      cylinder.SetCenter(C_Cent) 	  
      plane.SetNormal(normal)
      plane.SetOrigin(fend)
      Y_axis = np.array([0,1,0])
      sphere.SetCenter(fs)
      R_cylinder = 5   # Cylinder is SEt to be constant
      cylinder.SetRadius(R_cylinder)
      cylinder.SetHeight(H_cylinder_SEARCH)
      sphere.SetRadius(R_sphere)
      Transform = vtk.vtkTransform()
      Transform.PostMultiply()

      Transform = vtk.vtkTransform()
      Transform.PostMultiply()    
      Transform.Translate(-C_Cent[0],-C_Cent[1],-C_Cent[2])
      Transform.RotateWXYZ(-1*np.degrees(logic.py_ang(ft12,Y_axis)),np.cross(ft12,Y_axis))
      Transform.Translate(C_Cent[0],C_Cent[1],C_Cent[2])
      #transF = vtk.vtkTransformPolyDataFilter()
      #transF.SetInputConnection(cylinder.GetOutputPort())
      transF.SetTransform(Transform)
      transF.Update()

      pt = 3*[int()]
      Label_num_a =[]
      Label_l = []
      Num_OF_points_per_label=[]
      Num_of_points = 0
      inPtr =[]
      inPtr2 =[]
      Num_pt = int((H_cylinder_SEARCH))
      Cyl_pt = np.zeros((Num_pt,3))
      P_Source = vtk.vtkPointSource()
      P_Source.SetNumberOfPoints(10000)## CHANGED FROM 1000000 to 10000 to speed up the process
      P_Source.SetRadius(H_cylinder_SEARCH/2+20)
      P_Source.SetCenter(C_Cent)
      P_Source.Update()
      pointLocator = vtk.vtkPointLocator()
      pointLocator.SetDataSet(P_Source.GetOutput())
      pointLocator.AutomaticOn()
      pointLocator.SetNumberOfPointsPerBucket(1)
      pointLocator.BuildLocator()
      result = vtk.vtkIdList()
      pIJK = 3*[0]
      for i in range(0,Num_pt):
         Cyl_pt[i,:] = fstart+(i)*normal
         pointLocator.FindPointsWithinRadius(R_cylinder,fstart+(i)*normal,result)
         for j in range(0,result.GetNumberOfIds()): 
               point = P_Source.GetOutput().GetPoint(result.GetId(j))		
               trans.TransformPoint(point,pIJK) #RAS to IJK transformation 
               pt[0]= int(pIJK[0])
               pt[1]= int(pIJK[1])
               pt[2]= int(pIJK[2])
               if pt[0] < 0 or pt[1] < 0 or pt[2] < 0 or pt[0] >= labelDims[0] or pt[1] >= labelDims[1] or pt[2] >= labelDims[2]:
                     continue
               inPtr = imageRec_Main.GetOutput().GetScalarComponentAsDouble(pt[0],pt[1],pt[2],0)
               Label_num_a.append(inPtr)
               inPtr2 = imageCastLabel_A.GetOutput().GetScalarComponentAsDouble(pt[0],pt[1],pt[2],0)
               Label_l.append(inPtr2)
               Score[0] = Score[0]+int(Aparc2[Aparc.index(int(inPtr2))])
               Score[1] = Score[1]+int(Aparc3[Aparc.index(int(inPtr2))])
               
      Label_num_a = list(set(Label_num_a))# To remove redundancies
      #print(Label_num_a)
      length = Label_num_a.__len__()
      print('length: '+str(length))    
      addLine = 0*(Mem_array[:,0])
      for s in range(length):
        if Label_num_a[s] in Label_LUT:
            exist = Mem_array[:,Label_LUT.index(Label_num_a[s])]
            addLine = addLine+(exist)
      addLine = addLine>0
      numNewCells = int(sum(addLine))
      print('Number of lines: '+str(numNewCells))
      address = np.where(addLine==1)[:][0]
      numNewPts = int(sum(Mem_points[address]))
      print('Number of points: '+str(numNewPts))
      if numNewCells>0:
        newArrays =[]
        ptId = int()
        pcount = int()
        inLines.InitTraversal()
        address1 = np.ndarray.tolist(address)
        for inCellId in address1:
          npts = int(Mem_points[inCellId])
          if (oldTensors):
              for l in range(0,Num_arrays):
                  Score[l+2] = Score[l+2]+int(Arrays_list[l].GetTuple(ptId)[0])
          ptId = ptId+npts

        print('ptId: '+str(ptId))
      Label_l = list(set(Label_l)) # Remove duplicate list items
      print "Score: "+str((Score))
      Lg_Table[m,3:(3+Score.__len__())] = Score	
      Lg_Table[m,9] = H_cylinder_SEARCH
    Dir_name_save = qt.QFileDialog.getExistingDirectory() 
    os.chdir(Dir_name_save)
    np.savetxt("Lg_Table.csv",Lg_Table,delimiter=",")
     
		
 
  def run_7(self): #Run the computation for path refresh upon slider change
    ft1 = np.array(ftone)
    ft2 = np.array(fttwo)
    ft12  = ft1-ft2
    markupNode.GetNthFiducialPosition(0,ftone)
    markupNode.GetNthFiducialPosition(1,fttwo)
    markupNode.GetNthFiducialPosition(2,fss)
    ft1 = np.array(ftone)
    ft2 = np.array(fttwo)
    fs = np.array(fss)
    C_Cent = ((ft1+ft2)/2)
    cylinder.SetCenter(C_Cent)    
    ft12  = (ft2-ft1)
    normal = ft12/np.sqrt(sum(ft12*ft12))
    fend = C_Cent+(H_cylinder/2)*normal
    plane.SetNormal(normal)
    plane.SetOrigin(fend)
    Y_axis = np.array([0,1,0])
    sphere.SetCenter(fs)
    cylinder.SetRadius(R_cylinder)
    cylinder.SetHeight(H_cylinder)
    sphere.SetRadius(R_sphere)
    Transform = vtk.vtkTransform()
    Transform.PostMultiply()    
    Transform.Translate(-C_Cent[0],-C_Cent[1],-C_Cent[2])
    Transform.RotateWXYZ(-1*np.degrees(logic.py_ang(ft12,Y_axis)),np.cross(ft12,Y_axis))
    Transform.Translate(C_Cent[0],C_Cent[1],C_Cent[2])
    #transF = vtk.vtkTransformPolyDataFilter()
    #transF.SetInputConnection(cylinder.GetOutputPort())
    transF.SetTransform(Transform)
    transF.Update()
    
    

		
  def run_8(self):#Run the computation for collosion detection and score computation
    global Score, newArrays, Num_arrays
    Score=[]
    Score.append(0)#This first column is reserved
    GM_Counter = [0]*Aparc.__len__()   #This counter is for number of points inside of access path from each label region 
    markupNode.GetNthFiducialPosition(0,ftone)
    markupNode.GetNthFiducialPosition(1,fttwo)
    markupNode.GetNthFiducialPosition(2,fss)
    ft1 = np.array(ftone)
    ft2 = np.array(fttwo)
    fs = np.array(fss)
    C_Cent = ((ft1+ft2)/2)
    cylinder.SetCenter(C_Cent)    
    ft12  = (ft2-ft1)
    normal = ft12/np.sqrt(sum(ft12*ft12))
    fend = C_Cent+(H_cylinder/2)*normal
    fstart = C_Cent-(H_cylinder/2)*normal
    pt = 3*[int()]
    Label_num_a =[]
    Label_l = []
    Num_OF_points_per_label=[]
    Num_of_points = 0
    inPtr =[]
    inPtr2 =[]
    Num_pt = int((H_cylinder))
    Cyl_pt = np.zeros((Num_pt,3))
    P_Source = vtk.vtkPointSource()
    P_Source.SetNumberOfPoints(1000000)
    P_Source.SetRadius(H_cylinder/2+20)
    P_Source.SetCenter(C_Cent)
    P_Source.Update()
    pointLocator = vtk.vtkPointLocator()
    pointLocator.SetDataSet(P_Source.GetOutput())
    pointLocator.AutomaticOn()
    pointLocator.SetNumberOfPointsPerBucket(1)
    pointLocator.BuildLocator()
    P_Source_2 = vtk.vtkPointSource()
    P_Source_2.SetNumberOfPoints(500000)
    P_Source_2.SetRadius(R_sphere)
    P_Source_2.SetCenter(fs)
    P_Source_2.Update()
    pointLocator2 = vtk.vtkPointLocator()
    pointLocator2.SetDataSet(P_Source_2.GetOutput())
    pointLocator2.AutomaticOn()
    pointLocator2.SetNumberOfPointsPerBucket(1)
    pointLocator2.BuildLocator()
    result = vtk.vtkIdList()
    pIJK = 3*[0]
    for i in range(0,Num_pt):
        Cyl_pt[i,:] = fstart+(i)*normal
        pointLocator.FindPointsWithinRadius(R_cylinder,fstart+(i)*normal,result)
        for j in range(0,result.GetNumberOfIds()): 
               point = P_Source.GetOutput().GetPoint(result.GetId(j))		
               trans.TransformPoint(point,pIJK) #RAS to IJK transformation 
               pt[0]= int(pIJK[0])
               pt[1]= int(pIJK[1])
               pt[2]= int(pIJK[2])
               if pt[0] < 0 or pt[1] < 0 or pt[2] < 0 or pt[0] >= labelDims[0] or pt[1] >= labelDims[1] or pt[2] >= labelDims[2]:
                     continue
               inPtr = imageRec_Main.GetOutput().GetScalarComponentAsDouble(pt[0],pt[1],pt[2],0)
               Label_num_a.append(inPtr)
               inPtr2 = imageCastLabel_A.GetOutput().GetScalarComponentAsDouble(pt[0],pt[1],pt[2],0)
               Label_l.append(inPtr2)
               Score[0] = Score[0]+Aparc2[Aparc.index(int(inPtr2))]
               GM_Counter[Aparc.index(int(inPtr2))] = GM_Counter[Aparc.index(int(inPtr2))]+1
    pointLocator2.FindPointsWithinRadius(R_sphere,fs,result)
    for j in range(0,result.GetNumberOfIds()): 
               point = P_Source_2.GetOutput().GetPoint(result.GetId(j))		
               trans.TransformPoint(point,pIJK) #RAS to IJK transformation 
               pt[0]= int(pIJK[0])
               pt[1]= int(pIJK[1])
               pt[2]= int(pIJK[2])
               if pt[0] < 0 or pt[1] < 0 or pt[2] < 0 or pt[0] >= labelDims[0] or pt[1] >= labelDims[1] or pt[2] >= labelDims[2]:
                     continue
               inPtr = imageRec_Main.GetOutput().GetScalarComponentAsDouble(pt[0],pt[1],pt[2],0)
               Label_num_a.append(inPtr)
               inPtr2 = imageCastLabel_A.GetOutput().GetScalarComponentAsDouble(pt[0],pt[1],pt[2],0)
               Label_l.append(inPtr2)
               Score[0] = Score[0]+Aparc2[Aparc.index(int(inPtr2))]
               GM_Counter[Aparc.index(int(inPtr2))] = GM_Counter[Aparc.index(int(inPtr2))]+1			   
    Label_num_a = list(set(Label_num_a))# To remove redundancies
    #Label_num_a = Label_LUT
    print(Label_num_a)
    #print(Cyl_pt)
    length = Label_num_a.__len__()
    print('length: '+str(length))    
    addLine = 0*(Mem_array[:,0])
    for s in range(length):
        if Label_num_a[s] in Label_LUT:
            exist = Mem_array[:,Label_LUT.index(Label_num_a[s])]
            addLine = addLine+(exist)
    addLine = addLine>0
    numNewCells = int(sum(addLine))
    print('Number of lines: '+str(numNewCells))
    address = np.where(addLine==1)[:][0]
    numNewPts = int(sum(Mem_points[address]))
    print('Number of points: '+str(numNewPts))
    if numNewCells>0:
       outFibers = vtk.vtkPolyData()
       points = vtk.vtkPoints()
       points.Allocate(numNewPts)
       outFibers.SetPoints(points)
       outFibersCellArray = vtk.vtkCellArray()
       outFibersCellArray.Allocate(numNewPts+numNewCells)
       outFibers.SetLines(outFibersCellArray)
       newTensors = vtk.vtkFloatArray()
       newTensors.SetNumberOfComponents(9)
       newTensors.Allocate(9*numNewPts)
       outFibers.GetPointData().SetTensors(newTensors)
       newArrays =[]
       for j in range(0,Num_arrays): 
            newArrays.append(vtk.vtkFloatArray())
            newArrays[j].SetName(Arrays_list[j].GetName())
            newArrays[j].SetNumberOfComponents(1)
            newArrays[j].Allocate(numNewPts)
            outFibers.GetPointData().AddArray(newArrays[j])
            Score.append(0)
       ptId = int()
       pcount = int()
       inLines.InitTraversal()
       address1 = np.ndarray.tolist(address)
       killLoopFlag = False
       for inCellId in address1:
         npts = int(Mem_points[inCellId])
         outFibersCellArray.InsertNextCell(npts,range(ptId,ptId+npts))
         st_pt = int(sum(Mem_points[0:(inCellId)]))
         points.InsertPoints(ptId,npts,st_pt,inPts)
         if (oldTensors):
             newTensors.InsertTuples(ptId,npts,st_pt,oldTensors)
             for l in range(0,Num_arrays):
                  newArrays[l].InsertTuples(ptId,npts,st_pt,Arrays_list[l])
                  Score[l+1] = Score[l+1]+Arrays_list[l].GetTuple(ptId)[0]
         ptId = ptId+npts
         if (killLoopFlag):
            print('Canceled')
            break
       print('ptId: '+str(ptId))
    Selected_fibers_Node.SetAndObservePolyData(outFibers)
    # This part is for Model selectiond/deselction
    Label_l = list(set(Label_l)) # Remove duplicate list items
    print "Selected Labels: "+str(Label_l)
    print "Counter of Labels: "+str(GM_Counter)
    if Label_l.__len__()>0:
       self.Model_vis(Label_l)
 
  def py_ang(self,v1, v2):
   """ Returns the angle in radians between vectors 'v1' and 'v2'    """
   cosang = np.dot(v1, v2)
   sinang = la.norm(np.cross(v1, v2))
   return np.arctan2(sinang, cosang)
  
  def run_4(self): # Calculating and saving Mem_array, Mem_points and Label_LUT 
    global Mem_array, Mem_points, Label_LUT,trans
    Label_A_RASToIJK = vtk.vtkMatrix4x4()
    Rec_Buck_node.GetRASToIJKMatrix(Label_A_RASToIJK)
    trans = vtk.vtkTransform()
    trans.Identity();
    trans.PreMultiply();
    trans.SetMatrix(Label_A_RASToIJK)
    numPts = inPts.GetNumberOfPoints()
    numLines = inLines.GetNumberOfCells()
    pts = vtk.vtkIdList()
    addLines = []
    LabelScala_range = imageRec_Main.GetOutput().GetScalarRange()#Should be eaqual to a*b*c
    Mem_array = np.full((numLines,a*b*c), False, dtype=bool)
    Mem_points = np.zeros(numLines)
    Label_LUT = []
    inLines.InitTraversal()
    print(numLines)
    for inCellId in range(0,numLines):
        inLines.GetNextCell(pts);
        npts = int(pts.GetNumberOfIds())
        pIJK = 3*[0]
        pt = 3*[int()]
        p = 3*[0]
        inPtr =[]
        if npts <19:
           print('short line!!!')
           print(inCellId)
        else:
           for j in range(0,npts):
             p = inPts.GetPoint(pts.GetId(j))
             trans.TransformPoint(p,pIJK)
             pt[0]= int(pIJK[0])
             pt[1]= int(pIJK[1])
             pt[2]= int(pIJK[2])
             if pt[0] < 0 or pt[1] < 0 or pt[2] < 0 or pt[0] >= labelDims[0] or pt[1] >= labelDims[1] or pt[2] >= labelDims[2]:
                 continue
             inPtr = imageRec_Main.GetOutput().GetScalarComponentAsDouble(pt[0],pt[1],pt[2],0)
             if inPtr<1:
                continue
             if inPtr in Label_LUT:
                  Mem_array[inCellId,Label_LUT.index(inPtr)] = True
             else: 
                  Label_LUT.append(inPtr)
                  Mem_array[inCellId,Label_LUT.index(inPtr)] = True
            
        Mem_points[inCellId] = npts
	
    os.chdir(Dir_name_save)
    np.save('Mem_array.npy', Mem_array[:,0:Label_LUT.__len__()])
    np.save('Mem_points.npy',Mem_points)
    f = open("Label_LUT", "w")
    for item in Label_LUT:
         f.write("%s " % int(item))
    f.close()
  def run_5(self):
    global Mem_array, Mem_points, Label_LUT, FC_GM, SC_GM, Aparc, Aparc2, Aparc3
    Mem_array = np.load(Directory_name+'/Mem_array.npy')
    Mem_points = np.load(Directory_name+'/Mem_points.npy')
    f=open(Directory_name+'/Label_LUT', 'rb') 
    Label_LUT = map(int,f.read().split())
    f=open(Directory_name+'/FC_GM', 'rb')
    FC_GM = map(float,f.read().split())
    f=open(Directory_name+'/SC_GM', 'rb')
    SC_GM = map(float,f.read().split())
    f=open(Directory_name+'/Aparc', 'rb')  # All of the Labels of the Atlas regions
    Aparc = map(float,f.read().split())
    f=open(Directory_name+'/Aparc2', 'rb') # FC Absolute value of the Atlas regions
    Aparc2 = map(float,f.read().split())
    f=open(Directory_name+'/Aparc3', 'rb') # FC Absolute value of the Atlas regions
    Aparc3 = map(float,f.read().split())
  def Model_vis(self,Labels):
        len = Labels.__len__()# Labels are inpt data values
        print("Number of intersected Regions : "+str(len))
        for s in range(0,node_names.__len__()):
            child = Node_Models.GetNthChildNode(s)
            node = child.GetModelNode()
            node.SetDisplayVisibility(0)            
        for r in range(0,len):
            if Labels[r] in node_names:
               child = Node_Models.GetNthChildNode(node_names.index(Labels[r]))
               node = child.GetModelNode()
               d_node = node.GetDisplayNode()
               #print(node.GetName())
               node.SetDisplayVisibility(1)
               d_node.SetOpacity(0.75)    
	
  def run_9(self): # This part is for calculation of the damage or score and saving it in an array with the time in the file name
        
        		
		#Score =[4.5,5.5,6.5,4.2]
        global Dir_name_save  # Temporary , should be removed later
        Dir_name_save = qt.QFileDialog.getExistingDirectory() 
        os.chdir(Dir_name_save)
        Timenow = datetime.datetime.now()
        File_name = "Score_"+str(Timenow.hour)+"_"+str(Timenow.minute)+"_"+str(Timenow.second)+"_"+str(Timenow.microsecond)
        f = open(File_name, "w")
        for item in Score:
          f.write("%s " % float(item))
        f.close()
 
class NeuroPathTest(ScriptedLoadableModuleTest):
  """
  This is the test case for your scripted module.
  Uses ScriptedLoadableModuleTest base class, available at:
  https://github.com/Slicer/Slicer/blob/master/Base/Python/slicer/ScriptedLoadableModule.py
  """
 
  def setUp(self):
    """ Do whatever is needed to reset the state - typically a scene clear will be enough.
    """
    slicer.mrmlScene.Clear(0)
 
  def runTest(self):
    """Run as few or as many tests as needed here.
    """
    self.setUp()
    self.test_NeuroPath1()
 
  def test_NeuroPath1(self):
    """ Ideally you should have several levels of tests.  At the lowest level
    tests should exercise the functionality of the logic with different inputs
    (both valid and invalid).  At higher levels your tests should emulate the
    way the user would interact with your code and confirm that it still works
    the way you intended.
    One of the most important features of the tests is that it should alert other
    developers when their changes will have an impact on the behavior of your
    module.  For example, if a developer removes a feature that you depend on,
    your test should break so they know that the feature is needed.
    """
 
    self.delayDisplay("Starting the test")
    #
    # first, get some data
    #
    import urllib
    downloads = (
        ('http://slicer.kitware.com/midas3/download?items=5767', 'FA.nrrd', slicer.util.loadVolume),
        )
 
    for url,name,loader in downloads:
      filePath = slicer.app.temporaryPath + '/' + name
      if not os.path.exists(filePath) or os.stat(filePath).st_size == 0:
        logging.info('Requesting download %s from %s...\n' % (name, url))
        urllib.urlretrieve(url, filePath)
      if loader:
        logging.info('Loading %s...' % (name,))
        loader(filePath)
    self.delayDisplay('Finished with download and loading')
 
    volumeNode = slicer.util.getNode(pattern="FA")
    logic = NeuroPathLogic()
    self.assertTrue( logic.hasImageData(volumeNode) )
    self.delayDisplay('Test passed!')
