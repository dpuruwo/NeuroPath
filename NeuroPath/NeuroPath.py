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

    #
    # Parameters Area
    #
    parametersCollapsibleButton = ctk.ctkCollapsibleButton()
    parametersCollapsibleButton.text = "Parameters"
    parametersCollapsibleButton.collapsed = False
    self.parametersList = parametersCollapsibleButton
    self.layout.addWidget(parametersCollapsibleButton)
 
    # Layout within the dummy collapsible button
    parametersFormLayout = qt.QFormLayout(parametersCollapsibleButton)
    #
    # Fiducial  List
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
    parametersFormLayout.addRow("Fiducial List", self.ROISelector)
    self.ROISelector.connect("currentNodeChanged(vtkMRMLNode*)", self.OnSelect_ROI)

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
    self.SliderCylinder.setToolTip("Cylinder Radius")
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
    self.SliderCylinder_2.setToolTip("Cylinder Height")
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

    #
    # Input the Labeled Volume (???)
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
  

  def onSliderChange(self):
      global R_sphere,R_cylinder,H_cylinder
      R_sphere =  self.SliderSphere.value
      R_cylinder = self.SliderCylinder.value
      H_cylinder = self.SliderCylinder_2.value
      #logic.run_7()

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
      markupNode.GetNthFiducialPosition(0,fiducial1) #former ftone
      markupNode.GetNthFiducialPosition(1,fiducial2) #former fttwo
      markupNode.GetNthFiducialPosition(2,fiducial3) #former fss
      #print ftone 
      #print fttwo
      #print fss 
     
      ## updates sphere and cylinder values
      R_sphere =  self.SliderSphere.value
      R_cylinder = self.SliderCylinder.value
      H_cylinder = self.SliderCylinder_2.value
      
      
      #
      #creates an array
      #
      f1 = np.array(fiducial1)   #former ft1
      f2 = np.array(fiducial2)   #former ft2
      f3 = np.array(fiducial3)   #former fs
      deltaf1f2 = f1-f2  #position of fiducial 1 - fiducial 2. Why??
      Y_axis = np.array([0,1,0])
      C_Cent = ((f1+f2)/2) 
     
      #
      ##
      #

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

      ## not sure what these transformations are supposed to do
      Transform = vtk.vtkTransform()
      Transform.PostMultiply()
      Transform.Translate(-C_Cent[0],-C_Cent[1],-C_Cent[2])
      #np.degrees converts radians to degrees
      Transform.RotateWXYZ(-1*np.degrees(logic.py_ang(deltaf1f2,Y_axis)),np.cross(deltaf1f2,Y_axis))
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

  ### holds the ROI coordinates found in Annotation Properities
  def OnSelect_OptBox(self):
    global bounds,markupNode_2
    bounds = [0,0,0,0,0,0] # [L,R,P,A,I,S]
    markupNode_2 = self.ROISelector_2.currentNode()
    markupNode_2.GetRASBounds(bounds)
     #print bounds

  def OnSelect_label(self): # For loading input label, two main purposes generating grid and model handling
    # GetDimentions: number of points on each axis
    # imageCastLabel_A gets used in run_3, inPtr2  
    #labeldims also in run_3
    logic = NeuroPathLogic()
    self.logic.setAtlasImage(self.InputLabel_A.curretNode())
    
    global InputLabel_A,imageRec_Buck,imageCastLabel_A,imageRec_Box
    InputLabel_A = self.Label.currentNode()
    imageRec_Buck = vtk.vtkImageCast()
    imageRec_Buck.SetInputData(InputLabel_A.GetImageData())
    imageRec_Buck.Update()
    
    labelDims = imageRec_Buck.GetOutput().GetDimensions()
    print labelDims
    
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
    trans.Identity     
    trans.PreMultiply()
    trans.SetMatrix(Label_A_RASToIJK)
    space = 5
    a = int(labelDims[0]/space)
    b = int(labelDims[1]/space)
    c = int(labelDims[2]/space)

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
    # X_1 left corner coordinate of ROI Box (L)
    # X_2 right corner coordinate of ROI Box (R)
    # Y_1 posterior corner of ROI Box (P)
    # Y_2 anterior corner of ROI box (A)
    # Z_1 inferior corner of ROI box (I)
    # Z_2 superior corner of ROI box (S)

    X_1 = int(bounds[0]/10)  
    X_2 = int(bounds[1]/10)
    Y_1 = int(bounds[2]/10)
    Y_2 = int(bounds[3]/10)
    Z_1 = int(bounds[4]/10)
    Z_2 = int(bounds[5]/10)

    # (R - L + 1) * (A - P + 1) + 2*((A - P + 1 )*(S - I + 1) + (R - L + 1)*(S - I + 1)) 
    Num_sample = (X_2-X_1+1)*(Y_2-Y_1+1)+2*((Y_2-Y_1+1)*(Z_2-Z_1+1)+(X_2-X_1+1)*(Z_2-Z_1+1)) 
   

    ## Adds to Lg_Table the coordinates of the box
    # Left corner is stored in first row, 5th column and so on
    Lg_Table = np.zeros((Num_sample+1,10)) # 10 columns 
    Lg_Table[0,4] = X_1 
    Lg_Table[0,5] = X_2
    Lg_Table[0,6] = Y_1
    Lg_Table[0,7] = Y_2
    Lg_Table[0,8] = Z_1
    Lg_Table[0,9] = Z_2
    

    #loops i*j times 
    #
    k = 1 # counter
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
    markupNode.GetNthFiducialPosition(0,fiducial1)
    markupNode.GetNthFiducialPosition(2,fiducial3)
    f1 = np.array(fiducial1)
    f3 = np.array(fiducial3)
    Lg_Table[0,0:3] = f1 #location of target 1
    for m in range(1,Num_sample+1):

      Score = [0,0,0,0,0,0]
      f2=np.array[Lg_Table[m,0:3]]
      deltaf1f2 = (f2-f1)
      
      H_cylinder_SEARCH = np.sqrt(sum(deltaf1f2*deltaf1f2))
      normal = deltaf1f2/H_cylinder_SEARCH
      C_Cent = f1 + (H_cylinder_SEARCH/2)*normal
      fend = C_Cent+(H_cylinder_SEARCH/2)*normal
      fstart = C_Cent-(H_cylinder_SEARCH/2)*normal
      cylinder.SetCenter(C_Cent)
      plane.SetNormal(normal)
      plane.SetOrigin(fend)
      Y_axis = np.array([0,1,0])
      sphere.SetCenter(f3)
      R_cylinder = 5 #Cylinder is set to be constant
      cylinder.SetRadius(R_cylinder)
      cylinder.SetHeight(H_cylinder_SEARCH)
      sphere.SetRadius(R_sphere)
      

      Transform = vtk.vtkTransform()
      Transform.PostMultiply()
      Transform.Translate(-C_Cent[0],-C_Cent[1],-C_Cent[2])
      Transform.RotateWXYZ(-1*np.degrees(logic.py_ang(deltaf1f2,Y_axis)),np.cross(delt,Y_axis))
      Transform.Translate(C_Cent[0],C_Cent[1],C_Cent[2])
      transF.SetTransform(Transform)
      transF.Update()

      pt = 3*[int()] # [0,0,0]
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
     
      pIJK = 3*[0] #[0,0,0]
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

    ### where are these from?
    f=open(Directory_name+'/Aparc', 'rb')  # All of the Labels of the Atlas regions
    Aparc = map(float,f.read().split())
    f=open(Directory_name+'/Aparc2', 'rb') # FC Absolute value of the Atlas regions
    Aparc2 = map(float,f.read().split())
    f=open(Directory_name+'/Aparc3', 'rb') # FC Absolute value of the Atlas regions
    Aparc3 = map(float,f.read().split())

  def py_ang(self,v1, v2): # ask Dr. Eagleson
   # measure euclidian distance between y axis and deltaf1f2 (fiducial 1 - fiducial 2)
  # then cross product normalized 
  #the cross product can be thought of as a measure of perpendicularity in the same way that the dot product is a measure of parallelism.

  # v1 is (deltaf1f2)
  # v2 is Y_axis ->  Y_axis = np.array([0,1,0])

  # product of the Euclidean magnitudes of the two vectors;
  # v1 . v2 = v1v2cos(angle between) # measure "how much" one vector lies along another
    cosang = np.dot(v1, v2) 

   #normalize magnitude (or length) of deltaf1f2 and Y_axis vectors
    sinang = la.norm(np.cross(v1, v2)) 
    return np.arctan2(sinang, cosang)

  
### LEFT HERE
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



