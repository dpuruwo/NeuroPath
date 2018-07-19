import os
import vtk, qt, ctk, slicer
from slicer.ScriptedLoadableModule import *
import logging
import numpy as np
import numpy.linalg as la
import csv
 
#
# NeuroPath
#
class NeuroPath(ScriptedLoadableModule):
  def __init__(self, parent):
    ScriptedLoadableModule.__init__(self, parent)

    self.parent.title = "NeuroPath"
    self.parent.categories = ["Neurosurgery"]
    self.parent.dependencies = ["Markups"]
    self.parent.contributors = ["Daiana Pur (Western University), Adam Rankin (Robarts Research Institute)"]
    self.parent.helpText = """COMPLETE ME!"""
    self.parent.acknowledgementText = """COMPLETE ME!"""

    self.logic = NeuroPathLogic()

#
# NeuroPathWidget
#
class NeuroPathWidget(ScriptedLoadableModuleWidget):

  def setup(self):
    ScriptedLoadableModuleWidget.setup(self)

    # Inputs
    inputsCollapsibleButton = ctk.ctkCollapsibleButton()
    inputsCollapsibleButton.text = "Inputs"
    self.layout.addWidget(inputsCollapsibleButton)
    inputsFormLayout = qt.QFormLayout(inputsCollapsibleButton)

    self.modelPointsSelector = slicer.qMRMLNodeComboBox()
    self.modelPointsSelector.nodeTypes = ["vtkMRMLMarkupsFiducialNode"]
    self.modelPointsSelector.selectNodeUponCreation = True
    self.modelPointsSelector.addEnabled = True
    self.modelPointsSelector.removeEnabled = True
    self.modelPointsSelector.noneEnabled = True
    self.modelPointsSelector.showHidden = False
    self.modelPointsSelector.showChildNodeTypes = False
    self.modelPointsSelector.setMRMLScene(slicer.mrmlScene)
    self.modelPointsSelector.setToolTip("Pick the fiducials that define the cylinder top, bottom, and sphere centre.")
    inputsFormLayout.addRow("Model Fiducials: ", self.modelPointsSelector)
    self.modelPointsSelector.connect("currentNodeChanged(vtkMRMLNode*)", self.onModelPointsSelected)

    self.qSlider_SphereRadius = ctk.ctkSliderWidget()
    self.qSlider_SphereRadius.singleStep = 0.1
    self.qSlider_SphereRadius.minimum = 0
    self.qSlider_SphereRadius.maximum = 20
    self.qSlider_SphereRadius.value = 10
    self.qSlider_SphereRadius.setToolTip("Set sphere radius.")
    self.qSlider_SphereRadius.valueChanged.connect(self.onSliderChange)
    inputsFormLayout.addRow("Sphere Radius: ", self.qSlider_SphereRadius)

    self.qSlider_CylinderRadius = ctk.ctkSliderWidget()
    self.qSlider_CylinderRadius.singleStep = 0.1
    self.qSlider_CylinderRadius.minimum = 0
    self.qSlider_CylinderRadius.maximum = 20
    self.qSlider_CylinderRadius.value = 7.5
    self.qSlider_CylinderRadius.setToolTip("Set cylinder radius.")
    self.qSlider_CylinderRadius.valueChanged.connect(self.onSliderChange)
    inputsFormLayout.addRow("Cylinder Radius: ", self.qSlider_CylinderRadius)

    self.qSlider_CylinderHeight = ctk.ctkSliderWidget()
    self.qSlider_CylinderHeight.singleStep = 0.1
    self.qSlider_CylinderHeight.minimum = 5
    self.qSlider_CylinderHeight.maximum = 120
    self.qSlider_CylinderHeight.value = 30
    self.qSlider_CylinderHeight.setToolTip("Set cylinder height.")
    self.qSlider_CylinderHeight.valueChanged.connect(self.onSliderChange)
    inputsFormLayout.addRow("Cylinder Height", self.qSlider_CylinderHeight)

    self.updateFibersButton = qt.QPushButton("Update the fibers")
    self.updateFibersButton.toolTip = "Recalculate the fibers."
    inputsFormLayout.addRow("Fiber and Grey Matter:", self.updateFibersButton)
    self.updateFibersButton.stateChanged.connect(self.onUpdateFibersClicked)

    self.writeScoresButton = qt.QPushButton("Save the scores")
    self.writeScoresButton.toolTip = "Write the functional scores to disk."
    inputsFormLayout.addRow("Calculate and save scores:", self.writeScoresButton)
    self.writeScoresButton.connect("clicked(bool)", self.onSaveScoresClicked)

    self.labeledAtlasVolumeSelector = slicer.qMRMLNodeComboBox()
    self.labeledAtlasVolumeSelector.nodeTypes = ["vtkMRMLScalarVolumeNode"]
    self.labeledAtlasVolumeSelector.selectNodeUponCreation = True
    self.labeledAtlasVolumeSelector.addEnabled = False
    self.labeledAtlasVolumeSelector.removeEnabled = False
    self.labeledAtlasVolumeSelector.noneEnabled = True
    self.labeledAtlasVolumeSelector.showHidden = False
    self.labeledAtlasVolumeSelector.showChildNodeTypes = False
    self.labeledAtlasVolumeSelector.setMRMLScene(slicer.mrmlScene)
    self.labeledAtlasVolumeSelector.setToolTip("Pick the Atlas-based Segmented Brain Volume.")
    inputsFormLayout.addRow("Labeled Atlas Volume: ", self.labeledAtlasVolumeSelector)
    self.labeledAtlasVolumeSelector.connect("currentNodeChanged(vtkMRMLNode*)", self.onLabeledAtlasVolumeChanged)

    # Generate a rectangular Box for Score Guidance
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
    inputsFormLayout.addRow("Generate the Score Box: ", self.Rec_Box)
    self.Rec_Box.connect("currentNodeChanged(vtkMRMLNode*)", self.OnSelect_Rec_Box) # missing self.OnSelect_label function

    self.skullMaskVolumeSelector = slicer.qMRMLNodeComboBox()
    self.skullMaskVolumeSelector.nodeTypes = ["vtkMRMLScalarVolumeNode"]
    self.skullMaskVolumeSelector.selectNodeUponCreation = False
    self.skullMaskVolumeSelector.addEnabled = False
    self.skullMaskVolumeSelector.removeEnabled = False
    self.skullMaskVolumeSelector.noneEnabled = True
    self.skullMaskVolumeSelector.showHidden = False
    self.skullMaskVolumeSelector.showChildNodeTypes = False
    self.skullMaskVolumeSelector.setMRMLScene(slicer.mrmlScene)
    self.skullMaskVolumeSelector.setToolTip("Pick the skull mask volume.")
    inputsFormLayout.addRow("Skull Mask Volume: ", self.skullMaskVolumeSelector)
    self.skullMaskVolumeSelector.connect("currentNodeChanged(vtkMRMLNode*)", self.onSkullMaskChanged)

    self.tractographyFibreBundleSelector = slicer.qMRMLNodeComboBox()
    self.tractographyFibreBundleSelector.nodeTypes = ["vtkMRMLFiberBundleNode"]
    self.tractographyFibreBundleSelector.selectNodeUponCreation = True
    self.tractographyFibreBundleSelector.addEnabled = False
    self.tractographyFibreBundleSelector.removeEnabled = False
    self.tractographyFibreBundleSelector.noneEnabled = True
    self.tractographyFibreBundleSelector.showHidden = False
    self.tractographyFibreBundleSelector.showChildNodeTypes = False
    self.tractographyFibreBundleSelector.setMRMLScene(slicer.mrmlScene)
    self.tractographyFibreBundleSelector.setToolTip("Pick the full-brain Tractography file.")
    inputsFormLayout.addRow("Tractography Fiber Bundle: ", self.tractographyFibreBundleSelector)
    self.tractographyFibreBundleSelector.connect("currentNodeChanged(vtkMRMLNode*)", self.onFiberBundleChanged)

    self.roiSelector = slicer.qMRMLNodeComboBox()
    self.roiSelector.nodeTypes = ["vtkMRMLAnnotationROINode"]
    self.roiSelector.selectNodeUponCreation = True
    self.roiSelector.addEnabled = True
    self.roiSelector.removeEnabled = True
    self.roiSelector.noneEnabled = True
    self.roiSelector.showHidden = False
    self.roiSelector.showChildNodeTypes = False
    self.roiSelector.setMRMLScene(slicer.mrmlScene)
    self.roiSelector.setToolTip("Pick the ROI to define the region for labeling of fibers.")
    inputsFormLayout.addRow("Labeling ROI: ", self.roiSelector)
    self.roiSelector.connect("currentNodeChanged(vtkMRMLNode*)", self.onRoiChanged)

    self.saveArrayButton = qt.QPushButton("Generate and Save Arrays")
    self.saveArrayButton.toolTip = "Run the algorithm!"
    inputsFormLayout.addRow(self.saveArrayButton)
    self.saveArrayButton.connect('clicked(bool)', self.onSaveArrayClicked)

    self.loadArrayButton = qt.QPushButton("Load Saved Array")
    self.loadArrayButton.toolTip = "Run the algorithm!"
    inputsFormLayout.addRow(self.loadArrayButton)
    self.loadArrayButton.connect('clicked(bool)', self.onLoadArrayClicked)

    # output fiberbundle generator
    self.fiberOutputSelector = slicer.qMRMLNodeComboBox()
    self.fiberOutputSelector.nodeTypes = ["vtkMRMLFiberBundleNode"]
    self.fiberOutputSelector.selectNodeUponCreation = True
    self.fiberOutputSelector.addEnabled = True
    self.fiberOutputSelector.removeEnabled = True
    self.fiberOutputSelector.noneEnabled = True
    self.fiberOutputSelector.showHidden = False
    self.fiberOutputSelector.showChildNodeTypes = False
    self.fiberOutputSelector.setMRMLScene(slicer.mrmlScene)
    self.fiberOutputSelector.setToolTip("Create or pick the output fiber bundle node.")
    inputsFormLayout.addRow("Output fiber bundle: ", self.fiberOutputSelector)
    self.fiberOutputSelector.connect("currentNodeChanged(vtkMRMLNode*)", self.onOutputFiberChanged)

    # Input the Models Hierarchy Node
    self.brainSurfaceModelSelector = slicer.qMRMLNodeComboBox()
    self.brainSurfaceModelSelector.nodeTypes = ["vtkMRMLModelHierarchyNode"]
    self.brainSurfaceModelSelector.selectNodeUponCreation = True
    self.brainSurfaceModelSelector.addEnabled = False
    self.brainSurfaceModelSelector.removeEnabled = False
    self.brainSurfaceModelSelector.noneEnabled = True
    self.brainSurfaceModelSelector.showHidden = False
    self.brainSurfaceModelSelector.showChildNodeTypes = False
    self.brainSurfaceModelSelector.setMRMLScene(slicer.mrmlScene)
    self.brainSurfaceModelSelector.setToolTip("Pick the Hierarchy of Brain Surface Models.")
    inputsFormLayout.addRow("Brain surface models: ", self.brainSurfaceModelSelector)
    self.brainSurfaceModelSelector.connect("currentNodeChanged(vtkMRMLNode*)", self.brainSurfaceModelChanged)

    # Apply Button
    radioLayout = qt.QHBoxLayout()
    buttonGroup = qt.QButtonGroup()
    self.radioButtons=[]
    self.radioButtons.append(qt.QRadioButton('FC_GM'))
    self.radioButtons.append(qt.QRadioButton('SC_GM'))
    self.radioButtons.append(qt.QRadioButton('Opacity_reset'))
    buttonGroup.addButton(self.radisButtons[0],0)
    buttonGroup.addButton(self.radioButtons[1],1)
    buttonGroup.addButton(self.radioButtons[2],2)
    radioLayout.addWidget(self.radioButtons[0])
    radioLayout.addWidget(self.radioButtons[1])
    radioLayout.addWidget(self.radioButtons[2])
    inputsFormLayout.addRow("Select the Grey Matter Score to paint: ", radioLayout)
    self.radioButtons[0].connect('clicked()', self.functionConnectivityGreyMatterClicked)
    self.radioButtons[1].connect('clicked()', self.structuralConnectivityGreyMatterClicked)
    self.radioButtons[2].connect('clicked()', self.resetModelOpacities)

    self.runOptimizationButton = qt.QPushButton("Run the Optimization")
    self.runOptimizationButton.toolTip = "Run the algorithm!"
    inputsFormLayout.addRow("Run optimization: ", self.runOptimizationButton)
    self.runOptimizationButton.connect("clicked(bool)", self.onRunOptimizationClicked)

  def onModelPointsSelected(self):
    # Create cylinder and sphere from selected fiducials
    cylinderTop = [0,0,0]
    cylinderBottom = [0,0,0]
    sphereOrigin = [0,0,0]
    markupNode = self.modelPointsSelector.currentNode()
    if markupNode.GetNumberOfMarkups() < 3:
      print("Markup node is invalid. Needs at least 3 points.")
      return

    markupNode.GetNthFiducialPosition(0, cylinderTop)
    markupNode.GetNthFiducialPosition(1, cylinderBottom)
    markupNode.GetNthFiducialPosition(2, sphereOrigin)

    self.logic.computeModels(cylinderTop, cylinderBottom, self.qSlider_CylinderRadius.value, self.qSlider_CylinderHeight.value, sphereOrigin, self.qSlider_SphereRadius.value)

  def cleanup(self):
    pass

  def onSliderChange(self):
    self.logic.recomputePath(self.qSlider_SphereRadius.value, self.qSlider_CylinderRadius.value, self.qSlider_CylinderHeight.value)

  def onUpdateFibersClicked(self):
    if self.tractographyFibreBundleSelector.currentNode() is None:
      logging.error("Must select a tractography fiber bundle.")
      return

    if self.labeledAtlasVolumeSelector.currentNode() is None:
      logging.error("Must select labeled atlas.")
      return

    if self.skullMaskVolumeSelector.currentNode() is None:
      logging.error("Skull mask volume required.")
      return

    outFibers, Label_l = self.logic.updateFibers(self.labeledAtlasVolumeSelector.currentNode(),
                                                 self.tractographyFibreBundleSelector.currentNode(),
                                                 self.skullMaskVolumeSelector.currentNode(),
                                                 self.Aparc,
                                                 self.Aparc2,
                                                 self.arrays_list,
                                                 self.mem_points,
                                                 self.mem_array,
                                                 self.labelLUT)

    if self.fiberOutputSelector.currentNode() is not None:
      self.fiberOutputSelector.currentNode().SetAndObservePolyData(outFibers)

    if len(Label_l) > 0:
       self.changeModelVisibility(Label_l)

  # For loading input label, two main purposes generating grid and model handling
  def onLabeledAtlasVolumeChanged(self):
    pass

  # For generating grid image
  def OnSelect_Rec_Box(self):
    labeledAtlasNode = self.labeledAtlasVolumeSelector.currentNode()
    Rec_Box_node = self.Rec_Box.currentNode()

    labeledAtlasRASToIJK = vtk.vtkMatrix4x4()
    self.labeledAtlasVolumeSelector.currentNode().GetRASToIJKMatrix(labeledAtlasRASToIJK)
    Rec_Box_node.SetRASToIJKMatrix(labeledAtlasRASToIJK)
    origin = labeledAtlasNode.GetOrigin()
    Rec_Box_node.SetOrigin(origin)

    fileName = qt.QFileDialog.getOpenFileName()
    fileHandle = open(fileName,'rb')
    csv_f = csv.reader(fileHandle)
    ar = np.genfromtxt (fileName, delimiter=",")
    spaceCount = ar.shape[0]
    c = ar.shape[1]
    table = np.zeros((spaceCount,c))
    i = 0
    for row in csv_f:
      j=0
      for field in row:
        table[i,j] = float(field)
        j+=1
      i+=1

    skullTable = table
    trans = vtk.vtkTransform()
    trans.Identity()
    trans.PreMultiply()
    trans.SetMatrix(labeledAtlasRASToIJK)
    pIJK = [0,0,0]
    pIJK_far = [0,0,0]
    pt = [0,0,0]
    ftarget = np.array([int(table[0,0]),int(table[0,1]),int(table[0,2])])
    k=0
    for s in range(1, spaceCount):
      fstart_p = np.array([int(table[s,0]),int(table[s,1]),int(table[s,2])])
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
        if self.skullMaskVolumeSelector.currentNode().GetImageData().GetScalarComponentAsDouble(pt[0],pt[1],pt[2],0)>0:
          k+=1
          imageRec_Box.GetOutput().SetScalarComponentFromDouble(int(pIJK_far[0]),int(pIJK_far[1]),int(pIJK_far[2]),0,int(table[s,9]+1))
          skullTable[k,0:3] = point_p_far
          break

    Rec_Box_node.SetAndObserveImageData(imageRec_Box.GetOutput())
    os.chdir(qt.QFileDialog.getExistingDirectory())
    np.savetxt("skullTable.csv",skullTable,delimiter=",")

  def onSkullMaskChanged(self):
    skullMaskNode = self.skullMaskVolumeSelector.currentNode()
    labelDims = skullMaskNode.GetImageData().GetDimensions()
    space = 5
    a = int(labelDims[0]/space)
    b = int(labelDims[1]/space)
    c = int(labelDims[2]/space)

  def onFiberBundleChanged(self):
    tractographyNode = self.tractographyFibreBundleSelector.currentNode()
    tractographyPolyData = tractographyNode.GetPolyData()

    self.arrays_list =[]
    for i in range(1,tractographyPolyData.GetPointData().GetNumberOfArrays()): #Array 0 is the tensors
      self.arrays_list.append(tractographyPolyData.GetPointData().GetArray(i))

  def onSaveArrayClicked(self):
    saveDirectory = qt.QFileDialog.getExistingDirectory()
    if self.tractographyFibreBundleSelector.currentNode() is None or self.skullMaskVolumeSelector.currentNode() is None:
      logging.error("Tractography and skull mask required.")
      return

    self.mem_array, self.mem_points, self.labelLUT = self.logic.calculateArrayData(saveDirectory, self.tractographyFibreBundleSelector.currentNode(), self.skullMaskVolumeSelector.currentNode())

  def onLoadArrayClicked(self):
    loadDirectoryName = qt.QFileDialog.getExistingDirectory()

    self.mem_array = np.load(loadDirectoryName + '/Mem_array.npy')
    self.mem_points = np.load(loadDirectoryName + '/Mem_points.npy')
    f = open(loadDirectoryName + '/Label_LUT', 'rb')
    self.labelLUT = map(int, f.read().split())
    f = open(loadDirectoryName + '/FC_GM', 'rb')
    self.FC_GM = map(float, f.read().split())
    f = open(loadDirectoryName + '/SC_GM', 'rb')
    self.SC_GM = map(float, f.read().split())
    f = open(loadDirectoryName + '/Aparc', 'rb')  # All of the Labels of the Atlas regions
    self.Aparc = map(float, f.read().split())
    f = open(loadDirectoryName + '/Aparc2', 'rb')  # FC Absolute value of the Atlas regions
    self.Aparc2 = map(float, f.read().split())
    f = open(loadDirectoryName + '/Aparc3', 'rb')  # FC Absolute value of the Atlas regions
    self.Aparc3 = map(float, f.read().split())

  def onOutputFiberChanged(self):
    outputFiberNode = self.fiberOutputSelector.currentNode()
    if self.tractographyFibreBundleSelector.currentNode() is not None:
      outputFiberNode.SetAndObserveTransformNodeID(self.tractographyFibreBundleSelector.currentNode().GetTransformNodeID())
    outputFiberNode.CreateDefaultDisplayNodes()
    outputFiberNode.CreateDefaultStorageNode()
    outputFiberNode.SetName('Selected_Fibers')

  def brainSurfaceModelChanged(self):
    brainSurfaceNode = self.brainSurfaceModelSelector.currentNode()
    if brainSurfaceNode is not None:
      self.radioButtons.enabled = True
    else:
      self.radioButtons.enabled = False

    modelCount = brainSurfaceNode.GetNumberOfChildrenNodes()
    self.node_names = []
    for t in range(0, modelCount):
      child = brainSurfaceNode.GetNthChildNode(t)
      node = child.GetModelNode()
      name = node.GetName()
      name = name[6:name.__len__()]# Remove first part of the name which is usually Model_
      self.node_names.append(int(name[0:name.index('_')]))

  def functionConnectivityGreyMatterClicked(self):
    brainSurfaceNode = self.brainSurfaceModelSelector.currentNode()
    modelCount = brainSurfaceNode.GetNumberOfChildrenNodes()

    for t in range(0,modelCount):
      child = brainSurfaceNode.GetNthChildNode(t)
      modelNode = child.GetModelNode()
      name = modelNode.GetName()
      displayNode = modelNode.GetDisplayNode()
      displayNode.SetOpacity(self.FC_GM[t])
      if self.FC_GM[t]<0.5:
        displayNode.SetColor(2*self.FC_GM[t],1,0)
      else:
        displayNode.SetColor(1,2-2*self.FC_GM[t],0)

  def structuralConnectivityGreyMatterClicked(self):
    brainSurfaceNode = self.brainSurfaceModelSelector.currentNode()
    modelCount = brainSurfaceNode.GetNumberOfChildrenNodes()

    for t in range(0, modelCount):
      child = brainSurfaceNode.GetNthChildNode(t)
      node = child.GetModelNode()
      name = node.GetName()
      d_node =node.GetDisplayNode()
      d_node.SetOpacity(self.SC_GM[t])
      if self.SC_GM[t]<0.5:
        d_node.SetColor(2*self.SC_GM[t],0,1)
      else:
        d_node.SetColor(1,0,2-2*self.SC_GM[t])

  def resetModelOpacities(self):
    brainSurfaceNode = self.brainSurfaceModelSelector.currentNode()
    modelCount = brainSurfaceNode.GetNumberOfChildrenNodes()

    for t in range(0, modelCount):
      child = brainSurfaceNode.GetNthChildNode(t)
      node = child.GetModelNode()
      name = node.GetName()
      d_node =node.GetDisplayNode()
      d_node.SetOpacity(1)

  def onRunOptimizationClicked(self):
    if self.labeledAtlasVolumeSelector.currentNode() is None:
      logging.error("Labeled atlas required.")
      return

    bounds = [0,0,0,0,0,0]
    if self.roiSelector.currentNode() is None:
      if self.labeledAtlasVolumeSelector.currentNode() is None:
        logging.error("Input volume required to run the algorithm.")
        return
      else:
        bounds = self.labeledAtlasVolumeSelector.currentNode().GetImageData().GetBounds()
    else:
      self.roiSelector.currentNode().GetRASBounds(bounds)

    self.logic.runOptimization(self.labeledAtlasVolumeSelector.currentNode(), self.tractographyFibreBundleSelector.currentNode(), self.mem_points, self.mem_array, self.labelLUT, self.arrays_list, bounds)

  def changeModelVisibility(self, labels):
    if self.brainSurfaceModelSelector.currentNode() is None:
      print "Brain models required before visibility can be changed."
      return

    for s in range(0, len(self.node_names)):
      child = self.brainSurfaceModelSelector.currentNode().GetNthChildNode(s)
      modelNode = child.GetModelNode()
      modelNode.SetDisplayVisibility(0)
    for r in range(0, len(labels)):
      if labels[r] in self.node_names:
        child = self.brainSurfaceModelSelector.currentNode().GetNthChildNode(self.node_names.index(labels[r]))
        modelNode = child.GetModelNode()
        displayNode = modelNode.GetDisplayNode()
        modelNode.SetDisplayVisibility(1)
        displayNode.SetOpacity(0.75)

# NeuroPathLogic
class NeuroPathLogic(ScriptedLoadableModuleLogic):
  def __init__(self):
    self.score = None
    self.cylinderTop = None
    self.cylinderBottom = None
    self.sphereOrigin = None
    self.cylinderRadius = None
    self.sphereRadius = None
    self.cylinderHeight = None

    self.cylinderModelNode = None
    self.sphereModelNode = None

    self.cylinderPolyDataFilter = None

  def runOptimization(self, labeledAtlasNode, tractographyNode, mem_points, mem_array, label_LUT, arrays_list, bounds):
    if self.cylinderTop is None:
      print "Must set cylinderSource before calculating."
      return
    if tractographyNode is None:
      print "Must set tractography node before running algorithm"
      return

    xMin = int(bounds[0]/10)
    xMax = int(bounds[1]/10)
    yMin = int(bounds[2]/10)
    yMax = int(bounds[3]/10)
    zMin = int(bounds[4]/10)
    zMax = int(bounds[5]/10)
    sampleCount = (xMax-xMin+1)*(yMax-yMin+1)+2*((yMax-yMin+1)*(zMax-zMin+1)+(xMax-xMin+1)*(zMax-zMin+1))

    table = np.zeros((sampleCount+1,10))
    table[0,4] = xMin
    table[0,5] = xMax
    table[0,6] = yMin
    table[0,7] = yMax
    table[0,8] = zMin
    table[0,9] = zMax
    k = 1

    for i in range(0,(xMax-xMin)+1):
     for j in range(0,(yMax-yMin)+1):
       table[k,0:3] = [(xMin+i)*10,(yMin+j)*10,zMax*10]
       k+=1
      
    for i in range(0,(yMax-yMin)+1):
     for j in range(0,(zMax-zMin)+1):
       table[k,0:3] = [xMin*10,(yMin+i)*10,(zMin+j)*10]
       k+=1
    		
    for i in range(0,(yMax-yMin)+1):
     for j in range(0,(zMax-zMin)+1):
       table[k,0:3] = [xMax*10,(yMin+i)*10,(zMin+j)*10]
       k+=1
              
    for i in range(0,(xMax-xMin)+1):
     for j in range(0,(zMax-zMin)+1):
       table[k,0:3] = [(xMin+i)*10,yMin*10,(zMin+j)*10]
       k+=1
        
    for i in range(0,(xMax-xMin)+1):
     for j in range(0,(zMax-zMin)+1):
       table[k,0:3] = [(xMin+i)*10,yMax*10,(zMin+j)*10]
       k+=1

    fs = self.sphereOrigin
    table[0,0:3] = self.cylinderTop
    for m in range(1, sampleCount+1):
      self.score = [0,0,0,0,0,0]
      cylinderBottom = np.array(table[m,0:3])
      cylinderAxis = (cylinderBottom-self.cylinderTop)
      cylinderHeightSearch = np.sqrt(sum(cylinderAxis*cylinderAxis))
      normalizedCylinderAxis = cylinderAxis/cylinderHeightSearch
      cylinderCenter = self.cylinderTop + (cylinderHeightSearch/2)*normalizedCylinderAxis
      cylinderEnd = cylinderCenter+(cylinderHeightSearch/2)*normalizedCylinderAxis
      cylinderStart = cylinderCenter-(cylinderHeightSearch/2)*normalizedCylinderAxis
      cylinderSource = vtk.vtkCylinderSource()
      cylinderSource.SetCenter(cylinderCenter)
      plane = vtk.vtkPlane()
      plane.SetNormal(normalizedCylinderAxis)
      plane.SetOrigin(cylinderEnd)
      upAxis = np.array([0,1,0])
      cylinderSource.SetRadius(self.cylinderRadius)
      cylinderSource.SetHeight(cylinderHeightSearch)

      sphereSource = vtk.vtkSphereSource()
      sphereSource.SetCenter(self.sphereOrigin)
      sphereSource.SetRadius(self.sphereRadius)

      cylinderTransform = vtk.vtkTransform()
      cylinderTransform.PostMultiply()
      cylinderTransform.Translate(-cylinderCenter[0],-cylinderCenter[1],-cylinderCenter[2])
      cylinderTransform.RotateWXYZ(-1*np.degrees(self.py_ang(cylinderAxis,upAxis)),np.cross(cylinderAxis,upAxis))
      cylinderTransform.Translate(cylinderCenter[0],cylinderCenter[1],cylinderCenter[2])

      polyDataFilter = vtk.vtkTransformPolyDataFilter()
      polyDataFilter.SetInputConnection(cylinderSource.GetOutputPort())
      polyDataFilter.SetTransform(cylinderTransform)
      polyDataFilter.Update()

      pt = [0,0,0]
      Label_num_a =[]
      Label_l = []
      pointCount = int(cylinderHeightSearch)
      cylinderPoints = np.zeros((pointCount,3))
      pointSource = vtk.vtkPointSource()
      pointSource.SetNumberOfPoints(10000)## CHANGED FROM 1000000 to 10000 to speed up the process
      pointSource.SetRadius(cylinderHeightSearch/2+20)
      pointSource.SetCenter(cylinderCenter)
      pointSource.Update()

      pointLocator = vtk.vtkPointLocator()
      pointLocator.SetDataSet(pointSource.GetOutput())
      pointLocator.AutomaticOn()
      pointLocator.SetNumberOfPointsPerBucket(1)
      pointLocator.BuildLocator()
      result = vtk.vtkIdList()
      pIJK = [0,0,0]

      labelDimensions = imageRec_Main.GetOutput().GetDimensions()
      for i in range(0, pointCount):
        cylinderPoints[i,:] = cylinderStart+(i)*normalizedCylinderAxis
        pointLocator.FindPointsWithinRadius(self.cylinderRadius, cylinderStart+(i)*normalizedCylinderAxis,result)
        for j in range(0,result.GetNumberOfIds()):
          point = pointSource.GetOutput().GetPoint(result.GetId(j))
          cylinderTransform.TransformPoint(point,pIJK) #RAS to IJK transformation
          pt[0]= int(pIJK[0])
          pt[1]= int(pIJK[1])
          pt[2]= int(pIJK[2])
          if pt[0] < 0 or pt[1] < 0 or pt[2] < 0 or pt[0] >= labelDimensions[0] or pt[1] >= labelDimensions[1] or pt[2] >= labelDimensions[2]:
            continue
          inPtr = imageRec_Main.GetOutput().GetScalarComponentAsDouble(pt[0],pt[1],pt[2],0)
          Label_num_a.append(inPtr)
          inPtr2 = labeledAtlasNode.GetImageData().GetScalarComponentAsDouble(pt[0],pt[1],pt[2],0)
          Label_l.append(inPtr2)
          self.score[0] += int(self.Aparc2[self.Aparc.index(int(inPtr2))])
          self.score[1] += int(self.Aparc3[self.Aparc.index(int(inPtr2))])

      Label_num_a = list(set(Label_num_a))
      length = Label_num_a.__len__()
      addLine = 0*(mem_array[:,0])
      for s in range(length):
        if Label_num_a[s] in label_LUT:
          exist = mem_array[:,label_LUT.index(Label_num_a[s])]
          addLine += exist

      addLine = addLine>0
      numNewCells = int(sum(addLine))

      address = np.where(addLine==1)[:][0]
      tractographyPolyData = tractographyNode.GetPolyData()
      oldTensors = tractographyPolyData.GetPointData().GetTensors()

      if numNewCells>0:
        ptId = 0
        address1 = np.ndarray.tolist(address)
        for inCellId in address1:
          pointCountInCell = int(mem_points[inCellId])
          if oldTensors:
            for l in range(0, tractographyPolyData.GetPointData().GetNumberOfArrays()):
              self.score[3] += int(arrays_list[l].GetTuple(ptId)[0])
          ptId += pointCountInCell

      Label_l = list(set(Label_l)) # Remove duplicate list items
      table[m,3:(3+len(self.score))] = self.score
      table[m,9] = cylinderHeightSearch

    saveDirectory = qt.QFileDialog.getExistingDirectory()
    np.savetxt(saveDirectory + "/table.csv",table,delimiter=",")

  # Run the computation for path refresh
  def recomputePath(self):
    cylinderCenter = ((self.cylinderTop + self.cylinderBottom)/2)
    cylinderSource = vtk.vtkCylinderSource()
    cylinderSource.SetCenter(cylinderCenter)
    cylinderAxis = (self.cylinderTop-self.cylinderBottom)
    normalizedCylinderAxis = cylinderAxis/np.sqrt(sum(cylinderAxis*cylinderAxis))
    cylinderEnd = cylinderCenter+(self.cylinderHeight/2)*normalizedCylinderAxis
    plane = vtk.vtkPlane()
    plane.SetNormal(normalizedCylinderAxis)
    plane.SetOrigin(cylinderEnd)
    upAxis = np.array([0,1,0])
    sphereSource = vtk.vtkSphereSource()
    sphereSource.SetCenter(self.sphereOrigin)
    cylinderSource.SetRadius(self.cylinderRadius)
    cylinderSource.SetHeight(self.cylinderHeight)
    sphereSource.SetRadius(self.sphereRadius)

    cylinderTransform = vtk.vtkTransform()
    cylinderTransform.PostMultiply()
    cylinderTransform.Translate(-cylinderCenter[0],-cylinderCenter[1],-cylinderCenter[2])
    cylinderTransform.RotateWXYZ(-1*np.degrees(self.py_ang(cylinderAxis,upAxis)),np.cross(cylinderAxis,upAxis))
    cylinderTransform.Translate(cylinderCenter[0],cylinderCenter[1],cylinderCenter[2])

    self.cylinderPolyDataFilter.SetInputConnection(cylinderSource.GetOutputPort())
    self.cylinderPolyDataFilter.SetTransform(cylinderTransform)
    self.cylinderPolyDataFilter.Update()

  # Run the computation for collision detection and score computation
  def updateFibers(self, labeledAtlasNode, tractographyNode, skullMaskNode, Aparc, Aparc2, arrays_list, mem_points, mem_array, labelLUT):
    self.score=[]
    self.score.append(0) # This first column is reserved
    GM_Counter = [0]*len(Aparc) # This counter is for number of points inside of access path from each label region

    cylinderCenter = ((self.cylinderTop + self.cylinderBottom) / 2)
    cylinderSource = vtk.vtkCylinderSource()
    cylinderSource.SetCenter(cylinderCenter)
    cylinderAxis = (self.cylinderTop - self.cylinderBottom)
    normalizedCylinderAxis = cylinderAxis / np.sqrt(sum(cylinderAxis * cylinderAxis))
    cylinderStart = cylinderCenter-(self.cylinderHeight/2)*normalizedCylinderAxis
    pt = [0,0,0]
    Label_num_a =[]

    upAxis = np.array([0, 1, 0])
    cylinderTransform = vtk.vtkTransform()
    cylinderTransform.PostMultiply()
    cylinderTransform.Translate(-cylinderCenter[0], -cylinderCenter[1], -cylinderCenter[2])
    cylinderTransform.RotateWXYZ(-1 * np.degrees(self.py_ang(cylinderAxis, upAxis)), np.cross(cylinderAxis, upAxis))
    cylinderTransform.Translate(cylinderCenter[0], cylinderCenter[1], cylinderCenter[2])

    pointSource = vtk.vtkPointSource()
    pointSource.SetNumberOfPoints(1000000)
    pointSource.SetRadius(self.cylinderHeight2+20)
    pointSource.SetCenter(cylinderCenter)
    pointSource.Update()
    pointLocator = vtk.vtkPointLocator()
    pointLocator.SetDataSet(pointSource.GetOutput())
    pointLocator.AutomaticOn()
    pointLocator.SetNumberOfPointsPerBucket(1)
    pointLocator.BuildLocator()

    sphericalPointSource = vtk.vtkPointSource()
    sphericalPointSource.SetNumberOfPoints(500000)
    sphericalPointSource.SetRadius(self.sphereRadius)
    sphericalPointSource.SetCenter(self.sphereOrigin)
    sphericalPointSource.Update()
    sphericalPointLocator = vtk.vtkPointLocator()
    sphericalPointLocator.SetDataSet(sphericalPointSource.GetOutput())
    sphericalPointLocator.AutomaticOn()
    sphericalPointLocator.SetNumberOfPointsPerBucket(1)
    sphericalPointLocator.BuildLocator()

    result = vtk.vtkIdList()
    pIJK = [0,0,0]
    Label_l = []
    labelDimensions = skullMaskNode.GetImageData().GetDimensions()
    pointCount = int(self.cylinderHeight)
    cylinderPoints = np.zeros((pointCount, 3))
    for i in range(0,pointCount):
      cylinderPoints[i,:] = cylinderStart+(i)*normalizedCylinderAxis
      pointLocator.FindPointsWithinRadius(self.cylinderRadius, cylinderStart+(i)*normalizedCylinderAxis,result)
      for j in range(0,result.GetNumberOfIds()):
        point = pointSource.GetOutput().GetPoint(result.GetId(j))
        cylinderTransform.TransformPoint(point,pIJK) #RAS to IJK transformation
        pt[0]= int(pIJK[0])
        pt[1]= int(pIJK[1])
        pt[2]= int(pIJK[2])
        if pt[0] < 0 or pt[1] < 0 or pt[2] < 0 or pt[0] >= labelDimensions[0] or pt[1] >= labelDimensions[1] or pt[2] >= labelDimensions[2]:
          continue
        inPtr = skullMaskNode.GetImageData().GetScalarComponentAsDouble(pt[0],pt[1],pt[2],0)
        Label_num_a.append(inPtr)
        inPtr2 = labeledAtlasNode.GetImageData().GetScalarComponentAsDouble(pt[0],pt[1],pt[2],0)
        Label_l.append(inPtr2)
        self.score[0] += Aparc2[Aparc.index(int(inPtr2))]
        GM_Counter[Aparc.index(int(inPtr2))] += 1

    sphericalPointLocator.FindPointsWithinRadius(self.sphereRadius, self.sphereOrigin, result)
    for j in range(0,result.GetNumberOfIds()): 
      point = sphericalPointSource.GetOutput().GetPoint(result.GetId(j))
      cylinderTransform.TransformPoint(point,pIJK) #RAS to IJK transformation
      pt[0]= int(pIJK[0])
      pt[1]= int(pIJK[1])
      pt[2]= int(pIJK[2])
      if pt[0] < 0 or pt[1] < 0 or pt[2] < 0 or pt[0] >= labelDimensions[0] or pt[1] >= labelDimensions[1] or pt[2] >= labelDimensions[2]:
        continue
      inPtr = skullMaskNode.GetImageData().GetScalarComponentAsDouble(pt[0],pt[1],pt[2],0)
      Label_num_a.append(inPtr)
      inPtr2 = labeledAtlasNode.GetImageData().GetScalarComponentAsDouble(pt[0],pt[1],pt[2],0)
      Label_l.append(inPtr2)
      self.score[0] += Aparc2[Aparc.index(int(inPtr2))]
      GM_Counter[Aparc.index(int(inPtr2))] += 1

    Label_num_a = list(set(Label_num_a))# To remove redundancies

    length = len(Label_num_a)
    addLine = 0*(mem_array[:,0])
    for s in range(length):
      if Label_num_a[s] in labelLUT:
        exist = mem_array[:,labelLUT.index(Label_num_a[s])]
        addLine += exist
    addLine = addLine>0
    numNewCells = int(sum(addLine))

    address = np.where(addLine==1)[:][0]
    numNewPts = int(sum(mem_points[address]))

    outFibers = None
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
      for j in range(0, tractographyNode.GetPolyData().GetPointData().GetNumberOfArrays()-1):
        newArrays.append(vtk.vtkFloatArray())
        newArrays[j].SetName(arrays_list[j].GetName())
        newArrays[j].SetNumberOfComponents(1)
        newArrays[j].Allocate(numNewPts)
        outFibers.GetPointData().AddArray(newArrays[j])
        self.score.append(0)

      ptId = 0
      tractographyNode.GetPolyData().GetLines().InitTraversal()
      address1 = np.ndarray.tolist(address)
      killLoopFlag = False
      for inCellId in address1:
        pointCount = int(mem_points[inCellId])
        outFibersCellArray.InsertNextCell(pointCount,range(ptId,ptId+pointCount))
        st_pt = int(sum(mem_points[0:(inCellId)]))
        points.InsertPoints(ptId, pointCount, st_pt, tractographyNode.GetPolyData().GetPoints())
        if (tractographyNode.GetPolyData().GetPointData().GetTensors()):
          newTensors.InsertTuples(ptId,pointCount,st_pt,tractographyNode.GetPolyData().GetPointData().GetTensors())
          for l in range(0,tractographyNode.GetPolyData().GetPointData().GetNumberOfArrays()-1):
            newArrays[l].InsertTuples(ptId,pointCount,st_pt, arrays_list[l])
            self.score[l+1] += arrays_list[l].GetTuple(ptId)[0]
        ptId += pointCount
        if killLoopFlag:
          break
    return outFibers, list(set(Label_l))
 
  def py_ang(self, v1, v2):
   """ Returns the angle in radians between vectors 'v1' and 'v2' """
   return np.arctan2(la.norm(np.cross(v1, v2)), np.dot(v1, v2))

  # Calculating mem_array, mem_points and label_LUT
  def calculateArrayData(self, tractographyNode, skullMaskNode):
    Label_A_RASToIJK = vtk.vtkMatrix4x4()
    skullMaskNode.GetRASToIJKMatrix(Label_A_RASToIJK)
    trans = vtk.vtkTransform()
    trans.Identity()
    trans.PreMultiply()
    trans.SetMatrix(Label_A_RASToIJK)
    numLines = tractographyNode.GetPolyData().GetLines().GetNumberOfCells()
    pts = vtk.vtkIdList()
    _mem_array = np.full((numLines,a*b*c), False, dtype=bool)
    _mem_points = np.zeros(numLines)
    _label_LUT = []

    tractographyNode.GetPolyData().GetLines().InitTraversal()

    labelDimensions = skullMaskNode.GetImageData().GetDimensions()
    for inCellId in range(0,numLines):
      tractographyNode.GetPolyData().GetLines().GetNextCell(pts)
      npts = int(pts.GetNumberOfIds())
      pIJK = [0,0,0]
      pt = [0,0,0]

      if npts >= 19:
        for j in range(0,npts):
          p = tractographyNode.GetPolyData().GetPoints().GetPoint(pts.GetId(j))
          trans.TransformPoint(p,pIJK)
          pt[0]= int(pIJK[0])
          pt[1]= int(pIJK[1])
          pt[2]= int(pIJK[2])
          if pt[0] < 0 or pt[1] < 0 or pt[2] < 0 or pt[0] >= labelDimensions[0] or pt[1] >= labelDimensions[1] or pt[2] >= labelDimensions[2]:
            continue
          inPtr = skullMaskNode.GetImageData().GetScalarComponentAsDouble(pt[0],pt[1],pt[2],0)
          if inPtr<1:
            continue
          if inPtr in _label_LUT:
            _mem_array[inCellId,_label_LUT.index(inPtr)] = True
          else:
            _label_LUT.append(inPtr)
            _mem_array[inCellId,_label_LUT.index(inPtr)] = True
        _mem_points[inCellId] = npts
	
    os.chdir(qt.QFileDialog.getExistingDirectory())
    np.save('mem_array.npy', _mem_array[:,0:len(_label_LUT)])
    np.save('mem_points.npy',_mem_points)
    f = open("_label_LUT", "w")
    for item in _label_LUT:
      f.write("%s " % int(item))
    f.close()

    return _mem_array, _mem_points, _label_LUT

  def computeModels(self, cylinderTop, cylinderBottom, cylinderRadius, cylinderHeight, sphereCenter, sphereRadius):
    self.cylinderTop = np.array(cylinderTop)
    self.cylinderBottom = np.array(cylinderBottom)
    self.sphereOrigin = np.array(sphereCenter)
    self.cylinderRadius = cylinderRadius
    self.sphereRadius = sphereRadius
    self.cylinderHeight = cylinderHeight
    cylinderAxis = self.cylinderTop - self.cylinderBottom
    upAxis = np.array([0, 1, 0])
    cylinderCentre = ((self.cylinderTop + self.cylinderBottom) / 2)

    sphere = vtk.vtkSphereSource()
    sphere.SetCenter(self.sphereOrigin)
    sphere.SetThetaResolution(30)
    sphere.SetPhiResolution(30)
    sphere.SetRadius(self.sphereRadius)
    clip = vtk.vtkClipPolyData()
    clip.SetValue(0)
    clip.GenerateClippedOutputOn()
    clip.SetInputConnection(sphere.GetOutputPort())

    plane = vtk.vtkPlane()
    normalizedCylinderAxis = cylinderAxis / np.sqrt(sum(cylinderAxis * cylinderAxis))
    planeOrigin = cylinderCentre + (self.cylinderHeight / 2) * normalizedCylinderAxis
    plane.SetNormal(normalizedCylinderAxis)
    plane.SetOrigin(planeOrigin)
    clip.SetClipFunction(plane)

    cylinderSource = vtk.vtkCylinderSource()
    cylinderSource.SetRadius(self.cylinderRadius)
    cylinderSource.SetHeight(self.cylinderHeight)
    cylinderSource.SetCenter(cylinderCentre)
    cylinderSource.SetResolution(100)
    cylinderSource.CappingOff()

    cylinderTransform = vtk.vtkTransform()
    cylinderTransform.PostMultiply()
    cylinderTransform.Translate(-cylinderCentre[0], -cylinderCentre[1], -cylinderCentre[2])
    cylinderTransform.RotateWXYZ(-1 * np.degrees(self.logic.py_ang(cylinderAxis, upAxis)), np.cross(cylinderAxis, upAxis))
    cylinderTransform.Translate(cylinderCentre[0], cylinderCentre[1], cylinderCentre[2])
    self.cylinderPolyDataFilter = vtk.vtkTransformPolyDataFilter()
    self.cylinderPolyDataFilter.SetInputConnection(cylinderSource.GetOutputPort())
    self.cylinderPolyDataFilter.SetTransform(cylinderTransform)
    self.cylinderPolyDataFilter.Update()

    modelsLogic = slicer.modules.models.logic()
    self.cylinderModelNode = modelsLogic.AddModel(self.cylinderPolyDataFilter.GetOutputPort())
    self.cylinderModelNode.GetDisplayNode().SetSliceIntersectionVisibility(True)
    self.cylinderModelNode.GetDisplayNode().SetSliceIntersectionThickness(2)
    self.cylinderModelNode.GetDisplayNode().SetColor(0, 0, 1)
    self.cylinderModelNode.SetName('Cylinder')
    self.cylinderModelNode(self.cylinderModelNode)

    self.sphereModelNode = modelsLogic.AddModel(clip.GetOutputPort())
    self.sphereModelNode.GetDisplayNode().SetSliceIntersectionVisibility(True)
    self.sphereModelNode.GetDisplayNode().SetSliceIntersectionThickness(2)
    self.sphereModelNode.GetDisplayNode().SetColor(0, 1, 0)
    self.sphereModelNode.SetName('Sphere')
 
class NeuroPathTest(ScriptedLoadableModuleTest):
  def setUp(self):
    slicer.mrmlScene.Clear(0)
 
  def runTest(self):
    self.setUp()
    self.test_NeuroPath1()
 
  def test_NeuroPath1(self):
    self.delayDisplay("Starting the test")
    self.delayDisplay('Test passed!')
