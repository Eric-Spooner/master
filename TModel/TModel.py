import os
import unittest
import vtk, qt, ctk, slicer
from slicer.ScriptedLoadableModule import *
import logging


#
# TModel
#

class TModel(ScriptedLoadableModule):
    """Uses ScriptedLoadableModule base class, available at:
    https://github.com/Slicer/Slicer/blob/master/Base/Python/slicer/ScriptedLoadableModule.py
    """

    def __init__(self, parent):
        ScriptedLoadableModule.__init__(self, parent)
        self.parent.title = "TModel"  # TODO make this more human readable by adding spaces
        self.parent.categories = ["Examples"]
        self.parent.dependencies = []
        self.parent.contributors = ["John Doe (AnyWare Corp.)"]  # replace with "Firstname Lastname (Organization)"
        self.parent.helpText = """
This is an example of scripted loadable module bundled in an extension.
It performs a simple thresholding on the input volume and optionally captures a screenshot.
"""
        self.parent.helpText += self.getDefaultModuleDocumentationLink()
        self.parent.acknowledgementText = """
This file was originally developed by Jean-Christophe Fillion-Robin, Kitware Inc.
and Steve Pieper, Isomics, Inc. and was partially funded by NIH grant 3P41RR013218-12S1.
"""  # replace with organization, grant and thanks.


#
# TModelWidget
#

class TModelWidget(ScriptedLoadableModuleWidget):
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
        # input volume selector
        #
        self.inputSelector = slicer.qMRMLNodeComboBox()
        self.inputSelector.nodeTypes = ["vtkMRMLScalarVolumeNode"]
        self.inputSelector.selectNodeUponCreation = True
        self.inputSelector.addEnabled = False
        self.inputSelector.removeEnabled = False
        self.inputSelector.noneEnabled = False
        self.inputSelector.showHidden = False
        self.inputSelector.showChildNodeTypes = False
        self.inputSelector.setMRMLScene(slicer.mrmlScene)
        self.inputSelector.setToolTip("Pick the input volume to the algorithm.")
        parametersFormLayout.addRow("Input Volume: ", self.inputSelector)

        #
        # armature selector
        #
        self.armatureSelector = slicer.qMRMLNodeComboBox()
        self.armatureSelector.nodeTypes = ["vtkMRMLModelNode"]
        self.armatureSelector.selectNodeUponCreation = True
        self.armatureSelector.addEnabled = False
        self.armatureSelector.removeEnabled = False
        self.armatureSelector.noneEnabled = False
        self.armatureSelector.showHidden = False
        self.armatureSelector.showChildNodeTypes = False
        self.armatureSelector.setMRMLScene(slicer.mrmlScene)
        self.armatureSelector.setToolTip("Pick the input armature to the algorithm.")
        parametersFormLayout.addRow("Input Armature: ", self.armatureSelector)

        # output volume selector
        #
        self.outputSelector = slicer.qMRMLNodeComboBox()
        self.outputSelector.nodeTypes = ["vtkMRMLScalarVolumeNode"]
        self.outputSelector.selectNodeUponCreation = True
        self.outputSelector.addEnabled = True
        self.outputSelector.removeEnabled = True
        self.outputSelector.noneEnabled = True
        self.outputSelector.showHidden = False
        self.outputSelector.showChildNodeTypes = False
        self.outputSelector.setMRMLScene(slicer.mrmlScene)
        self.outputSelector.setToolTip("Pick the output to the algorithm.")
        parametersFormLayout.addRow("Output Volume: ", self.outputSelector)

        #
        # Apply Button
        #
        self.applyButton = qt.QPushButton("Apply")
        self.applyButton.toolTip = "Run the algorithm."
        self.applyButton.enabled = False
        parametersFormLayout.addRow(self.applyButton)

        # connections
        self.applyButton.connect('clicked(bool)', self.onApplyButton)
        self.inputSelector.connect("currentNodeChanged(vtkMRMLNode*)", self.onSelect)
        self.outputSelector.connect("currentNodeChanged(vtkMRMLNode*)", self.onSelect)

        # Add vertical spacer
        self.layout.addStretch(1)

        # Refresh Apply button state
        self.onSelect()

    def cleanup(self):
        pass

    def onSelect(self):
        self.applyButton.enabled = self.inputSelector.currentNode() and self.outputSelector.currentNode()

    def onApplyButton(self):
        logic = TModelLogic()
        logic.run(self.inputSelector.currentNode(), self.armatureSelector.currentNode(),
                  self.outputSelector.currentNode())


#
# TModelLogic
#

class TModelLogic(ScriptedLoadableModuleLogic):
    """This class should implement all the actual
    computation done by your module.  The interface
    should be such that other python code can import
    this class and make use of the functionality without
    requiring an instance of the Widget.
    Uses ScriptedLoadableModuleLogic base class, available at:
    https://github.com/Slicer/Slicer/blob/master/Base/Python/slicer/ScriptedLoadableModule.py
    """

    def hasImageData(self, volumeNode):
        """This is an example logic method that
        returns true if the passed in volume
        node has valid image data
        """
        if not volumeNode:
            logging.debug('hasImageData failed: no volume node')
            return False
        if volumeNode.GetImageData() is None:
            logging.debug('hasImageData failed: no image data in volume node')
            return False
        return True

    def isValidInputOutputData(self, inputVolumeNode, armatureModel, outputVolumeNode):
        """Validates if the output is not the same as input
        """
        if not inputVolumeNode:
            logging.debug('isValidInputOutputData failed: no input volume node defined')
            return False
        if not outputVolumeNode:
            logging.debug('isValidInputOutputData failed: no output volume node defined')
            return False
        if inputVolumeNode.GetID() == outputVolumeNode.GetID():
            logging.debug(
                'isValidInputOutputData failed: input and output volume is the same. Create a new volume for output to avoid this error.')
            return False
        return True

    def run(self, inputVolume, armatureModel, outputVolume):
        """
        Run the actual algorithm
        """
        """
        TODO: Validate the input data!
        """

        if not self.isValidInputOutputData(inputVolume, armatureModel, outputVolume):
            slicer.util.errorDisplay('Input volume is the same as output volume. Choose a different output volume.')
            return False

        logging.info('Processing started')

        expandedVolume = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLScalarVolumeNode", "expandedVolume");
        cliParams = {'inputVolume': inputVolume.GetID(), 'outputVolume': expandedVolume}
        cliNode = slicer.cli.run(slicer.modules.expand, None, cliParams, wait_for_completion=True)

        """
        Create TORSO
        """
        # extract the Torso and put it in the output
        cliParams = {'InputVolume': expandedVolume, 'OutputVolume': outputVolume.GetID(),
                     'Lower': 4, 'Upper': 9, 'ThresholdType': 'Outside'}
        cliNode = slicer.cli.run(slicer.modules.thresholdmaster, None, cliParams, wait_for_completion=True)

        """
               **********************************************
               **********************************************
               *                                            *
               *             LEFT ARM BEGIN                 *
               *                                            *
               **********************************************
               **********************************************
               """
        """
        Start with the left arm and hand
        """
        leftHand = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLScalarVolumeNode", "leftHand")

        # Get the left hand
        cliParams = {'InputVolume': expandedVolume, 'OutputVolume': leftHand,
                     'Lower': 21, 'Upper': 21, 'ThresholdType': 'Outside'}
        cliNode = slicer.cli.run(slicer.modules.thresholdscalarvolume, None, cliParams, wait_for_completion=True)

        leftHandRotateVolume = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLScalarVolumeNode", "leftHandRotateVolume")

        # Compute the thresholded output volume using the Threshold Scalar Volume CLI module
        cliParams = {'inputVolume': leftHand, 'ComponentToRotate': 21,
                     'ComponentFixed': 20, 'ArmaturePoly': armatureModel.GetID(),
                     'outputVolume': leftHandRotateVolume}
        cliNode = slicer.cli.run(slicer.modules.logic, None, cliParams, wait_for_completion=True)

        # Get the left under arm
        leftUnderArm = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLScalarVolumeNode", "leftUnderArm")

        # Compute the thresholded output volume using the Threshold Scalar Volume CLI module
        cliParams = {'InputVolume': expandedVolume, 'OutputVolume': leftUnderArm,
                     'Lower': 20, 'Upper': 20, 'ThresholdType': 'Outside'}
        cliNode = slicer.cli.run(slicer.modules.thresholdscalarvolume, None, cliParams, wait_for_completion=True)

        """ Combine hand with underarm
        """
        leftUnderArmAndHand = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLScalarVolumeNode", "leftUnderArmAndHand")

        cliParams = {'inputVolume1': leftUnderArm, 'inputVolume2': leftHandRotateVolume,
                     'outputVolume': leftUnderArmAndHand}
        cliNode = slicer.cli.run(slicer.modules.addvolume, None, cliParams, wait_for_completion=True)

        ####
        leftUnderArmAndHandRotate = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLScalarVolumeNode",
                                                                       "leftUnderArmAndHandRotate")

        # Compute the thresholded output volume using the Threshold Scalar Volume CLI module
        cliParams = {'inputVolume': leftUnderArmAndHand, 'ComponentToRotate': 20,
                     'ComponentFixed': 19, 'ArmaturePoly': armatureModel.GetID(),
                     'outputVolume': leftUnderArmAndHandRotate}
        cliNode = slicer.cli.run(slicer.modules.logic, None, cliParams, wait_for_completion=True)

        # Get the left upper arm
        leftUpperArm = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLScalarVolumeNode", "leftUpperArm")

        # Compute the thresholded output volume using the Threshold Scalar Volume CLI module
        cliParams = {'InputVolume': expandedVolume, 'OutputVolume': leftUpperArm,
                     'Lower': 19, 'Upper': 19, 'ThresholdType': 'Outside'}
        cliNode = slicer.cli.run(slicer.modules.thresholdscalarvolume, None, cliParams, wait_for_completion=True)

        """ Combine hand with underarm
        """
        leftUpperArmAndRest = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLScalarVolumeNode", "leftUpperArmAndRest")

        cliParams = {'inputVolume1': leftUpperArm, 'inputVolume2': leftUnderArmAndHandRotate,
                     'outputVolume': leftUpperArmAndRest}
        cliNode = slicer.cli.run(slicer.modules.addvolume, None, cliParams, wait_for_completion=True)

        ####
        leftUpperArmAndRestRotate = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLScalarVolumeNode",
                                                                       "leftUpperArmAndRestRotate")

        # Compute the thresholded output volume using the Threshold Scalar Volume CLI module
        cliParams = {'inputVolume': leftUpperArmAndRest, 'ComponentToRotate': 19,
                     'ComponentFixed': 4, 'ArmaturePoly': armatureModel.GetID(),
                     'outputVolume': leftUpperArmAndRestRotate}
        cliNode = slicer.cli.run(slicer.modules.logic, None, cliParams, wait_for_completion=True)

        """
          add the completed parts to the torso
        """
        cliParams = {'inputVolume1': outputVolume.GetID(), 'inputVolume2': leftUpperArmAndRestRotate,
                     'outputVolume': outputVolume.GetID()}
        cliNode = slicer.cli.run(slicer.modules.addvolume, None, cliParams, wait_for_completion=True)


        """
        **********************************************
        **********************************************
        *                                            *
        *              RIGHT LEG BEGIN                *
        *                                            *
        **********************************************
        **********************************************
        """
        rightFootAndFemure = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLScalarVolumeNode", "rightFootAndFemure")

        # Get the right foot and tibea
        cliParams = {'InputVolume': expandedVolume, 'OutputVolume': rightFootAndFemure,
                     'Lower': 11, 'Upper': 12, 'ThresholdType': 'Outside'}
        cliNode = slicer.cli.run(slicer.modules.thresholdscalarvolume, None, cliParams, wait_for_completion=True)

        rightFootAndFemureRotate = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLScalarVolumeNode",
                                                                      "rightFootAndFemureRotate")

        # Compute the thresholded output volume using the Threshold Scalar Volume CLI module
        cliParams = {'inputVolume': rightFootAndFemure, 'ComponentToRotate': 11,
                     'ComponentFixed': 10, 'ArmaturePoly': armatureModel.GetID(),
                     'outputVolume': rightFootAndFemureRotate}
        cliNode = slicer.cli.run(slicer.modules.logic, None, cliParams, wait_for_completion=True)

        # Get the right under arm
        rightTibea = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLScalarVolumeNode", "rightTibea")

        # Compute the thresholded output volume using the Threshold Scalar Volume CLI module
        cliParams = {'InputVolume': expandedVolume, 'OutputVolume': rightTibea,
                     'Lower': 10, 'Upper': 10, 'ThresholdType': 'Outside'}
        cliNode = slicer.cli.run(slicer.modules.thresholdscalarvolume, None, cliParams, wait_for_completion=True)

        """ Combine hand with underarm
        """
        rightTibeaAndRest = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLScalarVolumeNode", "rightTibeaAndRest")

        cliParams = {'inputVolume1': rightTibea, 'inputVolume2': rightFootAndFemureRotate,
                     'outputVolume': rightTibeaAndRest}
        cliNode = slicer.cli.run(slicer.modules.addvolume, None, cliParams, wait_for_completion=True)

        ####
        rightTibeaAndRestRotate = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLScalarVolumeNode",
                                                                     "rightTibeaAndRestRotate")

        # Compute the thresholded output volume using the Threshold Scalar Volume CLI module
        cliParams = {'inputVolume': rightTibeaAndRest, 'ComponentToRotate': 10,
                     'ComponentFixed': 9, 'ArmaturePoly': armatureModel.GetID(),
                     'outputVolume': rightTibeaAndRestRotate}
        cliNode = slicer.cli.run(slicer.modules.logic, None, cliParams, wait_for_completion=True)

        """
          add the completed parts to the torso
        """
        cliParams = {'inputVolume1': outputVolume.GetID(), 'inputVolume2': rightTibeaAndRestRotate,
                     'outputVolume': outputVolume.GetID()}
        cliNode = slicer.cli.run(slicer.modules.addvolume, None, cliParams, wait_for_completion=True)

        """
               **********************************************
               **********************************************
               *                                            *
               *             LEFT LEG BEGIN                 *
               *                                            *
               **********************************************
               **********************************************
               """
        """
            Start with the left arm and hand
            """
        leftFootAndFemure = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLScalarVolumeNode", "leftFootAndFemure")

        # Get the left foot and tibea
        cliParams = {'InputVolume': expandedVolume, 'OutputVolume': leftFootAndFemure,
                     'Lower': 14, 'Upper': 15, 'ThresholdType': 'Outside'}
        cliNode = slicer.cli.run(slicer.modules.thresholdscalarvolume, None, cliParams, wait_for_completion=True)

        leftFootAndFemureRotate = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLScalarVolumeNode",
                                                                     "leftFootAndFemureRotate")

        # Compute the thresholded output volume using the Threshold Scalar Volume CLI module
        cliParams = {'inputVolume': leftFootAndFemure, 'ComponentToRotate': 14,
                     'ComponentFixed': 13, 'ArmaturePoly': armatureModel.GetID(),
                     'outputVolume': leftFootAndFemureRotate}
        cliNode = slicer.cli.run(slicer.modules.logic, None, cliParams, wait_for_completion=True)

        # Get the left under arm
        leftTibea = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLScalarVolumeNode", "leftTibea")

        # Compute the thresholded output volume using the Threshold Scalar Volume CLI module
        cliParams = {'InputVolume': expandedVolume, 'OutputVolume': leftTibea,
                     'Lower': 13, 'Upper': 13, 'ThresholdType': 'Outside'}
        cliNode = slicer.cli.run(slicer.modules.thresholdscalarvolume, None, cliParams, wait_for_completion=True)

        """ Combine hand with underarm
        """
        leftTibeaAndRest = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLScalarVolumeNode", "leftTibeaAndRest")

        cliParams = {'inputVolume1': leftTibea, 'inputVolume2': leftFootAndFemureRotate,
                     'outputVolume': leftTibeaAndRest}
        cliNode = slicer.cli.run(slicer.modules.addvolume, None, cliParams, wait_for_completion=True)

        ####
        leftTibeaAndRestRotate = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLScalarVolumeNode",
                                                                    "leftTibeaAndRestRotate")

        # Compute the thresholded output volume using the Threshold Scalar Volume CLI module
        cliParams = {'inputVolume': leftTibeaAndRest, 'ComponentToRotate': 13,
                     'ComponentFixed': 8, 'ArmaturePoly': armatureModel.GetID(),
                     'outputVolume': leftTibeaAndRestRotate}
        cliNode = slicer.cli.run(slicer.modules.logic, None, cliParams, wait_for_completion=True)

        """
          add the completed parts to the torso
        """
        cliParams = {'inputVolume1': outputVolume.GetID(), 'inputVolume2': leftTibeaAndRestRotate,
                     'outputVolume': outputVolume.GetID()}
        cliNode = slicer.cli.run(slicer.modules.addvolume, None, cliParams, wait_for_completion=True)

        """
        **********************************************
        **********************************************
        *                                            *
        *             HEAD AND NACK BEGIN            *
        *                                            *
        **********************************************
        **********************************************
        """
        headVolume = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLScalarVolumeNode", "head")

        # Get the head
        cliParams = {'InputVolume': expandedVolume, 'OutputVolume': headVolume,
                     'Lower': 2, 'Upper': 2, 'ThresholdType': 'Outside'}
        cliNode = slicer.cli.run(slicer.modules.thresholdmaster, None, cliParams, wait_for_completion=True)

        headRotateVolume = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLScalarVolumeNode", "headRotateVolume")

        # Compute the thresholded output volume using the Threshold Scalar Volume CLI module
        cliParams = {'inputVolume': headVolume, 'ComponentToRotate': 2,
                     'ComponentFixed': 3, 'ArmaturePoly': armatureModel.GetID(),
                     'outputVolume': headRotateVolume}
        cliNode = slicer.cli.run(slicer.modules.logic, None, cliParams, wait_for_completion=True)

        """Get the nack
        """
        nackVolume = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLScalarVolumeNode", "nackVolume")

        # Compute the thresholded output volume using the Threshold Scalar Volume CLI module
        cliParams = {'InputVolume': expandedVolume, 'OutputVolume': nackVolume,
                     'Lower': 3, 'Upper': 3, 'ThresholdType': 'Outside'}
        cliNode = slicer.cli.run(slicer.modules.thresholdmaster, None, cliParams, wait_for_completion=True)

        """ Combine Nack with new Head
        """
        nackHead = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLScalarVolumeNode", "nackHead")

        cliParams = {'inputVolume1': nackVolume, 'inputVolume2': headRotateVolume, 'outputVolume': nackHead}
        cliNode = slicer.cli.run(slicer.modules.addvolume, None, cliParams, wait_for_completion=True)

        """ Rotate Nack Head against Torso
        """
        nackHeadRotate = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLScalarVolumeNode", "nackHeadRotate")

        # Compute the thresholded output volume using the Threshold Scalar Volume CLI module
        cliParams = {'inputVolume': nackHead, 'ComponentToRotate': 3,
                     'ComponentFixed': 6, 'ArmaturePoly': armatureModel.GetID(),
                     'outputVolume': nackHeadRotate}
        cliNode = slicer.cli.run(slicer.modules.logic, None, cliParams, wait_for_completion=True)

        """
          add the completed parts to the torso
        """
        cliParams = {'inputVolume1': outputVolume.GetID(), 'inputVolume2': nackHeadRotate,
                     'outputVolume': outputVolume.GetID()}
        cliNode = slicer.cli.run(slicer.modules.addvolume, None, cliParams, wait_for_completion=True)

        """
        **********************************************
        **********************************************
        *                                            *
        *             RIGHT ARM BEGIN                *
        *                                            *
        **********************************************
        **********************************************
        """

        """
          Start with the right arm and hand
        """
        rightHand = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLScalarVolumeNode", "rightHand")

        # Get the right hand
        cliParams = {'InputVolume': expandedVolume, 'OutputVolume': rightHand,
                     'Lower': 18, 'Upper': 18, 'ThresholdType': 'Outside'}
        cliNode = slicer.cli.run(slicer.modules.thresholdscalarvolume, None, cliParams, wait_for_completion=True)

        rightHandRotateVolume = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLScalarVolumeNode", "rightHandRotateVolume")

        # Compute the thresholded output volume using the Threshold Scalar Volume CLI module
        cliParams = {'inputVolume': rightHand, 'ComponentToRotate': 18,
                     'ComponentFixed': 17, 'ArmaturePoly': armatureModel.GetID(),
                     'outputVolume': rightHandRotateVolume}
        cliNode = slicer.cli.run(slicer.modules.logic, None, cliParams, wait_for_completion=True)

        # Get the right under arm
        rightUnderArm = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLScalarVolumeNode", "rightUnderArm")

        # Compute the thresholded output volume using the Threshold Scalar Volume CLI module
        cliParams = {'InputVolume': expandedVolume, 'OutputVolume': rightUnderArm,
                     'Lower': 17, 'Upper': 17, 'ThresholdType': 'Outside'}
        cliNode = slicer.cli.run(slicer.modules.thresholdscalarvolume, None, cliParams, wait_for_completion=True)

        """ Combine hand with underarm
        """
        rightUnderArmAndHand = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLScalarVolumeNode", "rightUnderArmAndHand")

        cliParams = {'inputVolume1': rightUnderArm, 'inputVolume2': rightHandRotateVolume,
                     'outputVolume': rightUnderArmAndHand}
        cliNode = slicer.cli.run(slicer.modules.addvolume, None, cliParams, wait_for_completion=True)

        ####
        rightUnderArmAndHandRotate = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLScalarVolumeNode",
                                                                        "rightUnderArmAndHandRotate")

        # Compute the thresholded output volume using the Threshold Scalar Volume CLI module
        cliParams = {'inputVolume': rightUnderArmAndHand, 'ComponentToRotate': 17,
                     'ComponentFixed': 16, 'ArmaturePoly': armatureModel.GetID(),
                     'outputVolume': rightUnderArmAndHandRotate}
        cliNode = slicer.cli.run(slicer.modules.logic, None, cliParams, wait_for_completion=True)

        # Get the right upper arm
        rightUpperArm = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLScalarVolumeNode", "rightUpperArm")

        # Compute the thresholded output volume using the Threshold Scalar Volume CLI module
        cliParams = {'InputVolume': expandedVolume, 'OutputVolume': rightUpperArm,
                     'Lower': 16, 'Upper': 16, 'ThresholdType': 'Outside'}
        cliNode = slicer.cli.run(slicer.modules.thresholdscalarvolume, None, cliParams, wait_for_completion=True)

        """ Combine hand with underarm
        """
        rightUpperArmAndRest = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLScalarVolumeNode", "rightUpperArmAndRest")

        cliParams = {'inputVolume1': rightUpperArm, 'inputVolume2': rightUnderArmAndHandRotate,
                     'outputVolume': rightUpperArmAndRest}
        cliNode = slicer.cli.run(slicer.modules.addvolume, None, cliParams, wait_for_completion=True)

        ####
        rightUpperArmAndRestRotate = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLScalarVolumeNode",
                                                                        "rightUpperArmAndRestRotate")

        # Compute the thresholded output volume using the Threshold Scalar Volume CLI module
        cliParams = {'inputVolume': rightUpperArmAndRest, 'ComponentToRotate': 16,
                     'ComponentFixed': 5, 'ArmaturePoly': armatureModel.GetID(),
                     'outputVolume': rightUpperArmAndRestRotate}
        cliNode = slicer.cli.run(slicer.modules.logic, None, cliParams, wait_for_completion=True)

        """
          add the completed parts to the torso
        """
        cliParams = {'inputVolume1': outputVolume.GetID(), 'inputVolume2': rightUpperArmAndRestRotate,
                     'outputVolume': outputVolume.GetID()}
        cliNode = slicer.cli.run(slicer.modules.addvolume, None, cliParams, wait_for_completion=True)

        logging.info('Processing completed')

        return True


class TModelTest(ScriptedLoadableModuleTest):
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
        self.test_TModel1()

    def test_TModel1(self):
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

        for url, name, loader in downloads:
            filePath = slicer.app.temporaryPath + '/' + name
            if not os.path.exists(filePath) or os.stat(filePath).st_size == 0:
                logging.info('Requesting download %s from %s...\n' % (name, url))
                urllib.urlretrieve(url, filePath)
            if loader:
                logging.info('Loading %s...' % (name,))
                loader(filePath)
        self.delayDisplay('Finished with download and loading')

        volumeNode = slicer.util.getNode(pattern="FA")
        logic = TModelLogic()
        self.assertIsNotNone(logic.hasImageData(volumeNode))
        self.delayDisplay('Test passed!')
