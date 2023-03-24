"""


:author: 
:contact: 
:email: 
:organization: 
:address: 
:copyright: 
:date: Mar 20 2023 12:01
:dragonflyVersion: 2022.2.0.1367
:UUID: 8833cb2cc73811edab7744032c94bd8e
"""

__version__ = '1.0.0'

import pandas as pd
import openpyxl

from PyQt5.QtWidgets import QFileDialog

import ORSModel
from ORSModel import MultiROI
from ORSModel import Managed
from ORSModel import ROI
from ORSModel import VisualRuler
from ORSModel import Vector3

from ORSServiceClass.actionAndMenu.menu import Menu
from ORSServiceClass.decorators.infrastructure import interfaceMethod
from ORSServiceClass.menuItems.contextualMenuItem import ContextualMenuItem
from ORSServiceClass.ORSWidget.chooseObjectAndNewName.chooseObjectAndNewName import ChooseObjectAndNewName
from ORSServiceClass.ORSWidget.SimpleEntryDialog.simpleEntryDialog import SimpleEntryDialog

from OrsHelpers.ListHelper import ListHelper
from OrsHelpers.multiroilabelhelper import MultiROILabelHelper
from OrsHelpers.primitivehelper import PrimitiveHelper

from OrsLibraries.workingcontext import WorkingContext


class DrawPairingsWithCOM_8833cb2cc73811edab7744032c94bd8e(ContextualMenuItem):

    @classmethod
    def getIsSelectionValid(cls, aCollectionOfObjects, implementation):
        """
        :param aCollectionOfObjects: a list of objects currently being selected, i.e. on which the menu item could be applied.
        :param implementation: a subclass of AbstractPlugin in the current context.
        :return: if True, the menu item will be displayed.
        """

        if aCollectionOfObjects is None:
            return False

        if len(aCollectionOfObjects) < 1:
            return False

        selectionIsOnlyMultiROIs = all([isinstance(obj, MultiROI) for obj in aCollectionOfObjects])
        return selectionIsOnlyMultiROIs

    @classmethod
    def getMenuItemForSelection(cls, aCollectionOfObjects, implementation):
        """
        Returns the menu item
        :param aCollectionOfObjects: a list of objects currently being selected, i.e. on which the menu item will be applied.
        :param implementation: a subclass of AbstractPlugin in the current context.
        :return: Menu
        """
        
        collectionString = ListHelper.asPythonCollectionString(aCollectionOfObjects)
        myMenu = Menu(title='Draw COM Pairing Rulers',
                      id_='DrawPairingsWithCOM_8833cb2cc73811edab7744032c94bd8e',
                      section='',
                      action=f'DrawPairingsWithCOM_8833cb2cc73811edab7744032c94bd8e.menuItemEntryPoint({collectionString}, {implementation.getVarName()})')
        return myMenu

    @classmethod
    def menuItemEntryPoint(cls, collectionString, implementation):
        """
        Will be executed when the menu item is selected.
        :param collectionString: a list of objects representation currently being selected, i.e. on which the menu item will be applied.
        :param implementation: a subclass of AbstractPlugin in the current context.
        """
        
        # aCollectionOfObjects is a Python list of objects currently being selected
        aCollectionOfObjects = ListHelper.fromPythonCollection(collectionString, asPythonList=True)

        synMROI = aCollectionOfObjects[0]
        vesMROI = cls._selectMultiROI()

        selectedFilename, filter = QFileDialog.getOpenFileName(caption="Please select a CSV with your COM information")

        if selectedFilename == '':
            return

        df = pd.read_excel(selectedFilename, engine='openpyxl')
        df1 = df[['synLabel', 'vesLabel']]
        print(df1.head())

        synIndex = cls._getIndexValue()
        vesIndex = int(df1.iloc[synIndex-1]['vesLabel'])

        synROI = cls._exportLabelToROI(synMROI, synIndex)
        vesROI = cls._exportLabelToROI(vesMROI, vesIndex)

        aLayoutName = 'toplayout\\scene_0\\0\\3D'
        associatedState = 'OrsStateRulerEdit'
        t1 = 0
        t2 = 0
        p1 = synROI.getCenterOfMass(pTimeStep=0)
        p2 = vesROI.getCenterOfMass(pTimeStep=0)
        cls._createPrimitive(aLayoutName, associatedState, t1, p1, t2, p2)

    @classmethod
    @interfaceMethod
    def _selectMultiROI(cls):
        class_name = "MultiROI"
        class_object = getattr(ORSModel, class_name)
        return ChooseObjectAndNewName.prompt([class_object], parent=WorkingContext.getCurrentContextWindow(), dialog_title="Select the MultiROI that pairs with your selected MultiROI", allowNone=True, getNewTitle=False)

    @classmethod
    @interfaceMethod
    def _getIndexValue(cls):
        value = SimpleEntryDialog.prompt(None, "Integer", "Select Index", 1)
        try:
            return int(value)
        except ValueError:
            return 1

    @classmethod
    @interfaceMethod
    def _exportLabelToROI(cls, multiROI, index):
        """
        Exports a selected label from a MultiROI as an ROI

        :param multiROI: multiROI to use
        :type multiROI: ORSModel.ors.MultiROI
        :param index: label index to export
        :type index: int
        """

        _idx = ORSModel.ors.ArrayUnsignedLong()
        _idx.insertAt(index=0, pValue=index)

        multiROI.setSelectedLabels(iTIndex=0, labels=_idx, selected=True)
        guid = MultiROILabelHelper.extractSelectedLabelsToROIs(multiroi=multiROI, tIndex=0)[0]
        guidROI = Managed.getObjectWithGUID(guid)
        guidROI.publish(logging=True)
        return(guidROI)

    @classmethod
    @interfaceMethod
    def _createPrimitive(cls, aLayoutName, associatedState, t1, p1, t2, p2):
        """
        Creates and publishes a Ruler

        :param aLayoutName: Name of Layout
        :type aLayoutName: str
        :param associatedState: State to be used
        :type associatedState: str
        :param t1: TimeStep for point 1
        :type t1: int
        :param p1: Position for point 1
        :type p1: ORSModel.ors.Vector3
        :param t2: TimeStep for point 2
        :type t2: int
        :param p2: Position for point 1
        :type p2: ORSModel.ors.Vector3
        """

        # Creating a blank ruler primitive
        newAnnotation = PrimitiveHelper.createPrimitive(primitiveClass=ORSModel.ors.VisualRuler, aLayoutName=aLayoutName,
                                                        associatedState=associatedState)

        # Create and move the ruler points
        success = PrimitiveHelper.addControlPoint(anAnnotation=newAnnotation,
                                                  timeStep=t1,
                                                  position=p1)

        success1 = PrimitiveHelper.addControlPoint(anAnnotation=newAnnotation,
                                                  timeStep=t2,
                                                  position=p2)

        # Publishing the ruler
        newAnnotation.setIsRepresentable(isRepresentable=True, logging=True)
        newAnnotation.publish(logging=True)
        newAnnotation.setDataDirty(logging=True)
