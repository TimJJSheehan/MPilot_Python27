from __future__ import division
from MPCore import MPilotEEMSFxnParent as mpefp
import numpy as np
import copy as cp

class Copy(mpefp._MPilotEEMSFxnParent):

    def __init__(self,mptCmdStruct=None):

        super(Copy, self).__init__(
            mptCmdStruct=mptCmdStruct,
            dataType=None,      # Set at runtime
            isDataLayer = None  # Set at runtime
            )

    # def __init__(self,mptCmdStruct=None):

    def _SetFxnDesc(self):
        
        self.fxnDesc['DisplayName'] = 'Copy'
        self.fxnDesc['ShortDesc'] = 'Copies the data from another field'
        self.fxnDesc['ReturnType'] = 'Any'
        
        self.fxnDesc['ReqArgs'] = {
            'InFieldName':'Field Name'
            }
        self.fxnDesc['OptArgs'] = {
            'OutFileName':'File Name',
            'Metadata':'Any',
            'PrecursorFieldNames':['Field Name','Field Name List']
            }
        
    # _SetFxnDesc(self):

    def DependencyNms(self):
        
        rtrn = self._ArgToList('InFieldName')
        rtrn += self._ArgToList('PrecursorFieldNames')
        return rtrn
    
    # def DependencyNms(self):
        
    def Exec(self,executedObjects):

        fldNmObj = executedObjects[self.ValFromArgByNm('InFieldName')]
        self.execRslt = cp.deepcopy(fldNmObj.ExecRslt())
        self.dataType = fldNmObj.DataType()
        self.isDataLayer = fldNmObj.IsDataLayer()

        executedObjects[self.RsltNm()] = self
        
    # def Exec(self,executedObjects):
    
# class Copy(mpefp._MPilotEEMSFxnParent):
        
class AMinusB(mpefp._MPilotEEMSFxnParent):
    
    def __init__(self,mptCmdStruct=None):

        super(AMinusB, self).__init__(
            mptCmdStruct=mptCmdStruct,
            dataType=None, # set at run time
            isDataLayer = True
            )

    # def __init__(self,mptCmdStruct=None):

    def _SetFxnDesc(self):
        
        self.fxnDesc['DisplayName'] = 'A Minus B'
        self.fxnDesc['ShortDesc'] = 'Performs A - B'
        self.fxnDesc['ReturnType'] = [
            'Integer',
            'Positive Integer',
            'Float',
            'Positive Float'
            ]                        
        
        self.fxnDesc['ReqArgs'] = {
            'A':'Field Name',
            'B':'Field Name'
            }
        self.fxnDesc['OptArgs'] = {
            'OutFileName':'File Name',
            'Metadata':'Any',
            'PrecursorFieldNames':['Field Name','Field Name List']
            }
        
    # _SetFxnDesc(self):

    def DependencyNms(self):
        
        rtrn = self._ArgToList('A')
        rtrn += self._ArgToList('B')
        rtrn += self._ArgToList('PrecursorFieldNames')
        return rtrn
    
    # def DependencyNms(self):
    
    def Exec(self,executedObjects):

        aObj = executedObjects[self.ValFromArgByNm('A')]
        bObj = executedObjects[self.ValFromArgByNm('B')]
        
        self._ValidateIsDataLayer(aObj)
        self._ValidateIsDataLayer(bObj)
        
        self._ValidateProgDataType(
            'A',
            aObj,
            ['Integer','Positive Integer','Float','Positive Float']
            )
        self._ValidateProgDataType(
            'B',
            bObj,
            ['Integer','Positive Integer','Float','Positive Float']
            )

        self.execRslt = aObj.ExecRslt() - bObj.ExecRslt()

        if aObj.DataType() in ['Integer','Positive Integer'] and \
            bObj.DataType() in ['Integer','Positive Integer']:
            if self.ExecRslt().min() > 0:
                self.dataType = 'Positive Integer'
            else:
                self.dataType = 'Integer'
        else:
            if self.ExecRslt().min() > 0:
                self.dataType = 'Positive Float'
            else:
                self.dataType = 'Float'
            
        executedObjects[self.RsltNm()] = self
        
    # def Exec(self,executedObjects):
    
# class AMinusB(mpefp._MPilotEEMSFxnParent):
        
class Sum(mpefp._MPilotEEMSFxnParent):

    def __init__(self,mptCmdStruct=None):

        super(Sum, self).__init__(
            mptCmdStruct=mptCmdStruct,
            dataType=None,      # Set at runtime
            isDataLayer = True
            )

    # def __init__(self,mptCmdStruct=None):
    
    def _SetFxnDesc(self):
        # description of command used for validation and
        # information display Each Pilot fxn command should have its
        # own description

        self.fxnDesc['DisplayName'] = 'Sum'
        self.fxnDesc['ShortDesc'] = 'Sums input variables'
        self.fxnDesc['ReturnType'] = [
            'Integer',
            'Positive Integer',
            'Float',
            'Positive Float'
            ]                        

        self.fxnDesc['ReqArgs'] = {
            'InFieldNames':['Field Name List']
            }
        self.fxnDesc['OptArgs'] = {
            'OutFileName':'File Name',
            'Metadata':'Any',
            'PrecursorFieldNames':['Field Name','Field Name List']
            }
        
    # _SetFxnDesc(self):

    def Exec(self,executedObjects):
        # executedObjects is a dictionary of executed MPilot function objects
        
        self._ValidateListLen('InFldNms',1)
        fldNms = self._ArgToList('InFieldNames')

        # Validate input data and gather types
        inDataTypes = {}
        for fldNm in fldNms:
            self._ValidateProgDataType(
                fldNm,
                executedObjects[fldNm],
                ['Integer','Positive Integer','Float','Positive Float']
                )
            self._ValidateIsDataLayer(executedObjects[fldNm])
            inDataTypes[executedObjects[fldNm].DataType()] = True
        
        # Result starts with a copy of the first input field.
        
        self.execRslt = cp.deepcopy(executedObjects[fldNms[0]].ExecRslt())
        
        for fldNm in fldNms[1:]:
            self.execRslt += executedObjects[fldNm].ExecRslt()

        # Determine type of output
        if 'Float' in inDataTypes.keys() or 'Positive Float' in inDataTypes.keys():
            if self.ExecRslt().min() > 0:
                self.dataType = 'Positive Float'
            else:
                self.dataType = 'Float'
        else:
            if self.ExecRslt().min() > 0:
                self.dataType = 'Positive Integer'
            else:
                self.dataType = 'Integer'        

        executedObjects[self.RsltNm()] = self
        
    # def Exec(self,executedObjects):
    
# class Sum(mpefp._MPilotEEMSFxnParent):
        
class WeightedSum(mpefp._MPilotEEMSFxnParent):

    def __init__(self,mptCmdStruct=None):

        super(WeightedSum, self).__init__(
            mptCmdStruct=mptCmdStruct,
            dataType=None,      # Set at runtime
            isDataLayer = True
            )

    # def __init__(self,mptCmdStruct=None):
    
    def _SetFxnDesc(self):
        # description of command used for validation and
        # information display Each Pilot fxn command should have its
        # own description

        self.fxnDesc['DisplayName'] = 'Weighted Sum'
        self.fxnDesc['ShortDesc'] = 'Takes the weighted sum of input variables'
        self.fxnDesc['ReturnType'] = [
            'Float',
            'Positive Float'
            ]                        

        self.fxnDesc['ReqArgs'] = {
            'InFieldNames':['Field Name List'],
            'Weights':['Float','Float List']            
            }
        self.fxnDesc['OptArgs'] = {
            'OutFileName':'File Name',
            'Metadata':'Any',
            'PrecursorFieldNames':['Field Name','Field Name List']
            }
        
    # _SetFxnDesc(self):

    def Exec(self,executedObjects):
        # executedObjects is a dictionary of executed MPilot function objects
        
        self._ValidateListLen('InFldNms',1)
        self._ValidateEqualListLens(['InFieldNames','Weights'])
        
        fldNms = self._ArgToList('InFieldNames')
        wts = self.ValFromArgByNm('Weights')

        # Validate input data and gather types
        for fldNm in fldNms:
            self._ValidateProgDataType(
                fldNm,
                executedObjects[fldNm],
                ['Integer','Positive Integer','Float','Positive Float']
                )
            self._ValidateIsDataLayer(executedObjects[fldNm])
        
        # Result starts with a copy of the first input field.
        self.execRslt = cp.deepcopy(executedObjects[fldNms[0]].ExecRslt()) * wts[0]
        
        for wt,fldNm in zip(wts,fldNms)[1:]:
            self.execRslt += executedObjects[fldNm].ExecRslt() * wt

        # Determine type of output
        if self.ExecRslt().min() > 0:
            self.dataType = 'Positive Float'
        else:
            self.dataType = 'Float'

        executedObjects[self.RsltNm()] = self
        
    # def Exec(self,executedObjects):

# class WeightedSum(mpefp._MPilotEEMSFxnParent):
        
class Multiply(mpefp._MPilotEEMSFxnParent):

    def __init__(self,mptCmdStruct=None):

        super(Multiply, self).__init__(
            mptCmdStruct=mptCmdStruct,
            dataType=None,      # Set at runtime
            isDataLayer = True
            )

    # def __init__(self,mptCmdStruct=None):
    
    def _SetFxnDesc(self):
        # description of command used for validation and
        # information display Each Pilot fxn command should have its
        # own description

        self.fxnDesc['DisplayName'] = 'Multiply'
        self.fxnDesc['ShortDesc'] = 'Multiplies input variables'
        self.fxnDesc['ReturnType'] = [
            'Integer',
            'Positive Integer',
            'Float',
            'Positive Float'
            ]                        

        self.fxnDesc['ReqArgs'] = {
            'InFieldNames':['Field Name List']
            }
        self.fxnDesc['OptArgs'] = {
            'OutFileName':'File Name',
            'Metadata':'Any',
            'PrecursorFieldNames':['Field Name','Field Name List']
            }
        
    # _SetFxnDesc(self):

    def Exec(self,executedObjects):
        # executedObjects is a dictionary of executed MPilot function objects
        
        self._ValidateListLen('InFldNms',1)
        fldNms = self._ArgToList('InFieldNames')

        # Validate input data and gather types
        inDataTypes = {}
        for fldNm in fldNms:
            self._ValidateProgDataType(
                fldNm,
                executedObjects[fldNm],
                ['Integer','Positive Integer','Float','Positive Float']
                )
            self._ValidateIsDataLayer(executedObjects[fldNm])
            inDataTypes[executedObjects[fldNm].DataType()] = True
        
        # Result starts with a copy of the first input field.
        
        self.execRslt = cp.deepcopy(executedObjects[fldNms[0]].ExecRslt())
        
        for fldNm in fldNms[1:]:
            self.execRslt *= executedObjects[fldNm].ExecRslt()

        # Determine type of output
        if 'Float' in inDataTypes.keys() or 'Positive Float' in inDataTypes.keys():
            if self.ExecRslt().min() > 0:
                self.dataType = 'Positive Float'
            else:
                self.dataType = 'Float'
        else:
            if self.ExecRslt().min() > 0:
                self.dataType = 'Positive Integer'
            else:
                self.dataType = 'Integer'        

        executedObjects[self.RsltNm()] = self
        
    # def Exec(self,executedObjects):
    
# class Multiply(mpefp._MPilotEEMSFxnParent):
        
class ADividedByB(mpefp._MPilotEEMSFxnParent):
    
    def __init__(self,mptCmdStruct=None):

        super(ADividedByB, self).__init__(
            mptCmdStruct=mptCmdStruct,
            dataType=None, # set at run time
            isDataLayer = True
            )

    # def __init__(self,mptCmdStruct=None):

    def _SetFxnDesc(self):
        
        self.fxnDesc['DisplayName'] = 'A Divided By B'
        self.fxnDesc['ShortDesc'] = 'Performs A / B'
        self.fxnDesc['ReturnType'] = [
            'Float',
            'Positive Float'
            ]                        
        
        self.fxnDesc['ReqArgs'] = {
            'A':'Field Name',
            'B':'Field Name'
            }
        self.fxnDesc['OptArgs'] = {
            'OutFileName':'File Name',
            'Metadata':'Any',
            'PrecursorFieldNames':['Field Name','Field Name List']
            }
        
    # _SetFxnDesc(self):

    def DependencyNms(self):
        
        rtrn = self._ArgToList('A')
        rtrn += self._ArgToList('B')
        rtrn += self._ArgToList('PrecursorFieldNames')
        return rtrn
    
    # def DependencyNms(self):
    
    def Exec(self,executedObjects):

        aObj = executedObjects[self.ValFromArgByNm('A')]
        bObj = executedObjects[self.ValFromArgByNm('B')]
        
        self._ValidateIsDataLayer(aObj)
        self._ValidateIsDataLayer(bObj)
        
        self._ValidateProgDataType(
            'A',
            aObj,
            ['Integer','Positive Integer','Float','Positive Float']
            )
        self._ValidateProgDataType(
            'B',
            bObj,
            ['Integer','Positive Integer','Float','Positive Float']
            )

        self.execRslt = aObj.ExecRslt() / bObj.ExecRslt()

        if self.ExecRslt().min() > 0:
            self.dataType = 'Positive Float'
        else:
            self.dataType = 'Float'
            
        executedObjects[self.RsltNm()] = self
        
    # def Exec(self,executedObjects):
    
# class Divide(mpefp._MPilotEEMSFxnParent):
        
class Minimum(mpefp._MPilotEEMSFxnParent):

    def __init__(self,mptCmdStruct=None):

        super(Minimum, self).__init__(
            mptCmdStruct=mptCmdStruct,
            dataType=None,   # Set at runtime
            isDataLayer=True
            )

    # def __init__(self,mptCmdStruct=None):

    def _SetFxnDesc(self):
        
        self.fxnDesc['DisplayName'] = 'Minimum'
        self.fxnDesc['ShortDesc'] = 'Takes the minimum input variables'
        self.fxnDesc['ReturnType'] = [
            'Integer',
            'Positive Integer',
            'Float',
            'Positive Float'
            ]
        
        self.fxnDesc['ReqArgs'] = {
            'InFieldNames':['Field Name','Field Name List']
            }
        self.fxnDesc['OptArgs'] = {
            'OutFileName':'File Name',
            'Metadata':'Any',
            'PrecursorFieldNames':['Field Name','Field Name List']
            }
        
    # _SetFxnDesc(self):

    def Exec(self,executedObjects):

        self._ValidateListLen('InFldNms',1)
        fldNms = self._ArgToList('InFieldNames')
        
        # Validate input data and gather types
        inDataTypes = {}
        for fldNm in fldNms:
            self._ValidateProgDataType(
                fldNm,
                executedObjects[fldNm],
                ['Integer','Positive Integer','Float','Positive Float']
                )
            self._ValidateIsDataLayer(executedObjects[fldNm])
            inDataTypes[executedObjects[fldNm].DataType()] = True
        
        # Traverse inputs and take minimum
        self.execRslt = cp.deepcopy(executedObjects[fldNms[0]].ExecRslt())

        for fldNm in fldNms[1:]:
            self.execRslt = np.ma.minimum(self.execRslt,executedObjects[fldNm].ExecRslt())

        # Determine type of output
        if 'Float' in inDataTypes.keys() or 'Positive Float' in inDataTypes.keys():
            if self.ExecRslt().min() > 0:
                self.dataType = 'Positive Float'
            else:
                self.dataType = 'Float'
        else:
            if self.ExecRslt().min() > 0:
                self.dataType = 'Positive Integer'
            else:
                self.dataType = 'Integer'        

        executedObjects[self.RsltNm()] = self
        
    # def Exec(self,executedObjects):
    
# class Minimum(mpefp._MPilotEEMSFxnParent):
        
class Maximum(mpefp._MPilotEEMSFxnParent):

    def __init__(self,mptCmdStruct=None):

        super(Maximum, self).__init__(
            mptCmdStruct=mptCmdStruct,
            dataType=None,   # Set at runtime
            isDataLayer=True
            )

    # def __init__(self,mptCmdStruct=None):

    def _SetFxnDesc(self):
        
        self.fxnDesc['DisplayName'] = 'Maximum'
        self.fxnDesc['ShortDesc'] = 'Takes the maximum input variables'
        self.fxnDesc['ReturnType'] = [
            'Integer',
            'Positive Integer',
            'Float',
            'Positive Float'
            ]
        
        self.fxnDesc['ReqArgs'] = {
            'InFieldNames':['Field Name','Field Name List']
            }
        self.fxnDesc['OptArgs'] = {
            'OutFileName':'File Name',
            'Metadata':'Any',
            'PrecursorFieldNames':['Field Name','Field Name List']
            }
        
    # _SetFxnDesc(self):

    def Exec(self,executedObjects):

        self._ValidateListLen('InFldNms',1)
        fldNms = self._ArgToList('InFieldNames')
        
        # Validate input data and gather types
        inDataTypes = {}
        for fldNm in fldNms:
            self._ValidateProgDataType(
                fldNm,
                executedObjects[fldNm],
                ['Integer','Positive Integer','Float','Positive Float']
                )
            self._ValidateIsDataLayer(executedObjects[fldNm])
            inDataTypes[executedObjects[fldNm].DataType()] = True
        
        # Traverse inputs and take minimum
        self.execRslt = cp.deepcopy(executedObjects[fldNms[0]].ExecRslt())

        for fldNm in fldNms[1:]:
            self.execRslt = np.ma.maximum(self.execRslt,executedObjects[fldNm].ExecRslt())

        # Determine type of output
        if 'Float' in inDataTypes.keys() or 'Positive Float' in inDataTypes.keys():
            if self.ExecRslt().min() > 0:
                self.dataType = 'Positive Float'
            else:
                self.dataType = 'Float'
        else:
            if self.ExecRslt().min() > 0:
                self.dataType = 'Positive Integer'
            else:
                self.dataType = 'Integer'        

        executedObjects[self.RsltNm()] = self
        
    # def Exec(self,executedObjects):
    
# class Maximum(mpefp._MPilotEEMSFxnParent):
        
class Mean(mpefp._MPilotEEMSFxnParent):

    def __init__(self,mptCmdStruct=None):

        super(Mean, self).__init__(
            mptCmdStruct=mptCmdStruct,
            dataType=None,
            isDataLayer = True
            )

    # def __init__(self,mptCmdStruct=None):
    
    def _SetFxnDesc(self):
        # description of command used for validation and
        # information display Each Pilot fxn command should have its
        # own description

        self.fxnDesc['DisplayName'] = 'Mean'
        self.fxnDesc['ShortDesc'] = 'Mean of input variables'
        self.fxnDesc['ReturnType'] = [
            'Float',
            'Positive Float'
            ]                        

        self.fxnDesc['ReqArgs'] = {
            'InFieldNames':['Field Name List']
            }
        self.fxnDesc['OptArgs'] = {
            'OutFileName':'File Name',
            'Metadata':'Any',
            'PrecursorFieldNames':['Field Name','Field Name List']
            }
        
    # _SetFxnDesc(self):

    def Exec(self,executedObjects):
        # executedObjects is a dictionary of executed MPilot function objects
        
        self._ValidateListLen('InFldNms',1)
        fldNms = self._ArgToList('InFieldNames')

        # Validate input data and gather types
        inDataTypes = {}
        for fldNm in fldNms:
            self._ValidateProgDataType(
                fldNm,
                executedObjects[fldNm],
                ['Integer','Positive Integer','Float','Positive Float']
                )
            self._ValidateIsDataLayer(executedObjects[fldNm])
            inDataTypes[executedObjects[fldNm].DataType()] = True
        
        # Result starts with a copy of the first input field.
        
        self.execRslt = cp.deepcopy(executedObjects[fldNms[0]].ExecRslt())
        
        for fldNm in fldNms[1:]:
            self.execRslt += executedObjects[fldNm].ExecRslt()

        self.execRslt = self.execRslt / len(fldNms)

        # Determine type of output
        if self.ExecRslt().min() > 0:
            self.dataType = 'Positive Float'
        else:
            self.dataType = 'Float'

        executedObjects[self.RsltNm()] = self
        
    # def Exec(self,executedObjects):
    
# class Mean(mpefp._MPilotEEMSFxnParent):
        
class WeightedMean(mpefp._MPilotEEMSFxnParent):

    def __init__(self,mptCmdStruct=None):

        super(WeightedMean, self).__init__(
            mptCmdStruct=mptCmdStruct,
            dataType=None,      # Set at runtime
            isDataLayer = True
            )

    # def __init__(self,mptCmdStruct=None):
    
    def _SetFxnDesc(self):
        # description of command used for validation and
        # information display Each Pilot fxn command should have its
        # own description

        self.fxnDesc['DisplayName'] = 'Weighted Mean'
        self.fxnDesc['ShortDesc'] = 'Takes the weighted mean of input variables'
        self.fxnDesc['ReturnType'] = [
            'Float',
            'Positive Float'
            ]                        

        self.fxnDesc['ReqArgs'] = {
            'InFieldNames':['Field Name List'],
            'Weights':['Float','Float List']            
            }
        self.fxnDesc['OptArgs'] = {
            'OutFileName':'File Name',
            'Metadata':'Any',
            'PrecursorFieldNames':['Field Name','Field Name List']
            }
        
    # _SetFxnDesc(self):

    def Exec(self,executedObjects):
        # executedObjects is a dictionary of executed MPilot function objects
        
        self._ValidateListLen('InFldNms',1)
        self._ValidateEqualListLens(['InFieldNames','Weights'])
        
        fldNms = self._ArgToList('InFieldNames')
        wts = self.ValFromArgByNm('Weights')

        # Validate input data and gather types
        for fldNm in fldNms:
            self._ValidateProgDataType(
                fldNm,
                executedObjects[fldNm],
                ['Integer','Positive Integer','Float','Positive Float']
                )
            self._ValidateIsDataLayer(executedObjects[fldNm])
        
        # Result starts with a copy of the first input field.
        self.execRslt = cp.deepcopy(executedObjects[fldNms[0]].ExecRslt()) * wts[0]
        
        for wt,fldNm in zip(wts,fldNms)[1:]:
            self.execRslt += executedObjects[fldNm].ExecRslt() * wt

        self.execRslt = self.execRslt / sum(wts)

        # Determine type of output
        if self.ExecRslt().min() > 0:
            self.dataType = 'Positive Float'
        else:
            self.dataType = 'Float'

        executedObjects[self.RsltNm()] = self
        
    # def Exec(self,executedObjects):

# class WeightedMean(mpefp._MPilotEEMSFxnParent):
        
class Normalize(mpefp._MPilotEEMSFxnParent):

    def __init__(self,mptCmdStruct=None):

        super(Normalize, self).__init__(
            mptCmdStruct=mptCmdStruct,
            dataType=None,      # Set at runtime
            isDataLayer = None  # Set at runtime
            )

    # def __init__(self,mptCmdStruct=None):

    def _SetFxnDesc(self):
        
        self.fxnDesc['DisplayName'] = 'Normalize'
        self.fxnDesc['ShortDesc'] = 'Normalizes the data from another field to range (default 0:1)'
        self.fxnDesc['ReturnType'] = [
            'Float',
            'Positive Float'
            ]                        
        
        self.fxnDesc['ReqArgs'] = {
            'InFieldName':'Field Name'
            }
        self.fxnDesc['OptArgs'] = {
            'StartVal':'Float',
            'EndVal':'Float',
            'OutFileName':'File Name',
            'Metadata':'Any',
            'PrecursorFieldNames':['Field Name','Field Name List']
            }
        
    # _SetFxnDesc(self):

    def DependencyNms(self):
        
        rtrn = self._ArgToList('InFieldName')
        rtrn += self._ArgToList('PrecursorFieldNames')
        return rtrn
    
    # def DependencyNms(self):
        
    def Exec(self,executedObjects):

        inObj = executedObjects[self.ValFromArgByNm('InFieldName')]
        self._ValidateIsDataLayer(inObj)
        
        startVal = self.ValFromArgByNm('StartVal')
        endVal = self.ValFromArgByNm('EndVal')
        if startVal is None: startVal = 0.0
        if endVal is None: endVal = 1.0

        inArr = inObj.ExecRslt()
        inMin = inArr.min()
        inMax = inArr.max()
        
        self.execRslt = (inArr-inMin) * (startVal-endVal) / (inMin-inMax) + startVal

        executedObjects[self.RsltNm()] = self
        
    # def Exec(self,executedObjects):

# class Normalize(mpefp._MPilotEEMSFxnParent):

class PrintVars(mpefp._MPilotEEMSFxnParent):

    def __init__(self,mptCmdStruct=None):

        super(PrintVars, self).__init__(
            mptCmdStruct=mptCmdStruct,
            dataType='None',      # Set at runtime
            isDataLayer = False    # Set at runtime
            )

    # def __init__(self,mptCmdStruct=None):

    def _SetFxnDesc(self):
        # description of command for validation and info display

        self.fxnDesc['DisplayName'] = 'Print variables(s) to screen or file'
        self.fxnDesc['ShortDesc'] = 'Prints each variable in a list of variable names.'
        self.fxnDesc['ReturnType'] = 'Bool'

        self.fxnDesc['ReqArgs'] = {
            'InFieldNames':['Field Name','Field Name List']
            }
        
        self.fxnDesc['OptArgs'] = {
            'Metadata':'Any',
            'PrecursorFieldNames':['Field Name','Field Name List'],
            'OutFileName':'File Name'
            }
        
    # _SetFxnDesc(self):

    def DependencyNms(self):
        
        rtrn = self._ArgToList('InFieldNames')
        rtrn += self._ArgToList('PrecursorFieldNames')
        return rtrn
    
    # def DependencyNms(self):

    def Exec(self,executedObjects):
        # executedObjects is a dictionary of executed MPilot function objects

        fldNms = self._ArgToList('InFieldNames')

        outFNm = self.ArgByNm('OutFileName')
        if outFNm is not None:
            with open(outFNm,'w') as outF:
                for fldNm in fldNms:
                    outF.write('{}: {}\n'.format(
                        fldNm,
                        executedObjects[fldNm].ExecRslt()
                        )
                    )
        else:
            for fldNm in fldNms:
                print '{}: {}'.format(
                    fldNm,
                    executedObjects[fldNm].ExecRslt()
                    )
        
        # if outFNm is not None:
        
        self.execRslt = True
        executedObjects[self.RsltNm()] = self
        
    # def Exec(self,executedObjects):
    
# class PrintVars(mpefp._MPilotEEMSFxnParent):
