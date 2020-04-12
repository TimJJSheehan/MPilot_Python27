# Stats functions for MPilot
#
# File Log:
# 2017.06.19 - tjs
#  Created for R2. This could use some really cool classes
#

from __future__ import division
from MPCore import MPilotEEMSFxnParent as mpefp
import numpy as np
from scipy import stats

class StatsR2(mpefp._MPilotEEMSFxnParent):

    def __init__(self,mptCmdStruct=None):
        
        super(type(self),self).__init__(
            mptCmdStruct=mptCmdStruct,
            isDataLayer = False,
            dataType='Float',
            )
        
    # def __init__(self,mptCmdStruct=None):

    def _SetFxnDesc(self):
        # description of command for validation and info display

        self.fxnDesc['DisplayName'] = 'Statistical R2'
        self.fxnDesc['ShortDesc'] = 'Calculates R2 for two arrays of values'
        self.fxnDesc['ReturnType'] = 'Float'

        self.fxnDesc['ReqArgs'] = {
            'XFieldName':'Field Name',
            'YFieldName':'Field Name'
            }
        
        self.fxnDesc['OptArgs'] = {
            'Metadata':'Any',
            'PrecursorFieldNames':['Field Name','Field Name List'],
            }
        
    # _SetFxnDesc(self):

    def DependencyNms(self):
        
        rtrn = self._ArgToList('XFieldName')
        rtrn += self._ArgToList('YFieldName')
        rtrn += self._ArgToList('PrecursorFieldNames')
        return rtrn
    
    # def DependencyNms(self):

    def Exec(self,executedObjects):
        # executedObjects is a dictionary of executed MPilot function objects

        xObj = executedObjects[self.ValFromArgByNm('XFieldName')]
        yObj = executedObjects[self.ValFromArgByNm('YFieldName')]

        self.execRslt = stats.linregress(xObj.ExecRslt(),yObj.ExecRslt()).rvalue ** 2
        executedObjects[self.RsltNm()] = self

    # def Exec(self,executedObjects):

# class StatsR2(mpefp._MPilotEEMSFxnParent):

class StatsPValue(mpefp._MPilotEEMSFxnParent):

    def __init__(self,mptCmdStruct=None):
        
        super(type(self), self).__init__(
            mptCmdStruct=mptCmdStruct,
            isDataLayer = False,
            dataType='Float',
            )
        
    # def __init__(self,mptCmdStruct=None):

    def _SetFxnDesc(self):
        # description of command for validation and info display

        self.fxnDesc['DisplayName'] = 'Statistical R2'
        self.fxnDesc['ShortDesc'] = 'Calculates R2 for two arrays of values'
        self.fxnDesc['ReturnType'] = 'Float'

        self.fxnDesc['ReqArgs'] = {
            'XFieldName':'Field Name',
            'YFieldName':'Field Name'
            }
        
        self.fxnDesc['OptArgs'] = {
            'Metadata':'Any',
            'PrecursorFieldNames':['Field Name','Field Name List'],
            }
        
    # _SetFxnDesc(self):

    def DependencyNms(self):
        
        rtrn = self._ArgToList('XFieldName')
        rtrn += self._ArgToList('YFieldName')
        rtrn += self._ArgToList('PrecursorFieldNames')
        return rtrn
    
    # def DependencyNms(self):

    def Exec(self,executedObjects):
        # executedObjects is a dictionary of executed MPilot function objects

        xObj = executedObjects[self.ValFromArgByNm('XFieldName')]
        yObj = executedObjects[self.ValFromArgByNm('YFieldName')]

        self.execRslt = stats.linregress(xObj.ExecRslt(),yObj.ExecRslt()).pvalue
        executedObjects[self.RsltNm()] = self

    # def Exec(self,executedObjects):
    
class StatsStdErr(mpefp._MPilotEEMSFxnParent):

    def __init__(self,mptCmdStruct=None):
        
        super(type(self), self).__init__(
            mptCmdStruct=mptCmdStruct,
            isDataLayer = False,
            dataType='Float',
            )
        
    # def __init__(self,mptCmdStruct=None):

    def _SetFxnDesc(self):
        # description of command for validation and info display

        self.fxnDesc['DisplayName'] = 'Statistical R2'
        self.fxnDesc['ShortDesc'] = 'Calculates R2 for two arrays of values'
        self.fxnDesc['ReturnType'] = 'Float'

        self.fxnDesc['ReqArgs'] = {
            'XFieldName':'Field Name',
            'YFieldName':'Field Name'
            }
        
        self.fxnDesc['OptArgs'] = {
            'Metadata':'Any',
            'PrecursorFieldNames':['Field Name','Field Name List'],
            }
        
    # _SetFxnDesc(self):

    def DependencyNms(self):
        
        rtrn = self._ArgToList('XFieldName')
        rtrn += self._ArgToList('YFieldName')
        rtrn += self._ArgToList('PrecursorFieldNames')
        return rtrn
    
    # def DependencyNms(self):

    def Exec(self,executedObjects):
        # executedObjects is a dictionary of executed MPilot function objects

        xObj = executedObjects[self.ValFromArgByNm('XFieldName')]
        yObj = executedObjects[self.ValFromArgByNm('YFieldName')]

        self.execRslt = stats.linregress(xObj.ExecRslt(),yObj.ExecRslt()).stderr
        executedObjects[self.RsltNm()] = self

    # def Exec(self,executedObjects):
    
# class StatsStdErr(mpefp._MPilotEEMSFxnParent):

class StatsSum(mpefp._MPilotEEMSFxnParent):

    def __init__(self,mptCmdStruct=None):
        
        super(type(self), self).__init__(
            mptCmdStruct=mptCmdStruct,
            isDataLayer = False,
            dataType='Float',
            )
        
    # def __init__(self,mptCmdStruct=None):

    def _SetFxnDesc(self):
        # description of command for validation and info display

        self.fxnDesc['DisplayName'] = 'Statistical Sum'
        self.fxnDesc['ShortDesc'] = 'The sum of values in a single array'
        self.fxnDesc['ReturnType'] = 'Float'

        self.fxnDesc['ReqArgs'] = {
            'InFieldName':'Field Name',
            }
        
        self.fxnDesc['OptArgs'] = {
            'Metadata':'Any',
            'PrecursorFieldNames':['Field Name','Field Name List'],
            }
        
    # _SetFxnDesc(self):

    def DependencyNms(self):
        
        rtrn = self._ArgToList('InFieldName')
        rtrn += self._ArgToList('PrecursorFieldNames')
        return rtrn
    
    # def DependencyNms(self):

    def Exec(self,executedObjects):
        # executedObjects is a dictionary of executed MPilot function objects
        
        self.execRslt = executedObjects[self.ValFromArgByNm('InFieldName')].ExecRslt().sum()
        executedObjects[self.RsltNm()] = self

    # def Exec(self,executedObjects):
    
# class StatsSum(mpefp._MPilotEEMSFxnParent):

class MeanPValStdErr(mpefp._MPilotEEMSFxnParent):

    def __init__(self,mptCmdStruct=None):
        
        super(type(self), self).__init__(
            mptCmdStruct=mptCmdStruct,
            isDataLayer = False,
            dataType='Float',
            )
        
    # def __init__(self,mptCmdStruct=None):

    def _SetFxnDesc(self):
        # description of command for validation and info display

        self.fxnDesc['DisplayName'] = 'Statistical R2, p Value, Standard Deviation'
        self.fxnDesc['ShortDesc'] = 'Calculates R2, p Value, Standard Deviation for two arrays of values'
        self.fxnDesc['ReturnType'] = 'String'

        self.fxnDesc['ReqArgs'] = {
            'XFieldName':'Field Name',
            'YFieldName':'Field Name'
            }
        
        self.fxnDesc['OptArgs'] = {
            'Metadata':'Any',
            'PrecursorFieldNames':['Field Name','Field Name List'],
            }
        
    # _SetFxnDesc(self):

    def DependencyNms(self):
        
        rtrn = self._ArgToList('XFieldName')
        rtrn += self._ArgToList('YFieldName')
        rtrn += self._ArgToList('PrecursorFieldNames')
        return rtrn
    
    # def DependencyNms(self):

    def Exec(self,executedObjects):
        # executedObjects is a dictionary of executed MPilot function objects

        xObj = executedObjects[self.ValFromArgByNm('XFieldName')]
        yObj = executedObjects[self.ValFromArgByNm('YFieldName')]

        slopb,intcpt,rVal,pVal,stdErr = stats.linregress(xObj.ExecRslt(),yObj.ExecRslt())
        self.execRslt = 'R2, p Value, Std Err:  {}   {}   {}'.format(rVal**2,pVal,stdErr)
        
        executedObjects[self.RsltNm()] = self

    # def Exec(self,executedObjects):
    
# class MeanPValStdErr(mpefp._MPilotEEMSFxnParent):
