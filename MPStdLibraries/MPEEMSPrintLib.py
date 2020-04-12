# Print functions for MPilot
#
# File Log:
# 2017.06.19 - tjs
#   Just a function for printing FieldName value pairs
#

import MPilotEEMSFxnParent as mpefp

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
    
