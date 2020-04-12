# EEMS I/O CSV Command classes for MPilot framework
#
# File Log:
# 2016.07.06 - tjs
#  Created _MPilotEEMSParent

from MPCore import MPilotEEMSFxnParent as mpefp
import numpy as np

class EEMSRead(mpefp._MPilotEEMSFxnParent):

    def __init__(
        self,
        mptCmdStruct=None
        ):

        self.mptCmdStruct = mptCmdStruct
        if not self.ArgExists('DataType'):
            dataType = 'Float'
        elif self.ArgByNm('DataType') == 'Integer':
                dataType = 'Integer'
        elif self.ArgByNm('DataType') == 'Float':
                dataType = 'Integer'
        else:
            raise Exception(
                '{}{}{}{}'.format(
                    '\n********************ERROR********************\n',
                    'DataType was: {}  Should be one of: Float Integer}\n'.format(
                        self.ArgByNm('DataType')
                        ),
                    'File: {}  Line number: {}\n'.format(
                        self.mptCmdStruct['cmdFileNm'],
                        self.mptCmdStruct['lineNo']
                        ),
                    'Full command:\n{}\n'.format(self.mptCmdStruct['rawCmdStr'])
                    ),
                )

        super(EEMSRead, self).__init__(
            mptCmdStruct=mptCmdStruct,
            dataType=dataType,
            isDataLayer=True)

    # def __init__(self,mptCmdStruct=None):
    
    def _SetFxnDesc(self):
        # description of command used for validation and
        # information display Each Pilot fxn command should have its
        # own description

        self.fxnDesc['DisplayName'] = 'EEMSREAD'
        self.fxnDesc['ShortDesc'] = 'Reads a variable from a file'
        self.fxnDesc['ReturnType'] = ['Float','Integer']

        self.fxnDesc['ReqArgs'] = {
            'InFileName':'File Name',
            'InFieldName':'Field Name'
            }
        self.fxnDesc['OptArgs'] = {
            'OutFileName':'File Name',
            'Metadata':'Any',
            'ReturnType':'Data Type Desc',
            'NewFieldName':'Field Name',
            'MissingVal':'Any',
            'PrecursorFieldNames':['Field Name','Field Name List'],
            'DataType':'Data Type Desc',
            }
        
    # _SetFxnDesc(self):

    def _SetRsltNm(self,nm):
        self.mptCmdStruct['parsedCmd']['rsltNm'] = nm

    def DependencyNms(self):
        
        rtrn = self._ArgToList('PrecursorFieldNames')
        return rtrn
    
    # def DependencyNms(self):
    
    def Exec(self,executedObjects):

        # Set the return name
        if self.RsltNm() is None:
            if ArgExists('NewFieldName'):
                _SetRsltNm(self.ArgByNm('NewFieldName'))
            else:
                _SetRsltNm(self.ArgByNm('InFieldName'))
        # if self.RsltNm() is None:
                           
        if self.ArgExists('ReturnType'):
            self.fxnDesc['ReturnType'] = self.ArgByNm('ReturnType')

        with open(self.ArgByNm('InFileName'),'rU') as inF:
            # Read first line, clean it, get index of field name
            try:
                colNdx = inF.readline().strip().replace('"','').split(',').index(self.ArgByNm('InFieldName'))
            except:
                raise Exception(
                    '{}{}{}{}{}{}'.format(
                        '\n********************ERROR********************\n',
                        'Invalid file:{}\n'.format(self.ArgByNm('InFileName')),
                        '  Problem may be an empty file, a file lacking a header line,\n',
                        '  or a field name not present in the file.\n',
                        'File: {}  Line number: {}\n'.format(
                            self.mptCmdStruct['cmdFileNm'],
                            self.mptCmdStruct['lineNo']
                            ),
                        'Full command:\n{}\n'.format(self.mptCmdStruct['rawCmdStr'])
                    ),
                )

            # Read in the data, line by line
            # Skip blank lines
            line = inF.readline()
            inVals = []
            while line != '':
                stripLine = line.strip().replace('"','')
                if stripLine != '':
                    charVal = stripLine.split(',')[colNdx]
                    inVals.append(float(charVal))
                line = inF.readline()

        # with open(self.ArgByNm('InFileName'),'r') as inF:

        # Heere Need to take care of return type and data type
        if self.ArgExists('MissingVal'):
            if self.DataType() == 'Integer':
                self.execRslt = np.ma.array(
                    inVals,mask=False,
                    dtype = np.int,
                    fill_value=int(self.ArgByNm('MissingVal'))
                    )
            else:
                self.execRslt = np.ma.array(
                    inVals,mask=False,
                    dtype = np.float,
                    fill_value=float(self.ArgByNm('MissingVal'))
                    )
        else:
            self.execRslt = np.ma.array(
                inVals,mask=False
                )
            
        # if self.ArgExists('MissingVal'):
        
        self.execRslt.soften_mask()

        if self.ArgExists('MissingVal'):
            self.execRslt.mask = np.ma.where(self.execRslt == float(self.ArgByNm('MissingVal')), True, False)


        # Check shape versus shape of any object with a data layer
        # (i.e. data layer is masked array)

        for exObjKey,exObjVal in executedObjects.items():
            # Skip results that are not data
            if not exObjVal.IsDataLayer():
                continue
            
            if  exObjVal.ExecRslt().shape != self.execRslt.shape:
                raise Exception(
                    '{}{}{}{}{}'.format(
                        '\n********************ERROR********************\n',
                        'Input variable shape did not match existing data:\n',
                        '  Input file name, field name, and shape: {}  {}'.format(
                            self.ArgByNm('InFileName'),
                            self.ArgByNm('InFieldName'),
                            self.execRslt.shape
                            ),
                        '  Existing field name and shape: {}  {}'.format(
                            exObjKey,
                            exObjVal.ExecRslt().shape
                            ),
                        'File: {}  Line number: {}\n'.format(
                            self.mptCmdStruct['cmdFileNm'],
                            self.mptCmdStruct['lineNo']
                            ),
                        'Full command:\n{}\n'.format(self.mptCmdStruct['rawCmdStr'])
                    ),
                )
            
        executedObjects[self.RsltNm()] = self
        
    # def Exec(self,executedObjects):

# class EEMSRead(mpefp._MPilotEEMSFxnParent)        

class EEMSWrite(mpefp._MPilotEEMSFxnParent):

    def __init__(
        self,
        mptCmdStruct=None
        ):

        super(EEMSWrite, self).__init__(
            mptCmdStruct=mptCmdStruct,
            dataType='Boolean',
            isDataLayer=False
            )

    # def __init__(self,mptCmdStruct=None):
    
    def _SetFxnDesc(self):
        # description of command used for validation and
        # information display Each Pilot fxn command should have its
        # own description

        self.fxnDesc['DisplayName'] = 'EEMSWrite'
        self.fxnDesc['ShortDesc'] = 'Writes one or more variables to a file'
        self.fxnDesc['ReturnType'] = 'Boolean'

        self.fxnDesc['ReqArgs'] = {
            'OutFileName':'File Name',
            'OutFieldNames':['Field Name','Field Name List']
            }
        self.fxnDesc['OptArgs'] = {
            'PrecursorFieldNames':['Field Name','Field Name List']
            }
        
    # _SetFxnDesc(self):

    def DependencyNms(self):
        
        # strip square brackets (if any) and split into list
        rtrn = self.ArgByNm('OutFieldNames').replace('[','').replace(']','').split(',')
        # optional dependencies
        if self.ArgExists('PrecursorFieldNames'):
            rtrn += self.ArgByNm('PrecursorFieldNames').replace('[','').replace(']','').split(',')
        return rtrn
    
    # def DependencyNms(self):

    def Exec(self,executedObjects):

        outFldNms = self.ArgByNm('OutFieldNames').replace('[','').replace(']','').split(',')

        # Check that all output variables are array
        for outFldNm in outFldNms:
            if not isinstance(executedObjects[outFldNm].ExecRslt(),np.ma.core.MaskedArray):
                raise Exception(
                    '{}{}{}{}'.format(
                        '\n********************ERROR********************\n',
                        'Trying to output non-array value: {}.\n'.format(outFldNm),
                        'File: {}  Line number: {}\n'.format(
                            self.mptCmdStruct['cmdFileNm'],
                            self.mptCmdStruct['lineNo']
                            ),
                        'Full command:\n{}\n'.format(self.mptCmdStruct['rawCmdStr'])
                    ),
                )

        # Check that all output variables have the same shape

        for outFldNm in outFldNms[1:]:
            if executedObjects[outFldNm].ExecRslt().shape != \
              executedObjects[outFldNms[0]].ExecRslt().shape:
                raise Exception(
                    '{}{}{}{}{}'.format(
                        '\n********************ERROR********************\n',
                        'Output variable shapes do not match.\n',
                        'Variables and shapes: {} {}, {} {}\n'.format(
                            outFldNms[0],
                            executedObjects[outFldNms[0]].ExecRslt().shape,
                            outFldNm,
                            executedObjects[outFldNm].ExecRslt().shape
                            ),
                        'File: {}  Line number: {}\n'.format(
                            self.mptCmdStruct['cmdFileNm'],
                            self.mptCmdStruct['lineNo']
                            ),
                        'Full command:\n{}\n'.format(self.mptCmdStruct['rawCmdStr'])
                    ),
                )

        # Now write those bad boys!
        with open(self.ArgByNm('OutFileName'),'w') as outF:

            outF.write('{}\n'.format(','.join(outFldNms)))

            # Assumption: 1d arrays
            outArr = np.ma.array([executedObjects[outFldNm].ExecRslt() for outFldNm in outFldNms])
            outArr = outArr.transpose([1,0])

            for rowNdx in range(outArr.shape[0]):
                outF.write('{}\n'.format(','.join([str(x) for x in outArr[rowNdx,:]])))
                                         
        # with open(self.ArgByNm('OutFileName'),'w') as outF:
        
        # self.execRslt = True
        # executedObjects[self.RsltNm()] = self
        
    # def Exec(self,executedObjects):

# class EEMSWrite(mpefp._MPilotEEMSFxnParent)        
