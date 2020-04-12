from __future__ import division
from MPCore import MPilotEEMSFxnParent as mpefp
import numpy as np
import copy as cp

# Conversion to fuzzy operations

class CvtToFuzzy(mpefp._MPilotEEMSFxnParent):
    
    def __init__(self,mptCmdStruct=None):
        
        super(CvtToFuzzy, self).__init__(
            mptCmdStruct=mptCmdStruct,
            dataType='Fuzzy',
            isDataLayer = True
            )

    # def __init__(self,mptCmdStruct=None):

    def _SetFxnDesc(self):
        # description of command for validation and info display

        self.fxnDesc['DisplayName'] = 'Convert to Fuzzy'
        self.fxnDesc['ShortDesc'] = 'Converts input values into fuzzy values using linear interpolation'
        self.fxnDesc['ReturnType'] = 'Fuzzy'

        self.fxnDesc['ReqArgs'] = {
            'InFieldName':'Field Name'
            }
        self.fxnDesc['OptArgs'] = {
            'TrueThreshold':'Float',
            'FalseThreshold':'Float',
            'OutFileName':'File Name',
            'Metadata':'Any',
            'PrecursorFieldNames':['Field Name','Field Name List'],
            'Direction':'Any',
            }
        
    # _SetFxnDesc(self):

    def DependencyNms(self):
        
        rtrn = self._ArgToList('InFieldName')
        rtrn += self._ArgToList('PrecursorFieldNames')
        return rtrn
    
    # def DependencyNms(self):

    def Exec(self,executedObjects):
        # executedObjects is a dictionary of executed MPilot function objects
        
        inFldNm = self.ArgByNm('InFieldName')
        inObj = executedObjects[inFldNm]
        
        self._ValidateProgDataType(
            inFldNm,
            inObj,
            [
                'Float',
                'Positive Float',
                'Integer',
                'Positive Integer',
                'Fuzzy'
                ]
            )

        self._ValidateIsDataLayer(inObj)
        
        inArr = inObj.ExecRslt()

        if 'Direction' in self.Args():
            if self.ArgByNm('Direction') not in ['LowToHigh','HighToLow']:
                raise Exception(
                    '{}{}{}{}{}{}{}'.format(
                        '\n********************ERROR********************\n',
                        'Invalid Direction specified in command:\n',
                        '  Direction must be one of:\n',
                        '    LowToHigh\n',
                        '    HighToLow\n',
                        'File: {}  Line number: {}\n'.format(
                            self.mptCmdStruct['cmdFileNm'],
                            self.mptCmdStruct['lineNo']
                            ),
                        'Full command:\n{}\n'.format(self.mptCmdStruct['rawCmdStr'])
                    ),
                )

        # if 'Direction' in self.Args():

        if 'FalseThreshold' in self.Args():
            falseThresh = self.ArgByNm('FalseThreshold')
        elif 'Direction' in self.Args():
            if self.ArgByNm('Direction') in ['LowToHigh']:
                falseThresh =  inArr.min()
            elif self.ArgByNm('Direction') in ['HighToLow']:
                falseThresh =  inArr.max()
                
        if 'TrueThreshold' in self.Args():
            trueThresh = self.ArgByNm('TrueThreshold')
        elif 'Direction' in self.Args():
            if self.ArgByNm('Direction') in ['LowToHigh']:
                trueThresh =  inArr.max()
            elif self.ArgByNm('Direction') in ['HighToLow']:
                trueThresh =  inArr.min()

        if trueThresh == falseThresh:
            raise Exception(
                '{}{}{}{}'.format(
                    '\n********************ERROR********************\n',
                    'True and False thresholds must not be equal:\n',
                    'File: {}  Line number: {}\n'.format(
                        self.mptCmdStruct['cmdFileNm'],
                        self.mptCmdStruct['lineNo']
                        ),
                    'Full command:\n{}\n'.format(self.mptCmdStruct['rawCmdStr'])
                ),
            )

        x1 = float(trueThresh)
        x2 = float(falseThresh)
        y1 = self.fuzzyMax
        y2 = self.fuzzyMin

        # linear conversion
        self.execRslt = (inArr - x1) * (y2-y1)/(x2-x1) + y1
        
        # bring back values to fuzzy limits
        self._InsureFuzzy(self.execRslt)

        executedObjects[self.RsltNm()] = self
        
    # def Exec(self,executedObjects):
    
# class CvtToFuzzy(mpefp._MPilotEEMSFxnParent):

class CvtToFuzzyZScore(mpefp._MPilotEEMSFxnParent):
    
    def __init__(self,mptCmdStruct=None):
        
        super(CvtToFuzzyZScore, self).__init__(
            mptCmdStruct=mptCmdStruct,
            dataType='Fuzzy',
            isDataLayer = True
            )

    # def __init__(self,mptCmdStruct=None):

    def _SetFxnDesc(self):
        # description of command for validation and info display

        self.fxnDesc['DisplayName'] = 'Convert to Fuzzy'
        self.fxnDesc['ShortDesc'] = 'Converts input values into fuzzy values using linear interpolation based on Z Score'
        self.fxnDesc['ReturnType'] = 'Fuzzy'

        self.fxnDesc['ReqArgs'] = {
            'InFieldName':'Field Name',
            'TrueThresholdZScore':'Float',
            'FalseThresholdZScore':'Float',
            }
        self.fxnDesc['OptArgs'] = {
            'OutFileName':'File Name',
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
        
        inFldNm = self.ArgByNm('InFieldName')
        inObj = executedObjects[inFldNm]
        
        self._ValidateProgDataType(
            inFldNm,
            inObj,
            [
                'Float',
                'Positive Float',
                'Integer',
                'Positive Integer',
                'Fuzzy'
                ]
            )

        self._ValidateIsDataLayer(inObj)
        
        inArr = inObj.ExecRslt()

        trueThreshZScore = float(self.ArgByNm('TrueThresholdZScore'))
        falseThreshZScore = float(self.ArgByNm('FalseThresholdZScore'))

        if trueThreshZScore == falseThreshZScore:
            raise Exception(
                '{}{}{}{}'.format(
                    '\n********************ERROR********************\n',
                    'True and False ZScore thresholds must not be equal:\n',
                    'File: {}  Line number: {}\n'.format(
                        self.mptCmdStruct['cmdFileNm'],
                        self.mptCmdStruct['lineNo']
                        ),
                    'Full command:\n{}\n'.format(self.mptCmdStruct['rawCmdStr'])
                ),
            )

        rawMean = np.ma.mean(inArr)
        rawStdDev = np.ma.std(inArr)
        
        x1 = rawMean + rawStdDev * trueThreshZScore
        x2 = rawMean + rawStdDev * falseThreshZScore
        
        y1 = self.fuzzyMax
        y2 = self.fuzzyMin

        # linear conversion
        self.execRslt = (inArr - x1) * (y2-y1)/(x2-x1) + y1
        
        # bring back values to fuzzy limits
        self._InsureFuzzy(self.execRslt)

        executedObjects[self.RsltNm()] = self
        
    # def Exec(self,executedObjects):
    
# class CvtToFuzzyZScore(mpefp._MPilotEEMSFxnParent):

class CvtToFuzzyCat(mpefp._MPilotEEMSFxnParent):
    
    def __init__(self,mptCmdStruct=None):
        
        super(CvtToFuzzyCat, self).__init__(
            mptCmdStruct=mptCmdStruct,
            dataType='Fuzzy',
            isDataLayer = True
            )

    # def __init__(self,mptCmdStruct=None):

    def _SetFxnDesc(self):
        # description of command for validation and info display

        self.fxnDesc['DisplayName'] = 'Convert to Fuzzy by Category'
        self.fxnDesc['ShortDesc'] = 'Converts integer input values into fuzzy based on user specification.'
        self.fxnDesc['ReturnType'] = 'Fuzzy'

        self.fxnDesc['ReqArgs'] = {
            'InFieldName':'Field Name',
            'RawValues':['Integer List'],
            'FuzzyValues':['Fuzzy Value List'],
            'DefaultFuzzyValue':'Fuzzy Value',
            }
        self.fxnDesc['OptArgs'] = {
            'OutFileName':'File Name',
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
        
        inFldNm = self.ArgByNm('InFieldName')
        inObj = executedObjects[inFldNm]
        
        # Result starts with a copy of the first input field, then add the rest
        # and divide by the number of them

        self._ValidateIsDataLayer(inObj)
        self._ValidateProgDataType(
            inFldNm,
            inObj,
            [
                'Integer',
                'Positive Integer',
                ]
            )

        self._ValidateEqualListLens(['FuzzyValues','RawValues'])
        self._ValidateArgListItemsUnique('RawValues')
            
        inArr = inObj.ExecRslt()
        
        self.execRslt = np.ma.empty(
            inArr.shape,
            dtype=float
            )

        self.execRslt.data[:] = float(self.ArgByNm('DefaultFuzzyValue'))

        fuzzyVals = self.ValFromArgByNm('FuzzyValues')
        rawVals = self.ValFromArgByNm('RawValues')
        
        for rawVal,fuzzyVal in zip(rawVals,fuzzyVals):
            np.place(self.execRslt.data,inArr.data == rawVal,fuzzyVal)

        self.execRslt.mask = cp.deepcopy(inArr.mask)

        # bring back values to fuzzy limits
        self._InsureFuzzy(self.execRslt)

        executedObjects[self.RsltNm()] = self
        
    # def Exec(self,executedObjects):
    
# class CvtToFuzzyCat(mpefp._MPilotEEMSFxnParent):

class CvtToFuzzyCurve(mpefp._MPilotEEMSFxnParent):
   
    def __init__(self,mptCmdStruct=None):
        
        super(CvtToFuzzyCurve, self).__init__(
            mptCmdStruct=mptCmdStruct,
            dataType='Fuzzy',
            isDataLayer = True
            )

    # def __init__(self,mptCmdStruct=None):

    def _SetFxnDesc(self):
        # description of command for validation and info display

        self.fxnDesc['DisplayName'] = 'Convert to Fuzzy Curve'
        self.fxnDesc['ShortDesc'] = 'Converts input values into fuzzy based on user-defined curve.'
        self.fxnDesc['ReturnType'] = 'Fuzzy'

        self.fxnDesc['ReqArgs'] = {
            'InFieldName':'Field Name',
            'RawValues':['Float List'],
            'FuzzyValues':['Fuzzy Value List'],
            }
        self.fxnDesc['OptArgs'] = {
            'OutFileName':'File Name',
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
        
        inFldNm = self.ArgByNm('InFieldName')
        inObj = executedObjects[inFldNm]
        
        # Result starts with a copy of the first input field, then add the rest
        # and divide by the number of them
        
        self._ValidateIsDataLayer(inObj)
        self._ValidateProgDataType(
            inFldNm,
            inObj,
            [
                'Float',
                'Positive Float',
                'Integer',
                'Positive Integer',
                'Fuzzy'
            ]
            )
        self._ValidateEqualListLens(['FuzzyValues','RawValues'])
        self._ValidateArgListItemsUnique('RawValues')
            
        inArr = inObj.ExecRslt()

        self.execRslt = np.ma.empty(
            inArr.shape,
            dtype=float
            )

        fuzzyVals = self.ValFromArgByNm('FuzzyValues')
        rawVals = self.ValFromArgByNm('RawValues')

        zippedPoints = sorted(zip(rawVals,fuzzyVals))

        # Set the fuzzy values corresponding to the raw values less
        # than the lowest raw value to the corresponding fuzzy value
        np.place(self.execRslt.data,inArr <= zippedPoints[0][0],zippedPoints[0][1])

        # Iterate over the line segments that approximate the curve
        # and assign fuzzy values.
        for ndx in range(1,len(zippedPoints)):

            # Linear equation formula for line segment
            m = (zippedPoints[ndx][1] - zippedPoints[ndx-1][1]) / \
            (zippedPoints[ndx][0] - zippedPoints[ndx-1][0])
            
            b = zippedPoints[ndx-1][1] - m * zippedPoints[ndx-1][0]

            whereNdxs = np.where(
                np.logical_and(
                    inArr.data > zippedPoints[ndx-1][0],
                    inArr.data <= zippedPoints[ndx][0]
                    )
                )

            self.execRslt.data[whereNdxs] = m * inArr.data[whereNdxs] + b
            
        # Set the fuzzy values corresponding to the raw values greater
        # than the highest raw value to the corresponding fuzzy value
        np.place(self.execRslt.data,inArr > zippedPoints[-1][0],zippedPoints[-1][1])
        
        self.execRslt.mask = cp.deepcopy(inArr.mask)

        # bring back values to fuzzy limits
        self._InsureFuzzy(self.execRslt)

        executedObjects[self.RsltNm()] = self
        
    # def Exec(self,executedObjects):

# class CvtToFuzzyCurve(mpefp._MPilotEEMSFxnParent):

class CvtToFuzzyCurveZScore(mpefp._MPilotEEMSFxnParent):
   
    def __init__(self,mptCmdStruct=None):
        
        super(CvtToFuzzyCurveZScore, self).__init__(
            mptCmdStruct=mptCmdStruct,
            dataType='Fuzzy',
            isDataLayer = True
            )

    # def __init__(self,mptCmdStruct=None):

    def _SetFxnDesc(self):
        # description of command for validation and info display

        self.fxnDesc['DisplayName'] = 'Convert to Fuzzy Curve'
        self.fxnDesc['ShortDesc'] = 'Converts input values into fuzzy based on user-defined curve.'
        self.fxnDesc['ReturnType'] = 'Fuzzy'

        self.fxnDesc['ReqArgs'] = {
            'InFieldName':'Field Name',
            'ZScoreValues':['Float List'],
            'FuzzyValues':['Fuzzy Value List'],
            }
        self.fxnDesc['OptArgs'] = {
            'OutFileName':'File Name',
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
        
        inFldNm = self.ArgByNm('InFieldName')
        inObj = executedObjects[inFldNm]
        
        # Result starts with a copy of the first input field, then add the rest
        # and divide by the number of them
        
        self._ValidateIsDataLayer(inObj)
        self._ValidateProgDataType(
            inFldNm,
            inObj,
            [
                'Float',
                'Positive Float',
                'Integer',
                'Positive Integer',
                'Fuzzy'
            ]
            )
        self._ValidateEqualListLens(['FuzzyValues','ZScoreValues'])
        self._ValidateArgListItemsUnique('RawValues')
            
        inArr = inObj.ExecRslt()

        self.execRslt = np.ma.empty(
            inArr.shape,
            dtype=float
            )

        fuzzyVals = self.ValFromArgByNm('FuzzyValues')
        zScoreVals = self.ValFromArgByNm('ZScoreValues')
        rawMean = np.ma.mean(inArr)
        rawStdDev = np.ma.std(inArr)
        
        rawVals = [rawMean + zScoreVal * rawStdDev for zScoreVal in zScoreVals]

        zippedPoints = sorted(zip(rawVals,fuzzyVals))

        # Set the fuzzy values corresponding to the raw values less
        # than the lowest raw value to the corresponding fuzzy value
        np.place(self.execRslt.data,inArr <= zippedPoints[0][0],zippedPoints[0][1])

        # Iterate over the line segments that approximate the curve
        # and assign fuzzy values.
        for ndx in range(1,len(zippedPoints)):

            # Linear equation formula for line segment
            m = (zippedPoints[ndx][1] - zippedPoints[ndx-1][1]) / \
            (zippedPoints[ndx][0] - zippedPoints[ndx-1][0])
            
            b = zippedPoints[ndx-1][1] - m * zippedPoints[ndx-1][0]

            whereNdxs = np.where(
                np.logical_and(
                    inArr.data > zippedPoints[ndx-1][0],
                    inArr.data <= zippedPoints[ndx][0]
                    )
                )

            self.execRslt.data[whereNdxs] = m * inArr.data[whereNdxs] + b
            
        # Set the fuzzy values corresponding to the raw values greater
        # than the highest raw value to the corresponding fuzzy value
        np.place(self.execRslt.data,inArr > zippedPoints[-1][0],zippedPoints[-1][1])
        
        self.execRslt.mask = cp.deepcopy(inArr.mask)

        # bring back values to fuzzy limits
        self._InsureFuzzy(self.execRslt)

        executedObjects[self.RsltNm()] = self
        
    # def Exec(self,executedObjects):

# class CvtToFuzzyCurveZScore(mpefp._MPilotEEMSFxnParent):

class CvtToBinary(mpefp._MPilotEEMSFxnParent):
    
    def __init__(self,mptCmdStruct=None):
        
        super(type(self), self).__init__(
            mptCmdStruct=mptCmdStruct,
            dataType='Fuzzy',
            isDataLayer = True
            )

    # def __init__(self,mptCmdStruct=None):

    def _SetFxnDesc(self):
        # description of command for validation and info display

        self.fxnDesc['DisplayName'] = 'Convert to Fuzzy'
        self.fxnDesc['ShortDesc'] = '''Converts input values into binary 0 or 1 based on threshold.
Direction = LowToHigh for values below threshold to be false and above to be true.
Direction = HighToLow for values below threshold to be true and above to be false.'''
        self.fxnDesc['ReturnType'] = 'Fuzzy'

        self.fxnDesc['ReqArgs'] = {
            'InFieldName':'Field Name',
            'Threshold':'Float',
            'Direction':'Any',
            }
        self.fxnDesc['OptArgs'] = {
            'OutFileName':'File Name',
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
        
        inFldNm = self.ArgByNm('InFieldName')
        inObj = executedObjects[inFldNm]
        
        self._ValidateProgDataType(
            inFldNm,
            inObj,
            [
                'Float',
                'Positive Float',
                'Integer',
                'Positive Integer',
                'Fuzzy'
                ]
            )

        self._ValidateIsDataLayer(inObj)
        
        inArr = inObj.ExecRslt()

        if self.ArgByNm('Direction') not in ['LowToHigh','HighToLow']:
            raise Exception(
                '{}{}{}{}{}{}{}'.format(
                    '\n********************ERROR********************\n',
                    'Invalid Direction specified in command:\n',
                    '  Direction must be one of:\n',
                    '    LowToHigh\n',
                    '    HighToLow\n',
                    'File: {}  Line number: {}\n'.format(
                        self.mptCmdStruct['cmdFileNm'],
                        self.mptCmdStruct['lineNo']
                        ),
                    'Full command:\n{}\n'.format(self.mptCmdStruct['rawCmdStr'])
                ),
            )

        # Convert to binary
        if self.ArgByNm('Direction') == 'LowToHigh':
            lowVal = 0.
            hiVal = 1.
        else:
            lowVal = 1.
            hiVal = 0.
            
        self.execRslt = np.ma.where(
            inArr < self.ValFromArgByNm('Threshold'),
            lowVal,
            hiVal
            )
            
        # bring back values to fuzzy limits
        self._InsureFuzzy(self.execRslt)

        executedObjects[self.RsltNm()] = self
        
    # def Exec(self,executedObjects):
    
# class CvtToBinary(mpefp._MPilotEEMSFxnParent):

# Fuzzy logic operations

class FuzzyUnion(mpefp._MPilotEEMSFxnParent):

    def __init__(self,mptCmdStruct=None):

        super(FuzzyUnion, self).__init__(
            mptCmdStruct=mptCmdStruct,
            dataType='Fuzzy',
            isDataLayer = True
            )

    # def __init__(self,mptCmdStruct=None):
    
    def _SetFxnDesc(self):
        # description of command used for validation and
        # information display Each Pilot fxn command should have its
        # own description

        self.fxnDesc['DisplayName'] = 'Fuzzy Union'
        self.fxnDesc['ShortDesc'] = 'Takes the fuzzy Union (mean) of fuzzy input variables'
        self.fxnDesc['ReturnType'] = 'Fuzzy'

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
        # executedObjects is a dictionary of executed MPilot function objects
        
        self._ValidateListLen('InFldNms',1)
        fldNms = self._ArgToList('InFieldNames')
        
        # Result starts with a copy of the first input field, then add the rest
        # and divide by the number of them
        
        self._ValidateProgDataType(fldNms[0],executedObjects[fldNms[0]],'Fuzzy')
        self.execRslt = cp.deepcopy(executedObjects[fldNms[0]].ExecRslt())
        
        for fldNm in fldNms[1:]:
            self._ValidateProgDataType(fldNm,executedObjects[fldNm],'Fuzzy')
            self.execRslt += executedObjects[fldNm].ExecRslt()

        self.execRslt = self.execRslt / float(len(fldNms))
        self._InsureFuzzy(self.execRslt)

        executedObjects[self.RsltNm()] = self
        
    # def Exec(self,executedObjects):
    
# class FuzzyUnion(mpefp._MPilotEEMSFxnParent):    

class FuzzyWeightedUnion(mpefp._MPilotEEMSFxnParent):

    def __init__(self,mptCmdStruct=None):

        super(FuzzyWeightedUnion, self).__init__(
            mptCmdStruct=mptCmdStruct,
            dataType='Fuzzy',
            isDataLayer = True
            )

    # def __init__(self,mptCmdStruct=None):
    
    def _SetFxnDesc(self):
        # description of command used for validation and
        # information display Each Pilot fxn command should have its
        # own description

        self.fxnDesc['DisplayName'] = 'Fuzzy Weighted Union'
        self.fxnDesc['ShortDesc'] = 'Takes the weighted fuzzy Union (mean) of fuzzy input variables'
        self.fxnDesc['ReturnType'] = 'Fuzzy'

        self.fxnDesc['ReqArgs'] = {
            'InFieldNames':['Field Name List'],
            'Weights':['Float List']            
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
        
        # Result starts with a copy of the first input field, then add the rest
        # and divide by the number of them
        
        self._ValidateProgDataType(fldNms[0],executedObjects[fldNms[0]],'Fuzzy')
        self.execRslt = cp.deepcopy(executedObjects[fldNms[0]].ExecRslt()) * wts[0]
        
        for wt,fldNm in zip(wts,fldNms)[1:]:
            self._ValidateProgDataType(fldNm,executedObjects[fldNm],'Fuzzy')
            self.execRslt += executedObjects[fldNm].ExecRslt() * wt

        self.execRslt = self.execRslt / sum(wts)
        self._InsureFuzzy(self.execRslt)

        executedObjects[self.RsltNm()] = self
        
    # def Exec(self,executedObjects):
    
# class FuzzyWeightedUnion(mpefp._MPilotEEMSFxnParent):    

class FuzzySelectedUnion(mpefp._MPilotEEMSFxnParent):

    def __init__(self,mptCmdStruct=None):

        super(FuzzySelectedUnion, self).__init__(
            mptCmdStruct=mptCmdStruct,
            dataType='Fuzzy',
            isDataLayer = True
            )

    # def __init__(self,mptCmdStruct=None):
    
    def _SetFxnDesc(self):
        # description of command used for validation and
        # information display Each Pilot fxn command should have its
        # own description

        self.fxnDesc['DisplayName'] = 'Fuzzy Selected Union'
        self.fxnDesc['ShortDesc'] = 'Takes the fuzzy Union (mean) of N Truest or Falsest fuzzy input variables'
        self.fxnDesc['ReturnType'] = 'Fuzzy'

        self.fxnDesc['ReqArgs'] = {
            'InFieldNames':['Field Name','Field Name List'],
            'TruestOrFalsest':'Truest Or Falsest',
            'NumberToConsider':'Positive Integer'
            }
        self.fxnDesc['OptArgs'] = {
            'OutFileName':'File Name',
            'Metadata':'Any',
            'PrecursorFieldNames':['Field Name','Field Name List']
            }
        
    # _SetFxnDesc(self):

    def Exec(self,executedObjects):
        # executedObjects is a dictionary of executed MPilot function objects
        
        fldNms = self._ArgToList('InFieldNames')
        numToCnsdr = self.ValFromArgByNm('NumberToConsider')
        
        if len(fldNms) < numToCnsdr:
            raise Exception(
                '{}{}{}{}{}'.format(
                    '\n********************ERROR********************\n',
                    'Number of InFieldNames must be greater than or equal to NumberToConsider:\n',
                    'Number of InFieldNames: {}  NumberToConsider: {}\n'.format(
                        len(fldNms),
                        numToCnsdr,
                        ),
                    'File: {}  Line number: {}\n'.format(
                        self.mptCmdStruct['cmdFileNm'],
                        self.mptCmdStruct['lineNo']
                        ),
                    'Full command:\n{}\n'.format(self.mptCmdStruct['rawCmdStr'])
                ),
            )

        # Create a stacked array with layers from input arrays, sort it, and use
        # that to calculate fuzzy xor. There is no np.ma.stacked, so the masks have
        # to be handled separately from the data. Note we are building the maximal
        # mask from all the inputs before broadcasting it to the size of the stacked
        # array. There are some issues with getting stackedArr to be writable. The
        # below code works.

        tmpMask = executedObjects[fldNms[0]].ExecRslt().mask
        for fldNm in fldNms[1:]:
            tmpMask = np.logical_or(tmpMask,executedObjects[fldNm].ExecRslt().mask)

        stackedArr = np.ma.array(
            np.stack([executedObjects[fldNm].ExecRslt().data for fldNm in fldNms]),
            mask = cp.deepcopy(
                np.broadcast_to(
                    tmpMask,
                    [len(fldNms)]+list(executedObjects[fldNm].ExecRslt().shape)
                    )
                )
            )
            
        stackedArr.sort(axis=0,kind='heapsort')

        if self.ValFromArgByNm('TruestOrFalsest') == 'Truest':
            self.execRslt = np.ma.mean(stackedArr[-numToCnsdr:],axis=0)
        else:
            self.execRslt = np.ma.mean(stackedArr[0:numToCnsdr],axis=0)
            
        executedObjects[self.RsltNm()] = self
        
    # def Exec(self,executedObjects):
    
# class FuzzySelectedUnion(mpefp._MPilotEEMSFxnParent):

class FuzzyOr(mpefp._MPilotEEMSFxnParent):

    def __init__(self,mptCmdStruct=None):

        super(FuzzyOr, self).__init__(
            mptCmdStruct=mptCmdStruct,
            dataType='Fuzzy',
            isDataLayer = True
            )

    # def __init__(self,mptCmdStruct=None):
    
    def _SetFxnDesc(self):

        self.fxnDesc['DisplayName'] = 'Fuzzy Or'
        self.fxnDesc['ShortDesc'] = 'Takes the fuzzy Or (maximum) of fuzzy input variables'
        self.fxnDesc['ReturnType'] = 'Fuzzy'
        
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
        
        # Traverse inputs and take maximum
        self._ValidateProgDataType(fldNms[0],executedObjects[fldNms[0]],'Fuzzy')
        self.execRslt = cp.deepcopy(executedObjects[fldNms[0]].ExecRslt())
        
        for fldNm in fldNms[1:]:
            self._ValidateProgDataType(fldNm,executedObjects[fldNm],'Fuzzy')
            self.execRslt = np.ma.maximum(self.execRslt,executedObjects[fldNm].ExecRslt())

        self.execRslt = self.execRslt
        self._InsureFuzzy(self.execRslt)

        executedObjects[self.RsltNm()] = self
        
    # def Exec(self,executedObjects):

# class FuzzyOr(mpefp._MPilotEEMSFxnParent):

class FuzzyAnd(mpefp._MPilotEEMSFxnParent):

    def __init__(self,mptCmdStruct=None):

        super(FuzzyAnd, self).__init__(
            mptCmdStruct=mptCmdStruct,
            dataType='Fuzzy',
            isDataLayer = True
            )

    # def __init__(self,mptCmdStruct=None):

    def _SetFxnDesc(self):
        
        self.fxnDesc['DisplayName'] = 'Fuzzy And'
        self.fxnDesc['ShortDesc'] = 'Takes the fuzzy And (minimum) of fuzzy input variables'
        self.fxnDesc['ReturnType'] = 'Fuzzy'
        
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
        
        # Traverse inputs and take maximum
        self._ValidateProgDataType(fldNms[0],executedObjects[fldNms[0]],'Fuzzy')
        self.execRslt = cp.deepcopy(executedObjects[fldNms[0]].ExecRslt())

        for fldNm in fldNms[1:]:
            self._ValidateProgDataType(fldNm,executedObjects[fldNm],'Fuzzy')
            self.execRslt = np.ma.minimum(self.execRslt,executedObjects[fldNm].ExecRslt())

        self._InsureFuzzy(self.execRslt)

        executedObjects[self.RsltNm()] = self
        
    # def Exec(self,executedObjects):
    
# class FuzzyAnd(mpefp._MPilotEEMSFxnParent):

class FuzzyXOr(mpefp._MPilotEEMSFxnParent):

    def __init__(self,mptCmdStruct=None):

        super(FuzzyXOr, self).__init__(
            mptCmdStruct=mptCmdStruct,
            dataType='Fuzzy',
            isDataLayer = True
            )

    # def __init__(self,mptCmdStruct=None):

    def _SetFxnDesc(self):
        
        self.fxnDesc['DisplayName'] = 'Fuzzy XOr'
        self.fxnDesc['ShortDesc'] = 'Computes Fuzzy XOr: Truest - (Truest - 2nd Truest) * (2nd Truest - full False)/(Truest - full False)'
        self.fxnDesc['ReturnType'] = 'Fuzzy'
        
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

        self._ValidateListLen('InFieldNames',2)
        fldNms = self._ArgToList('InFieldNames')
        
        # Validate inputs
        for fldNm in fldNms:
            self._ValidateProgDataType(fldNm,executedObjects[fldNm],'Fuzzy')
        
        # Create a stacked array with layers from input arrays, sort it, and use
        # that to calculate fuzzy xor. There is no np.ma.stacked, so the masks have
        # to be handled separately from the data. Note we are building the maximal
        # mask from all the inputs before broadcasting it to the size of the stacked
        # array. There are some issues with getting stackedArr to be writable. The
        # below code works.

        tmpMask = executedObjects[fldNms[0]].ExecRslt().mask
        for fldNm in fldNms[1:]:
            tmpMask = np.logical_or(tmpMask,executedObjects[fldNm].ExecRslt().mask)
            
        stackedArr = np.ma.array(
            np.stack([executedObjects[fldNm].ExecRslt().data for fldNm in fldNms]),
            mask = cp.deepcopy(
                np.broadcast_to(
                    tmpMask,
                    [len(fldNms)]+list(executedObjects[fldNm].ExecRslt().shape)
                    )
                )
            )
            
        stackedArr.sort(axis=0,kind='heapsort')
        
        self.execRslt = np.ma.where(
            stackedArr[-1] <= self.fuzzyMin,
            self.fuzzyMin,
            stackedArr[-1] - \
                (stackedArr[-1] - stackedArr[-2]) * \
                (stackedArr[-2] - self.fuzzyMin) / \
                (stackedArr[-1] - self.fuzzyMin)
            )

        self._InsureFuzzy(self.execRslt)

        executedObjects[self.RsltNm()] = self
        
    # def Exec(self,executedObjects):
    
# class FuzzyXOr(mpefp._MPilotEEMSFxnParent):

class FuzzyNot(mpefp._MPilotEEMSFxnParent):

    def __init__(self,mptCmdStruct=None):

        super(FuzzyNot, self).__init__(
            mptCmdStruct=mptCmdStruct,
            dataType='Fuzzy',
            isDataLayer = True
            )

    # def __init__(self,mptCmdStruct=None):

    def _SetFxnDesc(self):
        
        self.fxnDesc['DisplayName'] = 'Fuzzy Not'
        self.fxnDesc['ShortDesc'] = 'Takes the fuzzy And (minimum) of fuzzy input variables'
        self.fxnDesc['ReturnType'] = 'Fuzzy'
        
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

        inFldNm = self.ArgByNm('InFieldName')

        self._ValidateProgDataType(inFldNm,executedObjects[inFldNm],'Fuzzy')
        self.execRslt = -executedObjects[inFldNm].ExecRslt()

        self._InsureFuzzy(self.execRslt)

        executedObjects[self.RsltNm()] = self

    # def Exec(self,executedObjects):

# class FuzzyNot(mpefp._MPilotEEMSFxnParent):    
    
class CvtFromFuzzy(mpefp._MPilotEEMSFxnParent):
    
    def __init__(self,mptCmdStruct=None):
        
        super(type(self), self).__init__(
            mptCmdStruct=mptCmdStruct,
            dataType='Float',
            isDataLayer = True
            )

    # def __init__(self,mptCmdStruct=None):

    def _SetFxnDesc(self):
        # description of command for validation and info display

        self.fxnDesc['DisplayName'] = 'Convert from Fuzzy'
        self.fxnDesc['ShortDesc'] = 'Converts input fuzzy values into non-fuzzy values using linear interpolation'
        self.fxnDesc['ReturnType'] = 'Fuzzy'

        self.fxnDesc['ReqArgs'] = {
            'InFieldName':'Field Name',
            # Required for conversion
            'TrueThreshold':'Float',
            'FalseThreshold':'Float',
            }
        self.fxnDesc['OptArgs'] = {
            'OutFileName':'File Name',
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
        
        inFldNm = self.ArgByNm('InFieldName')
        inObj = executedObjects[inFldNm]
        
        self._ValidateProgDataType(
            inFldNm,
            inObj,
            [
                'Fuzzy'
                ]
            )

        self._ValidateIsDataLayer(inObj)
        
        inArr = inObj.ExecRslt()

        falseThresh = self.ArgByNm('FalseThreshold')
        trueThresh = self.ArgByNm('TrueThreshold')

        if trueThresh == falseThresh:
            raise Exception(
                '{}{}{}{}'.format(
                    '\n********************ERROR********************\n',
                    'True and False thresholds must not be equal:\n',
                    'File: {}  Line number: {}\n'.format(
                        self.mptCmdStruct['cmdFileNm'],
                        self.mptCmdStruct['lineNo']
                        ),
                    'Full command:\n{}\n'.format(self.mptCmdStruct['rawCmdStr'])
                ),
            )

        # Simply reverse the x and y used in CvtToFuzzy
        y1 = float(trueThresh)
        y2 = float(falseThresh)
        x1 = self.fuzzyMax
        x2 = self.fuzzyMin

        # linear conversion
        self.execRslt = (inArr - x1) * (y2-y1)/(x2-x1) + y1

        executedObjects[self.RsltNm()] = self
        
    # def Exec(self,executedObjects):
    
# class CvtFromFuzzy(mpefp._MPilotEEMSFxnParent):
