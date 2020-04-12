# MPilot is a framework designed to let programmers create tree-based
# workflow frameworks for users. This is the top level class used in
# the MPilot Framework.
#
# File Log:
# 2016.07.06 - tjs
#  Created _MPilotFxnParent

from collections import OrderedDict
import copy as cp
import MPilotParseMetadata as pmd

class _MPilotFxnParent(object):

    def __init__(self,mptCmdStruct = None):

        # self.mptCmdStruct contains information about the MPilot command
        # that is used to order commands, check command validity,
        # and provide detailed information in error messages.

        # The arguments that are used in executing the command
        
        self.execRslt = None
        self.fxnDesc = OrderedDict()
        self.fxnDesc['Name'] = self.__class__.__name__
        self._SetFxnDesc()

        if mptCmdStruct is None:
            
            self.mptCmdStruct = {
                'cmdFileNm':None,
                'lineNo':None,
                'rawCmdStr':None,
                'cleanCmdStr':None,
                'parsedCmd':None
                }

            self.mptCmdStruct['parsedCmd'] = {}
            self.mptCmdStruct['parsedCmd']['arguments'] = {}
                
        else:
            
            self.mptCmdStruct = cp.deepcopy(mptCmdStruct)
            self._ValidateStrCmd()
            
        # if mptCmdStruct is None:

    # def __init__(...)

    def _SetFxnDesc(self):
        # description of command used for validation and
        # information display Each Pilot fxn command should have its
        # own description

        raise Exception(
            '{}{}{}{}{}{}{}'.format(
                '\n********************ERROR********************\n',
                'Programming error:\n',
                '  Your program is using the inherited _MPilotFxnParent:_SetFxnDesc() method.',
                '  There should be a unique _SetFxnDesc() method for the defined MPilot command.',
                '  Check and correct the class definition of the MPilot command: {}\n'.format(self.fxnDesc['Name']),
                'File: {}  Line number: {}\n'.format(
                    self.mptCmdStruct['cmdFileNm'],self.mptCmdStruct['lineNo']
                    ),
                'Full command:\n{}\n'.format(self.mptCmdStruct['rawCmdStr'])
                ),
            )

        ### Below is provided as example of how to code a function
        ### description. You should do this for each MPilot Command
        ### you create.

        self.fxnDesc['Name'] = 'PilotCmd'
        self.fxnDesc['ReqArgs'] = {
            'InFileName':'File Name',   # for commands that read
            'InField':'Field Name',     # a field name or field name list
            'Other':'TypeOfArgumentVariable'
            }
        self.fxnDesc['OptArgs'] = {
            'OutFileName':'File Name',
            'PrecursorField':'Field Name', # a field name or field name list
            }
        
        self.fxnDesc['DisplayName'] = 'UserFriendlyName'
        self.fxnDesc['ShortDesc'] = 'QuickDescription'
        self.fxnDesc['ReturnType'] = 'TypeOfReturnVariable'

    # _SetFxnDesc(self):

    def _ValidateStrCmd(self):
        # Assumes syntax of parsed command is valid
        # Right now, does not check result name

        cmdArgs = self.mptCmdStruct['parsedCmd']['arguments']

        reqArgsSet = set(self.fxnDesc['ReqArgs'])
        optArgSet = set(self.fxnDesc['OptArgs'])
        inArgsSet = set(cmdArgs.keys())

        # are all required args there?
        if not reqArgsSet.issubset(inArgsSet):
            raise Exception(
                '{}{}{}{}'.format(
                    '\n********************ERROR********************\n',
                    'Missing input argument(s):\n  {}\n'.format(
                        ' '.join(reqArgsSet - inArgsSet)
                        ),
                    'File: {}  Line number: {}\n'.format(
                        self.mptCmdStruct['cmdFileNm'],self.mptCmdStruct['lineNo']
                        ),
                    'Full command:\n{}\n'.format(self.mptCmdStruct['rawCmdStr'])
                    ),
                )

        # are all of the args either optional or required?
        if not inArgsSet.issubset(set.union(reqArgsSet,optArgSet)):
            raise Exception(
                '{}{}{}{}'.format(
                    '\n********************ERROR********************\n',
                    'Invalid input argument(s):\n  {}\n'.format(
                        ' '.join(inArgsSet - set.union(reqArgsSet,optArgSet))
                        ),
                    'File: {}  Line number: {}\n'.format(
                        self.mptCmdStruct['cmdFileNm'],self.mptCmdStruct['lineNo']
                        ),
                    'Full command:\n{}\n'.format(self.mptCmdStruct['rawCmdStr'])
                    ),
                )

        # Is the result name a valid field name
        self._ValidateArgType(self.mptCmdStruct['parsedCmd']['rsltNm'],'Field Name')
        
        # Are the all the args of the correct type?
        # loop through the allowed argument
        descArgs = cp.deepcopy(self.fxnDesc['ReqArgs'])
        descArgs.update(self.fxnDesc['OptArgs'])

        for descArgNm,descArgType in descArgs.items():

            # if the arg was not in the input args,
            # don't worry about it
            if descArgNm not in inArgsSet: 
                continue
            
            # arg type spec can be list, if not make it a
            # list for ease of coding this
            if not isinstance(descArgType,list):
                descArgType = [descArgType]

            # loop through the valid arg types
            
            argIsValid = False
            for pType in descArgType:
                if  self._IsArgType(cmdArgs[descArgNm],pType):
                    argIsValid = True
                    break
            # for pType in argType:
                
            # determine if the arg is one the valid types

            if not argIsValid:
                raise Exception(
                    '{}{}{}{}{}'.format(
                        '\n********************ERROR********************\n',
                        'Invalid argument value:\n',
                        '  Argument name: {}\n  Should be one of:\n    {}\n  Value is:{}\n'.format(
                            descArgNm,
                            '\n    '.join(descArgType),
                            cmdArgs[descArgNm]
                            ),
                        'File: {}  Line number: {}\n'.format(
                            self.mptCmdStruct['cmdFileNm'],self.mptCmdStruct['lineNo']
                            ),
                        'Full command:\n{}\n'.format(self.mptCmdStruct['rawCmdStr'])
                        ),
                    )
                
    # def _ValidateStrCmd():
    # Functions to test individual file types

    def _ValidateArgExists(self,argNm):

        if not self.ArgExists(argNm):
            raise Exception(
                '{}{}{}{}{}'.format(
                    '\n********************ERROR********************\n',
                    'Trying to acces invalid argument:\n',
                    '  Argument name: {}\n'.format(
                        argNm
                        ),
                    'File: {}  Line number: {}\n'.format(
                        self.mptCmdStruct['cmdFileNm'],
                        self.mptCmdStruct['lineNo']
                        ),
                    'Full command:\n{}\n'.format(self.mptCmdStruct['rawCmdStr'])
                ),
            )

    # def _ValidateArgExists(self,argNm):

    def _ValidateArgType(self,arg,argType):

        if not self._IsArgType(arg,argType):

            raise Exception(
                '{}{}{}{}{}'.format(
                    '\n********************ERROR********************\n',
                    'Input argument has invalid type:\n',
                    'Value: {}  Is not: {}\n'.format(
                        arg,
                        argType
                        ),
                    
                    # 'Argument name: {}  Should be: {}  Is: {}\n'.format(
                    #     'InVals',
                    #     self.fxnDesc['ReqArgs']['InVals'],
                    #     pCmd['arguments']['InVals'],
                    #     ),
                    'File: {}  Line number: {}\n'.format(
                        self.mptCmdStruct['cmdFileNm'],
                        self.mptCmdStruct['lineNo']
                        ),
                    'Full command:\n{}\n'.format(self.mptCmdStruct['rawCmdStr'])
                ),
            )

    # def _ValidateArgType(self,arg,argType):
        
    def _ValidateListLen(self,argNm,minLen,maxLen=float('inf')):

        listLen = len(self.ArgByNm('InFieldNames').replace('[','').replace(']','').split(','))
        
        if not minLen <= listLen <= maxLen:
            raise Exception(
                '{}{}{}{}{}'.format(
                    '\n********************ERROR********************\n',
                    'Argument list has invalid length:\n',
                    'Argument name: {}  Should be: {} <= Length <= {}  Length is {}.\n'.format(
                        argNm,
                        minLen,
                        maxLen,
                        len(myLst)
                        ),
                    'File: {}  Line number: {}\n'.format(
                        self.mptCmdStruct['cmdFileNm'],
                        self.mptCmdStruct['lineNo']
                        ),
                    'Full command:\n{}\n'.format(self.mptCmdStruct['rawCmdStr'])
                ),
            )
            
    # def _ValidateListLen(self,myLst,minLen,maxLen=float('inf')):

    def _ValidateEqualListLens(self,argNms):
        
        listLen = -1
        for argNm in argNms:
            pListLen = len(self.ArgByNm(argNm).replace('[','').replace(']','').split(','))
            if listLen < 0:
                listLen = pListLen
            elif pListLen != listLen:
                raise Exception(
                    '{}{}{}{}{}'.format(
                        '\n********************ERROR********************\n',
                        'List arguments do not have identical lengths:\n',
                        '  These arguments must have identical lengths:\n    {}\n'.format(
                            '\n    '.join(argNms)
                            ),
                        'File: {}  Line number: {}\n'.format(
                            self.mptCmdStruct['cmdFileNm'],
                            self.mptCmdStruct['lineNo']
                            ),
                        'Full command:\n{}\n'.format(self.mptCmdStruct['rawCmdStr'])
                    ),
                )
            
    # def _ValidateEqualListLens(self,argNms):
    
    def _ValidateArgListItemsUnique(self,argNm):

        inLst = self.ValFromArgByNm('RawValues')
        
        if not isinstance(inLst,list): return
        
        if len(set(inLst)) != len(inLst):
            raise Exception(
                '{}{}{}{}{}{}'.format(
                    '\n********************ERROR********************\n',
                    'List argument does not have unique entries:\n',
                    '   Argument name: {}\n'.format(
                        argNm,
                        ),
                    '   Argument values:\n    {}\n'.format(
                        '\n    '.join([float(x) for x in inLst]),
                        ),
                    'File: {}  Line number: {}\n'.format(
                        self.mptCmdStruct['cmdFileNm'],
                        self.mptCmdStruct['lineNo']
                        ),
                    'Full command:\n{}\n'.format(self.mptCmdStruct['rawCmdStr'])
                ),
            )

        
    # def _ValidateEqualListLens(self,argNms):

    def _ArgToList(self,argNm):
        
        rtrn = [] # empty list if the arg does not exist in the cmd
        if self.ArgByNm(argNm) is not None:
            rtrn = self.ArgByNm(argNm).replace('[','').replace(']','').split(',')

        return rtrn
    
    # def _ArgToList(self,argNm):

    def _ScalarToList(self,in_scalar):

        # If input is a scalar put it in a list, otherwise leave unchanged
        if isinstance(in_scalar,(int,float)):
            return [in_scalar]
        else:
            return in_scalar

    # def _ScalarToList(self,in_scalar):
    
    def _IsArgType(self,inStr,inType):
        
        raise Exception(
            '{}{}{}{}{}{}{}'.format(
                '\n********************ERROR********************\n',
                'Programming error:\n',
                '  Your program is using the inherited _MPilotFxnParent:_IsArgType() method.',
                '  There should be a unique _IsArgType() method for the defined MPilot command.',
                '  Check and correct the class definition of the MPilot command: {}\n'.format(self.fxnDesc['Name']),
                'File: {}  Line number: {}\n'.format(
                    self.mptCmdStruct['cmdFileNm'],self.mptCmdStruct['lineNo']
                    ),
                'Full command:\n{}\n'.format(self.mptCmdStruct['rawCmdStr'])
                ),
            )
    
    # def _IsArgType(self,inStr,inType):
                    
    def __enter__(self):
        return(self)

    def __exit__(self,exc_type,exc_value,traceback):
        if exc_type is not None:
            print exc_type, exc_value, traceback

    ## Public methods

    def InitFromParsedCmd(self,parsedCmd):
        
        self.mptCmdStruct['parsedCmd'] = cp.deepcopy(parsedCmd)
        self.InitRawCmdFromParsedCmd()
        self.InitCleanCmdFromParsedCmd()
        
    # def InitFromParsedCmd(self,parsedCmd):
    
    def ArgTypesFromArg(self,argNm):

        if argNm in self.fxnDesc['ReqArgs']:
            return self.fxnDesc['ReqArgs'][argNm]
        elif argNm in self.fxnDesc['OptArgs']:
            return self.fxnDesc['OptArgs'][argNm]
        else:
            raise Exception(
                '{}{}{}{}'.format(
                    '\n********************ERROR********************\n',
                    'Illegal argument:  {}\n'.format(argNm),
                    'File: {}  Line number: {}\n'.format(
                        self.mptCmdStruct['cmdFileNm'],self.mptCmdStruct['lineNo']
                        ),
                    'Full command:\n{}\n'.format(self.mptCmdStruct['rawCmdStr'])
                    ),
                )

    # def ArgTypesFromArg(self,argNm):
    
    def SetRsltNm(self,nm): self.mptCmdStruct['parsedCmd']['rsltNm'] = nm

    def SetArg(self,argNm,value):
        self.mptCmdStruct['parsedCmd']['arguments'][argNm] = cp.deepcopy(value)

    def SetCmdFileNm(self,nm): self.mptCmdStruct['cmdFileNm'] = nm
        
    def SetLinNo(self,num): self.mptCmdStruct['LineNo'] = nm
        
    def SetRawCmdStr(self,cmdStr): self.mptCmdStruct['rawCmdStr'] = cmdStr
        
    def SetCleanCmdStr(self,cmdStr): self.mptCmdStruct['cleanCmdStr'] = cmdStr
        
    # MPilot Function description access
    def ArgExists(self,nm):
        
        rtrn = False
        if self.mptCmdStruct is not None:
            if 'parsedCmd' in self.mptCmdStruct:
                rtrn = nm in self.mptCmdStruct['parsedCmd']['arguments']
        return rtrn
    
    def RsltNm(self): return self.mptCmdStruct['parsedCmd']['rsltNm']

    def FxnNm(self): return self.fxnDesc['Name']

    def FxnReqArgs(self): return self.fxnDesc['ReqArgs']
    
    def FxnOptArgs(self): return self.fxnDesc['OptArgs']
        
    def FxnDisplayName(self): return self.fxnDesc['DisplayName']
        
    def FxnShortDesc(self): return self.fxnDesc['ShortDesc']
        
    def FxnReturnType(self): return self.fxnDesc['ReturnType']
        
    def FxnDesc(self): return self.fxnDesc
        
    def FormattedFxnDesc(self):

        rtrnStr = ''
        if 'Name' in self.fxnDesc:
            rtrnStr = '{}{}: {}\n'.format(rtrnStr,'Function Name',self.fxnDesc['Name'])
        if 'DisplayName' in self.fxnDesc:
            rtrnStr = '{}  {}: {}\n'.format(rtrnStr,'Display Name',self.fxnDesc['DisplayName'])
        if 'ShortDesc' in self.fxnDesc:
            rtrnStr = '{}  {}: {}\n'.format(rtrnStr,'Description',self.fxnDesc['ShortDesc'])
        if 'ReturnType' in self.fxnDesc:
            rtrnStr = '{}  {}: {}\n'.format(rtrnStr,'Return Type',self.fxnDesc['ReturnType'])
        if 'InputType' in self.fxnDesc:
            rtrnStr = '{}  {}: {}\n'.format(rtrnStr,'Input Type',self.fxnDesc['InputType'])
        if 'ReqArgs' in self.fxnDesc:
            rtrnStr = '{}  Required Arguments and data types:\n'.format(rtrnStr)
            for argKey,argVal in self.fxnDesc['ReqArgs'].items():
                rtrnStr = '{}    {}: {}\n'.format(rtrnStr,argKey,argVal)
        if 'OptArgs' in self.fxnDesc:
            rtrnStr = '{}  Optional Arguments and data types:\n'.format(rtrnStr)
            for argKey,argVal in self.fxnDesc['OptArgs'].items():
                rtrnStr = '{}    {}: {}\n'.format(rtrnStr,argKey,argVal)

        return rtrnStr
        
    # def FxnDescStr(self):

    # Object item access

    def StrCmd(self): return self.mptCmdStruct
        
    def CmdFileNm(self): return self.mptCmdStruct['cmdFileNm']
        
    def LineNo(self): return self.mptCmdStruct['lineNo']
        
    def RawCmdStr(self): return self.mptCmdStruct['rawCmdStr']
        
    def CleanCmdStr(self): return self.mptCmdStruct['cleanCmdStr']
        
    def ParsedCmd(self): return self.mptCmdStruct['parsedCmd']

    def FormattedCmd(self):
        
        rtrnStr = '{} = {}(\n'.format(
            self.ParsedCmd()['rsltNm'],
            self.ParsedCmd()['cmd']
            )
        rtrnStr = '{}{}\n)'.format(
            rtrnStr,
            ',\n'.join(
                ['    {} = {}'.format(pNm,pVal) for pNm,pVal in self.ParsedCmd()['arguments'].items()]
                )
            )

        return rtrnStr

    # def FormattedCmd(self):
        
    def InitRawCmdFromParsedCmd(self):
        newCmdStr = '{} = {}(\n  '.format(
            self.mptCmdStruct['parsedCmd']['rsltNm'],
            self.fxnDesc['Name']
            )

        argLines = ['{} = {}'.format(k,v) for k,v in self.mptCmdStruct['parsedCmd']['arguments'].items()]
        newCmdStr = '{}{}'.format(newCmdStr,',\n  '.join(argLines))
        newCmdStr = '{}\n)'.format(newCmdStr)
        self.SetRawCmdStr(newCmdStr)
        
    def InitCleanCmdFromParsedCmd(self):
        newCmdStr = '{} = {}(\n  '.format(
            self.mptCmdStruct['parsedCmd']['rsltNm'],
            self.fxnDesc['Name']
            )

        argLines = ['{} = {}'.format(k,v) for k,v in self.mptCmdStruct['parsedCmd']['arguments'].items()]
        newCmdStr = '{}{}'.format(newCmdStr,',\n  '.join(argLines))
        newCmdStr = '{}\n)'.format(newCmdStr)
        self.SetCleanCmdStr(newCmdStr)
        
    def RsltNm(self):
        if self.mptCmdStruct['parsedCmd'] is not None:
            return self.mptCmdStruct['parsedCmd']['rsltNm']
        else:
            return None
        
    def Args(self):
        if self.mptCmdStruct['parsedCmd'] is not None:
            return self.mptCmdStruct['parsedCmd']['arguments']
        else:
            return None
        
    def ArgByNm(self, argNm):
        rtrn = None
        if self.mptCmdStruct['parsedCmd'] is not None:
            if argNm in self.mptCmdStruct['parsedCmd']['arguments']:
                rtrn = self.mptCmdStruct['parsedCmd']['arguments'][argNm]
        return rtrn

    def MetadataByKey(self, metaDataKey):
        rtrn = None
        if self.mptCmdStruct['parsedCmd'] is not None:
            if 'Metadata' in self.mptCmdStruct['parsedCmd']['arguments']:
                parsedMetadata = pmd.ParseMetadata(self.mptCmdStruct['parsedCmd']['arguments']['Metadata'])
                if metaDataKey in parsedMetadata:
                    rtrn = parsedMetadata[metaDataKey]

        return rtrn

    # def MetadataByKey(self, metaDataKey):

# Get some other things
    def ExecRslt(self):
        return self.execRslt

    def DependencyNms(self):
        # Parse out the dependency names from the self.mptCmdStruct
        # and return them
        return None

    # Execute the command
    def Exec(self,executedObjects):
        # This is where the the execution string gets built
        # and executed. The results are added to executedObjects.
        # executedObjects is a dictionary, each entry is a data layer.
        # Depending on how the implementde MPilot framework is
        # designed, executedObjects may contain numpy arrays (masked
        # arrays recommended), as well as data descriptors.
        #
        # Any data checking needs to be done here.

        # Throw an error if this method is not implemented in descendent
        raise Exception(
            '{}{}{}{}{}{}{}'.format(
                '\n********************ERROR********************\n',
                'Programming error:\n',
                '  Your program is using the inherited _MPilotFxnParent:Exec() method.',
                '  There should be a unique Exec() method for the defined MPilot command.',
                '  Check and correct the class definition of the MPilot command: {}\n'.format(self.fxnDesc['Name']),
                'File: {}  Line number: {}\n'.format(
                    self.mptCmdStruct['cmdFileNm'],self.mptCmdStruct['lineNo']
                    ),
                'Full command:\n{}\n'.format(self.mptCmdStruct['rawCmdStr'])
                ),
            )

        ### Below is to give you an idea of what to do. i.e. calculate
        ### something and add it to the data.
        
        self.execRslt = None # Do a calculation
        if updateProgData:
            executedObjects[self.RsltNm()] = self.execRslt
            
    # def Exec(self,executedObjects,updateProgData=True):

# class _MPilotFxnParent(object):
    
