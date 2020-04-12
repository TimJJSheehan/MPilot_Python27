import MPilotParse as mpparse
from collections import OrderedDict
import os.path
import datetime

class MPilotProgram(object):
    
    def __init__(
        self,
        mpFramework,          # An MPilot Framework
        sourceProgFNm = None, # File name for MPilot script
        sourceProgStr = None, # MPilot script as a single string
        ):

        self.unorderedMPCmds = OrderedDict()
        self.orderedMPCmds = None
        self.mpFramework = mpFramework
        self.rslts = {}

        if sourceProgFNm is not None:
            
            if sourceProgStr is not None:
                raise Exception(
                    '{}{}{}{}'.format(
                        '\n********************ERROR********************\n',
                        'Illegal to specify both a sourceProgFNm and a sourceProgStr\n',
                        'sourceProgFNm: {}'.format(sourceProgFNm),
                        'sourceProgStr: {}...'.format(sourceProgStr[0:100]),
                    )
                )
                
            self.LoadMptFile(sourceProgFNm)

        elif sourceProgStr is not None:
            
            self.LoadMptStr(sourceProgStr)
            
        # if sourceProgFNm is not None:
        
    # def __init__(...)

    def __enter__(self):
        return(self)

    def __exit__(self,exc_type,exc_value,traceback):
        if exc_type is not None:
            print exc_type, exc_value, traceback

    def _OrderCmds(self):

        # If an item has a dependency that is not the result
        # of any command, then there is an unfulfillable
        # dependency

        for rsltNm,mpCmd in self.unorderedMPCmds.items():
        
            if mpCmd.DependencyNms() is not None:
                if not set(mpCmd.DependencyNms()).issubset(set(self.unorderedMPCmds)):
                    raise Exception(
                        '{}{}{}{}'.format(
                            '\n********************ERROR********************\n',
                            'Required variable(s) for computing a result is(are) missing:\n  {}\n'.format(
                                ' '.join(set(mpCmd.DependencyNms()) - set(self.unorderedMPCmds))
                            ),
                            'File: {}  Line number: {}\n'.format(
                                mpCmd.CmdFileNm(),
                                mpCmd.LineNo()
                            ),
                            'Full command:\n{}\n'.format(mpCmd.RawCmdStr())
                        )
                    )
                    
        # Now to order commands...    
        # Idea is to iteratively go through the unordered commands
        # list, moving items whose dependencies are in the ordered
        # commands to the end of the ordered commands.
        # The first items to go into the ordered commands will be
        # those with no dependencies.
        # If an interation moves no items to the ordered commands,
        # then there is a circular reference somewhere in the
        # unordered commands.

        self.orderedMPCmds = OrderedDict()
        unorderedCmdNms = self.unorderedMPCmds.keys()
        
        while len(unorderedCmdNms) > 0:
            startUnorderedLen = len(unorderedCmdNms)
            for rsltNm in unorderedCmdNms:
                if self.unorderedMPCmds[rsltNm].DependencyNms() is None:
                    # Command has no dependencies and thus is good to run
                    self.orderedMPCmds[rsltNm] = self.unorderedMPCmds[rsltNm]
                    unorderedCmdNms.remove(rsltNm)
                elif set(self.unorderedMPCmds[rsltNm].DependencyNms()).issubset(set(self.orderedMPCmds)):
                    # All dependencies are already in ordered commands
                    self.orderedMPCmds[rsltNm] = self.unorderedMPCmds[rsltNm]
                    unorderedCmdNms.remove(rsltNm)
        
            if startUnorderedLen == len(unorderedCmdNms):
                raise Exception(
                    '{}{}{}{}{}'.format(
                        '\n********************ERROR********************\n',
                        'Circular reference in script. One or more of the\n',
                        '  result variables in the script depends on itself.\n',
                        'File: {}\n'.format(mpCmd.CmdFileNm()),
                        'Cicular reference among commands for these variables:\n  {}'.format(
                            ' '.join(self.unorderedMPCmds)
                        )
                    )
                )
        
    # def _OrderCommands()

    def _ParseDict(self,rsltNm,treeImage,lvl):

        treeImage.append((rsltNm,lvl))
        lvl += 1

        if self.orderedMPCmds[rsltNm].DependencyNms() is not None:
            for subRsltNm in self.orderedMPCmds[rsltNm].DependencyNms():
                self._ParseDict(subRsltNm,treeImage,lvl)

        return treeImage

    # def _ParseDict(self,rsltNm,treeImage,lvl):

    def _ClearCmds(self):

        self.unorderedMPCmds = OrderedDict()
        self.orderedMPCmds = None
        
    # def _ClearCmds(self):
    
    # public classes

    def LoadMptFile(self,sourceProgFNm):
        # Replaces current set of commands with new commands
        # from a file.

        self._ClearCmds()
        
        if sourceProgFNm is not None:
            if not os.path.isfile(sourceProgFNm):
                raise Exception(
                    '{}{}'.format(
                        '\n********************ERROR********************\n',
                        'File does not exist: {}\n'.format(sourceProgFNm)
                        )
                    )

            self.mptCmdStructs = mpparse.ParseFileToCommands(sourceProgFNm)
                
            for rsltNm,mptCmdStruct in self.mptCmdStructs.items():
                self.unorderedMPCmds[rsltNm] = \
                  self.mpFramework.CreateFxnObject(mptCmdStruct['parsedCmd']['cmd'],mptCmdStruct)

            self._OrderCmds()

    # def LoadMptFile(self,sourceProgFNm):

    def LoadMptStr(self,inStr): # The string is a complete program
        # Replaces current set of commands with new commands
        # from a string.

        self._ClearCmds()
        
        # Parse the string into commands
        self.mptCmdStructs = mpparse.ParseStringToCommands(inStr,'Script string')

        # Create the commands and save them in this program object
        for rsltNm,mptCmdStruct in self.mptCmdStructs.items():
            self.unorderedMPCmds[rsltNm] = \
              self.mpFramework.CreateFxnObject(mptCmdStruct['parsedCmd']['cmd'],mptCmdStruct)

    # def LoadMptStr(self,sourceProgFNm):
        
    def MPFramework(self):
        # Access to the underlying mpFramework. Not necessarily
        # recommended.
        return self.mpFramework

    def AllMPFunctionInfo(self):
        return self.mpFramework.GetAllFxnClassInfo()

    def AllFormattedMPFunctionInfo(self):
        return self.mpFramework.GetAllFormattedFxnClassInfo()
        
    def AddCmd(self,mpCmd):
        
        if self.CmdExists(mpCmd):
            raise Exception(
                '{}{}{}{}'.format(
                    '\n********************ERROR********************\n',
                    'Attempt to duplicate result name:\n  {}\n'.format(mpCmd.RsltNm()),
                    'File: {}  Line number: {}\n'.format(
                        mpCmd.CmdFileNm(),
                        mpCmd.LineNo()
                    ),
                    'Full command:\n{}\n'.format(mpCmd.RawCmdStr())
                )
            )
                    
        self.unorderedMPCmds[mpCmd.RsltNm()] = mpCmd
        
        self.orderedMPCmds = None # invalidate ordering
        
    # def AddCmd(self,mpCmd):

    def CmdExists(self,mpCmd):
        return mpCmd.RsltNm() in self.unorderedMPCmds

    def RsltNmExists(self,rsltNm):
        return rsltNm in self.unorderedMPCmds
            
    def DelCmd(self,mpCmd):
        
        self.DelCmdByRsltNm(mpCmd.RsltNm())
        
    # def DelCmd(self,mpCmd):
        
    def DelCmdByRsltNm(self,rsltNm):
        
        if rsltNm in self.unorderedMPCmds:
            del self.unorderedMPCmds[rsltNm]

        self.orderedMPCmds = None # invalidate ordering
        
    # def DelCmd(self,mpCmd):
        
    def OrderedCmds(self):
        if self.orderedMPCmds is None:
            self._OrderCmds()
        return self.orderedMPCmds

    def UnorderedCmds(self):
        return self.unorderedMPCmds

    def ProgAsText(self):

        return '{}\n'.format(
            '\n\n'.join([c.FormattedCmd() for c in self.unorderedMPCmds.values()])
            )
        
    # def ProgAsText(self):

    def CmdTree(self):

        if self.orderedMPCmds is None:
            self._OrderCmds()
            
        # find the top node(s)
        
        # Create a dictionary of all the rsltNms that are
        # depended upon
        depends = {}
        for rsltNm,mpCmd in self.orderedMPCmds.items():
            if mpCmd.DependencyNms() is not None:
                for depRsltNm in mpCmd.DependencyNms():
                    depends[depRsltNm] = True

        # Compare depends{} contents with all commands to find those
        # upon which no others depend
        
        topNodeFlds = []
        for rsltNm in self.orderedMPCmds.keys():
            if rsltNm not in depends.keys():
                topNodeFlds.append(rsltNm)

        topNodeFlds = []
        for rsltNm in self.orderedMPCmds.keys():
            if rsltNm not in depends.keys():
                topNodeFlds.append(rsltNm)

        # order nodes by dependency, tracking depth

        treeImage = []
        for topNodeFld in topNodeFlds:
            self._ParseDict(topNodeFld,treeImage,0)

        return treeImage
    
    # def CmdTree(self):

    def CmdTreeWithLines(self):
        
        rtrn = ''
        for rsltNm,lvl in self.CmdTree():
            if lvl > 0:
                rtrn += '{}'.format((lvl-1) * '|   ')
                rtrn += '|---'
                rtrn += '{}\n'.format(rsltNm)
            else:
                rtrn += '{}\n'.format(rsltNm)

        return rtrn
                
    # def CmdTreeWithLines(self):

    def Run(self):
        
        if self.orderedMPCmds is None:
            self._OrderCmds()
            
        for mpCmdNm,mpCmd in self.orderedMPCmds.items():
            mpCmd.Exec(self.rslts)

    def Rslts(self):
        rtrn = OrderedDict()
        for rsltKey,rsltObj in self.rslts.items():
            rtrn[rsltKey] = rsltObj.ExecRslt()
        return rtrn

    def RsltByNm(self,nm):
        return self.rslts[nm].ExecRslt()

    def CmdByNm(self,nm):
        return self.unorderedMPCmds[nm]

    def ClearRslts(self): self.rslts = {}

# class MPilotProgram(object):
