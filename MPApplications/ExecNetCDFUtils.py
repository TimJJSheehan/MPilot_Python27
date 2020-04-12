#!/opt/local/bin/python

# Note the line above. If you are running on linux or OS X,
# you should have a similar line. This tells the operating
# system what command it should use to run this script. It
# should be pointing to version 2.7 of python

from MPCore import MPilotProgram as mpprog
from MPCore import MPilotFramework as mpf
from MPCore import MPilotParse as mpp
from collections import OrderedDict
import numpy as np
import sys
import os
import re

def UsageDie():
    print '{} ScriptFileName'.format(os.path.basename(sys.argv[0]))
    print
    print '  Execute an MPilot script file\n'
    print
    print 'OR'
    print
    print '{} -tree ScriptFileName'.format(os.path.basename(sys.argv[0]))
    print
    print '  Display the logic for an MPilot script\n' 
    print 
    print 'OR'
    print
    print '{} -macro ScriptFileName'.format(os.path.basename(sys.argv[0]))
    print
    print '  Display script with macro substitution completed\n' 
    print 
    print 'OR'
    print
    print '{} -list'.format(os.path.basename(sys.argv[0]))
    print 
    print '  List of available framework commands'
    print
    print 'OR'
    print
    print '\n{} -help [CommandName CommandName...]'.format(os.path.basename(sys.argv[0]))
    print 
    print '  Detailed description of specified CommandNames or all commands if no'
    print '  CommandNames are specified'
    print
    print 'OR'
    print
    print '\n{} -defmacros MacroNm:MacroVal[,MacroNm:MacroVal] ScriptFileName'.format(os.path.basename(sys.argv[0]))
    print 
    print '  for each MacroNm:MacroVal does a macro substitution on script'
    print '  CommandNames are specified'
    
    exit()
# def UsageDie():

def MacroSub(inStr):
    '''
    This reads in the file and does a macro substitution
    using lines like this:
    
    # MACRO SrchStr:RepStr

    to do substitutions in the script before feeding it to EEMS.
    '''
    macros = OrderedDict()
    rtrnStr = ''

    for inLine in inStr.splitlines():

        myMatch = re.match(r'^\s*#\s*MACRO\s+([^\s:]+)\s*:([^\s].*[^\s]*)\s*$',inLine)

        if myMatch is not None:
            macros[myMatch.groups()[0]] = myMatch.groups()[1]
        else:
            for srch,rep in macros.items():
                inLine = inLine.replace(srch,rep)

        rtrnStr = '{}{}\n'.format(rtrnStr,inLine)

    return rtrnStr

# def MacroSub(inFNm):

def ReadIt(inFNm):
    with open(inFNm,'r') as inF:
        rtrnStr = inF.read()
    return rtrnStr
# def ReadIt(inFNm):
        
def CreateFramework():

    return mpf.MPilotFramework([
        ('.MPNetCDFUtilsLib','MPStdLibraries')
        ])

# def CreateFramework():

def RunIt(framework,progStr):

    # This bit of code loads and runs the MPilot script
    # you specified on the command line
    with mpprog.MPilotProgram(
            framework,
            sourceProgStr = progStr
        ) as prog:

        # This runs the MPilot script you specified on the
        # command line.
        prog.Run()

# def RunIt(inFNm,outFNm)

def TreeIt(framework,inFNm):
    
    # This bit of code loads the MPilot script, generates
    # the logic tree, and prints the logic tree for the script
    # you specified on the command line
    with mpprog.MPilotProgram(
            framework,
            sourceProgFNm = inFNm
        ) as prog:

        # This runs the MPilot script you specified on the
        # command line.
        print prog.CmdTreeWithLines()

# def TreeIt(framework,inFNm):

def MacroIt(inStr,macroList):
    
    '''
    This reads in the file and does a macro substitution
    using macroList which is a list of entries like this:
    
    MacroNm:MacroVal
    '''
    rtrnStr = ''

    for inLine in inStr.splitlines():

        # don't do macros on commented lines
        if not inLine.startswith('#'):
            for macroDef in macroList:
                srch,repl = macroDef.split(':')
                inLine = inLine.replace(srch,repl)

        rtrnStr = '{}{}\n'.format(rtrnStr,inLine)
        
    # for inLine in instr.splitlines():

    return rtrnStr
    
# def MacroIt(progStr)
    
myFw = CreateFramework()

if len(sys.argv) < 2:

    UsageDie()
    
elif sys.argv[1] == '-tree' and len(sys.argv) == 3:

    TreeIt(myFw,sys.argv[2])

elif sys.argv[1] == '-macro' and len(sys.argv) == 3:

    print MacroSub(ReadIt(sys.argv[2]))
    
elif sys.argv[1] == '-list' and len(sys.argv) == 2:

    for fxnNm in myFw.FxnNames():
        print '{}'.format(fxnNm)
    print
    
elif sys.argv[1] == '-help':
    
    if len(sys.argv) == 2:
        print myFw.GetAllFormattedFxnClassInfo()
    else:
        for cmdNm in sys.argv[2:]:
            rslt = myFw.GetFxnFormattedClassInfo(cmdNm)
            if rslt is None:
                print '\n{} is not a valid command.'.format(cmdNm)
            else:
                print '\n{}'.format(rslt)

elif sys.argv[1] == '-defmacros':
    macroLst = sys.argv[2].split(',')
    progStr = MacroIt(ReadIt(sys.argv[3]),macroLst)
    RunIt(myFw,MacroSub(progStr))
    print '\nRun succeeded\n'


elif len(sys.argv) == 2: 

    RunIt(myFw,MacroSub(ReadIt(sys.argv[1])))
    print '\nRun succeeded\n'

else:

    UsageDie()

