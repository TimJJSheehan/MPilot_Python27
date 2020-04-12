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
    print '''
{} ScriptFileName

  Execute an MPilot script file

OR

{} -tree ScriptFileName

  Display the logic for an MPilot script

OR

{} -list

  List of available framework commands'

OR

{} -help [CommandName CommandName...]

  Detailed description of specified CommandNames or all commands if no
  CommandNames are specified
'''.format(
    os.path.basename(sys.argv[0]),
    os.path.basename(sys.argv[0]),
    os.path.basename(sys.argv[0]),
    os.path.basename(sys.argv[0])
    )
    exit()
# def UsageDie():

def MacroSub(inFNm):
    '''
    This reads in the file and does a macro substitution
    using lines like this:
    
    # MACRO SrchStr:RepStr

    to do substitutions in the script before feeding it to EEMS.
    '''
    with open(inFNm,'r') as inF:

        macros = OrderedDict()
        rtrnStr = ''

        for inLine in inF.readlines():

            myMatch = re.match(r'^\s*#\s*MACRO\s+([^\s:]+)\s*:([^\s].*[^\s])\s*$',inLine)

            if myMatch is not None:
                macros[myMatch.groups()[0]] = myMatch.groups()[1]
            else:
                for srch,rep in macros.items():
                    inLine = inLine.replace(srch,rep)

            rtrnStr = '{}{}'.format(rtrnStr,inLine)

        return rtrnStr

    # with open(sys.argv[1],'r') as inF, \...
    
# def MacroSub(inFNm):

def CreateFramework():

    return mpf.MPilotFramework([
        ('.MPEEMSBasicLib','MPStdLibraries'),
        ('.MPEEMSFuzzyLogicLib','MPStdLibraries'),
        ('.MPEEMSGraphLib','MPStdLibraries'),
        ('.MPEEMSCSVIO','MPStdLibraries'),
        ('.MPEEMSStatsLib','MPStdLibraries'),
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

# def CreateJSONFile(inFNm,outFNm)

def TreeIt(framework,inFNm):
    
    # This bit of code loads the MPilot script, generates
    # the logic tree, and prints the logic tree for the script
    # you specified on the command line
    with mpprog.MPilotProgram(
            framework,
            inFNm
        ) as prog:

        # This runs the MPilot script you specified on the
        # command line.
        print prog.CmdTreeWithLines()

# def TreeIt(framework,inFNm):

myFw = CreateFramework()

if len(sys.argv) < 2:

    UsageDie()
    
elif sys.argv[1] == '-tree' and len(sys.argv) == 3:

    TreeIt(myFw,sys.argv[2])
    
elif sys.argv[1] == '-list' and len(sys.argv) == 2:

    print
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

elif len(sys.argv) == 2: 

    RunIt(myFw,MacroSub(sys.argv[1]))
    print '\nRun succeeded\n'

else:

    UsageDie()

