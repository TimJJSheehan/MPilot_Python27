from collections import OrderedDict
import MPilotParseMetadata as pmd
import re
import copy as cp

'''
MPilot parser

Functions to:

 Parse a file into commands
 Parse a string into commands
 Parse a command into arguments

tjs 2016.06.20
  Started
tjs 2016.11.04
  Added parsing of a string
tjs 2017.05.11
  Added limited support for metadata. Metadata
  is of the form:
  
  Metadata = [key1:value1,...,keyn:valuen]

  Although much of what has been put into place in MPilotParseMetadata.py
  supports heirarchical metadata with quoted values, of the form:

  Metadata = [key1:value1,key2:[key21,value21,key22:"value22 with spaces, and chars#"],...,keyn:valuen]

  Time did not permit the modifications to functions in this file to implement that.

  In the next round of revisions, I will change the parsing to handle the nested/quoted
  metadata. This will likely involve a more tokenized approach to the command as a whole.
  Work towards this end in the now unused file NewParse.py

'''

def ParseStringToCommands(
    inStr,
    objNm=None
    ):
    # Each command must start on a new line
    # a mptCmdStruct is a list of dicts for the command string

    mptCmdStructs = OrderedDict()
    mptCmdStructstartLineNo = 0
    if objNm is None: objNm = 'Input string'

    cmdLine = ''      # buffer to build command from lines of input file
    inParens = False  # whether or not parsing is within parentheses
    parenCnt = 0      # count of parenthesis levels
    inLineCnt = 0     # line number of input file for error messages.

    rawCmd = ''
    cleanCmd = ''

    # Basically a finite state machine to build individual commands
    
    for inLine in inStr.split('\n'):

        inLineCnt +=1

        cleanLine = re.sub('#.*$','',inLine)
        cleanLine = re.sub('\s+','',cleanLine)
        
        # Only start gathering a command where
        # a command starts
        if rawCmd == '':
            if cleanLine == '':
                continue
            else:
                mptCmdStructstartLineNum = inLineCnt

        rawCmd = '{}{}\n'.format(rawCmd,inLine)

        for charNdx in range(len(cleanLine)):
            cleanCmd += cleanLine[charNdx]
            if cleanLine[charNdx] == '(':
                inParens = True
                parenCnt += 1
            elif cleanLine[charNdx] == ')':
                parenCnt -= 1

            if parenCnt < 0:
                raise Exception(
                    '{}{}{}{}'.format(
                        '\n********************ERROR********************\n',
                        'Unmatched right paren *)*\n',
                        '  input: {}, line {}:\n'.format(objNm,inLineCnt),
                        '  {}\n'.format(inLine),
                        )
                    )
            if inParens and parenCnt == 0:
                if charNdx < (len(cleanLine)-1):
                    raise Exception(
                        '{}{}{}{}'.format(
                            '\n********************ERROR********************\n',
                            'Extraneous characters beyond end of command\n',
                            '  Input: {}, line {}:\n'.format(objNm,inLineCnt),
                            '  {}\n'.format(inLine),
                            )
                        )

                else:

                    cleanCmd = re.sub(',+([\]\)])',r'\1',cleanCmd)
                    
                    mptCmdStructTmp = ({
                        'cmdFileNm':objNm,
                        'lineNo':mptCmdStructstartLineNum,
                        'rawCmdStr':rawCmd,
                        'cleanCmdStr':cleanCmd,
                        })

                    mptCmdStructTmp['parsedCmd'] = ParseCommandToArgs(mptCmdStructTmp)
                    
                    rsltNm = mptCmdStructTmp['parsedCmd']['rsltNm']

                    if mptCmdStructTmp['parsedCmd']['rsltNm'] in mptCmdStructs:
                        raise Exception(
                            '{}{}{}{}{}{}'.format(
                                '\n********************ERROR********************\n',
                                'Result appears more than once in input: {}.\n'.format(objNm),
                                'First occurence: line {}:\n'.format(mptCmdStructs[rsltNm]['lineNo']),
                                '  Command:\n{}\n'.format(mptCmdStructs[rsltNm]['rawCmdStr']),
                                'Second occurence: line {}:\n'.format(mptCmdStructTmp['lineNo']),
                                '  Command:\n{}\n'.format(mptCmdStructTmp['rawCmdStr']),
                                )
                            )

                    mptCmdStructs[mptCmdStructTmp['parsedCmd']['rsltNm']] = cp.deepcopy(mptCmdStructTmp)

                    rawCmd = ''
                    cleanCmd = ''
                    inParens = False
                    parenCnt = 0

        # for charNdx in range(len(cleanLine)):
    # for inLine in inStr.split('\n'):

    #    if cleanCmd != '':
    if rawCmd != '':
        raise Exception(
            '{}{}{}{}'.format(
                '\n********************ERROR********************\n',
                'Incomplete command in file:\n',
                '  input: {}, line {}:\n'.format(objNm,mptCmdStructstartLineNum),
                'Command:\n{}\n'.format(rawCmd),
                )
            )

    return mptCmdStructs

# def ParseStringToCommands(inFDef):

def ParseFileToCommands(inFDef):
    
    if isinstance(inFDef, basestring):
        fObj = open(inFDef, 'rU')
    else:
        fObj = inFDef

    with fObj as inFile:
        rtrn = ParseStringToCommands(fObj.read(),fObj.name)

    return rtrn

# def ParseFileToCommands(inFDef):

def ParseCommandToArgs(mptCmdStruct):
    # mptCmdStruct is a dict:
    # 'lineNo': the command's line within the input file
    # 'rawCmdStr': the command string as it appeared in the input
    #    file, comments, line breaks, and all
    # 'cleanCmdStr': the command string stripped of all its
    #    comments, line breaks, and all

    # strip white space
    cmdStr = re.sub('\s+','',mptCmdStruct['cleanCmdStr'])

    # parse the command string into result, command name, and arguments

    cmdMatchRE = re.compile(r'([^=]+)=([^\(]+)\(([^\)]+)\)')
    exprParse = cmdMatchRE.match(cmdStr)
    
    if not exprParse or len(exprParse.groups()) != 3:
        raise Exception(
            '{}{}{}{}'.format(
                '\n********************ERROR********************\n',
                'Invalid command format.\n',
                'File: {}  Line number: {}\n'.format(mptCmdStruct['cmdFileNm'],mptCmdStruct['lineNo']),
                'Full command:\n{}\n'.format(mptCmdStruct['rawCmdStr'])),
                )


    parsedCmd = OrderedDict()
    
    # Every command must have a result
    if exprParse.groups()[0] is None:
        raise Exception(
            '{}{}{}{}'.format(
                '\n********************ERROR********************\n',
                'Command is missing result name.\n',
                'File: {}  Line number: {}\n'.format(mptCmdStruct['cmdFileNm'],mptCmdStruct['lineNo']),
                'Full command:\n{}\n'.format(mptCmdStruct['rawCmdStr'])),
            )


    parsedCmd['rsltNm'] = exprParse.groups()[0]

    if (not re.match(r'^[a-zA-Z0-9\-\_]+$',parsedCmd['rsltNm']) and
        not re.match(r'^\[[a-zA-Z0-9\-\_,]+\]$',parsedCmd['rsltNm'])
        ):
        raise Exception(
            '{}{}{}{}'.format(
                '\n********************ERROR********************\n',
                'Invalid result name in command.\n',
                'File: {}  Line number: {}\n'.format(mptCmdStruct['cmdFileNm'],mptCmdStruct['lineNo']),
                'Full command:\n{}\n'.format(mptCmdStruct['rawCmdStr'])),
            )


    parsedCmd['cmd'] = exprParse.groups()[1]
    
    # Parse out the arguments
    argStr = exprParse.groups()[2]
    argPairs = []

    while argStr != '':
        argPairMatchObj = re.match(r'\s*([^=]*=\s*\[[^\[]*\])\s*,*\s*(.*)',argStr)
        if argPairMatchObj:
            argPairs.append(argPairMatchObj.groups()[0])
            argStr = argPairMatchObj.groups()[1]
        else:
            argPairMatchObj = re.match(r'\s*([^=,]*=\s*[^,]*)\s*,*\s*(.*)',argStr)
            if argPairMatchObj:
                argPairs.append(argPairMatchObj.groups()[0])
                argStr = argPairMatchObj.groups()[1]
            else:
                raise Exception(
                    '{}{}{}{}{}'.format(
                        '\n********************ERROR********************\n',
                        'Invalid argument specification.\n',
                        'File: {}  Line number: {}\n'.format(mptCmdStruct['cmdFileNm'],mptCmdStruct['lineNo']),
                        'Section of command: {}\n'.format(argStr),
                        'Full command:\n{}\n'.format(mptCmdStruct['rawCmdStr']),
                        )
                    )

        # if argPair:...else:

    # while argStr != '':

    parsedCmd['arguments'] = OrderedDict()
    for argPair in argPairs:

        argTokens = re.split(r'\s*=\s*',argPair,2)

        argTokens[0] = argTokens[0].strip()
        argTokens[1] = argTokens[1].strip()

        argTokens[1] = re.sub(r'\s*\[\s*','[',argTokens[1])
        argTokens[1] = re.sub(r'\s*\]\s*',']',argTokens[1])
        argTokens[1] = re.sub(r'\s*,\s*',',',argTokens[1])

        if (len(argTokens) != 2
            or argTokens[0] == ''
            or argTokens[1] == ''
            or argTokens[0] in parsedCmd
            ):

            raise Exception(
                '{}{}{}{}'.format(
                    '\n********************ERROR********************\n',
                    'Invalid argument specification. Line number {}\n'.format(mptCmdStruct['lineNo']),
                    'Argument specification: {}\n'.format(argPair),
                    'Full command:\n{}\n'.format(mptCmdStruct['rawCmdStr']),
                    )
                )

        parsedCmd['arguments'][argTokens[0]] = argTokens[1]
            
    # for argPair in argPairs:

    return parsedCmd

# def ParseCommandToArgs(inCmd):
    
def StringSpaceClean(inStr):
    '''
    When we split on quotes, we will get an odd number of
     strings if we have an even number of quotes, with odd
     indexed strings being quoted strings. We can strip
     spaces out of the strings that were not quoted (i.e,
     even indexed strings). Then we just join the split strings
     back together with quotes to get the stripping we were after.
     '''

    splitStr = inStr.split('"')
    if len(splitStr) %2 != 1:
        raise Exception('Unmatched quotes')
    for ndx in xrange(0, len(splitStr), 2):
        splitStr[ndx] = re.sub('\s+','',splitStr[ndx])

    # newlines to spaces
    return('"'.join(splitStr).replace('\n',' '))

# def StringSpaceCleaner(inStr):
