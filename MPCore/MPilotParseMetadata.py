#!/opt/local/bin/python

from collections import OrderedDict
import re

def ParseMetadata(mdStr):

    # These parsing functions will honor quoted strings
    # with spaces and special characters.
    
    mdTokens = StripTerminalCommas(TokenizeWithQuotes(mdStr))
    # for token in mdTokens:
    #     print token
    return ParseMetadataFromTokens(mdTokens,mdStr)[0]

# def ParseMetadata(mdStr):

def ParseMetadataFromTokens(mdTokens,mdStr,tokenNdx=0):

    # The idea here is that next token will be used to
    # determine if the current token is legal
    # if not, an error is thrown

    # each token is in a tuple with an index to mdStr, so
    # that when an error is thrown, the location in the mdStr
    # can be given. 

    # Check for balanced brackets
    if [x for x in mdTokens].count('[') != [x for x in mdTokens].count(']'):
        raise Exception(
            '{}{}{}{}{}\n'.format(
                '\n********************ERROR********************\n',
                'Illegal metadata:\n\n',
                '  Unbalanced brackets: number of "[" should match number of "]"\n\n',
                '  Full metadata definition:\n',
                '{}\n\n'.format(mdStr),
                )
            )
    
    mdDict = OrderedDict()

    if not IsTokenType(mdTokens[tokenNdx],'['):
        raise Exception(
            '{}{}{}{}{}{}\n'.format(
                '\n********************ERROR********************\n',
                'Illegal metadata:\n\n',
                '  Expected: "["\n',
                '  Got: "{}"\n\n'.format(mdTokens[tokenNdx][0]),
                '  Full metadata definition:\n',
                '{}\n\n'.format(mdStr),
                )
            )
    else:
        crntTokenType = '['
        crntTokenNdx = tokenNdx
        crntToken = mdTokens[crntTokenNdx]

    # if not IsTokenType(mdTokens[tokenNdx][0],'['):

    # Go through the tokens building the dict and trapping errors
    while crntTokenNdx < len(mdTokens):

        # is this the last token
        if crntTokenNdx + 1 < len(mdTokens):
            nextToken = mdTokens[crntTokenNdx+1]
        else:
            nextToken = None

        # check for unclosed metadata string
        if crntTokenNdx is None and not crntTokenType != ']':
            raise Exception(
                '{}{}{}{}{}\n'.format(
                    '\n********************ERROR********************\n',
                    'Illegal metadata:\n',
                    '  Metadata string not closed. May be missing "]"\n'
                    '  Full metadata definition:\n'
                    '{}\n\n'.format(mdStr),
                    )
                )
        
        if crntTokenType == '[':

            if IsTokenType(nextToken,'key'):
                crntTokenNdx = crntTokenNdx + 1
                crntTokenType = 'key'
                crntToken = mdTokens[crntTokenNdx]
                continue
            else:
                raise Exception(
                    '{}{}{}{}\n'.format(
                        '\n********************ERROR********************\n',
                        'Metadata error! Metadata key must follow "["\n',
                        '  Full metadata definition:\n\n',
                        '{}\n\n'.format(mdStr)
                        )
                    )
            
        elif crntTokenType == 'key':

            if IsTokenType(nextToken,':'):
                mdKey = crntToken
                crntTokenNdx = crntTokenNdx + 1
                crntTokenType = ':'
                crntToken = mdTokens[crntTokenNdx]
                continue
            else:
                raise Exception(
                    '{}{}{}{}\n'.format(
                        '\n********************ERROR********************\n',
                        'Metadata error! ":" must follow metadata key: "{}"\n\n'.format(crntToken) ,
                        '  Full metadata definition:\n\n',
                        '{}\n\n'.format(mdStr)
                        )
                    )

        elif crntTokenType == ':':
            if IsTokenType(nextToken,'value'):
                crntTokenNdx = crntTokenNdx + 1
                crntTokenType = 'value'
                crntToken = mdTokens[crntTokenNdx]
                mdDict[mdKey] = crntToken
                mdKey = None
                continue
            elif IsTokenType(nextToken,'['):
                mdDict[mdKey],crntTokenNdx = \
                  ParseMetadataFromTokens(mdTokens,mdStr,tokenNdx=crntTokenNdx+1)
                crntTokenType = 'value'
                crntToken = mdTokens[crntTokenNdx]
                mdKey = None
                continue
            else:
                raise Exception(
                    '{}{}{}{}\n'.format(
                        '\n********************ERROR********************\n',
                        'Metadata error! Metadata Value or "[" must follow : "{}"\n\n'.format(crntToken) ,
                        '  Full metadata definition:\n\n',
                        '{}\n\n'.format(mdStr)
                        )
                    )
            
        elif crntTokenType == 'value':

            if nextToken == None:
                return mdDict, None
            elif IsTokenType(nextToken,','):
                crntTokenNdx = crntTokenNdx + 1
                crntTokenType = ','
                crntToken = mdTokens[crntTokenNdx]
                continue
            elif IsTokenType(nextToken,']'):
                crntTokenNdx = crntTokenNdx + 1
                crntTokenType = ']'
                crntToken = mdTokens[crntTokenNdx]
                continue
            else:
                raise Exception(
                    '{}{}{}{}\n'.format(
                        '\n********************ERROR********************\n',
                        'Metadata error! "," or "]" must follow meta data value "{}"\n\n'.format(crntToken) ,
                        '  Full metadata definition:\n\n',
                        '{}\n\n'.format(mdStr)
                        )
                    )
            
        elif crntTokenType == ',':
            
            if IsTokenType(nextToken,'key'):
                crntTokenNdx = crntTokenNdx + 1
                crntTokenType = 'key'
                crntToken = mdTokens[crntTokenNdx]
                continue
            else:
                raise Exception(
                    '{}{}{}{}\n'.format(
                        '\n********************ERROR********************\n',
                        'Metadata error! Metadata key must follow ","\n',
                        '  Full metadata definition:\n\n',
                        '{}\n\n'.format(mdStr)
                        )
                    )
            
        elif crntTokenType == ']':

            if nextToken is None: # End of metadata
                return mdDict, None
            elif IsTokenType(nextToken,',') or IsTokenType(nextToken,']'):
                crntTokenNdx = crntTokenNdx + 1
                continue
            else:
                raise Exception(
                    '{}{}{}{}\n'.format(
                        '\n********************ERROR********************\n',
                        'Metadata error! "]" or "," key must follow "]"\n',
                        '  Full metadata definition:\n\n',
                        '{}\n\n'.format(mdStr)
                        )
                    )

    # while len(mdTokens) > 0:

    # Because the ']' causes a return from the function,
    # there is an error if we get here

    raise Exception(
        '{}{}{}\n'.format(
            '\n********************ERROR********************\n',
            'Metadata missing closing bracket: "]"\n',
            'Metadata string is {}\n'.format(mdStr)
            )
        )

# def ParseMetadataFromTokens(mdTokens,mdStr):

def IsTokenInTypes(token,types):

    rtrn = False
    for type in types:
        if IsTokenType(token,type):
            rtrn = True
    return rtrn

# def IsTokenInTypes(token,types)
    
def IsTokenType(token,type):

    reservedChars = ['[',':',',',']']
    rtrn = False
    
    if type in reservedChars:
        rtrn = token == type
    elif type in ['key','value']:
        # if it contains non whitespace, it is valid
        if token[0] not in reservedChars and \
          re.search('[^\s]',token) is not None:
            rtrn = True

    return rtrn
    
# def IsTokenType(token,types):

def TokenizeWithQuotes(inStr):

    '''
    This routine will preserve quoted strings with any characters in them.
    Note: Multiple white space in metadata will be converted into a
    single space.
    '''
    # Split out quoted strings. Odd # strings => even # quotes

    splitStr = inStr.split('"')

    if len(splitStr) %2 != 1:
        raise Exception('Unmatched quotes in command:\n{}\n'.format(mdStr))

    tokens = []

    for ndx in range(len(splitStr)):
        if ndx%2 == 0:
            tokens += TokenizeString(splitStr[ndx])
        else:
            tokens.append(re.sub(r'\s+',' ','"{}"'.format(splitStr[ndx])))

    return tokens

# def TokenizeMetadata(mdStr):

def TokenizeString(inStr):

    '''
    Because for the life of me, I cannot figure out how to
    include the newline properly in re.match, I'm preconverting
    newlines into spaces.

    Whitespace is used as a separator for tokens when it is not
    padding a special character. e.g. "This = x y" will be converted
    to ['This','=','x','y']
    '''
    inStr = re.sub(r'[\n\r]+',' ',inStr)
    
    singleCharTokens = ['\[','\]','\(','\)',',','=',':']
    singleCharTokensStr = ''.join(singleCharTokens)

    # Pull off leading characters
    
    tokens = []

    while len(inStr) > 0:
        
        for tok in singleCharTokens:
            # can we pull off a single character token
            tokMatch = re.match(r'\s*({})\s*(.*)$'.format(tok),inStr)
            if tokMatch is not None: break
            
        # can we pull off a string token
        if tokMatch is None:
            tokMatch = re.match(r'\s*([^\s{}]+)\s*(.*)$'.format(singleCharTokensStr),inStr)

        if tokMatch is not None:
            tokens.append(tokMatch.groups()[0])
            inStr = tokMatch.groups()[1]
            
        else:
            # Error check or return
            if re.match(r'^\s*$',inStr) is not None:
                raise Exception(
                    '{}{}{}{}'.format(
                        '\n********************ERROR********************\n',
                        'Programming error: cannot parse string:\n',
                        '{}\n'.format(inStr)
                        )
                    )

    # while len(inStr) > 0:

    return tokens

# def TokenizeString(inStr):

def StripTerminalCommas(inTokens):

    for ndx in list(reversed(range(len(inTokens)-1))):
        if inTokens[ndx] == ',' and inTokens[ndx+1] in [']',')']:
            inTokens.pop(ndx)

    return inTokens

# def StripTerminalCommas(inTokens)
