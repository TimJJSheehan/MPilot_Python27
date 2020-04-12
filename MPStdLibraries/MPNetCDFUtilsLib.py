'''
This is a big bunch of utilities designed to be used with NetCDF files.

The big idea here is that nc variables or global metadata can be read
in as a variable, worked on to reduce, filter, trim, combine, compare,
mask, and manipulate in other ways. Once the work is done, the
resulting variables can be graphed or written to file.

The work of manipulating data and maintaining indices is done by the
classes defined in this file and by the classes defined in
MPNetCDF4Variable.py.

If this is successful, users will have an MPilot workflow-based set
of utilities that will allow them to handle common NetCDF variable
manipulations without having to do a lot of the tedious work involved
with working with NetCDF variables and files using the python netCDF4
module directly.

The audience for this thing is scientists who understand their data
and need to do simple manipulations simply.

'''

from collections import OrderedDict
import netCDF4 as nc
import numpy as np
from scipy import stats
import copy as cp
from MPCore import MPilotFxnParent as mpfp
from MPUtilities import MPNetCDF4Variable as mpncv
import re
import os
import textwrap as tw
from MPUtilities import ArrayCongruence as ac
import matplotlib.pyplot as plt
from MPUtilities import Smooth as smooth
import random as rndm
from MPUtilities import NumpyToPng as nptopng
import csv

def _redimension_arr_for_reduction(in_arr,axes):

    '''
    Redimension an array so that the last dimension contains
    the axes over which a reduction is to take place. This is
    necessary for some specialized operations that can't easily
    be done over multiple axes. For instance, mode.
    
    For example say you have a 3d array with a shape (100,200,300),
    and you want to take the mode over the first two axes. The
    call you make to this function would be:

    _redimension_arr_for_reduction(my_arr,[0,1])

    This function would transpose and reshape the array and return an
    array with this shape:
    
    (300,20000)

    With the first axis corresponding to the third axis of the input
    array.

    You could then take the mode

    stats.mode(reshaped_arr,axis = len(reshape_dims)-1,nan_policy = 'propagate')

    Again, the remain dimension corresponds to the third dimension of your
    input array.

    '''

    # make axes iterable if they are not. e.g. if axes is a single value
    try:
        _ = (dim for dim in axes)
    except TypeError:
        axes = [axes]

    # See initial comments for general strategy
    
    transpose_axes = []   # these are preseved in both size and order
    reshape_dims = []     # these are combined and made the last axis

    # which axes to preserve
    for axis in range(len(in_arr.shape)):
        if axis not in axes:
            transpose_axes.append(axis)
            reshape_dims.append(in_arr.shape[axis])

    reshape_last_dim = 1
    # which axes to combine
    for axis in axes:
        transpose_axes.append(axis)
        reshape_last_dim *= (in_arr.shape[axis])

    reshape_dims.append(reshape_last_dim)

    # transpose, reshape, take mode, and reduce to preserved axes

    return in_arr.transpose(transpose_axes).reshape(reshape_dims)

# def _redimension_arr_for_reduction(in_arr,axes):

def _mode_by_axes(in_arr,axes):

    reshaped_arr = _redimension_arr_for_reduction(in_arr,axes)

    # mode unmasks the array. the masked_where statements remask it appropriately
    mode_arr = stats.mode(reshaped_arr,axis = reshaped_arr.ndim-1,nan_policy = 'propagate')[0]
    mode_arr = np.ma.masked_where(mode_arr == reshaped_arr.fill_value, mode_arr, copy = False)
    mode_arr = np.ma.masked_where(mode_arr == float('nan'), mode_arr, copy = False)

    # get rid of that pesky extra dimension
    sq_arr = np.ma.squeeze(mode_arr)

    return sq_arr

# def _mode_by_axes(in_arr,axes):
 
class _NetCDFUtilParent(mpfp._MPilotFxnParent):
    
    def __init__(
        self,
        mpt_cmd_struct=None,
        ):

        super(_NetCDFUtilParent,self).__init__(mpt_cmd_struct)
        
    # def __init__(...)

    def _IsArgType(self,inStr,inType):

        # If the string has a list of single word values it
        # can be.
        if inType.startswith('One of|'):
            if inStr in inType.split('|',2)[1].split():
                return True
            else:
                return False

        if inType.startswith('Tuple:'):
            
            types = inType.split(':',1)[1].split(':')
            types = [s.strip() for s in types]

            # strip out the whitespace and split the tup
            tup_entries = ''.join(inStr.split()).split(':')

            # check for length of tup
            if len(tup_entries) != len(types):
                return False

            # check type of each tup entry
            for tup_entry, type in zip(tup_entries,types):
                if not self._IsArgType(tup_entry,type):
                    return False
                
            # if we got to here, it's good!
            return True

        # if inType.startswith('Tuple:'):

        if inType.startswith('Tuple List:'):
            # a tuple list would look like this:
            # [lat:44.75:49.00, lon:-127.25:-129.25]
            types = inType.split(':',1)[1].split(':')
            types = [s.strip() for s in types]
            tups = inStr.replace('[','').replace(']','').split(',')
            
            for tup in tups:
                
                # strip out the whitespace and split the tup
                tup_entries = ''.join(tup.split()).split(':')
                
                # check for length of tup
                if len(tup_entries) != len(types):
                    return False
                
                # check type of each tup entry
                for tup_entry, type in zip(tup_entries,types):
                    if not self._IsArgType(tup_entry,type):
                        return False

            # if we got to here, it's good!
            return True

        # if inType.startswith('Tuple List:'):
            
        # Otherwise do the check by type.
        if re.match(r'.+ List$',inType):
            theType = inType.replace(' List','',1)
            if not re.match(r'\[.*\]',inStr): return False
            theStr = inStr.replace('[','').replace(']','')
        else:
            theStr = inStr
            theType = inType
            
        theStrs = theStr.split(',')
        if len(theStrs) == 0: return False

        for theStr in theStrs:
            if theType == 'Any':
                pass
            elif theType == 'String':
                pass
            elif theType == 'File Name':
                if not re.match(r'([a-zA-Z]:[\\/]){0,1}[\w\\/\.\- ~]*\w+\s*$',theStr):
                    return False
            elif theType == 'Field Name':
                if not re.match(r'^[a-zA-Z][-\w]*$',theStr):
                    return False
            elif theType == 'Integer':
                if not re.match(r'^[+-]{0,1}[0-9]+$',theStr):
                    return False
            elif theType == 'Binary':
                if not re.match(r'^[01]+$',theStr):
                    return False
            elif theType == 'Positive Integer':
                if not re.match(r'^[+]{0,1}^[0-9]$',theStr):
                    return False
                else:
                    if int(theStr) < 1:
                        return False
            elif theType == 'Float':
                if not re.match(r'^[+-]{0,1}([0-9]+\.*[0-9]*)(e[+-]{0,1}[0-9]+){0,1}$|(^[+-]{0,1}\.[0-9]+)(e[+-]{0,1}[0-9]+){0,1}$',theStr):
                    return False
                
            elif theType == 'Fuzzy Value':

                if not re.match(r'^[+-]{0,1}([0-9]+\.*[0-9]*)(e[+-]{0,1}[0-9]+){0,1}$|(^[+-]{0,1}\.[0-9]+)(e[+-]{0,1}[0-9]+){0,1}$',theStr):
                    return False
                if not self.fuzzyMin <= float(theStr) <= self.fuzzyMax:
                    return False

            elif theType == 'Positive Float':
                if not re.match(r'^[+]{0,1}([0-9]+\.*[0-9]*)(e[+-]{0,1}[0-9]+){0,1}$|(^[+-]{0,1}\.[0-9]+)(e[+-]{0,1}[0-9]+){0,1}$',theStr):
                    return False
            elif theType == 'Boolean':
                if not (inStr == 'True' or inStr == 'False'):
                    return False                
            else:
                raise Exception(
                    '{}{}{}'.format(
                        '\n********************ERROR********************\n',
                        'Class definition error.\n',
                        'Illegal argument type in function descriptions: {}'.format(inType)
                        )
                    )
        # for theStr in theStrs:

        return True
    
    # def _IsArgType(self,inStr,type):

    def ValFromArgByNm(self,argNm):

        # Return an actual value (be that a single value or a list)
        # of an argument. Arguments are read as strings, so this
        # converts that string into what it is supposed be. Right
        # now, it only does conversion on variables that are numeric.
        # If there is an option for what the argument can be, it
        # finds the most restrictive (e.g. integer instead of float)
        # definition for the argument and does the conversion based
        # on that.

        rtrn = None
        arg = self.ArgByNm(argNm)

        if arg is None: return None

        # Get the legal arg types for this argument
        if argNm in self.fxnDesc['ReqArgs']:
            legalArgTypes = self.fxnDesc['ReqArgs'][argNm]
        elif argNm in self.fxnDesc['OptArgs']:
            legalArgTypes = self.fxnDesc['OptArgs'][argNm]
        else:
            raise Exception(
                '{}{}{}{}{}{}'.format(
                    '\n********************ERROR********************\n',
                    'Programming error:\n',
                    '  Cannot find description of argter: {}'.format(argNm),
                    '  Check and correct the class definition of the MPilot command: {}\n'.format(
                        self.fxnDesc['Name']
                        ),
                    'File: {}  Line number: {}\n'.format(
                        self.mptCmdStruct['cmdFileNm'],self.mptCmdStruct['lineNo']
                        ),
                    'Full command:\n{}\n'.format(self.mptCmdStruct['rawCmdStr'])
                    ),
                )

        # if argNm self.fxnDesc['ReqArgs']:

        # From the choices it is allowed to be, find out
        # the most restrictive type that fits. Note,
        # order of checking matters here.

        argType = None
        
        if not isinstance(legalArgTypes, list):
            
            argType = legalArgTypes
            
        else:

            for simplestArgType in [
                'Positive Integer',
                'Integer',
                'Positive Integer List',
                'Integer List',
                'Float',
                'Fuzzy Value',
                'Positive Float',
                'Float List',
                'Fuzzy Value List',
                'Positive Float List',
                'Field Name',
                'Field Name List'
                ]:

                # Check the options for this argument
                if simplestArgType in legalArgTypes:
                    # is it this?
                    if self._IsArgType(arg,simplestArgType):
                        argType = simplestArgType
                        break

        # if not isinstance(argTypes, list):

        if argType is None:
            raise Exception(
                '{}{}{}{}{}{}{}'.format(
                    '\n********************ERROR********************\n',
                    'Illegal argument value(s): \n'.format(legalArgTypes),
                    '  Argument name: {}\n'.format(argNm),
                    '  Must be (one of): {}\n'.format(legalArgTypes),
                    '  Value is: {}\n'.format(arg),
                    'File: {}  Line number: {}\n'.format(
                        self.mptCmdStruct['cmdFileNm'],
                        self.mptCmdStruct['lineNo']
                        ),
                    'Full command:\n{}\n'.format(self.mptCmdStruct['rawCmdStr'])
                    )
                )

        # List of the argter values in string form
        argVals = arg.replace('[','').replace(']','').split(',')

        # Convert the string(s) into value(s)
        if argType in [
            'Integer',
            'Positive Integer',
            'Integer List',
            'Positive Integer List',
            ]:
            rtrn = [int(x) for x in argVals]
        elif argType in [
            'Float',
            'Fuzzy Value',
            'Positive Float',
            'Float List',
            'Fuzzy Value List',
            'Positive Float List',
            ]:
            rtrn = [float(x) for x in argVals]
        elif argType in [
            'Boolean'
            ]:
            if (arg == '0' or
                re.match(r'^[Ff][Aa][Ll][Ss][Ee]$',arg)
                ):
                rtrn = [False]
            elif (arg == '1' or
                re.match(r'^[Tt][Rr][Uu][Ee]$',arg)
                ):
                rtrn = [True]

        else: # No conversion required
            rtrn = argVals
        # if argType in [...] elif...else
            
        # Convert to single value is type is not list
        if argType.find('List') < 0:
            rtrn = rtrn[0]

        return rtrn

    # def ValFromArgByNm(self,argNm):

    def _MakeCongruentNCDimArr(self,ncdimarr_to_expand,ncdimarr_to_match):
        '''
        Creates a new NCDimArr by projecting or transposing
        arr_to_expand so that it matches arr_to_match in dimensions.

        Returns None if it can't be done.
        '''

        expand_mskdarr = ncdimarr_to_expand.data.data
        expand_dim_nms = ncdimarr_to_expand.dims.keys()
        
        match_mskdarr = ncdimarr_to_match.data.data
        match_dim_nms = ncdimarr_to_match.dims.keys()

        if ac.array_is_congruent_to_tgt(expand_mskdarr,expand_dim_nms,match_mskdarr,match_dim_nms):
            
            rtrn = cp.deepcopy(ncdimarr_to_expand)

        elif ac.can_expand_to_tgt(expand_mskdarr,expand_dim_nms,match_mskdarr,match_dim_nms):
            rtrn = mpncv.NCDimensionedVar(
                dims = cp.deepcopy(ncdimarr_to_match.dims),
                data = mpncv.NCVar(
                    data = ac.expand_arr_to_match_tgt(
                        expand_mskdarr,
                        expand_dim_nms,
                        match_mskdarr,
                        match_dim_nms
                        ),
                    dim_nms = ncdimarr_to_match.dims.keys()
                    ),
                name = ncdimarr_to_expand
                )

        elif ac.can_transpose_to_tgt(expand_mskdarr,expand_dim_nms,match_mskdarr,match_dim_nms):

            rtrn = mpncv.NCDimensionedVar(
                dims = cp.deepcopy(ncdimarr_to_match.dims),
                data = mpncv.NCVar(
                    data = ac.transpose_arr_to_match_tgt(
                        expand_mskdarr,
                        expand_dim_nms,
                        match_mskdarr,
                        match_dim_nms
                        ),
                    dim_nms = ncdimarr_to_match.dims.keys()
                    ),
                name = ncdimarr_to_expand
                )

        else:

            rtrn = None

        return rtrn
        
    # def _MakeCongruentNCDimArr(self,arr_to_expand,arr_to_match):

    def _NCDimArrsAreCongruent(self,arr_1,arr_2):
        
        arr_1_mskdarr = arr_1.data.data
        arr_1_dim_nms = arr_1.dims.keys()
        arr_2_mskdarr = arr_2.data.data
        arr_2_dim_nms = arr_2.dims.keys()

        return ac.array_is_congruent_to_tgt(
            arr_1_mskdarr,
            arr_1_dim_nms,
            arr_2_mskdarr,
            arr_2_dim_nms
            )

    # def _NCDimArrsAreCongruent(self,arr_1,arr_2):

    def _CreateNCVarName(self):

        # Gets the name for the NCVar. Default behavior
        # is the variable name used in the MPilot script

        rtrn = self.ValFromArgByNm('NewNCVarName')
        
        if self.ValFromArgByNm('NewNCVarName') is None:
            rtrn = self.RsltNm()

        return rtrn
            
    # def _CreateNCVarName(self):

    def _SetMetadata(self,ncdimvar=None):

        if ncdimvar is None:
            ncdimvar = self.execRslt
        
        metadata_defs = self._ArgToList('Metadata')
        
        if metadata_defs is not None:
            for key,val in [tup.split(':') for tup in metadata_defs]:
                val = val.replace('_',' ')
                ncdimvar.set_metadata_val(key,val)

    # def _SetMetadata(self,ncdimvar):
    
# class _NetCDFUtilParent(mpfp._MPilotFxnParent):

class ReadAscVariable(_NetCDFUtilParent):
    # Reads and arc ascii file. Assumption is one variable
    
    def __init__(
        self,
        mpt_cmd_struct=None,
        ):

        super(ReadAscVariable,self).__init__(mpt_cmd_struct)
        
    # def __init__(...)

    def _SetFxnDesc(self):
        
        self.fxnDesc['DisplayName'] = 'Read Arc ASCII Variable'
        self.fxnDesc['ShortDesc'] = 'Reads an Arc ASCII file as a single variable.'
        self.fxnDesc['ReturnType'] = 'NCVar'
        
        self.fxnDesc['ReqArgs'] = {
            'InFileName':'File Name',
            'InFieldName':'Field Name'
            }
        self.fxnDesc['OptArgs'] = {
            'PrecursorFieldNames':['Field Name','Field Name List'],
            'Metadata':'Tuple List: Field Name:Any',
            'NewNCVarName':'Field Name'
            }
        
    # def _SetFxnDesc(self):

    def DependencyNms(self):
        return self._ArgToList('PrecursorFieldNames')

    def Exec(self,executedObjects):
        
        if not os.path.isfile(self.ArgByNm('InFileName')):
            raise Exception(
                '{}{}{}{}'.format(
                    '\n********************ERROR********************\n',
                    'Read file does not exist: {}\n'.format(self.ArgByNm('InFileName')),
                    'Script File: {}  Line number: {}\n'.format(
                        self.mptCmdStruct['cmdFileNm'],
                        self.mptCmdStruct['lineNo']
                        ),
                    'Full command:\n{}\n'.format(self.mptCmdStruct['rawCmdStr'])
                ),
            )
        # if not os.path.isfile(self.ArgByNm('InFileName'),'r'):

        with open(self.ArgByNm('InFileName'), 'rU') as in_f:
            lines = in_f.readlines()

        # Grab the fields that define the file

        arcasc_params = OrderedDict()
        arcasc_params['ncols'] = None
        arcasc_params['nrows'] = None
        arcasc_params['xllcorner'] = None
        arcasc_params['yllcorner'] = None
        arcasc_params['cellsize'] = None
        arcasc_params['NODATA_value'] = None

        parm_nm,parm_val = lines.pop(0).strip().split()
        if parm_nm == 'ncols':
            arcasc_params['ncols'] = int(parm_val)
        parm_nm,parm_val = lines.pop(0).strip().split()
        if parm_nm == 'nrows':
            arcasc_params['nrows'] = int(parm_val)
        parm_nm,parm_val = lines.pop(0).strip().split()
        if parm_nm == 'xllcorner':
            arcasc_params['xllcorner'] = float(parm_val)
        parm_nm,parm_val = lines.pop(0).strip().split()
        if parm_nm == 'yllcorner':
            arcasc_params['yllcorner'] = float(parm_val)
        parm_nm,parm_val = lines.pop(0).strip().split()
        if parm_nm == 'cellsize':
            arcasc_params['cellsize'] = float(parm_val)
        parm_nm,parm_val = lines.pop(0).strip().split()
        if parm_nm == 'NODATA_value':
            arcasc_params['NODATA_value'] = float(parm_val)

        if None in arcasc_params.values():
            raise Exception(
                (5*'{}\n').format(
                    '\n********************ERROR********************\n',
                    'Header in Arc Ascii file missing field(s).',
                    '  Filename: {}'.format(self.ArgByNm('InFileName')),
                    '  Header values\n    {}'.format(
                       '\n    '.join([' '.join([k,str(arcasc_params[k])]) for k in arcasc_params.keys()])
                        ),
                    'Full command:\n{}'.format(self.mptCmdStruct['rawCmdStr'])
                ),
            )

        # Create the indices
        x_ndx = [
            (arcasc_params['xllcorner'] + arcasc_params['cellsize']/2.) + n * arcasc_params['cellsize']
            for n in range(arcasc_params['ncols'])
            ]

        y_ndx = [
            (arcasc_params['yllcorner'] + arcasc_params['cellsize']/2.) + n * arcasc_params['cellsize']
            for n in range(arcasc_params['nrows'])
            ]
        y_ndx.reverse()

        x_ndx = np.ma.masked_array(x_ndx)
        y_ndx = np.ma.masked_array(y_ndx)

        if len(lines) != y_ndx.shape[0]:
            raise Exception(
                (5*'{}\n').format(
                    '\n********************ERROR********************\n',
                    'Arc Ascii file, number of lines does not match expected number of rows.',
                    '  Filename: {}'.format(self.ArgByNm('InFileName')),
                    '  Expected number of rows:    {}'.format(y_ndx.shape[0]),
                    '  Actual number of rows:    {}'.format(len(lines)),
                    'Full command:\n{}'.format(self.mptCmdStruct['rawCmdStr'])
                    ),
                )

        # Create the array to hold the values and fill it from the file
        # filling in reverse order because y was reversed.

        # ma_arr = np.ma.masked_array(
        #     np.zeros((y_ndx.shape[0],x_ndx.shape[0])),
        #     mask=np.zeros((y_ndx.shape[0],x_ndx.shape[0]))
        #     )

        ma_arr = np.ma.zeros(
            (y_ndx.shape[0],x_ndx.shape[0])
            )

        row_ndx = 0

        for in_line in lines:
            
            vals = [float(x) for x in in_line.split()]
            
            if len(vals) != x_ndx.shape[0]:
                raise Exception(
                    (5*'{}\n').format(
                        '\n********************ERROR********************\n',
                        'Arc Ascii file, number of columns does not match expected number of cols.',
                        '  Filename: {}'.format(self.ArgByNm('InFileName')),
                        '  Expected number of columns:    {}'.format(x_ndx.shape[0]),
                        '  Actual number of columns:    {}'.format(len(vals)),
                        'Full command:\n{}'.format(self.mptCmdStruct['rawCmdStr'])
                        ),
                    )

            ma_arr[row_ndx,:] = vals
            row_ndx = row_ndx+1

        # while in_line = lines.pop(-1):

        # Make the NetCDF variable

        ma_arr.mask = np.where(ma_arr == arcasc_params['NODATA_value'], True,False)
        ma_arr.data[:] = np.where(ma_arr.mask, mpncv.lookup_fill_val('d'),ma_arr.data)
        
        dims = OrderedDict()
        dims['y'] = mpncv.NCVar(
            name = 'y',
            data = np.ma.masked_array(y_ndx, dtype = float),
            dim_nms = 'y'
            )
        dims['x'] = mpncv.NCVar(
            name = 'x',
            data = np.ma.masked_array(x_ndx, dtype = float),
            dim_nms = 'x'
            )

        self.execRslt = mpncv.NCDimensionedVar(
            dims = cp.deepcopy(dims),
            name = self._CreateNCVarName(),
            data = mpncv.NCVar(
                data = ma_arr,
                dim_nms = dims.keys()
                )
            )

        self._SetMetadata()
        self.execRslt.set_metadata_val('NodeCommand',self.mptCmdStruct['rawCmdStr'])
        executedObjects[self.RsltNm()] = self

    # def Exec(self,executedObjects):

# class ReadAscVariable(_NetCDFUtilParent):
        
class ReadNCVariable(_NetCDFUtilParent):

    def __init__(
        self,
        mpt_cmd_struct=None,
        ):

        super(ReadNCVariable,self).__init__(mpt_cmd_struct)
        
    # def __init__(...)

    def _SetFxnDesc(self):
        
        self.fxnDesc['DisplayName'] = 'Read NC Variable'
        self.fxnDesc['ShortDesc'] = 'Reads a NetCDF variable from a file'
        self.fxnDesc['ReturnType'] = 'NCVar'
        
        self.fxnDesc['ReqArgs'] = {
            'InFileName':'File Name',
            'InFieldName':'Field Name'
            }
        self.fxnDesc['OptArgs'] = {
            'PrecursorFieldNames':['Field Name','Field Name List'],
            'NewNCVarName':'Field Name',
            'DimensionValueOrders':'Tuple List: Field Name:One of| HighToLow LowToHigh',
            'Metadata':'Tuple List: Field Name:Any',
            }
                    
    # def _SetFxnDesc(self):

    def DependencyNms(self):
        return self._ArgToList('PrecursorFieldNames')

    def Exec(self,executedObjects):
        
        if not os.path.isfile(self.ArgByNm('InFileName')):
            raise Exception(
                '{}{}{}{}'.format(
                    '\n********************ERROR********************\n',
                    'Read file does not exist: {}\n'.format(self.ArgByNm('InFileName')),
                    'Script File: {}  Line number: {}\n'.format(
                        self.mptCmdStruct['cmdFileNm'],
                        self.mptCmdStruct['lineNo']
                        ),
                    'Full command:\n{}\n'.format(self.mptCmdStruct['rawCmdStr'])
                ),
            )
        # if not os.path.isfile(self.ArgByNm('InFileName'),'r'):

        with nc.Dataset(self.ArgByNm('InFileName'),'r') as in_ds:
            if self.ArgByNm('InFieldName') not in in_ds.variables:
                raise Exception(
                    '{}{}{}{}{}'.format(
                        '\n********************ERROR********************\n',
                        'Read failure for file: {}\n'.format(self.ArgByNm('InFileName')),
                        '  Variable not in file: {}\n'.format(self.ArgByNm('InFieldName')),
                        'Script File: {}  Line number: {}\n'.format(
                            self.mptCmdStruct['cmdFileNm'],
                            self.mptCmdStruct['lineNo']
                            ),
                        'Full command:\n{}\n'.format(self.mptCmdStruct['rawCmdStr'])
                    ),
                )
                
            self.execRslt = mpncv.NCDimensionedVar(
                nc_ds=in_ds,
                ncvar_nm=self.ArgByNm('InFieldName')
                )
            
        # with nc.Dataset(self.ArgByNm('InFileName'),'r') as in_ds:

        # Reorder based on dimensions?
        dim_val_orders = self.ValFromArgByNm('DimensionValueOrders')
        
        if dim_val_orders is not None:
            for dim_nm,dim_order in [tup.split(':') for tup in dim_val_orders]:
                if dim_order == 'LowToHigh':
                    self.execRslt.reorder_dim_lo_to_hi(dim_nm)
                elif dim_order == 'HighToLow':
                    self.execRslt.reorder_dim_hi_to_lo(dim_nm)
            
        self.execRslt.name = self._CreateNCVarName()
        self._SetMetadata()
        self.execRslt.set_metadata_val('NodeCommand',self.mptCmdStruct['rawCmdStr'])
        executedObjects[self.RsltNm()] = self

    # def Exec(self,executedObjects):
    
# class ReadNCVariable(mpfp._MPilotFxnParent):


class ReadAndStackNCVariables(_NetCDFUtilParent):
    
    def __init__(
        self,
        mpt_cmd_struct=None,
        ):

        super(ReadAndStackNCVariables,self).__init__(mpt_cmd_struct)
        
    # def __init__(...)

    def _SetFxnDesc(self):
        
        self.fxnDesc['DisplayName'] = 'Read and stack a set of congruent NC variables into a single variable'
        self.fxnDesc['ShortDesc'] = 'Reads and stacks multiple NetCDF variables'
        self.fxnDesc['ReturnType'] = 'NCVar'
        
        self.fxnDesc['ReqArgs'] = {
            'InFileNames':['File Name','File Name List'],
            'InFieldNames':['Field Name','Field Name List'],
            'NewDimension':'Tuple: Field Name:Float:Float'
            }
        self.fxnDesc['OptArgs'] = {
            'PrecursorFieldNames':['Field Name','Field Name List'],
            'NewNCVarName':'Field Name',
            'DimensionValueOrders':'Tuple List: Field Name:One of| HighToLow LowToHigh',
            'Metadata':'Tuple List: Field Name:Any',
            }
                    
    # def _SetFxnDesc(self):

    def DependencyNms(self):
        return self._ArgToList('PrecursorFieldNames')

    def Exec(self,executedObjects):

        infile_names = self._ArgToList('InFileNames')
        infield_names = self._ArgToList('InFieldNames')

        # Check for correct filename and fieldname list lengths
        if len(infile_names) > 1 and len(infield_names) > 1:
            if len(infile_names) != len(infield_names):
                raise Exception(
                    (5*'{}\n').format(
                        '\n********************ERROR********************',
                        'FileNames and FieldNames length mismatch:',
                        'Number FileNames: {}, number FieldNames: {}'.format(len(infile_names), len(infield_names)),
                        'File: {}  Line number: {}'.format(
                            self.mptCmdStruct['cmdFileNm'],
                            self.mptCmdStruct['lineNo']
                            ),
                        'Full command:{}'.format(self.mptCmdStruct['rawCmdStr'])
                    ),
                )
            # if len(infile_names) != len(infile_names):
        # if len(infile_names) > 1 and len(infile_names) > 1:

        # make the lists of file and field names match

        if len(infile_names) > 1 and len(infield_names) == 1:
            infield_names = len(infile_names) * infield_names
            
        if len(infield_names) > 1 and len(infile_names) == 1:
            infile_names = len(infield_names) * infile_names
        

        # Read the arrays, checking for congruence

        in_dim_vars = []
        for infile_name, infield_name in zip(infile_names,infield_names):
            
            if not os.path.isfile(infile_name):
                raise Exception(
                    '{}{}{}{}'.format(
                        '\n********************ERROR********************\n',
                        'Read file does not exist: {}\n'.format(infile_name),
                        'Script File: {}  Line number: {}\n'.format(
                            self.mptCmdStruct['cmdFileNm'],
                            self.mptCmdStruct['lineNo']
                            ),
                        'Full command:\n{}\n'.format(self.mptCmdStruct['rawCmdStr'])
                    ),
                )
            # if not os.path.isfile(self.ArgByNm('InFileName'),'r'):

            with nc.Dataset(infile_name,'r') as in_ds:
                if infield_name not in in_ds.variables:
                    raise Exception(
                        '{}{}{}{}{}'.format(
                            '\n********************ERROR********************\n',
                            'Read failure for file: {}\n'.format(infile_name),
                            '  Variable not in file: {}\n'.format(infield_name),
                            'Script File: {}  Line number: {}\n'.format(
                                self.mptCmdStruct['cmdFileNm'],
                                self.mptCmdStruct['lineNo']
                                ),
                            'Full command:\n{}\n'.format(self.mptCmdStruct['rawCmdStr'])
                        ),
                    )

                new_dim_var = mpncv.NCDimensionedVar(
                    nc_ds=in_ds,
                    ncvar_nm=infield_name
                    )

            # Reorder based on dimensions?
            dim_val_orders = self.ValFromArgByNm('DimensionValueOrders')

            if dim_val_orders is not None:
                for dim_nm,dim_order in [tup.split(':') for tup in dim_val_orders]:
                    if dim_order == 'LowToHigh':
                        new_dim_var.reorder_dim_lo_to_hi(dim_nm)
                    elif dim_order == 'HighToLow':
                        new_dim_var.reorder_dim_hi_to_lo(dim_nm)

            in_dim_vars.append(new_dim_var)

        # for infile_name, infield_name in zip(infile_names,infield_names):

        # Variables are all read in, now its time to combine them along new dimension

        # Check for matching shapes

        fld_info = ''
        for infile_name, infield_name, shape in zip(
            infile_names,
            infield_names,
            [tmp_in_dnm_var.data.data.shape for tmp_in_dnm_var in in_dim_vars]
            ):

            fld_info = '{}\nInFileName: {} InFieldName: {} Shape: {}'.format(fld_info, infile_name, infield_name, shape)

        for dim_var in in_dim_vars:
            
            if dim_var.data.data.shape != in_dim_vars[0].data.data.shape:

                raise Exception(
                    (5*'{}\n').format(
                        '\n********************ERROR********************',
                        'All fields must have same shape:',
                        'Fields and shapes:{}\n'.format(fld_info),
                        'File: {}  Line number: {}'.format(
                            self.mptCmdStruct['cmdFileNm'],
                            self.mptCmdStruct['lineNo']
                            ),
                        'Full command:{}'.format(self.mptCmdStruct['rawCmdStr'])
                    ),
                )

        # create the stacked array and its mask

        rslt_arr = np.ma.masked_array(
            np.stack([dim_var.data.data.data for dim_var in in_dim_vars],axis=-1),
            mask = np.stack([dim_var.data.data.mask for dim_var in in_dim_vars],axis=-1)
            )

        # create the new dimension's array
        new_dim_nm,new_dim_start,new_dim_step = self.ValFromArgByNm('NewDimension').split(':')
        new_dim_arr = np.ma.masked_array(
            [float(new_dim_start) + x * float(new_dim_step) for x in range(rslt_arr.shape[-1])],
            mask = False
            )

        new_dim_ncvar = mpncv.NCVar(
            data = new_dim_arr,
            name = new_dim_nm,
            dim_nms = new_dim_nm
            )

        # dimensions for the rslt_arr

        rslt_arr_dims = OrderedDict()
        for dim_nm in in_dim_vars[0].dims.keys():
            rslt_arr_dims[dim_nm] = cp.deepcopy(in_dim_vars[0].dims[dim_nm])
        rslt_arr_dims[new_dim_nm] = cp.deepcopy(new_dim_ncvar)

        # Create and populate new NCDimensionedVar()

        self.execRslt = mpncv.NCDimensionedVar(
            dims = rslt_arr_dims,
            name = self._CreateNCVarName(),
            data = mpncv.NCVar(
                data = rslt_arr,
                dim_nms = rslt_arr_dims.keys()
                ),
            )

        self._SetMetadata()
        self.execRslt.set_metadata_val('NodeCommand',self.mptCmdStruct['rawCmdStr'])
        executedObjects[self.RsltNm()] = self

    # def Exec(self,executedObjects):
    
# class ReadNCVariable(mpfp._MPilotFxnParent):







class WriteNCVariables(_NetCDFUtilParent):

    def __init__(
        self,
        mpt_cmd_struct=None,
        ):

        super(WriteNCVariables,self).__init__(mpt_cmd_struct)
        
    # def __init__(...)

    def _SetFxnDesc(self):
        
        self.fxnDesc['DisplayName'] = 'Write NC Variables'
        self.fxnDesc['ShortDesc'] = 'Writes a list of NCVars to a NetCDF File'
        self.fxnDesc['ReturnType'] = 'Boolean'
        
        self.fxnDesc['ReqArgs'] = {
            'OutFileName':'File Name',
            'OutFieldNames':['Field Name','Field Name List'],
            'Overwrite':'Boolean'
            }
        self.fxnDesc['OptArgs'] = {
            'DimensionTolerance':'Float',
            'PrecursorFieldNames':['Field Name','Field Name List']
            }
        
    # def _SetFxnDesc(self):

    def DependencyNms(self):

        rtrn = self._ArgToList('OutFieldNames') + self._ArgToList('PrecursorFieldNames')
        return rtrn
    
    # def DependencyNms(self):

    def Exec(self,executedObjects):

        if self.ValFromArgByNm('Overwrite'):
            write_or_append = 'w'
        else:
            write_or_append = 'r+'
            
        try: # check for writable file
            outF = open(self.ArgByNm('OutFileName'),write_or_append)
            outF.close()
            if self.ValFromArgByNm('Overwrite'):
                os.remove(self.ArgByNm('OutFileName'))
        except IOError:
            raise Exception(
                '{}{}{}{}'.format(
                    '\n********************ERROR********************\n',
                    'Unable to open file for writing: {}\n'.format(self.ArgByNm('OutFileName')),
                    'Script File: {}  Line number: {}\n'.format(
                        self.mptCmdStruct['cmdFileNm'],
                        self.mptCmdStruct['lineNo']
                        ),
                    'Full command:\n{}\n'.format(self.mptCmdStruct['rawCmdStr'])
                ),
            )

        # try: # check for writable file

        with nc.Dataset(self.ArgByNm('OutFileName'),write_or_append) as outDS:

            # Tolerance for checking matching dimensions
            tolerance = self.ValFromArgByNm('DimensionTolerance')
            if tolerance is None:
                tolerance = 0
        
            for outFldNm in self._ArgToList('OutFieldNames'):

                if not isinstance(executedObjects[outFldNm].ExecRslt(),mpncv.NCDimensionedVar):
                    raise Exception(
                        '{}{}{}{}{}'.format(
                            '\n********************ERROR********************\n',
                            'Write failure for file: {}\n'.format(self.ArgByNm('InFileName')),
                            '  Variable is not NCVar: {}\n'.format(outFldNm),
                            'Script File: {}  Line number: {}\n'.format(
                                self.mptCmdStruct['cmdFileNm'],
                                self.mptCmdStruct['lineNo']
                                ),
                            'Full command:\n{}\n'.format(self.mptCmdStruct['rawCmdStr'])
                        ),
                    )


                executedObjects[outFldNm].execRslt.add_to_nc_ds(outDS,tolerance)

        self.execRslt = True
        executedObjects[self.RsltNm] = self

    # def Exec(self,executedObjects):
    
# class WriteNCVariables(_NetCDFUtilParent):

class ReadNCGlobalMetadata(_NetCDFUtilParent):

    def __init__(
        self,
        mpt_cmd_struct=None,
        ):

        super(ReadNCGlobalMetadata,self).__init__(mpt_cmd_struct)
        
    # def __init__(...)

    def _SetFxnDesc(self):
        
        self.fxnDesc['DisplayName'] = 'Read NC Global Metadata'
        self.fxnDesc['ShortDesc'] = 'Reads global NetCDF metadata from a file'
        self.fxnDesc['ReturnType'] = 'GlobalMetadata'
        
        self.fxnDesc['ReqArgs'] = {
            'InFileName':'File Name',
            }
        self.fxnDesc['OptArgs'] = {
            'PrecursorFieldNames':['Field Name','Field Name List']
            }
        
    # def _SetFxnDesc(self):

    def DependencyNms(self):
        return self._ArgToList('PrecursorFieldNames')

    def Exec(self,executedObjects):
        
        if not os.path.isfile(self.ArgByNm('InFileName')):
            raise Exception(
                '{}{}{}{}'.format(
                    '\n********************ERROR********************\n',
                    'Read file does not exist: {}\n'.format(self.ArgByNm('InFileName')),
                    'Script File: {}  Line number: {}\n'.format(
                        self.mptCmdStruct['cmdFileNm'],
                        self.mptCmdStruct['lineNo']
                        ),
                    'Full command:\n{}\n'.format(self.mptCmdStruct['rawCmdStr'])
                ),
            )
        # if not os.path.isfile(self.ArgByNm('InFileName'),'r'):
            
        with nc.Dataset(self.ArgByNm('InFileName'),'r') as in_ds:
            self.execRslt = cp.deepcopy(in_ds.__dict__)
        executedObjects[self.RsltNm()] = self

    # def Exec(self,executedObjects):
        
# class ReadNCGlobalMetadata(_NetCDFUtilParent):

class WriteNCGlobalMetadata(_NetCDFUtilParent):

    def __init__(
        self,
        mpt_cmd_struct=None,
        ):

        super(WriteNCGlobalMetadata,self).__init__(mpt_cmd_struct)
        
    # def __init__(...)

    def _SetFxnDesc(self):
        
        self.fxnDesc['DisplayName'] = 'Write NC Global Metadata'
        self.fxnDesc['ShortDesc'] = 'Writes global NetCDF metadata to a file'
        self.fxnDesc['ReturnType'] = 'Boolean'
        
        self.fxnDesc['ReqArgs'] = {
            'OutFileName':'File Name',
            'InFieldName':'Field Name',
            'Overwrite':'Boolean'
            }
        self.fxnDesc['OptArgs'] = {
            'PrecursorFieldNames':['Field Name','Field Name List']
            }
        
    # def _SetFxnDesc(self):

    def DependencyNms(self):
        return self._ArgToList('PrecursorFieldNames') + self._ArgToList('InFieldName')

    def Exec(self,executedObjects):

        if self.ValFromArgByNm('Overwrite'):
            write_or_append = 'w'
        else:
            write_or_append = 'r+'
            
        try: # check for writable file
            outF = open(self.ArgByNm('OutFileName'),write_or_append)
            outF.close()
            if self.ValFromArgByNm('Overwrite'):
                os.remove(self.ArgByNm('OutFileName'))
        except IOError:
            raise Exception(
                '{}{}{}{}'.format(
                    '\n********************ERROR********************\n',
                    'Unable to open file for writing: {}\n'.format(self.ArgByNm('OutFileName')),
                    'Script File: {}  Line number: {}\n'.format(
                        self.mptCmdStruct['cmdFileNm'],
                        self.mptCmdStruct['lineNo']
                        ),
                    'Full command:\n{}\n'.format(self.mptCmdStruct['rawCmdStr'])
                ),
            )

        # try: # check for writable file

        nc_atts = OrderedDict()

        with nc.Dataset(self.ArgByNm('OutFileName'),write_or_append) as out_ds:
            out_ds.setncatts(
                executedObjects[self.ValFromArgByNm('InFieldName')].ExecRslt()
                )

        self.execRslt = True
        executedObjects[self.RsltNm] = self
                        
    # def Exec(self,executedObjects):
        
# class WriteNCGlobalMetadata(_NetCDFUtilParent):

class AddToNCGlobalMetadata(_NetCDFUtilParent):

    def __init__(
        self,
        mpt_cmd_struct=None,
        ):

        super(AddToNCGlobalMetadata,self).__init__(mpt_cmd_struct)
        
    # def __init__(...)

    def _SetFxnDesc(self):
        
        self.fxnDesc['DisplayName'] = 'Add fields NC Global Metadata'
        self.fxnDesc['ShortDesc'] = 'Writes additional global NetCDF metadata to a file'
        self.fxnDesc['ReturnType'] = 'Boolean'
        
        self.fxnDesc['ReqArgs'] = {
            'OutFileName':'File Name',
            'MetadataDefinitions':'Tuple List: Any:Any',
            'Overwrite':'Boolean'
            }
        self.fxnDesc['OptArgs'] = {
            'PrecursorFieldNames':['Field Name','Field Name List']
            }
        
    # def _SetFxnDesc(self):

    def DependencyNms(self):
        return self._ArgToList('PrecursorFieldNames') + self._ArgToList('InFieldName')

    def Exec(self,executedObjects):

        if self.ValFromArgByNm('Overwrite'):
            write_or_append = 'w'
        else:
            write_or_append = 'r+'
            
        try: # check for writable file
            outF = open(self.ArgByNm('OutFileName'),write_or_append)
            outF.close()
            if self.ValFromArgByNm('Overwrite'):
                os.remove(self.ArgByNm('OutFileName'))
        except IOError:
            raise Exception(
                '{}{}{}{}'.format(
                    '\n********************ERROR********************\n',
                    'Unable to open file for writing: {}\n'.format(self.ArgByNm('OutFileName')),
                    'Script File: {}  Line number: {}\n'.format(
                        self.mptCmdStruct['cmdFileNm'],
                        self.mptCmdStruct['lineNo']
                        ),
                    'Full command:\n{}\n'.format(self.mptCmdStruct['rawCmdStr'])
                ),
            )

        # try: # check for writable file

        nc_atts = OrderedDict()

        for md_def in self._ArgToList('MetadataDefinitions'):
            key,val = md_def.split(':',2)
            val = val.replace('_',' ')
            nc_atts[key] = val
            
        with nc.Dataset(self.ArgByNm('OutFileName'),write_or_append) as out_ds:
                out_ds.setncatts(nc_atts)

        self.execRslt = True
        executedObjects[self.RsltNm] = self
                        
    # def Exec(self,executedObjects):
        
# class AddToNCGlobalMetadata(_NetCDFUtilParent):

class CopyNCVar(_NetCDFUtilParent):

    '''
    Copies the variable.
    '''
    
    def __init__(
        self,
        mpt_cmd_struct=None,
        ):

        super(CopyNCVar,self).__init__(mpt_cmd_struct)
        
    # def __init__(...)

    def _SetFxnDesc(self):
        
        self.fxnDesc['DisplayName'] = 'CopyNCVar'
        self.fxnDesc['ShortDesc'] = 'Copies the variable.'
        
        self.fxnDesc['ReturnType'] = 'NCVar'
        
        self.fxnDesc['ReqArgs'] = {
            'InFieldName':'Field Name',
            }
        self.fxnDesc['OptArgs'] = {
            'PrecursorFieldNames':['Field Name','Field Name List'],
            'Metadata':'Tuple List: Field Name:Any',
            'NewNCVarName':'Field Name',
            }
        
    # def _SetFxnDesc(self):

    def DependencyNms(self):
        return self._ArgToList('InFieldName') + \
          self._ArgToList('PrecursorFieldNames')

    def Exec(self,executedObjects):

        new_ncdimvar = cp.deepcopy(
            executedObjects[self.ValFromArgByNm('InFieldName')].ExecRslt()
            )

        new_ncdimvar.name = self._CreateNCVarName()
            
        self.execRslt = new_ncdimvar
        self._SetMetadata()
        self.execRslt.set_metadata_val('NodeCommand',self.mptCmdStruct['rawCmdStr'])
        executedObjects[self.RsltNm()] = self
  
    # def Exec(self,executedObjects):
        
# class CopyNCVar(_NetCDFUtilParent):
    
class MaskNCVarByValues(_NetCDFUtilParent):

    def __init__(
        self,
        mpt_cmd_struct=None,
        ):

        super(MaskNCVarByValues,self).__init__(mpt_cmd_struct)
        
    # def __init__(...)

    def _SetFxnDesc(self):
        
        self.fxnDesc['DisplayName'] = 'MaskNCVarByValues'
        self.fxnDesc['ShortDesc'] = 'Masks NCVar using data values'
        self.fxnDesc['ReturnType'] = 'NCVar'
        
        self.fxnDesc['ReqArgs'] = {
            'InFieldName':'Field Name',
            'MaskInOrOut':'One of| In Out',
            'Values':['Float', 'Float List', 'Integer', 'Integer List']
            }
        self.fxnDesc['OptArgs'] = {
            'PrecursorFieldNames':['Field Name','Field Name List'],
            'Metadata':'Tuple List: Field Name:Any',
            'NewNCVarName':'Field Name'
            }
        
    # def _SetFxnDesc(self):

    def DependencyNms(self):
        return self._ArgToList('PrecursorFieldNames') + \
          self._ArgToList('InFieldName')

    def Exec(self,executedObjects):

        src_dimarr = executedObjects[self.ValFromArgByNm('InFieldName')].ExecRslt()
        src_marr = src_dimarr.data.data

        if self.ValFromArgByNm('MaskInOrOut') == 'In':

            # Start with everything false and mask in using values.

            new_mask = np.ma.make_mask_none(src_marr.shape)
            new_mask[:] = True

            for val in self._ScalarToList(self.ValFromArgByNm('Values')):

                # change mask array in place
                new_mask[np.ma.where(src_marr == int(val))] = False

        else: # mask_action == 'Out'    

            if isinstance(src_marr.mask,np.bool_):
                new_mask = np.ma.make_mask_none(src_marr.shape)
                new_mask[:] = src_marr.mask
            else:
                new_mask = cp.deepcopy(src_marr.mask)
            
            for val in self._ScalarToList(self.ValFromArgByNm('Values')):

                # change mask array in place
                new_mask[np.ma.where(src_marr == int(val))] = True

        # if self.ValFromArgByNm('MaskInOrOut') == 'In':

        new_marr = np.ma.masked_array(
            cp.deepcopy(src_marr.data),
            mask=new_mask
            )

        # Make the dimensioned array object
        self.execRslt = mpncv.NCDimensionedVar(
            dims = cp.deepcopy(src_dimarr.dims),
            name = self._CreateNCVarName(),
            data = mpncv.NCVar(
                data = new_marr,
                dim_nms = src_dimarr.dims.keys()
                )
            )
            
        self._SetMetadata()
        self.execRslt.set_metadata_val('NodeCommand',self.mptCmdStruct['rawCmdStr'])
        executedObjects[self.RsltNm()] = self

    # def Exec(self,executedObjects):
        
# class MaskNCVarByValues(_NetCDFUtilParent):

class MaskNCVarByValueRanges(_NetCDFUtilParent):

    def __init__(
        self,
        mpt_cmd_struct=None,
        ):

        super(MaskNCVarByValueRanges,self).__init__(mpt_cmd_struct)
        
    # def __init__(...)

    def _SetFxnDesc(self):
        
        self.fxnDesc['DisplayName'] = 'MaskNCVarByValueRanges'
        self.fxnDesc['ShortDesc'] = 'Masks NCVar using data value ranges'
        self.fxnDesc['ReturnType'] = 'NCVar'
        
        self.fxnDesc['ReqArgs'] = {
            'InFieldName':'Field Name',
            'MaskInOrOut':'One of| In Out',
            'ValueRanges':'Tuple List: Float:Float',
            }
        self.fxnDesc['OptArgs'] = {
            'PrecursorFieldNames':['Field Name','Field Name List'],
            'Metadata':'Tuple List: Field Name:Any',
            'NewNCVarName':'Field Name'
            }
        
    # def _SetFxnDesc(self):

    def DependencyNms(self):
        return self._ArgToList('PrecursorFieldNames') + \
          self._ArgToList('InFieldName')

    def Exec(self,executedObjects):

        src_dimarr = executedObjects[self.ValFromArgByNm('InFieldName')].ExecRslt()
        src_marr = src_dimarr.data.data
        
        if self.ValFromArgByNm('MaskInOrOut') == 'In':

            new_mask = np.ma.make_mask_none(src_marr.shape)
            new_mask[:] = True

            for tup in self._ArgToList('ValueRanges'):

                mask_val1,mask_val2 =''.join(tup.split()).split(':')
                mask_min_val = min(float(mask_val1),float(mask_val2))
                mask_max_val = max(float(mask_val1),float(mask_val2))

                new_mask[
                    np.ma.where(
                        np.ma.logical_and(src_marr >= mask_min_val, src_marr <= mask_max_val)
                        )
                    ] = False
        
        else: # masking out

            if isinstance(src_marr.mask,np.bool_):
                new_mask = np.ma.make_mask_none(src_marr.shape)
                new_mask[:] = src_marr.mask # gets the boolean value
            else:
                new_mask = cp.deepcopy(src_marr.mask) # gets the array
        
            for tup in self._ArgToList('ValueRanges'):

                mask_val1,mask_val2 =''.join(tup.split()).split(':')
                mask_min_val = min(float(mask_val1),float(mask_val2))
                mask_max_val = max(float(mask_val1),float(mask_val2))

                new_mask[
                    np.ma.where(
                        np.ma.logical_and(src_marr >= mask_min_val, src_marr <= mask_max_val)
                        )
                    ] = True

        # if self.ValFromArgByNm('MaskInOrOut') == 'In':

        new_marr = np.ma.masked_array(
            cp.deepcopy(src_marr.data),
            mask=new_mask
            )

        # Make the dimensioned array object
        self.execRslt = mpncv.NCDimensionedVar(
            dims = cp.deepcopy(src_dimarr.dims),
            name = self._CreateNCVarName(),
            data = mpncv.NCVar(
                data = new_marr,
                dim_nms = src_dimarr.dims.keys()
                )
            )
            
        self._SetMetadata()
        self.execRslt.set_metadata_val('NodeCommand',self.mptCmdStruct['rawCmdStr'])
        executedObjects[self.RsltNm()] = self

    # def Exec(self,executedObjects):
        
# class MaskNCVarByValueRanges(_NetCDFUtilParent):

class ClipExtentToUnmasked(_NetCDFUtilParent):

    '''
    Clips the variable to the extent of those cells
    that are not masked to True, i.e. cells that are
    not masked out.
    '''
    
    def __init__(
        self,
        mpt_cmd_struct=None,
        ):

        super(ClipExtentToUnmasked,self).__init__(mpt_cmd_struct)
        
    # def __init__(...)

    def _SetFxnDesc(self):
        
        self.fxnDesc['DisplayName'] = 'ClipExtentToUnmasked'
        self.fxnDesc['ShortDesc'] = 'Clips NCVar to extent of unmasked cells.'
        self.fxnDesc['ReturnType'] = 'NCVar'
        
        self.fxnDesc['ReqArgs'] = {
            'InFieldName':'Field Name',
            }
        self.fxnDesc['OptArgs'] = {
            'PrecursorFieldNames':['Field Name','Field Name List'],
            'Metadata':'Tuple List: Field Name:Any',
            'NewNCVarName':'Field Name'
            }
        
    # def _SetFxnDesc(self):

    def DependencyNms(self):
        return self._ArgToList('PrecursorFieldNames') + \
        self._ArgToList('InFieldName')

    def Exec(self,executedObjects):

        in_NC_dim_var = executedObjects[self.ValFromArgByNm('InFieldName')].ExecRslt()
        
        # Determine the extent slices of the dimension
        
        # tuple of indices of all unmasked cells. One tuple value per
        # index. Each tuple value is an array of index values along that
        # dimension.
        masked_in_ndxs = np.where(in_NC_dim_var.data.data.mask == 0)

        # This gets the extents for each dimension for the masked_in cells
        extent_slices = []
        for dim_ndx_arr in masked_in_ndxs:

            extent_slices.append(slice(dim_ndx_arr.min(),dim_ndx_arr.max()+1))
            
        # The NCVars for the indexes:
        new_dims = OrderedDict()
        # for dim_nm, dim_slice in zip(new_NC_dim_var.data.dim_nms, extent_slices):
        for dim_nm, dim_slice in zip(in_NC_dim_var.data.dim_nms, extent_slices):
            in_dim = in_NC_dim_var.dims[dim_nm]
            new_dims[dim_nm] = mpncv.NCVar(
            data = cp.deepcopy(in_dim.data[dim_slice]),
            dim_nms = cp.deepcopy(in_dim.dim_nms),
            name = dim_nm
            )

        # Create the new NCDimensionedVar for the clipped array

        self.execRslt = mpncv.NCDimensionedVar(
            dims = new_dims,
            name = self._CreateNCVarName(),
            data = mpncv.NCVar(
                data = cp.deepcopy(in_NC_dim_var.data.data[extent_slices]),
                dim_nms = new_dims.keys(),
                )
            )

        self._SetMetadata()
        self.execRslt.set_metadata_val('NodeCommand',self.mptCmdStruct['rawCmdStr'])
        executedObjects[self.RsltNm()] = self
        
    # def Exec(self,executedObjects):
        
# class ClipExtentToUnmasked(_NetCDFUtilParent):


# class MinimumExtentFromNCVars(_NetCDFUtilParent):

#     '''
#     Process a set of NCVar or NCDimensionedVar objects and
#     finds the minimum extent across them. Returns None if
#     there is no overlap in one of the dimensions considered.
#     '''
    
#     def __init__(
#         self,
#         mpt_cmd_struct=None,
#         ):

#         super(ClipToMatchTemplate,self).__init__(mpt_cmd_struct)
        
#     # def __init__(...)

#     def _SetFxnDesc(self):
        
#         self.fxnDesc['DisplayName'] = 'MinimumExtentFromNCVars'
#         self.fxnDesc['ShortDesc'] = 'Gets the minimum extent from a set of NCVars'
        
#         self.fxnDesc['ReturnType'] = 'Extents'
        
#         self.fxnDesc['ReqArgs'] = {
#             'InFieldNames':'Field Name List',
#             'DimensionNames':'Any'
#             }
#         self.fxnDesc['OptArgs'] = {
#             'PrecursorFieldNames':['Field Name','Field Name List'],
#             }
        
#     # def _SetFxnDesc(self):

#     def DependencyNms(self):
#         return self._ArgToList('PrecursorFieldNames') + \
#           self._ArgToList('InFieldNames')

#     def Exec(self,executedObjects):

#         # ncdimvar is a NCDimensionedVar
#         # dim_ncvar is the ncvar of a dimension

#         extents = mpnc.extents()
#         for fld_nm in self._ArgToList('InFieldNames'):
#             fld_extents = executedObjects[fld_nm].extents
#             for dim_nm in self._ArgToList('DimensionNames'):
#                 if dim_nm in extents.extents.keys():
#                     extents.set_extent(dim_nm,(float('inf'),float(-inf)
                                               
#                 extent_min = min(extents.extents[dim_nm][0],fld_extents[dim_nm][0])
#                 extent_max = max(extents.extents[dim_nm][1],fld_extents[dim_nm][1])
#                 extents.set_extent(dim_nm,extent_min,extent_max)

#                 # heeere
                    
                    

#         comp_tolerance = self.ValFromArgByNm('ComparisonTolerance')
#         if comp_tolerance is None: comp_tolerance = 0
            
#         to_clip_ncdimvar = executedObjects[self.ValFromArgByNm('FieldNameToClip')].ExecRslt()
#         tmplt_ncdimvar = executedObjects[self.ValFromArgByNm('TemplateFieldName')].ExecRslt()
        
#         # Get the extents from the template NCVar
#         clip_dim_extents = OrderedDict()
#         for clip_dim_nm,dim_ncvar in  tmplt_ncdimvar.dims.items():
#             clip_dim_extents[clip_dim_nm] = {
#                 'min':dim_ncvar.data.min() - comp_tolerance,
#                 'max':dim_ncvar.data.max() + comp_tolerance
#                 }

#         # Get the slices for trimming the NCVar that will be clipped
#         to_clip_slices = []
#         for to_clip_dim_nm, to_clip_dim_ncvar in to_clip_ncdimvar.dims.items():

#             if to_clip_dim_nm in clip_dim_extents.keys():

#                 dim_ndx_arr = np.where(
#                     np.logical_and(
#                         to_clip_dim_ncvar.data >= clip_dim_extents[to_clip_dim_nm]['min'],
#                         to_clip_dim_ncvar.data <= clip_dim_extents[to_clip_dim_nm]['max'],
#                         )
#                     )[0]

#                 to_clip_slices.append(slice(dim_ndx_arr.min(),dim_ndx_arr.max()+1))

#                 slice_len = to_clip_slices[-1].stop - to_clip_slices[-1].start

#                 if slice_len != tmplt_ncdimvar.dims[to_clip_dim_nm].shape[0]:
#                     raise Exception(
#                         (3*'{}\n').format(
#                             'Grid spacing mismatch when clipping to template',
#                             '  TemplateFieldName, NCVar name: {}  {}'.format(
#                                 self.ValFromArgByNm('TemplateFieldName'),
#                                 tmplt_ncdimvar.name
#                                 ),
#                             '  FieldNameToClip, NCVar name:   {}  {}'.format(
#                                 self.ValFromArgByNm('FieldNameToClip'),
#                                 to_clip_ncdimvar.name
#                                 )
#                             )
#                         ) # raise Exception(...)

#             else:
            
#                 to_clip_slices.append(slice(0,len(to_clip_dim_ncvar.data)))

#             # if to_clip_dim_nm in clip_dim_extents.keys():...else...
                
#         # for to_clip_dim_nm, to_clip_dim_ncvar in to_clip_ncdimvar.dims.items():

#         # Get name of new NCDimensionedVar
#         if self.ValFromArgByNm('NewNCVarName') is not None:
#             new_ncdimvar_nm = self.ValFromArgByNm('NewNCVarName')
#         else:
#             new_ncdimvar_nm = to_clip_ncdimvar.name

#         # Create the new clipped NCDimensionedVar

#         new_ncdimvar = mpncv.NCDimensionedVar(name=new_ncdimvar_nm)
        
#         # The data NCVar:
#         new_ncdimvar.data = mpncv.NCVar(
#             data = cp.deepcopy(to_clip_ncdimvar.data.data[to_clip_slices]),
#             dim_nms = cp.deepcopy(to_clip_ncdimvar.data.dim_nms),
#             name = new_ncdimvar_nm
#             )

#         # The NCVars for the indexes:
#         new_dims = OrderedDict()
#         do_dim_copy_from_tmplt = self.ValFromArgByNm('CopyTemplateDimensions') is None or \
#             self.ValFromArgByNm('CopyTemplateDimensions') == True
            
#         for dim_nm, dim_slice in zip(new_ncdimvar.data.dim_nms, to_clip_slices):
#             # Use the template dimensions for the dimensions of the new NCVar
#             # this is the recommended default
#             if dim_nm in tmplt_ncdimvar.dims and do_dim_copy_from_tmplt:
                
#                 new_dims[dim_nm] = cp.deepcopy(tmplt_ncdimvar.dims[dim_nm])
                
#             else: # clip the dimension of the toclip NCVar
                
#                 in_dim = to_clip_ncdimvar.dims[dim_nm]
#                 new_dims[dim_nm] = mpncv.NCVar(
#                 data = cp.deepcopy(in_dim.data[dim_slice]),
#                 dim_nms = cp.deepcopy(in_dim.dim_nms),
#                 name = dim_nm
#                 )
                
#         # for dim_nm, dim_slice in zip(new_ncdimvar.data.dim_nms, to_clip_slices):

#         new_ncdimvar.dims = new_dims

#         self.execRslt = new_ncdimvar
#         executedObjects[self.RsltNm()] = self
        
#     # def Exec(self,executedObjects):
    
# # class MinExtentFromNCVars(_NetCDFUtilParent):

class ClipToMatchTemplate(_NetCDFUtilParent):

    '''
    Clips the variable to the extent of those cells
    that are not masked to True, i.e. cells that are
    not masked out.

    Uses index values in the template NCVar to select valid
    index values in the NCVar to clip. Checks to make sure
    that the size of the resulting indices match those of the
    template.
    '''
    
    def __init__(
        self,
        mpt_cmd_struct=None,
        ):

        super(ClipToMatchTemplate,self).__init__(mpt_cmd_struct)
        
    # def __init__(...)

    def _SetFxnDesc(self):
        
        self.fxnDesc['DisplayName'] = 'ClipToMatchTemplate'
        self.fxnDesc['ShortDesc'] = (5*'\n{}').format(
            '    - Clips extents of FieldNameToClip to match those of TemplateFieldName,',
            '      doing value comparison with optional comparison tolerance.',
            '    - Insures that index sizes match.',
            '    - Replaces indices of clipped data with those of the template data unless',
            '      otherwise specified.'
            )
        
        self.fxnDesc['ReturnType'] = 'NCVar'
        
        self.fxnDesc['ReqArgs'] = {
            'TemplateFieldName':'Field Name',
            'FieldNameToClip':'Field Name'
            }
        self.fxnDesc['OptArgs'] = {
            'PrecursorFieldNames':['Field Name','Field Name List'],
            'Metadata':'Tuple List: Field Name:Any',
            'NewNCVarName':'Field Name',
            'CopyTemplateDimensions':'Boolean',
            'ComparisonTolerance':'Float'
            }
        
    # def _SetFxnDesc(self):

    def DependencyNms(self):
        return self._ArgToList('PrecursorFieldNames') + \
          self._ArgToList('TemplateFieldName') + \
          self._ArgToList('FieldNameToClip') 

    def Exec(self,executedObjects):

        # ncdimvar is a NCDimensionedVar
        # dim_ncvar is the ncvar of a dimension

        comp_tolerance = self.ValFromArgByNm('ComparisonTolerance')
        if comp_tolerance is None: comp_tolerance = 0
            
        to_clip_ncdimvar = executedObjects[self.ValFromArgByNm('FieldNameToClip')].ExecRslt()
        tmplt_ncdimvar = executedObjects[self.ValFromArgByNm('TemplateFieldName')].ExecRslt()
        
        # Get the extents from the template NCVar
        clip_dim_extents = OrderedDict()
        for clip_dim_nm,dim_ncvar in  tmplt_ncdimvar.dims.items():
            clip_dim_extents[clip_dim_nm] = {
                'min':dim_ncvar.data.min() - comp_tolerance,
                'max':dim_ncvar.data.max() + comp_tolerance
                }

        # Get the slices for trimming the NCVar that will be clipped
        to_clip_slices = []
        for to_clip_dim_nm, to_clip_dim_ncvar in to_clip_ncdimvar.dims.items():

            if to_clip_dim_nm in clip_dim_extents.keys():

                dim_ndx_arr = np.where(
                    np.logical_and(
                        to_clip_dim_ncvar.data >= clip_dim_extents[to_clip_dim_nm]['min'],
                        to_clip_dim_ncvar.data <= clip_dim_extents[to_clip_dim_nm]['max'],
                        )
                    )[0]

                to_clip_slices.append(slice(dim_ndx_arr.min(),dim_ndx_arr.max()+1))

                slice_len = to_clip_slices[-1].stop - to_clip_slices[-1].start

                if slice_len != tmplt_ncdimvar.dims[to_clip_dim_nm].shape[0]:
                    raise Exception(
                        (3*'{}\n').format(
                            'Grid spacing mismatch when clipping to template',
                            '  TemplateFieldName, NCVar name: {}  {}'.format(
                                self.ValFromArgByNm('TemplateFieldName'),
                                tmplt_ncdimvar.name
                                ),
                            '  FieldNameToClip, NCVar name:   {}  {}'.format(
                                self.ValFromArgByNm('FieldNameToClip'),
                                to_clip_ncdimvar.name
                                )
                            )
                        ) # raise Exception(...)

            else:
            
                to_clip_slices.append(slice(0,len(to_clip_dim_ncvar.data)))

            # if to_clip_dim_nm in clip_dim_extents.keys():...else...
                
        # for to_clip_dim_nm, to_clip_dim_ncvar in to_clip_ncdimvar.dims.items():

        # The NCVars for the indexes:
        new_dims = OrderedDict()
        do_dim_copy_from_tmplt = self.ValFromArgByNm('CopyTemplateDimensions') is None or \
            self.ValFromArgByNm('CopyTemplateDimensions') == True
            
        for dim_nm, dim_slice in zip(to_clip_ncdimvar.data.dim_nms, to_clip_slices):
            # Use the template dimensions for the dimensions of the new NCVar
            # this is the recommended default
            if dim_nm in tmplt_ncdimvar.dims and do_dim_copy_from_tmplt:
                
                new_dims[dim_nm] = cp.deepcopy(tmplt_ncdimvar.dims[dim_nm])
                
            else: # clip the dimension of the toclip NCVar
                
                in_dim = to_clip_ncdimvar.dims[dim_nm]
                new_dims[dim_nm] = mpncv.NCVar(
                data = cp.deepcopy(in_dim.data[dim_slice]),
                dim_nms = cp.deepcopy(in_dim.dim_nms),
                name = dim_nm
                )
                
        # for dim_nm, dim_slice in zip(to_clip_ncdimvar.data.dim_nms, to_clip_slices):

        # Create the new NCDimensionedVar for the clipped array

        self.execRslt = mpncv.NCDimensionedVar(
            dims = new_dims,
            name = self._CreateNCVarName(),
            data = mpncv.NCVar(
                data = cp.deepcopy(to_clip_ncdimvar.data.data[to_clip_slices]),
                dim_nms = new_dims.keys(),
                )
            )

        self._SetMetadata()
        self.execRslt.set_metadata_val('NodeCommand',self.mptCmdStruct['rawCmdStr'])
        executedObjects[self.RsltNm()] = self

    # def Exec(self,executedObjects):
        
# class ClipToMatchTemplate(_NetCDFUtilParent):

class ReplaceMask(_NetCDFUtilParent):

    '''
    Replaces the mask in one NCVar with a copy of the mask
    from another.
    '''
    
    def __init__(
        self,
        mpt_cmd_struct=None,
        ):

        super(ReplaceMask,self).__init__(mpt_cmd_struct)
        
    # def __init__(...)

    def _SetFxnDesc(self):
        
        self.fxnDesc['DisplayName'] = 'ReplaceMask'
        self.fxnDesc['ShortDesc'] = 'Replaces the mask in one NCVar with a copy of the mask from another.'        
        self.fxnDesc['ReturnType'] = 'NCVar'
        
        self.fxnDesc['ReqArgs'] = {
            'TargetFieldName':'Field Name',
            'FromFieldName':'Field Name'
            }
        self.fxnDesc['OptArgs'] = {
            'PrecursorFieldNames':['Field Name','Field Name List'],
            'Metadata':'Tuple List: Field Name:Any',
            'NewNCVarName':'Field Name',
            }
        
    # def _SetFxnDesc(self):

    def DependencyNms(self):
        return self._ArgToList('PrecursorFieldNames') + \
          self._ArgToList('TargetFieldName') + \
          self._ArgToList('FromFieldName') 

    def Exec(self,executedObjects):
        
        tgt_ncdimvar = executedObjects[self.ValFromArgByNm('TargetFieldName')].ExecRslt()
        new_ncdimvar = cp.deepcopy(tgt_ncdimvar)
        new_ncdimvar.name = self._CreateNCVarName()
        from_ncdimvar = executedObjects[self.ValFromArgByNm('FromFieldName')].ExecRslt()

        tgt_mask_arr = tgt_ncdimvar.data.data.mask
        tgt_dim_nms = tgt_ncdimvar.dims.keys()
        from_mask_arr = from_ncdimvar.data.data.mask
        from_dim_nms = from_ncdimvar.dims.keys()
        
        # Checks and manipulations to make this work

        if ac.array_is_congruent_to_tgt(
            from_mask_arr,
            from_dim_nms,
            tgt_mask_arr,
            tgt_dim_nms
            ):

            new_ncdimvar.data.data.mask = cp.deepcopy(from_mask_arr)

        elif ac.can_expand_to_tgt(
            from_mask_arr,
            from_dim_nms,
            tgt_mask_arr,
            tgt_dim_nms
            ):

            new_ncdimvar.data.data.mask = ac.expand_arr_to_match_tgt(
                from_mask_arr,
                from_dim_nms,
                tgt_mask_arr,
                tgt_dim_nms
                )

        elif ac.can_transpose_to_tgt(
            from_mask_arr,
            from_dim_nms,
            tgt_mask_arr,
            tgt_dim_nms
            ):

            new_ncdimvar.data.data.mask = ac.transpose_arr_to_match_tgt(
                from_mask_arr,
                from_dim_nms,
                tgt_mask_arr,
                tgt_dim_nms
                )

        else:
            
            raise Exception(
                (3*'{}\n').format(
                    'Mask from FromFieldName cannot be made to match TargetFieldName mask shape.',
                    '  FromFieldName, NCVar name: {}  {}'.format(
                        self.ValFromArgByNm('FromFieldName'),
                        from_ncdimvar.name
                        ),
                    '  FieldNameToClip, NCVar name:   {}  {}'.format(
                        self.ValFromArgByNm('TargetFieldName'),
                        tgt_ncdimvar.name
                        )
                    )
                ) # raise Exception(...)

        # Rename it?
                
        self.execRslt = new_ncdimvar
        self._SetMetadata()
        self.execRslt.set_metadata_val('NodeCommand',self.mptCmdStruct['rawCmdStr'])
        executedObjects[self.RsltNm()] = self
        
    # def Exec(self,executedObjects):
        
# class ReplaceMask(_NetCDFUtilParent):

class AddToMask(_NetCDFUtilParent):

    '''
    Replaces the mask in one NCVar with a copy of the mask
    from another.
    '''
    
    def __init__(
        self,
        mpt_cmd_struct=None,
        ):

        super(AddToMask,self).__init__(mpt_cmd_struct)
        
    # def __init__(...)

    def _SetFxnDesc(self):
        
        self.fxnDesc['DisplayName'] = 'AddToMask'
        self.fxnDesc['ShortDesc'] = 'Combines from masks with the target mask in the target NCVar.'
        self.fxnDesc['ReturnType'] = 'NCVar'
        
        self.fxnDesc['ReqArgs'] = {
            'TargetFieldName':'Field Name',
            'FromFieldNames':['Field Name','Field Name List']
            }
        self.fxnDesc['OptArgs'] = {
            'PrecursorFieldNames':['Field Name','Field Name List'],
            'Metadata':'Tuple List: Field Name:Any',
            'NewNCVarName':'Field Name',
            }
        
    # def _SetFxnDesc(self):

    def DependencyNms(self):
        return self._ArgToList('PrecursorFieldNames') + \
          self._ArgToList('TargetFieldName') + \
          self._ArgToList('FromFieldNames') 

    def Exec(self,executedObjects):

        tgt_ncdimvar = executedObjects[self.ValFromArgByNm('TargetFieldName')].ExecRslt()
        new_ncdimvar = cp.deepcopy(tgt_ncdimvar)
        new_ncdimvar.name = self._CreateNCVarName()        
        tgt_mask_arr = tgt_ncdimvar.data.data.mask
        tgt_dim_nms = tgt_ncdimvar.dims.keys()
            
        for from_fld_nm in self._ArgToList('FromFieldNames'):

            from_ncdimvar = executedObjects[from_fld_nm].ExecRslt()

            from_mask_arr = from_ncdimvar.data.data.mask
            from_dim_nms = from_ncdimvar.dims.keys()

            # Checks and manipulations to make this work

            if ac.array_is_congruent_to_tgt(
                from_mask_arr,
                from_dim_nms,
                tgt_mask_arr,
                tgt_dim_nms
                ):

                tmp_mask = from_mask_arr

            elif ac.can_expand_to_tgt(
                from_mask_arr,
                from_dim_nms,
                tgt_mask_arr,
                tgt_dim_nms
                ):

                tmp_mask = ac.expand_arr_to_match_tgt(
                    from_mask_arr,
                    from_dim_nms,
                    tgt_mask_arr,
                    tgt_dim_nms
                    )

            elif ac.can_transpose_to_tgt(
                from_mask_arr,
                from_dim_nms,
                tgt_mask_arr,
                tgt_dim_nms
                ):

                tmp_mask = ac.transpose_arr_to_match_tgt(
                    from_mask_arr,
                    from_dim_nms,
                    tgt_mask_arr,
                    tgt_dim_nms
                    )

            else:

                raise Exception(
                    (3*'{}\n').format(
                        'Mask from FromFieldNames cannot be made to match TargetFieldName mask shape.',
                        '  FromFieldName, NCVar name, shape: {}  {}  {}'.format(
                            self.ValFromArgByNm('FromFieldNames'),
                            from_ncdimvar.name,
                            from_ncdimvar.data.shape,
                            ),
                        '  TargetFieldName, NCVar name:   {}  {}  {}'.format(
                            self.ValFromArgByNm('TargetFieldName'),
                            tgt_ncdimvar.name,
                            tgt_ncdimvar.data.shape
                            )
                        )
                    ) # raise Exception(...)

            new_ncdimvar.data.data.mask = np.ma.mask_or(new_ncdimvar.data.data.mask,tmp_mask)

        # for from_fld_nm in self._ArgToList('FromFieldNames'):

        self.execRslt = new_ncdimvar
        self._SetMetadata()
        self.execRslt.set_metadata_val('NodeCommand',self.mptCmdStruct['rawCmdStr'])
        executedObjects[self.RsltNm()] = self
        
    # def Exec(self,executedObjects):
        
# class AddToMask(_NetCDFUtilParent):

class MaskByTemplate(_NetCDFUtilParent):

    '''Adds to the mask of a variable using another variable as a template.'''

    def __init__(
        self,
        mpt_cmd_struct=None,
        ):

        super(MaskByTemplate,self).__init__(mpt_cmd_struct)
        
    # def __init__(...)

    def _SetFxnDesc(self):
        
        self.fxnDesc['DisplayName'] = 'Mask By Template'
        self.fxnDesc['ShortDesc'] = 'Adds to the mask of a variable using one or more variables as templates.'
        self.fxnDesc['ReturnType'] = 'NCVar'
        
        self.fxnDesc['ReqArgs'] = {
            'TargetFieldName':'Field Name',
            'FromFieldNames':['Field Name','Field Name List'],
            'MaskInOrOut':'One of| In Out'
            }
        self.fxnDesc['OptArgs'] = {
            'MaskValue':'Float',
            'PrecursorFieldNames':['Field Name','Field Name List'],
            'Metadata':'Tuple List: Field Name:Any',
            'NewNCVarName':'Field Name',
            }
        
    # def _SetFxnDesc(self):

    def DependencyNms(self):
        return self._ArgToList('PrecursorFieldNames') + \
          self._ArgToList('TargetFieldName') + \
          self._ArgToList('FromFieldNames') 

    def Exec(self,executedObjects):

        tgt_ncdimvar = executedObjects[self.ValFromArgByNm('TargetFieldName')].ExecRslt()
        new_ncdimvar = cp.deepcopy(tgt_ncdimvar)
        from_dim_nms = new_ncdimvar.dims.keys()
        new_ncdimvar.name = self._CreateNCVarName()        
        tgt_mask_arr = tgt_ncdimvar.data.data.mask
        tgt_dim_nms = tgt_ncdimvar.dims.keys()

        for from_fld_nm in self._ArgToList('FromFieldNames'):

            from_arr = executedObjects[from_fld_nm].ExecRslt().data.data

            # Create the mask to add

            if self.ValFromArgByNm('MaskValue') is not None:
                add_mask = np.ma.where(from_arr == self.ValFromArgByNm('MaskValue'), False, True)
            else:
                add_mask = cp.deepcopy(from_arr.mask)

            if self.ValFromArgByNm('MaskInOrOut') == 'Out':
                add_mask = np.where(add_mask == True, False, True)

            if ac.array_is_congruent_to_tgt(
                add_mask,
                from_dim_nms,
                tgt_mask_arr,
                tgt_dim_nms
                ):

                pass # add mask is good

            elif ac.can_expand_to_tgt(
                add_mask,
                from_dim_nms,
                tgt_mask_arr,
                tgt_dim_nms
                ):

                add_mask = ac.expand_arr_to_match_tgt(
                    add_mask,
                    from_dim_nms,
                    tgt_mask_arr,
                    tgt_dim_nms
                    )

            elif ac.can_transpose_to_tgt(
                add_mask,
                from_dim_nms,
                tgt_mask_arr,
                tgt_dim_nms
                ):

                tmp_mask = ac.transpose_arr_to_match_tgt(
                    add_mask,
                    from_dim_nms,
                    tgt_mask_arr,
                    tgt_dim_nms
                    )

            else:

                raise Exception(
                    (3*'{}\n').format(
                        'Mask from FromFieldName cannot be made to match TargetFieldName mask shape.',
                        '  FromFieldName, NCVar name, shape: {}  {}  {}'.format(
                            self.ValFromArgByNm('FromFieldNames'),
                            from_ncdimvar.name,
                            from_ncdimvar.data.shape,
                            ),
                        '  TargetFieldName, NCVar name:   {}  {}  {}'.format(
                            self.ValFromArgByNm('TargetFieldName'),
                            tgt_ncdimvar.name,
                            tgt_ncdimvar.data.shape
                            )
                        )
                    ) # raise Exception(...)

            new_ncdimvar.data.data.mask = np.ma.mask_or(new_ncdimvar.data.data.mask,add_mask)

        # for from_fld_nm in self._ArgToList('FromFieldNames'):
        
        self.execRslt = new_ncdimvar
        self._SetMetadata()
        self.execRslt.set_metadata_val('NodeCommand',self.mptCmdStruct['rawCmdStr'])
        executedObjects[self.RsltNm()] = self
        
    # def Exec(self,executedObjects):

# class MaskByTemplate(_NetCDFUtilParent):

class MaskByMatchingValues(_NetCDFUtilParent):

    '''Adds to the mask of a variable using another variable as a template.'''

    def __init__(
        self,
        mpt_cmd_struct=None,
        ):

        super(MaskByMatchingValues,self).__init__(mpt_cmd_struct)
        
    # def __init__(...)

    def _SetFxnDesc(self):
        
        self.fxnDesc['DisplayName'] = 'Mask By Matching Values'
        self.fxnDesc['ShortDesc'] = 'Masks an array where values match another array. Also combines masks.'
        self.fxnDesc['ReturnType'] = 'NCVar'
        
        self.fxnDesc['ReqArgs'] = {
            'InFieldName':'Field Name',
            'MatchFieldName':'Field Name',
            'MaskInOrOut':'One of| In Out'
            }
        self.fxnDesc['OptArgs'] = {
            'PrecursorFieldNames':['Field Name','Field Name List'],
            'Metadata':'Tuple List: Field Name:Any',
            'NewNCVarName':'Field Name',
            }
        
    # def _SetFxnDesc(self):

    def DependencyNms(self):
        return self._ArgToList('PrecursorFieldNames') + \
          self._ArgToList('InFieldName') + \
          self._ArgToList('MatchFieldName') 

    def Exec(self,executedObjects):

        new_ncdimvar = cp.deepcopy(executedObjects[self.ValFromArgByNm('InFieldName')].ExecRslt())
        new_ncdimvar.name = self._CreateNCVarName()        

        value_arr = self._MakeCongruentNCDimArr(
            executedObjects[self.ValFromArgByNm('MatchFieldName')].ExecRslt(),
            new_ncdimvar
            ).data.data

        if self.ValFromArgByNm('MaskInOrOut') == 'In':
            np.ma.masked_where(
                value_arr != new_ncdimvar.data.data,
                new_ncdimvar.data.data,
                copy = False
                )
        else:
            np.ma.masked_where(
                value_arr == new_ncdimvar.data.data,
                new_ncdimvar.data.data,
                copy = False
                )

        self.execRslt = new_ncdimvar
        self._SetMetadata()
        self.execRslt.set_metadata_val('NodeCommand',self.mptCmdStruct['rawCmdStr'])
        executedObjects[self.RsltNm()] = self
        
    # def Exec(self,executedObjects):
    
# class MaskByMatchingValues(_NetCDFUtilParent):

class ClipByDimensionValueRanges(_NetCDFUtilParent):

    '''
    Clips the variable to the extent of the specified
    index values. e.g. lat:44.75:49.00 would trim the
    extent to cells 49.00 <= lat index value <= 49.00.
    Optional padding value can be used to insure that
    rounding doesn't place desired results out of range.
    '''
    
    def __init__(
        self,
        mpt_cmd_struct=None,
        ):

        super(ClipByDimensionValueRanges,self).__init__(mpt_cmd_struct)
        
    # def __init__(...)

    def _SetFxnDesc(self):
        
        self.fxnDesc['DisplayName'] = 'ClipByIndexValueRanges'
        self.fxnDesc['ShortDesc'] = 'Clips array to extent calculated for value ranges (e.g. lat:44.75:49.00).'
        
        self.fxnDesc['ReturnType'] = 'NCVar'
        
        self.fxnDesc['ReqArgs'] = {
            'InFieldName':'Field Name',
            'Extents':'Tuple List: Field Name:Float:Float'
            }
        self.fxnDesc['OptArgs'] = {
            'PrecursorFieldNames':['Field Name','Field Name List'],
            'Metadata':'Tuple List: Field Name:Any',
            'NewNCVarName':'Field Name',
            'ComparisonTolerance':'Float'
            }
        
    # def _SetFxnDesc(self):

    def DependencyNms(self):
        return self._ArgToList('InFieldName') + \
          self._ArgToList('PrecursorFieldNames')

    def Exec(self,executedObjects):

        comp_tolerance = self.ValFromArgByNm('ComparisonTolerance')
        if comp_tolerance is None: comp_tolerance = 0
            
        to_clip_ncdimvar = executedObjects[self.ValFromArgByNm('InFieldName')].ExecRslt()

        # Get the ranges for the dimensions
        clip_dim_extents = OrderedDict()
        for tup in self._ArgToList('Extents'):
            
            # strip white space and split
            dim_nm,dim_val1,dim_val2 =''.join(tup.split()).split(':')
            dim_min_val = min(float(dim_val1),float(dim_val2)) - comp_tolerance
            dim_max_val = max(float(dim_val1),float(dim_val2)) + comp_tolerance

            if dim_nm not in to_clip_ncdimvar.dims:
                raise Exception(
                    (3*'{}\n').format(
                        'Error: Non-existent dimension specified: {}'.format(dim_nm),
                        '  Result name: {}'.format(self.RsltNm()),
                        '  Dimension name: {}'.format(dim_nm)
                        )
                    )

            # Error if dim specified twice
            if dim_nm in clip_dim_extents.keys():

                raise Exception(
                    (3*'{}\n').format(
                        'Error: Each dimension may only be specified one time',
                        '  Result name: {}'.format(self.RsltNm()),
                        '  Dimension name: {}'.format(dim_nm)
                        )
                    )

            clip_dim_extents[dim_nm] = {
                'min':dim_min_val,
                'max':dim_max_val
                }

        # Get the slices for trimming the NCVar that will be clipped
        to_clip_slices = []
        for to_clip_dim_nm, to_clip_dim_ncvar in to_clip_ncdimvar.dims.items():

            if to_clip_dim_nm in clip_dim_extents.keys():

                dim_ndx_arr = np.where(
                    np.logical_and(
                        to_clip_dim_ncvar.data >= clip_dim_extents[to_clip_dim_nm]['min'],
                        to_clip_dim_ncvar.data <= clip_dim_extents[to_clip_dim_nm]['max'],
                        )
                    )[0]

                to_clip_slices.append(slice(dim_ndx_arr.min(),dim_ndx_arr.max()+1))

            else:
            
                to_clip_slices.append(slice(0,len(to_clip_dim_ncvar.data)))

            # if to_clip_dim_nm in clip_dim_extents.keys():...else...

        # for to_clip_dim_nm, to_clip_dim_ncvar in to_clip_ncdimvar.dims.items():

        # The NCVars for the indexes:
        new_dims = OrderedDict()
        do_dim_copy_from_tmplt = self.ValFromArgByNm('CopyTemplateDimensions') is None or \
            self.ValFromArgByNm('CopyTemplateDimensions') == True
            
        for dim_nm, dim_slice in zip(to_clip_ncdimvar.data.dim_nms, to_clip_slices):
            
            in_dim = to_clip_ncdimvar.dims[dim_nm]
            new_dims[dim_nm] = mpncv.NCVar(
            data = cp.deepcopy(in_dim.data[dim_slice]),
            dim_nms = cp.deepcopy(in_dim.dim_nms),
            name = dim_nm
            )
 
        # Create the new NCDimensionedVar for the clipped array

        self.execRslt = mpncv.NCDimensionedVar(
            dims = new_dims,
            name = self._CreateNCVarName(),
            data = mpncv.NCVar(
                data = cp.deepcopy(to_clip_ncdimvar.data.data[to_clip_slices]),
                dim_nms = new_dims.keys(),
                )
            )

        self._SetMetadata()
        self.execRslt.set_metadata_val('NodeCommand',self.mptCmdStruct['rawCmdStr'])
        executedObjects[self.RsltNm()] = self
        
    # def Exec(self,executedObjects):

# class ClipByDimensionValueRanges(_NetCDFUtilParent):

class ClipByDimensionIndexRanges(_NetCDFUtilParent):

    '''
    Clips the variable to the extent of the specified
    index ranges. e.g. lat:11:15 would trim the
    extent to cells to the lat[11:16] (Note: it is inclusive)
    '''
    
    def __init__(
        self,
        mpt_cmd_struct=None,
        ):

        super(ClipByDimensionIndexRanges,self).__init__(mpt_cmd_struct)
        
    # def __init__(...)

    def _SetFxnDesc(self):
        
        self.fxnDesc['DisplayName'] = 'ClipByIndexRanges'
        self.fxnDesc['ShortDesc'] = 'Clips array to extent indexes (not values).'
        
        self.fxnDesc['ReturnType'] = 'NCVar'
        
        self.fxnDesc['ReqArgs'] = {
            'InFieldName':'Field Name',
            'Extents':'Tuple List: Field Name:Integer:Integer'
            }
        self.fxnDesc['OptArgs'] = {
            'PrecursorFieldNames':['Field Name','Field Name List'],
            'Metadata':'Tuple List: Field Name:Any',
            'NewNCVarName':'Field Name',
            }
        
    # def _SetFxnDesc(self):

    def DependencyNms(self):
        return self._ArgToList('InFieldName') + \
          self._ArgToList('PrecursorFieldNames')

    def Exec(self,executedObjects):

        comp_tolerance = self.ValFromArgByNm('ComparisonTolerance')
        if comp_tolerance is None: comp_tolerance = 0
            
        to_clip_ncdimvar = executedObjects[self.ValFromArgByNm('InFieldName')].ExecRslt()

        # Get the ranges for the dimensions
        # used to order slices to match order of ncvar indices
        to_clip_slices = OrderedDict()
        for dim_nm,dim_len in zip(to_clip_ncdimvar.dim_nms,to_clip_ncdimvar.shape):
            to_clip_slices[dim_nm] = slice(0,dim_len)

        used_dim_nms = []
        for tup in self._ArgToList('Extents'):
            
            # strip white space and split
            dim_nm,ndx_val1,ndx_val2 =''.join(tup.split()).split(':')
            dim_min_ndx_val = min(int(ndx_val1),int(ndx_val2))
            dim_max_ndx_val = max(int(ndx_val1),int(ndx_val2)) + 1

            if dim_nm not in to_clip_ncdimvar.dims:
                raise Exception(
                    (3*'{}\n').format(
                        'Error: Non-existent dimension specified: {}'.format(dim_nm),
                        '  Result name: {}'.format(self.RsltNm()),
                        '  Dimension name: {}'.format(dim_nm)
                        )
                    )

            # Error if dim specified twice
            if dim_nm in used_dim_nms:

                raise Exception(
                    (3*'{}\n').format(
                        'Error: Each dimension may only be specified one time',
                        '  Result name: {}'.format(self.RsltNm()),
                        '  Dimension name: {}'.format(dim_nm)
                        )
                    )
            else:
                used_dim_nms.append(dim_nm)

            to_clip_slices[dim_nm] = slice(dim_min_ndx_val,dim_max_ndx_val)

        # for tup in self._ArgToList('Extents'):

        # The NCVars for the indexes:
        new_dims = OrderedDict()
            
        for dim_nm, dim_slice in to_clip_slices.items():
            
            in_dim = to_clip_ncdimvar.dims[dim_nm]
            new_dims[dim_nm] = mpncv.NCVar(
            data = cp.deepcopy(in_dim.data[dim_slice]),
            dim_nms = cp.deepcopy(in_dim.dim_nms),
            name = dim_nm
            )

        # Create the new NCDimensionedVar for the clipped array

        self.execRslt = mpncv.NCDimensionedVar(
            dims = new_dims,
            name = self._CreateNCVarName(),
            data = mpncv.NCVar(
                data = cp.deepcopy(to_clip_ncdimvar.data.data[to_clip_slices.values()]),
                dim_nms = new_dims.keys(),
                )
            )

        self._SetMetadata()
        self.execRslt.set_metadata_val('NodeCommand',self.mptCmdStruct['rawCmdStr'])
        executedObjects[self.RsltNm()] = self

    # def Exec(self,executedObjects):

# class ClipByDimensionIndexRanges(_NetCDFUtilParent):

class MaskByDimensionPatterns(_NetCDFUtilParent):

    '''
    Masks the variable using by repeating specified mask
    over a dimension. For example if you wanted to get
    only the April,May,June values for a variable with a
    month dimension, and the variable starts with January
    you would specify

    Patterns = [month:111000111111]

    This pattern will be replicated across the month
    dimension. Where there are 1's the data will be masked
    out.
    '''
    
    def __init__(
        self,
        mpt_cmd_struct=None,
        ):

        super(MaskByDimensionPatterns,self).__init__(mpt_cmd_struct)
        
    # def __init__(...)

    def _SetFxnDesc(self):
        
        self.fxnDesc['DisplayName'] = 'MaskByDimensionPatterns'
        self.fxnDesc['ShortDesc'] = 'Masks out data based on where 1s are in the mask pattern.'
        
        self.fxnDesc['ReturnType'] = 'NCVar'
        
        self.fxnDesc['ReqArgs'] = {
            'InFieldName':'Field Name',
            'Patterns':'Tuple List: Field Name:Binary'
            }
        self.fxnDesc['OptArgs'] = {
            'PrecursorFieldNames':['Field Name','Field Name List'],
            'Metadata':'Tuple List: Field Name:Any',
            'NewNCVarName':'Field Name',
            }
        
    # def _SetFxnDesc(self):

    def DependencyNms(self):
        return self._ArgToList('InFieldName') + \
          self._ArgToList('PrecursorFieldNames')

    def Exec(self,executedObjects):

        new_ncdimvar = cp.deepcopy(
            executedObjects[self.ValFromArgByNm('InFieldName')].ExecRslt()
            )
        new_ncdimvar.name = self._CreateNCVarName()

        new_ncdimvar_mask = new_ncdimvar.data.data.mask
        
        for pattern_spec in self._ArgToList('Patterns'):

            '''
            For each pattern spec:
            - Create a 1-d pattern array with the same length as the dimension
            - Generate a transpose of the mask array so that the dimension
              of the pattern spec is the last dimension. This allows for array
              broadcast.
            - Do a mask_or using the 1-d pattern array
            '''

            dim_nm, mask_pattern = pattern_spec.split(':')

            if dim_nm not in new_ncdimvar.dim_nms:
                raise Exception(
                    (3*'{}\n').format(
                        'Error: Non-existent dimension specified: {}'.format(dim_nm),
                        '  Result name: {}'.format(self.RsltNm()),
                        '  Dimension name: {}'.format(dim_nm)
                        )
                    )

            # Mask repeated to length of dimension
            mask_pattern = [int(d) for d in mask_pattern]
            dim = new_ncdimvar.dims[dim_nm]
            # add 1 because of integer arithmetic, then trim to dimension length
            dim_mask = mask_pattern * (len(dim.data) / len(mask_pattern) + 1)
            dim_mask = np.array(dim_mask[:len(dim.data)],dtype=bool)

            # Do the transpose and add to the mask

            axis = new_ncdimvar.dim_nms.index(dim_nm)

            # create the new dimension index order 
            forward_transpose = range(new_ncdimvar_mask.ndim)
            forward_transpose[axis:-1] = forward_transpose[axis+1:]
            forward_transpose[-1] = axis

            backward_transpose = range(new_ncdimvar_mask.ndim)
            backward_transpose[axis+1:] = backward_transpose[axis:-1]
            backward_transpose[axis] = len(backward_transpose) - 1
            
            # transpose forward, do masking, then transpose back to
            # original order
            trans_mask = new_ncdimvar_mask.transpose(forward_transpose)
            trans_mask = np.ma.mask_or(trans_mask,dim_mask)
            new_ncdimvar.data.data.mask = trans_mask.transpose(backward_transpose)
            
        # for pattern_spec in self._ArgToLists('Patterns'):

        # if self.ValFromArgByNm('NewNCVarName') is not None:
        #     new_ncdimvar.name = self.ValFromArgByNm('NewNCVarName')
            
        self.execRslt = new_ncdimvar
        self._SetMetadata()
        self.execRslt.set_metadata_val('NodeCommand',self.mptCmdStruct['rawCmdStr'])
        executedObjects[self.RsltNm()] = self
        
# class MaskByDimensionPatterns(_NetCDFUtilParent):

class UnmaskNCVar(_NetCDFUtilParent):

    '''
    Unmasks the variable. Actually, it sets the mask to false
    so that no cells in the variable are masked out.
    '''
    
    def __init__(
        self,
        mpt_cmd_struct=None,
        ):

        super(UnmaskNCVar,self).__init__(mpt_cmd_struct)
        
    # def __init__(...)

    def _SetFxnDesc(self):
        
        self.fxnDesc['DisplayName'] = 'UnmaskNCVar'
        self.fxnDesc['ShortDesc'] = 'Masks in all data in the variable.'
        
        self.fxnDesc['ReturnType'] = 'NCVar'
        
        self.fxnDesc['ReqArgs'] = {
            'InFieldName':'Field Name',
            }
        self.fxnDesc['OptArgs'] = {
            'PrecursorFieldNames':['Field Name','Field Name List'],
            'Metadata':'Tuple List: Field Name:Any',
            'NewNCVarName':'Field Name',
            }
        
    # def _SetFxnDesc(self):

    def DependencyNms(self):
        return self._ArgToList('InFieldName') + \
          self._ArgToList('PrecursorFieldNames')

    def Exec(self,executedObjects):

        new_ncdimvar = cp.deepcopy(
            executedObjects[self.ValFromArgByNm('InFieldName')].ExecRslt()
            )

        new_ncdimvar.data.data.mask[:] = False
        new_ncdimvar.name = self._CreateNCVarName()
            
        self.execRslt = new_ncdimvar
        self._SetMetadata()
        self.execRslt.set_metadata_val('NodeCommand',self.mptCmdStruct['rawCmdStr'])
        executedObjects[self.RsltNm()] = self
  
# class UnmaskNCVar(_NetCDFUtilParent):

class UnmaskNCVarWithValue(_NetCDFUtilParent):

    '''
    Unmasks the variable. Actually, it sets the mask to false
    so that no cells in the variable are masked out.

    Sets the value of the formerly masked cells to ReplacementValue.
    
    '''
    
    def __init__(
        self,
        mpt_cmd_struct=None,
        ):

        super(UnmaskNCVarWithValue,self).__init__(mpt_cmd_struct)
        
    # def __init__(...)

    def _SetFxnDesc(self):
        
        self.fxnDesc['DisplayName'] = 'UnmaskNCVarWithValue'
        self.fxnDesc['ShortDesc'] = 'Masks in all data in the variable and sets formerly masked cells to ReplacementValue.'
        
        self.fxnDesc['ReturnType'] = 'NCVar'
        
        self.fxnDesc['ReqArgs'] = {
            'InFieldName':'Field Name',
            'ReplacementValue':'Float'
            }
        self.fxnDesc['OptArgs'] = {
            'PrecursorFieldNames':['Field Name','Field Name List'],
            'Metadata':'Tuple List: Field Name:Any',
            'NewNCVarName':'Field Name',
            }
        
    # def _SetFxnDesc(self):

    def DependencyNms(self):
        return self._ArgToList('InFieldName') + \
          self._ArgToList('PrecursorFieldNames')

    def Exec(self,executedObjects):

        new_ncdimvar = cp.deepcopy(
            executedObjects[self.ValFromArgByNm('InFieldName')].ExecRslt()
            )

        new_ncdimvar.data.data.data[:] = np.where(
            new_ncdimvar.data.data.mask,
            self.ValFromArgByNm('ReplacementValue'),
            new_ncdimvar.data.data.data
            )

        new_ncdimvar.data.data.mask[:] = False

        new_ncdimvar.name = self._CreateNCVarName()
        
        self.execRslt = new_ncdimvar
        self._SetMetadata()
        self.execRslt.set_metadata_val('NodeCommand',self.mptCmdStruct['rawCmdStr'])
        executedObjects[self.RsltNm()] = self
  
# class UnmaskNCVarWithValue(_NetCDFUtilParent):

################################################################################
# Summary over dimensions operations
################################################################################

class _SummaryOverDimensions(_NetCDFUtilParent):

    def __init__(
        self,
        mpt_cmd_struct=None,
        ):

        super(_SummaryOverDimensions,self).__init__(mpt_cmd_struct)
        
    '''
    This should only be used as a parent class for classes
    that do a summary over one or more dimensions.
    '''

    def _Exec(self,executedObjects,summary_type):

        parent_ncdimvar = executedObjects[self.ValFromArgByNm('InFieldName')].ExecRslt()

        # Insure specified dimensions exist
        # and get dimension indices to operate over

        valid_dim_nms = parent_ncdimvar.dims.keys()

        # Some classes using this parent only allow one dimension
        if 'Dimensions' in self.FxnReqArgs():
            dim_nms = self._ArgToList('Dimensions')
        elif 'Dimension' in self.FxnReqArgs():
            dim_nms = self._ArgToList('Dimension')            
            
        dim_ndxs = []
        for dim_nm in dim_nms:
            
            if dim_nm in valid_dim_nms:
                dim_ndxs.append(valid_dim_nms.index(dim_nm))
            else:
                raise Exception(
                    (3*'{}\n').format(
                        'Error: Invalid dimension specified.',
                        '  Result name: {}'.format(self.RsltNm()),
                        '  Dimension name: {}'.format(dim_nm)
                        )
                    )

        # for dim_nm in self._ArgToList('Dimensions'):

        # Create dimensions for new variable. Keep correct
        # dimension order
        new_dim_ncvars = OrderedDict()
        for dim_nm in parent_ncdimvar.dims.keys():
            if dim_nm not in dim_nms:
                new_dim_ncvars[dim_nm] = cp.deepcopy(parent_ncdimvar.dims[dim_nm])

        # Create and populate new NCDimensionedVar()
        new_ncdimvar = mpncv.NCDimensionedVar(
            dims = new_dim_ncvars,
            name = self._CreateNCVarName()
            )

        if summary_type == 'sum':
            data_nparr = parent_ncdimvar.data.data.sum(axis=tuple(dim_ndxs))
        elif summary_type == 'mean':
            data_nparr = parent_ncdimvar.data.data.mean(axis=tuple(dim_ndxs))
        elif summary_type == 'max':
            data_nparr = parent_ncdimvar.data.data.max(axis=tuple(dim_ndxs))
        elif summary_type == 'min':
            data_nparr = parent_ncdimvar.data.data.min(axis=tuple(dim_ndxs))
        elif summary_type == 'std':
            data_nparr = parent_ncdimvar.data.data.std(axis=tuple(dim_ndxs))
        elif summary_type == 'var': # variance
            data_nparr = parent_ncdimvar.data.data.var(axis=tuple(dim_ndxs))
        elif summary_type == 'count':
            data_nparr = parent_ncdimvar.data.data.count(axis=tuple(dim_ndxs))
        elif summary_type == 'mode': # mode
            data_nparr = _mode_by_axes(parent_ncdimvar.data.data,tuple(dim_ndxs))
        elif summary_type == 'median': # median
            data_nparr = np.ma.median(parent_ncdimvar.data.data,axis=tuple(dim_ndxs))
            
        elif summary_type == 'dimension_val_of_max':
            
            data_nparr = np.ma.masked_array(
                parent_ncdimvar.data.data.argmax(axis=dim_ndxs[0]),
                parent_ncdimvar.data.data.max(axis=dim_ndxs[0]).mask, #cheat to get mask
                dtype = float
                )

            dim_arr = parent_ncdimvar.dims.values()[dim_ndxs[0]].data.data
            
            for arr_val_ndx in range(dim_arr.shape[0]):
                data_nparr[np.ma.where(data_nparr == arr_val_ndx)] = dim_arr[arr_val_ndx]

        elif summary_type == 'dimension_val_of_min':
            
            data_nparr = np.ma.masked_array(
                parent_ncdimvar.data.data.argmin(axis=dim_ndxs[0]),
                parent_ncdimvar.data.data.min(axis=dim_ndxs[0]).mask, #cheat to get mask
                dtype = float
                )

            dim_arr = parent_ncdimvar.dims.values()[dim_ndxs[0]].data.data
            
            for arr_val_ndx in range(dim_arr.shape[0]):
                data_nparr[np.ma.where(data_nparr == arr_val_ndx)] = dim_arr[arr_val_ndx]

        else:
            raise Exception(
                (3*'{}\n').format(
                    'Error: Illegal summary operation.',
                    '  Summary operation: {}'.format(summary_type)
                    )
                ) # raise Exception(...)
            
        new_ncdimvar.data = mpncv.NCVar(
            data = data_nparr,
            dim_nms = new_dim_ncvars.keys(),
            name = self._CreateNCVarName()
            )
            
        self.execRslt = new_ncdimvar
        self._SetMetadata()
        self.execRslt.set_metadata_val('NodeCommand',self.mptCmdStruct['rawCmdStr'])
        executedObjects[self.RsltNm()] = self
        
    # def _Exec(self,executedObjects,summary_type):
        
# class _SummaryOverDimensions(_NetCDFUtilParent)
        
class MeanOverDimensions(_SummaryOverDimensions):

    def __init__(
        self,
        mpt_cmd_struct=None,
        ):

        super(MeanOverDimensions,self).__init__(mpt_cmd_struct)
        
    # def __init__(...)

    def _SetFxnDesc(self):
        
        self.fxnDesc['DisplayName'] = 'MeanOverDimensions'
        self.fxnDesc['ShortDesc'] = 'Takes the mean of an NCVar over one or more dimensions.'
        
        self.fxnDesc['ReturnType'] = 'NCVar, Float, or Integer'
        
        self.fxnDesc['ReqArgs'] = {
            'InFieldName':'Field Name',
            'Dimensions':['Field Name','Field Name List']
            }
        self.fxnDesc['OptArgs'] = {
            'PrecursorFieldNames':['Field Name','Field Name List'],
            'Metadata':'Tuple List: Field Name:Any',
            'NewNCVarName':'Field Name',
            }
        
    # def _SetFxnDesc(self):

    def DependencyNms(self):
        return self._ArgToList('InFieldName') + \
          self._ArgToList('PrecursorFieldNames')

    def Exec(self,executedObjects):

        self._Exec(executedObjects,'mean')
        
    # def Exec(self,executedObjects):

# class MeanOverDimensions(_NetCDFUtilParent):
    
class SumOverDimensions(_SummaryOverDimensions):

    def __init__(
        self,
        mpt_cmd_struct=None,
        ):

        super(SumOverDimensions,self).__init__(mpt_cmd_struct)
        
    # def __init__(...)

    def _SetFxnDesc(self):
        
        self.fxnDesc['DisplayName'] = 'SumOverDimensions'
        self.fxnDesc['ShortDesc'] = 'Sums an NCVar over one or more dimensions.'
        
        self.fxnDesc['ReturnType'] = 'NCVar, Float, or Integer'
        
        self.fxnDesc['ReqArgs'] = {
            'InFieldName':'Field Name',
            'Dimensions':['Field Name','Field Name List']
            }
        self.fxnDesc['OptArgs'] = {
            'PrecursorFieldNames':['Field Name','Field Name List'],
            'Metadata':'Tuple List: Field Name:Any',
            'NewNCVarName':'Field Name',
            }
        
    # def _SetFxnDesc(self):

    def DependencyNms(self):
        return self._ArgToList('InFieldName') + \
          self._ArgToList('PrecursorFieldNames')

    def Exec(self,executedObjects):

        self._Exec(executedObjects,'sum')
        
    # def Exec(self,executedObjects):

# class SumOverDimensions(_NetCDFUtilParent):
    
class MinOverDimensions(_SummaryOverDimensions):

    def __init__(
        self,
        mpt_cmd_struct=None,
        ):

        super(MinOverDimensions,self).__init__(mpt_cmd_struct)
        
    # def __init__(...)

    def _SetFxnDesc(self):
        
        self.fxnDesc['DisplayName'] = 'MinOverDimensions'
        self.fxnDesc['ShortDesc'] = 'Minimum of an NCVar over one or more dimensions.'
        
        self.fxnDesc['ReturnType'] = 'NCVar, Float, or Integer'
        
        self.fxnDesc['ReqArgs'] = {
            'InFieldName':'Field Name',
            'Dimensions':['Field Name','Field Name List']
            }
        self.fxnDesc['OptArgs'] = {
            'PrecursorFieldNames':['Field Name','Field Name List'],
            'Metadata':'Tuple List: Field Name:Any',
            'NewNCVarName':'Field Name',
            }
        
    # def _SetFxnDesc(self):

    def DependencyNms(self):
        return self._ArgToList('InFieldName') + \
          self._ArgToList('PrecursorFieldNames')

    def Exec(self,executedObjects):

        self._Exec(executedObjects,'min')
        
    # def Exec(self,executedObjects):

# class MinOverDimensions(_NetCDFUtilParent):
    
class MaxOverDimensions(_SummaryOverDimensions):

    def __init__(
        self,
        mpt_cmd_struct=None,
        ):

        super(MaxOverDimensions,self).__init__(mpt_cmd_struct)
        
    # def __init__(...)

    def _SetFxnDesc(self):
        
        self.fxnDesc['DisplayName'] = 'MaxOverDimensions'
        self.fxnDesc['ShortDesc'] = 'Maximum of an NCVar over one or more dimensions.'
        
        self.fxnDesc['ReturnType'] = 'NCVar, Float, or Integer'
        
        self.fxnDesc['ReqArgs'] = {
            'InFieldName':'Field Name',
            'Dimensions':['Field Name','Field Name List']
            }
        self.fxnDesc['OptArgs'] = {
            'PrecursorFieldNames':['Field Name','Field Name List'],
            'Metadata':'Tuple List: Field Name:Any',
            'NewNCVarName':'Field Name',
            }
        
    # def _SetFxnDesc(self):

    def DependencyNms(self):
        return self._ArgToList('InFieldName') + \
          self._ArgToList('PrecursorFieldNames')

    def Exec(self,executedObjects):

        self._Exec(executedObjects,'max')
        
    # def Exec(self,executedObjects):

# class MaxOverDimensions(_NetCDFUtilParent):
    
class StdDevOverDimensions(_SummaryOverDimensions):

    def __init__(
        self,
        mpt_cmd_struct=None,
        ):

        super(StdDevOverDimensions,self).__init__(mpt_cmd_struct)
        
    # def __init__(...)

    def _SetFxnDesc(self):
        
        self.fxnDesc['DisplayName'] = 'StdDevOverDimensions'
        self.fxnDesc['ShortDesc'] = 'Standard deviation of an NCVar over one or more dimensions.'
        
        self.fxnDesc['ReturnType'] = 'NCVar, Float, or Integer'
        
        self.fxnDesc['ReqArgs'] = {
            'InFieldName':'Field Name',
            'Dimensions':['Field Name','Field Name List']
            }
        self.fxnDesc['OptArgs'] = {
            'PrecursorFieldNames':['Field Name','Field Name List'],
            'Metadata':'Tuple List: Field Name:Any',
            'NewNCVarName':'Field Name',
            }
        
    # def _SetFxnDesc(self):

    def DependencyNms(self):
        return self._ArgToList('InFieldName') + \
          self._ArgToList('PrecursorFieldNames')

    def Exec(self,executedObjects):

        self._Exec(executedObjects,'std')
        
    # def Exec(self,executedObjects):

# class StdDevOverDimensions(_NetCDFUtilParent):
    
class VarianceOverDimensions(_SummaryOverDimensions):

    def __init__(
        self,
        mpt_cmd_struct=None,
        ):

        super(VarianceOverDimensions,self).__init__(mpt_cmd_struct)
        
    # def __init__(...)

    def _SetFxnDesc(self):
        
        self.fxnDesc['DisplayName'] = 'VarianceOverDimensions'
        self.fxnDesc['ShortDesc'] = 'Variance of an NCVar over one or more dimensions.'
        
        self.fxnDesc['ReturnType'] = 'NCVar, Float, or Integer'
        
        self.fxnDesc['ReqArgs'] = {
            'InFieldName':'Field Name',
            'Dimensions':['Field Name','Field Name List']
            }
        self.fxnDesc['OptArgs'] = {
            'PrecursorFieldNames':['Field Name','Field Name List'],
            'Metadata':'Tuple List: Field Name:Any',
            'NewNCVarName':'Field Name',
            }
        
    # def _SetFxnDesc(self):

    def DependencyNms(self):
        return self._ArgToList('InFieldName') + \
          self._ArgToList('PrecursorFieldNames')

    def Exec(self,executedObjects):

        self._Exec(executedObjects,'var')
        
    # def Exec(self,executedObjects):

# class VarianceOverDimensions(_NetCDFUtilParent):

class CountOverDimensions(_SummaryOverDimensions):

    def __init__(
        self,
        mpt_cmd_struct=None,
        ):

        super(CountOverDimensions,self).__init__(mpt_cmd_struct)
        
    # def __init__(...)

    def _SetFxnDesc(self):
        
        self.fxnDesc['DisplayName'] = 'CountOverDimensions'
        self.fxnDesc['ShortDesc'] = 'Count of an NCVar over one or more dimensions.'
        
        self.fxnDesc['ReturnType'] = 'NCVar, or Integer'
        
        self.fxnDesc['ReqArgs'] = {
            'InFieldName':'Field Name',
            'Dimensions':['Field Name','Field Name List']
            }
        self.fxnDesc['OptArgs'] = {
            'PrecursorFieldNames':['Field Name','Field Name List'],
            'Metadata':'Tuple List: Field Name:Any',
            'NewNCVarName':'Field Name',
            }
        
    # def _SetFxnDesc(self):

    def DependencyNms(self):
        return self._ArgToList('InFieldName') + \
          self._ArgToList('PrecursorFieldNames')

    def Exec(self,executedObjects):

        self._Exec(executedObjects,'count')
        
    # def Exec(self,executedObjects):

# class StdDevOverDimensions(_NetCDFUtilParent):

class ModeOverDimensions(_SummaryOverDimensions):

    def __init__(
        self,
        mpt_cmd_struct=None,
        ):

        super(ModeOverDimensions,self).__init__(mpt_cmd_struct)
        
    # def __init__(...)

    def _SetFxnDesc(self):
        
        self.fxnDesc['DisplayName'] = 'ModeOverDimensions'
        self.fxnDesc['ShortDesc'] = 'Mode of an NCVar over one or more dimensions.'
        
        self.fxnDesc['ReturnType'] = 'NCVar, Float, or Integer'
        
        self.fxnDesc['ReqArgs'] = {
            'InFieldName':'Field Name',
            'Dimensions':['Field Name','Field Name List']
            }
        self.fxnDesc['OptArgs'] = {
            'PrecursorFieldNames':['Field Name','Field Name List'],
            'Metadata':'Tuple List: Field Name:Any',
            'NewNCVarName':'Field Name',
            }
        
    # def _SetFxnDesc(self):

    def DependencyNms(self):
        return self._ArgToList('InFieldName') + \
          self._ArgToList('PrecursorFieldNames')

    def Exec(self,executedObjects):

        self._Exec(executedObjects,'mode')
        
    # def Exec(self,executedObjects):

# class ModeOverDimensions(_NetCDFUtilParent):

class MedianOverDimensions(_SummaryOverDimensions):

    def __init__(
        self,
        mpt_cmd_struct=None,
        ):

        super(MedianOverDimensions,self).__init__(mpt_cmd_struct)
        
    # def __init__(...)

    def _SetFxnDesc(self):
        
        self.fxnDesc['DisplayName'] = 'MedianOverDimensions'
        self.fxnDesc['ShortDesc'] = 'Takes the median of an NCVar over one or more dimensions.'
        
        self.fxnDesc['ReturnType'] = 'NCVar, Float, or Integer'
        
        self.fxnDesc['ReqArgs'] = {
            'InFieldName':'Field Name',
            'Dimensions':['Field Name','Field Name List']
            }
        self.fxnDesc['OptArgs'] = {
            'PrecursorFieldNames':['Field Name','Field Name List'],
            'Metadata':'Tuple List: Field Name:Any',
            'NewNCVarName':'Field Name',
            }
        
    # def _SetFxnDesc(self):

    def DependencyNms(self):
        return self._ArgToList('InFieldName') + \
          self._ArgToList('PrecursorFieldNames')

    def Exec(self,executedObjects):

        self._Exec(executedObjects,'median')
        
    # def Exec(self,executedObjects):

# class MedianOverDimensions(_NetCDFUtilParent):

class DimensionValueOfMaxOverDimension(_SummaryOverDimensions):
    
    def __init__(
        self,
        mpt_cmd_struct=None,
        ):

        super(DimensionValueOfMaxOverDimension,self).__init__(mpt_cmd_struct)
        
    # def __init__(...)

    def _SetFxnDesc(self):
        
        self.fxnDesc['DisplayName'] = 'Dimension Value of Max Over Dimension'
        self.fxnDesc['ShortDesc'] = 'Dimension value of maximum of an NCVar over a dimension  (e.g. which layer had the max).'
        
        self.fxnDesc['ReturnType'] = 'NCVar, Float, or Integer'
        
        self.fxnDesc['ReqArgs'] = {
            'InFieldName':'Field Name',
            'Dimension':'Field Name',
            }
        self.fxnDesc['OptArgs'] = {
            'PrecursorFieldNames':['Field Name','Field Name List'],
            'Metadata':'Tuple List: Field Name:Any',
            'NewNCVarName':'Field Name',
            }
        
    # def _SetFxnDesc(self):

    def DependencyNms(self):
        return self._ArgToList('InFieldName') + \
          self._ArgToList('PrecursorFieldNames')

    def Exec(self,executedObjects):

        self._Exec(executedObjects,'dimension_val_of_max')
        
    # def Exec(self,executedObjects):

#class DimensionValueOfMaxOverDimension(_SummaryOverDimensions):

class DimensionValueOfMinOverDimension(_SummaryOverDimensions):
    
    def __init__(
        self,
        mpt_cmd_struct=None,
        ):

        super(DimensionValueOfMinOverDimension,self).__init__(mpt_cmd_struct)
        
    # def __init__(...)

    def _SetFxnDesc(self):
        
        self.fxnDesc['DisplayName'] = 'Dimension Value of Min Over Dimension'
        self.fxnDesc['ShortDesc'] = 'Dimension value of minimum of an NCVar over a dimension  (e.g. which layer had the min).'
        
        self.fxnDesc['ReturnType'] = 'NCVar, Float, or Integer'
        
        self.fxnDesc['ReqArgs'] = {
            'InFieldName':'Field Name',
            'Dimension':'Field Name',
            }
        self.fxnDesc['OptArgs'] = {
            'PrecursorFieldNames':['Field Name','Field Name List'],
            'Metadata':'Tuple List: Field Name:Any',
            'NewNCVarName':'Field Name',
            }
        
    # def _SetFxnDesc(self):

    def DependencyNms(self):
        return self._ArgToList('InFieldName') + \
          self._ArgToList('PrecursorFieldNames')

    def Exec(self,executedObjects):

        self._Exec(executedObjects,'dimension_val_of_min')
        
    # def Exec(self,executedObjects):

#class DimensionValueOfMinOverDimension(_SummaryOverDimensions):



    

################################################################################
# Stats between arrays
################################################################################

class R2(_NetCDFUtilParent):

    def __init__(
        self,
        mpt_cmd_struct=None,
        ):

        super(R2,self).__init__(mpt_cmd_struct)
        
    # def __init__(...)

    def _SetFxnDesc(self):
        # description of command for validation and info display

        self.fxnDesc['DisplayName'] = 'Statistical R2'
        self.fxnDesc['ShortDesc'] = 'Calculates R2 between two arrays of values'
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

        x_arr = executedObjects[self.ValFromArgByNm('XFieldName')].ExecRslt().data.data.compressed()
        y_arr = executedObjects[self.ValFromArgByNm('YFieldName')].ExecRslt().data.data.compressed()

        self.execRslt = stats.linregress(x_arr,y_arr).rvalue ** 2
        executedObjects[self.RsltNm()] = self

    # def Exec(self,executedObjects):

# class R2(_NetCDFUtilParent):

class PValue(_NetCDFUtilParent):

    def __init__(
        self,
        mpt_cmd_struct=None,
        ):

        super(PValue,self).__init__(mpt_cmd_struct)
        
    # def __init__(...)

    def _SetFxnDesc(self):
        # description of command for validation and info display

        self.fxnDesc['DisplayName'] = 'Statistical P Value'
        self.fxnDesc['ShortDesc'] = 'Calculates P Value for two arrays of values'
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

        x_arr = executedObjects[self.ValFromArgByNm('XFieldName')].ExecRslt().data.data.compressed()
        y_arr = executedObjects[self.ValFromArgByNm('YFieldName')].ExecRslt().data.data.compressed()

        self.execRslt = stats.linregress(x_arr,y_arr).pvalue
        executedObjects[self.RsltNm()] = self

    # def Exec(self,executedObjects):

# class PValue():
    
class StdErr(_NetCDFUtilParent):

    def __init__(
        self,
        mpt_cmd_struct=None,
        ):

        super(StdErr,self).__init__(mpt_cmd_struct)
        
    # def __init__(...)


    def _SetFxnDesc(self):
        # description of command for validation and info display

        self.fxnDesc['DisplayName'] = 'Statistical Standard Error'
        self.fxnDesc['ShortDesc'] = 'Calculates Standard Error for two arrays of values'
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

        x_arr = executedObjects[self.ValFromArgByNm('XFieldName')].ExecRslt().data.data.compressed()
        y_arr = executedObjects[self.ValFromArgByNm('YFieldName')].ExecRslt().data.data.compressed()

        self.execRslt = stats.linregress(x_arr,y_arr).stderr
        executedObjects[self.RsltNm()] = self

    # def Exec(self,executedObjects):
    
# class StdErr(_NetCDFUtilParent):

class MeanSquaredError(_NetCDFUtilParent):

    def __init__(
        self,
        mpt_cmd_struct=None,
        ):

        super(MeanSquaredError,self).__init__(mpt_cmd_struct)
        
    # def __init__(...)


    def _SetFxnDesc(self):
        # description of command for validation and info display

        self.fxnDesc['DisplayName'] = 'Statistical Mean Squared Error'
        self.fxnDesc['ShortDesc'] = 'Calculates Mean Squared Error for two arrays of values'
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

        x_arr = executedObjects[self.ValFromArgByNm('XFieldName')].ExecRslt().data.data.compressed()
        y_arr = executedObjects[self.ValFromArgByNm('YFieldName')].ExecRslt().data.data.compressed()

        self.execRslt = ((x_arr - y_arr) ** 2).sum() / float(x_arr.shape[0])
        executedObjects[self.RsltNm()] = self

    # def Exec(self,executedObjects):
    
# class MeanSquaredError(_NetCDFUtilParent):

################################################################################
# Operations involving multiple arrays
################################################################################

class ZScore(_SummaryOverDimensions):

    def __init__(
        self,
        mpt_cmd_struct=None,
        ):

        super(ZScore,self).__init__(mpt_cmd_struct)
        
    # def __init__(...)

    def _SetFxnDesc(self):
        
        self.fxnDesc['DisplayName'] = 'ZScore'
        self.fxnDesc['ShortDesc'] = 'ZScore of array. (Value - Mean) / Std Dev. Mean and Std Dev arrays must have matching dimensions and be projectable to Value array.'
        
        self.fxnDesc['ReturnType'] = 'NCVar'
        
        self.fxnDesc['ReqArgs'] = {
            'ValueFieldName':'Field Name',
            'MeanFieldName':'Field Name',
            'StdDevFieldName':'Field Name',
            }
        
        self.fxnDesc['OptArgs'] = {
            'Metadata':'Tuple List: Field Name:Any',
            'NewNCVarName':'Field Name',
            'PrecursorFieldNames':['Field Name','Field Name List'],
            }
        
    # def _SetFxnDesc(self):

    def DependencyNms(self):
        return self._ArgToList('ValueFieldName') + \
            self._ArgToList('MeanFieldName') + \
            self._ArgToList('StdDevFieldName') + \
            self._ArgToList('PrecursorFieldNames')

    def Exec(self,executedObjects):
        
        std_ncdimvar = executedObjects[self.ValFromArgByNm('StdDevFieldName')].ExecRslt()
        std_mskdarr = std_ncdimvar.data.data
        std_dim_nms = std_ncdimvar.dims.keys()
        
        mean_ncdimvar = executedObjects[self.ValFromArgByNm('MeanFieldName')].ExecRslt()
        mean_mskdarr = mean_ncdimvar.data.data
        mean_dim_nms = mean_ncdimvar.dims.keys()
        
        value_ncdimvar = executedObjects[self.ValFromArgByNm('ValueFieldName')].ExecRslt()
        value_mskdarr = value_ncdimvar.data.data
        value_dim_nms = value_ncdimvar.dims.keys()

        # StdDev and Mean arrays must be congruent
        if not ac.array_is_congruent_to_tgt(
            std_mskdarr,
            std_dim_nms,
            mean_mskdarr,
            mean_dim_nms
            ):

            raise Exception(
                (3*'{}\n').format(
                    'Error: Std Dev and Mean arrays must be congruent.',
                    '  Std Dev NCVar name, MPilot variable name: {}  {}'.format(
                        std_ncdimvar.name,
                        self.ArgByNm('StdDevFieldName')
                        ),
                    '  Mean NCVar name, MPilot variable name: {}  {}'.format(
                        mean_ncdimvar.name,
                        self.ArgByNm('MeanFieldName')
                        ),
                    )
                ) # raise Exception(...)

        # if not ac.array_is_congruent_to_tgt(

        # If need be, expand or transpose the std dev and man arrays to match the value array
        if ac.array_is_congruent_to_tgt(
            std_mskdarr,
            std_dim_nms,
            value_mskdarr,
            value_dim_nms
            ):

            # We don't need to change the std dev and mean arrays
            pass

        elif ac.can_expand_to_tgt(
            std_mskdarr,
            std_dim_nms,
            value_mskdarr,
            value_dim_nms
            ):

            std_mskdarr = ac.expand_arr_to_match_tgt(
                std_mskdarr,
                std_dim_nms,
                value_mskdarr,
                value_dim_nms
                )

            mean_mskdarr = ac.expand_arr_to_match_tgt(
                mean_mskdarr,
                mean_dim_nms,
                value_mskdarr,
                value_dim_nms
                )

        elif ac.can_transpose_to_tgt(
                std_mskdarr,
                std_dim_nms,
                value_mskdarr,
                value_dim_nms
            ):

            std_mskdarr = ac.transpose_arr_to_match_tgt(
                std_mskdarr,
                std_dim_nms,
                value_mskdarr,
                value_dim_nms
                )

            mean_mskdarr = ac.transpose_arr_to_match_tgt(
                mean_mskdarr,
                mean_dim_nms,
                value_mskdarr,
                value_dim_nms
                )

        else:
            
            raise Exception(
                (4*'{}\n').format(
                    'Std Dev and Mean arrays cannot be made to match Value array shape.',
                    '  ValueFieldName, NCVar name: {}  {}'.format(
                        self.ValFromArgByNm('ValueFieldName'),
                        value_ncdimvar.name
                        ),
                    '  StdDevFieldName, NCVar name:   {}  {}'.format(
                        self.ValFromArgByNm('StdDevFieldName'),
                        std_ncdimvar.name
                        ),
                    '  MeanFieldName, NCVar name:   {}  {}'.format(
                        self.ValFromArgByNm('MeanFieldName'),
                        mean_ncdimvar.name
                        ),
                    )
                ) # raise Exception(...)


        # if ac.array_is_congruent_to_tgt(...)...elif...elif...else

        self.execRslt = mpncv.NCDimensionedVar(
            dims = cp.deepcopy(value_ncdimvar.dims),
            data = mpncv.NCVar(
                data = (value_mskdarr - mean_mskdarr) / std_mskdarr,
                dim_nms = value_ncdimvar.dims.keys()
                ),
            name = self._CreateNCVarName()
            )

        self._SetMetadata()
        self.execRslt.set_metadata_val('NodeCommand',self.mptCmdStruct['rawCmdStr'])
        executedObjects[self.RsltNm()] = self
        
    # def Exec(self,executedObjects):

# class ZScore(_SummaryOverDimensions):

################################################################################
# Scalar arithmetic
################################################################################

class _ApplyScalarOperation(_NetCDFUtilParent):
    
    def __init__(
        self,
        mpt_cmd_struct=None,
        ):

        super(_ApplyScalarOperation,self).__init__(mpt_cmd_struct)
        
    # def __init__(...)
    
    def _Exec(self,executedObjects,operation):

        new_ncdimvar = cp.deepcopy(executedObjects[self.ValFromArgByNm('InFieldName')].ExecRslt())
        new_ncdimvar.name = self._CreateNCVarName()

        if operation == 'mult':
            new_ncdimvar.data.data = new_ncdimvar.data.data * self.ValFromArgByNm('Scalar')
        elif operation == 'div':
            new_ncdimvar.data.data = new_ncdimvar.data.data / self.ValFromArgByNm('Scalar')
        elif operation == 'add':
            new_ncdimvar.data.data = new_ncdimvar.data.data + self.ValFromArgByNm('Scalar')
        elif operation == 'sub':
            new_ncdimvar.data.data = new_ncdimvar.data.data - self.ValFromArgByNm('Scalar')
        else:
            raise Exception(
                (2*'{}\n').format(
                    'Error: Illegal summary operation.',
                    '  Summary operation: {}'.format(summary_type)
                    )
                ) # raise Exception(...)
            
        # Get name of new NCDimensionedVar

        self.execRslt = new_ncdimvar
        self._SetMetadata()
        self.execRslt.set_metadata_val('NodeCommand',self.mptCmdStruct['rawCmdStr'])
        executedObjects[self.RsltNm()] = self
        
    # def _Exec(self,executedObjects,operation):
        
# class _ApplyScalarOperation(_NetCDFUtilParent):

class MultiplyByScalar(_ApplyScalarOperation):
    
    def __init__(
        self,
        mpt_cmd_struct=None,
        ):

        super(MultiplyByScalar,self).__init__(mpt_cmd_struct)
        
    # def __init__(...)

    def _SetFxnDesc(self):
        
        self.fxnDesc['DisplayName'] = 'MultiplyByScalar'
        self.fxnDesc['ShortDesc'] = 'Multipy NCVar data array by a scalar'
        
        self.fxnDesc['ReturnType'] = 'NCVar'
        
        self.fxnDesc['ReqArgs'] = {
            'InFieldName':'Field Name',
            'Scalar':'Float'
            }
        self.fxnDesc['OptArgs'] = {
            'PrecursorFieldNames':['Field Name','Field Name List'],
            'Metadata':'Tuple List: Field Name:Any',
            'NewNCVarName':'Field Name',
            }
        
    # def _SetFxnDesc(self):

    def DependencyNms(self):
        return self._ArgToList('InFieldName') + \
          self._ArgToList('PrecursorFieldNames')

    def Exec(self,executedObjects):

        self._Exec(executedObjects,'mult')
        
    # def Exec(self,executedObjects):

# class MultiplyByScalar(_ApplyScalarOperation):

class DivideByScalar(_ApplyScalarOperation):
    
    def __init__(
        self,
        mpt_cmd_struct=None,
        ):

        super(DivideByScalar,self).__init__(mpt_cmd_struct)
        
    # def __init__(...)

    def _SetFxnDesc(self):
        
        self.fxnDesc['DisplayName'] = 'DivideByScalar'
        self.fxnDesc['ShortDesc'] = 'Divide NCVar data array by a scalar'
        
        self.fxnDesc['ReturnType'] = 'NCVar'
        
        self.fxnDesc['ReqArgs'] = {
            'InFieldName':'Field Name',
            'Scalar':'Float'
            }
        self.fxnDesc['OptArgs'] = {
            'PrecursorFieldNames':['Field Name','Field Name List'],
            'Metadata':'Tuple List: Field Name:Any',
            'NewNCVarName':'Field Name',
            }
        
    # def _SetFxnDesc(self):

    def DependencyNms(self):
        return self._ArgToList('InFieldName') + \
          self._ArgToList('PrecursorFieldNames')

    def Exec(self,executedObjects):

        self._Exec(executedObjects,'div')
        
    # def Exec(self,executedObjects):

# class DivideByScalar(_ApplyScalarOperation):

class AddScalar(_ApplyScalarOperation):
    
    def __init__(
        self,
        mpt_cmd_struct=None,
        ):

        super(AddScalar,self).__init__(mpt_cmd_struct)
        
    # def __init__(...)

    def _SetFxnDesc(self):
        
        self.fxnDesc['DisplayName'] = 'AddScalar'
        self.fxnDesc['ShortDesc'] = 'Add a scalar to NCVar data array'
        
        self.fxnDesc['ReturnType'] = 'NCVar'
        
        self.fxnDesc['ReqArgs'] = {
            'InFieldName':'Field Name',
            'Scalar':'Float'
            }
        self.fxnDesc['OptArgs'] = {
            'PrecursorFieldNames':['Field Name','Field Name List'],
            'Metadata':'Tuple List: Field Name:Any',
            'NewNCVarName':'Field Name',
            }
        
    # def _SetFxnDesc(self):

    def DependencyNms(self):
        return self._ArgToList('InFieldName') + \
          self._ArgToList('PrecursorFieldNames')

    def Exec(self,executedObjects):

        self._Exec(executedObjects,'add')
        
    # def Exec(self,executedObjects):

# class AddScalar(_ApplyScalarOperation):

class SubtractScalar(_ApplyScalarOperation):
    
    def __init__(
        self,
        mpt_cmd_struct=None,
        ):

        super(SubtractScalar,self).__init__(mpt_cmd_struct)
        
    # def __init__(...)

    def _SetFxnDesc(self):
        
        self.fxnDesc['DisplayName'] = 'SubtractScalar'
        self.fxnDesc['ShortDesc'] = 'Subtract a scalar from NCVar data array'
        
        self.fxnDesc['ReturnType'] = 'NCVar'
        
        self.fxnDesc['ReqArgs'] = {
            'InFieldName':'Field Name',
            'Scalar':'Float'
            }
        self.fxnDesc['OptArgs'] = {
            'PrecursorFieldNames':['Field Name','Field Name List'],
            'Metadata':'Tuple List: Field Name:Any',
            'NewNCVarName':'Field Name',
            }
        
    # def _SetFxnDesc(self):

    def DependencyNms(self):
        return self._ArgToList('InFieldName') + \
          self._ArgToList('PrecursorFieldNames')

    def Exec(self,executedObjects):

        self._Exec(executedObjects,'sub')
        
    # def Exec(self,executedObjects):

# class SubtractScalar(_ApplyScalarOperation):

class Normalize(_NetCDFUtilParent):
    
    def __init__(
        self,
        mpt_cmd_struct=None,
        ):

        super(Normalize,self).__init__(mpt_cmd_struct)
        
    # def __init__(...)

    def _SetFxnDesc(self):
        
        self.fxnDesc['DisplayName'] = 'Normalize'
        self.fxnDesc['ShortDesc'] = 'Normalizes an NCVar data array to a mean value of 1.'
        
        self.fxnDesc['ReturnType'] = 'NCVar'
        
        self.fxnDesc['ReqArgs'] = {
            'InFieldName':'Field Name',
            }
        self.fxnDesc['OptArgs'] = {
            'PrecursorFieldNames':['Field Name','Field Name List'],
            'Metadata':'Tuple List: Field Name:Any',
            'NewNCVarName':'Field Name',
            }
        
    # def _SetFxnDesc(self):

    def DependencyNms(self):
        return self._ArgToList('InFieldName') + \
          self._ArgToList('PrecursorFieldNames')

    def Exec(self,executedObjects):

        new_ncdimvar = cp.deepcopy(executedObjects[self.ValFromArgByNm('InFieldName')].ExecRslt())
        new_ncdimvar.name = self._CreateNCVarName()
        new_ncdimvar.data.data /= new_ncdimvar.data.data.mean()
        self.execRslt = new_ncdimvar
        self._SetMetadata()
        self.execRslt.set_metadata_val('NodeCommand',self.mptCmdStruct['rawCmdStr'])
        executedObjects[self.RsltNm()] = self
        
    # def Exec(self,executedObjects):

# class SubtractScalar(_ApplyScalarOperation):

################################################################################
# Array arithmetic
################################################################################

class AbsoluteValue(_NetCDFUtilParent):
    
    def __init__(
        self,
        mpt_cmd_struct=None,
        ):

        super(AbsoluteValue,self).__init__(mpt_cmd_struct)
        
    # def __init__(...)

    def _SetFxnDesc(self):
        
        self.fxnDesc['DisplayName'] = 'Absolute Value'
        self.fxnDesc['ShortDesc'] = 'Takes absolute value of array values'
        
        self.fxnDesc['ReturnType'] = 'NCVar'
        
        self.fxnDesc['ReqArgs'] = {
            'InFieldName':'Field Name',
            }
        self.fxnDesc['OptArgs'] = {
            'PrecursorFieldNames':['Field Name','Field Name List'],
            'Metadata':'Tuple List: Field Name:Any',
            'NewNCVarName':'Field Name',
            }
        
    # def _SetFxnDesc(self):

    def DependencyNms(self):
        return self._ArgToList('InFieldName') + \
          self._ArgToList('PrecursorFieldNames')

    def Exec(self,executedObjects):

        new_ncdimvar = cp.deepcopy(executedObjects[self.ValFromArgByNm('InFieldName')].ExecRslt())
        new_ncdimvar.data.data = np.ma.absolute(new_ncdimvar.data.data)
        self.execRslt = new_ncdimvar
        self.execRslt.name = self._CreateNCVarName()        
        self._SetMetadata()
        self.execRslt.set_metadata_val('NodeCommand',self.mptCmdStruct['rawCmdStr'])
        executedObjects[self.RsltNm()] = self
        
    # def Exec(self,executedObjects):

# class AddScalar(_ApplyScalarOperation):

class _DoArrayArithmetic(_NetCDFUtilParent):

    def __init__(
        self,
        mpt_cmd_struct=None,
        ):

        super(_DoArrayArithmetic,self).__init__(mpt_cmd_struct)
        
    # def __init__(...)

    def DependencyNms(self):
        return self._ArgToList('InFieldNames') + \
          self._ArgToList('PrecursorFieldNames')

    def _Exec(self,executedObjects,operation):

        fldnm_list = self._ArgToList('InFieldNames')

        if len(fldnm_list) < 2:
            raise Exception(
                (2*'{}\n').format(
                    'Error. Must have 2 or more Field Names',
                    '  Result Name, Field Name List: {}  [{}]'.format(
                        self.RsltNm(),
                        ','.join(fldnm_list)
                        ),
                    )
                ) # raise Exception(...)

        # if len(fldnm_list) < 2:

        new_ncdimvar = cp.deepcopy(executedObjects[fldnm_list[0]].ExecRslt())
        new_ncdimvar.name = self._CreateNCVarName()
        
        for fldnm in fldnm_list[1:]:
            
            ncdimvar_for_op = executedObjects[fldnm].ExecRslt()

            # If necessary, expand the array being applied to new_ncdimvar

            if self._NCDimArrsAreCongruent(ncdimvar_for_op,new_ncdimvar):

                delete_ncdimvar_for_op = False
                
            else:

                ncdimvar_for_op = self._MakeCongruentNCDimArr(ncdimvar_for_op,new_ncdimvar)

                if ncdimvar_for_op is not None:
                    delete_ncdimvar_for_op = True
                else:
                    raise Exception(
                        (4*'{}\n').format(
                            'Error. Cannot make array congruent to first array in arithmetic array operation.',
                            '  Result name: {}'.format(self.RsltNm()),
                            '  First array: Field Name, NCVar name, dim names, shape: {}  {}  {}  {}'.format(
                                fldnm_list[0],
                                executedObjects[fldnm_list[0]].ExecRslt().name,
                                executedObjects[fldnm_list[0]].ExecRslt().dims.keys(),
                                executedObjects[fldnm_list[0]].ExecRslt().data.data.shape,
                                ),
                            '  Incongruent array: Field Name, NCVar name, dim names, shape: {}  {}  {}  {}'.format(
                                fldnm_list[0],
                                executedObjects[fldnm].ExecRslt().name,
                                executedObjects[fldnm].ExecRslt().dims.keys(),
                                executedObjects[fldnm].ExecRslt().data.data.shape,
                                ),
                            )
                        ) # raise Exception(...)
                        
            # if self._NCDimArrsAreCongruent(ncdimvar_for_op,new_ncdimvar):...else...

            # do the operation
            if operation == 'add':
                new_ncdimvar.data.data = new_ncdimvar.data.data + ncdimvar_for_op.data.data
            elif operation == 'subtract':
                new_ncdimvar.data.data = new_ncdimvar.data.data - ncdimvar_for_op.data.data
            elif operation == 'multiply':
                new_ncdimvar.data.data = new_ncdimvar.data.data * ncdimvar_for_op.data.data
            elif operation == 'divide':
                new_ncdimvar.data.data = new_ncdimvar.data.data / ncdimvar_for_op.data.data
            else:
                raise Exception(
                    (2*'{}\n').format(
                        'Error: Illegal array arithmetic operation.',
                        '  Operation: {}'.format(summary_type)
                        )
                    ) # raise Exception(...)                    

            if delete_ncdimvar_for_op:
                del ncdimvar_for_op
                    
        # for fldnm in fldnm_list[1:]:

        self.execRslt = new_ncdimvar
        self._SetMetadata()
        self.execRslt.set_metadata_val('NodeCommand',self.mptCmdStruct['rawCmdStr'])
        executedObjects[self.RsltNm()] = self
        
    # def _Exec(self,executedObjects,operation):
    
# class _DoArrayArithmetic(_NetCDFUtilParent):
    
class SumArrays(_DoArrayArithmetic):
    
    def __init__(
        self,
        mpt_cmd_struct=None,
        ):

        super(SumArrays,self).__init__(mpt_cmd_struct)

    # def __init__(...)

    def _SetFxnDesc(self):
        
        self.fxnDesc['DisplayName'] = 'SumArrays'
        self.fxnDesc['ShortDesc'] = 'Sum a list of arrays. Default NCVar name and size is that of first array in list.'
        
        self.fxnDesc['ReturnType'] = 'NCVar'
        
        self.fxnDesc['ReqArgs'] = {
            'InFieldNames':'Field Name List',
            }
        self.fxnDesc['OptArgs'] = {
            'PrecursorFieldNames':['Field Name','Field Name List'],
            'Metadata':'Tuple List: Field Name:Any',
            'NewNCVarName':'Field Name',
            }
        
    # def _SetFxnDesc(self):

    def Exec(self,executedObjects):

        self._Exec(executedObjects,'add')
        
    # def _Exec(self,executedObjects):
     
# class SumArrays(_DoArrayArithmetic):

class SubtractArrays(_DoArrayArithmetic):
    
    def __init__(
        self,
        mpt_cmd_struct=None,
        ):

        super(SubtractArrays,self).__init__(mpt_cmd_struct)

    # def __init__(...)

    def _SetFxnDesc(self):
        
        self.fxnDesc['DisplayName'] = 'SubtractArrays'
        self.fxnDesc['ShortDesc'] = 'Subtract arrays from the first array in the list. Default NCVar name and size is that of first array in list.'
        
        self.fxnDesc['ReturnType'] = 'NCVar'
        
        self.fxnDesc['ReqArgs'] = {
            'InFieldNames':'Field Name List',
            }
        self.fxnDesc['OptArgs'] = {
            'PrecursorFieldNames':['Field Name','Field Name List'],
            'Metadata':'Tuple List: Field Name:Any',
            'NewNCVarName':'Field Name',
            }
        
    # def _SetFxnDesc(self):

    def Exec(self,executedObjects):

        self._Exec(executedObjects,'subtract')
        
    # def _Exec(self,executedObjects):
     
# class SubtractArrays(_DoArrayArithmetic):

class MultiplyArrays(_DoArrayArithmetic):
    
    def __init__(
        self,
        mpt_cmd_struct=None,
        ):

        super(MultiplyArrays,self).__init__(mpt_cmd_struct)

    # def __init__(...)

    def _SetFxnDesc(self):
        
        self.fxnDesc['DisplayName'] = 'MultiplyArrays'
        self.fxnDesc['ShortDesc'] = 'Multiply arrays in a list of arrays list. Default NCVar name and size is that of first array in list.'        
        self.fxnDesc['ReturnType'] = 'NCVar'
        
        self.fxnDesc['ReqArgs'] = {
            'InFieldNames':'Field Name List',
            }
        self.fxnDesc['OptArgs'] = {
            'PrecursorFieldNames':['Field Name','Field Name List'],
            'Metadata':'Tuple List: Field Name:Any',
            'NewNCVarName':'Field Name',
            }
        
    # def _SetFxnDesc(self):

    def Exec(self,executedObjects):

        self._Exec(executedObjects,'multiply')
        
    # def _Exec(self,executedObjects):

# class MultiplyArrays(_DoArrayArithmetic):

class DivideArrays(_DoArrayArithmetic):
    
    def __init__(
        self,
        mpt_cmd_struct=None,
        ):

        super(DivideArrays,self).__init__(mpt_cmd_struct)

    # def __init__(...)

    def _SetFxnDesc(self):
        
        self.fxnDesc['DisplayName'] = 'MultiplyArrays'
        self.fxnDesc['ShortDesc'] = 'Divide the first array in the list by each subsequent array. Default NCVar name and size is that of first array in list.'
        
        self.fxnDesc['ReturnType'] = 'NCVar'
        
        self.fxnDesc['ReqArgs'] = {
            'InFieldNames':'Field Name List',
            }
        self.fxnDesc['OptArgs'] = {
            'PrecursorFieldNames':['Field Name','Field Name List'],
            'Metadata':'Tuple List: Field Name:Any',
            'NewNCVarName':'Field Name',
            }
        
    # def _SetFxnDesc(self):

    def Exec(self,executedObjects):

        self._Exec(executedObjects,'divide')
        
    # def _Exec(self,executedObjects):
# class DivideArrays(_DoArrayArithmetic):

class ConcatenateArrays(_NetCDFUtilParent):
    
    def __init__(
        self,
        mpt_cmd_struct=None,
        ):

        super(ConcatenateArrays,self).__init__(mpt_cmd_struct)

    # def __init__(...)

    def _SetFxnDesc(self):
        
        self.fxnDesc['DisplayName'] = 'ConcatenateArrays'
        self.fxnDesc['ShortDesc'] = 'Concatenates in order listed.'
        
        self.fxnDesc['ReturnType'] = 'NCVar'
        
        self.fxnDesc['ReqArgs'] = {
            'InFieldNames':'Field Name List',
            'DimensionName':'Field Name'
            }
        self.fxnDesc['OptArgs'] = {
            'PrecursorFieldNames':['Field Name','Field Name List'],
            'Metadata':'Tuple List: Field Name:Any',
            'NewNCVarName':'Field Name',
            }
        
    # def _SetFxnDesc(self):

    def DependencyNms(self):
        return self._ArgToList('InFieldNames') + \
          self._ArgToList('PrecursorFieldNames')

    def Exec(self,executedObjects):

        in_fld_nms = self._ArgToList('InFieldNames')
        cat_dim_nm = self.ValFromArgByNm('DimensionName')

        new_ncdimvar = cp.deepcopy(executedObjects[in_fld_nms[0]].ExecRslt())
        new_ncdimvar_fldnm = in_fld_nms[0]
        cat_ncdimvar_fldms = in_fld_nms[1:]
        cat_ncdimvars = []

        for ndx in range(1,len(in_fld_nms)):
            cat_ncdimvars.append(executedObjects[in_fld_nms[ndx]].ExecRslt())

        # Check for cat dimension presence
        if self.ValFromArgByNm('DimensionName') not in new_ncdimvar.dim_nms:
            raise Exception(
                (1*'{}\n').format(
                    'Error: Invalid Dimension Name:  {}'.format(cat_dim_nm),
                    )
                ) # raise Exception(...)      

        # Check that every ncdimvar has the right dimensions and sizes
        new_ncdimvar_dim_szs = [len(new_ncdimvar.dims[nm].data) for nm in new_ncdimvar.dim_nms]

        for cat_ncdimvar,cat_ncdimvar_fldnm in zip(cat_ncdimvars,cat_ncdimvar_fldms):

            if cat_ncdimvar.dim_nms != new_ncdimvar.dim_nms:
                raise Exception(
                    (3*'{}\n').format(
                        'Error: Dimension names of fields do not match.',
                        '  Field Name, NCVar Name, Dimension Names: {}  {}  {}'.format(
                            new_ncdimvar_fldnm,
                            new_ncdimvar.name,
                            new_ncdimvar.dim_nms
                            ),
                        '  Field Name, NCVar Name, Dimension Names: {}  {}  {}'.format(
                            cat_ncdimvar_fldnm,
                            cat_ncdimvar.name,
                            cat_ncdimvar.dim_nms
                            )
                        )
                    )
                # raise Exception(

            new_ncdimvar_dim_szs = [len(new_ncdimvar.dims[nm].data) for nm in new_ncdimvar.dim_nms]
            
            for dim_nm,dim_sz in zip(new_ncdimvar.dim_nms,new_ncdimvar_dim_szs):
                
                if dim_nm == cat_dim_nm: continue

                if len(cat_ncdimvar.dims[dim_nm].data) != dim_sz:
                    raise Exception(
                        (3*'{}\n').format(
                            'Error: Dimension sizes of non-concatenated dimension do not match.',
                            '  Field Name, NCVar Name, Dimension Name, Dimension Size: {}  {}  {}  {}'.format(
                                new_ncdimvar_fldnm,
                                new_ncdimvar.name,
                                dim_nm,
                                dim_sz
                                ),
                            '  Field Name, NCVar Name, Dimension Name, Dimension Size: {}  {}  {}  {}'.format(
                                cat_ncdimvar_fldnm,
                                cat_ncdimvar.name,
                                dim_nm,
                                len(cat_ncdimvar.dims[dim_nm].data)
                                )
                            )
                        )
                    # raise Exception(

            # for dim_nm,dim_sz in zip(new_ncdimvar.dim_nms,new_ncdimvar_dim_szs)

        # for cat_ncdimvar,cat_ncdimvar_fldnm in zip(cat_ncdimvars,cat_ncdimvar_fldnms):
        
        # Now do the concatenating

        cat_axis = new_ncdimvar.dim_nms.index(cat_dim_nm)

        for cat_ncdimvar in cat_ncdimvars:

            # cat the data
            new_ncdimvar.data.data = np.ma.concatenate(
                [new_ncdimvar.data.data,cat_ncdimvar.data.data],
                axis = cat_axis
                )

            # cat the axis dimension
            new_ncdimvar.dims[cat_dim_nm].data = np.ma.concatenate(
                [new_ncdimvar.dims[cat_dim_nm].data,cat_ncdimvar.dims[cat_dim_nm].data]
                )

        new_ncdimvar.name = self._CreateNCVarName()

        self.execRslt = new_ncdimvar        
        self._SetMetadata()
        self.execRslt.set_metadata_val('NodeCommand',self.mptCmdStruct['rawCmdStr'])
        executedObjects[self.RsltNm()] = self
        
    # def _Exec(self,executedObjects):
# class ConcatenateArrays():

################################################################################
# Class Composition
################################################################################

class ClassFractionComposition(_NetCDFUtilParent):

    def __init__(
        self,
        mpt_cmd_struct=None,
        ):

        super(ClassFractionComposition,self).__init__(mpt_cmd_struct)

    # def __init__(...)

    def _SetFxnDesc(self):
        
        self.fxnDesc['DisplayName'] = 'Class Fraction Composition'
        self.fxnDesc['ShortDesc'] = 'Calculates fraction of each class in array, with optional weighting.'
        self.fxnDesc['ReturnType'] = 'NCVar'
        
        self.fxnDesc['ReqArgs'] = {
            'InFieldName':'Field Name'
            }
        self.fxnDesc['OptArgs'] = {
            'PrintResults':'Boolean',
            'PrecursorFieldNames':['Field Name','Field Name List'],
            'Metadata':'Tuple List: Field Name:Any',
            'NewNCVarName':'Field Name'
            }
        
    # def _SetFxnDesc(self):

    def DependencyNms(self):
        
        rtrn = self._ArgToList('InFieldName')
        if self.ArgByNm('WeightFieldName') is not None:
            rtrn += self._ArgToList('WeightFieldName')
        rtrn += self._ArgToList('PrecursorFieldNames')
        return rtrn
    
    # def DependencyNms(self):

    def Exec(self,executedObjects):

        class_ncdimvar = executedObjects[self.ValFromArgByNm('InFieldName')].ExecRslt()
        class_vals, class_cnts = np.unique(class_ncdimvar.data.data.flatten().compressed(),return_counts = True)
        class_fracs = class_cnts.astype(np.float) / class_cnts.sum()

        # Create the new NCDimensionedVar
        new_ncdimvar = mpncv.NCDimensionedVar(name=self._CreateNCVarName())

        # The data NCVar:
        new_ncdimvar.data = mpncv.NCVar(
            data = class_fracs,
            dim_nms = ['ClassValuesFor{}'.format(new_ncdimvar.name)],
            name = new_ncdimvar.name
            )

        # The NCVars for the indexes:
        new_dims = OrderedDict()
        new_dims['ClassValuesFor{}'.format(new_ncdimvar.name)] = mpncv.NCVar(
            data = class_vals,
            dim_nms = ['ClassValuesFor{}'.format(new_ncdimvar.name)],
            name = 'ClassValuesFor{}'.format(new_ncdimvar.name)
            )

        new_ncdimvar.dims = new_dims

        self.execRslt = new_ncdimvar
        self._SetMetadata()
        self.execRslt.set_metadata_val('NodeCommand',self.mptCmdStruct['rawCmdStr'])
        self.execRslt.set_metadata_val('Description','Class values for new_ncdimvar')
        executedObjects[self.RsltNm()] = self

        if self.ValFromArgByNm('PrintResults') == True:

            print 'Class composition for {} (NCVariableName {}):'.format(self.RsltNm(),new_ncdimvar.name)
            print 'Class Val, Class Fraction'
            for class_val,class_frac in zip(class_vals,class_fracs):
                print '{}, {}'.format(class_val,class_frac)
        
    # def Exec(self)

# class ClassFractionComposition(_NetCDFUtilParent):
    
################################################################################
# Graphics
################################################################################

class _GraphicsParent(_NetCDFUtilParent):

    def _InsureDirExists(self,dir_nm):
        
        if os.path.isdir(dir_nm):
            # Dir exists and it is a dir. We're good to go
            pass
        elif os.path.exists(dir_nm):
            # A file with dir_nm exists. No can do
            raise Exception(
                (3*'{}\n').format(
                    'Error: Cannot create a directory because a file with that name exists.',
                    '  Result name: {}'.format(self.RsltNm),
                    '  Directory name: {}'.format(dir_nm)
                    )
                ) # raise Exception(...)                    

        else:
            os.makedirs(dir_nm)

        return True
        
    # def _InsureOutDirExists(self,dir_nm):

    def _RenderLayer(
        self,
        data_obj,
        min_val,
        max_val,
        in_fldnm,
        out_fnm,
        ):

        # Here is the numpy array to render
        tmp_arr = data_obj.ExecRslt().data.data

        if tmp_arr.ndim != 2:
            print (3*'{}\n').format(
                'Warning: Trying to render non 2-D layer.',
                '  Variable name: {}'.format(self.RsltNm()),
                '  Skipping...'
                )
            
            return False

        # Flip axes if needed
        if self.ValFromArgByNm('FlipX'):
            tmp_arr = np.flipud(tmp_arr)
            
        if self.ValFromArgByNm('FlipY'):
            tmp_arr = np.fliplr(tmp_arr)
        
        # Build the graphic
        # Note x is associated with second dimension and y with first
        x_size = data_obj.ExecRslt().shape[1]/100.
        y_size = data_obj.ExecRslt().shape[0]/100.

        if not self.ValFromArgByNm('PlainImageOnly'):
            # Additions are for stuff around the edges
            x_size += 2
            y_size += 1.5
            
        fig = plt.figure(
            figsize = (
                x_size,
                y_size
                ),
            dpi = 600
            )

        face_color = self.ValFromArgByNm('FaceColor')
        if face_color is None:
            face_color = 'grey'
        else:
            rgb_r, rgb_g, rgb_b = [float(x) for x in face_color.split(':')]
            face_color = (rgb_r/255., rgb_g/255., rgb_b/255.)
            
        ax1 = fig.add_subplot(111,facecolor=face_color)
        if self.ValFromArgByNm('PlainImageOnly'):
            ax1.axis('off')

        ax1.tick_params(
                axis='both',
                which='both',
                bottom=False,
                top=False,
                left=False,
                right=False,
                labelbottom=False,
                labeltop=False,
                labelleft=False,
                labelright=False,
                )

        origin = self.ArgByNm('Origin') if self.ArgByNm('Origin') is not None else 'lower'

        my_img = ax1.imshow(
            tmp_arr,
            aspect='auto',
            interpolation='none',
            origin=origin
            )
        
        my_img.set_clim(min_val,max_val)

        c_map = self.ValFromArgByNm('ColorMap')
        if c_map is None:
            c_map = 'RdYlBu'

        my_img.set_cmap(c_map)

        if self.ValFromArgByNm('PlainImageOnly'):
            pass
        else:
            title_txt = self.ValFromArgByNm('Title')
            if title_txt is None:
                title_txt = in_fldnm
            else:
                title_txt = title_txt.replace('_',' ')

            title_size = self.ValFromArgByNm('TitleSize')
            if title_size is None:
                title_size = 10

            title_offset = self.ValFromArgByNm('TitleOffset')
            if title_offset is None:
                title_offset = 1

            ax1.set_title(
                title_txt,
                fontsize=title_size,
                fontweight='bold',
                y = title_offset
                )

            tick_label_size = self.ValFromArgByNm('TickLabelSize')
            if tick_label_size is None:
                tick_label_size = 10

            cbar = fig.colorbar(my_img)

            tick_width = self.ValFromArgByNm('TickWidth')
            if tick_width is None:
                tick_width = 1

            tick_length = self.ValFromArgByNm('TickLength')
            if tick_length is None:
                tick_length = 1

            cbar.ax.tick_params(
                labelsize=tick_label_size,
                width=tick_width,
                length=tick_length,
                )
            
        # if self.ValFromArgByNm('PlainImageOnly'):...else...
        
        if out_fnm is None:
            out_fnm = tf.mktemp(suffix='.png')
            plt.savefig(out_fnm)
            os.system('open -a Preview {}'.format(out_fnm))

        else:
            #print 'saving',out_fnm
            plt.savefig(out_fnm)
        
        plt.close(fig) # fig must be closed to get it out of memory

        return True

    # def _RenderLayer(self,executedObjects,in_fldnm,out_fnm):

    def _XYLineGraph(
            self,
            x_vals,
            y_vals,
            out_fnm,
            title = None,
            x_lbl = None,
            y_lbl = None,
            x_min = None,
            x_max = None,
            y_min = None,
            y_max = None,
            color = None,
            line_style = None,            
            ):

        if len(x_vals) != len(y_vals):
            raise Exception(
                (4*'{}\n').format(
                    'Error: Number of x values and y values are not the same.',
                    '  Result name: {}'.format(self.RsltNm),
                    '  Number of x vals, y vals: {}, {}'.format(len(x_vals),len(y_vals)),
                    '  Out file name: {}'.format(out_fnm)
                    )
                ) # raise Exception(...)

        
         # Create the graph
        fig = plt.figure()
        ax1 = fig.add_subplot(111,facecolor='white')
        
        if x_lbl is not None: ax1.set_xlabel(x_lbl.replace('_',' '))
        if y_lbl is not None: ax1.set_ylabel(y_lbl.replace('_',' '))
        if title is not None: ax1.set_title(title.replace('_',' '))
        if x_min is None: x_min = min(x_vals)
        if x_max is None: x_max = max(x_vals)
        if y_min is None: y_min = min(y_vals)
        if y_max is None: y_max = max(y_vals)
        ax1.set_xlim(x_min,x_max)
        ax1.set_ylim(y_min,y_max)
        if color is None: color = 'black'
        if line_style is None: line_style = '-'

        ax1.plot(
            x_vals,
            y_vals,
            color=color,
            linestyle = line_style
            )
            
        # for dataKey,dataObj in dataObjs.items():
        
        if out_fnm is None:
            out_fnm = tf.mktemp(suffix='.png')
            plt.savefig(out_fnm)
            os.system('open -a Preview {}'.format(out_fnm))

        else:
            # print 'saving',out_fnm
            plt.savefig(out_fnm)
          
        plt.close(fig) # fig must be closed to get it out of memory

        return True

    # def _XYLineGraph(...)
    
    def _XYLinesOnOneGraph(
            self,
            x_arrs,
            y_arrs,
            title,
            x_lbl,
            y_lbl,
            line_lbls,
            out_fnm,
            fill_between = False
            ):

        # All on one graph

        if not isinstance(y_arrs,list): y_arrs = [y_arrs]
        if not isinstance(x_arrs,list): x_arrs = [x_arrs]
        
        # Right number of arrays and things?

        if len(x_arrs) != len(y_arrs) and len(x_arrs) != 1:
            raise Exception(
                (3*'{}\n').format(
                    'Error: Number of x arrays must be 1 or the same number as y arrays.',
                    '  Result name: {}'.format(self.RsltNm),
                    '  Number of x arrays, y arrays: {}, {}'.format(len(x_arrs),len(y_arrs))
                    )
                ) # raise Exception(...)                    

        if line_lbls is not None:
            if len(y_arrs) != len(line_lbls):
                raise Exception(
                    (3*'{}\n').format(
                        'Error: Number of y arrays must be the same as number y line labels.',
                        '  Result name: {}'.format(self.RsltNm),
                        '  Number of x arrays, y line labels: {}'.format(len(x_arrs),len(line_lbls))
                        )
                    ) # raise Exception(...)



        fig = plt.figure()
        ax1 = fig.add_subplot(111)
        
        # Determine mins and maxes
        x_min = self.ValFromArgByNm('XMin')
        if x_min is None:
            x_min = float('inf')
            for x_arr in x_arrs:
                if isinstance(x_arr,mpncv.NCDimensionedVar):
                    x_min = min(x_min,x_arr.data.data.min())
                elif isinstance(x_arr,mpncv.NCVar):
                    x_min = min(x_min,x_arr.data.min())
                else:
                    x_min = min(x_min,x_arr.min())
                    
        x_max = self.ValFromArgByNm('XMax')
        if x_max is None:
            x_max = float('-inf')
            for x_arr in x_arrs:
                if isinstance(x_arr,mpncv.NCDimensionedVar):
                    x_max = max(x_max,x_arr.data.data.max())
                elif isinstance(x_arr,mpncv.NCVar):
                    x_max = max(x_max,x_arr.data.max())
                else:
                    x_max = max(x_max,x_arr.max())
                    
        y_min = self.ValFromArgByNm('YMin')
        if y_min is None:
            y_min = float('inf')
            for y_arr in y_arrs:
                if isinstance(y_arr,mpncv.NCDimensionedVar):
                    y_min = min(y_min,y_arr.data.data.min())
                elif isinstance(y_arr,mpncv.NCVar):
                    y_min = min(y_min,y_arr.data.min())
                else:
                    y_min = min(y_min,y_arr.min())
                    
        y_max = self.ValFromArgByNm('YMax')
        if y_max is None:
            y_max = float('-inf')
            for y_arr in y_arrs:
                if isinstance(y_arr,mpncv.NCDimensionedVar):
                    y_max = max(y_max,y_arr.data.data.max())
                elif isinstance(y_arr,mpncv.NCVar):
                    y_max = max(y_max,y_arr.data.max())
                else:
                    y_max = max(y_max,y_arr.max())

        # Arrange the colors for the lines
        if self.ValFromArgByNm('Colors') is not None:
            
            color_list = self._ArgToList('Colors')
            
            if len(color_list) != len(x_arrs):
                raise Exception(
                    (3*'{}\n').format(
                        'Error: Colors must match x field names 1 to 1.',
                        '  Result name: {}'.format(self.RsltNm),
                        '  Color list: {}'.format(color_list),
                        )
                    ) # raise Exception(...)

            # Convert to rgb if colors are specified that way
            for color_ndx in range(len(color_list)):
                if re.match('^rgb[0-9a-fA-F]{6}$',color_list[color_ndx]):
                    color_list[color_ndx] = color_list[color_ndx].replace('rgb','#')
            # for color_ndx in len(color_list):
            
        else:

            color_list = ['black'] * len(x_arrs)
            
        # if self.ValFromArgByNm('Colors') is not None:

        # Arrange the linestyles for the lines
        if self.ValFromArgByNm('LineStyles') is not None:
            
            linestyle_list = self._ArgToList('LineStyles')
            
            if len(linestyle_list) != len(x_arrs):
                raise Exception(
                    (3*'{}\n').format(
                        'Error: LineStyles must match x field names 1 to 1.',
                        '  Result name: {}'.format(self.RsltNm),
                        '  LineStyle list: {}'.format(linestyle_list),
                        )
                    ) # raise Exception(...)
        else:

            linestyle_list = ['-'] * len(x_arrs)
            
        # if self.ValFromArgByNm('LineStyles') is not None:

        if self.ValFromArgByNm('LineTitles') is not None:
            
            line_lbls = self._ArgToList('LineTitles')
            if len(line_lbls) != len(x_arrs):
                raise Exception(
                    (4*'{}\n').format(
                        'Error: LineTitles must match x field names 1 to 1.',
                        '  Result name: {}'.format(self.RsltNm),
                        '  Line titles: {}'.format(line_lbls),
                        )
                    ) # raise Exception(...)

            # Replace underscores with space
            line_lbls = [s.replace('_',' ') for s in line_lbls]

        # if self.ValFromArgByNm('LineTitles') is not None:
        
        # Create the graph    
        ax1.set_xlabel(x_lbl)
        ax1.set_ylabel(y_lbl)
        ax1.set_title(title.replace('_',' '))
        ax1.set_xlim(x_min,x_max)
        ax1.set_ylim(y_min,y_max)

        # Now loop through the data,
        # and plot the lines

        # create these before using them so we can fill area if needed
        x_plot_arrs = []
        y_plot_arrs = []
        
        for ndx in range(len(y_arrs)):

            # Get x and y arrays. Pull them out of other class type
            # if needed.

            if len(x_arrs) == 1:
                x_arr = x_arrs[0]
            else:
                x_arr = x_arrs[ndx]

            if isinstance(x_arr,mpncv.NCDimensionedVar):
                x_arr = x_arr.data.data
            elif isinstance(x_arr,mpncv.NCVar):
                x_arr = x_arr.data
            
            y_arr = y_arrs[ndx]

            if isinstance(y_arr,mpncv.NCDimensionedVar):
                y_arr = y_arr.data.data
            elif isinstance(y_arr,mpncv.NCVar):
                y_arr = y_arr.data

            x_plot_arr = np.ma.MaskedArray(x_arr,mask=np.ma.mask_or(x_arr.mask,y_arr.mask)).flatten().compressed()
            y_plot_arr = np.ma.MaskedArray(y_arr,mask=np.ma.mask_or(x_arr.mask,y_arr.mask)).flatten().compressed()

            if self.ValFromArgByNm('TriangleSmooth') is not None:
                y_plot_arr = smooth.SmoothTriangle(y_plot_arr,self.ValFromArgByNm('TriangleSmooth'))

            ax1.plot(
                x_plot_arr,
                y_plot_arr,
                color=color_list[ndx],
                linestyle = linestyle_list[ndx],
                label=line_lbls[ndx].replace('_',' ')
                )

            if fill_between:
                # if filling between we will need these
                x_plot_arrs.append(x_plot_arr)
                y_plot_arrs.append(y_plot_arr)
            
        # for ndx in range(len(y_arrs)):

        if fill_between:
            for ndx in range(len(y_arrs)):
                if ndx == 0:
                    ax1.fill_between(
                        x_plot_arrs[ndx],
                        y_plot_arrs[ndx],
                        np.full(y_plot_arrs[0].shape,y_min,dtype=y_plot_arrs[0].dtype),
                        color = color_list[ndx]
                        )
                else:
                    ax1.fill_between(
                        x_plot_arrs[ndx],
                        y_plot_arrs[ndx],
                        y_plot_arrs[ndx-1],
                        )
                    
        ax1.legend()
        
        if out_fnm is None:
            out_fnm = tf.mktemp(suffix='.png')
            plt.savefig(out_fnm)
            os.system('open -a Preview {}'.format(out_fnm))

        else:
            # print 'saving',out_fnm
            plt.savefig(out_fnm)
          
        plt.close(fig) # fig must be closed to get it out of memory

        return True
        
    #     def _XYLinesOnOneGraph(self,x_arrs,y_arrs,line_lbls,out_fnm):
    
    def _BoxAndWhiskerPlot(
        self,
        box_data,
        box_lbls,
        out_fnm,
        ):

        title_size = 14
        x_lbl_size = 12
        y_lbl_size = 12
        tick_lbl_size = 10
        fig_width = 5.1
        fig_height = 4.2

        fig = plt.figure(
            figsize = (
                fig_width,
                fig_height,
                ),
            dpi = 600
            )
        
        ax1 = fig.add_subplot(111)

        ax1.tick_params(
            labelsize = tick_lbl_size,
            )
        
        title = self.ValFromArgByNm('Title')
        if title is None: title = 'Box Plot'

        x_lbl = self.ValFromArgByNm('XLabel')
        if x_lbl is not None:
            ax1.set_xlabel(
                x_lbl.replace('_',' '),
                size = x_lbl_size,
                fontweight = 'bold'
                )

        y_lbl = self.ValFromArgByNm('YLabel')
        if y_lbl is not None:
            ax1.set_ylabel(
                y_lbl.replace('_',' '),
                size = y_lbl_size,
                fontweight = 'bold'
                )

        ax1.set_title(
            title.replace('_',' '),
            fontsize=title_size,
            fontweight='bold'
            )

        y_max = self.ValFromArgByNm('YMax')
        if y_max is None:
            y_max = float('-inf')
            for box in box_data:
                y_max = max(y_max,box.max())

        y_min = self.ValFromArgByNm('YMin')
        if y_min is None:
            y_min = float('inf')
            for box in box_data:
                y_min = min(y_min,box.min())

        ax1.set_ylim(y_min,y_max)
 
        show_fliers = self.ValFromArgByNm('ShowOutliers')
        if show_fliers is None:
            show_fliers = True # default if not specified
        
        show_mean_line = self.ValFromArgByNm('ShowMeanLine')
        if show_mean_line is None:
            show_mean_line = False # default if not specified
        
        # make the actual boxplot
        ax1.boxplot(
            box_data,
            labels = box_lbls,
            showfliers = show_fliers,
            showmeans = show_mean_line,
            meanline = show_mean_line,
            )

        if self.ValFromArgByNm('ShowPoints') == True:
            ngroup = len(box_data)
            
            for ndx in range(len(box_data)):
                box_datum = box_data[ndx]
                xs = [1+ndx + rndm.triangular(-0.05,0.05,0.0) for i in range(len(box_datum))]
                plt.scatter(xs, box_datum, color = 'black',marker = 'x',s=3,alpha=0.4)

        # if self.ValFromArgByNm('ShowPoints') == True:
        
        if out_fnm is None:
            out_fnm = tf.mktemp(suffix='.png')
            plt.savefig(out_fnm)
            os.system('open -a Preview {}'.format(out_fnm))

        else:
            # print 'saving',out_fnm
            plt.savefig(out_fnm)
          
        plt.close(fig) # fig must be closed to get it out of memory

        return True
    
    # def _BoxAndWhiskerPlot(...):

    def _HorizontalBarChart(
        self,
        bar_vals,
        bar_lbls,
        out_fnm
        ):

        title_size = 14
        x_lbl_size = 12
        y_lbl_size = 12
        tick_lbl_size = 8
        fig_width = 6
        fig_height = 4.2

        fig = plt.figure(
            figsize = (
                fig_width,
                fig_height,
                ),
            dpi = 600
            )

        bar_colors = len(bar_vals)/2 * ['black','grey']
        if (len(bar_vals) / 2) * 2 != len(bar_vals):
            bar_colors += ['black']
        
        plt.gcf().subplots_adjust(left = 0.2)
        ax1 = fig.add_subplot(111)

        ax1.tick_params(
            labelsize = tick_lbl_size,
            )
        
        title = self.ValFromArgByNm('Title')
        if title is None: title = 'Bar Chart'

        x_lbl = self.ValFromArgByNm('XLabel')
        if x_lbl is not None:
            ax1.set_xlabel(
                x_lbl.replace('_',' '),
                size = x_lbl_size,
                fontweight = 'bold'
                )

        y_lbl = self.ValFromArgByNm('YLabel')
        if y_lbl is not None:
            ax1.set_ylabel(
                y_lbl.replace('_',' '),
                size = y_lbl_size,
                fontweight = 'bold'
                )

        ax1.set_title(
            title.replace('_',' '),
            fontsize=title_size,
            fontweight='bold'
            )

        x_max = self.ValFromArgByNm('XMax')
        if x_max is None:
            x_max = float('-inf')
            for x_val in bar_vals:
                x_max = max(x_max,x_val)

        x_min = self.ValFromArgByNm('XMin')
        if x_min is None:
            x_min = 0

        ax1.barh(
            np.arange(len(bar_lbls)),
            bar_vals,
            align='center',
            color=bar_colors
            )

        ax1.set_xlim(x_min,x_max)
        ax1.set_yticks(np.arange(len(bar_vals)))
        ax1.set_yticklabels(bar_lbls)
        ax1.invert_yaxis() # labels top to bottom
            
        if out_fnm is None:
            out_fnm = tf.mktemp(suffix='.png')
            plt.savefig(out_fnm)
            os.system('open -a Preview {}'.format(out_fnm))

        else:
            # print 'saving',out_fnm
            plt.savefig(out_fnm)
          
        plt.close(fig) # fig must be closed to get it out of memory

        return True
        
    # def _HorizontalBarChart(...)
    
    def _ScatterPlot(
        self,
        x_vals,
        y_vals,
        pt_lbls,
        out_fnm,
        ):

        title_sz = 14
        title_wt = 'bold'
        pt_lbl_size = 9
        x_lbl_size = 12
        x_lbl_weight = 'bold'
        y_lbl_size = 12
        y_lbl_weight = 'bold'
        tick_lbl_size = 12
        title = self.ValFromArgByNm('Title')
        if title is None: title = 'Scatter Plot'

        x_lbl = self.ValFromArgByNm('XLabel')
        if x_lbl is None:
            x_lbl = 'X Value'
        else:
            x_lbl = x_lbl.replace('_',' ')

        y_lbl = self.ValFromArgByNm('YLabel')
        if y_lbl is None:
            y_lbl = 'Y Value'
        else:
            y_lbl = y_lbl.replace('_',' ')

        # Determine mins and maxes
        x_min = self.ValFromArgByNm('XMin')
        x_line = self.ValFromArgByNm('XLineVal')
        if x_min is None:
            x_min = x_vals.min()
            if x_line is not None:
                x_min = min(x_min, x_line)
                
        x_max = self.ValFromArgByNm('XMax')
        if x_max is None:
            x_max = x_vals.max()
            if x_line is not None:
                x_max = max(x_max, x_line)

        y_min = self.ValFromArgByNm('YMin')
        y_line = self.ValFromArgByNm('YLineVal')
        if y_min is None:
            y_min = y_vals.min()
            if y_line is not None:
                y_min = min(y_min, y_line)
                
        y_max = self.ValFromArgByNm('YMax')
        if y_max is None:
            y_max = y_vals.max()
            if y_line is not None:
                y_max = max(y_max, y_line)


        same_xy_scale = self.ValFromArgByNm('SameXYScale')
        if same_xy_scale:
            x_max = y_max = max(x_max,y_max)
            x_min = y_min = min(x_min,y_min)

        # Create the graph

        plot_wid = self.ValFromArgByNm('PlotWidth')
        if plot_wid is None:
            plot_wid = 5

        plot_ht = self.ValFromArgByNm('PlotHeight')
        if plot_ht is None:
            plot_ht = 5
            
        plt.rcParams["figure.figsize"] = (plot_wid,plot_ht)
        
        fig = plt.figure()
        ax1 = fig.add_subplot(111)
        
        ax1.set_xlabel(
            x_lbl,
            size = x_lbl_size,
            fontweight = x_lbl_weight
            )
        ax1.set_ylabel(
            y_lbl,
            size = y_lbl_size,
            fontweight = y_lbl_weight
            )
        ax1.set_title(
            title.replace('_',' '),
            size = title_sz,
            fontweight = title_wt
            )
        ax1.set_xlim(x_min,x_max)
        ax1.set_ylim(y_min,y_max)
        ax1.tick_params(
            labelsize = tick_lbl_size
            )

        # points
        ax1.scatter(
            x_vals,
            y_vals,
            marker='.',
            s=2,
            color='black'
            )

        # labels

        if pt_lbls is not None:
            for x_val, y_val, pt_lbl in zip (x_vals,y_vals,pt_lbls):
                ax1.text(
                    x_val,
                    y_val,
                    pt_lbl,
                    size = pt_lbl_size,
                    color='black'
                    )
            
        # x and y lines
        if y_line is not None:
            ax1.axvline(
                x=y_line,
                color='black',
                linewidth=0.5
                )

        if x_line is not None:
            ax1.axhline(
                y=x_line,
                color='black',
                linewidth=0.5
                )

        if self.ValFromArgByNm('DoXYLine') is not None:
            xymin = max(x_min,y_min)
            xymax = min(x_max,y_max)
            ax1.plot(
                [xymin,xymax],
                [xymin,xymax],
                linewidth=1,
                color='black'
                )

        left_padding = self.ValFromArgByNm('LeftPadding')
        if left_padding is None:
            left_padding = 0.15

        right_padding = self.ValFromArgByNm('RightPadding')
        if right_padding is None:
            right_padding = 0.05

        bottom_padding = self.ValFromArgByNm('BottomPadding')
        if bottom_padding is None:
            bottom_padding = 0.1

        top_padding = self.ValFromArgByNm('TopPadding')
        if top_padding is None:
            top_padding = 0.05

        plt.subplots_adjust(
            left=left_padding,
            bottom=bottom_padding,
            right=(1.0-right_padding),
            top=(1.0-top_padding),
            wspace=0,
            hspace=0
            ) 

        if out_fnm is None:
            out_fnm = tf.mktemp(suffix='.png')
            plt.savefig(out_fnm)
            os.system('open -a Preview {}'.format(out_fnm))

        else:
            # print 'saving',out_fnm
            plt.savefig(out_fnm)

        plt.close(fig) # fig must be closed to get it out of memory

        return True

     # def _ScatterPlot(...)

    def _RenderLayerAtResolution(
        self,
        data_obj,
        min_val,
        max_val,
        in_fldnm,
        out_fnm,
        ):

        # Renders a layer with one pixel per cell
        
        # Here is the numpy array to render
        tmp_arr = data_obj.ExecRslt().data.data

        if tmp_arr.ndim != 2:
            print (3*'{}\n').format(
                'Warning: Trying to render non 2-D layer.',
                '  Variable name: {}'.format(self.RsltNm()),
                '  Skipping...'
                )
            
            return False
        
        if out_fnm is None:
            out_fnm = tf.mktemp(suffix='.png')
        # else:
            #print 'saving',out_fnm

        reverse_ramp = self.ValFromArgByNm('ReverseColorMap')
        if reverse_ramp != True:
            reverse_ramp = False

        color_ramp_nm = self.ValFromArgByNm('ColorMap')
        if color_ramp_nm  is not None:
            color_ramp = nptopng.color_ramp(
                color_ramp_nm,
                reverse_ramp = reverse_ramp
                )
        else:
            color_ramp = None

        flip_up_down = self.ValFromArgByNm('FlipX')
        if flip_up_down != True:
            flip_up_down = False
            
        flip_left_right = self.ValFromArgByNm('FlipY')
            
        if flip_left_right != True:
            flip_left_right = False

        missing_color = self.ValFromArgByNm('FaceColor')
        if missing_color is not None:
            missing_color = self.ValFromArgByNm('FaceColor').split(':')

        nptopng.render_array_to_png(
            tmp_arr,
            out_fnm,
            value_range = [min_val,max_val],
            clr_ramp = color_ramp,
            missing_color = missing_color,
            flip_up_down = flip_up_down,
            flip_left_right = flip_left_right,
            )

        # RenderKey: Make and save an array for a key (i.e. bar with colors),
        # and a file giving values for n divisions of the bar, so when you
        # manually add divisions to the bar you can just cut and paste

        render_key = self.ValFromArgByNm('RenderKey')

        if render_key is not None:

            color_key_short_dim_sz = 51
            color_key_long_dim_sz = 512

            color_key_arr = np.zeros([color_key_short_dim_sz,color_key_long_dim_sz])
            color_key_arr[0:color_key_short_dim_sz,:] = [ \
                min_val + float ( x ) * ( max_val-min_val ) / \
                float ( color_key_arr.shape[1] - 1) \
                for x in range ( color_key_arr.shape[1] ) \
            ]

            if render_key == 'Horizontal':
                pass
            elif render_key == 'Vertical':
                color_key_arr = np.flip(np.transpose(color_key_arr),0)
                
            nptopng.render_array_to_png(
                color_key_arr,
                out_fnm.replace('.png','_key.png'),
                value_range = [min_val,max_val],
                clr_ramp = color_ramp,
                )

            with open(out_fnm.replace('.png','_key_divisions.txt'),'w') as out_f:
                
                for num_divisions in range(1,7):
                    step_val = (max_val - min_val) / float(num_divisions)
                    out_f.write('{} divisions\n'.format(num_divisions))
                    out_f.write('{}\n'.format(','.join([str(min_val + x * step_val) for x in range(num_divisions+1)])))

            # if render_key is not None:
            
        return True

    # def _RenderLayerAtResolution(self,...):
     
# class _GraphicsParent(_NetCDFUtilParent):

class RenderLayers(_GraphicsParent):
       
    def __init__(
        self,
        mpt_cmd_struct=None,
        ):

        super(RenderLayers,self).__init__(mpt_cmd_struct)

    # def __init__(...)

    def _SetFxnDesc(self):
        
        self.fxnDesc['DisplayName'] = 'RenderLayers'
        self.fxnDesc['ShortDesc'] = 'Renders layer'
        self.fxnDesc['ReturnType'] = 'Boolean'
        
        self.fxnDesc['ReqArgs'] = {
            'InFieldNames':['Field Name','Field Name List'],
            'OutDirectoryName':'File Name'
            }
        self.fxnDesc['OptArgs'] = {
            'PrecursorFieldNames':['Field Name','Field Name List'],
            'PlainImageOnly':'Boolean',
            'FileNamePrepend':'Any',
            'FileNameAppend':'Any',
            'Title':'Any',
            'ColorMap':'Any',
            'FaceColor':'Tuple: Integer:Integer:Integer',
            'FaceColor':'Any',
            'Origin':'Any',
            'FlipX':'Boolean',
            'FlipY':'Boolean',
            'MinVal':'Float',
            'MaxVal':'Float',
            'TitleSize':'Float',
            'TitleOffset':'Float',
            'TickLabelSize':'Float',
            'TickWidth':'Float',
            'TickLength':'Float',
            }
        
    # _SetFxnDesc(self):
    
    def DependencyNms(self):
        
        rtrn = self._ArgToList('InFieldNames')
        rtrn += self._ArgToList('PrecursorFieldNames')
        return rtrn
    
    # def DependencyNms(self):
        
    def Exec(self,executedObjects):

        self._InsureDirExists(self.ArgByNm('OutDirectoryName'))

        # Get the common min and max

        min_val = self.ValFromArgByNm('MinVal')
        if min_val is None:
            min_val = float('inf')
            for in_fldnm in self._ArgToList('InFieldNames'):
                min_val = min(min_val,executedObjects[in_fldnm].ExecRslt().data.data.min())
                
        max_val = self.ValFromArgByNm('MaxVal')
        if max_val is None:
            max_val = float('-inf')
            for in_fldnm in self._ArgToList('InFieldNames'):
                max_val = max(max_val,executedObjects[in_fldnm].ExecRslt().data.data.max())
        
        # For each variable render the layer and write it
        for in_fldnm in self._ArgToList('InFieldNames'):

            prepend_str = self.ArgByNm('FileNamePrepend')
            if prepend_str is None: prepend_str = ''
            append_str = self.ArgByNm('FileNameAppend')
            if append_str is None: append_str = ''

            file_nm = '{}{}{}'.format(prepend_str,in_fldnm,append_str)
            
            # print 'Rendering {}'.format(file_nm)
            self._RenderLayer(
                executedObjects[in_fldnm],
                min_val,
                max_val,
                in_fldnm,
                '{}/{}.png'.format(self.ArgByNm('OutDirectoryName'),file_nm)
                )

        self.execRslt = True
        executedObjects[self.RsltNm()] = self

    # def Exec(self,executedObjects):
        
# class RenderLayers(_GraphicsParent):

class RenderLayersAtResolution(_GraphicsParent):
       
    def __init__(
        self,
        mpt_cmd_struct=None,
        ):

        super(RenderLayersAtResolution,self).__init__(mpt_cmd_struct)

    # def __init__(...)

    def _SetFxnDesc(self):
        
        self.fxnDesc['DisplayName'] = 'RenderLayersKeepResolution'
        self.fxnDesc['ShortDesc'] = 'Renders layers at their data resolution'
        self.fxnDesc['ReturnType'] = 'Boolean'
        
        self.fxnDesc['ReqArgs'] = {
            'InFieldNames':['Field Name','Field Name List'],
            'OutDirectoryName':'File Name'
            }
        self.fxnDesc['OptArgs'] = {
            'PrecursorFieldNames':['Field Name','Field Name List'],
            'FileNamePrepend':'Any',
            'FileNameAppend':'Any',
            'ColorMap':'Any',
            'ReverseColorMap':'Boolean',
            'FaceColor':'Tuple: Integer:Integer:Integer',
            'FlipX':'Boolean',
            'FlipY':'Boolean',
            'MinVal':'Float',
            'MaxVal':'Float',
            'RenderKey':'One of| Horizontal Vertical',
            }
        
    # _SetFxnDesc(self):
    
    def DependencyNms(self):
        
        rtrn = self._ArgToList('InFieldNames')
        rtrn += self._ArgToList('PrecursorFieldNames')
        return rtrn
    
    # def DependencyNms(self):
        
    def Exec(self,executedObjects):

        self._InsureDirExists(self.ArgByNm('OutDirectoryName'))

        # Get the common min and max

        min_val = self.ValFromArgByNm('MinVal')
        if min_val is None:
            min_val = float('inf')
            for in_fldnm in self._ArgToList('InFieldNames'):
                min_val = min(min_val,executedObjects[in_fldnm].ExecRslt().data.data.min())
                
        max_val = self.ValFromArgByNm('MaxVal')
        if max_val is None:
            max_val = float('-inf')
            for in_fldnm in self._ArgToList('InFieldNames'):
                max_val = max(max_val,executedObjects[in_fldnm].ExecRslt().data.data.max())
                
        # For each variable render the layer and write it
        for in_fldnm in self._ArgToList('InFieldNames'):
            
            prepend_str = self.ArgByNm('FileNamePrepend')
            if prepend_str is None: prepend_str = ''
            append_str = self.ArgByNm('FileNameAppend')
            if append_str is None: append_str = ''

            file_nm = '{}{}{}'.format(prepend_str,in_fldnm,append_str)
            
            # print 'Rendering {}'.format(file_nm)

            
            self._RenderLayerAtResolution(
                executedObjects[in_fldnm],
                min_val,
                max_val,
                in_fldnm,
                '{}/{}.png'.format(self.ArgByNm('OutDirectoryName'),file_nm)
                )

        self.execRslt = True
        executedObjects[self.RsltNm()] = self

    # def Exec(self,executedObjects):
        
# class RenderLayersAtResolution(_GraphicsParent):

class Graph1DVsDim(_GraphicsParent):

    def __init__(
        self,
        mpt_cmd_struct=None,
        ):

        super(Graph1DVsDim,self).__init__(mpt_cmd_struct)

    # def __init__(...)

    def _SetFxnDesc(self):
        
        self.fxnDesc['DisplayName'] = 'Graph1DVsDim'
        self.fxnDesc['ShortDesc'] = 'Line graph of variable(s) versus a dimension.'
        self.fxnDesc['ReturnType'] = 'Boolean'
        
        self.fxnDesc['ReqArgs'] = {
            'InFieldNames':['Field Name','Field Name List'],
            'OutFileName':'File Name',
            }
        self.fxnDesc['OptArgs'] = {
            'Metadata':'Any',
            'PrecursorFieldNames':['Field Name','Field Name List'],
            'Title':'Any',
            'XLabel':'Any',
            'YLabel':'Any',
            'XMin':'Float',
            'XMax':'Float',
            'YMin':'Float',
            'YMax':'Float',
            'Colors':'Any',
            'LineStyles':'Any',
            'LineTitles':'Any',
            'TriangleSmooth':'Integer'
            }
        
    # _SetFxnDesc(self):

    def DependencyNms(self):
        
        rtrn = self._ArgToList('InFieldNames')
        rtrn += self._ArgToList('PrecursorFieldNames')
        return rtrn
    
    # def DependencyNms(self):
        
    def Exec(self,executedObjects):

        x_arrs = []
        y_arrs = []
        line_titles = []

        for in_fldnm in self._ArgToList('InFieldNames'):

            line_titles.append(in_fldnm)

            in_marr = cp.deepcopy(executedObjects[in_fldnm].ExecRslt().data.data)

            if in_marr.ndim != 1:

                raise Exception(
                    (4*'{}\n').format(
                        'Array is not 1D',
                        '  Field name, variable name, shape: {}  {}  {}'.format(
                            in_fldnm,
                            executedObjects[in_fldnm].ExecRslt().name,
                            in_marr.shape
                            ),
                        'File: {}  Line number: {}'.format(
                            self.mptCmdStruct['cmdFileNm'],self.mptCmdStruct['lineNo']
                            ),
                        'Full command:\n{}'.format(self.mptCmdStruct['rawCmdStr'])
                        ),
                    )

            # if in_marr.shape != 1:
            
            dim_marr = cp.deepcopy(executedObjects[in_fldnm].ExecRslt().dims.values()[0].data)
            tmp_mask = np.ma.mask_or(dim_marr.mask,in_marr.mask)

            dim_marr.mask = tmp_mask
            in_marr.mask = tmp_mask

            # x_arrs.append(in_marr.compressed())
            # y_arrs.append(dim_marr.compressed())

            y_arrs.append(in_marr)
            x_arrs.append(dim_marr)

        # for dataKey,dataObj in dataObjs.items():

        # Set some parameters for graph:
        title = self.ValFromArgByNm('Title')
        if title is None:
            title = 'Y vs X'
        else:
            title = title.replace('_',' ')
            
        x_lbl = self.ValFromArgByNm('XLabel')
        if x_lbl is None:
            x_lbl = 'X'
        else:
            x_lbl = x_lbl.replace('_',' ')

        y_lbl = self.ValFromArgByNm('YLabel')
        if y_lbl is None:
            y_lbl = 'Y'
        else:
            y_lbl = y_lbl.replace('_',' ')

        self.execRslt = self._XYLinesOnOneGraph(
            x_arrs,
            y_arrs,
            title,
            x_lbl,
            y_lbl,
            line_titles,
            self.ArgByNm('OutFileName')
            )

        executedObjects[self.RsltNm()] = self
    
    # def Exec(self,executedObjects):
    
# class Graph1DVsDim(_GraphicsParent):

class GraphLinesTogether(_GraphicsParent):

    def __init__(
        self,
        mpt_cmd_struct=None,
        ):

        super(GraphLinesTogether,self).__init__(mpt_cmd_struct)

    # def __init__(...)

    def _SetFxnDesc(self):
        
        self.fxnDesc['DisplayName'] = 'GraphLinesTogether'
        self.fxnDesc['ShortDesc'] = 'Line graph of variable(s) versus other variable(s).'
        self.fxnDesc['ReturnType'] = 'Boolean'
        
        self.fxnDesc['ReqArgs'] = {
            'XFieldNames':['Field Name','Field Name List'],
            'YFieldNames':['Field Name','Field Name List'],
            'OutFileName':'File Name',
            }
        self.fxnDesc['OptArgs'] = {
            'Metadata':'Any',
            'PrecursorFieldNames':['Field Name','Field Name List'],
            'Title':'Any',
            'XLabel':'Any',
            'YLabel':'Any',
            'XMin':'Float',
            'XMax':'Float',
            'YMin':'Float',
            'YMax':'Float',
            'Colors':'Any',
            'LineStyles':'Any',
            'LineTitles':'Any',
            'TriangleSmooth':'Integer',
            'DoStacked':'Boolean',
            'DoFillBetween':'Boolean',
            }
        
    # _SetFxnDesc(self):

    def DependencyNms(self):
        
        rtrn = self._ArgToList('XFieldNames')
        rtrn += self._ArgToList('YFieldNames')
        rtrn += self._ArgToList('PrecursorFieldNames')
        return rtrn
    
    # def DependencyNms(self):
        
    def Exec(self,executedObjects):

        x_arrs = []
        y_arrs = []

        if len(self._ArgToList('XFieldNames')) != len(self._ArgToList('YFieldNames')):

            raise Exception(
                (5*'{}\n').format(
                    'XFieldNames and YFieldNames must correspond 1 to 1',
                    'XFieldNames: {}'.format(', '.join(self._ArgToList('XFieldNames'))),
                    'YFieldNames: {}'.format(', '.join(self._ArgToList('YFieldNames'))),
                    'File: {}  Line number: {}'.format(
                        self.mptCmdStruct['cmdFileNm'],self.mptCmdStruct['lineNo']
                        ),
                    'Full command:\n{}'.format(self.mptCmdStruct['rawCmdStr'])
                    ),
                )

        line_titles = self._ArgToList('LineTitles')
        if line_titles == None:
            line_titles = self._ArgToList('XFieldNames')

        if len(self._ArgToList('XFieldNames')) != len(self._ArgToList('LineTitles')):
            raise Exception(
                (5*'{}\n').format(
                    'XFieldNames and LineTitles must correspond 1 to 1',
                    'XFieldNames: {}'.format(', '.join(self._ArgToList('XFieldNames'))),
                    'LineTitles: {}'.format(', '.join(self._ArgToList('LineTitles'))),
                    'File: {}  Line number: {}'.format(
                        self.mptCmdStruct['cmdFileNm'],self.mptCmdStruct['lineNo']
                        ),
                    'Full command:\n{}'.format(self.mptCmdStruct['rawCmdStr'])
                    ),
                )
        

        for x_fldnm,y_fldnm in zip(self._ArgToList('XFieldNames'), self._ArgToList('YFieldNames')):

            x_marr = executedObjects[x_fldnm].ExecRslt().data.data
            y_marr = executedObjects[y_fldnm].ExecRslt().data.data

            if x_marr.shape != y_marr.shape or x_marr.ndim != 1:
                raise Exception(
                    (4*'{}\n').format(
                        'X and Y fields must be 1D and have identical shapes',
                        '  XFieldName, shape: {}  {}  {}'.format(
                            x_fldnm,
                            x_marr.shape
                            ),
                        '  YFieldName, shape: {}  {}  {}'.format(
                            y_fldnm,
                            y_marr.shape
                            ),
                        'File: {}  Line number: {}'.format(
                            self.mptCmdStruct['cmdFileNm'],self.mptCmdStruct['lineNo']
                            ),
                        'Full command:\n{}'.format(self.mptCmdStruct['rawCmdStr'])
                        ),
                    )

            # if x_marr.shape != y_marr.shape or if x_marr.ndim != 1:
            
            # Combine masks
            tmp_mask = np.ma.mask_or(x_marr.mask,y_marr.mask)

            x_arrs.append(np.ma.array(x_marr,mask=tmp_mask,copy=True))
            y_arrs.append(np.ma.array(y_marr,mask=tmp_mask,copy=True))

        # for dataKey,dataObj in dataObjs.items():

        if self.ValFromArgByNm('DoStacked'):
            for ndx in range(1,len(y_arrs)):
                y_arrs[ndx] += y_arrs[ndx-1]

            fill_between = False
            if self.ValFromArgByNm('DoFillBetween'):
                fill_between = True
                
        # if self.ValFromArgByNm('DoStacked'):

        # Set some parameters for graph:
        title = self.ValFromArgByNm('Title')
        if title is None:
            title = 'Y vs X'
        else:
            title = title.replace('_',' ')
            
        x_lbl = self.ValFromArgByNm('XLabel')
        if x_lbl is None:
            x_lbl = 'X'
        else:
            x_lbl = x_lbl.replace('_',' ')

        y_lbl = self.ValFromArgByNm('YLabel')
        if y_lbl is None:
            y_lbl = 'Y'
        else:
            y_lbl = y_lbl.replace('_',' ')

        self.execRslt = self._XYLinesOnOneGraph(
            x_arrs,
            y_arrs,
            title,
            x_lbl,
            y_lbl,
            line_titles,
            self.ArgByNm('OutFileName'),
            fill_between = fill_between,
            )

        executedObjects[self.RsltNm()] = self
    
    # def Exec(self,executedObjects):
    
# class GraphLinesTogether(_GraphicsParent):

class GraphDistributionsTogether(_GraphicsParent):

    def __init__(
        self,
        mpt_cmd_struct=None,
        ):

        super(GraphDistributionsTogether,self).__init__(mpt_cmd_struct)

    # def __init__(...)

    def _SetFxnDesc(self):
        
        self.fxnDesc['DisplayName'] = 'GraphDistributionsTogether'
        self.fxnDesc['ShortDesc'] = 'Creates a set of connected line distributions for one or more variables on one graph.'
        self.fxnDesc['ReturnType'] = 'Boolean'
        
        self.fxnDesc['ReqArgs'] = {
            'InFieldNames':['Field Name','Field Name List'],
            'OutFileName':'File Name',
            }
        self.fxnDesc['OptArgs'] = {
            'Metadata':'Any',
            'PrecursorFieldNames':['Field Name','Field Name List'],
            'Normalize':'Boolean',
            'Title':'Any',
            'XLabel':'Any',
            'YLabel':'Any',
            'XMin':'Float',
            'XMax':'Float',
            'YMin':'Float',
            'YMax':'Float',
            'Colors':'Any',
            'LineStyles':'Any',
            'LineTitles':'Any',
            'NumberOfBins':'Integer'
            }
        
    # _SetFxnDesc(self):
 
    def DependencyNms(self):
        
        rtrn = self._ArgToList('InFieldNames')
        rtrn += self._ArgToList('PrecursorFieldNames')
        return rtrn
    
    # def DependencyNms(self):
        
    def Exec(self,executedObjects):

        x_arrs = []
        y_arrs = []
        line_titles = []

        num_bins = self.ValFromArgByNm('NumberOfBins')
        if num_bins is None:
            num_bins = 50

        # Prep arrays for plotting and get the min and max for the distribution
        arrays_to_plot = OrderedDict()
        bin_min = float('inf')
        bin_max = float('-inf')
        
        for in_fldnm in self._ArgToList('InFieldNames'):

            line_titles.append(in_fldnm)
            
            in_arr = executedObjects[in_fldnm].ExecRslt().data.data
            
            # Flatten and compress the arrays
            flat_arr = in_arr.flatten().compressed()

            # Save array for plotting
            arrays_to_plot[in_fldnm] = flat_arr

            # Update min and max for plots
            bin_min = min(bin_min,in_arr.min())
            bin_max = max(bin_max,in_arr.max())
             
        # for in_fldnm in self._ArgToList('InFieldNames'):

        # Get the edges defining the divisions for the distribution
        bin_wid = (bin_max - bin_min) / num_bins
        
        edges = np.array([bin_min + tmp_ndx * bin_wid for tmp_ndx in range(num_bins+1)])
        ctrs = 0.5 * (edges[1:] + edges[:-1])
        edges[0] -= bin_wid/10000. # Makes the loop for binning work right

        # Make the quantity arrays for the plot
        for in_fldnm in self._ArgToList('InFieldNames'):

            quantities = np.zeros(num_bins)
            # Calculate the quantities
            for bin_ndx in range(num_bins):
                quantities[bin_ndx] = np.where(
                    np.logical_and(edges[bin_ndx] < arrays_to_plot[in_fldnm], arrays_to_plot[in_fldnm] <= edges[bin_ndx+1]),
                    1,
                    0
                    ).sum()

            # Normalize if specified

            if self.ValFromArgByNm('Normalize'):
                quantities /= quantities.sum()
                

            # Add the x and y arrays to the arrays that will be plotted
            x_arrs.append(np.ma.MaskedArray(ctrs,mask=False))
            y_arrs.append(np.ma.MaskedArray(quantities,mask=False))

        # for in_fldnm in self._ArgToList('InFieldNames'):

        # for ndx in range(len(x_arrs)):
            # print line_titles[ndx]

            # print '  x',x_arrs[ndx].min(),x_arrs[ndx].max()
            # print '  y',y_arrs[ndx].min(),y_arrs[ndx].max()

        # Set some parameters for graph:
        title = self.ValFromArgByNm('Title')
        if title is None:
            title = 'Distribution'
        else:
            title = title.replace('_',' ')
           
        x_lbl = self.ValFromArgByNm('XLabel')
        if x_lbl is None:
            x_lbl = 'Value'
        else:
            x_lbl = x_lbl.replace('_',' ')

        y_lbl = self.ValFromArgByNm('YLabel')
        if y_lbl is None:
            if self.ValFromArgByNm('Normalize') == True:
                y_lbl = 'Fraction'
            else:
                y_lbl = 'Count'
        else:
            y_lbl = y_lbl.replace('_',' ')

        self.execRslt = self._XYLinesOnOneGraph(
            x_arrs,
            y_arrs,
            title,
            x_lbl,
            y_lbl,
            line_titles,
            self.ArgByNm('OutFileName')
            )

        executedObjects[self.RsltNm()] = self

    # def Exec(self,executedObjects):
            
# class GraphDistributionsTogether(_GraphicsParent):

class GraphWeightedDistributions(_GraphicsParent):

    def __init__(
        self,
        mpt_cmd_struct=None,
        ):

        super(GraphWeightedDistributions,self).__init__(mpt_cmd_struct)

    # def __init__(...)

    def _SetFxnDesc(self):
        
        self.fxnDesc['DisplayName'] = 'GraphWeightedDistributions'
        self.fxnDesc['ShortDesc'] = 'Creates a set of connected line distributions for one or more variables on one graph, weighted by a weighting array.'
        self.fxnDesc['ReturnType'] = 'Boolean'
        
        self.fxnDesc['ReqArgs'] = {
            'InFieldNames':['Field Name','Field Name List'],
            'WeightFieldName':'Field Name',
            'OutFileName':'File Name',
            }
        self.fxnDesc['OptArgs'] = {
            'Metadata':'Any',
            'PrecursorFieldNames':['Field Name','Field Name List'],
            'Normalize':'Boolean',
            'Title':'Any',
            'XLabel':'Any',
            'YLabel':'Any',
            'XMin':'Float',
            'XMax':'Float',
            'YMin':'Float',
            'YMax':'Float',
            'Colors':'Any',
            'LineStyles':'Any',
            'LineTitles':'Any',
            'NumberOfBins':'Integer'
            }
        
    # _SetFxnDesc(self):

    def DependencyNms(self):
        
        rtrn = self._ArgToList('InFieldNames')
        rtrn += self._ArgToList('WeightFieldName')
        rtrn += self._ArgToList('PrecursorFieldNames')
        return rtrn
    
    # def DependencyNms(self):
        
    def Exec(self,executedObjects):

        x_arrs = []
        y_arrs = []
        line_titles = []

        weight_ncdimvar = executedObjects[self.ValFromArgByNm('WeightFieldName')].ExecRslt()
        
        num_bins = self.ValFromArgByNm('NumberOfBins')
        if num_bins is None:
            num_bins = 50

        # Prep arrays for plotting and get the min and max for the distribution
        arrays_to_plot = OrderedDict()
        bin_min = float('inf')
        bin_max = float('-inf')
        for in_fldnm in self._ArgToList('InFieldNames'):

            line_titles.append(in_fldnm)
            
            in_ncdimvar = cp.deepcopy(executedObjects[in_fldnm].ExecRslt())
            in_marr = in_ncdimvar.data.data
            
            # make a weight array that is congruent to in_arr

            if ac.array_is_congruent_to_tgt(
                weight_ncdimvar.data.data,
                weight_ncdimvar.dim_nms,
                in_ncdimvar.data.data,
                in_ncdimvar.dim_nms
                ):

                tmp_wt_marr = cp.deepcopy(weight_ncdimvar.data.data)

            elif ac.can_expand_to_tgt(
                weight_ncdimvar.data.data,
                weight_ncdimvar.dim_nms,
                in_ncdimvar.data.data,
                in_ncdimvar.dim_nms
                ):

                tmp_wt_marr =  ac.expand_arr_to_match_tgt(
                    weight_ncdimvar.data.data,
                    weight_ncdimvar.dim_nms,
                    in_ncdimvar.data.data,
                    in_ncdimvar.dim_nms
                    )

            elif ac.can_transpose_to_tgt(
                weight_ncdimvar.data.data,
                weight_ncdimvar.dim_nms,
                in_ncdimvar.data.data,
                in_ncdimvar.dim_nms
                ):
            
                tmp_wt_marr  = ac.transpose_arr_to_match_tgt(
                    weight_ncdimvar.data.data,
                    weight_ncdimvar.dim_nms,
                    in_ncdimvar.data.data,
                    in_ncdimvar.dim_nms
                    )

            else:

                raise Exception(
                    (4*'{}\n').format(
                        'Cannot make weight congruent to variable to be weighted:'.format(dim_nm),
                        '  Weight name, variable name: {}  {}'.format(
                            self.ValFromArgByNm('WeightFieldName'),
                            in_fldnm
                            ),
                        'File: {}  Line number: {}'.format(
                            self.mptCmdStruct['cmdFileNm'],self.mptCmdStruct['lineNo']
                            ),
                        'Full command:\n{}'.format(self.mptCmdStruct['rawCmdStr'])
                        ),
                    )

            # if ac.array_is_congruent_to_tgt(...)...elif...else
                
            # Make masks congruent

            new_mask = np.ma.mask_or(
                tmp_wt_marr.mask,
                in_marr.mask
                )

            tmp_wt_marr.mask[:] = new_mask[:]
            in_marr.mask[:] = new_mask[:]

            # Flatten and compress the arrays
            tmp_wt_arr = tmp_wt_marr.flatten().compressed()
            in_arr = in_marr.flatten().compressed()

            # Save array for plotting
            arrays_to_plot[in_fldnm] = in_arr

            # Update min and max for plots
            bin_min = min(bin_min,in_arr.min())
            bin_max = max(bin_max,in_arr.max())
             
        # for in_fldnm in self._ArgToList('InFieldNames'):

        # Get the edges defining the divisions for the distribution
        bin_wid = (bin_max - bin_min) / num_bins
        
        edges = np.array([bin_min + tmp_ndx * bin_wid for tmp_ndx in range(num_bins+1)])
        ctrs = 0.5 * (edges[1:] + edges[:-1])
        edges[0] -= bin_wid/10000. # Makes the loop for binning work right

        # Make the quantity arrays for the plot
        for in_fldnm in self._ArgToList('InFieldNames'):

            quantities = np.zeros(num_bins)
            # Calculate the quantities
            for bin_ndx in range(num_bins):
                quantities[bin_ndx] = np.where(
                    np.logical_and(edges[bin_ndx] < arrays_to_plot[in_fldnm], arrays_to_plot[in_fldnm] <= edges[bin_ndx+1]),
                    tmp_wt_arr,
                    0
                    ).sum()

            # Normalize if specified

            if self.ValFromArgByNm('Normalize'):
                quantities /= quantities.sum()

            # Add the x and y arrays to the arrays that will be plotted
            x_arrs.append(np.ma.MaskedArray(ctrs,mask=False))
            y_arrs.append(np.ma.MaskedArray(quantities,mask=False))

        # for in_fldnm in self._ArgToList('InFieldNames'):

        # Set some parameters for graph:
        title = self.ValFromArgByNm('Title')
        if title is None:
            title = 'Distribution'
        else:
            title = title.replace('_',' ')
           
        x_lbl = self.ValFromArgByNm('XLabel')
        if x_lbl is None:
            x_lbl = 'Value'
        else:
            x_lbl = x_lbl.replace('_',' ')

        y_lbl = self.ValFromArgByNm('YLabel')
        if y_lbl is None:
            if self.ValFromArgByNm('Normalize') == True:
                y_lbl = 'Fraction'
            else:
                y_lbl = 'Count'
        else:
            y_lbl = y_lbl.replace('_',' ')
            
        self.execRslt = self._XYLinesOnOneGraph(
            x_arrs,
            y_arrs,
            title,
            x_lbl,
            y_lbl,
            line_titles,
            self.ArgByNm('OutFileName')
            )

        executedObjects[self.RsltNm()] = self
        
# class GraphWeightedDistributions(_GraphicsParent):

class ScatterPlot(_GraphicsParent):

    def __init__(
        self,
        mpt_cmd_struct=None,
        ):

        super(ScatterPlot,self).__init__(mpt_cmd_struct)

    # def __init__(...)

    def _SetFxnDesc(self):
        
        self.fxnDesc['DisplayName'] = '2-d Scatter Plot'
        self.fxnDesc['ShortDesc'] = 'Creates a scatter plot. Has option for point labels.'
        self.fxnDesc['ReturnType'] = 'Boolean'
        
        self.fxnDesc['ReqArgs'] = {
            'XFieldName':'Field Name',
            'YFieldName':'Field Name',
            'OutFileName':'File Name',
            }
        self.fxnDesc['OptArgs'] = {
            'PrecursorFieldNames':['Field Name','Field Name List'],
            'Title':'Any',
            'PointLabels':'Any List',
            'XLabel':'Any',
            'YLabel':'Any',
            'XMin':'Float',
            'XMax':'Float',
            'YMin':'Float',
            'YMax':'Float',
            'SameXYScale':'Boolean',
            'XLineVal':'Float',
            'YLineVal':'Float',
            'LeftPadding':'Float',
            'RightPadding':'Float',
            'BottomPadding':'Float',
            'TopPadding':'Float',
            'DoXYLine':'Boolean',
            'PlotWidth':'Float',
            'PlotHeight':'Float',
            }
        
    # _SetFxnDesc(self):

    def DependencyNms(self):
        
        rtrn = self._ArgToList('XFieldName')
        rtrn += self._ArgToList('YFieldName')
        rtrn += self._ArgToList('PrecursorFieldNames')
        return rtrn
    
    # def DependencyNms(self):

    def Exec(self,executedObjects):

        x_arr = executedObjects[self.ValFromArgByNm('XFieldName')].ExecRslt().data.data.compressed()
        y_arr = executedObjects[self.ValFromArgByNm('YFieldName')].ExecRslt().data.data.compressed()

        if x_arr.shape != y_arr.shape:
            raise Exception(
                (5*'{}\n').format(
                    'X and Y arrays must have identical numbers of unmasked elements'.format(dim_nm),
                    '  X field name and number unmasked elements: {}  {}'.format(
                        self.ValFromArgByNm('XFieldName'),
                        x_arr.shape
                        ),
                    '  Y field name and  number unmasked elements: {}  {}'.format(
                        self.ValFromArgByNm('YFieldName'),
                        y_arr.shape
                        ),
                    'File: {}  Line number: {}'.format(
                        self.mptCmdStruct['cmdFileNm'],self.mptCmdStruct['lineNo']
                        ),
                    'Full command:\n{}'.format(self.mptCmdStruct['rawCmdStr'])
                    ),
                )

        if self.ArgByNm('PointLabels') is not None:
            point_lbls = self._ArgToList('PointLabels')
        else:
            point_lbls = None
            
        if point_lbls is not None:
            if len(point_lbls) != x_arr.size:
                raise Exception(                
                    (4*'{}\n').format(
                        'Arrays and Point lables must have the same number of valid elements.',
                        '  Number of valid elements in arrays, labels: {}  {}'.format(
                            x_arr.size,
                            len(point_lbls)
                            ),
                        'File: {}  Line number: {}'.format(
                            self.mptCmdStruct['cmdFileNm'],self.mptCmdStruct['lineNo']
                            ),
                        'Full command:\n{}'.format(self.mptCmdStruct['rawCmdStr'])
                        ),
                    )
                
        # if point_lbls is not None:

        self.execRslt = self._ScatterPlot(
            x_arr,
            y_arr,
            point_lbls,
            self.ArgByNm('OutFileName')
            )

        executedObjects[self.RsltNm()] = self

    # def Exec(self,executedObjects):

# class ScatterPlot(_GraphicsParent):

class BoxAndWhiskerPlot(_GraphicsParent):

    def __init__(
        self,
        mpt_cmd_struct=None,
        ):

        super(BoxAndWhiskerPlot,self).__init__(mpt_cmd_struct)

    # def __init__(...)

    def _SetFxnDesc(self):
        
        self.fxnDesc['DisplayName'] = 'Box and Whisker Plot'
        self.fxnDesc['ShortDesc'] = 'Creates a box and whisker for each array specified.'
        self.fxnDesc['ReturnType'] = 'Boolean'
        
        self.fxnDesc['ReqArgs'] = {
            'InFieldNames':['Field Name','Field Name List'],
            'OutFileName':'File Name',
            }
        self.fxnDesc['OptArgs'] = {
            'Metadata':'Any',
            'PrecursorFieldNames':['Field Name','Field Name List'],
            'Title':'Any',
            'BoxLabels':'Any List',
            'XLabel':'Any',
            'YLabel':'Any',
            'YMin':'Float',
            'YMax':'Float',
            'ShowOutliers':'Boolean',
            'ShowMeanLine':'Boolean',
            'ShowPoints':'Boolean',
            }
        
    # _SetFxnDesc(self):
 
    def DependencyNms(self):
        
        rtrn = self._ArgToList('InFieldNames')
        rtrn += self._ArgToList('PrecursorFieldNames')
        return rtrn
    
    # def DependencyNms(self):
        
    def Exec(self,executedObjects):

        box_lbls = self._ArgToList('BoxLabels')
        in_fldnms = self._ArgToList('InFieldNames')

        if len(box_lbls) is not 0:
            if len(box_lbls) != len(in_fldnms):
                raise Exception(
                    (5*'{}\n').format(
                        '\n********************ERROR********************\n',
                        'Number of box labels does not match number of fields.',
                        '  InFieldNames: {}'.format(' '.join(in_fldnms)),
                        '  BoxLabels: {}'.format(' '.join(box_lbls)),
                        'Full command:\n{}'.format(self.mptCmdStruct['rawCmdStr'])
                    ),
                )

            box_lbls = [box_lbl.replace('_',' ') for box_lbl in box_lbls]

            # if len(box_lbls) != len(in_fldnms):
            
        else:
              box_lbls = in_fldnms
                        
        # if box_lbls is not None:...else...

        box_data = [executedObjects[in_fldnm].ExecRslt().data.data.compressed() for in_fldnm in in_fldnms]

        self.execRslt = self._BoxAndWhiskerPlot(
            box_data,
            box_lbls,
            self.ValFromArgByNm('OutFileName'),
            )

        executedObjects[self.RsltNm()] = self

    # def Exec(self,executedObjects):
            
# class BoxAndWhiskerPlot(_GraphicsParent):

class HorizontalBarChart(_GraphicsParent):

    def __init__(
        self,
        mpt_cmd_struct=None,
        ):

        super(HorizontalBarChart,self).__init__(mpt_cmd_struct)

    # def __init__(...)

    def _SetFxnDesc(self):
        
        self.fxnDesc['DisplayName'] = 'Horizontal Bar Chart'
        self.fxnDesc['ShortDesc'] = 'Creates horizontal bar chart.'
        self.fxnDesc['ReturnType'] = 'Boolean'
        
        self.fxnDesc['ReqArgs'] = {
            'InFieldName':'Field Name',
            'BarLabels':'Any List',
            'OutFileName':'File Name',
            }
        self.fxnDesc['OptArgs'] = {
            'Metadata':'Any',
            'PrecursorFieldNames':['Field Name','Field Name List'],
            'Title':'Any',
            'XLabel':'Any',
            'XMax':'Float',
            'XMin':'Float',
            }
        
    # _SetFxnDesc(self):
 
    def DependencyNms(self):
        
        rtrn = self._ArgToList('InFieldName')
        rtrn += self._ArgToList('PrecursorFieldNames')
        return rtrn
    
    # def DependencyNms(self):
        
    def Exec(self,executedObjects):

        bar_lbls = self._ArgToList('BarLabels')
        in_fldnm = self.ValFromArgByNm('InFieldName')

        bar_vals = executedObjects[in_fldnm].ExecRslt().data.data.compressed()
        
        if len(bar_lbls) != len(bar_vals):
            raise Exception(
                (5*'{}\n').format(
                    '\n********************ERROR********************\n',
                    'Number of bar labels does not match number of values.',
                    '  InFieldName: {}'.format(in_fldnm),
                    '  BarLabels: {}'.format(' '.join(bar_lbls)),
                    'Full command:\n{}'.format(self.mptCmdStruct['rawCmdStr'])
                ),
            )

            bar_lbls = [bar_lbl.replace('_',' ') for bar_lbl in bar_lbls]

        # if len(bar_lbls) != len(in_fldnms):

        self.execRslt = self._HorizontalBarChart(
            bar_vals,
            bar_lbls,
            self.ValFromArgByNm('OutFileName'),
            )

        executedObjects[self.RsltNm()] = self

    # def Exec(self,executedObjects):

# class HorizontalBarChart(_GraphicsParent):
    
################################################################################
# Strided and Chunked statistics
################################################################################

'''
Chunking: taking a summary over each n entries in a dimension:
Take a dimension with a length of 12

Chunk 4 woud produce a result with a length of 3
0 1 2 3   4 5 6 7   8 9 10 11  # index within 1D arra
0 0 0 0   1 1 1 1   2 2  2  2  # contributes to chunk index

Striding takes uses every nth entry
Stride 4 would produce a result with a lenght of 4
0 1 2 3   4 5 6 7   8 9 10 11  # 
0 1 2 3   0 1 2 3   0 1  2  3  # contributes to stride index

For the 0-11 array above
Chunk 4 sum would be:
[6, 22, 38]

For the 0-11 array above
striding 4 sum would be:
[12, 15, 18, 21]

The general idea here is to split the axis being chunked or strided
into blocks of the chunk or stride size. Summarizing over the first
dimension of that split does chunking; summarizing over the second
does splitting.
'''

class _SummaryOverChunkOrStride(_NetCDFUtilParent):

    def __init__(
        self,
        mpt_cmd_struct=None,
        ):

        super(_SummaryOverChunkOrStride,self).__init__(mpt_cmd_struct)
        
    '''
    This should only be used as a parent class for classes
    that do a summary over one or more dimensions.
    '''

    # def __init__(...)

    def DependencyNms(self):
        return self._ArgToList('InFieldName') + \
          self._ArgToList('PrecursorFieldNames')

    def _reshape_arr_for_chunk_or_split(self,in_arr,axis,division_size):

        # Reshape the array. 
        new_shape = []
        new_shape[0:axis] = in_arr.shape[0:axis]
        new_shape += [in_arr.shape[axis]/division_size,division_size]
        new_shape += list(in_arr.shape[axis+1:])

        # return the summed strided values
        return in_arr.reshape(new_shape)

    # def _reshape_arr_for_chunk_or_split(in_arr,axis,division_size):

    def _Exec(self,executedObjects,summary_type,chunk_or_stride):

        parent_ncdimvar = executedObjects[self.ValFromArgByNm('InFieldName')].ExecRslt()
        dim_nm = self.ValFromArgByNm('DimensionName')
        
        if chunk_or_stride == 'Chunk':
            division_sz = int(self.ValFromArgByNm('ChunkSize'))
        elif chunk_or_stride == 'Stride':
            division_sz = int(self.ValFromArgByNm('StrideSize'))
        else:
            raise Exception(
                (2*'{}\n').format(
                    'Error: chunk_or_stride must be one of Chunk, Stride.',
                    '  Result name: {}'.format(self.RsltNm()),
                    )
                )
            
        # Insure specified dimension exists
        if dim_nm not in parent_ncdimvar.data.dim_nms:
            raise Exception(
                (3*'{}\n').format(
                    'Error: Invalid dimension specified.',
                    '  Result name: {}'.format(self.RsltNm()),
                    '  Dimension name: {}'.format(dim_nm)
                    )
                )

        # Insure even division of axis by division_sz

        if parent_ncdimvar.dims[dim_nm].shape[0] % division_sz != 0:
            raise Exception(
                (4*'{}\n').format(
                    'Error: Dimension length must be evenly divisible by {}Size.'.format(chunk_or_stride),
                    '  Result name: {}'.format(self.RsltNm()),
                    '  Dimension name, length: {}  {}'.format(dim_nm,parent_ncdimvar.dims[dim_nm].shape[0]),
                    '  {}Size: {}'.format(chunk_or_stride,division_sz)
                    )
                )

        # Now get the view of the data array with the target dimension split
        dimension_axis = parent_ncdimvar.dim_nms.index(dim_nm)
        array_with_split_dim = self._reshape_arr_for_chunk_or_split(
            parent_ncdimvar.data.data,
            dimension_axis,
            division_sz
            )

        # Now do the chunk or stride summary
        # operation_axis is the axis of the array_with_split_dim on which to summarize
        # first dim is stride, second is chunk
        if chunk_or_stride == 'Chunk':
            operation_axis = dimension_axis+1
            remaining_axis = dimension_axis
        elif chunk_or_stride == 'Stride':
            operation_axis = dimension_axis
            remaining_axis = dimension_axis+1
        else:
            raise Exception('Programming error: invalid chunk_or_stride')

        if summary_type == 'sum':
            data_nparr = array_with_split_dim.sum(axis=operation_axis)
        elif summary_type == 'mean':
            data_nparr = array_with_split_dim.mean(axis=operation_axis)
        elif summary_type == 'max':
            data_nparr = array_with_split_dim.max(axis=operation_axis)
        elif summary_type == 'min':
            data_nparr = array_with_split_dim.min(axis=operation_axis)
        elif summary_type == 'std':
            data_nparr = array_with_split_dim.std(axis=operation_axis)
        elif summary_type == 'var': # variance
            data_nparr = array_with_split_dim.var(axis=operation_axis)
        else:
            raise Exception(
                (3*'{}\n').format(
                    'Error: Illegal summary operation.',
                    '  Summary operation: {}'.format(summary_type)
                    )
                ) # raise Exception(...)

        # if summary_type == 'sum':...elif...else

        # Now let's build that new ncdimvar

        # First, take care of the dimension replacing the axis dimension

        new_dim_nm,new_dim_start,new_dim_step = self.ValFromArgByNm('NewDimension').split(':')
        new_dim_start = float(new_dim_start)
        new_dim_step = float(new_dim_step)

        new_dim_marr = np.ma.masked_array(
            [new_dim_start + x * new_dim_step for x in range(array_with_split_dim.shape[remaining_axis])],
            mask = False
            )

        new_dim_ncvar = mpncv.NCVar(
            data = new_dim_marr,
            name = new_dim_nm,
            dim_nms = new_dim_nm
            )

        # Then build the new set of dimensions

        key_to_replace = parent_ncdimvar.dims.keys()[dimension_axis]
        new_dims = OrderedDict([
            (new_dim_nm,new_dim_ncvar) if k == key_to_replace else (k,v) \
            for k,v in parent_ncdimvar.dims.items()
            ])

        # Then build the new ncvar

        new_ncvar = mpncv.NCVar(
            data = data_nparr,
            dim_nms = new_dims.keys(),
            name = self._CreateNCVarName()
            )
        
        # Then build the new ncdimvar

        new_ncdimvar = mpncv.NCDimensionedVar(
            dims = new_dims,
            data = new_ncvar,
            name = self._CreateNCVarName()
            )

        self.execRslt = new_ncdimvar
        self._SetMetadata()
        self.execRslt.set_metadata_val('NodeCommand',self.mptCmdStruct['rawCmdStr'])
        executedObjects[self.RsltNm()] = self
        
# class _SummaryOverChunkOrStride(_NetCDFUtilParent):

class MeanByChunk(_SummaryOverChunkOrStride):

    def __init__(
        self,
        mpt_cmd_struct=None,
        ):

        super(MeanByChunk,self).__init__(mpt_cmd_struct)
        
    # def __init__(...)

    def _SetFxnDesc(self):
        
        self.fxnDesc['DisplayName'] = 'MeanByChunk'
        self.fxnDesc['ShortDesc'] = 'Takes the mean of an NCVar by chunks over a dimension.'
        
        self.fxnDesc['ReturnType'] = 'NCVar',
        
        self.fxnDesc['ReqArgs'] = {
            'InFieldName':'Field Name',
            'DimensionName':'Field Name',
            'ChunkSize':'Integer',
            'NewDimension':'Tuple: Field Name:Float:Float'
            }
        self.fxnDesc['OptArgs'] = {
            'PrecursorFieldNames':['Field Name','Field Name List'],
            'Metadata':'Tuple List: Field Name:Any',
            'NewNCVarName':'Field Name',
            }
        
    # def _SetFxnDesc(self):

    def Exec(self,executedObjects):

        self._Exec(executedObjects,'mean','Chunk')
        
    # def Exec(self,executedObjects):

# class MeanByChunk(_NetCDFUtilParent):
    
class SumByChunk(_SummaryOverChunkOrStride):

    def __init__(
        self,
        mpt_cmd_struct=None,
        ):

        super(SumByChunk,self).__init__(mpt_cmd_struct)
        
    # def __init__(...)

    def _SetFxnDesc(self):
        
        self.fxnDesc['DisplayName'] = 'SumByChunk'
        self.fxnDesc['ShortDesc'] = 'Takes the sum of an NCVar by chunks over a dimension.'
        
        self.fxnDesc['ReturnType'] = 'NCVar, Float, or Integer'
        
        self.fxnDesc['ReqArgs'] = {
            'InFieldName':'Field Name',
            'DimensionName':'Field Name',
            'ChunkSize':'Integer',
            'NewDimension':'Tuple: Field Name:Float:Float'
            }
        self.fxnDesc['OptArgs'] = {
            'PrecursorFieldNames':['Field Name','Field Name List'],
            'Metadata':'Tuple List: Field Name:Any',
            'NewNCVarName':'Field Name',
            }
        
    # def _SetFxnDesc(self):

    def Exec(self,executedObjects):

        self._Exec(executedObjects,'sum','Chunk')
        
    # def Exec(self,executedObjects):

# class SumByChunk(_NetCDFUtilParent):
    
class MinByChunk(_SummaryOverChunkOrStride):

    def __init__(
        self,
        mpt_cmd_struct=None,
        ):

        super(MinByChunk,self).__init__(mpt_cmd_struct)
        
    # def __init__(...)

    def _SetFxnDesc(self):
        
        self.fxnDesc['DisplayName'] = 'MinByChunk'
        self.fxnDesc['ShortDesc'] = 'Takes the minimum of an NCVar by chunks over a dimension.'
        
        self.fxnDesc['ReturnType'] = 'NCVar, Float, or Integer'
        
        self.fxnDesc['ReqArgs'] = {
            'InFieldName':'Field Name',
            'DimensionName':'Field Name',
            'ChunkSize':'Integer',
            'NewDimension':'Tuple: Field Name:Float:Float'
            }
        self.fxnDesc['OptArgs'] = {
            'PrecursorFieldNames':['Field Name','Field Name List'],
            'Metadata':'Tuple List: Field Name:Any',
            'NewNCVarName':'Field Name',
            }
        
    # def _SetFxnDesc(self):

    def Exec(self,executedObjects):

        self._Exec(executedObjects,'min','Chunk')
        
    # def Exec(self,executedObjects):

# class MinByChunk(_NetCDFUtilParent):
    
class MaxByChunk(_SummaryOverChunkOrStride):

    def __init__(
        self,
        mpt_cmd_struct=None,
        ):

        super(MaxByChunk,self).__init__(mpt_cmd_struct)
        
    # def __init__(...)

    def _SetFxnDesc(self):
        
        self.fxnDesc['DisplayName'] = 'MaxByChunk'
        self.fxnDesc['ShortDesc'] = 'Takes the maximum of an NCVar by chunks over a dimension.'
        
        self.fxnDesc['ReturnType'] = 'NCVar, Float, or Integer'
        
        self.fxnDesc['ReqArgs'] = {
            'InFieldName':'Field Name',
            'DimensionName':'Field Name',
            'ChunkSize':'Integer',
            'NewDimension':'Tuple: Field Name:Float:Float'
            }
        self.fxnDesc['OptArgs'] = {
            'PrecursorFieldNames':['Field Name','Field Name List'],
            'Metadata':'Tuple List: Field Name:Any',
            'NewNCVarName':'Field Name',
            }
        
    # def _SetFxnDesc(self):

    def Exec(self,executedObjects):

        self._Exec(executedObjects,'max','Chunk')
        
    # def Exec(self,executedObjects):

# class MaxByChunk(_NetCDFUtilParent):
    
class StdDevByChunk(_SummaryOverChunkOrStride):

    def __init__(
        self,
        mpt_cmd_struct=None,
        ):

        super(StdDevByChunk,self).__init__(mpt_cmd_struct)
        
    # def __init__(...)

    def _SetFxnDesc(self):
        
        self.fxnDesc['DisplayName'] = 'StdDevByChunk'
        self.fxnDesc['ShortDesc'] = 'Takes the standard deviation of an NCVar by chunks over a dimension.'
        
        self.fxnDesc['ReturnType'] = 'NCVar, Float, or Integer'
        
        self.fxnDesc['ReqArgs'] = {
            'InFieldName':'Field Name',
            'DimensionName':'Field Name',
            'ChunkSize':'Integer',
            'NewDimension':'Tuple: Field Name:Float:Float'
            }
        self.fxnDesc['OptArgs'] = {
            'PrecursorFieldNames':['Field Name','Field Name List'],
            'Metadata':'Tuple List: Field Name:Any',
            'NewNCVarName':'Field Name',
            }
        
    # def _SetFxnDesc(self):

    def Exec(self,executedObjects):

        self._Exec(executedObjects,'std','Chunk')
        
    # def Exec(self,executedObjects):

# class StdDevByChunk(_NetCDFUtilParent):
    
class VarianceByChunk(_SummaryOverChunkOrStride):

    def __init__(
        self,
        mpt_cmd_struct=None,
        ):

        super(VarianceByChunk,self).__init__(mpt_cmd_struct)
        
    # def __init__(...)

    def _SetFxnDesc(self):
        
        self.fxnDesc['DisplayName'] = 'VarianceByChunk'
        self.fxnDesc['ShortDesc'] = 'Takes the variance of an NCVar by chunks over a dimension.'
        
        self.fxnDesc['ReturnType'] = 'NCVar, Float, or Integer'
        
        self.fxnDesc['ReqArgs'] = {
            'InFieldName':'Field Name',
            'DimensionName':'Field Name',
            'ChunkSize':'Integer',
            'NewDimension':'Tuple: Field Name:Float:Float'
            }
        self.fxnDesc['OptArgs'] = {
            'PrecursorFieldNames':['Field Name','Field Name List'],
            'Metadata':'Tuple List: Field Name:Any',
            'NewNCVarName':'Field Name',
            }
        
    # def _SetFxnDesc(self):

    def Exec(self,executedObjects):

        self._Exec(executedObjects,'var','Chunk')
        
    # def Exec(self,executedObjects):

# class VarianceByChunk(_NetCDFUtilParent):

class MeanByStride(_SummaryOverChunkOrStride):

    def __init__(
        self,
        mpt_cmd_struct=None,
        ):

        super(MeanByStride,self).__init__(mpt_cmd_struct)
        
    # def __init__(...)

    def _SetFxnDesc(self):
        
        self.fxnDesc['DisplayName'] = 'MeanByStride'
        self.fxnDesc['ShortDesc'] = 'Takes the mean of an NCVar by strides over a dimension.'
        
        self.fxnDesc['ReturnType'] = 'NCVar',
        
        self.fxnDesc['ReqArgs'] = {
            'InFieldName':'Field Name',
            'DimensionName':'Field Name',
            'StrideSize':'Integer',
            'NewDimension':'Tuple: Field Name:Float:Float'
            }
        self.fxnDesc['OptArgs'] = {
            'PrecursorFieldNames':['Field Name','Field Name List'],
            'Metadata':'Tuple List: Field Name:Any',
            'NewNCVarName':'Field Name',
            }
        
    # def _SetFxnDesc(self):

    def Exec(self,executedObjects):

        self._Exec(executedObjects,'mean','Stride')
        
    # def Exec(self,executedObjects):

# class MeanByStride(_NetCDFUtilParent):
    
class SumByStride(_SummaryOverChunkOrStride):

    def __init__(
        self,
        mpt_cmd_struct=None,
        ):

        super(SumByStride,self).__init__(mpt_cmd_struct)
        
    # def __init__(...)

    def _SetFxnDesc(self):
        
        self.fxnDesc['DisplayName'] = 'SumByStride'
        self.fxnDesc['ShortDesc'] = 'Takes the sum of an NCVar by striding over a dimension.'
        
        self.fxnDesc['ReturnType'] = 'NCVar, Float, or Integer'
        
        self.fxnDesc['ReqArgs'] = {
            'InFieldName':'Field Name',
            'DimensionName':'Field Name',
            'StrideSize':'Integer',
            'NewDimension':'Tuple: Field Name:Float:Float'
            }
        self.fxnDesc['OptArgs'] = {
            'PrecursorFieldNames':['Field Name','Field Name List'],
            'Metadata':'Tuple List: Field Name:Any',
            'NewNCVarName':'Field Name',
            }
        
    # def _SetFxnDesc(self):

    def Exec(self,executedObjects):

        self._Exec(executedObjects,'sum','Stride')
        
    # def Exec(self,executedObjects):

# class SumByStride(_NetCDFUtilParent):
    
class MinByStride(_SummaryOverChunkOrStride):

    def __init__(
        self,
        mpt_cmd_struct=None,
        ):

        super(MinByStride,self).__init__(mpt_cmd_struct)
        
    # def __init__(...)

    def _SetFxnDesc(self):
        
        self.fxnDesc['DisplayName'] = 'MinByStride'
        self.fxnDesc['ShortDesc'] = 'Takes the minimum of an NCVar by striding over a dimension.'
        
        self.fxnDesc['ReturnType'] = 'NCVar, Float, or Integer'
        
        self.fxnDesc['ReqArgs'] = {
            'InFieldName':'Field Name',
            'DimensionName':'Field Name',
            'StrideSize':'Integer',
            'NewDimension':'Tuple: Field Name:Float:Float'
            }
        self.fxnDesc['OptArgs'] = {
            'PrecursorFieldNames':['Field Name','Field Name List'],
            'Metadata':'Tuple List: Field Name:Any',
            'NewNCVarName':'Field Name',
            }
        
    # def _SetFxnDesc(self):

    def Exec(self,executedObjects):

        self._Exec(executedObjects,'min','Stride')
        
    # def Exec(self,executedObjects):

# class MinByStride(_NetCDFUtilParent):
    
class MaxByStride(_SummaryOverChunkOrStride):

    def __init__(
        self,
        mpt_cmd_struct=None,
        ):

        super(MaxByStride,self).__init__(mpt_cmd_struct)
        
    # def __init__(...)

    def _SetFxnDesc(self):
        
        self.fxnDesc['DisplayName'] = 'MaxByStride'
        self.fxnDesc['ShortDesc'] = 'Takes the maximum of an NCVar by striding over a dimension.'
        
        self.fxnDesc['ReturnType'] = 'NCVar, Float, or Integer'
        
        self.fxnDesc['ReqArgs'] = {
            'InFieldName':'Field Name',
            'DimensionName':'Field Name',
            'StrideSize':'Integer',
            'NewDimension':'Tuple: Field Name:Float:Float'
            }
        self.fxnDesc['OptArgs'] = {
            'PrecursorFieldNames':['Field Name','Field Name List'],
            'Metadata':'Tuple List: Field Name:Any',
            'NewNCVarName':'Field Name',
            }
        
    # def _SetFxnDesc(self):

    def Exec(self,executedObjects):

        self._Exec(executedObjects,'max','Stride')
        
    # def Exec(self,executedObjects):

# class MaxByStride(_NetCDFUtilParent):
    
class StdDevByStride(_SummaryOverChunkOrStride):

    def __init__(
        self,
        mpt_cmd_struct=None,
        ):

        super(StdDevByStride,self).__init__(mpt_cmd_struct)
        
    # def __init__(...)

    def _SetFxnDesc(self):
        
        self.fxnDesc['DisplayName'] = 'StdDevByStride'
        self.fxnDesc['ShortDesc'] = 'Takes the standard deviation of an NCVar by striding over a dimension.'
        
        self.fxnDesc['ReturnType'] = 'NCVar, Float, or Integer'
        
        self.fxnDesc['ReqArgs'] = {
            'InFieldName':'Field Name',
            'DimensionName':'Field Name',
            'StrideSize':'Integer',
            'NewDimension':'Tuple: Field Name:Float:Float'
            }
        self.fxnDesc['OptArgs'] = {
            'PrecursorFieldNames':['Field Name','Field Name List'],
            'Metadata':'Tuple List: Field Name:Any',
            'NewNCVarName':'Field Name',
            }
        
    # def _SetFxnDesc(self):

    def Exec(self,executedObjects):

        self._Exec(executedObjects,'std','Stride')
        
    # def Exec(self,executedObjects):

# class StdDevByStride(_NetCDFUtilParent):
    
class VarianceByStride(_SummaryOverChunkOrStride):

    def __init__(
        self,
        mpt_cmd_struct=None,
        ):

        super(VarianceByStride,self).__init__(mpt_cmd_struct)
        
    # def __init__(...)

    def _SetFxnDesc(self):
        
        self.fxnDesc['DisplayName'] = 'VarianceByStride'
        self.fxnDesc['ShortDesc'] = 'Takes the variance of an NCVar by striding over a dimension.'
        
        self.fxnDesc['ReturnType'] = 'NCVar, Float, or Integer'
        
        self.fxnDesc['ReqArgs'] = {
            'InFieldName':'Field Name',
            'DimensionName':'Field Name',
            'StrideSize':'Integer',
            'NewDimension':'Tuple: Field Name:Float:Float'
            }
        self.fxnDesc['OptArgs'] = {
            'PrecursorFieldNames':['Field Name','Field Name List'],
            'Metadata':'Tuple List: Field Name:Any',
            'NewNCVarName':'Field Name',
            }
        
    # def _SetFxnDesc(self):

    def Exec(self,executedObjects):

        self._Exec(executedObjects,'var','Stride')
        
    # def Exec(self,executedObjects):

# class VarianceByStride(_NetCDFUtilParent):

################################################################################
# Special functions
################################################################################

class MakeArrayWithValue(_NetCDFUtilParent):
    
    def __init__(
        self,
        mpt_cmd_struct=None,
        ):

        super(MakeArrayWithValue,self).__init__(mpt_cmd_struct)

    # def __init__(...)

    def _SetFxnDesc(self):
        
        self.fxnDesc['DisplayName'] = 'MakeArrayWithValue'
        self.fxnDesc['ShortDesc'] = 'Makes an array with all cells having specified values.'
        
        self.fxnDesc['ReturnType'] = 'NCVar'
        
        self.fxnDesc['ReqArgs'] = {
            'InFieldName':'Field Name',
            'Value':'Float',
            }
        self.fxnDesc['OptArgs'] = {
            'PrecursorFieldNames':['Field Name','Field Name List'],
            'Metadata':'Tuple List: Field Name:Any',
            'NewNCVarName':'Field Name',
            }
        
    # def _SetFxnDesc(self):
    
    def DependencyNms(self):
        return self._ArgToList('InFieldName') + \
          self._ArgToList('PrecursorFieldNames')

    def Exec(self,executedObjects):

        new_ncdimvar = cp.deepcopy(executedObjects[self.ValFromArgByNm('InFieldName')].ExecRslt())
        new_ncdimvar.data.data.data[:] = self.ValFromArgByNm('Value') 
        new_ncdimvar.name = self._CreateNCVarName()
            
        self.execRslt = new_ncdimvar
        self._SetMetadata()
        self.execRslt.set_metadata_val('NodeCommand',self.mptCmdStruct['rawCmdStr'])
        executedObjects[self.RsltNm()] = self
        
    # def _Exec(self,executedObjects):
     
# class MakeArrayWithValue(_NetCDFUtilParent):

class MakeNumberedCellsArray(_NetCDFUtilParent):
    
    def __init__(
        self,
        mpt_cmd_struct=None,
        ):

        super(MakeNumberedCellsArray,self).__init__(mpt_cmd_struct)

    # def __init__(...)

    def _SetFxnDesc(self):
        
        self.fxnDesc['DisplayName'] = 'MakeNumberedCellsArray'
        self.fxnDesc['ShortDesc'] = 'Makes an array of integers with cells numbered starting at 0.'
        
        self.fxnDesc['ReturnType'] = 'NCVar'
        
        self.fxnDesc['ReqArgs'] = {
            'InFieldName':'Field Name',
            }
        self.fxnDesc['OptArgs'] = {
            'PrecursorFieldNames':['Field Name','Field Name List'],
            'Metadata':'Tuple List: Field Name:Any',
            'NewNCVarName':'Field Name',
            }
        
    # def _SetFxnDesc(self):

    def Exec(self,executedObjects):

        in_ncdimvar = executedObjects[self.ArgByNm('InFieldName')].ExecRslt()
        
        rslt_arr = cp.deepcopy(in_ncdimvar.data.data.data)
        rslt_arr_ravel = rslt_arr.ravel()
        rslt_arr_ravel[:] = np.arange(len(rslt_arr_ravel)) # uses same data array as rslt_arr

        self.execRslt = mpncv.NCDimensionedVar(
            name = self._CreateNCVarName(),
            dims = cp.deepcopy(in_ncdimvar.dims),
            data = mpncv.NCVar(
                data = nc.ma.masked_array(rslt_arr,cp.deepcopy(in_ncdimvar.data.data.mask)),
                dim_nms = in_ncdimvar.dims.keys()
                ),
            )
        
        self._SetMetadata()
        self.execRslt.set_metadata_val('NodeCommand',self.mptCmdStruct['rawCmdStr'])
        executedObjects[self.RsltNm()] = self
        
    # def _Exec(self,executedObjects):
     
# class MakeNumberedCellsArray(_NetCDFUtilParent):

class ReclassifyByValueRanges(_NetCDFUtilParent):

    '''
    Reclassifies the input data based on value ranges.
    e.g. a tuple of 4:27:3 would cause any cell with
    4 <= value <= 27 to have value 3 in the resulting
    array.
    '''
    
    def __init__(
        self,
        mpt_cmd_struct=None,
        ):

        super(ReclassifyByValueRanges,self).__init__(mpt_cmd_struct)
        
    # def __init__(...)

    def _SetFxnDesc(self):
        
        self.fxnDesc['DisplayName'] = 'ReclassifyByValueRanges'
        self.fxnDesc['ShortDesc'] = 'Reclassifies input array using value ranges. Tuple is low value:high value:new value.'
        
        self.fxnDesc['ReturnType'] = 'NCVar'
        
        self.fxnDesc['ReqArgs'] = {
            'InFieldName':'Field Name',
            'ReclassSpecifications':'Tuple List: Integer:Integer:Integer'
            }
        self.fxnDesc['OptArgs'] = {
            'PrecursorFieldNames':['Field Name','Field Name List'],
            'Metadata':'Tuple List: Field Name:Any',
            'NewNCVarName':'Field Name',
            'DefaultValue':'Integer'
            }
        
    # def _SetFxnDesc(self):

    def DependencyNms(self):
        return self._ArgToList('InFieldName') + \
          self._ArgToList('PrecursorFieldNames')

    def Exec(self,executedObjects):

        in_ncdimvar = executedObjects[self.ValFromArgByNm('InFieldName')].ExecRslt()
        in_marr = in_ncdimvar.data.data

        rslt_marr = cp.deepcopy(in_marr)
        dflt_val = self.ValFromArgByNm('DefaultValue')
        if dflt_val is not None:
            rslt_marr[~rslt_marr.mask] = dflt_val

        for tup in self._ArgToList('ReclassSpecifications'):

            class_val1,class_val2,new_val =''.join(tup.split()).split(':')
            class_min_val = min(float(class_val1),float(class_val2))
            class_max_val = max(float(class_val1),float(class_val2))
            new_val = float(new_val)

            rslt_marr = np.ma.where(
                np.ma.logical_and(in_marr >= class_min_val, in_marr <= class_max_val),
                new_val,
                rslt_marr
                )

        # Create new dimensioned array
        dims = cp.deepcopy(executedObjects[self.ArgByNm('InFieldName')].ExecRslt().dims)

        self.execRslt = mpncv.NCDimensionedVar(
            name = self._CreateNCVarName(),
            dims = dims,
            data = mpncv.NCVar(
                data=rslt_marr,
                dim_nms=dims.keys()
                )
            )

        self._SetMetadata()
        self.execRslt.set_metadata_val('NodeCommand',self.mptCmdStruct['rawCmdStr'])
        executedObjects[self.RsltNm()] = self
        
    # def Exec(self,executedObjects):

# class ReclassifyByValueRanges(_NetCDFUtilParent):

class QuantifyByCrossReference(_NetCDFUtilParent):

    # This MPilot function is for looking up a table value associated with
    # two value classes. For instance, a vegetation difference table might
    # look like this:
    #     1  2  3  4  5
    # 1   0  1  2  3  3
    # 2   1  0  1  2  3
    # 3   2  1  0  1  2
    # 4   3  2  1  0  1
    # 5   3  3  2  1  0
    #
    # So the difference between veg type 2 and veg type 4 yields a value
    # of 2.
    #
    # The reason for creating this is so that vegetation types from different
    # time periods and scenarios can be compared and difference levels
    # generated. I'm sure there are other uses
    #
    # 2018.07.18 - tjs
    # created

    def __init__(
        self,
        mpt_cmd_struct=None,
        ):

        super(QuantifyByCrossReference,self).__init__(mpt_cmd_struct)
        
    # def __init__(...)
    
    def _SetFxnDesc(self):
        
        self.fxnDesc['DisplayName'] = 'Classify by Cross Reference'
        self.fxnDesc['ShortDesc'] = 'Using a 2d table gets a lookup value based on 2 inputs'
        self.fxnDesc['ReturnType'] = 'NCVar'
        
        self.fxnDesc['ReqArgs'] = {
            'InFieldNames':'Field Name List',
            'CSVTableFileName':'File Name',
            }
        self.fxnDesc['OptArgs'] = {
            'PrecursorFieldNames':['Field Name','Field Name List'],
            'Metadata':'Tuple List: Field Name:Any',
            'NewNCVarName':'Field Name',
            }
                    
    # def _SetFxnDesc(self):

    def DependencyNms(self):
        return self._ArgToList('InFieldNames') + self._ArgToList('PrecursorFieldNames')

    def _classify_by_x_reference(
        row_val,
        col_val,
        row_lu_lst,
        col_lu_lst,
        lu_val_arr
        ):

        nan_val = float('nan')

        try:
            row_ndx = row_lu_lst.index(row_val)
        except ValueError:
            return nan_val

        try:
            col_ndx = col_lu_lst.index(col_val)
        except ValueError:
            return nan_val

        return lu_val_arr[row_ndx,col_ndx]

    # def _classify_by_x_reference(...)

    def Exec(self,executedObjects):

        fldnm_list = self._ArgToList('InFieldNames')

        if len(fldnm_list) != 2:
            raise Exception(
                (3*'{}\n').format(
                    'Error. Must have exactly 2 In Field Names',
                    '  Result Name, Field Name List: {}  [{}]'.format(
                        self.RsltNm(),
                        ','.join(fldnm_list)
                        ),
                    'Full command:\n{}'.format(self.mptCmdStruct['rawCmdStr'])
                    )
                ) # raise Exception(...)

        # if len(fldnm_list) < 2:
        
        row_ncdimvar = executedObjects[fldnm_list[0]].ExecRslt()
        col_ncdimvar = executedObjects[fldnm_list[1]].ExecRslt()

        # Make congruent arrays used in the lookup

        if self._NCDimArrsAreCongruent(row_ncdimvar,col_ncdimvar):

            delete_tmp_row_ncdimvar = delete_tmp_col_ncdimvar = False
            tmp_row_ncdimvar = row_ncdimvar
            tmp_col_ncdimvar = col_ncdimvar
                
        else:

            # try to expand row_ncdimvar to match
            tmp_row_ncdimvar = self._MakeCongruentNCDimArr(row_ncdimvar,col_ncdimvar)
            if tmp_row_ncdimvar is not None:
                tmp_col_ncdimvar = col_ncdimvar
                delete_tmp_row_ncdimvar = True
                delete_tmp_col_ncdimvar = False
            else:
                tmp_col_ncdimvar = self._MakeCongruentNCDimArr(col_ncdimvar,row_ncdimvar)
                if tmp_col_ncdimvar is not None:
                    tmp_row_ncdimvar = row_ncdimvar
                    delete_tmp_row_ncdimvar = False
                    delete_tmp_col_ncdimvar = True
                else:
                    raise Exception(
                        (5*'{}\n').format(
                            'Error. Cannot make arrays congruent for operation.',
                            '  Result name: {}'.format(self.RsltNm()),
                            '  First array: Field Name, NCVar name, dim names, shape: {}  {}  {}  {}'.format(
                                fldnm_list[0],
                                executedObjects[fldnm_list[0]].ExecRslt().name,
                                executedObjects[fldnm_list[0]].ExecRslt().dims.keys(),
                                executedObjects[fldnm_list[0]].ExecRslt().data.data.shape,
                                ),
                            '  Second array: Field Name, NCVar name, dim names, shape: {}  {}  {}  {}'.format(
                                fldnm_list[1],
                                executedObjects[fldnm_list[1]].ExecRslt().name,
                                executedObjects[fldnm_list[1]].ExecRslt().dims.keys(),
                                executedObjects[fldnm_list[1]].ExecRslt().data.data.shape,
                                ),
                            'Full command:\n{}'.format(self.mptCmdStruct['rawCmdStr'])
                            )
                        ) # raise Exception(...)
                        
            # if self._NCDimArrsAreCongruent(ncdimvar_for_op,new_ncdimvar):...else...

        if not os.path.isfile(self.ArgByNm('CSVTableFileName')):
            raise Exception(
                4*'{}'.format(
                    '\n********************ERROR********************\n',
                    'CSV table file does not exist: {}\n'.format(self.ArgByNm('CSVTableFileName')),
                    'Script File: {}  Line number: {}\n'.format(
                        self.mptCmdStruct['cmdFileNm'],
                        self.mptCmdStruct['lineNo']
                        ),
                    'Full command:\n{}\n'.format(self.mptCmdStruct['rawCmdStr'])
                ),
            )
        # if not os.path.isfile(self.ArgByNm('CSVTableFileName'),'r'):

        # Read in the csv table
        raw_lookup = np.genfromtxt(self.ArgByNm('CSVTableFileName'),delimiter=',')
        lu_col_keys = raw_lookup[0,1:].tolist()
        lu_row_keys = raw_lookup[1:,0].tolist()
        lu_values = raw_lookup[1:,1:]

        # The data arrays corresponding to the row heads and col heads in the lu array
        row_key_data_arr = tmp_row_ncdimvar.data.data
        col_key_data_arr = tmp_col_ncdimvar.data.data

        rslt_arr = np.full_like(
            row_key_data_arr,
            np.ma.default_fill_value (np.ma.default_fill_value(np.dtype(float))),
            dtype = np.float
            )

        # Go through all the possible combinations of row/col classes and assign lu value
        # done in place for efficiency
        for lu_row_key in lu_row_keys:
            for lu_col_key in lu_col_keys:

                # get the lookup value associated with the keys
                lookup_val = lu_values[lu_row_keys.index(lu_row_key),lu_col_keys.index(lu_col_key)]

                # fill in entries associated with the row/col combination
                np.place (
                    rslt_arr,
                    np.logical_and (row_key_data_arr == lu_row_key, col_key_data_arr == lu_col_key),
                    float(lookup_val)
                    )
                
            # for lu_col_key in lu_col_keys:
        # for lu_row_key in lu_row_keys:

        rslt_arr = np.ma.masked_where(rslt_arr == np.ma.default_fill_value(np.dtype(float)), rslt_arr)

        # The data NCVar:
        self.execRslt = new_ncdimvar = mpncv.NCDimensionedVar(
            dims = cp.deepcopy(tmp_row_ncdimvar.dims),
            name = self._CreateNCVarName(),
            data = mpncv.NCVar(
                data = rslt_arr,
                dim_nms = tmp_row_ncdimvar.dims.keys(),
                )
            )

        self._SetMetadata()
        self.execRslt.set_metadata_val('NodeCommand',self.mptCmdStruct['rawCmdStr'])
        self.execRslt.set_metadata_val(
            'Description',
            'Lookup results for {} x {} using table {}'.format(
                row_ncdimvar.name,
                col_ncdimvar.name,
                self.ArgByNm('CSVTableFileName')
                )
            )

        executedObjects[self.RsltNm()] = self

        # cleanup
        if delete_tmp_row_ncdimvar == True:
            del tmp_row_ncdimvar
        if delete_tmp_col_ncdimvar == True:
            del tmp_col_ncdimvar

    # def Exec(self,executedObjects):

# class QuantifyByCrossReference(_NetCDFUtilParent):
    
################################################################################

class FlipDimensions(_NetCDFUtilParent):

    '''
    Reverses the order of dimension values in an ncvar
    '''
    def __init__(
        self,
        mpt_cmd_struct=None,
        ):

        super(FlipDimensions,self).__init__(mpt_cmd_struct)
        
    # def __init__(...)

    def _SetFxnDesc(self):
        
        self.fxnDesc['DisplayName'] = 'FlipDimensions'
        self.fxnDesc['ShortDesc'] = 'Reverses the order of dimension values in a variable'
        
        self.fxnDesc['ReturnType'] = 'NCVar'
        
        self.fxnDesc['ReqArgs'] = {
            'InFieldName':'Field Name',
            'DimensionNames':['Field Name','Field Name List']
            }
        self.fxnDesc['OptArgs'] = {
            'PrecursorFieldNames':['Field Name','Field Name List'],
            'Metadata':'Tuple List: Field Name:Any',
            'NewNCVarName':'Field Name',
            }
        
    # def _SetFxnDesc(self):

    def DependencyNms(self):
        return self._ArgToList('InFieldName') + \
          self._ArgToList('PrecursorFieldNames')

    def Exec(self,executedObjects):
        
        new_ncdimvar = cp.deepcopy(executedObjects[self.ValFromArgByNm('InFieldName')].ExecRslt())
        new_ncdimvar.name = self._CreateNCVarName()
        
        for dim_nm in self._ArgToList('DimensionNames'):
            
            if dim_nm not in new_ncdimvar.dim_nms:
                raise Exception(
                    '{}{}{}'.format(
                        'Illegal dimension specified: {}'.format(dim_nm),
                        'File: {}  Line number: {}\n'.format(
                            self.mptCmdStruct['cmdFileNm'],self.mptCmdStruct['lineNo']
                            ),
                        'Full command:\n{}\n'.format(self.mptCmdStruct['rawCmdStr'])
                        ),
                    )

            axis = new_ncdimvar.dim_nms.index(dim_nm)

            new_ncdimvar.data.data = np.flip(new_ncdimvar.data.data,axis)
            new_ncdimvar.dims[dim_nm].data = np.flip(new_ncdimvar.dims[dim_nm].data,0)
            
        # for dim_nm in self._ArgToList('DimensionNames'):
        
        self.execRslt = new_ncdimvar
        self._SetMetadata()
        self.execRslt.set_metadata_val('NodeCommand',self.mptCmdStruct['rawCmdStr'])
        executedObjects[self.RsltNm()] = self
        
    # def Exec(self,executedObjects):

# class FlipDimensions(_NetCDFUtilParent):

class FlattenNCVar(_NetCDFUtilParent):

    '''
    Turns a variable into a 1d array
    '''
    def __init__(
        self,
        mpt_cmd_struct=None,
        ):

        super(FlattenNCVar,self).__init__(mpt_cmd_struct)
        
    # def __init__(...)

    def _SetFxnDesc(self):
        
        self.fxnDesc['DisplayName'] = 'FlattenNCVar'
        self.fxnDesc['ShortDesc'] = 'Converts a variable array into a 1d array'
        
        self.fxnDesc['ReturnType'] = 'NCVar'
        
        self.fxnDesc['ReqArgs'] = {
            'InFieldName':'Field Name',
            'NewDimension':'Tuple: Field Name:Float:Float'
            }
        self.fxnDesc['OptArgs'] = {
            'DeleteMissing':'Boolean',
            'PrecursorFieldNames':['Field Name','Field Name List'],
            'Metadata':'Tuple List: Field Name:Any',
            'NewNCVarName':'Field Name',
            }
        
    # def _SetFxnDesc(self):

    def DependencyNms(self):
        return self._ArgToList('InFieldName') + \
          self._ArgToList('PrecursorFieldNames')

    def Exec(self,executedObjects):

        # Make the ncdimvar array
        if self.ValFromArgByNm('DeleteMissing') == True:
            rslt_marr = executedObjects[self.ValFromArgByNm('InFieldName')].ExecRslt().data.data.compressed()
        else:
            rslt_marr = executedObjects[self.ValFromArgByNm('InFieldName')].ExecRslt().data.data.flatten()

        # Make the dimension array
        new_dim_nm,new_dim_start,new_dim_step = self.ValFromArgByNm('NewDimension').split(':')
        new_dim_start = float(new_dim_start)
        new_dim_step = float(new_dim_step)

        new_dim_marr = np.ma.masked_array(
            [new_dim_start + x * new_dim_step for x in range(rslt_marr.shape[0])],
            mask = False
            )

        new_dim_ncvar = mpncv.NCVar(
            data = new_dim_marr,
            name = new_dim_nm,
            dim_nms = new_dim_nm
            )

        # Then build the new ncdimvar

        new_ncdimvar = mpncv.NCDimensionedVar(
            dims = OrderedDict([(new_dim_nm,new_dim_ncvar)]),
            data = mpncv.NCVar(
                data = rslt_marr,
                dim_nms = [new_dim_nm],
                ),
            name = self._CreateNCVarName()
            )

        self.execRslt = new_ncdimvar
        self._SetMetadata()
        self.execRslt.set_metadata_val('NodeCommand',self.mptCmdStruct['rawCmdStr'])
        executedObjects[self.RsltNm()] = self
        
    # def Exec(self,executedObjects):

# class FlattenNCVar(_NetCDFUtilParent):

class ProjectToMatchTemplate(_NetCDFUtilParent):

    '''
    Projects (expands) one array to match dimensions of another
    '''
    def __init__(
        self,
        mpt_cmd_struct=None,
        ):

        super(ProjectToMatchTemplate,self).__init__(mpt_cmd_struct)
        
    # def __init__(...)

    def _SetFxnDesc(self):
        
        self.fxnDesc['DisplayName'] = 'ProjectDimsToMatchTemplate'
        self.fxnDesc['ShortDesc'] = 'Projects a FieldNameToProject to match those of TemplateFieldName'
        
        self.fxnDesc['ReturnType'] = 'NCVar'
        
        self.fxnDesc['ReqArgs'] = {
            'FieldNameToProject':'Field Name',
            'TemplateFieldName':'Field Name',
            }
        self.fxnDesc['OptArgs'] = {
            'PrecursorFieldNames':['Field Name','Field Name List'],
            'Metadata':'Tuple List: Field Name:Any',
            'NewNCVarName':'Field Name',
            }
        
    # def _SetFxnDesc(self):

    def DependencyNms(self):
        return self._ArgToList('FieldNameToProject') + \
          self._ArgToList('TemplateFieldName') + \
          self._ArgToList('PrecursorFieldNames')

    def Exec(self,executedObjects):

        to_project_ncdimvar = executedObjects[self.ValFromArgByNm('FieldNameToProject')].ExecRslt()
        tmplt_ncdimvar = executedObjects[self.ValFromArgByNm('TemplateFieldName')].ExecRslt()

        self.execRslt = self._MakeCongruentNCDimArr(to_project_ncdimvar,tmplt_ncdimvar)

        if self.execRslt is None:
            raise Exception(
                (5*'{}\n').format(
                    '\n********************ERROR********************\n',
                    'FieldNameToProject cannot be expanded to match TemplateFieldName',
                    '  FieldNameToProject dimensions: {}'.format(','.join(to_project_ncdimvar.dim_nms)),
                    '  TemplateFieldName dimensions: {}'.format(','.join(tmplt_ncdimvar.dim_nms)),
                    'Full command:\n{}'.format(self.mptCmdStruct['rawCmdStr'])
                ),
            )

        self.execRslt.name = self._CreateNCVarName()
        self._SetMetadata()
        self.execRslt.set_metadata_val('NodeCommand',self.mptCmdStruct['rawCmdStr'])
        executedObjects[self.RsltNm()] = self
        
    # def Exec(self,executedObjects):

# class ProjectToMatchTemplate(_NetCDFUtilParent):

################################################################################

class CopyMetadata(_NetCDFUtilParent):

    '''
    Copies metadata from one variable to another
    '''
    
    def __init__(
        self,
        mpt_cmd_struct=None,
        ):

        super(CopyMetadata,self).__init__(mpt_cmd_struct)
        
    # def __init__(...)

    def _SetFxnDesc(self):
        
        self.fxnDesc['DisplayName'] = 'CopyMetadata'
        self.fxnDesc['ShortDesc'] = 'Copies metadata from one field to another'
        
        self.fxnDesc['ReturnType'] = 'NCVar'
        
        self.fxnDesc['ReqArgs'] = {
            'InFieldName':'Field Name',
            'InMetadataFieldName':'Field Name',
            }
        self.fxnDesc['OptArgs'] = {
            'PrecursorFieldNames':['Field Name','Field Name List'],
            'Metadata':'Tuple List: Field Name:Any',
            'NewNCVarName':'Field Name',
            }
        
    # def _SetFxnDesc(self):

    def DependencyNms(self):
        return self._ArgToList('InFieldName') + \
          self._ArgToList('InMetadataFieldName') + \
          self._ArgToList('PrecursorFieldNames')

    def Exec(self,executedObjects):

        new_ncdimvar = cp.deepcopy(executedObjects[self.ValFromArgByNm('InFieldName')].ExecRslt())
        new_ncdimvar.name = self._CreateNCVarName()
        md_src_dict = executedObjects[self.ValFromArgByNm('InMetadataFieldName')].ExecRslt().metadata
        
        md_nms = self._ArgToList('MetaDataNames')
        if md_nms is None:
            new_ncdimvar.metadata = md_src_dict
        else:
            for md_nm in md_nms:
                if md_nm in md_src_dict:
                    new_ncdimvar.set_metadata_val(md_nm,md_src_dict[md_nm])
                else:
                    # ignore absent requested metadata items
                    pass
        # if md_nms is None:

        self.execRslt = new_ncdimvar
        self.execRslt.set_metadata_val('NodeCommand',self.mptCmdStruct['rawCmdStr'])
        executedObjects[self.RsltNm()] = self
        
    # def Exec(self,executedObjects):

# class CopyMetadata(_NetCDFUtilParent):

class SetMetadata(_NetCDFUtilParent):

    '''
    Sets metadata.
    '''
    
    def __init__(
        self,
        mpt_cmd_struct=None,
        ):

        super(SetMetadata,self).__init__(mpt_cmd_struct)
        
    # def __init__(...)

    def _SetFxnDesc(self):
        
        self.fxnDesc['DisplayName'] = 'SetMetadata'
        self.fxnDesc['ShortDesc'] = 'Sets metadata field value. Underscore becomes space in the value field of Name:Value tuple.'
        
        self.fxnDesc['ReturnType'] = 'NCVar'
        
        self.fxnDesc['ReqArgs'] = {
            'InFieldName':'Field Name',
            'MetadataDefinitions':'Tuple List: Any:Any',
            }
        self.fxnDesc['OptArgs'] = {
            'PrecursorFieldNames':['Field Name','Field Name List'],
            'Metadata':'Tuple List: Field Name:Any',
            'NewNCVarName':'Field Name',
            }
        
    # def _SetFxnDesc(self):

    def DependencyNms(self):
        return self._ArgToList('InFieldName') + \
          self._ArgToList('PrecursorFieldNames')

    def Exec(self,executedObjects):

        new_ncdimvar = cp.deepcopy(executedObjects[self.ValFromArgByNm('InFieldName')].ExecRslt())
        new_ncdimvar.name = self._CreateNCVarName()

        for md_def in self._ArgToList('MetaDataDefinitions'):
            key,val = md_def.split(':',2)
            val = val.replace('_',' ')
            new_ncdimvar.set_metadata_val(key,val)
                    
        self.execRslt = new_ncdimvar
        self.execRslt.set_metadata_val('NodeCommand',self.mptCmdStruct['rawCmdStr'])
        executedObjects[self.RsltNm()] = self
        
    # def Exec(self,executedObjects):

# class SetMetadata(_NetCDFUtilParent):

################################################################################

class PrintVariables(_NetCDFUtilParent):

    def __init__(
        self,
        mpt_cmd_struct=None,
        ):

        super(PrintVariables,self).__init__(mpt_cmd_struct)
        
    # def __init__(...)

    def _SetFxnDesc(self):
        
        self.fxnDesc['DisplayName'] = 'PrintVariables'
        self.fxnDesc['ShortDesc'] = 'Prints some NCVars.'
        self.fxnDesc['ReturnType'] = 'Boolean'
        
        self.fxnDesc['ReqArgs'] = {
            'InFieldNames':['Field Name','Field Name List'],
            }
        self.fxnDesc['OptArgs'] = {
            'PrecursorFieldNames':['Field Name','Field Name List'],
            'AsCSV':'Boolean',
            'SuppressCSVHeaders':'Boolean'
            }
        
    # def _SetFxnDesc(self):

    def DependencyNms(self):
        rtrn = self._ArgToList('InFieldNames') + \
        self._ArgToList('PrecursorFieldNames')
        return rtrn

    def Exec(self,executedObjects):
        
        print_objs = OrderedDict()
        for in_fldnm in self._ArgToList('InFieldNames'):
            if isinstance(executedObjects[in_fldnm].ExecRslt(), mpncv.NCDimensionedVar):
                print_objs[in_fldnm] = executedObjects[in_fldnm].ExecRslt().data.data                
            else:
                print_objs[in_fldnm] = executedObjects[in_fldnm].ExecRslt()

        if self.ValFromArgByNm('AsCSV') == True:

            in_fldnms = self._ArgToList('InFieldNames')
            
            # Check that all variables are 1d with the same index
            dim_nms_for_checking = executedObjects[in_fldnms[0]].ExecRslt().dim_nms
            if len(dim_nms_for_checking) != 1:
                raise Exception(
                    (4*'{}\n').format(
                        '\n********************ERROR********************',
                        'Writing variables as CSV only valid for 1d variables with matching index.',
                        '  Variable: {} has indexes {}'.format(in_fldnms[0],dim_nms_for_checking),
                        'Full command:\n{}'.format(self.mptCmdStruct['rawCmdStr'])
                    ),
                )
            
            for in_fldnm in in_fldnms:
                if executedObjects[in_fldnm].ExecRslt().dim_nms != dim_nms_for_checking:
                    
                    raise Exception(
                        (4*'{}\n').format(
                            '\n********************ERROR********************',
                            'Writing variables as CSV only valid for 1d variables with matching index.',
                            '  Variable: {} has index {}'.format(in_fldnms[0],dim_nms_for_checking),
                            '  Variable: {} has index(es) {}'.format(
                                in_fldnm,dim_nms_for_checking,
                                executedObjects[in_fldnm].ExecRslt().dim_nms
                                ),
                            'Full command:\n{}'.format(self.mptCmdStruct['rawCmdStr'])
                        ),
                    )
            
            # for in_fldnm in in_fldnms:

            dimension_ncvar = executedObjects[in_fldnms[0]].ExecRslt().dims.values()[0]

            if self.ValFromArgByNm('SuppressCSVHeaders') != True:
                print '{},{}'.format(dimension_ncvar.name,','.join(in_fldnms))

            for ndx in range(len(dimension_ncvar.data.data)):
                outline = '{}'.format(dimension_ncvar.data.data[ndx])
                
                for in_fldnm in in_fldnms:
                    outline = '{},{}'.format(outline,executedObjects[in_fldnm].ExecRslt().data.data[ndx])

                print outline

        else:
            for in_fldnm in self._ArgToList('InFieldNames'):
            
                print 'InFieldName: {}'.format(in_fldnm)
                print print_objs[in_fldnm]
        
        self.execRslt = True
        executedObjects[self.RsltNm] = self
            
    # def Exec(self,executedObjects):
        
# class PrintVariables(_NetCDFUtilParent):

class PrintString(_NetCDFUtilParent):

    def __init__(
        self,
        mpt_cmd_struct=None,
        ):

        super(PrintString,self).__init__(mpt_cmd_struct)
        
    # def __init__(...)

    def _SetFxnDesc(self):
        
        self.fxnDesc['DisplayName'] = 'PrintString'
        self.fxnDesc['ShortDesc'] = 'Prints a string. Use an underscore to represent a space.'
        self.fxnDesc['ReturnType'] = 'Boolean'
        
        self.fxnDesc['ReqArgs'] = {
            'String':'String'
            }
        self.fxnDesc['OptArgs'] = {
            'PrecursorFieldNames':['Field Name','Field Name List']
            }
        
    # def _SetFxnDesc(self):

    def DependencyNms(self):
        rtrn = self._ArgToList('PrecursorFieldNames')
        return rtrn

    def Exec(self,executedObjects):
        print '{}\n'.format(self.ArgByNm('String').replace('_',' '))

        self.execRslt = True
        executedObjects[self.RsltNm] = self

    # def Exec(self,executedObjects):
        
# class PrintVariables(_NetCDFUtilParent):

################################################################################
# Fuzzy Logic Functions
################################################################################

class _FuzzyParent(_NetCDFUtilParent):
    
    fuzzyMax = 1.0
    #fuzzyMin = -1.0
    fuzzyMin = 0.0

    def __init__(
        self,
        mpt_cmd_struct=None,
        ):

        super(_FuzzyParent,self).__init__(mpt_cmd_struct)
        
    def _InsureFuzzy(self,arr):
        
        arr[np.ma.where(arr > self.fuzzyMax)] = self.fuzzyMax
        arr[np.ma.where(arr < self.fuzzyMin)] = self.fuzzyMin
                
    # def _InsureFuzzy(self,arr):
       
# class _FuzzyParent(_NetCDFUtilParent):

class _CvtToFuzzyCurve(_FuzzyParent):
   
    def __init__(self,mptCmdStruct=None):
        
        super(_CvtToFuzzyCurve, self).__init__(mptCmdStruct)

        self.raw_vals = None
        self.fuzzy_vals = None

    # def __init__(self,mptCmdStruct=None):

    def DependencyNms(self):
        
        rtrn = self._ArgToList('InFieldName')
        rtrn += self._ArgToList('PrecursorFieldNames')
        return rtrn
    
    # def DependencyNms(self):

    def _Exec(self,executedObjects):
        # executedObjects is a dictionary of executed MPilot function objects

        # This assumes self.raw_vals and self.fuzzy_vals have been set by the
        # descendent function. These two object variables will also be used
        # by GraphConversionCurves
        
        in_fld_nm = self.ArgByNm('InFieldName')
        in_NC_dim_var = executedObjects[in_fld_nm].ExecRslt()
        in_arr = executedObjects[in_fld_nm].ExecRslt().data.data

        rslt_arr = np.ma.empty(
            in_arr.shape,
            dtype=float
            )

        # Sort the raw values and associated fuzzy values
        self.raw_vals, self.fuzzy_vals = map(list,zip(*sorted(zip(self.raw_vals,self.fuzzy_vals))))
        
        # Set the fuzzy values corresponding to the raw values less
        # than the lowest raw value to the corresponding fuzzy value
        np.place(rslt_arr.data,in_arr <= self.raw_vals[0],self.fuzzy_vals[0])

        # Iterate over the line segments that approximate the curve
        # and assign fuzzy values.
        for ndx in range(1,len(self.raw_vals)):

            # Linear equation formula for line segment
            m = (self.fuzzy_vals[ndx] - self.fuzzy_vals[ndx-1]) / \
            (self.raw_vals[ndx] - self.raw_vals[ndx-1])
            
            b = self.fuzzy_vals[ndx-1] - m * self.raw_vals[ndx-1]

            whereNdxs = np.where(
                np.logical_and(
                    in_arr.data > self.raw_vals[ndx-1],
                    in_arr.data <= self.raw_vals[ndx]
                    )
                )

            rslt_arr.data[whereNdxs] = m * in_arr.data[whereNdxs] + b

        # Set the fuzzy values corresponding to the raw values greater
        # than the highest raw value to the corresponding fuzzy value
        np.place(rslt_arr.data,in_arr > self.raw_vals[-1],self.fuzzy_vals[-1])
        
        rslt_arr.mask = cp.deepcopy(in_arr.mask)

        # bring back values to fuzzy limits
        self._InsureFuzzy(rslt_arr)

        self.execRslt = mpncv.NCDimensionedVar(
            dims = cp.deepcopy(in_NC_dim_var.dims),
            name = self._CreateNCVarName(),
            data = mpncv.NCVar(
                data = rslt_arr,
                dim_nms = in_NC_dim_var.dims.keys(),
                )
            )
        
        self._SetMetadata()
        self.execRslt.set_metadata_val('NodeCommand',self.mptCmdStruct['rawCmdStr'])
        executedObjects[self.RsltNm()] = self
        
    # def _Exec(self,executedObjects):

# class _CvtToFuzzyCurve(_FuzzyParent):

class CvtToFuzzyZScores(_CvtToFuzzyCurve):
    
    def __init__(self,mptCmdStruct=None):
        
        super(CvtToFuzzyZScores, self).__init__(mptCmdStruct)

    # def __init__(self,mptCmdStruct=None):

    def _SetFxnDesc(self):
        # description of command for validation and info display

        self.fxnDesc['DisplayName'] = 'Convert to Fuzzy using Z Score'
        self.fxnDesc['ShortDesc'] = 'Converts input values into fuzzy values using linear interpolation based on Z Scores'
        self.fxnDesc['ReturnType'] = 'Fuzzy'

        self.fxnDesc['ReqArgs'] = {
            'InFieldName':'Field Name',
            'ZScores':['Float List'],
            'FuzzyValues':['Fuzzy Value List'],
            }
        self.fxnDesc['OptArgs'] = {
            'Metadata':'Tuple List: Field Name:Any',
            'NewNCVarName':'Field Name',            
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
        
        in_fld_nm = self.ArgByNm('InFieldName')
        in_NC_dim_var = executedObjects[in_fld_nm].ExecRslt()
        in_arr = executedObjects[in_fld_nm].ExecRslt().data.data

        self.fuzzy_vals = self.ValFromArgByNm('FuzzyValues')
        zscores = self.ValFromArgByNm('ZScores')

        if len(self.fuzzy_vals) != len(zscores):
            raise Exception(
                (3*'{}\n').format(
                    'Error: Number of zscores and raw values must match.',
                    '  Result name: {}'.format(self.RsltNm),
                    '  Number of ZScores, fuzzy values: {}'.format(len(zscores),len(self.fuzzy_vals))
                    )
                ) # raise Exception(...) 

        raw_mean = np.ma.mean(in_arr)
        raw_stddev = np.ma.std(in_arr)
                
        self.raw_vals = []
        for zscore in zscores:
            self.raw_vals.append(raw_mean + raw_stddev * zscore)
        
        self._Exec(executedObjects)

    # def Exec(self,executedObjects):
    
# class CvtToFuzzyZScore(_CvtToFuzzyCurve):

class CvtToFuzzyCat(_FuzzyParent):
    
    def __init__(self,mptCmdStruct=None):
        
        super(CvtToFuzzyCat, self).__init__(mptCmdStruct)

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
            'Metadata':'Tuple List: Field Name:Any',
            'NewNCVarName':'Field Name',            
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
        
        in_fld_nm = self.ArgByNm('InFieldName')
        in_NC_dim_var = executedObjects[in_fld_nm].ExecRslt()
        in_arr = executedObjects[in_fld_nm].ExecRslt().data.data
        
        # Result starts with a copy of the first input field, then add the rest
        # and divide by the number of them
           
        rslt_arr = np.ma.empty(
            in_arr.shape,
            dtype=float
            )

        rslt_arr.data[:] = float(self.ArgByNm('DefaultFuzzyValue'))

        fuzzyVals = self.ValFromArgByNm('FuzzyValues')
        rawVals = self.ValFromArgByNm('RawValues')
        
        for rawVal,fuzzyVal in zip(rawVals,fuzzyVals):
            np.place(rslt_arr.data,in_arr.data == rawVal,fuzzyVal)

        rslt_arr.mask = cp.deepcopy(in_arr.mask)

        # bring back values to fuzzy limits
        self._InsureFuzzy(rslt_arr)

        self.execRslt = mpncv.NCDimensionedVar(
            dims = cp.deepcopy(in_NC_dim_var.dims),
            name = self._CreateNCVarName(),
            data = mpncv.NCVar(
                data = rslt_arr,
                dim_nms = in_NC_dim_var.dims.keys(),
                )
            )
        
        self._SetMetadata()
        self.execRslt.set_metadata_val('NodeCommand',self.mptCmdStruct['rawCmdStr'])
        executedObjects[self.RsltNm()] = self
        
    # def Exec(self,executedObjects):
    
# class CvtToFuzzyCat(_FuzzyParent):

class CvtToFuzzy(_CvtToFuzzyCurve):

    def __init__(
        self,
        mpt_cmd_struct=None,
        ):

        super(CvtToFuzzy,self).__init__(mpt_cmd_struct)

    # def __init__(...)

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
            'Metadata':'Tuple List: Field Name:Any',
            'NewNCVarName':'Field Name',            
            'Direction':'One of| LowToHigh HighToLow',
            'PrecursorFieldNames':['Field Name','Field Name List'],
            }
        
    # _SetFxnDesc(self):
    
    def Exec(self,executedObjects):
        # executedObjects is a dictionary of executed MPilot function objects

        in_fld_nm = self.ArgByNm('InFieldName')
        in_NC_dim_var = executedObjects[in_fld_nm].ExecRslt()
        in_arr = executedObjects[in_fld_nm].ExecRslt().data.data

        direction = self.ArgByNm('Direction')
        if direction is None:
            direction = 'LowToHigh'  # default direction
        
        if 'FalseThreshold' in self.Args():
            false_thresh = self.ValFromArgByNm('FalseThreshold')
        elif direction in ['HighToLow']:
            false_thresh = in_arr.max()
        else:
            false_thresh =  in_arr.min()
            
        if 'TrueThreshold' in self.Args():
            true_thresh = self.ValFromArgByNm('TrueThreshold')
        elif direction in ['HighToLow']:
            true_thresh =  in_arr.min()
        else:
            true_thresh =  in_arr.max()

        if true_thresh == false_thresh:
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

        self.raw_vals = [false_thresh,true_thresh]
        self.fuzzy_vals = [self.fuzzyMin,self.fuzzyMax]
            
        self._Exec(executedObjects)
    # def Exec(self,executedObjects):
    
# class CvtToFuzzy(_CvtToFuzzyCurve):

class CvtToFuzzyCurve(_CvtToFuzzyCurve):
   
    def __init__(self,mptCmdStruct=None):
        
        super(CvtToFuzzyCurve, self).__init__(mptCmdStruct)

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
            'Metadata':'Tuple List: Field Name:Any',
            'NewNCVarName':'Field Name',            
            'PrecursorFieldNames':['Field Name','Field Name List'],
            }
        
    # _SetFxnDesc(self):

    def DependencyNms(self):
        
        rtrn = self._ArgToList('InFieldName')
        rtrn += self._ArgToList('PrecursorFieldNames')
        return rtrn
    
    # def DependencyNms(self):

    def Exec(self,executedObjects):
        
        self.fuzzy_vals = self.ValFromArgByNm('FuzzyValues')
        self.raw_vals = self.ValFromArgByNm('RawValues')

        self._Exec(executedObjects)
        
    # def Exec(self,executedObjects):

# class CvtToFuzzyCurve(_FuzzyParent):

class CvtToFuzzyPercentiles(_CvtToFuzzyCurve):
    
    def __init__(self,mptCmdStruct=None):
        
        super(CvtToFuzzyPercentiles, self).__init__(mptCmdStruct)

    # def __init__(self,mptCmdStruct=None):

    def _SetFxnDesc(self):
        # description of command for validation and info display

        self.fxnDesc['DisplayName'] = 'Convert to Fuzzy using Percentile'
        self.fxnDesc['ShortDesc'] = 'Converts input values into fuzzy values using linear interpolation based on percentiles.'
        self.fxnDesc['ReturnType'] = 'Fuzzy'

        self.fxnDesc['ReqArgs'] = {
            'InFieldName':'Field Name',
            'Percentiles':['Float List'],
            'FuzzyValues':['Fuzzy Value List'],
            }
        self.fxnDesc['OptArgs'] = {
            'Metadata':'Tuple List: Field Name:Any',
            'NewNCVarName':'Field Name',            
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
        
        in_fld_nm = self.ArgByNm('InFieldName')
        in_NC_dim_var = executedObjects[in_fld_nm].ExecRslt()
        in_arr = executedObjects[in_fld_nm].ExecRslt().data.data

        self.fuzzy_vals = self.ValFromArgByNm('FuzzyValues')
        percentiles = self.ValFromArgByNm('Percentiles')

        if len(self.fuzzy_vals) != len(percentiles):
            raise Exception(
                (3*'{}\n').format(
                    'Error: Number of percentiles and raw values must match.',
                    '  Result name: {}'.format(self.RsltNm),
                    '  Number of percentiles, fuzzy values: {}'.format(len(percentiles),len(self.fuzzy_vals))
                    )
                ) # raise Exception(...) 

        compressed_arr = in_arr.compressed()
        self.raw_vals = []
        for percentile in percentiles:
            self.raw_vals.append(np.percentile(compressed_arr,percentile))

        self._Exec(executedObjects)

    # def Exec(self,executedObjects):
    
# class CvtToFuzzyPercentile(_CvtToFuzzyCurve):

class CvtToBinary(_NetCDFUtilParent):
    
    def __init__(self,mptCmdStruct=None):
        
        super(type(self), self).__init__(mptCmdStruct)

    # def __init__(self,mptCmdStruct=None):

    def _SetFxnDesc(self):
        # description of command for validation and info display

        self.fxnDesc['DisplayName'] = 'Convert to Binary'
        self.fxnDesc['ShortDesc'] = '''Converts input values into binary 0 or 1 based on threshold.
Direction = LowToLow for values below threshold to be 0 and above to be 1.
Direction = HighToLow for values below threshold to be 1 and above to be 0.'''
        self.fxnDesc['ReturnType'] = 'Float'

        self.fxnDesc['ReqArgs'] = {
            'InFieldName':'Field Name',
            'Threshold':'Float',
            'Direction':'One of| LowToLow HighToLow',
            }
        self.fxnDesc['OptArgs'] = {
            'NewNCVarName':'Field Name',
            'Metadata':'Tuple List: Field Name:Any',
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
        
        in_fld_nm = self.ArgByNm('InFieldName')
        in_NC_dim_var = executedObjects[in_fld_nm].ExecRslt()
        in_arr = executedObjects[in_fld_nm].ExecRslt().data.data

        if self.ArgByNm('Direction') not in ['LowToLow','HighToLow']:
            raise Exception(
                '{}{}{}{}{}{}{}'.format(
                    '\n********************ERROR********************\n',
                    'Invalid Direction specified in command:\n',
                    '  Direction must be one of:\n',
                    '    LowToLow\n',
                    '    HighToLow\n',
                    'File: {}  Line number: {}\n'.format(
                        self.mptCmdStruct['cmdFileNm'],
                        self.mptCmdStruct['lineNo']
                        ),
                    'Full command:\n{}\n'.format(self.mptCmdStruct['rawCmdStr'])
                ),
            )

        # Convert to binary
        if self.ArgByNm('Direction') == 'LowToLow':
            lowVal = 0.
            hiVal = 1.
        else:
            lowVal = 1.
            hiVal = 0.
            
        rslt_arr = np.ma.where(
            in_arr < self.ValFromArgByNm('Threshold'),
            lowVal,
            hiVal
            )
            
        self.execRslt = mpncv.NCDimensionedVar(
            dims = cp.deepcopy(in_NC_dim_var.dims),
            name = self._CreateNCVarName(),
            data = mpncv.NCVar(
                data = rslt_arr,
                dim_nms = in_NC_dim_var.dims.keys(),
                )
            )

        self._SetMetadata()
        self.execRslt.set_metadata_val('NodeCommand',self.mptCmdStruct['rawCmdStr'])
        executedObjects[self.RsltNm()] = self
     
    # def Exec(self,executedObjects):
    
# class CvtToBinary(_FuzzyParent):

# Fuzzy logic operations

class FuzzyUnion(_FuzzyParent):

    def __init__(self,mptCmdStruct=None):

        super(FuzzyUnion, self).__init__(mptCmdStruct)

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
            'NewNCVarName':'Field Name',
            'Metadata':'Tuple List: Field Name:Any',
            'PrecursorFieldNames':['Field Name','Field Name List']
            }
        
    # _SetFxnDesc(self):

    def DependencyNms(self):
        
        rtrn = self._ArgToList('InFieldNames')
        rtrn += self._ArgToList('PrecursorFieldNames')
        return rtrn
    
    # def DependencyNms(self):

    def Exec(self,executedObjects):
        # executedObjects is a dictionary of executed MPilot function objects
        
        self._ValidateListLen('InFldNms',1)
        fldnm_list = self._ArgToList('InFieldNames')
        
        # Result starts with a copy of the first input field, then add the rest
        # and divide by the number of them
        
        rslt_arr = cp.deepcopy(executedObjects[fldnm_list[0]].ExecRslt().data.data)
        
        for fld_nm in fldnm_list[1:]:
            rslt_arr += executedObjects[fld_nm].ExecRslt().data.data

        rslt_arr = rslt_arr / float(len(fldnm_list))

        # for dimensioning result variable
        in_NC_dim_var = executedObjects[fldnm_list[0]].ExecRslt()

        self.execRslt = mpncv.NCDimensionedVar(
            dims = cp.deepcopy(in_NC_dim_var.dims),
            name = self._CreateNCVarName(),
            data = mpncv.NCVar(
                data = rslt_arr,
                dim_nms = in_NC_dim_var.dims.keys(),
                )
            )

        self._SetMetadata()
        self.execRslt.set_metadata_val('NodeCommand',self.mptCmdStruct['rawCmdStr'])
        executedObjects[self.RsltNm()] = self
        
    # def Exec(self,executedObjects):
    
# class FuzzyUnion(_FuzzyParent):    

class FuzzyWeightedUnion(_FuzzyParent):

    def __init__(self,mptCmdStruct=None):

        super(FuzzyWeightedUnion, self).__init__(mptCmdStruct)

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
            'NewNCVarName':'Field Name',
            'Metadata':'Tuple List: Field Name:Any',
            'PrecursorFieldNames':['Field Name','Field Name List']
            }
        
    # _SetFxnDesc(self):
    
    def DependencyNms(self):
        
        rtrn = self._ArgToList('InFieldNames')
        rtrn += self._ArgToList('PrecursorFieldNames')
        return rtrn
    
    # def DependencyNms(self):

    def Exec(self,executedObjects):
        # executedObjects is a dictionary of executed MPilot function objects
        
        self._ValidateListLen('InFldNms',1)
        self._ValidateEqualListLens(['InFieldNames','Weights'])
        
        fldnm_list = self._ArgToList('InFieldNames')
        wts = self.ValFromArgByNm('Weights')
        
        # Result starts with a copy of the first input field, then add the rest
        # and divide by the number of them
        
        rslt_arr = cp.deepcopy(executedObjects[fldnm_list[0]].ExecRslt().data.data) * wts[0]
        
        for wt,fldnm in zip(wts,fldnm_list)[1:]:
            rslt_arr += executedObjects[fldnm].ExecRslt().data.data * wt

        rslt_arr = rslt_arr / sum(wts)
        self._InsureFuzzy(rslt_arr)

        # for dimensioning result variable
        in_NC_dim_var = executedObjects[fldnm_list[0]].ExecRslt()

        self.execRslt = mpncv.NCDimensionedVar(
            dims = cp.deepcopy(in_NC_dim_var.dims),
            name = self._CreateNCVarName(),
            data = mpncv.NCVar(
                data = rslt_arr,
                dim_nms = in_NC_dim_var.dims.keys(),
                )
            )

        self._SetMetadata()
        self.execRslt.set_metadata_val('NodeCommand',self.mptCmdStruct['rawCmdStr'])
        executedObjects[self.RsltNm()] = self
        
    # def Exec(self,executedObjects):
    
# class FuzzyWeightedUnion(_FuzzyParent):    

class FuzzySelectedUnion(_FuzzyParent):

    def __init__(self,mptCmdStruct=None):

        super(FuzzySelectedUnion, self).__init__(mptCmdStruct)

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
            'TruestOrFalsest':'One of| Truest Falsest',
            'NumberToConsider':'Positive Integer'
            }
        self.fxnDesc['OptArgs'] = {
            'NewNCVarName':'Field Name',
            'Metadata':'Tuple List: Field Name:Any',
            'PrecursorFieldNames':['Field Name','Field Name List']
            }
        
    # _SetFxnDesc(self):

    def DependencyNms(self):
        
        rtrn = self._ArgToList('InFieldNames')
        rtrn += self._ArgToList('PrecursorFieldNames')
        return rtrn
    
    # def DependencyNms(self):

    def Exec(self,executedObjects):
        # executedObjects is a dictionary of executed MPilot function objects
        fldnm_list = self._ArgToList('InFieldNames')
        numToCnsdr = self.ValFromArgByNm('NumberToConsider')
        
        if len(fldnm_list) < numToCnsdr:
            raise Exception(
                '{}{}{}{}{}'.format(
                    '\n********************ERROR********************\n',
                    'Number of InFieldNames must be greater than or equal to NumberToConsider:\n',
                    'Number of InFieldNames: {}  NumberToConsider: {}\n'.format(
                        len(fldnm_list),
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

        tmpMask = executedObjects[fldnm_list[0]].ExecRslt().data.data.mask
        for fldnm in fldnm_list[1:]:
            tmpMask = np.logical_or(tmpMask,executedObjects[fldnm].ExecRslt().data.data.mask)

        stackedArr = np.ma.array(
            np.stack([executedObjects[fldnm].ExecRslt().data.data.data for fldnm in fldnm_list]),
            mask = cp.deepcopy(
                np.broadcast_to(
                    tmpMask,
                    [len(fldnm_list)]+list(executedObjects[fldnm].ExecRslt().data.data.shape)
                    )
                )
            )
            
        stackedArr.sort(axis=0,kind='heapsort')

        if self.ValFromArgByNm('TruestOrFalsest') == 'Truest':
            rslt_arr = np.ma.mean(stackedArr[-numToCnsdr:],axis=0)
        else:
            rslt_arr = np.ma.mean(stackedArr[0:numToCnsdr],axis=0)
            
        # for dimensioning result variable
        in_NC_dim_var = executedObjects[fldnm_list[0]].ExecRslt()

        self.execRslt = mpncv.NCDimensionedVar(
            dims = cp.deepcopy(in_NC_dim_var.dims),
            name = self._CreateNCVarName(),
            data = mpncv.NCVar(
                data = rslt_arr,
                dim_nms = in_NC_dim_var.dims.keys(),
                )
            )

        self._SetMetadata()
        self.execRslt.set_metadata_val('NodeCommand',self.mptCmdStruct['rawCmdStr'])
        executedObjects[self.RsltNm()] = self
        
    # def Exec(self,executedObjects):
    
# class FuzzySelectedUnion(_FuzzyParent):

class FuzzyOr(_FuzzyParent):

    def __init__(self,mptCmdStruct=None):

        super(FuzzyOr, self).__init__(mptCmdStruct)

    # def __init__(self,mptCmdStruct=None):
    
    def _SetFxnDesc(self):

        self.fxnDesc['DisplayName'] = 'Fuzzy Or'
        self.fxnDesc['ShortDesc'] = 'Takes the fuzzy Or (maximum) of fuzzy input variables'
        self.fxnDesc['ReturnType'] = 'Fuzzy'
        
        self.fxnDesc['ReqArgs'] = {
            'InFieldNames':['Field Name','Field Name List']
            }
        self.fxnDesc['OptArgs'] = {
            'NewNCVarName':'Field Name',
            'Metadata':'Tuple List: Field Name:Any',
            'PrecursorFieldNames':['Field Name','Field Name List']
            }
        
    # _SetFxnDesc(self):

    def DependencyNms(self):
        
        rtrn = self._ArgToList('InFieldNames')
        rtrn += self._ArgToList('PrecursorFieldNames')
        return rtrn
    
    # def DependencyNms(self):

    def Exec(self,executedObjects):

        self._ValidateListLen('InFldNms',1)
        fldnm_list = self._ArgToList('InFieldNames')

        # Traverse inputs and take maximum
        rslt_arr = cp.deepcopy(executedObjects[fldnm_list[0]].ExecRslt().data.data)
        
        for fldnm in fldnm_list[1:]:
            rslt_arr = np.ma.maximum(rslt_arr,executedObjects[fldnm].ExecRslt().data.data)

        self._InsureFuzzy(rslt_arr)

        # for dimensioning result variable
        in_NC_dim_var = executedObjects[fldnm_list[0]].ExecRslt()

        self.execRslt = mpncv.NCDimensionedVar(
            dims = cp.deepcopy(in_NC_dim_var.dims),
            name = self._CreateNCVarName(),
            data = mpncv.NCVar(
                data = rslt_arr,
                dim_nms = in_NC_dim_var.dims.keys(),
                )
            )

        self._SetMetadata()
        self.execRslt.set_metadata_val('NodeCommand',self.mptCmdStruct['rawCmdStr'])
        executedObjects[self.RsltNm()] = self
        
    # def Exec(self,executedObjects):

# class FuzzyOr(_FuzzyParent):

class FuzzyAnd(_FuzzyParent):

    def __init__(self,mptCmdStruct=None):

        super(FuzzyAnd, self).__init__(mptCmdStruct)

    # def __init__(self,mptCmdStruct=None):

    def _SetFxnDesc(self):
        
        self.fxnDesc['DisplayName'] = 'Fuzzy And'
        self.fxnDesc['ShortDesc'] = 'Takes the fuzzy And (minimum) of fuzzy input variables'
        self.fxnDesc['ReturnType'] = 'Fuzzy'
        
        self.fxnDesc['ReqArgs'] = {
            'InFieldNames':['Field Name','Field Name List']
            }
        self.fxnDesc['OptArgs'] = {
            'NewNCVarName':'Field Name',
            'Metadata':'Tuple List: Field Name:Any',
            'PrecursorFieldNames':['Field Name','Field Name List']
            }
        
    # _SetFxnDesc(self):

    def DependencyNms(self):
        
        rtrn = self._ArgToList('InFieldNames')
        rtrn += self._ArgToList('PrecursorFieldNames')
        return rtrn
    
    # def DependencyNms(self):

    def Exec(self,executedObjects):

        self._ValidateListLen('InFldNms',1)
        fldnm_list = self._ArgToList('InFieldNames')
        
        # Traverse inputs and take maximum
        rslt_arr = cp.deepcopy(executedObjects[fldnm_list[0]].ExecRslt().data.data)

        for fldnm in fldnm_list[1:]:
            rslt_arr = np.ma.minimum(rslt_arr,executedObjects[fldnm].ExecRslt().data.data)

        self._InsureFuzzy(rslt_arr)

        # for dimensioning result variable
        in_NC_dim_var = executedObjects[fldnm_list[0]].ExecRslt()

        self.execRslt = mpncv.NCDimensionedVar(
            dims = cp.deepcopy(in_NC_dim_var.dims),
            name = self._CreateNCVarName(),
            data = mpncv.NCVar(
                data = rslt_arr,
                dim_nms = in_NC_dim_var.dims.keys(),
                )
            )

        self._SetMetadata()
        self.execRslt.set_metadata_val('NodeCommand',self.mptCmdStruct['rawCmdStr'])
        executedObjects[self.RsltNm()] = self
                
    # def Exec(self,executedObjects):
    
# class FuzzyAnd(_FuzzyParent):

class FuzzyXOr(_FuzzyParent):

    def __init__(self,mptCmdStruct=None):

        super(FuzzyXOr, self).__init__(mptCmdStruct)

    # def __init__(self,mptCmdStruct=None):

    def _SetFxnDesc(self):
        
        self.fxnDesc['DisplayName'] = 'Fuzzy XOr'
        self.fxnDesc['ShortDesc'] = 'Computes Fuzzy XOr: Truest - (Truest - 2nd Truest) * (2nd Truest - full False)/(Truest - full False)'
        self.fxnDesc['ReturnType'] = 'Fuzzy'
        
        self.fxnDesc['ReqArgs'] = {
            'InFieldNames':['Field Name List']
            }
        self.fxnDesc['OptArgs'] = {
            'NewNCVarName':'Field Name',
            'Metadata':'Tuple List: Field Name:Any',
            'PrecursorFieldNames':['Field Name','Field Name List']
            }
        
    # _SetFxnDesc(self):

    def DependencyNms(self):
        
        rtrn = self._ArgToList('InFieldNames')
        rtrn += self._ArgToList('PrecursorFieldNames')
        return rtrn
    
    # def DependencyNms(self):

    def Exec(self,executedObjects):

        self._ValidateListLen('InFieldNames',2)
        fldnm_list = self._ArgToList('InFieldNames')
        
        # Create a stacked array with layers from input arrays, sort it, and use
        # that to calculate fuzzy xor. There is no np.ma.stacked, so the masks have
        # to be handled separately from the data. Note we are building the maximal
        # mask from all the inputs before broadcasting it to the size of the stacked
        # array. There are some issues with getting stackedArr to be writable. The
        # below code works.

        tmpMask = executedObjects[fldnm_list[0]].ExecRslt().data.data.mask
        for fldnm in fldnm_list[1:]:
            tmpMask = np.logical_or(tmpMask,executedObjects[fldnm].ExecRslt().data.data.mask)
            
        stackedArr = np.ma.array(
            np.stack([executedObjects[fldnm].ExecRslt().data.data.data for fldnm in fldnm_list]),
            mask = cp.deepcopy(
                np.broadcast_to(
                    tmpMask,
                    [len(fldnm_list)]+list(executedObjects[fldnm].ExecRslt().data.data.shape)
                    )
                )
            )
            
        stackedArr.sort(axis=0,kind='heapsort')

        # the array multiply here triggers a multiplication overflow warning.
        # does not seem to affect results
        rslt_arr = np.ma.where(
            stackedArr[-1] <= self.fuzzyMin,
            self.fuzzyMin,
            stackedArr[-1] - \
                (stackedArr[-1] - stackedArr[-2]) * \
                (stackedArr[-2] - self.fuzzyMin) / \
                (stackedArr[-1] - self.fuzzyMin)
            )

        self._InsureFuzzy(rslt_arr)

        # for dimensioning result variable
        in_NC_dim_var = executedObjects[fldnm_list[0]].ExecRslt()

        self.execRslt = mpncv.NCDimensionedVar(
            dims = cp.deepcopy(in_NC_dim_var.dims),
            name = self._CreateNCVarName(),
            data = mpncv.NCVar(
                data = rslt_arr,
                dim_nms = in_NC_dim_var.dims.keys(),
                )
            )

        self._SetMetadata()
        self.execRslt.set_metadata_val('NodeCommand',self.mptCmdStruct['rawCmdStr'])
        executedObjects[self.RsltNm()] = self
                
    # def Exec(self,executedObjects):
    
# class FuzzyXOr(_FuzzyParent):

class FuzzyNot(_FuzzyParent):

    def __init__(self,mptCmdStruct=None):

        super(FuzzyNot, self).__init__(mptCmdStruct)

    # def __init__(self,mptCmdStruct=None):

    def _SetFxnDesc(self):
        
        self.fxnDesc['DisplayName'] = 'Fuzzy Not'
        self.fxnDesc['ShortDesc'] = 'Takes the fuzzy And (minimum) of fuzzy input variables'
        self.fxnDesc['ReturnType'] = 'Fuzzy'
        
        self.fxnDesc['ReqArgs'] = {
            'InFieldName':'Field Name'
            }
        self.fxnDesc['OptArgs'] = {
            'NewNCVarName':'Field Name',
            'Metadata':'Tuple List: Field Name:Any',
            'PrecursorFieldNames':['Field Name','Field Name List']
            }
        
    # _SetFxnDesc(self):

    def DependencyNms(self):
        
        rtrn = self._ArgToList('InFieldName')
        rtrn += self._ArgToList('PrecursorFieldNames')
        return rtrn
    
    # def DependencyNms(self):

    def Exec(self,executedObjects):

        in_fld_nm = self.ArgByNm('InFieldName')
        rslt_arr = -executedObjects[in_fld_nm].ExecRslt().data.data

        self._InsureFuzzy(rslt_arr)

        # for dimensioning result variable
        in_NC_dim_var = executedObjects[in_fld_nm].ExecRslt()

        self.execRslt = mpncv.NCDimensionedVar(
            dims = cp.deepcopy(in_NC_dim_var.dims),
            name = self._CreateNCVarName(),
            data = mpncv.NCVar(
                data = rslt_arr,
                dim_nms = in_NC_dim_var.dims.keys(),
                )
            )

        self._SetMetadata()
        self.execRslt.set_metadata_val('NodeCommand',self.mptCmdStruct['rawCmdStr'])
        executedObjects[self.RsltNm()] = self
        
    # def Exec(self,executedObjects):

# class FuzzyNot(_FuzzyParent):    
    
class FuzzyMultiply(_FuzzyParent):

    def __init__(self,mptCmdStruct=None):

        super(FuzzyMultiply, self).__init__(mptCmdStruct)

    # def __init__(self,mptCmdStruct=None):

    def _SetFxnDesc(self):
        
        self.fxnDesc['DisplayName'] = 'Fuzzy Product'
        self.fxnDesc['ShortDesc'] = 'Takes the procduct of fuzzy input variables. Suitable for risk (likelihood * value).'
        self.fxnDesc['ReturnType'] = 'Fuzzy'
        
        self.fxnDesc['ReqArgs'] = {
            'InFieldNames':['Field Name','Field Name List']
            }
        self.fxnDesc['OptArgs'] = {
            'NewNCVarName':'Field Name',
            'Metadata':'Tuple List: Field Name:Any',
            'PrecursorFieldNames':['Field Name','Field Name List']
            }
        
    # _SetFxnDesc(self):

    def DependencyNms(self):
        
        rtrn = self._ArgToList('InFieldNames')
        rtrn += self._ArgToList('PrecursorFieldNames')
        return rtrn
    
    # def DependencyNms(self):

    def Exec(self,executedObjects):

        self._ValidateListLen('InFldNms',1)
        fldnm_list = self._ArgToList('InFieldNames')

        if self.fuzzyMin != 0. or self.fuzzyMax != 1. :
            raise Exception(
                '{}{}{}'.format(
                    '\n********************ERROR********************\n',
                    'Fuzzy Product only suitable when Fuzzy values range from 0 to 1.\n',
                    'Fuzzy range is {} to {}'.format(self.fuzzyMin,self.fuzzyMax)
                    )
                )
            
        fldnm_list = self._ArgToList('InFieldNames')

        if len(fldnm_list) < 2:
            raise Exception(
                (2*'{}\n').format(
                    'Error. Must have 2 or more Field Names',
                    '  Result Name, Field Name List: {}  [{}]'.format(
                        self.RsltNm(),
                        ','.join(fldnm_list)
                        ),
                    )
                ) # raise Exception(...)

        # if len(fldnm_list) < 2:

        new_ncdimvar = cp.deepcopy(executedObjects[fldnm_list[0]].ExecRslt())
        new_ncdimvar.name = self._CreateNCVarName()
        
        for fldnm in fldnm_list[1:]:
            
            ncdimvar_for_op = executedObjects[fldnm].ExecRslt()

            # If necessary, expand the array being applied to new_ncdimvar

            if self._NCDimArrsAreCongruent(ncdimvar_for_op,new_ncdimvar):

                delete_ncdimvar_for_op = False
                
            else:

                ncdimvar_for_op = self._MakeCongruentNCDimArr(ncdimvar_for_op,new_ncdimvar)

                if ncdimvar_for_op is not None:
                    delete_ncdimvar_for_op = True
                else:
                    raise Exception(
                        (4*'{}\n').format(
                            'Error. Cannot make array congruent to first array in arithmetic array operation.',
                            '  Result name: {}'.format(self.RsltNm()),
                            '  First array: Field Name, NCVar name, dim names, shape: {}  {}  {}  {}'.format(
                                fldnm_list[0],
                                executedObjects[fldnm_list[0]].ExecRslt().name,
                                executedObjects[fldnm_list[0]].ExecRslt().dims.keys(),
                                executedObjects[fldnm_list[0]].ExecRslt().data.data.shape,
                                ),
                            '  Incongruent array: Field Name, NCVar name, dim names, shape: {}  {}  {}  {}'.format(
                                fldnm_list[0],
                                executedObjects[fldnm].ExecRslt().name,
                                executedObjects[fldnm].ExecRslt().dims.keys(),
                                executedObjects[fldnm].ExecRslt().data.data.shape,
                                ),
                            )
                        ) # raise Exception(...)
                        
            # if self._NCDimArrsAreCongruent(ncdimvar_for_op,new_ncdimvar):...else...

            # do the operation
            new_ncdimvar.data.data = new_ncdimvar.data.data * ncdimvar_for_op.data.data

            if delete_ncdimvar_for_op:
                del ncdimvar_for_op
                    
        # for fldnm in fldnm_list[1:]:

        self.execRslt = new_ncdimvar
        self._SetMetadata()
        self.execRslt.set_metadata_val('NodeCommand',self.mptCmdStruct['rawCmdStr'])
        executedObjects[self.RsltNm()] = self
        
    # def Exec(self,executedObjects):

#class FuzzyMultiply(_DoArrayArithmetic):

class CvtFromFuzzy(_FuzzyParent):
    
    def __init__(self,mptCmdStruct=None):
        
        super(type(self), self).__init__(mptCmdStruct)

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
            'NewNCVarName':'Field Name',
            'Metadata':'Tuple List: Field Name:Any',
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
        
        in_fld_nm = self.ArgByNm('InFieldName')
        in_arr = executedObjects[in_fld_nm].ExecRslt().data.data
        
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
        rslt_arr = (in_arr - x1) * (y2-y1)/(x2-x1) + y1

        # for dimensioning result variable
        in_NC_dim_var = executedObjects[self.ValFromArgByNm('InFieldName')].ExecRslt()

        self.execRslt = mpncv.NCDimensionedVar(
            dims = cp.deepcopy(in_NC_dim_var.dims),
            name = self._CreateNCVarName(),
            data = mpncv.NCVar(
                data = rslt_arr,
                dim_nms = in_NC_dim_var.dims.keys(),
                )
            )

        self._SetMetadata()
        self.execRslt.set_metadata_val('NodeCommand',self.mptCmdStruct['rawCmdStr'])
        executedObjects[self.RsltNm()] = self
                
    # def Exec(self,executedObjects):
    
# class CvtFromFuzzy(_FuzzyParent):

class FuzzyUnionOverDimensions(_SummaryOverDimensions):

    def __init__(self,mptCmdStruct=None):

        super(FuzzyUnionOverDimensions, self).__init__(mptCmdStruct)

    # def __init__(self,mptCmdStruct=None):
    
    def _SetFxnDesc(self):
        # description of command used for validation and
        # information display Each Pilot fxn command should have its
        # own description

        self.fxnDesc['DisplayName'] = 'Fuzzy Union Over Dimensions'
        self.fxnDesc['ShortDesc'] = 'Takes the fuzzy Union (mean) of a fuzzy input'
        self.fxnDesc['ReturnType'] = 'NCVar'

        self.fxnDesc['ReqArgs'] = {
            'InFieldName':'Field Name',
            'Dimensions':['Field Name','Field Name List']
            }
        self.fxnDesc['OptArgs'] = {
            'NewNCVarName':'Field Name',
            'Metadata':'Tuple List: Field Name:Any',
            'PrecursorFieldNames':['Field Name','Field Name List']
            }
        
    # _SetFxnDesc(self):

    def DependencyNms(self):
        
        rtrn = self._ArgToList('InFieldName')
        rtrn += self._ArgToList('PrecursorFieldNames')
        return rtrn
    
    # def DependencyNms(self):

    def Exec(self,executedObjects):
        
        self._Exec(executedObjects,'mean')
        
    # def Exec(self,executedObjects):
    
# class FuzzyUnionOverDimensions(_SummaryOverDimensions):

class FuzzyAndOverDimensions(_SummaryOverDimensions):

    def __init__(self,mptCmdStruct=None):

        super(FuzzyAndOverDimensions, self).__init__(mptCmdStruct)

    # def __init__(self,mptCmdStruct=None):
    
    def _SetFxnDesc(self):
        # description of command used for validation and
        # information display Each Pilot fxn command should have its
        # own description

        self.fxnDesc['DisplayName'] = 'Fuzzy And Over Dimensions'
        self.fxnDesc['ShortDesc'] = 'Takes the fuzzy And (minimum) of a fuzzy input'
        self.fxnDesc['ReturnType'] = 'NCVar'

        self.fxnDesc['ReqArgs'] = {
            'InFieldName':'Field Name',
            'Dimensions':['Field Name','Field Name List']
            }
        self.fxnDesc['OptArgs'] = {
            'NewNCVarName':'Field Name',
            'Metadata':'Tuple List: Field Name:Any',
            'PrecursorFieldNames':['Field Name','Field Name List']
            }
        
    # _SetFxnDesc(self):

    def DependencyNms(self):
        
        rtrn = self._ArgToList('InFieldName')
        rtrn += self._ArgToList('PrecursorFieldNames')
        return rtrn
    
    # def DependencyNms(self):

    def Exec(self,executedObjects):
        
        self._Exec(executedObjects,'min')
        
    # def Exec(self,executedObjects):
    
# class FuzzyAndOverDimensions(_SummaryOverDimensions):

class FuzzyOrOverDimensions(_SummaryOverDimensions):

    def __init__(self,mptCmdStruct=None):

        super(FuzzyOrOverDimensions, self).__init__(mptCmdStruct)

    # def __init__(self,mptCmdStruct=None):
    
    def _SetFxnDesc(self):
        # description of command used for validation and
        # information display Each Pilot fxn command should have its
        # own description

        self.fxnDesc['DisplayName'] = 'Fuzzy And Over Dimensions'
        self.fxnDesc['ShortDesc'] = 'Takes the fuzzy Or (maximum) of a fuzzy input'
        self.fxnDesc['ReturnType'] = 'NCVar'

        self.fxnDesc['ReqArgs'] = {
            'InFieldName':'Field Name',
            'Dimensions':['Field Name','Field Name List']
            }
        self.fxnDesc['OptArgs'] = {
            'NewNCVarName':'Field Name',
            'Metadata':'Tuple List: Field Name:Any',
            'PrecursorFieldNames':['Field Name','Field Name List']
            }
        
    # _SetFxnDesc(self):

    def DependencyNms(self):
        
        rtrn = self._ArgToList('InFieldName')
        rtrn += self._ArgToList('PrecursorFieldNames')
        return rtrn
    
    # def DependencyNms(self):

    def Exec(self,executedObjects):
        
        self._Exec(executedObjects,'max')
        
    # def Exec(self,executedObjects):
    
# class FuzzyOrOverDimensions(_SummaryOverDimensions):

class FuzzySelectedUnionOverDimensions(_FuzzyParent):

    def __init__(self,mptCmdStruct=None):

        super(FuzzySelectedUnionOverDimensions, self).__init__(mptCmdStruct)

    # def __init__(self,mptCmdStruct=None):
    
    def _SetFxnDesc(self):
        # description of command used for validation and
        # information display Each Pilot fxn command should have its
        # own description

        self.fxnDesc['DisplayName'] = 'Fuzzy Selected Union Over Dimensions'
        self.fxnDesc['ShortDesc'] = 'Takes the fuzzy Union (mean) of N Truest or Falsest fuzzy values over dimension(s)'
        self.fxnDesc['ReturnType'] = 'NCVar'

        self.fxnDesc['ReqArgs'] = {
            'InFieldName':'Field Name',
            'Dimensions':['Field Name','Field Name List'],
            'TruestOrFalsest':'One of| Truest Falsest',
            'NumberToConsider':'Positive Integer'
            }
        self.fxnDesc['OptArgs'] = {
            'NewNCVarName':'Field Name',
            'Metadata':'Tuple List: Field Name:Any',
            'PrecursorFieldNames':['Field Name','Field Name List']
            }
        
    # _SetFxnDesc(self):

    def DependencyNms(self):
        
        rtrn = self._ArgToList('InFieldName')
        rtrn += self._ArgToList('PrecursorFieldNames')
        return rtrn
    
    # def DependencyNms(self):

    def Exec(self,executedObjects):
        # executedObjects is a dictionary of executed MPilot function objects
        
        in_fld_nm = self._ArgToList('InFieldName')
        numToCnsdr = self.ValFromArgByNm('NumberToConsider')
        parent_ncdimvar = executedObjects[self.ValFromArgByNm('InFieldName')].ExecRslt()

        # Insure specified dimensions exist
        # and get dimension indices to operate over
        valid_dim_nms = parent_ncdimvar.dims.keys()
        dim_ndxs = []
        for dim_nm in self._ArgToList('Dimensions'):
            
            if dim_nm in valid_dim_nms:
                dim_ndxs.append(valid_dim_nms.index(dim_nm))
            else:
                raise Exception(
                    (3*'{}\n').format(
                        'Error: Invalid dimension specified.',
                        '  Result name: {}'.format(self.RsltNm()),
                        '  Dimension name: {}'.format(dim_nm)
                        )
                    )

        # for dim_nm in self._ArgToList('Dimensions'):
        
        reshaped_arr = cp.deepcopy(_redimension_arr_for_reduction(parent_ncdimvar.data.data,tuple(dim_ndxs)))
        reshaped_arr.sort(axis=reshaped_arr.ndim-1,kind='heapsort')

        if reshaped_arr.shape[-1] < numToCnsdr:
            raise Exception(
                5*'{}\n'.format(
                    '\n********************ERROR********************\n',
                    'Product of specified dimension sizes must be greater than or equal to NumberToConsider:',
                    'Product of specified dimension sizes: {}  NumberToConsider: {}'.format(
                        len(reshaped_arr.shape[-1]),
                        numToCnsdr,
                        ),
                    'File: {}  Line number: {}'.format(
                        self.mptCmdStruct['cmdFileNm'],
                        self.mptCmdStruct['lineNo']
                        ),
                    'Full command:\n{}'.format(self.mptCmdStruct['rawCmdStr'])
                ),
            )

        dim_slices = []
        for dim_sz in reshaped_arr.shape:
            dim_slices.append(slice(0,dim_sz))

        # Take the mean of the truest or falsest numToCnsdr
        if self.ValFromArgByNm('TruestOrFalsest') == 'Truest':            
            dim_slices[-1] = slice(reshaped_arr.shape[-1]-numToCnsdr,reshaped_arr.shape[-1])
        else:
            dim_slices[-1] = slice(0,numToCnsdr)

        rslt_arr = np.ma.mean(reshaped_arr[dim_slices],axis=reshaped_arr.ndim-1)

        # if self.ValFromArgByNm('TruestOrFalsest') == 'Truest':
        #     rslt_arr = np.ma.mean(reshaped_arr[-numToCnsdr:],axis=reshaped_arr.ndim-1)
        # else:
        #     rslt_arr = np.ma.mean(reshaped_arr[0:numToCnsdr],axis=reshaped_arr.ndim-1)

        # Create dimensions for new variable. Keep correct
        # dimension order
        new_dim_ncvars = OrderedDict()
        for dim_nm in parent_ncdimvar.dims.keys():
            if dim_nm not in self._ArgToList('Dimensions'):
                new_dim_ncvars[dim_nm] = cp.deepcopy(parent_ncdimvar.dims[dim_nm])

        # Create and populate new NCDimensionedVar()

        self.execRslt = mpncv.NCDimensionedVar(
            dims = new_dim_ncvars,
            name = self._CreateNCVarName(),
            data = mpncv.NCVar(
                data = rslt_arr,
                dim_nms = new_dim_ncvars.keys()
                ),
            )

        self._SetMetadata()
        self.execRslt.set_metadata_val('NodeCommand',self.mptCmdStruct['rawCmdStr'])
        executedObjects[self.RsltNm()] = self
    # def Exec(self,executedObjects):
    
# class FuzzySelectedUnion(_FuzzyParent):

################################################################################

class StackNCVars(_NetCDFUtilParent):

    def __init__(
        self,
        mpt_cmd_struct=None,
        ):

        super(StackNCVars,self).__init__(mpt_cmd_struct)

    def _SetFxnDesc(self):
        
        self.fxnDesc['DisplayName'] = 'Stack NCVars'
        self.fxnDesc['ShortDesc'] = 'Takes a set of NCVars and combines them along a new dimension.'
        
        self.fxnDesc['ReturnType'] = 'NCVar'
        
        self.fxnDesc['ReqArgs'] = {
            'InFieldNames':['Field Name List'],
            'NewDimension':'Tuple: Field Name:Float:Float'
            }
        self.fxnDesc['OptArgs'] = {
            'PrecursorFieldNames':['Field Name','Field Name List'],
            'Metadata':'Tuple List: Field Name:Any',
            'NewNCVarName':'Field Name',
            }
        
    # def _SetFxnDesc(self):

    def DependencyNms(self):
        
        rtrn = self._ArgToList('InFieldNames')
        rtrn += self._ArgToList('PrecursorFieldNames')
        return rtrn
    
    # def DependencyNms(self):

    def Exec(self,executedObjects):
        # executedObjects is a dictionary of executed MPilot function objects
        
        fldnm_list = self._ArgToList('InFieldNames')
        parent_ncdimvar = executedObjects[fldnm_list[0]].ExecRslt()
        
        for fldnm in fldnm_list[1:]:
            if executedObjects[fldnm].ExecRslt().data.data.shape != parent_ncdimvar.data.data.shape:
                fldnms_and_shapes = ''.join(
                    ['  {} {}\n'.format(fldnm_tmp,executedObjects[fldnm_tmp].ExecRslt().data.data.shape) for fldnm_tmp in fldnm_list]
                    )
                    
                raise Exception(
                    (5*'{}\n').format(
                        '\n********************ERROR********************',
                        'All fields must have same shape:',
                        'Fields and shapes: {}'.format(fldnms_and_shapes),
                        'File: {}  Line number: {}'.format(
                            self.mptCmdStruct['cmdFileNm'],
                            self.mptCmdStruct['lineNo']
                            ),
                        'Full command:{}'.format(self.mptCmdStruct['rawCmdStr'])
                    ),
                )
        
        # create the stacked array and its mask
        rslt_arr = np.ma.masked_array(
            np.stack([executedObjects[fldnm].ExecRslt().data.data.data for fldnm in fldnm_list],axis=-1),
            mask = np.stack([executedObjects[fldnm].ExecRslt().data.data.mask for fldnm in fldnm_list],axis=-1)
            )

        # create the new dimension's array
        new_dim_nm,new_dim_start,new_dim_step = self.ValFromArgByNm('NewDimension').split(':')
        new_dim_arr = np.ma.masked_array(
            [float(new_dim_start) + x * float(new_dim_step) for x in range(rslt_arr.shape[-1])],
            mask = False
            )

        new_dim_ncvar = mpncv.NCVar(
            data = new_dim_arr,
            name = new_dim_nm,
            dim_nms = new_dim_nm
            )

        # dimensions for the rslt_arr

        rslt_arr_dims = OrderedDict()
        for dim_nm in parent_ncdimvar.dims.keys():
            rslt_arr_dims[dim_nm] = cp.deepcopy(parent_ncdimvar.dims[dim_nm])
        rslt_arr_dims[new_dim_nm] = cp.deepcopy(new_dim_ncvar)

        # Create and populate new NCDimensionedVar()

        self.execRslt = mpncv.NCDimensionedVar(
            dims = rslt_arr_dims,
            name = self._CreateNCVarName(),
            data = mpncv.NCVar(
                data = rslt_arr,
                dim_nms = rslt_arr_dims.keys()
                ),
            )

        self._SetMetadata()
        self.execRslt.set_metadata_val('NodeCommand',self.mptCmdStruct['rawCmdStr'])
        executedObjects[self.RsltNm()] = self
        
    # def Exec(self,executedObjects):
    
# class StackNCVars(_NetCDFUtilParent)

class MaskRandomSample(_NetCDFUtilParent):

    # Take a random sample from a netcdf variable. If only one sample is taken,
    # it returns a netcdf variable with the same dimensions as the input variable.
    # If more than one sample is needed, then it returns a netcdf array with an
    # addional dimension (last dimension) which corresponds to each of the samples.
    # This can be used to get mutually exclusive samples, as you would do in order
    # to have separate training and evaluation data sets, or to get non mutually
    # exclusive samples, as you would do to create a lot of training data sets from
    # a small parent data set.

    def __init__(
        self,
        mpt_cmd_struct=None,
        ):

        super(MaskRandomSample,self).__init__(mpt_cmd_struct)
        
    # def __init__(...)

    def _SetFxnDesc(self):
        
        self.fxnDesc['DisplayName'] = 'Mask Random Sample'
        self.fxnDesc['ShortDesc'] = 'Creates an array masked for random sample from parent. For more than one sample creates an array with additional dimension named SampleNumber.'
        
        self.fxnDesc['ReturnType'] = 'NCVar'
        
        self.fxnDesc['ReqArgs'] = {
            'InFieldName':'Field Name',
            'SampleSize':'Integer',
            }
        self.fxnDesc['OptArgs'] = {
            'PrecursorFieldNames':['Field Name','Field Name List'],
            'Metadata':'Tuple List: Field Name:Any',
            'NewNCVarName':'Field Name',
            'MutuallyExclusiveSamples':'Boolean',
            'NumberSamples':'Integer',
            'MaskInOrOut':'One of| In Out',
            }
        
    # def _SetFxnDesc(self):

    def DependencyNms(self):
        
        rtrn = self._ArgToList('InFieldName')
        rtrn += self._ArgToList('PrecursorFieldNames')
        return rtrn
    
    # def DependencyNms(self):

    def Exec(self,executedObjects):
    
        src_dimarr = executedObjects[self.ValFromArgByNm('InFieldName')].ExecRslt()
        src_marr = src_dimarr.data.data

        sample_sz = self.ValFromArgByNm('SampleSize')
        
        num_samples = self.ValFromArgByNm('NumberSamples')
        if num_samples is None:
            num_samples = 1
            
        mutex_samples = self.ValFromArgByNm('MutuallyExclusiveSamples')
        if mutex_samples is None:
            mutex_samples = True

        mask_in = self.ValFromArgByNm('MaskInOrOut')
        if mask_in is None or mask_in == 'In':
            mask_in = True
        else:
            mask_in = False

        # Check to make sure there are enough data points to create samples
        data_points_req = sample_sz
        if mutex_samples:
            data_points_req = sample_sz * num_samples

        if src_marr.count() < data_points_req:
            raise Exception(
                (4*'{}\n').format(
                    '\n********************ERROR********************\n',
                    'You are requesting sample(s) with more data points than available.',
                    '  Data points available: {}  Data points requested: {}'.format(
                        src_marr.count(),data_points_req
                        ),
                    'Full command:\n{}'.format(self.mptCmdStruct['rawCmdStr'])
                    ),
                )
            
        # Create the array for the output mask and populate it

        if mask_in:
            mask_start_val = True
        else:
            mask_start_val = False
            
        if num_samples == 1:

            out_mask = np.full(src_marr.shape,mask_start_val)

            # insure original mask is in place
            if not mask_in:
                out_mask = np.ma.mask_or(out_mask,src_marr.mask)
                
            mask_element_positions = np.random.choice(
                src_marr.count(),
                size = sample_sz,
                replace = False
                )

            # (~marr.mask).nonzero() are the ndxs for the masked in elements
            # of the starting array.
            
            selected_ndxs = tuple(np.take((~src_marr.mask).nonzero(), mask_element_positions, axis=1))
            
            out_mask[selected_ndxs] = not mask_start_val

            # create the resulting array

            new_marr = np.ma.masked_array(
                cp.deepcopy(src_marr.data),
                mask=out_mask
                )

            # Make the dimensioned array object
            self.execRslt = mpncv.NCDimensionedVar(
                dims = cp.deepcopy(src_dimarr.dims),
                name = self._CreateNCVarName(),
                data = mpncv.NCVar(
                    data = new_marr,
                    dim_nms = src_dimarr.dims.keys()
                    )
                )
            
        else:
                
            out_mask = np.full(list(src_marr.shape)+[num_samples],mask_start_val)
            
            # if masking out variables, need to carry along the original mask
            if not mask_in:
                for ndx in range(num_samples):
                    out_mask[...,ndx] = np.ma.mask_or(out_mask[...,ndx],src_marr.mask)
                
            out_data = np.array(
                np.transpose(
                    np.broadcast_to(
                        src_marr.data,
                        [num_samples] + list(src_marr.data.shape)
                        ),
                    range(1,len(src_marr.shape)+1) +[0] # transpose masking
                    )
                )

            if mutex_samples:
                
                # this guarantees samples are independent between layers
                mask_element_positions = np.random.choice(
                    src_marr.count(),
                    size = sample_sz * num_samples,
                    replace = False
                    )

            else:

                # creating samples for layers independently
                mask_element_positions = []
                for sample_ndx in range(num_samples):
                    mask_element_positions += np.random.choice(
                    src_marr.count(),
                    size = sample_sz,
                    replace = True
                    ).tolist()

                mask_element_positions = np.array(mask_element_positions)
                
            # if mutex_samples:
                
            selected_ndxs = tuple(np.take(
                (~src_marr.mask).nonzero(),
                mask_element_positions,
                axis=1
                ))

            # selected_ndxs have the dimensions for src_marr, so we need to add
            # the dimensions for the sample layer dimension
            layer_ndx_lst = []
            for layer_ndx in range(num_samples):
                layer_ndx_lst += [layer_ndx] * sample_sz
                    
            selected_ndxs = tuple(list(selected_ndxs) + [np.array(layer_ndx_lst)])

            # apply to mask
            out_mask[selected_ndxs] = not mask_start_val
            
            new_marr = np.ma.masked_array(
                out_data,
                mask=out_mask
                )

            dims = cp.deepcopy(src_dimarr.dims)

            dims['SampleNumber'] = mpncv.NCVar(
                data = np.arange(num_samples),
                dim_nms = ('SampleNumber'),
                name = 'SampleNumber'
                )

            # Make the dimensioned array object
            self.execRslt = mpncv.NCDimensionedVar(
                dims = dims,
                name = self._CreateNCVarName(),
                data = mpncv.NCVar(
                    data = new_marr,
                    dim_nms = dims.keys()
                    )
                )
            
        # if num_samples == 1:...else...

        self._SetMetadata()
        self.execRslt.set_metadata_val('NodeCommand',self.mptCmdStruct['rawCmdStr'])
        executedObjects[self.RsltNm()] = self

    # def Exec(self,executedObjects):
        
# class MaskRandomSample(_NetCDFUtilParent):

class InvertMask(_NetCDFUtilParent):
    
    def __init__(
        self,
        mpt_cmd_struct=None,
        ):

        super(InvertMask,self).__init__(mpt_cmd_struct)
        
    # def __init__(...)

    def _SetFxnDesc(self):
        
        self.fxnDesc['DisplayName'] = 'Invert Mask'
        self.fxnDesc['ShortDesc'] = 'Creates an array whose mask is the inverse of the input array.'
        
        self.fxnDesc['ReturnType'] = 'NCVar'
        
        self.fxnDesc['ReqArgs'] = {
            'InFieldName':'Field Name',
            }
        self.fxnDesc['OptArgs'] = {
            'PrecursorFieldNames':['Field Name','Field Name List'],
            'Metadata':'Tuple List: Field Name:Any',
            'NewNCVarName':'Field Name',
            }
        
    # def _SetFxnDesc(self):

    def DependencyNms(self):
        
        rtrn = self._ArgToList('InFieldName')
        rtrn += self._ArgToList('PrecursorFieldNames')
        return rtrn
    
    # def DependencyNms(self):

    def Exec(self,executedObjects):

        src_dimarr = executedObjects[self.ValFromArgByNm('InFieldName')].ExecRslt()
        src_marr = src_dimarr.data.data
        
        new_marr = np.ma.masked_array(
            cp.deepcopy(src_marr.data),
            mask= ~ src_marr.mask
            )

        # Make the dimensioned array object
        self.execRslt = mpncv.NCDimensionedVar(
            dims = cp.deepcopy(src_dimarr.dims),
            name = self._CreateNCVarName(),
            data = mpncv.NCVar(
                data = new_marr,
                dim_nms = src_dimarr.dims.keys()
                )
            )
            
        self._SetMetadata()
        self.execRslt.set_metadata_val('NodeCommand',self.mptCmdStruct['rawCmdStr'])
        executedObjects[self.RsltNm()] = self

    # def Exec(self,executedObjects):

#class InvertMask(_NetCDFUtilParent):

#### Genetic Fuzzy logic operators

class CvtToFuzzyCurveGenetic(CvtToFuzzyCurve):
   
    def __init__(self,mptCmdStruct=None):
        
        super(CvtToFuzzyCurveGenetic, self).__init__(mptCmdStruct)

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
            'RawMin':['Float'],
            'RawMax':['Float'],
            }
        self.fxnDesc['OptArgs'] = {
            'Metadata':'Any',
            'PrecursorFieldNames':['Field Name','Field Name List'],
            }
        
    # _SetFxnDesc(self):

    def Mutate(
        self,
        rawMin=None,
        rawMax=None,
        fzMin=None,
        fzMax=None,
        degree='Limited',
        prob=0.1
        ):

        if rawMin is None: rawMin = self.ValFromArgByNm('RawMin')
        if rawMax is None: rawMax = self.ValFromArgByNm('RawMax')
        if fzMin is None: fzMin = self.fuzzyMin
        if fzMax is None: fzMax = self.fuzzyMax
            
        # one draw to determine if there is a mutation
        if rndm.random() < prob:

            numVals = len(self.ValFromArgByNm('FuzzyValues'))

            if degree == 'Full':  # full mutation of all elements

                rawVals = [ rndm.uniform(rawMin,rawMax) for x in range(numVals)]
                fzVals = [ rndm.uniform(fzMin,fzMax) for x in range(numVals)]

            elif degree == 'Limited':

                # Mutate one value and limit the severity of the
                # mutation by combining a new value with the
                # current value

                ndxToMutate = rndm.randint(0,numVals-1)
                mutationWt = rndm.random()
                
                rawMutation = rndm.uniform(rawMin,rawMax)
                fzMutation = rndm.uniform(fzMin,fzMax)
                              
                rawVals = self.ValFromArgByNm('RawValues')
                fzVals = self.ValFromArgByNm('FuzzyValues')
            
                rawVals[ndxToMutate] = mutationWt * rawMutation + (1-mutationWt) * rawVals[ndxToMutate]
                fzVals[ndxToMutate] = mutationWt * fzMutation + (1-mutationWt) * fzVals[ndxToMutate]

            else:

                raise Exception(
                    '{}{}{}{}{}'.format(
                        '\n********************ERROR********************\n',
                        'Invalid specification in command:\n',
                        '  Degree must be one of:\n',
                        '    Limited\n',
                        '    Full\n'
                    ),
                )

            # This sorts the raw values and fuzzy values so that they are in low to high
            # order in the created command. This makes for much easier reading of the
            # command produced. These are short lists, so overhead should be very low

            sortedVals = map(list,zip(*sorted(zip(rawVals,fzVals))))

            self.SetArg('RawValues','[{}]'.format(','.join([str(x) for x in sortedVals[0]])))
            self.SetArg('FuzzyValues','[{}]'.format(','.join([str(x) for x in sortedVals[1]])))

        # if rndm.random() < prob:...elif...else
    # def Mutate(...)

# class CvtToFuzzyCurveGenetic(CvtToFuzzyCurve):

class FuzzyWeightedUnionGenetic(FuzzyWeightedUnion):

    def __init__(self,mptCmdStruct=None):
        
        super(FuzzyWeightedUnionGenetic, self).__init__(mptCmdStruct)

    # def __init__(self,mptCmdStruct=None):
    
    def Mutate(
        self,
        degree='Limited',
        prob=0.1
        ):

        # one draw to determine if there is a mutation
        if rndm.random() < prob:
            
            numVals = len(self.ValFromArgByNm('Weights'))

            if degree == 'Full':  # full mutation of all elements
                
                wts = [ rndm.random() for x in range(numVals)]
                
            elif degree == 'Limited':

                # Mutate one value and limit the severity of the
                # mutation by combining a new value with the
                # current value

                ndxToMutate = rndm.randint(0,numVals-1)
                mutationWt = rndm.random()
                mutation = rndm.random()
                
                wts = self.ValFromArgByNm('Weights')
                              
                wts[ndxToMutate] = mutationWt * mutation + (1-mutationWt) * wts[ndxToMutate]

            else:

                raise Exception(
                    '{}{}{}{}{}'.format(
                        '\n********************ERROR********************\n',
                        'Invalid specification in command:\n',
                        '  Degree must be one of:\n',
                        '    Limited\n',
                        '    Full\n'
                    ),
                )
            
            self.SetArg('Weights','[{}]'.format(','.join([str(x/max(wts)) for x in wts])))

        # if rndm.random() < prob:...elif...else

    # def Mutate(...)
        
# class FuzzyWeightedUnionGenetic(FuzzyWeightedUnion):

class FuzzyWeightedUnionWith0Genetic(FuzzyWeightedUnion):

    """
    Created to add a likelihood of introducing a 0 weight. This is
    useful for scenarios in which fewer values considered adds value
    to the fitness function.
    """

    def __init__(self,mptCmdStruct=None):
        
        super(FuzzyWeightedUnionWith0Genetic, self).__init__(mptCmdStruct)

    # def __init__(self,mptCmdStruct=None):

    def Mutate(
        self,
        degree='Limited',
        prob=0.1
        ):

        # one draw to determine if there is a mutation
        if rndm.random() < prob:
            
            numVals = len(self.ValFromArgByNm('Weights'))

            if degree == 'Full':  # full mutation of all elements
                
                wts = [ rndm.random() for x in range(numVals)]
                
            elif degree == 'Limited':

                # Mutate one value and limit the severity of the
                # mutation by combining a new value with the
                # current value

                ndxToMutate = rndm.randint(0,numVals-1)
                wts = self.ValFromArgByNm('Weights')

                # draw to determine if weight to apply is 0

                if rndm.random() < 0.5: # 1 in 20 chance of getting a 0
                    wts[ndxToMutate] = 0.
                    
                else:
                    mutationWt = rndm.random()
                    mutation = rndm.random()

                    wts[ndxToMutate] = mutationWt * mutation + (1-mutationWt) * wts[ndxToMutate]

            else:

                raise Exception(
                    '{}{}{}{}{}'.format(
                        '\n********************ERROR********************\n',
                        'Invalid specification in command:\n',
                        '  Degree must be one of:\n',
                        '    Limited\n',
                        '    Full\n'
                    ),
                )
            
            self.SetArg('Weights','[{}]'.format(','.join([str(x/max(wts)) for x in wts])))

        # if rndm.random() < prob:...elif...else

    # def Mutate(...)
        
# class FuzzyWeightedUnionWith0Genetic(FuzzyWeightedUnion):

class FuzzySelectedUnionGenetic(FuzzySelectedUnion):

    def __init__(self,mptCmdStruct=None):
        
        super(FuzzySelectedUnionGenetic, self).__init__(mptCmdStruct)

    # def __init__(self,mptCmdStruct=None):
    
    def Mutate(
        self,
        degree='Limited',
        prob=0.1
        ):

        # one draw to determine if there is a mutation
        if rndm.random() < prob:

            numVals = len(self.ValFromArgByNm('InFieldNames'))
            tOrF = self.ValFromArgByNm('TruestOrFalsest')

            if degree == 'Full':  # full mutation of all elements

                torF = 'Truest' if rndm.randint(0,1) else 'Falsest'
                numToConsider = rndm.randint(1,numVals)
                
            elif degree == 'Limited':

                # Shift from Truest from Falsest more likely
                # the greater fraction of all inputs are used
                # in the operation. The change in NumberToConsider
                # is based on the total number.

                numToConsider = self.ValFromArgByNm('NumberToConsider')

                if rndm.random() < 0.4 * numToConsider / numVals:
                    
                    # Mutate TruestOrFalsest
                    torF = 'Truest' if self.ValFromArgByNm('TruestOrFalsest') == 'Truest' else 'Falsest'

                else:

                    torF = self.ValFromArgByNm('TruestOrFalsest')
                    
                    # Mutate NumberToConsider
                    moreOrLess = rndm.sample({-1,1},1)[0]
                    if moreOrLess == -1:
                        # Reduce the number to consider
                        if numToConsider > 1:
                            # Change by at least 1
                            randChg = rndm.randint(0,numToConsider-2) + 1
                            numToConsider = numToConsider - randChg
                    else:
                        if numToConsider < numVals:
                            randChg = rndm.randint(1,numVals-numToConsider)
                            numToConsider = numToConsider + randChg

                # if rndm.random() < numToConsider / numVals:
                
            else:

                raise Exception(
                    '{}{}{}{}{}'.format(
                        '\n********************ERROR********************\n',
                        'Invalid specification in command:\n',
                        '  Degree must be one of:\n',
                        '    Limited\n',
                        '    Full\n'
                    ),
                )
            
            # if degree == 'Full':...elif...else
            
            self.SetArg('TruestOrFalsest',tOrF)
            self.SetArg('NumberToConsider',str(numToConsider))
            
        # if rndm.random() < prob:

    # def Mutate(...)
        
# class FuzzySelectedUnionGenetic(FuzzySelectedUnion):

class GraphConversionCurves(_GraphicsParent):
    ''' Creates a graph of the conversion curve used by a
    CvtToFuzzyCurve command.'''

    def __init__(
        self,
        mpt_cmd_struct=None,
        ):

        super(GraphConversionCurves,self).__init__(mpt_cmd_struct)

    # def __init__(...)

    def _SetFxnDesc(self):
        
        self.fxnDesc['DisplayName'] = 'Graph Conversion Curves'
        self.fxnDesc['ShortDesc'] = 'Line graph of fuzzy conversion curve(s).'
        self.fxnDesc['ReturnType'] = 'Boolean'
        
        self.fxnDesc['ReqArgs'] = {
            'InFieldNames':['Field Name','Field Name List'],
            'OutDirectoryName':'File Name'
            }
        self.fxnDesc['OptArgs'] = {
            'PrecursorFieldNames':['Field Name','Field Name List'],
            'FileNamePrepend':'Any',
            'FileNameAppend':'Any',
            'Title':'Any',
            'XLabel':'Any',
            'YLabel':'Any',
            'XMin':'Float',
            'XMax':'Float',
            'YMin':'Float',
            'YMax':'Float',
            'Color':'Any',
            'LineStyle':'Any',
            }
        
    # _SetFxnDesc(self):

    def DependencyNms(self):
        
        rtrn = self._ArgToList('InFieldNames')
        rtrn += self._ArgToList('PrecursorFieldNames')
        return rtrn
    
    # def DependencyNms(self):

    def Exec(self,executedObjects):

        self._InsureDirExists(self.ArgByNm('OutDirectoryName'))
        
        prepend_str = self.ArgByNm('FileNamePrepend')
        if prepend_str is None: prepend_str = ''
        append_str = self.ArgByNm('FileNameAppend')
        if append_str is None: append_str = ''

        # For each variable render the layer and write it
        for in_fldnm in self._ArgToList('InFieldNames'):

            in_exec_rslt = executedObjects[in_fldnm]

            raw_padding = abs(
                (in_exec_rslt.raw_vals[-1] - in_exec_rslt.raw_vals[0]) / 10
                )
            raw_vals = \
              [in_exec_rslt.raw_vals[0] - raw_padding] + \
              in_exec_rslt.raw_vals + \
              [in_exec_rslt.raw_vals[-1] + raw_padding]
              
            fuzzy_vals = \
                [in_exec_rslt.fuzzy_vals[0]] + \
                in_exec_rslt.fuzzy_vals + \
                [in_exec_rslt.fuzzy_vals[-1]]

            fuzzy_padding = abs(in_exec_rslt.fuzzyMax - in_exec_rslt.fuzzyMin) / 50
            raw_padding = abs(
                (in_exec_rslt.raw_vals[-1] - in_exec_rslt.raw_vals[0]) / 50
                )
            
            if not isinstance (in_exec_rslt,_CvtToFuzzyCurve):
                raise Exception(
                    (3*'{}\n').format(
                        'Error: Input field {} is not the result of a plottable convert to fuzzy function.'.format(
                            in_fldnm
                            ),
                    'File: {}  Line number: {}\n'.format(
                        self.mptCmdStruct['cmdFileNm'],
                        self.mptCmdStruct['lineNo']
                        ),
                    'Full command:\n{}\n'.format(self.mptCmdStruct['rawCmdStr'])
                ),
                    )
                    
            out_file_nm = '{}{}{}.png'.format(self.ValFromArgByNm('OutDirectoryName'),prepend_str,in_fldnm,append_str)

            title = self.ValFromArgByNm('Title')
            if title is None:
                title = 'Fuzzy Conversion Curve for {}'.format(in_fldnm)
            x_label = self.ValFromArgByNm('XLabel')
            if x_label is None:
                x_label = 'Raw value'
            y_label = self.ValFromArgByNm('YLabel')
            if y_label is None:
                y_label = 'Fuzzy value'
            x_min = self.ValFromArgByNm('XMin')
            if x_min is None:
                x_min = raw_vals[0] - raw_padding
            x_max = self.ValFromArgByNm('XMax')
            if x_max is None:
                x_max = raw_vals[-1] + raw_padding
            y_min = self.ValFromArgByNm('YMin')
            if y_min is None:
                y_min = in_exec_rslt.fuzzyMin - fuzzy_padding
            y_max = self.ValFromArgByNm('YMax')
            if y_max is None:
                y_max = in_exec_rslt.fuzzyMax + fuzzy_padding
            color = self.ValFromArgByNm('Color')
            if color is None:
                color = 'black'
            line_style = self.ValFromArgByNm('LineStyle')
            if line_style is None:
                line_style = '-'

            self._XYLineGraph(
                raw_vals,
                fuzzy_vals,
                out_file_nm,
                title = title,
                x_lbl = x_label,
                y_lbl = y_label,
                x_min = x_min,
                x_max = x_max,
                y_min = y_min,
                y_max = y_max,
                color = color,
                line_style = line_style
                )

        # for in_fldnm in self._ArgToList('InFieldNames'):
        self.execRslt = True
        executedObjects[self.RsltNm] = self
        
    # def Exec(self,executedObjects):
        
# class GraphConversionCurves(_GraphicsParent):

################################################################################

class MaskNCVarByDimensionIndexRanges(_NetCDFUtilParent):

    '''
    Masks the variable to the extent of the specified
    index ranges. e.g. lat:11:15 would trim the
    extent to cells to the lat[11:16] (Note: it is inclusive)
    '''
    
    def __init__(
        self,
        mpt_cmd_struct=None,
        ):

        super(MaskNCVarByDimensionIndexRanges,self).__init__(mpt_cmd_struct)
        
    # def __init__(...)

    def _SetFxnDesc(self):
        
        self.fxnDesc['DisplayName'] = 'MaskByIndexRanges'
        self.fxnDesc['ShortDesc'] = 'Masks array to extent indexes (not values).'
        
        self.fxnDesc['ReturnType'] = 'NCVar'
        
        self.fxnDesc['ReqArgs'] = {
            'InFieldName':'Field Name',
            'Extents':'Tuple List: Field Name:Integer:Integer',
            'MaskInOrOut':'One of| In Out',
            }
        self.fxnDesc['OptArgs'] = {
            'PrecursorFieldNames':['Field Name','Field Name List'],
            'Metadata':'Tuple List: Field Name:Any',
            'NewNCVarName':'Field Name',
            }
        
    # def _SetFxnDesc(self):

    def DependencyNms(self):
        return self._ArgToList('InFieldName') + \
          self._ArgToList('PrecursorFieldNames')


    def Exec(self,executedObjects):

        comp_tolerance = self.ValFromArgByNm('ComparisonTolerance')
        if comp_tolerance is None: comp_tolerance = 0
            
        to_mask_ncdimvar = executedObjects[self.ValFromArgByNm('InFieldName')].ExecRslt()
        to_mask_marr = to_mask_ncdimvar.data.data

        # Get the ranges for the dimensions
        # used to order slices to match order of ncvar indices
        add_mask_slices = OrderedDict()
        for dim_nm,dim_len in zip(to_mask_ncdimvar.dim_nms,to_mask_ncdimvar.shape):
            add_mask_slices[dim_nm] = slice(0,dim_len)

        used_dim_nms = []
        for tup in self._ArgToList('Extents'):
            
            # strip white space and split
            dim_nm,ndx_val1,ndx_val2 =''.join(tup.split()).split(':')
            dim_min_ndx_val = min(int(ndx_val1),int(ndx_val2))
            dim_max_ndx_val = max(int(ndx_val1),int(ndx_val2)) + 1

            if dim_nm not in to_mask_ncdimvar.dims:
                raise Exception(
                    (3*'{}\n').format(
                        'Error: Non-existent dimension specified: {}'.format(dim_nm),
                        '  Result name: {}'.format(self.RsltNm()),
                        '  Dimension name: {}'.format(dim_nm)
                        )
                    )

            # Error if dim specified twice
            if dim_nm in used_dim_nms:

                raise Exception(
                    (3*'{}\n').format(
                        'Error: Each dimension may only be specified one time',
                        '  Result name: {}'.format(self.RsltNm()),
                        '  Dimension name: {}'.format(dim_nm)
                        )
                    )
            else:
                used_dim_nms.append(dim_nm)

            add_mask_slices[dim_nm] = slice(dim_min_ndx_val,dim_max_ndx_val)

        # for tup in self._ArgToList('Extents'):

        add_mask = np.full_like(to_mask_marr.mask,True)
        add_mask[add_mask_slices.values()] = False
        if self.ValFromArgByNm('MaskInOrOut') == 'Out':
            add_mask = np.where(add_mask == True, False, True)

        self.execRslt = cp.deepcopy(to_mask_ncdimvar)
        self.execRslt.name = self.RsltNm()
        self.execRslt.data.data.mask = np.ma.mask_or(
            self.execRslt.data.data.mask,
            add_mask
            )

        self._SetMetadata()
        self.execRslt.set_metadata_val('NodeCommand',self.mptCmdStruct['rawCmdStr'])
        executedObjects[self.RsltNm()] = self
            
    # def Exec(self,executedObjects):

# class MaskNCVarByDimensionIndexRanges(_NetCDFUtilParent):

class CompositeArrays(_NetCDFUtilParent):

    ''' For an ordered set of arrays, create a composite of what is
    masked in. For a given cell, the value of the first masked in
    of the input arrays is used. For example:
    
    x = masked out, 1 = value in first array, 2 = value in second array

    First array, second array, composite result:

    xxxx     xxxx             xxxx
    x11x     xxxx             x11x
    x11x     x22x             x11x
    xxxx     x22x             x22x
    xxxx     xxxx             xxxx
    '''
    
    def __init__(
        self,
        mpt_cmd_struct=None,
        ):

        super(CompositeArrays,self).__init__(mpt_cmd_struct)

    def _SetFxnDesc(self):
        
        self.fxnDesc['DisplayName'] = 'Composite Arrays'
        self.fxnDesc['ShortDesc'] = 'Takes an ordered set of NCVars and composites them into one array.'
        
        self.fxnDesc['ReturnType'] = 'NCVar'
        
        self.fxnDesc['ReqArgs'] = {
            'InFieldNames':['Field Name List'],
            }
        self.fxnDesc['OptArgs'] = {
            'PrecursorFieldNames':['Field Name','Field Name List'],
            'Metadata':'Tuple List: Field Name:Any',
            }
        
    # def _SetFxnDesc(self):

    def DependencyNms(self):
        
        rtrn = self._ArgToList('InFieldNames')
        rtrn += self._ArgToList('PrecursorFieldNames')
        return rtrn
    
    # def DependencyNms(self):

    def Exec(self,executedObjects):
        # executedObjects is a dictionary of executed MPilot function objects
        
        fldnm_list = self._ArgToList('InFieldNames')
        new_ncdimvar = cp.deepcopy(executedObjects[fldnm_list[0]].ExecRslt())
        new_marr = new_ncdimvar.data.data
        new_ncdimvar.name = self.RsltNm()
                
        for fldnm in fldnm_list[1:]:

            marr_to_add = executedObjects[fldnm].ExecRslt().data.data
            
            if marr_to_add.shape != new_marr.shape:
                fldnms_and_shapes = ''.join(
                    ['  {} {}\n'.format(fldnm_tmp,executedObjects[fldnm_tmp].ExecRslt().data.data.shape) for fldnm_tmp in fldnm_list]
                    )
                    
                raise Exception(
                    (5*'{}\n').format(
                        '\n********************ERROR********************',
                        'All fields must have same shape:',
                        'Fields and shapes: {}'.format(fldnms_and_shapes),
                        'File: {}  Line number: {}'.format(
                            self.mptCmdStruct['cmdFileNm'],
                            self.mptCmdStruct['lineNo']
                            ),
                        'Full command:{}'.format(self.mptCmdStruct['rawCmdStr'])
                    ),
                )


            # Add to the result array only those cells that are unmasked in the current array
            # but not in the result array.

            tmp_where = np.logical_and(new_marr.mask,np.logical_not(marr_to_add.mask))
            np.copyto(new_marr.mask,False, where = tmp_where)
            np.copyto(new_marr.data,marr_to_add.data, where = tmp_where)

        self._SetMetadata()
        self.execRslt = new_ncdimvar
        self.execRslt.set_metadata_val('NodeCommand',self.mptCmdStruct['rawCmdStr'])
        executedObjects[self.RsltNm()] = self
        
    # def Exec(self,executedObjects):

# class CompositeArrays(_NetCDFUtilParent)

class ReadCSVColumnAsNCVariable(_NetCDFUtilParent):

    def __init__(
        self,
        mpt_cmd_struct=None,
        ):

        super(ReadCSVColumnAsNCVariable,self).__init__(mpt_cmd_struct)
        
    # def __init__(...)
    
    def _SetFxnDesc(self):
        
        self.fxnDesc['DisplayName'] = 'Read a csv column into a netCDF variable.'
        self.fxnDesc['ShortDesc'] = 'Read a column from a csv file into a netCDF variable.'
        self.fxnDesc['ReturnType'] = 'NCVar'
        
        self.fxnDesc['ReqArgs'] = {
            'InFieldName':'Field Name',
            'IndexFieldName':'Field Name',
            'InFileName':'File Name',
            }
        self.fxnDesc['OptArgs'] = {
            'MissingValue':'Float',
            'PrecursorFieldNames':['Field Name','Field Name List'],
            'Metadata':'Tuple List: Field Name:Any',
            'NewNCVarName':'Field Name',
            }
                    
    # def _SetFxnDesc(self):

    def DependencyNms(self):
        return self._ArgToList('PrecursorFieldNames')

    def Exec(self,executedObjects):

        if not os.path.isfile(self.ValFromArgByNm('InFileName')):
            raise Exception(
                (4*'{}\n').format(
                    '\n********************ERROR********************',
                    'Input file does not exist: {}'.format(self.ArgByNm('InFileName')),
                    'Script File: {}  Line number: {}'.format(
                        self.mptCmdStruct['cmdFileNm'],
                        self.mptCmdStruct['lineNo']
                        ),
                    'Full command:\n{}'.format(self.mptCmdStruct['rawCmdStr'])
                ),
            )
        # if not os.path.isfile(self.ValFromArgByNm('InFileName'),'r'):

        field_name = self.ValFromArgByNm('InFieldName')
        index_name = self.ValFromArgByNm('IndexFieldName')

        with open(self.ValFromArgByNm('InFileName'),'rb') as in_f:

            rdr = csv.reader(in_f)
            heads = rdr.next()

            if field_name not in heads:
                raise Exception(
                    (5*'{}\n').format(
                        '\n********************ERROR********************',
                        'Read failure for file: {}'.format(self.ArgByNm('InFileName')),
                        '  Variable not in file: {}'.format(self.ArgByNm('InFieldName')),
                        'Script File: {}  Line number: {}'.format(
                            self.mptCmdStruct['cmdFileNm'],
                            self.mptCmdStruct['lineNo']
                            ),
                        'Full command:\n{}'.format(self.mptCmdStruct['rawCmdStr'])
                    ),
                )

            if index_name not in heads:
                raise Exception(
                    (5*'{}\n').format(
                        '\n********************ERROR********************',
                        'Read failure for file: {}'.format(self.ArgByNm('InFileName')),
                        '  Index variable not in file: {}'.format(self.ArgByNm('IndexFieldName')),
                        'Script File: {}  Line number: {}'.format(
                            self.mptCmdStruct['cmdFileNm'],
                            self.mptCmdStruct['lineNo']
                            ),
                        'Full command:\n{}'.format(self.mptCmdStruct['rawCmdStr'])
                    ),
                )

            new_var_val_col_ndx = heads.index(field_name)
            new_var_ndx_col_ndx = heads.index(index_name)

            # Note: this should have checking for bad values, but doesn't
            
            new_var_val = []
            new_var_ndx = []
            for row in rdr:
                if ''.join(row) == '':
                    break
                new_var_val.append(float(row[new_var_val_col_ndx]))
                new_var_ndx.append(float(row[new_var_ndx_col_ndx]))

        # with open(self._ArgToList('InFieldName'),'rb') as in_f:

        new_var_ndx = np.ma.masked_array(new_var_ndx,mask = False)

        missing_val = self.ValFromArgByNm('MissingValue')
        if  missing_val is not None:
            new_var_arr = np.ma.array(new_var_val,dtype = float,fill_value=missing_val)
            new_var_arr = np.ma.masked_where(new_var_arr == missing_val,new_var_val,copy=False)
        else:
            new_var_arr = np.ma.masked_array(new_var_val,mask = False)
            
        dims = OrderedDict()
        dims[index_name] = mpncv.NCVar(
            name = index_name,
            data = np.ma.masked_array(new_var_ndx, dtype = float),
            dim_nms = index_name
            )

        self.execRslt = mpncv.NCDimensionedVar(
            dims = cp.deepcopy(dims),
            name = self._CreateNCVarName(),
            data = mpncv.NCVar(
                data = new_var_arr,
                dim_nms = dims.keys()
                )
            )

        self._SetMetadata()
        self.execRslt.set_metadata_val('NodeCommand',self.mptCmdStruct['rawCmdStr'])
        executedObjects[self.RsltNm()] = self

    # def Exec(self,executedObjects):
    
# class ReadCSVColumnAsNCVariable(_NetCDFUtilParent):

class WriteNCVariablesToCSV(_NetCDFUtilParent):

    # Does not write index

    def __init__(
        self,
        mpt_cmd_struct=None,
        ):

        super(WriteNCVariablesToCSV,self).__init__(mpt_cmd_struct)
        
    # def __init__(...)
    
    def _SetFxnDesc(self):
        
        self.fxnDesc['DisplayName'] = 'Write one or more netCDF variables to a CSV file.'
        self.fxnDesc['ShortDesc'] = 'Write one or more netCDF variables to a CSV file.'
        self.fxnDesc['ReturnType'] = 'Boolean'
        
        self.fxnDesc['ReqArgs'] = {
            'OutFieldNames':['Field Name','Field Name List'],
            'OutFileName':'File Name',
            }
        self.fxnDesc['OptArgs'] = {
            'PrecursorFieldNames':['Field Name','Field Name List'],
            'Metadata':'Tuple List: Field Name:Any',
            }
                    
    # def _SetFxnDesc(self):

    def DependencyNms(self):
        return self._ArgToList('OutFieldNames') + self._ArgToList('PrecursorFieldNames')

    def Exec(self,executedObjects):

        try: # check for writable file
            outF = open(self.ArgByNm('OutFileName'),'w')
            outF.close()
        except IOError:
            raise Exception(
                (4*'{}\n').format(
                    '\n********************ERROR********************',
                    'Unable to open file for writing: {}'.format(self.ArgByNm('OutFileName')),
                    'Script File: {}  Line number: {}'.format(
                        self.mptCmdStruct['cmdFileNm'],
                        self.mptCmdStruct['lineNo']
                        ),
                    'Full command:\n{}'.format(self.mptCmdStruct['rawCmdStr'])
                ),
            )

        # try: # check for writable file

        # assure congruent arrays
        var_0_nm = self._ArgToList('OutFieldNames')[0]
        for outfldnm in self._ArgToList('OutFieldNames')[1:]:
            if not self._NCDimArrsAreCongruent(executedObjects[var_0_nm].ExecRslt(),executedObjects[outfldnm].ExecRslt()):
                raise Exception(
                    (4*'{}\n').format(
                        '\n********************ERROR********************',
                        'Variables do not have the same shape: {}, {}'.format(var_0_nm,outfldnm),
                        'Script File: {}  Line number: {}'.format(
                            self.mptCmdStruct['cmdFileNm'],
                            self.mptCmdStruct['lineNo']
                            ),
                        'Full command:\n{}'.format(self.mptCmdStruct['rawCmdStr'])
                    ),
                )
                
        with open(self.ValFromArgByNm('OutFileName'),'w') as out_f:

            # write the field names
            out_f.write('{}\n'.format(','.join(self._ArgToList('OutFieldNames'))))

            # get the output data arrays and make a common mask
            out_arrs = []
            for outfldnm in self._ArgToList('OutFieldNames'):
                out_arrs.append(executedObjects[outfldnm].ExecRslt().data.data)

            # combine masks
            out_mask = out_arrs[0].mask
            for out_arr in out_arrs[1:]:
                out_mask = np.ma.mask_or(out_mask,out_arr.mask)

            for ndx in range(len(out_arrs)):
                out_arrs[ndx] = np.ma.array(out_arrs[ndx].data,mask=out_mask).compressed()
                
            for row_ndx in range(len(out_arrs[0])):
                out_f.write('{}\n'.format(','.join([str(out_arrs[x][row_ndx]) for x in range(len(out_arrs))])))
                
        self.execRslt = True
        executedObjects[self.RsltNm()] = self

    # def Exec(self,executedObjects):
    
# class WriteNCVariablesToCSV(_NetCDFUtilParent):

class Write1DNCVariablesToCSV(_NetCDFUtilParent):

    def __init__(
        self,
        mpt_cmd_struct=None,
        ):

        super(Write1DNCVariablesToCSV,self).__init__(mpt_cmd_struct)
        
    # def __init__(...)
    
    def _SetFxnDesc(self):
        
        self.fxnDesc['DisplayName'] = 'Write one or more one dimensional netCDF variables to a CSV file.'
        self.fxnDesc['ShortDesc'] = 'Write one or more one dimensional netCDF variables to a CSV file.'
        self.fxnDesc['ReturnType'] = 'Boolean'
        
        self.fxnDesc['ReqArgs'] = {
            'OutFieldNames':['Field Name','Field Name List'],
            'OutFileName':'File Name',
            }
        self.fxnDesc['OptArgs'] = {
            'MissingValue':'Any',
            'PrecursorFieldNames':['Field Name','Field Name List'],
            'Metadata':'Tuple List: Field Name:Any',
            }
                    
    # def _SetFxnDesc(self):

    def DependencyNms(self):
        return self._ArgToList('OutFieldNames') + self._ArgToList('PrecursorFieldNames')

    def Exec(self,executedObjects):

        try: # check for writable file
            outF = open(self.ArgByNm('OutFileName'),'w')
            outF.close()
        except IOError:
            raise Exception(
                (4*'{}\n').format(
                    '\n********************ERROR********************',
                    'Unable to open file for writing: {}'.format(self.ArgByNm('OutFileName')),
                    'Script File: {}  Line number: {}'.format(
                        self.mptCmdStruct['cmdFileNm'],
                        self.mptCmdStruct['lineNo']
                        ),
                    'Full command:\n{}'.format(self.mptCmdStruct['rawCmdStr'])
                ),
            )

        # try: # check for writable file

        # assure 1d congruent arrays
        var_0_nm = self._ArgToList('OutFieldNames')[0]
        var_0 = executedObjects[var_0_nm].ExecRslt()

        if len(var_0.data.data.shape) != 1:
            raise Exception(
                (4*'{}\n').format(
                    '\n********************ERROR********************',
                    'Variables must have only one dimension: {}'.format(var_0_nm,outfldnm),
                    'Script File: {}  Line number: {}'.format(
                        self.mptCmdStruct['cmdFileNm'],
                        self.mptCmdStruct['lineNo']
                        ),
                    'Full command:\n{}'.format(self.mptCmdStruct['rawCmdStr'])
                ),
            )
            
        
        for outfldnm in self._ArgToList('OutFieldNames')[1:]:
            if not self._NCDimArrsAreCongruent(executedObjects[var_0_nm].ExecRslt(),executedObjects[outfldnm].ExecRslt()):
                raise Exception(
                    (4*'{}\n').format(
                        '\n********************ERROR********************',
                        'Variables do not have the same shape: {}, {}'.format(var_0_nm,outfldnm),
                        'Script File: {}  Line number: {}'.format(
                            self.mptCmdStruct['cmdFileNm'],
                            self.mptCmdStruct['lineNo']
                            ),
                        'Full command:\n{}'.format(self.mptCmdStruct['rawCmdStr'])
                    ),
                )

        # Get the dimension array
        ndx_var = var_0.dims.values()[0]
        ndx_arr = ndx_var.data
        
        # print type(var_0.dims.values()[0].data.data)
        # print dir(var_0.dims.values()[0].data.data)
        # print var_0.dims.values()[0].data.data
        # print var_0.dims.values()[0].data.mask
        # print var_0.dims.values()[0].name
        # exit('testing')

        
            
        with open(self.ValFromArgByNm('OutFileName'),'w') as out_f:

            # write the field names
            out_f.write(
                '{},{}\n'.format(
                   ndx_var.name,
                    ','.join(self._ArgToList('OutFieldNames'))
                    )
                )

            # get the output data arrays and make a common mask
            out_arrs = []
            for outfldnm in self._ArgToList('OutFieldNames'):
                out_arrs.append(executedObjects[outfldnm].ExecRslt().data.data)

            # combine masks
            out_mask = out_arrs[0].mask
            for out_arr in out_arrs[1:]:
                out_mask = np.ma.mask_or(out_mask,out_arr.mask)

            for ndx in range(len(out_arrs)):
                out_arrs[ndx] = np.ma.array(out_arrs[ndx].data,mask=out_mask)

            missing_val = self.ArgByNm('MissingValue')

            for row_ndx in range(len(out_arrs[0])):
                outline =  '{},{}\n'.format(
                    ndx_arr[row_ndx],
                    ','.join([str(out_arrs[x][row_ndx]) for x in range(len(out_arrs))])
                    )

                if missing_val is not None:
                    outline = outline.replace('--',missing_val)

                out_f.write(outline)
                
        self.execRslt = True
        executedObjects[self.RsltNm()] = self

    # def Exec(self,executedObjects):
    
# class Write1DNCVariablesToCSV(_NetCDFUtilParent):
