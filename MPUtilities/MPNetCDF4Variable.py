'''
MPNetCDF4Variable.py

This file contains the classes NCVar, which corresponds to
a NetCDF dataset variable and NCDimensionedVar which contains a
data variable along with variables for its dimensions.

The strategy for this is to use an NCDimensionedVar for each data layer.
Grouping the dimensions with the variable data allows for dimension
trimming and resizing when a variable is trimmed. It also allows for
writing a variable to a NetCDF file along with its dimensions.

The downside, of course is the overhead of having identical dimensions
tied to multiple data variables. But it is worth the overhead to make
the whole system run smoothly when combining and outputting variables.

History:

2017.10.11 - tjs - original version working
2017.10.12 - Split into NCVar and NCIndexedVar
'''

from collections import OrderedDict
import numpy as np
import netCDF4 as nc
import copy as cp

# Utility to find fill value
def lookup_fill_val(ndx_val):

    # fill_vars_lu = {
    #     'NC_FILL_BYTE':-127,
    #     'NC_FILL_CHAR':0,
    #     'NC_FILL_SHORT':-32767,
    #     'NC_FILL_INT':-2147483647L,
    #     'NC_FILL_FLOAT':9.9692099683868690e+36,
    #     'NC_FILL_DOUBLE':9.9692099683868690e+36,
    #     'b':-127,
    #     'c':0,
    #     's':-32767,
    #     'h':-32767,
    #     'f':9.9692099683868690e+36,
    #     'd':9.9692099683868690e+36,
    #     'i':-2147483647L,
    #     'l':-2147483647L,
    #     'int8':-127,
    #     'int16':-32767,
    #     'float32':9.9692099683868690e+36,
    #     'float64':9.9692099683868690e+36,
    #     'int32':-2147483647L
    #     }

    fill_vars_lu = {
        'NC_FILL_BYTE':nc.default_fillvals['S1'],
        'NC_FILL_CHAR':nc.default_fillvals['S1'],
        'NC_FILL_SHORT':nc.default_fillvals['i2'],
        'NC_FILL_INT':nc.default_fillvals['i4'],
        'NC_FILL_FLOAT':nc.default_fillvals['f4'],
        'NC_FILL_DOUBLE':nc.default_fillvals['f8'],

        'l':nc.default_fillvals['i8'],  # int64
        'f':nc.default_fillvals['f4'],  # float32
        'L':nc.default_fillvals['u8'],  # uint64
        'b':nc.default_fillvals['i1'],  # int8
        'I':nc.default_fillvals['u4'],  # uint32
        'S':nc.default_fillvals['S1'],  # |S1
        'h':nc.default_fillvals['i2'],  # int16
        'B':nc.default_fillvals['u1'],  # uint8
        'i':nc.default_fillvals['i4'],  # int32
        'H':nc.default_fillvals['u2'],  # uint16
        'd':nc.default_fillvals['f8'],  # float64

        'int64':nc.default_fillvals['i8'],
        'float32':nc.default_fillvals['f4'],
        'uint64':nc.default_fillvals['u8'],
        'int8':nc.default_fillvals['i1'],
        'uint32':nc.default_fillvals['u4'],
        '|S1':nc.default_fillvals['S1'],
        'int16':nc.default_fillvals['i2'],
        'uint8':nc.default_fillvals['u1'],
        'int32':nc.default_fillvals['i4'],
        'uint16':nc.default_fillvals['u2'],
        'float64':nc.default_fillvals['f8'],        
        }
        
    try:
        return fill_vars_lu[ndx_val]
    except KeyError:
        return None

# def lookup_fill_val(ndx_val):

def masked_arr_from_arr(in_arr):

    # Makes a masked array, sets the fill_value,
    # and sets the actual value of masked array elements
    # to the fill_value

    # Set the data
    if isinstance(in_arr, np.ma.MaskedArray):
        rtrn_arr = np.ma.MaskedArray(
            cp.deepcopy(in_arr),
            fill_value = lookup_fill_val(in_arr.dtype.char)
            )
    else:
        rtrn_arr = np.ma.MaskedArray(
            cp.deepcopy(in_arr),
            mask = False,
            fill_value = lookup_fill_val(in_arr.dtype.char)
            )

    # make sure masked values have the fill_value
    if rtrn_arr.shape != ():
        if isinstance(rtrn_arr.mask,np.bool_):
            if rtrn_arr.mask == True:
                rtrn_arr.data[:] = rtrn_arr.fill_value
        # elif rtrn_arr.mask != False:
        else:
            rtrn_arr.data[:] = np.where(
                rtrn_arr.mask[:] == True,
                rtrn_arr.fill_value,
                rtrn_arr.data[:]
                )
            
    return rtrn_arr

# def masked_arr_from_arr(in_arr):

def _arrays_are_similar(arr_1,arr_2,tolerance=0,check_mask=True):

    if arr_1.shape != arr_2.shape:
        return False

    if abs((arr_1 - arr_2).max()) > tolerance:
        return False

    if check_mask:
        if isinstance(arr_1,np.ma.core.MaskedArray) != isinstance(arr_2,np.ma.core.MaskedArray):
            return False
        if isinstance(arr_1,np.ma.core.MaskedArray):
            if arr_1.mask != arr_2.mask:
                return False

    return True

# def _arrays_are_similar(arr_1,arr_2,tolerance=0,check_mask=True):
    
class NCVar(object):

    def __init__(
        self,
        
        # Use these to create from scratch:
        data = None,      # numpy array
        metadata = None,  # OrderedDict
        dim_nms = None,   # list of index names
        
        # Use this to create from NetCDF Dataset
        ncvar = None,  # A netcdf variable
        
        # Name is optional with ncvar, required with numpy array
        name = None,      # string
        ):

        # If building from variable in a dataset
        if ncvar != None:
                        
            self.init_from_ncvar(ncvar)
            
            if name == None:
                self.name = ncvar.name
            else:
                self.name = name
                
        else:

            # Building from scratch
            self.dim_nms = dim_nms
            self.metadata = metadata
            self.data = data
            self.name = name

        # if ncvar != None:

    # def __init__(self):

    def __enter__(self):
        return(self)

    def __exit__(self,exc_type,exc_value,traceback):
        if exc_type is not None:
            print exc_type, exc_value, traceback
            
    @property
    def dim_nms(self):
        return self._dim_nms

    @dim_nms.setter
    def dim_nms(self,nms):
        if nms == None:
            self._dim_nms = None
        elif isinstance(nms,basestring):
            self._dim_nms = [nms]
        elif isinstance(nms,(list,tuple)):
            self._dim_nms = []
            for nm in nms:
                if isinstance(nm,basestring):
                    self._dim_nms.append(nm)
                else:
                    raise Exception('string or list of strings required required')
        else:
            raise Exception('string or list of strings required required')

    # def dim_nms(self,dim_nms):

    @property
    def metadata(self):
        return self._metadata

    @metadata.setter
    def metadata(self,md_ord_dict):
        if md_ord_dict is None:
            self._metadata = OrderedDict()
        elif isinstance(md_ord_dict,OrderedDict):
            self._metadata = cp.deepcopy(md_ord_dict)
        else:
            raise Exception('OrderedDict required')

    @property
    def name(self):
        return self._name

    @name.setter
    def name(self,new_nm):
        self._name = new_nm

    def set_metadata_val(self,var_nm,var_val):
        if var_val is None:
            if var_nm in self._metadata.keys():
                del self._metadata[var_nm]
        else:
            self._metadata[var_nm] = var_val
        
    def del_metadata_val(self,var_nm):
        self.set_metadata_val(var_nm,None)
            
    @property
    def data(self):
        return self._data

    @property
    def shape(self):
        if isinstance(self.data,np.ndarray):
            return self.data.shape
        else:
            return ()

    @data.setter
    def data(self,nc_data):

        if isinstance(nc_data,np.ndarray):
            self._data = masked_arr_from_arr(nc_data)
        else:
            self._data = cp.deepcopy(nc_data)

    # def data(self,nc_arr):

    def init_from_ncvar(self,ncvar):

        # Set the dimensions
        self._dim_nms = cp.deepcopy(ncvar.dimensions)

        # Set the data
        self.data = ncvar[:]

        # Set the metadata
        self.metadata = None
        self.set_metadata_val('missing_value',self.data.fill_value)
        self.set_metadata_val('_FillValue',self.data.fill_value)

        for att_nm in ncvar.ncattrs():
            # keeping scale factor effectively doubles the scaling
            # so discard it.
            if att_nm == 'scale_factor': 
                continue
            elif att_nm not in self.metadata:
                self.set_metadata_val(att_nm,ncvar.getncattr(att_nm))

    # def init_from_ncvar(self,ncvar,nc_ds):

    def add_to_nc_ds(self,nc_ds):

        if self.name in nc_ds.variables:
            raise Exception(
                '{}{}{}'.format(
                    'Error: Variable already in output dataset\n',
                    '  Output file name:{}\n'.format(nc_ds.filepath()),
                    '  Output variable name: {}\n'.format(self.name),
                    )
                )

        # Creates an NetCDF variable in the NetCDF dataset nc_ds
        if '_FillValue' in self.metadata.keys():

            new_ncvar = nc_ds.createVariable(
                self.name,
                self.data.dtype.char,
                self.dim_nms,
                fill_value = self.metadata['_FillValue']
                )
            
        else:

            new_ncvar = nc_ds.createVariable(
                self.name,
                self.data.dtype.char,
                self.dim_nms,
                fill_value = lookup_fill_val(self.data.dtype.char)
                )
            
        # if '_FillValue' in self.metadata.keys():

        setattr(new_ncvar,'missing_value',new_ncvar._FillValue)
        
        new_ncvar[:] = cp.deepcopy(self.data)

        # Copy the attributes over
        for attr_nm in self.metadata.keys():
            if attr_nm not in new_ncvar.ncattrs():
                setattr(new_ncvar,attr_nm,self.metadata[attr_nm])
        # for attr_nm in self.metadata.keys():

    # def add_to_nc_ds(self,nc_ds):

    def set_masked_data_values(self):
        
        if self.data.shape != ():
            self.data.data[:] = np.where(
                self.data.mask[:] == True,
                self.data.fill_value,
                self.data.data[:]
                )

    # def set_masked_data_values(self):

# class NCVar(object):

class NCDimensionedVar(object):

    def __init__(
        self,

        # Use these to create from scratch
        dims = None,       # OrderedDict of 1-D NCVar objects
        data = None,       # NCVar object

        # Use these to create from a variable in a NetCDF dataset
        nc_ds = None,      # NetCDF Dataset object
        ncvar_nm = None,   # Name of variable in NetCDF Dataset object

        # name required to create from scratch, option if not
        name = None
        ):

        if (nc_ds != None and ncvar_nm == None) or \
          (nc_ds == None and ncvar_nm != None):

            raise Exception('Must provide both nc_ds and ncvar_nm to initialize from NetCDF')

        if nc_ds != None:

            # Get the data for the variable data
            if ncvar_nm not in nc_ds.variables:
                raise Exception('Variable {} not in NetCDF Dataset'.format(nc_var_nm))                

            self._data = NCVar(ncvar=nc_ds.variables[ncvar_nm])

            if name == None:
                self.name = ncvar_nm
            else:
                self.name = name

            # Get the data for the variable's dimensions

            self._dims = OrderedDict()
            for dim_nm in self.data.dim_nms:
                
                if dim_nm not in nc_ds.variables.keys():
                    self._dims[dim_nm] = NCVar(
                        np.arange(nc_ds.dimensions[dim_nm].size,dtype=float),
                        metadata = OrderedDict([('units','index')]),
                        dim_nms = [dim_nm]
                        )
                else:
                    self._dims[dim_nm] = NCVar(ncvar=nc_ds.variables[dim_nm])
        
            # for dim_nm in self.data.dim_nms:
            
        else:
            
            self.dims = dims
            self.data = data
            self.name = name

        # if ncvar != None:
            
    # def __init__(self):

    def __enter__(self):
        return(self)

    def __exit__(self,exc_type,exc_value,traceback):
        if exc_type is not None:
            print exc_type, exc_value, traceback

    @property
    def dim_nms(self):
        return self._dims.keys()

    @property
    def dims(self):
        return self._dims

    @dims.setter
    def dims(self,dim_ord_dict):

        if dim_ord_dict is None:
            self._dims = dim_ord_dict
        elif isinstance(dim_ord_dict,OrderedDict):
            for dim_nm,nc_var in dim_ord_dict.items():
                if len(nc_var.shape) != 1:
                    raise Exception('OrderedDict of 1 D NCVar objects required')
            self._dims = dim_ord_dict
        else:
            raise Exception('OrderedDict of 1 D NCVar objects required')
            
    # def dims(self,dim_ord_dict):
    
    @property
    def metadata(self):
        return self.data._metadata

    @metadata.setter
    def metadata(self,md_ord_dict):
        if isinstance(md_ord_dict,OrderedDict) or md_ord_dict is None:
            self.data._metadata = md_ord_dict
        else:
            raise Exception('OrderedDict or None required')
    
    @property
    def data(self):
        return self._data

    @data.setter
    def data(self,nc_var):
        
        if isinstance(nc_var,NCVar) or nc_var is None:
            self._data = nc_var
        else:
            raise Exception('NCVar or None required')

    @property
    def name(self):
        return self._name
                
    @name.setter
    def name(self,new_nm):
        self._name = new_nm
        if self.data is not None:
            self.data.name = new_nm

    @property
    def shape(self):
        rtrn = None
        if self.data is not None:
            rtrn = self.data.shape
        return rtrn
    # def shape(self):

    @property
    def extents(self):
        
        rtrn = OrderedDict()
        for dim_nm,dim_data in self.dims.items():
            rtrn[dim_nm] = (dim_data.data.min(),dim_data.data.max())

        return rtrn

    def get_extent(self,dim_nm):
        
        return (dim_data.data.min(),dim_data.data.max())
    
    def dim_axis_from_name(self,dim_nm):
        
        if dim_nm not in self.dim_nms:
            return None
        else:
            return self.dim_nms.index(dim_nm)
        
    # def dim_ndx_from_name(self,dim_nm):

    def set_dim(self,dim_nm,dim_ncvar):

        if dim_val is None:
            self._dims[dim_nm] = dim_ncvar
        elif not isinstance(dim_ncvar,NCVar):
            raise Exception('OrderedDict of 1 D NCVar objects required')
        
        if len(dim_val.shape != 1):
            raise Exception('OrderedDict of 1 D NCVar objects required')
        self._dims[dim_nm] = dim_val
            
    # def set_dim(self,dim_nm,dim_val):

    def set_metadata_val(self,var_nm,var_val):
        self.data.set_metadata_val(var_nm,var_val)
        
    def dims_are_valid(self):
        
        rtrn = True
        
        if len(self.data.shape) != len(self.dims.keys()):
            
            rtrn = False
            
        else:
            
            for ndx in range(len(self.data.shape)):
                dim_key = self.dims.keys()[ndx]
                if self.data.shape[ndx] != self.dims[dim_key].data.shape[0]:
                    rtrn = False
                    break

        # if len(data.shape) != len(dims.keys()):

        return rtrn

    # def dims_are_valid(self):

    def flip_dim(self,dim_nm):

        # flips the order of values in an axis and the associated dimension
        dim_axis = self.dim_axis_from_name(dim_nm)
        
        if dim_axis is None:
            raise Exception(
                '{}{}{}{}'.format(
                    'Illegal dimension name\n',
                    'Dimension name: {}\n'.format(dim_nm),
                    'Variable name:  {}\n'.format(self.name),
                    )
                )

        self.dims[dim_nm].data = np.flip(self.dims[dim_nm].data,axis=0)
        self.data.data = np.flip(self.data.data,axis=dim_axis)

    # def flip_dim(self,dim_nm):

    def reorder_dim_lo_to_hi(self,dim_nm):

        if dim_nm not in self.dim_nms:
            raise Exception(
                '{}{}{}'.format(
                    'Illegal dimension name\n',
                    'Dimension name: {}\n'.format(dim_nm),
                    'Variable name:  {}\n'.format(self.name),
                    )
                )

        if self.dims[dim_nm].data[0] > self.dims[dim_nm].data[-1]:
            self.flip_dim(dim_nm)

    # def reorder_dim_lo_to_hi(self,dim_nm):

    def reorder_dim_hi_to_lo(self,dim_nm):

        if dim_nm not in self.dim_nms:
            raise Exception(
                '{}{}{}{}'.format(
                    'Illegal dimension name\n',
                    'Dimension name: {}\n'.format(dim_nm),
                    'Variable name:  {}\n'.format(self.name),
                    )
                )

        if self.dims[dim_nm].data[0] < self.dims[dim_nm].data[-1]:
            self.flip_dim(dim_nm)

    # def reorder_dim_hi_to_lo(self,dim_nm):
        
    def add_to_nc_ds(self,nc_ds,tolerance):

        # Creates an NetCDF variable in the NetCDF dataset nc_ds

        # Check the dimensions in nc_ds. If they exist they must be
        # congruent with the variable. If they don't exist,
        # they must be added

        for dim_nm in self.dims:

            # does the dim variable exist in nc_ds
            if dim_nm in nc_ds.variables:

                if not _arrays_are_similar(
                    nc_ds.variables[dim_nm],
                    self.dims[dim_nm].data,
                    tolerance = tolerance,
                    check_mask=False
                    ):
                    # _arrays_are_similar(arr_1,arr_2,tolerance=0,check_mask=True):
                    
                    raise Exception(
                        '{}{}{}{}'.format(
                            'Dimension conflict\n',
                            'Dimension name: {}\n'.format(dim_nm),
                            'Variable name:  {}\n'.format(self.name),
                            'Output dataset: {}\n'.format(nc_ds.filepath())
                            )
                        )


                # ds_dim_var = nc_ds.variables[dim_nm]

                # dim_err = False
                # if isinstance(ds_dim_var,np.ma.MaskedArray):
                    
                #     if not np.array_equal(ds_dim_var[:].data, self.dims[dim_nm].data.data) or \
                #       not np.array_equal(ds_dim_var[:].mask, self.dims[dim_nm].data.mask):

                #         dim_err = True
                      
                # elif not np.array_equal(ds_dim_var[:],self.dims[dim_nm].data):

                #     dim_err = True

                # if dim_err:
                #     raise Exception(
                #         '{}{}{}{}'.format(
                #             'Dimension conflict\n',
                #             'Dimension name: {}\n'.format(dim_nm),
                #             'Variable name:  {}\n'.format(self.name),
                #             'Output dataset: {}\n'.format(nc_ds.filepath())
                #             )
                #         )

                # # if isinstance(ds_dim_var,np.ma.MaskedArray):

                # Add the dimension to dimensions?
                if dim_nm not in nc_ds.dimensions:
                    nc_ds.createDimension(dim_nm, self.dims[dimNm].shape[0])

            else: # dim_nm is not in variables

                # Add the dimension to dimensions?
                if dim_nm not in nc_ds.dimensions:
                    nc_ds.createDimension(dim_nm, self.dims[dim_nm].shape[0])

                # Add the dimension to the variables

                self.dims[dim_nm].add_to_nc_ds(nc_ds)

            # if dim_nm in nc_ds.variables:
        # for dim_nm in self.dims.keys():

        # Dimensions are taken care of. Now its time to add the data.

        self.data.add_to_nc_ds(nc_ds)

    # def create_ncvar(self,nc_ds):

# class NCDimensionedVar(object):

class Extents(OrderedDict):

    '''
    A set of dimensional extents not associated with a particular variable.
    Useful for defining a minimal or maximal set of extents over multiple
    variables.

    ExtentName : (min,max)
    '''

    def __init__(
        self,
        dim_extents=None,
        ):

        if dim_extents is None:
            
            # Empty is o.k.
            super(Extent,self).__init__()
            
        else:
            
            try:
                self._extents = cp.deepcopy(dim_extents)
            except:
                raise Exception(
                    (2*'{}\n').format(
                        'Error: Extent specification must be an OrderedDict of the form:',
                        '  ExtentName : (min,max)'
                        )
                    ) # raise Exception(...)
                    
       # if dim_extents is None:

    # def __init__(...)

    def __enter__(self):
        return(self)

    def __exit__(self,exc_type,exc_value,traceback):
        if exc_type is not None:
            print exc_type, exc_value, traceback

    def set_extent(self,dim_nm,min,max):
        self[dim_nm] = (min,max)

    def get_extent(self,dim_nm):
        if dim_nm in self:
            return self[dim_nm]
        else:
            return None

    def del_extent(self,dim_nm):
        if dim_nm in self:
            del self[dim_nm]

    def extent_exists(self,dim_nm):
        return dim_nm in self

    @property
    def extents(self):
        return self

# class Extents(object):
    
