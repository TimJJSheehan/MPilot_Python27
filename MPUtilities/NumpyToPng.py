import png
import numpy as np

class color_ramp(object):

    def __init__(
        self,
        ramp = None,
        reverse_ramp = None,
        r_vals = None,
        g_vals = None,
        b_vals = None,
        bit_depth = 8
        ):

        if ramp is not None:
            self.init_ramp(ramp)
            
        else:
            
            if r_vals == None:
                self.r_vals = range(2**bit_depth)
            else:
                self.r_vals = r_vals

            if g_vals == None:
                self.g_vals = range(2**bit_depth)
            else:
                self.g_vals = g_vals

            if b_vals == None:
                self.b_vals = range(2**bit_depth)
            else:
                self.b_vals = b_vals

            self.__bit_depth = bit_depth
        
        # if color_ramp is not None:...else:...

        if reverse_ramp:
            self.reverse_ramp()

    # def __init__(...)

    def init_ramp(self, ramp_nm):

        if ramp_nm == 'red_to_yellow_01':
            self.r_vals = 256 * [255]
            self.g_vals = range(0,256)
            self.b_vals = 256 * [0]
            self.__bit_depth = 8

        elif ramp_nm == 'red_to_green_01':
            self.r_vals = range(255,-1,-1)
            self.g_vals = range(0,256)
            self.b_vals = 256 * [0]
            self.__bit_depth = 8

        elif ramp_nm == 'yellow_to_green_01':
            self.r_vals = range(255,-1,-1)
            self.g_vals = 256 * [255]
            self.b_vals = 256 * [0]
            self.__bit_depth = 8

        elif ramp_nm == 'red_to_yellow_to_green_01':
            self.r_vals = 128 * [255] + range(255,-1,-2)
            self.g_vals = range(0,256,2) + 128 * [255]
            self.b_vals = 256 * [0]
            self.__bit_depth = 8

        elif ramp_nm == 'red_to_yellow_to_green_02':
            self.r_vals = 128 * [255] + range(255,-1,-2)
            self.g_vals = range(0,256,2) + 128 * [255]
            self.b_vals = 256 * [32]
            self.__bit_depth = 8

        elif ramp_nm == 'red_to_yellow_to_green_03':
            self.r_vals = 128 * [255] + range(255,-1,-2)
            self.g_vals = range(0,256,2) + 128 * [255]
            self.b_vals = 256 * [64]
            self.__bit_depth = 8

        elif ramp_nm == 'magenta_to_lemon_to_lime_01':
            self.r_vals = 128 * [255] + range(255,-1,-2)
            self.g_vals = range(0,256,2) + 128 * [255]
            self.b_vals = 256 * [128]
            self.__bit_depth = 8

            
        elif ramp_nm == 'jet':

            self.r_vals = 96 * [0] +  range(0,256,4) + 64 * [255] + range(255,128,-4)
            self.g_vals = 32 * [0] + range(0,256,4) + 64 * [255] + range(255,-1, -4) + 32 *[0]
            self.b_vals = range(128,255,4) + 64 * [255] + range(255,-1, -4) + 96 * [0]
            self.__bit_depth = 8
            
        elif ramp_nm == 'bwr':
            self.r_vals = range(0,256,2) + 128 * [255]
            self.g_vals = range(0,256,2) + range(255,-1,-2)
            self.b_vals = 128 * [255] + range(255,-1,-2)
            self.__bit_depth = 8
            
        elif ramp_nm == 'black_to_white':
            self.r_vals = range(0,256,1)
            self.g_vals = range(0,256,1)
            self.b_vals = range(0,256,1)
            self.__bit_depth = 8
            
        elif ramp_nm == 'white_to_black':
            self.r_vals = range(255,-1,-1)
            self.g_vals = range(255,-1,-1)
            self.b_vals = range(255,-1,-1)
            self.__bit_depth = 8
            
        else:
            raise Exception(
                (2*'{}\n').format(
                    '********************ERROR********************',
                    'Nonexistent color_ramp specified: {}'.format(ramp_nm)
                    )
                )

    # def init_ramp(self, ramp_nm):

    def reverse_ramp(self):
        self.r_vals.reverse()
        self.g_vals.reverse()
        self.b_vals.reverse()

    @property
    def bit_depth(self):
        return self.__bit_depth

    def r_from_ndx(self,ndx):
        return self.r_vals[ndx]
        
    def g_from_ndx(self,ndx):
        return self.g_vals[ndx]
        
    def b_from_ndx(self,ndx):
        return self.b_vals[ndx]

# class color_ramp(object)

def render_array_to_png(
    arr,
    out_fnm,
    clr_ramp = None,
    value_range = None,
    missing_color = [255,255,255],
    flip_up_down = False,
    flip_left_right = False
    ):

    if clr_ramp is None:
        clr_ramp = color_ramp()  # 8 bit grey scale

    if value_range is None:
        value_range = [arr.min(),arr.max()]

    if missing_color is None:
        missing_color = [255,255,255]

    arr = np.ma.where(arr < value_range[0], value_range[0], arr)
    arr = np.ma.where(arr > value_range[1], value_range[1], arr)

    arr_color_ndx = ((2**clr_ramp.bit_depth-1)*((arr - float(value_range[0]))/(max(value_range) - min(value_range)))).astype(np.uint16)

    # create color array and initialize to background color    
    arr_colored = np.array(np.broadcast_to(arr_color_ndx,(3,arr_color_ndx.shape[0],arr_color_ndx.shape[1])))

    for rgb_ndx in [0,1,2]:
        arr_colored[rgb_ndx,:,:] = missing_color[rgb_ndx]

    if isinstance(arr_color_ndx, np.ma.MaskedArray):

        if arr_color_ndx.mask is True:
            
            # single value True mask. completely masked out
            pass # no valid data to render

        elif arr_color_ndx.mask is False:

            # single value True mask. completely masked in
            
            for row_ndx in range(arr_color_ndx.shape[0]):
                for col_ndx in range(arr_color_ndx.shape[0]):
                    arr_colored[0,row_ndx,col_ndx] = clr_ramp.r_from_ndx(arr_color_ndx[row_ndx,col_ndx])
                    arr_colored[1,row_ndx,col_ndx] = clr_ramp.g_from_ndx(arr_color_ndx[row_ndx,col_ndx])
                    arr_colored[2,row_ndx,col_ndx] = clr_ramp.b_from_ndx(arr_color_ndx[row_ndx,col_ndx])

        else:
            
            # partially masked
            row_ndxs,col_ndxs = np.ma.where(arr_color_ndx.mask == False)
            point_ndxs = zip(row_ndxs,col_ndxs)

            for row_ndx, col_ndx in point_ndxs:
                arr_colored[0,row_ndx,col_ndx] = clr_ramp.r_from_ndx(arr_color_ndx[row_ndx,col_ndx])
                arr_colored[1,row_ndx,col_ndx] = clr_ramp.g_from_ndx(arr_color_ndx[row_ndx,col_ndx])
                arr_colored[2,row_ndx,col_ndx] = clr_ramp.b_from_ndx(arr_color_ndx[row_ndx,col_ndx])

    else:

        # not masked

        for row_ndx in range(arr_color_ndx.shape[0]):
            for col_ndx in range(arr_color_ndx.shape[0]):
                arr_colored[0,row_ndx,col_ndx] = clr_ramp.r_from_ndx(arr_color_ndx[row_ndx,col_ndx])
                arr_colored[1,row_ndx,col_ndx] = clr_ramp.g_from_ndx(arr_color_ndx[row_ndx,col_ndx])
                arr_colored[2,row_ndx,col_ndx] = clr_ramp.b_from_ndx(arr_color_ndx[row_ndx,col_ndx])
        
    # if isinstance(arr_color_ndx, np.ma.MaskedArray):

    if flip_up_down:
        arr_colored = np.flip(arr_colored,axis=1)
    
    if flip_left_right:
        arr_colored = np.flip(arr_colored,axis=2)
    
    with open(out_fnm, 'wb') as f:
        
        writer = png.Writer(height=arr_colored.shape[1], width=arr_colored.shape[2], bitdepth=clr_ramp.bit_depth)
        # Convert array to the Python list of lists expected by the png writer.
        arr_as_list = np.transpose(arr_colored,axes=(1,2,0)).flatten().reshape(arr_colored.shape[1],3*arr_colored.shape[2])

        writer.write(f, arr_as_list)

# def render_array(...)
