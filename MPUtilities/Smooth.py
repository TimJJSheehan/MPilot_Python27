# from pylab import *
# import sys
import numpy as np
# from scipy import stats
# from netCDF4 import *
# from scipy.io import netcdf
# import os.path
# from math import *
# import datetime
# import decimal
# import re
# import matplotlib.pyplot as plt


def SmoothTriangle(data,degree):
    """
    performs moving triangle smoothing with a variable degree.
    at ends, smooths asymmetrically using only as much of the triangle
    numpy.array as possible
    """
    # create the triangle [1,2...degree...2,1]
    triangle=np.array(range(degree+1) + list(reversed(range(degree))))
    smoothed=[]
    for i in range(len(data)):
        triangleRangeStart = max(degree-i, 0)
        triangleRangeEnd = min(degree+len(data)-i,degree*2+1)
        dataRangeStart = max(0,i-degree)
        dataRangeEnd = min(len(data),i+degree+1)
        
        point=data[dataRangeStart:dataRangeEnd]*triangle[triangleRangeStart:triangleRangeEnd]
        smoothed.append(sum(point)/sum(triangle[triangleRangeStart:triangleRangeEnd]))

    return np.array(smoothed)

