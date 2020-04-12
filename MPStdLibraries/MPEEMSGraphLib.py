from __future__ import division
from MPCore import MPilotEEMSFxnParent as mpefp
import numpy as np
import copy as cp
import matplotlib.pyplot as plt
import tempfile as tf
import os
from collections import OrderedDict

class HistoDist(mpefp._MPilotEEMSFxnParent):

    def __init__(self,mptCmdStruct=None):

        super(HistoDist, self).__init__(
            mptCmdStruct=mptCmdStruct,
            dataType='Boolean',
            isDataLayer = False
            )

    # def __init__(self,mptCmdStruct=None):

    def _SetFxnDesc(self):
        
        self.fxnDesc['DisplayName'] = 'HistoDist'
        self.fxnDesc['ShortDesc'] = 'Creates a histogram for a variable\'s distribution'
        self.fxnDesc['ReturnType'] = 'Boolean'
        
        self.fxnDesc['ReqArgs'] = {
            'InFieldName':'Field Name'
            }
        self.fxnDesc['OptArgs'] = {
            'OutFileName':'File Name',
            'Metadata':'Any',
            'PrecursorFieldNames':['Field Name','Field Name List'],
            'XMin':'Float',
            'XMax':'Float',
            'YMin':'Float',
            'YMax':'Float',
            'Color':'Any'
            }
        
    # _SetFxnDesc(self):

    def DependencyNms(self):
        
        rtrn = self._ArgToList('InFieldName')
        rtrn += self._ArgToList('PrecursorFieldNames')
        return rtrn
    
    # def DependencyNms(self):
        
    def Exec(self,executedObjects):
        
        dataObj = executedObjects[self.ValFromArgByNm('InFieldName')]
        self._ValidateIsDataLayer(dataObj)

        outFNm = self.ArgByNm('OutFileName')

        fig = plt.figure()
        ax1 = fig.add_subplot(111)

        xMin = self.ValFromArgByNm('XMin')
        if xMin is None:
            if dataObj.DataType() == 'Fuzzy':
                xMin = self.fuzzyMin
            else:
                xMin = dataObj.ExecRslt().min()

        xMax = self.ValFromArgByNm('XMax')
        if xMax is None:
            if dataObj.DataType() == 'Fuzzy':
                xMax = self.fuzzyMax
            else:
                xMax = dataObj.ExecRslt().max()

        # YMax must be set if we are going to use y limits
        yMax = self.ValFromArgByNm('YMax')
        yMin = self.ValFromArgByNm('yMin')
        if yMax is not None:
            if yMin is None: yMin = 0.
            ax1.set_ylim(yMin,yMax)

        if self.ValFromArgByNm('Color') is not None:
            faceColor = self.ValFromArgByNm('Color')
        else:
            faceColor = 'grey'
                
        hist = plt.hist(
            dataObj.ExecRslt().ravel(),
            bins=50,
            range=(xMin,xMax),
            facecolor = faceColor,
            )
        
        ax1.set_xlabel('Value')
        ax1.set_ylabel('Count')
        ax1.set_title('Distribution for {}'.format(self.ValFromArgByNm('InFieldName')))
        
        if outFNm is None:
            outFNm = tf.mktemp(suffix='.png')
            plt.savefig(outFNm)
            os.system('open -a Preview {}'.format(outFNm))

        else:
            plt.savefig(outFNm)
        
        self.execRslt = True
        
        executedObjects[self.RsltNm()] = self
        
        plt.close(fig) # fig must be closed to get it out of memory
        
    # def Exec(self,executedObjects):
    
# class HistoDist(mpefp._MPilotEEMSFxnParent):

class LineDist(mpefp._MPilotEEMSFxnParent):

    def __init__(self,mptCmdStruct=None):

        super(LineDist, self).__init__(
            mptCmdStruct=mptCmdStruct,
            dataType='Boolean',
            isDataLayer = False
            )

    # def __init__(self,mptCmdStruct=None):

    def _SetFxnDesc(self):
        
        self.fxnDesc['DisplayName'] = 'LineDist'
        self.fxnDesc['ShortDesc'] = 'Creates a connected line segment graph distribution for one or more variables'
        self.fxnDesc['ReturnType'] = 'Boolean'
        
        self.fxnDesc['ReqArgs'] = {
            'InFieldNames':['Field Name','Field Name List']
            }
        self.fxnDesc['OptArgs'] = {
            'OutFileName':'File Name',
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
            'Bins':'Integer'
            }
        
    # _SetFxnDesc(self):

    def DependencyNms(self):
        
        rtrn = self._ArgToList('InFieldNames')
        rtrn += self._ArgToList('PrecursorFieldNames')
        return rtrn
    
    # def DependencyNms(self):
        
    def Exec(self,executedObjects):

        print 'LineDist.Exec()',self.ArgByNm('OutFileName')


        title = self.ValFromArgByNm('Title')
        if title is None: title = 'Distribution'

        xLabel = self.ValFromArgByNm('XLabel')
        if xLabel is None: xLabel = 'Value'

        yLabel = self.ValFromArgByNm('YLabel')
        if yLabel is None: yLabel = 'Count'

        # Get the data objects used to calculate the distrubutions
        dataObjs = OrderedDict()
        for argNm in self._ArgToList('InFieldNames'):
            
            dataObjs[argNm] = executedObjects[argNm]
            self._ValidateIsDataLayer(executedObjects[argNm])

        outFNm = self.ArgByNm('OutFileName')

        fig = plt.figure()
        ax1 = fig.add_subplot(111)

        # Determine mins and maxes
        xMin = self.ValFromArgByNm('XMin')
        if xMin is None:
            xType = 'Fuzzy'
            for argNm,dataObj in dataObjs.items():
                if dataObj.DataType() != 'Fuzzy':
                    xType = 'Float'
                    
            if xType == 'Fuzzy':
                xMin = self.fuzzyMin
            else:
                xMin = float('inf')
                for argNm,dataObj in dataObjs.items():
                    xMin = min(xMin,dataObj.ExecRslt().min())

        # if xMin is None:
        
        xMax = self.ValFromArgByNm('XMax')
        if xMax is None:
            xType = 'Fuzzy'
            for argNm,dataObj in dataObjs.items():
                if dataObj.DataType() != 'Fuzzy':
                    xType = 'Float'
                    
            if xType == 'Fuzzy':
                xMax = self.fuzzyMax
            else:
                xMax = float('-inf')
                for argNm,dataObj in dataObjs.items():
                    xMax = max(xMax,dataObj.ExecRslt().max())

                    
        # YMax must be set if we are going to use y limits
        yMax = self.ValFromArgByNm('YMax')
        yMin = self.ValFromArgByNm('YMin')
        if yMin is None: yMin = 0.

        # Arrange the colors for the lines
        if self.ValFromArgByNm('Colors') is not None:
            
            colorLst = self._ArgToList('Colors')
            
            if len(colorLst) != len(dataObjs):
                raise Exception(
                    '{}{}{}{}{}{}{}{}'.format(
                        '\n********************ERROR********************\n',
                        'Colors must match field names 1 to 1:\n',
                        'FieldNames:\n  ',
                        '{}\n'.format('\n  '.join(dataObjs.keys())),
                        'Colors:\n  ',
                        '{}\n'.format('\n  '.join(colorLst)),
                        'File: {}  Line number: {}\n'.format(
                            self.mptCmdStruct['cmdFileNm'],
                            self.mptCmdStruct['lineNo']
                            ),
                        'Full command:\n{}\n'.format(self.mptCmdStruct['rawCmdStr'])
                    ),
                )

            colors = OrderedDict(zip(dataObjs.keys(),colorLst))

        else:
            colors = OrderedDict(zip(dataObjs.keys(),['black'] * len(dataObjs)))
                
        # if self.ValFromArgByNm('Colors') is not None:

        # Create the graph    
        ax1.set_xlabel(xLabel)
        ax1.set_ylabel(yLabel)
        ax1.set_title(title)
        ax1.set_xlim(xMin,xMax)

        # Now loop throgh the data,
        # make the y points for the graph,
        # and plot the line

        bins = self.ValFromArgByNm('Bins')
        if bins is None:
            bins = 100
        else:
            bins = int(bins)
            
        for dataKey,dataObj in dataObjs.items():

            # first flatten the array, then take only those values that are
            # not masked off, find the centers of the bins for plotting.
            # There is one less bin center than bin edges
            y,edges = np.histogram(dataObj.ExecRslt().flatten().compressed(),bins=bins)
            ctrs = 0.5 * (edges[1:] + edges[:-1])
            ax1.plot(ctrs,y,color=colors[dataKey],label=dataKey)
            
        # for dataKey,dataObj in dataObjs.items():


        if yMax is None:
            yMax = y.max()
        print 'line 10',y.max()
        
        ax1.set_ylim(yMin,yMax)
        
        ax1.legend()
        
        if outFNm is None:
            outFNm = tf.mktemp(suffix='.png')
            plt.savefig(outFNm)
            os.system('open -a Preview {}'.format(outFNm))

        else:
            plt.savefig(outFNm)
        
        self.execRslt = True
        
        executedObjects[self.RsltNm()] = self
        
        plt.close(fig) # fig must be closed to get it out of memory
        
    # def Exec(self,executedObjects):
    
# class HistoDist(mpefp._MPilotEEMSFxnParent):

class RenderLayer(mpefp._MPilotEEMSFxnParent):

    def __init__(self,mptCmdStruct=None):

        super(RenderLayer, self).__init__(
            mptCmdStruct=mptCmdStruct,
            dataType='Boolean',
            isDataLayer = False
            )

    # def __init__(self,mptCmdStruct=None):

    def _SetFxnDesc(self):
        
        self.fxnDesc['DisplayName'] = 'RenderLayer'
        self.fxnDesc['ShortDesc'] = 'Renders a layer'
        self.fxnDesc['ReturnType'] = 'Boolean'
        
        self.fxnDesc['ReqArgs'] = {
            'InFieldName':'Field Name'
            }
        self.fxnDesc['OptArgs'] = {
            'OutFileName':'File Name',
            'Metadata':'Any',
            'PrecursorFieldNames':['Field Name','Field Name List'],
            'ColorMap':'Any',
            'Origin':'Any',
            'FlipX':'Boolean',
            'FlipY':'Boolean',
            'MinVal':'Float',
            'MaxVal':'Float'
            }
        
    # _SetFxnDesc(self):

    def DependencyNms(self):
        
        rtrn = self._ArgToList('InFieldName')
        rtrn += self._ArgToList('PrecursorFieldNames')
        return rtrn
    
    # def DependencyNms(self):
        
    def Exec(self,executedObjects):
        
        dataObj = executedObjects[self.ValFromArgByNm('InFieldName')]
        self._ValidateIsDataLayer(dataObj)

        outFNm = self.ArgByNm('OutFileName')
        
        if self.ArgByNm('FlipX') is not None:
            flipX = self.ValFromArgByNm('FlipX')
        else:
            flipX = False

        if self.ArgByNm('FlipY') is not None:
            flipY = self.ValFromArgByNm('FlipY')
        else:
            flipY = False

        fig = plt.figure(
            # Additions are for stuff around the edges
            figsize = (
                dataObj.ExecRslt().shape[0]/100. + 5.50,
                dataObj.ExecRslt().shape[1]/100. + 0.71
                ),
            dpi = 100
            )
        ax1 = fig.add_subplot(111,facecolor='grey')


        if self.ArgByNm('MinVal') is not None:
            minVal = self.ValFromArgByNm('MinVal')
        elif dataObj.DataType() == 'Fuzzy':
            minVal = self.fuzzyMin
        else:
            minVal = dataObj.ExecRslt().min()
        
        if self.ArgByNm('MaxVal') is not None:
            maxVal = self.ValFromArgByNm('MaxVal')
        elif dataObj.DataType() == 'Fuzzy':
            maxVal = self.fuzzyMax
        else:
            maxVal = dataObj.ExecRslt().max()

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
                labelright=False
                )

        origin = self.ArgByNm('Origin') if self.ArgByNm('Origin') is not None else 'lower'

        tmpRslt = dataObj.ExecRslt()
        if flipX == True:
            tmpRslt = np.flipud(tmpRslt)
        if flipY == True:
            tmpRslt = np.fliplr(tmpRslt)

        myImg = ax1.imshow(
            tmpRslt,
            aspect='auto',
            interpolation='nearest',
            origin=origin
            )
        
        myImg.set_clim(minVal,maxVal)

        cmap = self.ValFromArgByNm('ColorMap')
        if cmap is None: cmap = 'RdYlBu'

        myImg.set_cmap(cmap)

        ax1.set_title(
            '{}'.format(self.ValFromArgByNm('InFieldName')),
            fontsize=18,
            fontweight='bold'
            )

        cbar = fig.colorbar(myImg)
        cbar.ax.tick_params(labelsize=14)
        
        if outFNm is None:
            outFNm = tf.mktemp(suffix='.png')
            plt.savefig(outFNm)
            os.system('open -a Preview {}'.format(outFNm))

        else:
            plt.savefig(outFNm)
        
        self.execRslt = True
        
        executedObjects[self.RsltNm()] = self

        plt.close(fig) # fig must be closed to get it out of memory

    # def Exec(self,executedObjects):
    
# class RenderLayer(mpefp._MPilotEEMSFxnParent):

class ScatterXY(mpefp._MPilotEEMSFxnParent):

    def __init__(self,mptCmdStruct=None):

        super(ScatterXY, self).__init__(
            mptCmdStruct=mptCmdStruct,
            dataType='Boolean',
            isDataLayer = False
            )

    # def __init__(self,mptCmdStruct=None):

    def _SetFxnDesc(self):
        
        self.fxnDesc['DisplayName'] = 'ScatterXY'
        self.fxnDesc['ShortDesc'] = 'Creates a scatter plot for two fields'
        self.fxnDesc['ReturnType'] = 'Boolean'
        
        self.fxnDesc['ReqArgs'] = {
            'XFieldName':'Field Name',
            'YFieldName':'Field Name'
            }
        self.fxnDesc['OptArgs'] = {
            'OutFileName':'File Name',
            'Metadata':'Any',
            'PrecursorFieldNames':['Field Name','Field Name List'],
            }
        
    # _SetFxnDesc(self):

    def DependencyNms(self):
        
        rtrn = self._ArgToList('XFieldName') + self._ArgToList('YFieldName')
        rtrn += self._ArgToList('PrecursorFieldNames')
        return rtrn
    
    # def DependencyNms(self):
        
    def Exec(self,executedObjects):
        
        xDataObj = executedObjects[self.ValFromArgByNm('XFieldName')]
        yDataObj = executedObjects[self.ValFromArgByNm('YFieldName')]
        self._ValidateIsDataLayer(xDataObj)
        self._ValidateIsDataLayer(yDataObj)

        outFNm = self.ArgByNm('OutFileName')

        fig = plt.figure()
        ax1 = fig.add_subplot(111)

        if xDataObj.DataType() == 'Fuzzy':
            xMinVal = self.fuzzyMin
            xMaxVal = self.fuzzyMax
        else:
            xMinVal = xDataObj.ExecRslt().min()
            xMaxVal = xDataObj.ExecRslt().max()

        if yDataObj.DataType() == 'Fuzzy':
            yMinVal = self.fuzzyMin
            yMaxVal = self.fuzzyMax
        else:
            yMinVal = yDataObj.ExecRslt().min()
            yMaxVal = yDataObj.ExecRslt().max()

        fig = plt.figure()
        ax1 = fig.add_subplot(111)

        ax1.set_xlabel(self.ValFromArgByNm('XFieldName'))
        ax1.set_ylabel(self.ValFromArgByNm('YFieldName'))

        ax1.set_xlim(xMinVal,xMaxVal)
        ax1.set_ylim(yMinVal,yMaxVal)
        
        ax1.set_title('{} vs {}'.format(
            self.ValFromArgByNm('YFieldName'),
            self.ValFromArgByNm('XFieldName')
            )
            )

        myWhere = np.logical_or(
            xDataObj.ExecRslt().mask == False,
            yDataObj.ExecRslt().mask == False
            )

        # make the scatter plot
        plt.scatter(
            xDataObj.ExecRslt()[myWhere].ravel(),
            yDataObj.ExecRslt()[myWhere].ravel(),
            c='r',
            lw=0,
            marker='o',
            s=1,
            alpha=0.1
            )
                    
        if outFNm is None:
            outFNm = tf.mktemp(suffix='.png')
            plt.savefig(outFNm)
            os.system('open -a Preview {}'.format(outFNm))

        else:
            plt.savefig(outFNm)
        
        self.execRslt = True
        
        executedObjects[self.RsltNm()] = self
        
        plt.close(fig) # fig must be closed to get it out of memory
        
    # def Exec(self,executedObjects):
    
# class ScatterXY(mpefp._MPilotEEMSFxnParent):
