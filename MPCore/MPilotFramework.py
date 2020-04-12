import re
from MPCore import MPilotFxnParent as mpfp
import importlib
import inspect
from collections import OrderedDict

class MPilotFramework(object):

    # The MPilot framework manages MPilot libraries
    # and makes the classes representing the MPilot functions
    # avaiable to other programs. It makes available only
    # classes that are descended from _MPilotFxnDefParent
    # and that do not begin with an underscore character.
    #
    # Note: a module specification is either a module directly
    # accessible by python, or a specification of a module's
    # location and the package that contains it, for example, 
    # if there is a package call MathLibs that contains a module
    # called Addition, you would pass in
    # ('.Addition','MathLibs')
    
    def __init__(self, moduleSpecLst):

        if not isinstance(moduleSpecLst,list):
            moduleNmLst = [moduleSpecLst]
        self.moduleSpecLst = moduleSpecLst
        self.pilotFxnNmsByModule = {}
        self.pilotFxnClasses = {}
        for modSpec in self.moduleSpecLst:
            self.AddModule(modSpec)

    # def __init__(self, moduleNmLst):

    def __enter__(self):
        return self

    def __exit__(self,exc_type,exc_value,traceback):
        if exc_type is not None:
            print exc_type, exc_value, traceback

    # def _LoadModules(self):

    ################################################################################
    # Public classes
    ################################################################################
            
    def AddModule(self,modSpec):
            
        # Import each module

        # If the module is part of a package, the pass in a tuple with the
        # module name relative to the package and the package name. For
        # example if there is a package call MathLibs that contains a module
        # called Addition, you would pass in
        # ('.Addition','MathLibs')
        # Add the classes that are PilotFxnDefParent descendants
        # to the dict of PilotFxnClasses unless there is a name
        # conflict


        # print 'importing', modSpec

        if isinstance(modSpec,tuple):
            mod = importlib.import_module(modSpec[0],package = modSpec[1])
            modNm = '{}{}'.format(modSpec[1],modSpec[0])
        else:
            mod = importlib.import_module(modSpec)
            modNm = modSpec

        self.pilotFxnNmsByModule[modNm] = []
            
        for attNm in dir(mod):

            attr = getattr(mod,attNm)

            # We only want:
            #   intended to be public
            #   is a class
            #   descended from mpfp._MPilotFxnParent
            #   defined in the current module
            if attNm.find('_') == 0: continue
            if not inspect.isclass(attr): continue
            if not issubclass(attr,mpfp._MPilotFxnParent): continue
            if not getattr(mod,attNm).__module__ == modNm: continue
            
            # Is it not a duplicate of loaded classes
            if attNm in self.pilotFxnClasses:
                raise Exception(
                    '{}{}{}{}'.format(
                        '\n********************ERROR********************\n',
                        'Pilot Function Class defined in two modules.\n',
                        '  Class name: {}\n'.format(attNm),
                        '  Modules: {}, {}\n'.format(
                            self.pilotFxnClasses[attNm]['modNm'],modNm
                            )
                        )
                    )
            self.pilotFxnClasses[attNm] = {'modNm':modNm,'mod':mod}
            self.pilotFxnNmsByModule[modNm].append(attNm)

        # for attNm in dir(mod):
    # def AddModule(self,modNm):
        
    def CreateFxnObject(self,pilotFxnNm,mptCmdStruct=None):

        if pilotFxnNm not in self.pilotFxnClasses:
            raise Exception(
                '{}{}{}'.format(
                    '\n********************ERROR********************\n',
                    'Function not in MPilot framework.\n',
                    '  Function class name: {}\n'.format(pilotFxnNm),
                    )
                )
        
        if mptCmdStruct is None:
            rtrn = getattr(self.pilotFxnClasses[pilotFxnNm]['mod'],pilotFxnNm)()
        else:
            rtrn = getattr(self.pilotFxnClasses[pilotFxnNm]['mod'],pilotFxnNm)(mptCmdStruct)

        return rtrn
    # def CreateFxnObject(self,pilotFxnNm,pilotFxnArgs):

    def GetAllFxnClassInfo(self):
        
        rtrn = []
        for pilotFxnNm in sorted(self.pilotFxnClasses.keys()):
            with getattr(self.pilotFxnClasses[pilotFxnNm]['mod'],pilotFxnNm)() as fxn:
                rtrn.append(fxn.FxnDesc())

        return rtrn
    
    # def GetAllFxnClassInfo(self):

    def GetAllFormattedFxnClassInfo(self):
        
        rtrn = ''
        for pilotFxnNm in sorted(self.pilotFxnClasses.keys()):
            with getattr(self.pilotFxnClasses[pilotFxnNm]['mod'],pilotFxnNm)() as fxn:
                rtrn = '{}\n{}\n'.format(rtrn,fxn.FormattedFxnDesc())

        return rtrn
    
    # def GetAllFxnClassInfo(self,verbose=False):

    def FxnNames(self): return sorted(self.pilotFxnClasses.keys())

    def GetFxnFormattedClassInfo(self,pilotFxnNm):

        rtrn = None
        
        if pilotFxnNm in self.pilotFxnClasses.keys():
            with getattr(self.pilotFxnClasses[pilotFxnNm]['mod'],pilotFxnNm)() as fxn:
                rtrn = fxn.FormattedFxnDesc()

        return rtrn
    
# class MPilotFramework(object):
