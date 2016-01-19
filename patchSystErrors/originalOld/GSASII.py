#!/usr/bin/env python
# -*- coding: utf-8 -*-
#GSASII
########### SVN repository information ###################
# $Date: 2015-12-21 12:03:08 -0500 (Mon, 21 Dec 2015) $
# $Author: toby $
# $Revision: 2101 $
# $URL: https://subversion.xray.aps.anl.gov/pyGSAS/trunk/GSASII.py $
# $Id: GSASII.py 2101 2015-12-21 17:03:08Z toby $
########### SVN repository information ###################
'''
*GSAS-II Main Module*
=====================

Main routines for the GSAS-II program
'''

import os
import sys
import math
import copy
import random as ran
import time
import copy
import glob
import imp
import inspect
import numpy as np
import scipy as sp
import wx
import matplotlib as mpl
try:
    import OpenGL as ogl
except ImportError:
    print('*******************************************************')
    print('PyOpenGL is missing from your python installation')
    print('     - we will try to install it')
    print('*******************************************************')
    def install_with_easyinstall(package):
        try: 
            print "trying a system-wide PyOpenGl install"
            easy_install.main(['-f',os.path.split(__file__)[0],package])
            return
        except:
            pass
        try: 
            print "trying a user level PyOpenGl install"
            easy_install.main(['-f',os.path.split(__file__)[0],'--user',package])
            return
        except:
            print "Install of '+package+' failed. Please report this information:"
            import traceback
            print traceback.format_exc()
            sys.exit()
    from setuptools.command import easy_install
    install_with_easyinstall('PyOpenGl')
    print('*******************************************************')         
    print('OpenGL has been installed. Restarting GSAS-II')
    print('*******************************************************')         
    loc = os.path.dirname(__file__)
    import subprocess
    subprocess.Popen([sys.executable,os.path.join(loc,'GSASII.py')])
    sys.exit()
    
# load the GSAS routines
import GSASIIpath
GSASIIpath.SetVersionNumber("$Revision: 2101 $")
import GSASIIIO as G2IO
import GSASIIgrid as G2gd
import GSASIIctrls as G2G
import GSASIIplot as G2plt
import GSASIIpwd as G2pwd
import GSASIIpwdGUI as G2pdG
import GSASIIphsGUI as G2phsG
import GSASIIimgGUI as G2imG
import GSASIIspc as G2spc
import GSASIIstrMain as G2stMn
import GSASIIstrIO as G2stIO
import GSASIImath as G2mth
import GSASIImapvars as G2mv
import GSASIIobj as G2obj
import GSASIIlattice as G2lat
import GSASIIlog as log

__version__ = '0.2.0'

# PATCH: for Mavericks (OS X 10.9.x), wx produces an annoying warning about LucidaGrandeUI.
# In case stderr has been suppressed there, redirect python error output to stdout. Nobody
# else should care much about this. 
sys.stderr = sys.stdout

def create(parent):
    return GSASII(parent)

def SetDefaultDData(dType,histoName,NShkl=0,NDij=0):
    if dType in ['SXC','SNC']:
        return {'Histogram':histoName,'Show':False,'Scale':[1.0,True],
            'Babinet':{'BabA':[0.0,False],'BabU':[0.0,False]},
            'Extinction':['Lorentzian','None', {'Tbar':0.1,'Cos2TM':0.955,
            'Eg':[1.e-10,False],'Es':[1.e-10,False],'Ep':[1.e-10,False]}],
            'Flack':[0.0,False]}
    elif dType == 'SNT':
        return {'Histogram':histoName,'Show':False,'Scale':[1.0,True],
            'Babinet':{'BabA':[0.0,False],'BabU':[0.0,False]},
            'Extinction':['Lorentzian','None', {
            'Eg':[1.e-10,False],'Es':[1.e-10,False],'Ep':[1.e-10,False]}]}
    elif 'P' in dType:
        return {'Histogram':histoName,'Show':False,'Scale':[1.0,False],
            'Pref.Ori.':['MD',1.0,False,[0,0,1],0,{},[],0.1],
            'Size':['isotropic',[1.,1.,1.],[False,False,False],[0,0,1],
                [1.,1.,1.,0.,0.,0.],6*[False,]],
            'Mustrain':['isotropic',[1000.0,1000.0,1.0],[False,False,False],[0,0,1],
                NShkl*[0.01,],NShkl*[False,]],
            'HStrain':[NDij*[0.0,],NDij*[False,]],                          
            'Extinction':[0.0,False],'Babinet':{'BabA':[0.0,False],'BabU':[0.0,False]}}

class GSASII(wx.Frame):
    '''Define the main GSAS-II frame and its associated menu items
    '''
    def MenuBinding(self,event):
        '''Called when a menu is clicked upon; looks up the binding in table
        '''
        log.InvokeMenuCommand(event.GetId(),self,event)
            
    def Bind(self,eventtype,handler,*args,**kwargs):
        '''Override the Bind function so that we can wrap calls that will be logged.
        
        N.B. This is a bit kludgy. Menu bindings with an id are wrapped and
        menu bindings with an object and no id are not. 
        '''
        if eventtype == wx.EVT_MENU and 'id' in kwargs:
            menulabels = log.SaveMenuCommand(kwargs['id'],self,handler)
            if menulabels:
                wx.Frame.Bind(self,eventtype,self.MenuBinding,*args,**kwargs)
                return
        wx.Frame.Bind(self,eventtype,handler,*args,**kwargs)      
    
    def _Add_FileMenuItems(self, parent):
        item = parent.Append(
            help='Open a GSAS-II project file (*.gpx)', id=wx.ID_ANY,
            kind=wx.ITEM_NORMAL,text='&Open project...')
        self.Bind(wx.EVT_MENU, self.OnFileOpen, id=item.GetId())
        item = parent.Append(
            help='Save project under current name', id=wx.ID_ANY,
            kind=wx.ITEM_NORMAL,text='&Save project')
        self.Bind(wx.EVT_MENU, self.OnFileSave, id=item.GetId())
        item = parent.Append(
            help='Save current project to new file', id=wx.ID_ANY,
            kind=wx.ITEM_NORMAL,text='Save project as...')
        self.Bind(wx.EVT_MENU, self.OnFileSaveas, id=item.GetId())
        item = parent.Append(
            help='Create empty new project, saving current is optional', id=wx.ID_ANY,
            kind=wx.ITEM_NORMAL,text='&New project')
        self.Bind(wx.EVT_MENU, self.OnFileClose, id=item.GetId())
        item = parent.Append(
            help='Exit from GSAS-II', id=wx.ID_ANY,
            kind=wx.ITEM_NORMAL,text='&Exit')
        self.Bind(wx.EVT_MENU, self.OnFileExit, id=item.GetId())
        
    def _Add_DataMenuItems(self,parent):
        item = parent.Append(
            help='',id=wx.ID_ANY,
            kind=wx.ITEM_NORMAL,
            text='Read image data...')
        self.Bind(wx.EVT_MENU, self.OnImageRead, id=item.GetId())
        item = parent.Append(
            help='',id=wx.ID_ANY,
            kind=wx.ITEM_NORMAL,
            text='Read Powder Pattern Peaks...')
        self.Bind(wx.EVT_MENU, self.OnReadPowderPeaks, id=item.GetId())
        item = parent.Append(
            help='',id=wx.ID_ANY,
            kind=wx.ITEM_NORMAL,
            text='Sum powder data')
        self.Bind(wx.EVT_MENU, self.OnPwdrSum, id=item.GetId())
        item = parent.Append(
            help='',id=wx.ID_ANY,
            kind=wx.ITEM_NORMAL,
            text='Sum image data')
        self.Bind(wx.EVT_MENU, self.OnImageSum, id=item.GetId())
        item = parent.Append(
            help='',id=wx.ID_ANY,
            kind=wx.ITEM_NORMAL,
            text='Add new phase')
        self.Bind(wx.EVT_MENU, self.OnAddPhase, id=item.GetId())
        item = parent.Append(
            help='',id=wx.ID_ANY,
            kind=wx.ITEM_NORMAL,
            text='Delete phase')
        self.Bind(wx.EVT_MENU, self.OnDeletePhase, id=item.GetId())
        item = parent.Append(
            help='',id=wx.ID_ANY,
            kind=wx.ITEM_NORMAL,
            text='Rename data') 
        self.Bind(wx.EVT_MENU, self.OnRenameData, id=item.GetId())
        item = parent.Append(
            help='',id=wx.ID_ANY,
            kind=wx.ITEM_NORMAL,
            text='Delete data')
        self.Bind(wx.EVT_MENU, self.OnDataDelete, id=item.GetId())
                
    def _Add_CalculateMenuItems(self,parent):
        item = parent.Append(help='Make new PDFs from selected powder patterns', 
            id=wx.ID_ANY, kind=wx.ITEM_NORMAL,text='Make new PDFs')
        self.MakePDF.append(item)
#        item.Enable(False)
        self.Bind(wx.EVT_MENU, self.OnMakePDFs, id=item.GetId())
        
        item = parent.Append(help='View least squares parameters', 
            id=wx.ID_ANY, kind=wx.ITEM_NORMAL,text='&View LS parms')
        self.Bind(wx.EVT_MENU, self.ShowLSParms, id=item.GetId())
        
        item = parent.Append(help='', id=wx.ID_ANY, kind=wx.ITEM_NORMAL,
            text='&Refine')
        if len(self.Refine): # extend state for new menus to match main (on mac)
            state = self.Refine[0].IsEnabled()
        else:
            state = False
        item.Enable(state)
        self.Refine.append(item)
        self.Bind(wx.EVT_MENU, self.OnRefine, id=item.GetId())
        
        item = parent.Append(help='', id=wx.ID_ANY, kind=wx.ITEM_NORMAL,
            text='Sequential refine')
        if len(self.SeqRefine): # extend state for new menus to match main (on mac)
            state = self.SeqRefine[0].IsEnabled()
        else:
            state = False
        item.Enable(state)
        self.SeqRefine.append(item) # save menu obj for use in self.EnableSeqRefineMenu
        self.Bind(wx.EVT_MENU, self.OnSeqRefine, id=item.GetId())
        
    def _init_Imports(self):
        '''import all the G2phase*.py & G2sfact*.py & G2pwd*.py files that 
        are found in the path
        '''

        self.ImportPhaseReaderlist = []
        self._init_Import_routines('phase',self.ImportPhaseReaderlist,'Phase')
        self.ImportSfactReaderlist = []
        self._init_Import_routines('sfact',self.ImportSfactReaderlist,'Struct_Factor')
        self.ImportPowderReaderlist = []
        self._init_Import_routines('pwd',self.ImportPowderReaderlist,'Powder_Data')
        self.ImportSmallAngleReaderlist = []
        self._init_Import_routines('sad',self.ImportSmallAngleReaderlist,'SmallAngle_Data')
        self.ImportImageReaderlist = []
        self._init_Import_routines('img',self.ImportImageReaderlist,'Images')
        self.ImportMenuId = {}

    def _init_Import_routines(self,prefix,readerlist,errprefix):
        '''import all the import readers matching the prefix
        '''
        #path2GSAS2 = os.path.dirname(os.path.realpath(__file__)) # location of this file
        #pathlist = sys.path[:]
        #if path2GSAS2 not in pathlist: pathlist.append(path2GSAS2)
        #path2GSAS2 = os.path.join(
        #    os.path.dirname(os.path.realpath(__file__)), # location of this file
        #    'imports')
        pathlist = sys.path[:]
        #if path2GSAS2 not in pathlist: pathlist.append(path2GSAS2)
        if '.' not in pathlist: pathlist.append('.') # insert the directory where G2 is started

        filelist = []
        for path in pathlist:
            for filename in glob.iglob(os.path.join(
                path,
                "G2"+prefix+"*.py")):
                filelist.append(filename)    
                #print 'debug: found',filename
        filelist = sorted(list(set(filelist))) # remove duplicates
        for filename in filelist:
            path,rootname = os.path.split(filename)
            pkg = os.path.splitext(rootname)[0]
            try:
                fp = None
                fp, fppath,desc = imp.find_module(pkg,[path,])
                pkg = imp.load_module(pkg,fp,fppath,desc)
                for clss in inspect.getmembers(pkg): # find classes defined in package
                    if clss[0].startswith('_'): continue
                    if inspect.isclass(clss[1]):
                        # check if we have the required methods
                        for m in 'Reader','ExtensionValidator','ContentsValidator':
                            if not hasattr(clss[1],m): break
                            if not callable(getattr(clss[1],m)): break
                        else:
                            reader = clss[1]() # create an import instance
                            if reader.UseReader:
                                readerlist.append(reader)
            except AttributeError:
                print 'Import_'+errprefix+': Attribute Error '+str(filename)
            #except ImportError:
            #    print 'Import_'+errprefix+': Error importing file '+str(filename)
            except Exception,errmsg:
                print('\nImport_'+errprefix+': Error importing file '+str(filename))
                print('Error message: '+str(errmsg)+'\n')
            if fp: fp.close()

    def EnableSeqRefineMenu(self):
        '''Enable or disable the sequential refinement menu items based on the
        contents of the Controls 'Seq Data' item (if present)
        '''
        controls = self.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(self,self.root, 'Controls'))
        if controls.get('Seq Data'):
            for i in self.SeqRefine: i.Enable(True)
        else:
            for i in self.SeqRefine: i.Enable(False)

    def PreviewFile(self,filename,fp):
        'confirm we have the right file'
        rdmsg = 'File '+str(filename)+' begins:\n\n'
        for i in range(3):
            rdmsg += fp.readline()
        rdmsg += '\n\nDo you want to read this file?'
        if not all([ord(c) < 128 and ord(c) != 0 for c in rdmsg]): # show only if ASCII
            rdmsg = 'File '+str(
                filename)+' is a binary file. Do you want to read this file?'
        # it would be better to use something that
        # would resize better, but this will do for now
        dlg = wx.MessageDialog(
            self, rdmsg,
            'Is this the file you want?', 
            wx.YES_NO | wx.ICON_QUESTION,
            )
        dlg.SetSize((700,300)) # does not resize on Mac
        result = wx.ID_NO
        try:
            result = dlg.ShowModal()
        finally:
            dlg.Destroy()
        if result == wx.ID_NO: return True
        return False
    
    def OnImportGeneric(self,reader,readerlist,label,multiple=False,
                        usedRanIdList=[],Preview=True):
        '''Used for all imports, including Phases, datasets, images...

        Called from :meth:`GSASII.OnImportPhase`, :meth:`GSASII.OnImportImage`,
        :meth:`GSASII.OnImportSfact`, :meth:`GSASII.OnImportPowder` and
        :meth:`GSASII.OnImportSmallAngle`

        Uses reader_objects subclassed from :class:`GSASIIIO.ImportPhase`,
        :class:`GSASIIIO.ImportStructFactor`,
        :class:`GSASIIIO.ImportPowderData`,
        :class:`GSASIIIO.ImportSmallAngleData` or
        :class:`GSASIIIO.ImportImage`.
        If a specific reader is specified, only that method will be called,
        but if no reader is specified, every one that is potentially
        compatible (by file extension) will be tried on the file(s)
        selected in the Open File dialog.

        :param reader_object reader: This will be a reference to
          a particular object to be used to read a file or None,
          if every appropriate reader should be used.

        :param list readerlist: a list of reader objects appropriate for
          the current read attempt. At present, this will be either
          self.ImportPhaseReaderlist, self.ImportSfactReaderlist
          self.ImportPowderReaderlist or self.ImportImageReaderlist
          (defined in _init_Imports from the files found in the path),
          but in theory this list could be tailored.
          Used only when reader is None.

        :param str label: string to place on the open file dialog:
          Open `label` input file

        :param bool multiple: True if multiple files can be selected
          in the file dialog. False is default. At present True is used
          only for reading of powder data.

        :param list usedRanIdList: an optional list of random Ids that 
          have been used and should not be reused

        :param bool Preview: indicates if a preview of the file should
          be shown. Default is True, but set to False for image files
          which are all binary. 

        :returns: a list of reader objects (rd_list) that were able
          to read the specified file(s). This list may be empty. 
        '''
        self.lastimport = ''
        self.zipfile = None
        singlereader = True
        if reader is None:
            singlereader = False
            multiple = False
            #print "use all formats"
            choices = "any file (*.*)|*.*"
            choices += "|zip archive (.zip)|*.zip"
            extdict = {}
            # compile a list of allowed extensions
            for rd in readerlist:
                fmt = rd.formatName
                for extn in rd.extensionlist:
                    if not extdict.get(extn): extdict[extn] = []
                    extdict[extn] += [fmt,]
            for extn in sorted(extdict.keys(),cmp=lambda x,y: cmp(x.lower(), y.lower())):
                fmt = ''
                for f in extdict[extn]:
                    if fmt != "": fmt += ', '
                    fmt += f
                choices += "|" + fmt + " file (*" + extn + ")|*" + extn
        else:
            readerlist = [reader,]
            # compile a list of allowed extensions
            choices = reader.formatName + " file ("
            w = ""
            for extn in reader.extensionlist:
                if w != "": w += ";"
                w += "*" + extn
            choices += w + ")|" + w
            choices += "|zip archive (.zip)|*.zip"
            if not reader.strictExtension:
                choices += "|any file (*.*)|*.*"
        # get the file(s) (probably should remove wx.CHANGE_DIR here)
        if multiple:
            mode = style=wx.OPEN | wx.CHANGE_DIR | wx.MULTIPLE
        else:
            mode = style=wx.OPEN | wx.CHANGE_DIR
        filelist = self.GetImportFile(message="Choose "+label+" input file",
                    defaultFile="",wildcard=choices, style=mode)
        rd_list = []
        filelist1 = []
        for filename in filelist:
            # is this a zip file?
            if os.path.splitext(filename)[1].lower() == '.zip':
                extractedfiles = G2IO.ExtractFileFromZip(
                    filename,parent=self,
                    multipleselect=True)
                if extractedfiles is None: continue # error or Cancel 
                if extractedfiles != filename:
                    self.zipfile = filename # save zip name
                    filelist1 += extractedfiles
                    continue
            filelist1.append(filename)
        filelist = filelist1
        for filename in filelist:
            # is this a zip file?
            if os.path.splitext(filename)[1].lower() == '.zip':
                extractedfile = G2IO.ExtractFileFromZip(filename,parent=self)
                if extractedfile is None: continue # error or Cancel 
                if extractedfile != filename:
                    filename,self.zipfile = extractedfile,filename # now use the file that was created
            # determine which formats are compatible with this file
            primaryReaders = []
            secondaryReaders = []
            for rd in readerlist:
                flag = rd.ExtensionValidator(filename)
                if flag is None: 
                    secondaryReaders.append(rd)
                elif flag:
                    primaryReaders.append(rd)
            if len(secondaryReaders) + len(primaryReaders) == 0 and reader:
                self.ErrorDialog('Not supported','The selected reader cannot read file '+filename)
                return []
            elif len(secondaryReaders) + len(primaryReaders) == 0:
                self.ErrorDialog('No Format','No matching format for file '+filename)
                return []

            fp = None
            msg = ''
            fp = open(filename,'Ur')
            if len(filelist) == 1 and Preview:
                if self.PreviewFile(filename,fp): return []
            self.lastimport = filename # this is probably not what I want to do -- it saves only the
            # last name in a series. See rd.readfilename for a better name.

            # try the file first with Readers that specify the
            # file's extension and later with ones that merely allow it
            errorReport = ''
            for rd in primaryReaders+secondaryReaders:
                rd.ReInitialize() # purge anything from a previous read
                fp.seek(0)  # rewind
                rd.errors = "" # clear out any old errors
                if not rd.ContentsValidator(fp): # rejected on cursory check
                    errorReport += "\n  "+rd.formatName + ' validator error'
                    if rd.errors: 
                        errorReport += ': '+rd.errors
                    continue 
                repeat = True
                rdbuffer = {} # create temporary storage for file reader
                block = 0
                while repeat: # loop if the reader asks for another pass on the file
                    block += 1
                    repeat = False
                    fp.seek(0)  # rewind
                    rd.objname = os.path.basename(filename)
                    flag = False
                    if GSASIIpath.GetConfigValue('debug'):
                        flag = rd.Reader(filename,fp,self,
                                         buffer=rdbuffer,
                                         blocknum=block,
                                         usedRanIdList=usedRanIdList,
                                         )
                    else:
                        try:
                            flag = rd.Reader(filename,fp,self,
                                             buffer=rdbuffer,
                                             blocknum=block,
                                             usedRanIdList=usedRanIdList,
                                             )
                        except rd.ImportException as detail:
                            rd.errors += "\n  Read exception: "+str(detail)
                        except Exception as detail:
                            import traceback
                            rd.errors += "\n  Unhandled read exception: "+str(detail)
                            rd.errors += "\n  Traceback info:\n"+str(traceback.format_exc())
                    if flag: # this read succeeded
                        rd.readfilename = filename
                        rd_list.append(copy.deepcopy(rd)) # save the result before it is written over
                        if rd.repeat:
                            repeat = True
                        continue
                    errorReport += '\n'+rd.formatName + ' read error'
                    if rd.errors:
                        errorReport += ': '+rd.errors
                if rd_list: # read succeeded, was there a warning or any errors? 
                    if rd.warnings:
                        self.ErrorDialog('Read Warning','The '+ rd.formatName+
                                         ' reader reported a warning message:\n\n'+
                                         rd.warnings)
                    break # success in reading, try no further
            else:
                if singlereader:
                    self.ErrorDialog('Read Error','The '+ rd.formatName+
                                     ' reader was not able to read file '+filename+msg+
                                     '\n\nError message(s):\n'+errorReport
                                     )
                else:
                    self.ErrorDialog('Read Error','No reader was able to read file '+filename+msg+
                                     '\n\nError messages:\n'+errorReport
                                     )
            if fp: fp.close()
        return rd_list

    def _Add_ImportMenu_Phase(self,parent):
        '''configure the Import Phase menus accord to the readers found in _init_Imports
        '''
        submenu = wx.Menu()
        item = parent.AppendMenu(wx.ID_ANY, 'Phase',
            submenu, help='Import phase data')
        for reader in self.ImportPhaseReaderlist:
            item = submenu.Append(wx.ID_ANY,help=reader.longFormatName,
                kind=wx.ITEM_NORMAL,text='from '+reader.formatName+' file')
            self.ImportMenuId[item.GetId()] = reader
            self.Bind(wx.EVT_MENU, self.OnImportPhase, id=item.GetId())
        item = submenu.Append(wx.ID_ANY,
                              help='Import phase data, use file to try to determine format',
                              kind=wx.ITEM_NORMAL,
                              text='guess format from file')
        self.Bind(wx.EVT_MENU, self.OnImportPhase, id=item.GetId())

    def GetImportFile(self, message, defaultDir="", defaultFile="", style=wx.OPEN,
                      *args, **kwargs):
        '''Defines a customized dialog that gets gets files from the appropriate
        Import directory (unless overridden with defaultDir). The default
        import directory is found in self.ImportDir, which will be set to
        the location of the (first) file entered in this dialog.
           
        :returns: a list of files or an empty list
        '''
        dlg = wx.FileDialog(self, message, defaultDir, defaultFile, *args,
                            style=style, **kwargs)
        if not defaultDir and self.ImportDir:
            dlg.SetDirectory(self.ImportDir)
        try:
            if dlg.ShowModal() == wx.ID_OK:
                if style & wx.MULTIPLE:
                    filelist = dlg.GetPaths()
                    if len(filelist) == 0: return []
                else:
                    filelist = [dlg.GetPath(),]
                # not sure if we want to do this (why use wx.CHANGE_DIR?)
                if style & wx.CHANGE_DIR: # to get Mac/Linux to change directory like windows!
                    os.chdir(dlg.GetDirectory())
                if not defaultDir and self.ImportDir: # if we are using a default and
                    # user moved, save the new location
                    self.ImportDir = os.path.split(filelist[0])[0]
                    #print 'default moved to',self.ImportDir               
            else: # cancel was pressed
                return []
        finally:
            dlg.Destroy()
        return filelist            
        
    def OnImportPhase(self,event):
        '''Called in response to an Import/Phase/... menu item
        to read phase information.
        dict self.ImportMenuId is used to look up the specific
        reader item associated with the menu item, which will be
        None for the last menu item, which is the "guess" option
        where all appropriate formats will be tried.
        '''
        # look up which format was requested
        reqrdr = self.ImportMenuId.get(event.GetId())
        
        # make a list of phase names, ranId's and the histograms used in those phases
        phaseRIdList,usedHistograms = self.GetPhaseInfofromTree()
        phaseNameList = usedHistograms.keys() # phase names in use
        usedHKLFhists = [] # used single-crystal histograms
        for p in usedHistograms:
            for h in usedHistograms[p]:
                if h.startswith('HKLF ') and h not in usedHKLFhists:
                    usedHKLFhists.append(h)
        rdlist = self.OnImportGeneric(reqrdr,
                                  self.ImportPhaseReaderlist,
                                  'phase',usedRanIdList=phaseRIdList)
        if len(rdlist) == 0: return
        # for now rdlist is only expected to have one element
        # but below will allow multiple phases to be imported
        # if ever the import routines ever implement multiple phase reads. 
        self.CheckNotebook()
        newPhaseList = []
        for rd in rdlist:
            PhaseName = ''
            dlg = wx.TextEntryDialog( # allow editing of phase name
                                    self, 'Enter the name for the new phase',
                                    'Edit phase name', rd.Phase['General']['Name'],
                                    style=wx.OK)
            while PhaseName == '':
                dlg.CenterOnParent()
                if dlg.ShowModal() == wx.ID_OK:
                    PhaseName = dlg.GetValue().strip()
                else:
                    dlg.Destroy()
                    return
            dlg.Destroy()
            # make new phase names unique
            rd.Phase['General']['Name'] = G2obj.MakeUniqueLabel(PhaseName,phaseNameList)
            PhaseName = rd.Phase['General']['Name'][:]
            newPhaseList.append(PhaseName)
            print 'Read phase '+str(PhaseName)+' from file '+str(self.lastimport)
            if not G2gd.GetPatternTreeItemId(self,self.root,'Phases'):
                sub = self.PatternTree.AppendItem(parent=self.root,text='Phases')
            else:
                sub = G2gd.GetPatternTreeItemId(self,self.root,'Phases')
            psub = self.PatternTree.AppendItem(parent=sub,text=PhaseName)
            self.PatternTree.SetItemPyData(psub,rd.Phase)
            self.PatternTree.Expand(self.root) # make sure phases are seen
            self.PatternTree.Expand(sub) 
            self.PatternTree.Expand(psub)
            self.PatternTree.SelectItem(psub) # show the page to complete the initialization (yuk!)
            wx.Yield() # make sure call of GSASII.OnPatternTreeSelChanged happens before we go on

            if rd.Constraints:
                sub = G2gd.GetPatternTreeItemId(self,self.root,'Constraints') # was created in CheckNotebook if needed
                Constraints = self.PatternTree.GetItemPyData(sub)                
                # TODO: make sure that NEWVAR names are unique here?
                for i in rd.Constraints:
                    if type(i) is dict:
                        #for j in i: print j,' --> ',i[j]
                        if '_Explain' not in Constraints: Constraints['_Explain'] = {}
                        Constraints['_Explain'].update(i)
                        continue
                    Constraints['Phase'].append(i)
        if not newPhaseList: return # somehow, no new phases
        # get a list of existing histograms
        PWDRlist = []
        HKLFlist = []
        if self.PatternTree.GetCount():
            item, cookie = self.PatternTree.GetFirstChild(self.root)
            while item:
                name = self.PatternTree.GetItemText(item)
                if name.startswith('PWDR ') and name not in PWDRlist:
                    PWDRlist.append(name)
                if name.startswith('HKLF ') and name not in HKLFlist:
                    HKLFlist.append(name)
                item, cookie = self.PatternTree.GetNextChild(self.root, cookie)
        TextList = PWDRlist + HKLFlist
        if not len(TextList): return # no data loaded yet
        header = 'Select histogram(s) to add to new phase(s):'
        for phaseName in newPhaseList:
            header += '\n  '+str(phaseName)

        notOK = True
        while notOK:
            result = G2G.ItemSelector(TextList,self,header,header='Add histogram(s)',multiple=True)
            if not result: return
            # check that selected single crystal histograms are not already in use!
            used = [TextList[i] for i in result if TextList[i] in usedHKLFhists]
            #for i in result:
            #    if TextList[i] in usedHKLFhists: used.append(TextList[i])
            if used:
                msg = 'The following single crystal histogram(s) are already in use'
                for i in used:
                    msg += '\n  '+str(i)
                msg += '\nAre you sure you want to add them to this phase? '
                msg += 'Associating a single crystal dataset to >1 histogram is usually an error, '
                msg += 'so No is suggested here.'
                if self.ErrorDialog('Likely error',msg,self,wtype=wx.YES_NO) == wx.ID_YES: notOK = False
            else:
                notOK = False
        # connect new phases to histograms
        sub = G2gd.GetPatternTreeItemId(self,self.root,'Phases')
        if not sub:
            raise Exception('ERROR -- why are there no phases here?')
        wx.BeginBusyCursor()
        item, cookie = self.PatternTree.GetFirstChild(sub)
        while item: # loop over (new) phases
            phaseName = self.PatternTree.GetItemText(item)
            data = self.PatternTree.GetItemPyData(item)
            item, cookie = self.PatternTree.GetNextChild(sub, cookie)
            if phaseName not in newPhaseList: continue
            generalData = data['General']
            SGData = generalData['SGData']
            Super = generalData.get('Super',0)
            SuperVec = []
            if Super:
                SuperVec = np.array(generalData['SuperVec'][0])
            UseList = data['Histograms']
            NShkl = len(G2spc.MustrainNames(SGData))
            NDij = len(G2spc.HStrainNames(SGData))
            for i in result:
                histoName = TextList[i]
                if histoName in HKLFlist:
                    #redo UpdateHKLFdata(histoName) here:
                    Id = G2gd.GetPatternTreeItemId(self,self.root,histoName)
                    refDict,reflData = self.PatternTree.GetItemPyData(Id)
                    G,g = G2lat.cell2Gmat(generalData['Cell'][1:7])
                    Super = reflData.get('Super',0)
                    for iref,ref in enumerate(reflData['RefList']):
                        hkl = ref[:3]
                        if Super:
                            H = list(hkl+SuperVec*ref[3])
                        else:
                            H = hkl
                        ref[4+Super] = np.sqrt(1./G2lat.calc_rDsq2(H,G))
                        iabsnt = G2spc.GenHKLf(H,SGData)[0]
                        if iabsnt:  #flag space gp. absences
                            if Super: 
                                if not ref[2+Super]:
                                    ref[3+Super] = 0
                                else:
                                    ref[3+Super] = 1    #twin id?
                            else:
                                ref[3] = 0
                    UseList[histoName] = SetDefaultDData(reflData['Type'],histoName)
                elif histoName in PWDRlist:
                    Id = G2gd.GetPatternTreeItemId(self,self.root,histoName)
                    refList = self.PatternTree.GetItemPyData(
                        G2gd.GetPatternTreeItemId(self,Id,'Reflection Lists'))
                    refList[generalData['Name']] = {}
                    UseList[histoName] = SetDefaultDData('PWDR',histoName,NShkl=NShkl,NDij=NDij)
                else:
                    raise Exception('Unexpected histogram '+str(histoName))
        wx.EndBusyCursor()
        return # success
        
    def _Add_ImportMenu_Image(self,parent):
        '''configure the Import Image menus accord to the readers found in _init_Imports
        '''
        submenu = wx.Menu()
        item = parent.AppendMenu(wx.ID_ANY, 'Image',
            submenu, help='Import image file')
        for reader in self.ImportImageReaderlist:
            item = submenu.Append(wx.ID_ANY,help=reader.longFormatName,
                kind=wx.ITEM_NORMAL,text='from '+reader.formatName+' file')
            self.ImportMenuId[item.GetId()] = reader
            self.Bind(wx.EVT_MENU, self.OnImportImage, id=item.GetId())
        item = submenu.Append(wx.ID_ANY,
                              help='Import image data, use file to try to determine format',
                              kind=wx.ITEM_NORMAL,
                              text='guess format from file')
        self.Bind(wx.EVT_MENU, self.OnImportImage, id=item.GetId())
        
    def OnImportImage(self,event):
        '''Called in response to an Import/Image/... menu item
        to read an image from a file. Like all the other imports,
        dict self.ImportMenuId is used to look up the specific
        reader item associated with the menu item, which will be
        None for the last menu item, which is the "guess" option
        where all appropriate formats will be tried.

        A reader object is filled each time an image is read. 
        '''
        # look up which format was requested
        reqrdr = self.ImportMenuId.get(event.GetId())
        rdlist = self.OnImportGeneric(reqrdr,
                                  self.ImportImageReaderlist,
                                  'image',multiple=True,Preview=False)
        if not rdlist: return
        first = True
        for rd in rdlist:
            if first:
                first = False
                self.CheckNotebook()
            if rd.repeatcount == 1 and not rd.repeat: # skip image number if only one in set
                rd.Data['ImageTag'] = None
            else:
                rd.Data['ImageTag'] = rd.repeatcount
            G2IO.LoadImage2Tree(rd.readfilename,self,rd.Comments,rd.Data,rd.Npix,rd.Image)
        self.PatternTree.SelectItem(G2gd.GetPatternTreeItemId(self,self.Image,'Image Controls'))             #show last image to have beeen read
                    
    def _Add_ImportMenu_Sfact(self,parent):
        '''configure the Import Structure Factor menus accord to the readers found in _init_Imports
        '''
        submenu = wx.Menu()
        item = parent.AppendMenu(wx.ID_ANY, 'Structure Factor',
            submenu, help='Import Structure Factor data')
        for reader in self.ImportSfactReaderlist:
            item = submenu.Append(wx.ID_ANY,help=reader.longFormatName,                
                kind=wx.ITEM_NORMAL,text='from '+reader.formatName+' file')
            self.ImportMenuId[item.GetId()] = reader
            self.Bind(wx.EVT_MENU, self.OnImportSfact, id=item.GetId())
        item = submenu.Append(wx.ID_ANY,
            help='Import Structure Factor, use file to try to determine format',
            kind=wx.ITEM_NORMAL,
            text='guess format from file')
        self.Bind(wx.EVT_MENU, self.OnImportSfact, id=item.GetId())

    def OnImportSfact(self,event):
        '''Called in response to an Import/Structure Factor/... menu item
        to read single crystal datasets. 
        dict self.ImportMenuId is used to look up the specific
        reader item associated with the menu item, which will be
        None for the last menu item, which is the "guess" option
        where all appropriate formats will be tried.
        '''
        # get a list of existing histograms
        HKLFlist = []
        if self.PatternTree.GetCount():
            item, cookie = self.PatternTree.GetFirstChild(self.root)
            while item:
                name = self.PatternTree.GetItemText(item)
                if name.startswith('HKLF ') and name not in HKLFlist:
                    HKLFlist.append(name)
                item, cookie = self.PatternTree.GetNextChild(self.root, cookie)
        # look up which format was requested
        reqrdr = self.ImportMenuId.get(event.GetId())
        rdlist = self.OnImportGeneric(reqrdr,self.ImportSfactReaderlist,
            'Structure Factor',multiple=True)
        if len(rdlist) == 0: return
        self.CheckNotebook()
        newHistList = []
        for rd in rdlist:
            HistName = rd.objname
            if len(rdlist) <= 2: 
                dlg = wx.TextEntryDialog( # allow editing of Structure Factor name
                    self, 'Enter the name for the new Structure Factor',
                    'Edit Structure Factor name', HistName,
                    style=wx.OK)
                dlg.CenterOnParent()
                if dlg.ShowModal() == wx.ID_OK:
                    HistName = dlg.GetValue()
                dlg.Destroy()
            HistName = 'HKLF '+HistName
            # make new histogram names unique
            if len(rd.Banks):
                for Bank in rd.Banks:
                    valuesdict = {'wtFactor':1.0,'Dummy':False,'ranId':ran.randint(0,sys.maxint),}
                    HistName = G2obj.MakeUniqueLabel(HistName,HKLFlist)
                    print 'Read structure factor table '+str(HistName)+' from file '+str(self.lastimport)
                    Id = self.PatternTree.AppendItem(parent=self.root,text=HistName)
                    if not Bank['RefDict'].get('FF'):
                        Bank['RefDict']['FF'] = {}
                    self.PatternTree.SetItemPyData(Id,[valuesdict,Bank['RefDict']])
                    Sub = self.PatternTree.AppendItem(Id,text='Instrument Parameters')
                    self.PatternTree.SetItemPyData(Sub,copy.copy(rd.Parameters))
                    self.PatternTree.SetItemPyData(
                        self.PatternTree.AppendItem(Id,text='Reflection List'),{})  #dummy entry for GUI use
                    newHistList.append(HistName)
            else:
                valuesdict = {'wtFactor':1.0,'Dummy':False,'ranId':ran.randint(0,sys.maxint),}
                HistName = G2obj.MakeUniqueLabel(HistName,HKLFlist)
                print 'Read structure factor table '+str(HistName)+' from file '+str(self.lastimport)
                if not rd.RefDict.get('FF'):
                    rd.RefDict['FF'] = {}
                Id = self.PatternTree.AppendItem(parent=self.root,text=HistName)
                self.PatternTree.SetItemPyData(Id,[valuesdict,rd.RefDict])
                Sub = self.PatternTree.AppendItem(Id,text='Instrument Parameters')
                self.PatternTree.SetItemPyData(Sub,rd.Parameters)
                self.PatternTree.SetItemPyData(
                    self.PatternTree.AppendItem(Id,text='Reflection List'),{})  #dummy entry for GUI use
                newHistList.append(HistName)
                
            self.PatternTree.SelectItem(Id)
            self.PatternTree.Expand(Id)
            self.Sngl = True

        if not newHistList: return # somehow, no new histograms
        # make a list of phase names
        phaseRIdList,usedHistograms = self.GetPhaseInfofromTree()
        phaseNameList = usedHistograms.keys() # phase names in use
        if not phaseNameList: return # no phases yet, nothing to do
        header = 'Select phase(s) to add the new\nsingle crystal dataset(s) to:'
        for Name in newHistList:
            header += '\n  '+str(Name)
        result = G2G.ItemSelector(phaseNameList,self,header,header='Add to phase(s)',multiple=True)
        if not result: return
        # connect new phases to histograms
        sub = G2gd.GetPatternTreeItemId(self,self.root,'Phases')
        if not sub:
            raise Exception('ERROR -- why are there no phases here?')
        wx.BeginBusyCursor()
        item, cookie = self.PatternTree.GetFirstChild(sub)
        iph = -1
        while item: # loop over (new) phases
            iph += 1
            phaseName = self.PatternTree.GetItemText(item)
            data = self.PatternTree.GetItemPyData(item)
            item, cookie = self.PatternTree.GetNextChild(sub, cookie)
            if iph not in result: continue
            generalData = data['General']
            SGData = generalData['SGData']
            Super = generalData.get('Super',0)
            SuperVec = []
            if Super:
                SuperVec = np.array(generalData['SuperVec'][0])
            UseList = data['Histograms']
            for histoName in newHistList:
                #redo UpdateHKLFdata(histoName) here:
                Id = G2gd.GetPatternTreeItemId(self,self.root,histoName)
                refDict,reflData = self.PatternTree.GetItemPyData(Id)
                UseList[histoName] = SetDefaultDData(reflData['Type'],histoName)
                G,g = G2lat.cell2Gmat(generalData['Cell'][1:7])
                if 'TwMax' in reflData:     #nonmerohedral twins present
                    UseList[histoName]['Twins'] = []
                    for iT in range(reflData['TwMax'][0]+1):
                        if iT in reflData['TwMax'][1]:
                            UseList[histoName]['Twins'].append([False,0.0])
                        else:
                            UseList[histoName]['Twins'].append([np.array([[1,0,0],[0,1,0],[0,0,1]]),[1.0,False,reflData['TwMax'][0]]])
                else:   #no nonmerohedral twins
                    UseList[histoName]['Twins'] = [[np.array([[1,0,0],[0,1,0],[0,0,1]]),[1.0,False,0]],]
                for iref,ref in enumerate(reflData['RefList']):
                    hkl = ref[:3]
                    if Super:
                        H = list(hkl+SuperVec*ref[3])
                    else:
                        H = hkl
                    ref[4+Super] = np.sqrt(1./G2lat.calc_rDsq2(H,G))
                    iabsnt,mul,Uniq,phi = G2spc.GenHKLf(H,SGData)
                    if iabsnt:  #flag space gp. absences
                        if Super: 
                            if not ref[2+Super]:
                                ref[3+Super] = 0
                            else:
                                ref[3+Super] = 1    #twin id?
                        else:
                            ref[3] = 0
        wx.EndBusyCursor()
        
        return # success

    def _Add_ImportMenu_powder(self,parent):
        '''configure the Powder Data menus accord to the readers found in _init_Imports
        '''
        submenu = wx.Menu()
        item = parent.AppendMenu(wx.ID_ANY, 'Powder Data',
            submenu, help='Import Powder data')
        for reader in self.ImportPowderReaderlist:
            item = submenu.Append(wx.ID_ANY,help=reader.longFormatName,
                kind=wx.ITEM_NORMAL,text='from '+reader.formatName+' file')
            self.ImportMenuId[item.GetId()] = reader
            self.Bind(wx.EVT_MENU, self.OnImportPowder, id=item.GetId())
        item = submenu.Append(wx.ID_ANY,
            help='Import powder data, use file to try to determine format',
            kind=wx.ITEM_NORMAL,text='guess format from file')
        self.Bind(wx.EVT_MENU, self.OnImportPowder, id=item.GetId())
        submenu.AppendSeparator()
        item = submenu.Append(wx.ID_ANY,
            help='Create a powder data set entry that will be simulated',
            kind=wx.ITEM_NORMAL,text='Simulate a dataset')
        self.Bind(wx.EVT_MENU, self.OnDummyPowder, id=item.GetId())
        
    def OpenPowderInstprm(self,instfile):
        '''Read a GSAS-II (new) instrument parameter file

        :param str instfile: name of instrument parameter file

        '''
        if os.path.splitext(instfile)[1].lower() != '.instprm': # invalid file
            return None            
        if not os.path.exists(instfile): # no such file
            return None
        File = open(instfile,'r')
        lines = File.readlines()
        File.close()
        return lines        
           
    def ReadPowderInstprm(self,instLines):
        '''Read lines from a GSAS-II (new) instrument parameter file

        :param list instLines: strings from GSAS-II parameter file; can be concatenated with ';'

        '''
        if not instLines[0].startswith('#GSAS-II'): # not a valid file
            return None
        newItems = []
        newVals = []
        for S in instLines:
            if S[0] == '#':
                continue
            S = S.replace(' ','')
            SS = S[:-1].split(';')
            for s in SS:
                [item,val] = s.split(':')
                newItems.append(item)
                try:
                    newVals.append(float(val))
                except ValueError:
                    newVals.append(val)                        
        return G2IO.makeInstDict(newItems,newVals,len(newVals)*[False,]),{}
        
    def ReadPowderIparm(self,instfile,bank,databanks,rd):
        '''Read a GSAS (old) instrument parameter file

        :param str instfile: name of instrument parameter file

        :param int bank: the bank number read in the raw data file

        :param int databanks: the number of banks in the raw data file.
          If the number of banks in the data and instrument parameter files
          agree, then the sets of banks are assumed to match up and bank
          is used to select the instrument parameter file. If not, the user
          is asked to make a selection.

        :param obj rd: the raw data (histogram) data object. This
          sets rd.instbank.

        '''
        if not os.path.exists(instfile): # no such file
            return {}
        fp = 0
        try:
            fp = open(instfile,'Ur')
            Iparm = {}
            for S in fp:
                if '#' in S[0]:
                    continue
                Iparm[S[:12]] = S[12:-1]
        except IOError:
            print('Error reading file:'+str(instfile))
        if fp:        
            fp.close()

        ibanks = int(Iparm.get('INS   BANK  ','1').strip())
        hType = Iparm['INS   HTYPE '].strip()
        if ibanks == 1: # there is only one bank here, return it
            rd.instbank = 1
            return Iparm
        if 'PNT' in hType:
            rd.instbank = bank
        elif ibanks != databanks:
            # number of banks in data and prm file not not agree, need a
            # choice from a human here
            choices = []
            for i in range(1,1+ibanks):
                choices.append('Bank '+str(i))
            bank = rd.BlockSelector(
                choices, self,
                title='Select an instrument parameter bank for '+
                os.path.split(rd.powderentry[0])[1]+' BANK '+str(bank)+
                '\nOr use Cancel to select from the default parameter sets',
                header='Block Selector')
        if bank is None: return {}
        # pull out requested bank # bank from the data, and change the bank to 1
        IparmS = {}
        for key in Iparm:
            if 'INS' in key[:3]:    #skip around rubbish lines in some old iparm files
                if key[4:6] == "  ":
                    IparmS[key] = Iparm[key]
                elif int(key[4:6].strip()) == bank:
                    IparmS[key[:4]+' 1'+key[6:]] = Iparm[key]
        rd.instbank = bank
        return IparmS
                        
    def GetPowderIparm(self,rd, prevIparm, lastIparmfile, lastdatafile):
        '''Open and read an instrument parameter file for a data file
        Returns the list of parameters used in the data tree

        :param obj rd: the raw data (histogram) data object.

        :param str prevIparm: not used

        :param str lastIparmfile: Name of last instrument parameter
          file that was read, or a empty string. 

        :param str lastdatafile: Name of last data file that was read.

        :returns: a list of two dicts, the first containing instrument parameters
          and the second used for TOf lookup tables for profile coeff.

        '''
        def SetPowderInstParms(Iparm, rd):
            '''extracts values from instrument parameters in rd.instdict
            or in array Iparm.
            Create and return the contents of the instrument parameter tree entry.
            '''
            Irads = {0:' ',1:'CrKa',2:'FeKa',3:'CuKa',4:'MoKa',5:'AgKa',6:'TiKa',7:'CoKa'}
            DataType = Iparm['INS   HTYPE '].strip()[:3]  # take 1st 3 chars
            # override inst values with values read from data file
            Bank = rd.powderentry[2]    #should be used in multibank iparm files
            if rd.instdict.get('type'):
                DataType = rd.instdict.get('type')
            data = [DataType,]
            instname = Iparm.get('INS  1INAME ')
            irad = int(Iparm.get('INS  1 IRAD ','0'))
            if instname:
                rd.Sample['InstrName'] = instname.strip()
            if 'C' in DataType:
                wave1 = None
                wave2 = 0.0
                if rd.instdict.get('wave'):
                    wl = rd.instdict.get('wave')
                    wave1 = wl[0]
                    if len(wl) > 1: wave2 = wl[1]
                s = Iparm['INS  1 ICONS']
                if not wave1:
                    wave1 = G2IO.sfloat(s[:10])
                    wave2 = G2IO.sfloat(s[10:20])
                v = (wave1,wave2,
                     G2IO.sfloat(s[20:30]),G2IO.sfloat(s[55:65]),G2IO.sfloat(s[40:50])) #get lam1, lam2, zero, pola & ratio
                if not v[1]:
                    names = ['Type','Lam','Zero','Polariz.','U','V','W','X','Y','SH/L','Azimuth'] 
                    v = (v[0],v[2],v[4])
                    codes = [0,0,0,0]
                else:
                    names = ['Type','Lam1','Lam2','Zero','I(L2)/I(L1)','Polariz.','U','V','W','X','Y','SH/L','Azimuth']
                    codes = [0,0,0,0,0,0]
                data.extend(v)
                if 'INS  1PRCF  ' in Iparm:
                    v1 = Iparm['INS  1PRCF  '].split()                                                  
                    v = Iparm['INS  1PRCF 1'].split()
                    data.extend([float(v[0]),float(v[1]),float(v[2])])                  #get GU, GV & GW - always here
                    azm = float(Iparm.get('INS  1DETAZM','0.0'))
                    v = Iparm['INS  1PRCF 2'].split()
                    if v1[0] == 3:
                        data.extend([float(v[0]),float(v[1]),float(v[2])+float(v[3],azm)])  #get LX, LY, S+H/L & azimuth
                    else:
                        data.extend([0.0,0.0,0.002,azm])                                      #OK defaults if fxn #3 not 1st in iprm file                    
                else:
                    v1 = Iparm['INS  1PRCF1 '].split()                                                  
                    v = Iparm['INS  1PRCF11'].split()
                    data.extend([float(v[0]),float(v[1]),float(v[2])])                  #get GU, GV & GW - always here
                    azm = float(Iparm.get('INS  1DETAZM','0.0'))
                    v = Iparm['INS  1PRCF12'].split()
                    if v1[0] == 3:
                        data.extend([float(v[0]),float(v[1]),float(v[2])+float(v[3],azm)])  #get LX, LY, S+H/L & azimuth
                    else:
                        data.extend([0.0,0.0,0.002,azm])                                      #OK defaults if fxn #3 not 1st in iprm file
                codes.extend([0,0,0,0,0,0,0])
                Iparm1 = G2IO.makeInstDict(names,data,codes)
                Iparm1['Source'] = [Irads[irad],Irads[irad]]
                Iparm1['Bank'] = [Bank,Bank,0]
                return [Iparm1,{}]
            elif 'T' in DataType:
                names = ['Type','fltPath','2-theta','difC','difA', 'difB','Zero','alpha','beta-0','beta-1',
                    'beta-q','sig-0','sig-1','sig-2','sig-q', 'X','Y','Azimuth',]
                codes = [0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0, 0,0,0,]
                azm = 0.
                if 'INS  1DETAZM' in Iparm:
                    azm = float(Iparm['INS  1DETAZM'])
                s = Iparm['INS   FPATH1'].split()
                fltPath0 = G2IO.sfloat(s[0])
                s = Iparm['INS  1BNKPAR'].split()
                fltPath1 = G2IO.sfloat(s[0])
                data.extend([fltPath0+fltPath1,])               #Flight path source-sample-detector
                data.extend([G2IO.sfloat(s[1]),])               #2-theta for bank
                s = Iparm['INS  1 ICONS'].split()
                data.extend([G2IO.sfloat(s[0]),G2IO.sfloat(s[1]),0.0,G2IO.sfloat(s[2])])    #difC,difA,difB,Zero
                if 'INS  1PRCF  ' in Iparm:
                    s = Iparm['INS  1PRCF  '].split()
                    pfType = int(s[0])
                    s = Iparm['INS  1PRCF 1'].split()
                    if abs(pfType) == 1:
                        data.extend([G2IO.sfloat(s[1]),G2IO.sfloat(s[2]),G2IO.sfloat(s[3])]) #alpha, beta-0, beta-1
                        s = Iparm['INS  1PRCF 2'].split()
                        data.extend([0.0,0.0,G2IO.sfloat(s[1]),G2IO.sfloat(s[2]),0.0,0.0,0.0,azm])    #beta-q, sig-0, sig-1, sig-2, sig-q, X, Y
                    elif abs(pfType) in [3,4,5]:
                        data.extend([G2IO.sfloat(s[0]),G2IO.sfloat(s[1]),G2IO.sfloat(s[2])]) #alpha, beta-0, beta-1
                        if abs(pfType) == 4:
                            data.extend([0.0,0.0,G2IO.sfloat(s[3]),0.0,0.0,0.0,0.0,azm])    #beta-q, sig-0, sig-1, sig-2, sig-q, X, Y
                        else:
                            s = Iparm['INS  1PRCF 2'].split()
                            data.extend([0.0,0.0,G2IO.sfloat(s[0]),G2IO.sfloat(s[1]),0.0,0.0,0.0,azm])    #beta-q, sig-0, sig-1, sig-2, sig-q, X, Y                       
                else:
                    s = Iparm['INS  1PRCF1 '].split()
                    pfType = int(s[0])
                    s = Iparm['INS  1PRCF11'].split()
                    if abs(pfType) == 1:
                        data.extend([G2IO.sfloat(s[1]),G2IO.sfloat(s[2]),G2IO.sfloat(s[3])]) #alpha, beta-0, beta-1
                        s = Iparm['INS  1PRCF12'].split()
                        data.extend([0.0,0.0,G2IO.sfloat(s[1]),G2IO.sfloat(s[2]),0.0,0.0,0.0,azm])    #beta-q, sig-0, sig-1, sig-2, sig-q, X, Y
                    elif abs(pfType) in [3,4,5]:
                        data.extend([G2IO.sfloat(s[0]),G2IO.sfloat(s[1]),G2IO.sfloat(s[2])]) #alpha, beta-0, beta-1
                        if abs(pfType) == 4:
                            data.extend([0.0,0.0,G2IO.sfloat(s[3]),0.0,0.0,0.0,0.0,azm])    #beta-q, sig-0, sig-1, sig-2, sig-q, X, Y
                        else:
                            s = Iparm['INS  1PRCF12'].split()
                            data.extend([0.0,0.0,G2IO.sfloat(s[0]),G2IO.sfloat(s[1]),0.0,0.0,0.0,azm])    #beta-q, sig-0, sig-1, sig-2, sig-q, X, Y                       
                Inst1 = G2IO.makeInstDict(names,data,codes)
                Inst1['Bank'] = [Bank,Bank,0]
                Inst2 = {}
                if pfType < 0:
                    Ipab = 'INS  1PAB'+str(-pfType)
                    Npab = int(Iparm[Ipab+'  '].strip())
                    Inst2['Pdabc'] = []
                    for i in range(Npab):
                        k = Ipab+str(i+1).rjust(2)
                        s = Iparm[k].split()
                        Inst2['Pdabc'].append([float(t) for t in s])
                    Inst2['Pdabc'] = np.array(Inst2['Pdabc'])
                    Inst2['Pdabc'].T[3] += Inst2['Pdabc'].T[0]*Inst1['difC'][0] #turn 3rd col into TOF
                if 'INS  1I ITYP' in Iparm:
                    s = Iparm['INS  1I ITYP'].split()
                    Ityp = int(s[0])
                    Tminmax = [float(s[1])*1000.,float(s[2])*1000.]
                    Itypes = ['Exponential','Maxwell/Exponential','','Maxwell/Chebyschev','']
                    if Ityp in [1,2,4]:
                        Inst2['Itype'] = Itypes[Ityp-1]
                        Inst2['Tminmax'] = Tminmax
                        Icoeff = []
                        Iesd = []
                        Icovar = []                   
                        for i in range(3):
                            s = Iparm['INS  1ICOFF'+str(i+1)].split()
                            Icoeff += [float(S) for S in s]
                            s = Iparm['INS  1IECOF'+str(i+1)].split()
                            Iesd += [float(S) for S in s]
                        NT = 10
                        for i in range(8):
                            s = Iparm['INS  1IECOR'+str(i+1)]
                            if i == 7:
                                NT = 8
                            Icovar += [float(s[6*j:6*j+6]) for j in range(NT)]
                        Inst2['Icoeff'] = Icoeff
                        Inst2['Iesd'] = Iesd
                        Inst2['Icovar'] = Icovar
                return [Inst1,Inst2]

        # stuff we might need from the reader
        filename = rd.powderentry[0]
        bank = rd.powderentry[2]
        numbanks = rd.numbanks
        # is there an instrument parameter file defined for the current data set?
        # or if this is a read on a set of set of files, use the last one again
        #if rd.instparm or (lastdatafile == filename and lastIparmfile):
        if rd.instparm or lastIparmfile:
            if rd.instparm:
                instfile = os.path.join(os.path.split(filename)[0],
                                    rd.instparm)
            else:
                # for multiple reads of one data file, reuse the inst parm file
                instfile = lastIparmfile
            if os.path.exists(instfile):
                #print 'debug: try read',instfile
                Lines = self.OpenPowderInstprm(instfile)
                instParmList = None
                if Lines is not None:
                    instParmList = self.ReadPowderInstprm(Lines)
                if instParmList is not None:
                    rd.instfile = instfile
                    rd.instmsg = 'GSAS-II file '+instfile
                    return instParmList
                Iparm = self.ReadPowderIparm(instfile,bank,numbanks,rd)
                if Iparm:
                    #print 'debug: success'
                    rd.instfile = instfile
                    rd.instmsg = instfile + ' bank ' + str(rd.instbank)
                    return SetPowderInstParms(Iparm,rd)
            else:
                self.ErrorDialog('Open Error','Error opening instrument parameter file '
                    +str(instfile)+' requested by file '+ filename)
        # is there an instrument parameter file matching the current file
        # with extension .inst or .prm? If so read it
        basename = os.path.splitext(filename)[0]
        for ext in '.instprm','.prm','.inst','.ins':
            instfile = basename + ext
            Lines = self.OpenPowderInstprm(instfile)
            instParmList = None
            if Lines is not None:
                instParmList = self.ReadPowderInstprm(Lines)
            if instParmList is not None:
                rd.instfile = instfile
                rd.instmsg = 'GSAS-II file '+instfile
                return instParmList
            Iparm = self.ReadPowderIparm(instfile,bank,numbanks,rd)
            if Iparm:
                #print 'debug: success'
                rd.instfile = instfile
                rd.instmsg = instfile + ' bank ' + str(rd.instbank)
                return SetPowderInstParms(Iparm,rd)
            else:
                #print 'debug: open/read failed',instfile
                pass # fail silently

        # did we read the data file from a zip? If so, look there for a
        # instrument parameter file
        if self.zipfile:
            for ext in '.instprm','.prm','.inst','.ins':
                instfile = G2IO.ExtractFileFromZip(
                    self.zipfile,
                    selection=os.path.split(basename + ext)[1],
                    parent=self)
                if instfile is not None and instfile != self.zipfile:
                    print 'debug:',instfile,'created from ',self.zipfile
                    Lines = self.OpenPowderInstprm(instfile)
                    instParmList = None
                    if Lines is not None:
                        instParmList = self.ReadPowderInstprm(Lines)
                    if instParmList is not None:
                        rd.instfile = instfile
                        rd.instmsg = 'GSAS-II file '+instfile
                        return instParmList
                    Iparm = self.ReadPowderIparm(instfile,bank,numbanks,rd)
                    if Iparm:
                        rd.instfile = instfile
                        rd.instmsg = instfile + ' bank ' + str(rd.instbank)
                        return SetPowderInstParms(Iparm,rd)
                    else:
                        #print 'debug: open/read for',instfile,'from',self.zipfile,'failed'
                        pass # fail silently

        while True: # loop until we get a file that works or we get a cancel
            instfile = ''
            dlg = wx.FileDialog(
                self,
                'Choose inst. param file for "'
                +rd.idstring
                +'" (or Cancel for default)',
                '.', '',
                'GSAS iparm file (*.prm,*.inst,*.ins)|*.prm;*.inst;*.ins|'
                'GSAS-II iparm file (*.instprm)|*.instprm|'
                'All files (*.*)|*.*', 
                wx.OPEN|wx.CHANGE_DIR)
            if os.path.exists(lastIparmfile):
                dlg.SetFilename(lastIparmfile)
            if dlg.ShowModal() == wx.ID_OK:
                instfile = dlg.GetPath()
            dlg.Destroy()
            if not instfile: break
            Lines = self.OpenPowderInstprm(instfile)
            instParmList = None
            if Lines is not None:
                instParmList = self.ReadPowderInstprm(Lines)
            if instParmList is not None:
                rd.instfile = instfile
                rd.instmsg = 'GSAS-II file '+instfile
                return instParmList
            Iparm = self.ReadPowderIparm(instfile,bank,numbanks,rd)
            if Iparm:
                #print 'debug: success with',instfile
                rd.instfile = instfile
                rd.instmsg = instfile + ' bank ' + str(rd.instbank)
                return SetPowderInstParms(Iparm,rd)
            else:
                self.ErrorDialog('Read Error',
                                 'Error opening/reading file '+str(instfile))
        
        # still no success: offer user choice of defaults
        import defaultIparms as dI
        while True: # loop until we get a choice
            choices = []
            head = 'Select from default instrument parameters for '+rd.idstring

            for l in dI.defaultIparm_lbl:
                choices.append('Defaults for '+l)
            res = rd.BlockSelector(
                choices,
                ParentFrame=self,
                title=head,
                header='Select default inst parms',
                useCancel=False)
            if res is None: continue
            rd.instfile = ''
            rd.instmsg = 'default: '+dI.defaultIparm_lbl[res]
            return self.ReadPowderInstprm(dI.defaultIparms[res])

    def OnImportPowder(self,event):
        '''Called in response to an Import/Powder Data/... menu item
        to read a powder diffraction data set. 
        dict self.ImportMenuId is used to look up the specific
        reader item associated with the menu item, which will be
        None for the last menu item, which is the "guess" option
        where all appropriate formats will be tried.

        Also reads an instrument parameter file for each dataset.
        '''
        # get a list of existing histograms
        PWDRlist = []
        if self.PatternTree.GetCount():
            item, cookie = self.PatternTree.GetFirstChild(self.root)
            while item:
                name = self.PatternTree.GetItemText(item)
                if name.startswith('PWDR ') and name not in PWDRlist:
                    PWDRlist.append(name)
                item, cookie = self.PatternTree.GetNextChild(self.root, cookie)
        # look up which format was requested
        reqrdr = self.ImportMenuId.get(event.GetId())  
        rdlist = self.OnImportGeneric(
            reqrdr,self.ImportPowderReaderlist,'Powder Data',multiple=True)
        if len(rdlist) == 0: return
        self.CheckNotebook()
        Iparm = None
        lastIparmfile = ''
        lastdatafile = ''
        newHistList = []
        self.EnablePlot = False
        for rd in rdlist:
            if 'Instrument Parameters' not in rd.pwdparms:
                # get instrument parameters for each dataset, unless already set 
                Iparm1,Iparm2 = self.GetPowderIparm(rd, Iparm, lastIparmfile, lastdatafile)
                if rd.repeat_instparm: 
                    lastIparmfile = rd.instfile
                # override any keys in read instrument parameters with ones set in import
                for key in Iparm1:  
                    if key in rd.instdict:
                        Iparm1[key] = rd.instdict[key]
            else:
                Iparm1,Iparm2 = rd.pwdparms['Instrument Parameters']
            lastdatafile = rd.powderentry[0]
            HistName = rd.idstring
            HistName = 'PWDR '+HistName
            # make new histogram names unique
            HistName = G2obj.MakeUniqueLabel(HistName,PWDRlist)
            print 'Read powder data '+str(HistName)+ \
                ' from file '+str(rd.readfilename) + \
                ' with parameters from '+str(rd.instmsg)
            # data are read, now store them in the tree
            Id = self.PatternTree.AppendItem(parent=self.root,text=HistName)
            if 'T' in Iparm1['Type'][0]:
                if not rd.clockWd and rd.GSAS:
                    rd.powderdata[0] *= 100.        #put back the CW centideg correction
                cw = np.diff(rd.powderdata[0])
                rd.powderdata[0] = rd.powderdata[0][:-1]+cw/2.
                if rd.GSAS:     #NB: old GSAS wanted intensities*CW even if normalized!
                    npts = min(len(rd.powderdata[0]),len(rd.powderdata[1]),len(cw))
                    rd.powderdata[1] = rd.powderdata[1][:npts]/cw[:npts]
                    rd.powderdata[2] = rd.powderdata[2][:npts]*cw[:npts]**2  #1/var=w at this point
                else:       #NB: from topas/fullprof type files
                    rd.powderdata[1] = rd.powderdata[1][:-1]
                    rd.powderdata[2] = rd.powderdata[2][:-1]
                if 'Itype' in Iparm2:
                    Ibeg = np.searchsorted(rd.powderdata[0],Iparm2['Tminmax'][0])
                    Ifin = np.searchsorted(rd.powderdata[0],Iparm2['Tminmax'][1])
                    rd.powderdata[0] = rd.powderdata[0][Ibeg:Ifin]
                    YI,WYI = G2pwd.calcIncident(Iparm2,rd.powderdata[0])
                    rd.powderdata[1] = rd.powderdata[1][Ibeg:Ifin]/YI
                    var = 1./rd.powderdata[2][Ibeg:Ifin]
                    var += WYI*rd.powderdata[1]**2
                    var /= YI**2
                    rd.powderdata[2] = 1./var
                rd.powderdata[3] = np.zeros_like(rd.powderdata[0])                                        
                rd.powderdata[4] = np.zeros_like(rd.powderdata[0])                                        
                rd.powderdata[5] = np.zeros_like(rd.powderdata[0])                                        
            valuesdict = {
                'wtFactor':1.0,
                'Dummy':False,
                'ranId':ran.randint(0,sys.maxint),
                'Offset':[0.0,0.0],'delOffset':0.02,'refOffset':-1.0,'refDelt':0.01,
                'qPlot':False,'dPlot':False,'sqrtPlot':False
                }
            rd.Sample['ranId'] = valuesdict['ranId'] # this should be removed someday
            self.PatternTree.SetItemPyData(Id,[valuesdict,rd.powderdata])
            self.PatternTree.SetItemPyData(
                self.PatternTree.AppendItem(Id,text='Comments'),
                rd.comments)
            Tmin = min(rd.powderdata[0])
            Tmax = max(rd.powderdata[0])
            self.PatternTree.SetItemPyData(
                self.PatternTree.AppendItem(Id,text='Limits'),
                rd.pwdparms.get('Limits',[(Tmin,Tmax),[Tmin,Tmax]])
                )
            self.PatternId = G2gd.GetPatternTreeItemId(self,Id,'Limits')
            self.PatternTree.SetItemPyData(
                self.PatternTree.AppendItem(Id,text='Background'),
                rd.pwdparms.get('Background',
                    [['chebyschev',True,3,1.0,0.0,0.0],{'nDebye':0,'debyeTerms':[],'nPeaks':0,'peaksList':[]}])
                    )
            self.PatternTree.SetItemPyData(
                self.PatternTree.AppendItem(Id,text='Instrument Parameters'),
                [Iparm1,Iparm2])
            self.PatternTree.SetItemPyData(
                self.PatternTree.AppendItem(Id,text='Sample Parameters'),
                rd.Sample)
            self.PatternTree.SetItemPyData(
                self.PatternTree.AppendItem(Id,text='Peak List')
                ,{'peaks':[],'sigDict':{}})
            self.PatternTree.SetItemPyData(
                self.PatternTree.AppendItem(Id,text='Index Peak List'),
                [[],[]])
            self.PatternTree.SetItemPyData(
                self.PatternTree.AppendItem(Id,text='Unit Cells List'),
                [])
            self.PatternTree.SetItemPyData(
                self.PatternTree.AppendItem(Id,text='Reflection Lists'),
                {})
            # if any Control values have been set, move them into tree
            Controls = self.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(self,self.root, 'Controls'))
            Controls.update(rd.Controls)
            newHistList.append(HistName)
        else:
            self.EnablePlot = True
            self.PatternTree.Expand(Id)
            self.PatternTree.SelectItem(Id)
            
        if not newHistList: return # somehow, no new histograms
        # make a list of phase names
        phaseRIdList,usedHistograms = self.GetPhaseInfofromTree()
        phaseNameList = usedHistograms.keys() # phase names in use
        if not phaseNameList: return # no phases yet, nothing to do
        header = 'Select phase(s) to add the new\npowder dataset(s) to:'
        for Name in newHistList:
            header += '\n  '+str(Name)

        result = G2G.ItemSelector(phaseNameList,self,header,header='Add to phase(s)',multiple=True)
        if not result: return
        # connect new phases to histograms
        sub = G2gd.GetPatternTreeItemId(self,self.root,'Phases')
        if not sub:
            raise Exception('ERROR -- why are there no phases here?')
        item, cookie = self.PatternTree.GetFirstChild(sub)
        iph = -1
        while item: # loop over (new) phases
            iph += 1
            phaseName = self.PatternTree.GetItemText(item)
            data = self.PatternTree.GetItemPyData(item)
            item, cookie = self.PatternTree.GetNextChild(sub, cookie)
            if iph not in result: continue
            generalData = data['General']
            SGData = generalData['SGData']
            UseList = data['Histograms']
            NShkl = len(G2spc.MustrainNames(SGData))
            NDij = len(G2spc.HStrainNames(SGData))
            for histoName in newHistList:
                UseList[histoName] = SetDefaultDData('PWDR',histoName,NShkl=NShkl,NDij=NDij)
                Id = G2gd.GetPatternTreeItemId(self,self.root,histoName)
                refList = self.PatternTree.GetItemPyData(
                    G2gd.GetPatternTreeItemId(self,Id,'Reflection Lists'))
                refList[generalData['Name']] = []
        return # success

    def OnDummyPowder(self,event):
        '''Called in response to Import/Powder Data/Simulate menu item
        to create a Dummy powder diffraction data set. 

        Reads an instrument parameter file and then gets input from the user
        '''
        # get a list of existing histograms
        PWDRlist = []
        if self.PatternTree.GetCount():
            item, cookie = self.PatternTree.GetFirstChild(self.root)
            while item:
                name = self.PatternTree.GetItemText(item)
                if name.startswith('PWDR ') and name not in PWDRlist:
                    PWDRlist.append(name)
                item, cookie = self.PatternTree.GetNextChild(self.root, cookie)
        # Initialize a base class reader
        rd = G2IO.ImportPowderData(
            extensionlist=tuple(),
            strictExtension=False,
            formatName = 'Simulate dataset',
            longFormatName = 'Compute a simulated pattern')
        rd.powderentry[0] = '' # no filename
        # #self.powderentry[1] = pos # bank offset (N/A here)
        rd.powderentry[2] = 1 # only one bank
        rd.comments.append('This is a dummy dataset for powder pattern simulation')
        self.CheckNotebook()
        Iparm = None
        lastIparmfile = ''
        lastdatafile = ''
        self.zipfile = None
        # get instrument parameters for it
        Iparm1,Iparm2 = self.GetPowderIparm(rd, Iparm, lastIparmfile, lastdatafile)
        if 'T' in Iparm1['Type'][0]:
            print('TOF simulation not supported yet')
            return False
        else:
            # need to get name, 2theta start, end, step
            rd.idstring = ' CW'
            if 'X' in Iparm1['Type'][0]:
                rd.idstring = 'CW x-ray simulation'
            else:
                rd.idstring = 'CW neutron simulation'
            # base initial range on wavelength
            wave = Iparm1.get('Lam')
            if wave:
                wave = wave[0]
            else:
                wave = Iparm1.get('Lam1')
                if wave:
                    wave = wave[0]
        N = 0
        while (N < 3): # insist on a dataset with a few points
            names = ('dataset name', 'start angle', 'end angle', 'step size')
            if not wave or wave < 1.0:
                inp = [rd.idstring, 10.,40.,0.005] # see names for what's what
            else:
                inp = [rd.idstring, 10.,80.,0.01] # see names for what's what
            dlg = G2G.ScrolledMultiEditor(
                self,[inp] * len(inp),range(len(inp)),names,
                header='Enter simulation name and range',
                minvals=(None,0.001,0.001,0.0001),
                maxvals=(None,180.,180.,.1),
                sizevals=((225,-1),)
                )
            dlg.CenterOnParent()
            if dlg.ShowModal() == wx.ID_OK:
                if inp[1] > inp[2]:
                    end,start,step = inp[1:]
                else:                
                    start,end,step = inp[1:]
                step = abs(step)
            else:
                return False
            N = int((end-start)/step)+1
            x = np.linspace(start,end,N,True)
            N = len(x)
        rd.powderdata = [
            np.array(x), # x-axis values
            np.zeros_like(x), # powder pattern intensities
            np.ones_like(x), # 1/sig(intensity)^2 values (weights)
            np.zeros_like(x), # calc. intensities (zero)
            np.zeros_like(x), # calc. background (zero)
            np.zeros_like(x), # obs-calc profiles
            ]
        Tmin = rd.powderdata[0][0]
        Tmax = rd.powderdata[0][-1]
        # data are read, now store them in the tree
        HistName = inp[0]
        HistName = 'PWDR '+HistName
        HistName = G2obj.MakeUniqueLabel(HistName,PWDRlist)  # make new histogram names unique
        Id = self.PatternTree.AppendItem(parent=self.root,text=HistName)
        valuesdict = {
            'wtFactor':1.0,
            'Dummy':True,
            'ranId':ran.randint(0,sys.maxint),
            'Offset':[0.0,0.0],'delOffset':0.02,'refOffset':-1.0,'refDelt':0.01,
            'qPlot':False,'dPlot':False,'sqrtPlot':False
            }
        self.PatternTree.SetItemPyData(Id,[valuesdict,rd.powderdata])
        self.PatternTree.SetItemPyData(
            self.PatternTree.AppendItem(Id,text='Comments'),
            rd.comments)
        self.PatternTree.SetItemPyData(
            self.PatternTree.AppendItem(Id,text='Limits'),
            [(Tmin,Tmax),[Tmin,Tmax]])
        self.PatternId = G2gd.GetPatternTreeItemId(self,Id,'Limits')
        self.PatternTree.SetItemPyData(
            self.PatternTree.AppendItem(Id,text='Background'),
            [['chebyschev',True,3,1.0,0.0,0.0],
             {'nDebye':0,'debyeTerms':[],'nPeaks':0,'peaksList':[]}])
        self.PatternTree.SetItemPyData(
            self.PatternTree.AppendItem(Id,text='Instrument Parameters'),
            [Iparm1,Iparm2])
        self.PatternTree.SetItemPyData(
            self.PatternTree.AppendItem(Id,text='Sample Parameters'),
            rd.Sample)
        self.PatternTree.SetItemPyData(
            self.PatternTree.AppendItem(Id,text='Peak List')
            ,{'peaks':[],'sigDict':{}})
        self.PatternTree.SetItemPyData(
            self.PatternTree.AppendItem(Id,text='Index Peak List'),
            [[],[]])
        self.PatternTree.SetItemPyData(
            self.PatternTree.AppendItem(Id,text='Unit Cells List'),
            [])
        self.PatternTree.SetItemPyData(
            self.PatternTree.AppendItem(Id,text='Reflection Lists'),
            {})
        self.PatternTree.Expand(Id)
        self.PatternTree.SelectItem(Id)
        print('Added simulation powder data '+str(HistName)+ 
              ' with parameters from '+str(rd.instmsg))

        # make a list of phase names
        phaseRIdList,usedHistograms = self.GetPhaseInfofromTree()
        phaseNameList = usedHistograms.keys() # phase names in use
        if not phaseNameList: return # no phases yet, nothing to do
        header = 'Select phase(s) to add the new\npowder simulation (dummy) dataset to:'
        result = G2G.ItemSelector(phaseNameList,self,header,header='Add to phase(s)',multiple=True)
        if not result: return
        # connect new phases to histograms
        sub = G2gd.GetPatternTreeItemId(self,self.root,'Phases')
        if not sub:
            raise Exception('ERROR -- why are there no phases here?')
        item, cookie = self.PatternTree.GetFirstChild(sub)
        iph = -1
        while item: # loop over (new) phases
            iph += 1
            phaseName = self.PatternTree.GetItemText(item)
            data = self.PatternTree.GetItemPyData(item)
            item, cookie = self.PatternTree.GetNextChild(sub, cookie)
            if iph not in result: continue
            generalData = data['General']
            SGData = generalData['SGData']
            UseList = data['Histograms']
            NShkl = len(G2spc.MustrainNames(SGData))
            NDij = len(G2spc.HStrainNames(SGData))
            UseList[HistName] = SetDefaultDData('PWDR',HistName,NShkl=NShkl,NDij=NDij)
            Id = G2gd.GetPatternTreeItemId(self,self.root,HistName)
            refList = self.PatternTree.GetItemPyData(
                G2gd.GetPatternTreeItemId(self,Id,'Reflection Lists'))
            refList[generalData['Name']] = []
        return # success
        
    def OnPreferences(self,event):
        'Edit the GSAS-II configuration variables'
        dlg = G2G.SelectConfigSetting(self)
        dlg.ShowModal() == wx.ID_OK
        dlg.Destroy()

    def _Add_ImportMenu_smallangle(self,parent):
        '''configure the Small Angle Data menus accord to the readers found in _init_Imports
        '''
        submenu = wx.Menu()
        item = parent.AppendMenu(wx.ID_ANY, 'Small Angle Data',
            submenu, help='Import small angle data')
        for reader in self.ImportSmallAngleReaderlist:
            item = submenu.Append(wx.ID_ANY,help=reader.longFormatName,
                kind=wx.ITEM_NORMAL,text='from '+reader.formatName+' file')
            self.ImportMenuId[item.GetId()] = reader
            self.Bind(wx.EVT_MENU, self.OnImportSmallAngle, id=item.GetId())
        # item = submenu.Append(wx.ID_ANY,
        #     help='Import small angle data, use file to try to determine format',
        #     kind=wx.ITEM_NORMAL,text='guess format from file')
        # self.Bind(wx.EVT_MENU, self.OnImportSmallAngle, id=item.GetId())

    def OnImportSmallAngle(self,event):
        '''Called in response to an Import/Small Angle Data/... menu item
        to read a small angle diffraction data set. 
        dict self.ImportMenuId is used to look up the specific
        reader item associated with the menu item, which will be
        None for the last menu item, which is the "guess" option
        where all appropriate formats will be tried.

        '''
        
        def GetSASDIparm(reader):
            parm = reader.instdict
            Iparm = {'Type':[parm['type'],parm['type'],0],'Lam':[parm['wave'],
                parm['wave'],0],'Azimuth':[0.,0.,0]}           
            return Iparm,{}
            
        # get a list of existing histograms
        SASDlist = []
        if self.PatternTree.GetCount():
            item, cookie = self.PatternTree.GetFirstChild(self.root)
            while item:
                name = self.PatternTree.GetItemText(item)
                if name.startswith('SASD ') and name not in SASDlist:
                    SASDlist.append(name)
                item, cookie = self.PatternTree.GetNextChild(self.root, cookie)
        # look up which format was requested
        reqrdr = self.ImportMenuId.get(event.GetId())  
        rdlist = self.OnImportGeneric(
            reqrdr,self.ImportSmallAngleReaderlist,'Small Angle Data',multiple=True)
        if len(rdlist) == 0: return
        self.CheckNotebook()
        Iparm = None
        lastdatafile = ''
        newHistList = []
        self.EnablePlot = False
        for rd in rdlist:
            lastdatafile = rd.smallangleentry[0]
            HistName = rd.idstring
            HistName = 'SASD '+HistName
            # make new histogram names unique
            HistName = G2obj.MakeUniqueLabel(HistName,SASDlist)
            print 'Read small angle data '+str(HistName)+ \
                ' from file '+str(self.lastimport)
            # data are read, now store them in the tree
            Id = self.PatternTree.AppendItem(parent=self.root,text=HistName)
            Iparm1,Iparm2 = GetSASDIparm(rd)
#            if 'T' in Iparm1['Type'][0]:
#                if not rd.clockWd and rd.GSAS:
#                    rd.powderdata[0] *= 100.        #put back the CW centideg correction
#                cw = np.diff(rd.powderdata[0])
#                rd.powderdata[0] = rd.powderdata[0][:-1]+cw/2.
#                rd.powderdata[1] = rd.powderdata[1][:-1]/cw
#                rd.powderdata[2] = rd.powderdata[2][:-1]*cw**2  #1/var=w at this point
#                if 'Itype' in Iparm2:
#                    Ibeg = np.searchsorted(rd.powderdata[0],Iparm2['Tminmax'][0])
#                    Ifin = np.searchsorted(rd.powderdata[0],Iparm2['Tminmax'][1])
#                    rd.powderdata[0] = rd.powderdata[0][Ibeg:Ifin]
#                    YI,WYI = G2pwd.calcIncident(Iparm2,rd.powderdata[0])
#                    rd.powderdata[1] = rd.powderdata[1][Ibeg:Ifin]/YI
#                    var = 1./rd.powderdata[2][Ibeg:Ifin]
#                    var += WYI*rd.powderdata[1]**2
#                    var /= YI**2
#                    rd.powderdata[2] = 1./var
#                rd.powderdata[3] = np.zeros_like(rd.powderdata[0])                                        
#                rd.powderdata[4] = np.zeros_like(rd.powderdata[0])                                        
#                rd.powderdata[5] = np.zeros_like(rd.powderdata[0])                                        
            Tmin = min(rd.smallangledata[0])
            Tmax = max(rd.smallangledata[0])
            valuesdict = {
                'wtFactor':1.0,
                'Dummy':False,
                'ranId':ran.randint(0,sys.maxint),
                }
            rd.Sample['ranId'] = valuesdict['ranId'] # this should be removed someday
            self.PatternTree.SetItemPyData(Id,[valuesdict,rd.smallangledata])
            self.PatternTree.SetItemPyData(
                self.PatternTree.AppendItem(Id,text='Comments'),
                rd.comments)
            self.PatternTree.SetItemPyData(
                self.PatternTree.AppendItem(Id,text='Limits'),
                [(Tmin,Tmax),[Tmin,Tmax]])
            self.PatternId = G2gd.GetPatternTreeItemId(self,Id,'Limits')
            self.PatternTree.SetItemPyData(
                self.PatternTree.AppendItem(Id,text='Instrument Parameters'),
                [Iparm1,Iparm2])
            self.PatternTree.SetItemPyData(
                self.PatternTree.AppendItem(Id,text='Substances'),G2pdG.SetDefaultSubstances())
            self.PatternTree.SetItemPyData(
                self.PatternTree.AppendItem(Id,text='Sample Parameters'),
                rd.Sample)
            self.PatternTree.SetItemPyData(
                self.PatternTree.AppendItem(Id,text='Models'),G2pdG.SetDefaultSASDModel())
            newHistList.append(HistName)
        else:
            self.EnablePlot = True
            self.PatternTree.Expand(Id)
            self.PatternTree.SelectItem(Id)
            
        if not newHistList: return # somehow, no new histograms
        return # success

    def OnMacroRecordStatus(self,event,setvalue=None):
        '''Called when the record macro menu item is used which toggles the
        value. Alternately a value to be set can be provided. Note that this
        routine is made more complex because on the Mac there are lots of menu
        items (listed in self.MacroStatusList) and this loops over all of them.
        '''
        nextvalue = log.ShowLogStatus() != True
        if setvalue is not None:
            nextvalue = setvalue 
        if nextvalue:
            log.LogOn()
            set2 = True
        else:
            log.LogOff()
            set2 = False
        for menuitem in self.MacroStatusList:
            menuitem.Check(set2)

    def _init_Macro(self):
        '''Define the items in the macro menu.
        '''
        menu = self.MacroMenu
        item = menu.Append(
                help='Start or stop recording of menu actions, etc.', id=wx.ID_ANY,
                kind=wx.ITEM_CHECK,text='Record actions')
        self.MacroStatusList.append(item)
        item.Check(log.ShowLogStatus())
        self.Bind(wx.EVT_MENU, self.OnMacroRecordStatus, item)

        # this may only be of value for development work
        item = menu.Append(
            help='Show logged commands', id=wx.ID_ANY,
            kind=wx.ITEM_NORMAL,text='Show log')
        def OnShowLog(event):
            print 70*'='
            print 'List of logged actions'
            for i,line in enumerate(log.G2logList):
                if line: print i,line
            print 70*'='
        self.Bind(wx.EVT_MENU, OnShowLog, item)

        item = menu.Append(
            help='Clear logged commands', id=wx.ID_ANY,
            kind=wx.ITEM_NORMAL,text='Clear log')
        def OnClearLog(event): log.G2logList=[None]
        self.Bind(wx.EVT_MENU, OnClearLog, item)
        
        item = menu.Append(
            help='Save logged commands to file', id=wx.ID_ANY,
            kind=wx.ITEM_NORMAL,text='Save log')
        def OnSaveLog(event):
            import cPickle
            defnam = os.path.splitext(
                os.path.split(self.GSASprojectfile)[1]
                )[0]+'.gcmd'
            dlg = wx.FileDialog(self,
                'Choose an file to save past actions', '.', defnam, 
                'GSAS-II cmd output (*.gcmd)|*.gcmd',
                wx.FD_SAVE|wx.FD_OVERWRITE_PROMPT|wx.CHANGE_DIR)
            dlg.CenterOnParent()
            try:
                if dlg.ShowModal() == wx.ID_OK:
                    filename = dlg.GetPath()
                    # make sure extension is correct
                    filename = os.path.splitext(filename)[0]+'.gcmd'
                else:
                    filename = None
            finally:
                dlg.Destroy()
            if filename:
                fp = open(filename,'wb')
                fp.write(str(len(log.G2logList)-1)+'\n')
                for item in log.G2logList:
                    if item: cPickle.dump(item,fp)
                fp.close()
        self.Bind(wx.EVT_MENU, OnSaveLog, item)

        item = menu.Append(
            help='Load logged commands from file', id=wx.ID_ANY,
            kind=wx.ITEM_NORMAL,text='Load log')
        def OnLoadLog(event):
            # this appends. Perhaps we should ask to clear? 
            import cPickle
            defnam = os.path.splitext(
                os.path.split(self.GSASprojectfile)[1]
                )[0]+'.gcmd'
            dlg = wx.FileDialog(self,
                'Choose an file to read saved actions', '.', defnam, 
                'GSAS-II cmd output (*.gcmd)|*.gcmd',
                wx.OPEN|wx.CHANGE_DIR)
            dlg.CenterOnParent()
            try:
                if dlg.ShowModal() == wx.ID_OK:
                    filename = dlg.GetPath()
                    # make sure extension is correct
                    filename = os.path.splitext(filename)[0]+'.gcmd'
                else:
                    filename = None
            finally:
                dlg.Destroy()
            if filename and os.path.exists(filename):
                fp = open(filename,'rb')
                lines = fp.readline()
                for i in range(int(lines)):
                    log.G2logList.append(cPickle.load(fp))
                fp.close()
        self.Bind(wx.EVT_MENU, OnLoadLog, item)

        item = menu.Append(
            help='Replay saved commands', id=wx.ID_ANY,
            kind=wx.ITEM_NORMAL,text='Replay log')
        self.Bind(wx.EVT_MENU, log.ReplayLog, item)

    def _init_Exports(self,menu):
        '''Find exporter routines and add them into menus
        '''
        # set up the top-level menus
        projectmenu = wx.Menu()
        item = menu.AppendMenu(
            wx.ID_ANY, 'Entire project as',
            projectmenu, help='Export entire project')

        phasemenu = wx.Menu()
        item = menu.AppendMenu(
            wx.ID_ANY, 'Phase as',
            phasemenu, help='Export phase or sometimes phases')

        powdermenu = wx.Menu()
        item = menu.AppendMenu(
            wx.ID_ANY, 'Powder data as',
            powdermenu, help='Export powder diffraction histogram(s)')

        singlemenu = wx.Menu()
        item = menu.AppendMenu(
            wx.ID_ANY, 'Single crystal data as',
            singlemenu, help='Export single crystal histogram(s)')

        imagemenu = wx.Menu()
        item = menu.AppendMenu(
            wx.ID_ANY, 'Image data as',
            imagemenu, help='Export powder image(s) data')

        mapmenu = wx.Menu()
        item = menu.AppendMenu(
            wx.ID_ANY, 'Maps as',
            mapmenu, help='Export density map(s)')

        # pdfmenu = wx.Menu()
        # item = menu.AppendMenu(
        #     wx.ID_ANY, 'PDFs as',
        #     pdfmenu, help='Export pair distribution function(s)')

        # find all the exporter files
        pathlist = sys.path
        filelist = []
        for path in pathlist:
            for filename in glob.iglob(os.path.join(path,"G2export*.py")):
                filelist.append(filename)    
        filelist = sorted(list(set(filelist))) # remove duplicates
        self.exporterlist = []
        # go through the routines and import them, saving objects that
        # have export routines (method Exporter)
        for filename in filelist:
            path,rootname = os.path.split(filename)
            pkg = os.path.splitext(rootname)[0]
            try:
                fp = None
                fp, fppath,desc = imp.find_module(pkg,[path,])
                pkg = imp.load_module(pkg,fp,fppath,desc)
                for clss in inspect.getmembers(pkg): # find classes defined in package
                    if clss[0].startswith('_'): continue
                    if inspect.isclass(clss[1]):
                        # check if we have the required methods
                        for m in 'Exporter','loadParmDict':
                            if not hasattr(clss[1],m): break
                            if not callable(getattr(clss[1],m)): break
                        else:
                            exporter = clss[1](self) # create an export instance
                            self.exporterlist.append(exporter)
            except AttributeError:
                print 'Import_'+errprefix+': Attribute Error'+str(filename)
                pass
            except ImportError:
                print 'Import_'+errprefix+': Error importing file'+str(filename)
                pass
            if fp: fp.close()
        # Add submenu item(s) for each Exporter by its self-declared type (can be more than one)
        for obj in self.exporterlist:
            #print 'exporter',obj
            for typ in obj.exporttype:
                if typ == "project":
                    submenu = projectmenu
                elif typ == "phase":
                    submenu = phasemenu
                elif typ == "powder":
                    submenu = powdermenu
                elif typ == "single":
                    submenu = singlemenu
                elif typ == "image":
                    submenu = imagemenu
                elif typ == "map":
                    submenu = mapmenu
                # elif typ == "pdf":
                #     submenu = pdfmenu
                else:
                    print("Error, unknown type in "+str(obj))
                    break
                item = submenu.Append(
                    wx.ID_ANY,
                    help=obj.longFormatName,
                    kind=wx.ITEM_NORMAL,
                    text=obj.formatName)
                self.Bind(wx.EVT_MENU, obj.Exporter, id=item.GetId())
                self.ExportLookup[item.GetId()] = typ # lookup table for submenu item
        item = imagemenu.Append(wx.ID_ANY,
                        help='Export image controls and masks for multiple images',
                        kind=wx.ITEM_NORMAL,
                        text='Multiple image controls and masks')
        self.Bind(wx.EVT_MENU, self.OnSaveMultipleImg, id=item.GetId())
        #code to debug an Exporter. hard-code the routine below, to allow a reload before use
        # def DebugExport(event):
        #      print 'start reload'
        #      reload(G2IO)
        #      import G2export_pwdr as dev
        #      reload(dev)
        #      dev.ExportPowderFXYE(self).Exporter(event)
        # item = menu.Append(
        #     wx.ID_ANY,kind=wx.ITEM_NORMAL,
        #     help="debug exporter",text="test Export FXYE")
        # self.Bind(wx.EVT_MENU, DebugExport, id=item.GetId())
        # # #self.ExportLookup[item.GetId()] = 'image'
        # self.ExportLookup[item.GetId()] = 'powder'

    def _Add_ExportMenuItems(self,parent):
        # item = parent.Append(
        #     help='Select PWDR item to enable',id=wx.ID_ANY,
        #     kind=wx.ITEM_NORMAL,
        #     text='Export Powder Patterns...')
        # self.ExportPattern.append(item)
        # item.Enable(False)
        # self.Bind(wx.EVT_MENU, self.OnExportPatterns, id=item.GetId())

        item = parent.Append(
            help='',id=wx.ID_ANY,
            kind=wx.ITEM_NORMAL,
            text='Export All Peak Lists...')
        self.ExportPeakList.append(item)
        item.Enable(True)
        self.Bind(wx.EVT_MENU, self.OnExportPeakList, id=item.GetId())

        item = parent.Append(
            help='',id=wx.ID_ANY,
            kind=wx.ITEM_NORMAL,
            text='Export HKLs...')
        self.ExportHKL.append(item)
        self.Bind(wx.EVT_MENU, self.OnExportHKL, id=item.GetId())

        item = parent.Append(
            help='Select PDF item to enable',
            id=wx.ID_ANY,
            kind=wx.ITEM_NORMAL,
            text='Export PDF...')
        self.ExportPDF.append(item)
        item.Enable(False)
        self.Bind(wx.EVT_MENU, self.OnExportPDF, id=item.GetId())

    def FillMainMenu(self,menubar):
        '''Define contents of the main GSAS-II menu for the (main) data tree window.
        For the mac, this is also called for the data item windows as well so that
        the main menu items are data menu as well.
        '''
        File = wx.Menu(title='')
        menubar.Append(menu=File, title='&File')
        self._Add_FileMenuItems(File)
        Data = wx.Menu(title='')
        menubar.Append(menu=Data, title='Data')
        self._Add_DataMenuItems(Data)
        Calculate = wx.Menu(title='')        
        menubar.Append(menu=Calculate, title='&Calculate')
        self._Add_CalculateMenuItems(Calculate)
        Import = wx.Menu(title='')        
        menubar.Append(menu=Import, title='Import')
        self._Add_ImportMenu_Image(Import)
        self._Add_ImportMenu_Phase(Import)
        self._Add_ImportMenu_powder(Import)
        self._Add_ImportMenu_Sfact(Import)
        self._Add_ImportMenu_smallangle(Import)
        item = File.Append(wx.ID_PREFERENCES, text = "&Preferences")
        self.Bind(wx.EVT_MENU, self.OnPreferences, item)

        #======================================================================
        # Code to help develop/debug an importer, much is hard-coded below
        # but module is reloaded before each use, allowing faster testing
        # def DebugImport(event):
        #     print 'start reload'
        #     import G2phase_ISO as dev
        #     reload(dev)
        #     rd = dev.ISODISTORTPhaseReader()
        #     self.ImportMenuId[event.GetId()] = rd
        #     self.OnImportPhase(event)
            # or ----------------------------------------------------------------------
            #self.OnImportGeneric(rd,[],'test of ISODISTORTPhaseReader')
            # special debug code
            # or ----------------------------------------------------------------------
            # filename = '/Users/toby/projects/branton/subgroup_cif.txt'
            # fp = open(filename,'Ur')
            # if not rd.ContentsValidator(fp):
            #     print 'not validated'
            #     # make a list of used phase ranId's
            # phaseRIdList = []
            # sub = G2gd.GetPatternTreeItemId(self,self.root,'Phases')
            # if sub:
            #     item, cookie = self.PatternTree.GetFirstChild(sub)
            #     while item:
            #         phaseName = self.PatternTree.GetItemText(item)
            #         ranId = self.PatternTree.GetItemPyData(item).get('ranId')
            #         if ranId: phaseRIdList.append(ranId)
            #         item, cookie = self.PatternTree.GetNextChild(sub, cookie)
            # if rd.Reader(filename,fp,usedRanIdList=phaseRIdList):
            #     print 'read OK'
        # item = Import.Append(
        #     wx.ID_ANY,kind=wx.ITEM_NORMAL,
        #     help="debug importer",text="test importer")
        # self.Bind(wx.EVT_MENU, DebugImport, id=item.GetId())
        #======================================================================
        self.ExportMenu = wx.Menu(title='')
        menubar.Append(menu=self.ExportMenu, title='Export')
        self._init_Exports(self.ExportMenu)
        self._Add_ExportMenuItems(self.ExportMenu)
        if GSASIIpath.GetConfigValue('Enable_logging'):
            self.MacroMenu = wx.Menu(title='')
            menubar.Append(menu=self.MacroMenu, title='Macro')
            self._init_Macro()
        HelpMenu=G2G.MyHelp(self,helpType='Data tree',
            morehelpitems=[
                           ('&Tutorials','Tutorials'), 
                           ])
        menubar.Append(menu=HelpMenu,title='&Help')

    def _init_ctrls(self, parent):
        wx.Frame.__init__(self, name='GSASII', parent=parent,
            size=wx.Size(400, 250),style=wx.DEFAULT_FRAME_STYLE, title='GSAS-II data tree')
        clientSize = wx.ClientDisplayRect()
        Size = self.GetSize()
        xPos = clientSize[2]-Size[0]
        self.SetPosition(wx.Point(xPos,clientSize[1]))
        self._init_Imports()
        #initialize Menu item objects (these contain lists of menu items that are enabled or disabled)
        self.MakePDF = []
        self.Refine = []
        self.SeqRefine = [] # pointer(s) to Sequential Refinement menu objects
        #self.ExportPattern = []
        self.ExportPeakList = []
        self.ExportHKL = []
        self.ExportPDF = []
        self.ExportPhase = []
        self.ExportCIF = []
        #
        self.GSASIIMenu = wx.MenuBar()
        # create a list of all dataframe menus (appended in PrefillDataMenu)
        self.dataMenuBars = [self.GSASIIMenu]
        self.MacroStatusList = []
        self.FillMainMenu(self.GSASIIMenu)
        self.SetMenuBar(self.GSASIIMenu)
        self.Bind(wx.EVT_SIZE, self.OnSize)
        self.Status = self.CreateStatusBar()
        self.mainPanel = wx.Panel(self,-1)
        
        wxID_PATTERNTREE = wx.NewId()
        #self.PatternTree = wx.TreeCtrl(id=wxID_PATTERNTREE, # replaced for logging
        self.PatternTree = G2G.G2TreeCtrl(id=wxID_PATTERNTREE,
            parent=self.mainPanel, pos=wx.Point(0, 0),style=wx.TR_DEFAULT_STYLE )
        self.PatternTree.Bind(wx.EVT_TREE_SEL_CHANGED, self.OnPatternTreeSelChanged)
        self.PatternTree.Bind(wx.EVT_TREE_ITEM_COLLAPSED,
            self.OnPatternTreeItemCollapsed, id=wxID_PATTERNTREE)
        self.PatternTree.Bind(wx.EVT_TREE_ITEM_EXPANDED,
            self.OnPatternTreeItemExpanded, id=wxID_PATTERNTREE)
        self.PatternTree.Bind(wx.EVT_TREE_DELETE_ITEM,
            self.OnPatternTreeItemDelete, id=wxID_PATTERNTREE)
        self.PatternTree.Bind(wx.EVT_TREE_KEY_DOWN,
            self.OnPatternTreeKeyDown, id=wxID_PATTERNTREE)
        self.PatternTree.Bind(wx.EVT_TREE_BEGIN_RDRAG,
            self.OnPatternTreeBeginRDrag, id=wxID_PATTERNTREE)        
        self.PatternTree.Bind(wx.EVT_TREE_END_DRAG,
            self.OnPatternTreeEndDrag, id=wxID_PATTERNTREE)        
        #self.root = self.PatternTree.AddRoot('Loaded Data: ')
        self.root = self.PatternTree.root
        plotFrame = wx.Frame(None,-1,'GSASII Plots',size=wx.Size(700,600), \
            style=wx.DEFAULT_FRAME_STYLE ^ wx.CLOSE_BOX)
        self.G2plotNB = G2plt.G2PlotNoteBook(plotFrame)
        plotFrame.Show()
        
        self.dataDisplay = None
        
    def __init__(self, parent):
        self.ExportLookup = {}
        self._init_ctrls(parent)
        self.Image = wx.Image(
            os.path.join(GSASIIpath.path2GSAS2,'gsas2.ico'),
            wx.BITMAP_TYPE_ICO)
        if "wxMSW" in wx.PlatformInfo:
            img = self.Image.Scale(16, 16).ConvertToBitmap()
        elif "wxGTK" in wx.PlatformInfo:
            img = self.Image.Scale(22, 22).ConvertToBitmap()
        else:
            img = self.Image.ConvertToBitmap()
        self.SetIcon(wx.IconFromBitmap(img))
        self.Bind(wx.EVT_CLOSE, self.ExitMain)
        # various defaults
        self.oldFocus = None
        self.GSASprojectfile = ''
        self.dirname = os.path.expanduser('~')       #start in the users home directory by default; may be meaningless
        self.exportDir = None  # the last directory used for exports, if any.
        self.undofile = ''
        self.TreeItemDelete = False
        self.plotStyle = {'qPlot':False,'dPlot':False,'sqrtPlot':False}
        self.Weight = False
        self.IfPlot = False
        self.DDShowAll = False
        self.atmSel = ''
        self.PatternId = 0
        self.PickId = 0
        self.PickIdText = None
        self.PeakTable = []
        self.LimitsTable = []
        self.ifX20 = True   #use M20 /= (1+X20) in powder indexing, etc.
        self.HKL = []
        self.Lines = []
        self.itemPicked = None
        self.dataFrame = None
        self.Interpolate = 'nearest'
        self.ContourColor = 'Paired'
        self.VcovColor = 'RdYlGn'
        self.RamaColor = 'Blues'
        self.Projection = 'equal area'
        self.logPlot = False
        self.sqPlot = False
        self.ErrorBars = False
        self.Contour = False
        self.Legend = False
        self.SinglePlot = True
        self.SubBack = False
        self.seqReverse = False
        self.seqLines = True #draw lines between points
        self.plotView = 0
        self.Image = 0
        self.oldImagefile = '' # the name of the last image file read
        self.oldImageTag = None # the name of the tag for multi-image files
        self.ImageZ = []  # this contains the image plotted and used for integration
        # self.ImageZ and self.oldImagefile are set in GSASIIplot.PlotImage
        # and GSASIIIO.GetImageData
        # any changes to self.ImageZ should initialize self.oldImagefile to force a reread
        self.Integrate = 0
        self.imageDefault = {}
        self.IntgOutList = [] # list of integration tree item Ids created in G2IO.SaveIntegration
        self.AutointPWDRnames = [] # list of autoint created PWDR tree item names (to be deleted on a reset)
        self.autoIntFrame = None
        self.IntegratedList = [] # list of already integrated IMG tree items
        self.Sngl = False
        self.ifGetRing = False
        self.MaskKey = ''           #trigger for making image masks
        self.StrainKey = ''         #ditto for new strain d-zeros
        self.EnablePlot = True
        self.hist = ''              # selected histogram in Phase/Data tab
        if GSASIIpath.GetConfigValue('Starting_directory'):
            try: 
                os.chdir(GSASIIpath.GetConfigValue('Starting_directory'))
            except:
                print('Ignoring Config Starting_directory value: '+
                      GSASIIpath.GetConfigValue('Starting_directory'))
        arg = sys.argv
        if len(arg) > 1 and arg[1]:
            self.GSASprojectfile = os.path.splitext(arg[1])[0]+'.gpx'
            self.dirname = os.path.dirname(arg[1])
            if self.dirname: os.chdir(self.dirname)
            try:
                self.StartProject()         #open the file if possible
            except Exception:
                print 'Error opening or reading file',arg[1]
                import traceback
                print traceback.format_exc()
              
        self.ImportDir = os.path.normpath(os.getcwd()) # specifies a default path to be used for imports
        if GSASIIpath.GetConfigValue('Import_directory'):
            self.ImportDir = GSASIIpath.GetConfigValue('Import_directory')
            
    def GetTreeItemsList(self,item):
        return self.PatternTree._getTreeItemsList(item)

    def OnSize(self,event):
        'Called when the main window is resized. Not sure why'
        w,h = self.GetClientSizeTuple()
        self.mainPanel.SetSize(wx.Size(w,h))
        self.PatternTree.SetSize(wx.Size(w,h))
                        
    def OnPatternTreeSelChanged(self, event):
        '''Called when a data tree item is selected'''
        if self.TreeItemDelete:
            self.TreeItemDelete = False
        else:
            pltNum = self.G2plotNB.nb.GetSelection()
            if pltNum >= 0:                         #to avoid the startup with no plot!
                pltPage = self.G2plotNB.nb.GetPage(pltNum)
                pltPlot = pltPage.figure
            item = event.GetItem()
            G2gd.MovePatternTreeToGrid(self,item)
            if self.oldFocus:
                self.oldFocus.SetFocus()
        
    def OnPatternTreeItemCollapsed(self, event):
        'Called when a tree item is collapsed - all children will be collapsed'
        self.PatternTree.CollapseAllChildren(event.GetItem())

    def OnPatternTreeItemExpanded(self, event):
        'Called when a tree item is expanded'
        self.OnPatternTreeSelChanged(event)
        event.Skip()
        
    def OnPatternTreeItemDelete(self, event):
        'Called when a tree item is deleted -- not sure what this does'
        self.TreeItemDelete = True

    def OnPatternTreeItemActivated(self, event):
        'Called when a tree item is activated'
        event.Skip()
        
    def OnPatternTreeBeginRDrag(self,event):
        event.Allow()
        self.BeginDragId = event.GetItem()
        self.ParentId = self.PatternTree.GetItemParent(self.BeginDragId)
        DragText = self.PatternTree.GetItemText(self.BeginDragId)
        self.DragData = [[DragText,self.PatternTree.GetItemPyData(self.BeginDragId)],]
        item, cookie = self.PatternTree.GetFirstChild(self.BeginDragId)
        while item:     #G2 data tree has no sub children under a child of a tree item
            name = self.PatternTree.GetItemText(item)
            self.DragData.append([name,self.PatternTree.GetItemPyData(item)])
            item, cookie = self.PatternTree.GetNextChild(self.BeginDragId, cookie)                            
        
    def OnPatternTreeEndDrag(self,event):
        event.Allow()
        self.EndDragId = event.GetItem()
        try:
            NewParent = self.PatternTree.GetItemParent(self.EndDragId)
        except:
            self.EndDragId = self.PatternTree.GetLastChild(self.root)
            NewParent = self.root
        if self.ParentId != NewParent:
            self.ErrorDialog('Drag not allowed','Wrong parent for item dragged')
        else:
            Name,Item = self.DragData[0]
            NewId = self.PatternTree.InsertItem(self.ParentId,self.EndDragId,Name,data=None)
            self.PatternTree.SetItemPyData(NewId,Item)
            for name,item in self.DragData[1:]:     #loop over children
                Id = self.PatternTree.AppendItem(parent=NewId,text=name)
                self.PatternTree.SetItemPyData(Id,item)
            self.PatternTree.Delete(self.BeginDragId)
            G2gd.MovePatternTreeToGrid(self,NewId)
        
    def OnPatternTreeKeyDown(self,event): #doesn't exactly work right with Shift key down
        'Allows stepping through the tree with the up/down arrow keys'
        self.oldFocus = wx.Window.FindFocus()
        keyevt = event.GetKeyEvent()
        key = event.GetKeyCode()
        item = self.PatternTree.GetSelection()
        if type(item) is int: return # is this the toplevel in tree?
        name = self.PatternTree.GetItemText(item)
        parent = self.PatternTree.GetItemParent(item)
        if key == wx.WXK_UP:
            if keyevt.GetModifiers() == wx.MOD_SHIFT and parent != self.root:
                if type(parent) is int: return # is this the toplevel in tree?
                prev = self.PatternTree.GetPrevSibling(parent)
                NewId = G2gd.GetPatternTreeItemId(self,prev,name)
                if NewId:
                    self.PatternTree.Collapse(parent)
                    self.PatternTree.Expand(prev)
                    self.oldFocus = wx.Window.FindFocus()
                    wx.CallAfter(self.PatternTree.SelectItem,NewId)
                else:
                    wx.CallAfter(self.PatternTree.SelectItem,item)
            else:    
                self.PatternTree.GetPrevSibling(item)
                self.PatternTree.SelectItem(item)
        elif key == wx.WXK_DOWN:
            if keyevt.GetModifiers() == wx.MOD_SHIFT and parent != self.root:
                next = self.PatternTree.GetNextSibling(parent)
                NewId = G2gd.GetPatternTreeItemId(self,next,name)
                if NewId:
                    self.PatternTree.Collapse(parent)
                    self.PatternTree.Expand(next)
                    self.oldFocus = wx.Window.FindFocus()
                    wx.CallAfter(self.PatternTree.SelectItem,NewId)
                else:
                    wx.CallAfter(self.PatternTree.SelectItem,item)
            else:    
                self.PatternTree.GetNextSibling(item)
                self.PatternTree.SelectItem(item)
                
    def OnReadPowderPeaks(self,event):
        'Bound to menu Data/Read Powder Peaks'
        Cuka = 1.54052
        self.CheckNotebook()
        dlg = wx.FileDialog(self, 'Choose file with peak list', '.', '', 
            'peak files (*.txt)|*.txt|All files (*.*)|*.*',wx.OPEN|wx.CHANGE_DIR)
        try:
            if dlg.ShowModal() == wx.ID_OK:
                self.HKL = []
                self.powderfile = dlg.GetPath()
                comments,peaks,limits,wave = G2IO.GetPowderPeaks(self.powderfile)
                Id = self.PatternTree.AppendItem(parent=self.root,text='PKS '+os.path.basename(self.powderfile))
                data = ['PKS',wave,0.0]
                names = ['Type','Lam','Zero'] 
                codes = [0,0,0]
                inst = [G2IO.makeInstDict(names,data,codes),{}]
                self.PatternTree.SetItemPyData(self.PatternTree.AppendItem(Id,text='Instrument Parameters'),inst)
                self.PatternTree.SetItemPyData(self.PatternTree.AppendItem(Id,text='Comments'),comments)
                self.PatternTree.SetItemPyData(self.PatternTree.AppendItem(Id,text='Limits'),[tuple(limits),limits])
                self.PatternTree.SetItemPyData(self.PatternTree.AppendItem(Id,text='Index Peak List'),[peaks,[]])
                self.PatternTree.SetItemPyData(self.PatternTree.AppendItem(Id,text='Unit Cells List'),[])             
                self.PatternTree.Expand(Id)
                self.PatternTree.SelectItem(Id)
                os.chdir(dlg.GetDirectory())           # to get Mac/Linux to change directory!
        finally:
            dlg.Destroy()
                        
    def OnImageRead(self,event):
        '''Called to read in an image in any known format. *** Depreciated. ***
        Note: When removed, G2IO.ReadLoadImage can also be removed
        '''
        G2G.G2MessageBox(self,'Please use the Import/Image/... menu item rather than this','depreciating menu item')
        self.CheckNotebook()
        dlg = wx.FileDialog(
            self, 'Choose image files', '.', '',
            'Any supported image file (*.edf;*.tif;*.tiff;*.mar*;*.ge*;*.avg;*.sum;*.img;*.stl;*.G2img;*.png)|'
            '*.edf;*.tif;*.tiff;*.mar*;*.ge*;*.avg;*.sum;*.img;*.stl;*.G2img;*.png;*.zip|'
            'European detector file (*.edf)|*.edf|'
            'Any detector tif (*.tif;*.tiff)|*.tif;*.tiff|'
            'MAR file (*.mar*)|*.mar*|'
            'GE Image (*.ge*;*.avg;*.sum)|*.ge*;*.avg;*.sum|'
            'ADSC Image (*.img)|*.img|'
            'Rigaku R-Axis IV (*.stl)|*.stl|'
            'GSAS-II Image (*.G2img)|*.G2img|'
            'Portable Network Graphics image (*.png)|*.png|'
            'Zip archive (*.zip)|*.zip|'
            'All files (*.*)|*.*',
            wx.OPEN | wx.MULTIPLE|wx.CHANGE_DIR)
        try:
            if dlg.ShowModal() == wx.ID_OK:
                imagefiles = dlg.GetPaths()
                imagefiles.sort()
                for imagefile in imagefiles:
                    G2IO.ReadLoadImage(imagefile,self)
                os.chdir(dlg.GetDirectory())           # to get Mac/Linux to change directory!
                self.PatternTree.SelectItem(G2gd.GetPatternTreeItemId(self,self.Image,'Image Controls'))             #show last image to be read
        finally:
            path = dlg.GetDirectory()           # to get Mac/Linux to change directory!
            os.chdir(path)
            dlg.Destroy()

    def CheckNotebook(self):
        '''Make sure the data tree has the minimally expected controls.
        '''
        if not G2gd.GetPatternTreeItemId(self,self.root,'Notebook'):
            sub = self.PatternTree.AppendItem(parent=self.root,text='Notebook')
            self.PatternTree.SetItemPyData(sub,[''])
        if not G2gd.GetPatternTreeItemId(self,self.root,'Controls'):
            sub = self.PatternTree.AppendItem(parent=self.root,text='Controls')
            self.PatternTree.SetItemPyData(sub,copy.copy(G2obj.DefaultControls))
        if not G2gd.GetPatternTreeItemId(self,self.root,'Covariance'):
            sub = self.PatternTree.AppendItem(parent=self.root,text='Covariance')
            self.PatternTree.SetItemPyData(sub,{})
        if not G2gd.GetPatternTreeItemId(self,self.root,'Constraints'):
            sub = self.PatternTree.AppendItem(parent=self.root,text='Constraints')
            self.PatternTree.SetItemPyData(sub,{'Hist':[],'HAP':[],'Phase':[]})
        if not G2gd.GetPatternTreeItemId(self,self.root,'Restraints'):
            sub = self.PatternTree.AppendItem(parent=self.root,text='Restraints')
            self.PatternTree.SetItemPyData(sub,{})
        if not G2gd.GetPatternTreeItemId(self,self.root,'Rigid bodies'):
            sub = self.PatternTree.AppendItem(parent=self.root,text='Rigid bodies')
            self.PatternTree.SetItemPyData(sub,{'Vector':{'AtInfo':{}},
                'Residue':{'AtInfo':{}},'RBIds':{'Vector':[],'Residue':[]}})
                
    class CopyDialog(wx.Dialog):
        '''Creates a dialog for copying control settings between
        data tree items'''
        def __init__(self,parent,title,text,data):
            wx.Dialog.__init__(self,parent,-1,title, 
                pos=wx.DefaultPosition,style=wx.DEFAULT_DIALOG_STYLE)
            self.data = data
            panel = wx.Panel(self)
            mainSizer = wx.BoxSizer(wx.VERTICAL)
            topLabl = wx.StaticText(panel,-1,text)
            mainSizer.Add((10,10),1)
            mainSizer.Add(topLabl,0,wx.ALIGN_CENTER_VERTICAL|wx.LEFT,10)
            mainSizer.Add((10,10),1)
            ncols = len(data)/40+1
            dataGridSizer = wx.FlexGridSizer(cols=ncols,hgap=2,vgap=2)
            for id,item in enumerate(self.data):
                ckbox = wx.CheckBox(panel,id,item[1])
                ckbox.Bind(wx.EVT_CHECKBOX,self.OnCopyChange)                    
                dataGridSizer.Add(ckbox,0,wx.LEFT,10)
            mainSizer.Add(dataGridSizer,0,wx.EXPAND)
            OkBtn = wx.Button(panel,-1,"Ok")
            OkBtn.Bind(wx.EVT_BUTTON, self.OnOk)
            cancelBtn = wx.Button(panel,-1,"Cancel")
            cancelBtn.Bind(wx.EVT_BUTTON, self.OnCancel)
            btnSizer = wx.BoxSizer(wx.HORIZONTAL)
            btnSizer.Add((20,20),1)
            btnSizer.Add(OkBtn)
            btnSizer.Add((20,20),1)
            btnSizer.Add(cancelBtn)
            btnSizer.Add((20,20),1)
            
            mainSizer.Add(btnSizer,0,wx.EXPAND|wx.BOTTOM|wx.TOP, 10)
            panel.SetSizer(mainSizer)
            panel.Fit()
            self.Fit()
        
        def OnCopyChange(self,event):
            id = event.GetId()
            self.data[id][0] = self.FindWindowById(id).GetValue()        
            
        def OnOk(self,event):
            parent = self.GetParent()
            parent.Raise()
            self.EndModal(wx.ID_OK)              
            
        def OnCancel(self,event):
            parent = self.GetParent()
            parent.Raise()
            self.EndModal(wx.ID_CANCEL)              
            
        def GetData(self):
            return self.data
        
    class SumDialog(wx.Dialog):
        'Allows user to supply scale factor(s) when summing data'
        def __init__(self,parent,title,text,dataType,data):
            wx.Dialog.__init__(self,parent,-1,title, 
                pos=wx.DefaultPosition,style=wx.DEFAULT_DIALOG_STYLE)
            self.data = data
            panel = wx.Panel(self)
            mainSizer = wx.BoxSizer(wx.VERTICAL)
            topLabl = wx.StaticText(panel,-1,text)
            mainSizer.Add((10,10),1)
            mainSizer.Add(topLabl,0,wx.ALIGN_CENTER_VERTICAL|wx.LEFT,10)
            mainSizer.Add((10,10),1)
            dataGridSizer = wx.FlexGridSizer(cols=2,hgap=2,vgap=2)
            for id,item in enumerate(self.data[:-1]):
                name = wx.TextCtrl(panel,-1,item[1],size=wx.Size(200,20))
                name.SetEditable(False)
                scale = wx.TextCtrl(panel,id,'%.3f'%(item[0]),style=wx.TE_PROCESS_ENTER)
                scale.Bind(wx.EVT_TEXT_ENTER,self.OnScaleChange)
                scale.Bind(wx.EVT_KILL_FOCUS,self.OnScaleChange)
                dataGridSizer.Add(scale,0,wx.LEFT,10)
                dataGridSizer.Add(name,0,wx.RIGHT,10)
            if dataType:
                dataGridSizer.Add(wx.StaticText(panel,-1,'Sum result name: '+dataType),0, \
                    wx.LEFT|wx.TOP|wx.ALIGN_CENTER_VERTICAL,10)
                self.name = wx.TextCtrl(panel,-1,self.data[-1],size=wx.Size(200,20),style=wx.TE_PROCESS_ENTER)
                self.name.Bind(wx.EVT_TEXT_ENTER,self.OnNameChange)
                self.name.Bind(wx.EVT_KILL_FOCUS,self.OnNameChange)
                dataGridSizer.Add(self.name,0,wx.RIGHT|wx.TOP,10)
            mainSizer.Add(dataGridSizer,0,wx.EXPAND)
            OkBtn = wx.Button(panel,-1,"Ok")
            OkBtn.Bind(wx.EVT_BUTTON, self.OnOk)
            cancelBtn = wx.Button(panel,-1,"Cancel")
            cancelBtn.Bind(wx.EVT_BUTTON, self.OnCancel)
            btnSizer = wx.BoxSizer(wx.HORIZONTAL)
            btnSizer.Add((20,20),1)
            btnSizer.Add(OkBtn)
            btnSizer.Add((20,20),1)
            btnSizer.Add(cancelBtn)
            btnSizer.Add((20,20),1)
            
            mainSizer.Add(btnSizer,0,wx.EXPAND|wx.BOTTOM|wx.TOP, 10)
            panel.SetSizer(mainSizer)
            panel.Fit()
            self.Fit()

        def OnScaleChange(self,event):
            id = event.GetId()
            value = self.FindWindowById(id).GetValue()
            try:
                self.data[id][0] = float(value)
                self.FindWindowById(id).SetValue('%.3f'%(self.data[id][0]))
            except ValueError:
                if value and '-' not in value[0]:
                    print 'bad input - numbers only'
                    self.FindWindowById(id).SetValue('0.000')
            
        def OnNameChange(self,event):
            self.data[-1] = self.name.GetValue() 
            
        def OnOk(self,event):
            parent = self.GetParent()
            parent.Raise()
            self.EndModal(wx.ID_OK)              
            
        def OnCancel(self,event):
            parent = self.GetParent()
            parent.Raise()
            self.EndModal(wx.ID_CANCEL)              
            
        def GetData(self):
            return self.data
                        
    def OnPwdrSum(self,event):
        'Sum together powder data(?)'
        TextList = []
        DataList = []
        SumList = []
        Names = []
        Inst = None
        SumItemList = []
        Comments = ['Sum equals: \n']
        if self.PatternTree.GetCount():
            item, cookie = self.PatternTree.GetFirstChild(self.root)
            while item:
                name = self.PatternTree.GetItemText(item)
                Names.append(name)
                if 'PWDR' in name:
                    TextList.append([0.0,name])
                    DataList.append(self.PatternTree.GetItemPyData(item)[1])    # (x,y,w,yc,yb,yd)
                    if not Inst:
                        Inst = self.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(self,item, 'Instrument Parameters'))
                item, cookie = self.PatternTree.GetNextChild(self.root, cookie)
            if len(TextList) < 2:
                self.ErrorDialog('Not enough data to sum','There must be more than one "PWDR" pattern')
                return
            TextList.append('default_sum_name')                
            dlg = self.SumDialog(self,'Sum data','Enter scale for each pattern in summation','PWDR',TextList)
            try:
                if dlg.ShowModal() == wx.ID_OK:
                    lenX = 0
                    Xminmax = [0,0]
                    Xsum = []
                    Ysum = []
                    Vsum = []
                    result = dlg.GetData()
                    for i,item in enumerate(result[:-1]):
                        scale,name = item
                        data = DataList[i]
                        if scale:
                            Comments.append("%10.3f %s" % (scale,' * '+name))
                            x,y,w,yc,yb,yd = data   #numpy arrays!
                            v = 1./w
                            if lenX:
                                if lenX != len(x):
                                    self.ErrorDialog('Data length error','Data to be summed must have same number of points'+ \
                                        '\nExpected:'+str(lenX)+ \
                                        '\nFound:   '+str(len(x))+'\nfor '+name)
                                    return
                            else:
                                lenX = len(x)
                            if Xminmax[1]:
                                if Xminmax != [x[0],x[-1]]:
                                    self.ErrorDialog('Data range error','Data to be summed must span same range'+ \
                                        '\nExpected:'+str(Xminmax[0])+' '+str(Xminmax[1])+ \
                                        '\nFound:   '+str(x[0])+' '+str(x[-1])+'\nfor '+name)
                                    return
                                else:
                                    for j,yi in enumerate(y):
                                         Ysum[j] += scale*yi
                                         Vsum[j] += abs(scale)*v[j]
                            else:
                                Xminmax = [x[0],x[-1]]
                                YCsum = YBsum = YDsum = [0.0 for i in range(lenX)]
                                for j,yi in enumerate(y):
                                    Xsum.append(x[j])
                                    Ysum.append(scale*yi)
                                    Vsum.append(abs(scale*v[j]))
                    Wsum = 1./np.array(Vsum)
                    outname = 'PWDR '+result[-1]
                    Id = 0
                    if outname in Names:
                        dlg2 = wx.MessageDialog(self,'Overwrite data?','Duplicate data name',wx.OK|wx.CANCEL)
                        try:
                            if dlg2.ShowModal() == wx.ID_OK:
                                Id = G2gd.GetPatternTreeItemId(self,self.root,name)
                                self.PatternTree.Delete(Id)
                        finally:
                            dlg2.Destroy()
                    Id = self.PatternTree.AppendItem(parent=self.root,text=outname)
                    if Id:
                        Sample = G2pdG.SetDefaultSample()
                        valuesdict = {
                            'wtFactor':1.0,
                            'Dummy':False,
                            'ranId':ran.randint(0,sys.maxint),
                            'Offset':[0.0,0.0],'delOffset':0.02,'refOffset':-1.0,'refDelt':0.01,
                            'qPlot':False,'dPlot':False,'sqrtPlot':False
                            }
                        self.PatternTree.SetItemPyData(Id,[valuesdict,[np.array(Xsum),np.array(Ysum),np.array(Wsum),
                            np.array(YCsum),np.array(YBsum),np.array(YDsum)]])
                        self.PatternTree.SetItemPyData(self.PatternTree.AppendItem(Id,text='Comments'),Comments)                    
                        self.PatternTree.SetItemPyData(self.PatternTree.AppendItem(Id,text='Limits'),[tuple(Xminmax),Xminmax])
                        self.PatternTree.SetItemPyData(self.PatternTree.AppendItem(Id,text='Background'),[['chebyschev',True,3,1.0,0.0,0.0],
                            {'nDebye':0,'debyeTerms':[],'nPeaks':0,'peaksList':[]}])
                        self.PatternTree.SetItemPyData(self.PatternTree.AppendItem(Id,text='Instrument Parameters'),Inst)
                        self.PatternTree.SetItemPyData(self.PatternTree.AppendItem(Id,text='Sample Parameters'),Sample)
                        self.PatternTree.SetItemPyData(self.PatternTree.AppendItem(Id,text='Peak List'),{'peaks':[],'sigDict':{}})
                        self.PatternTree.SetItemPyData(self.PatternTree.AppendItem(Id,text='Index Peak List'),[[],[]])
                        self.PatternTree.SetItemPyData(self.PatternTree.AppendItem(Id,text='Unit Cells List'),[])             
                        self.PatternTree.SetItemPyData(self.PatternTree.AppendItem(Id,text='Reflection Lists'),{})             
                        self.PatternTree.SelectItem(Id)
                        self.PatternTree.Expand(Id)
            finally:
                dlg.Destroy()

    def OnImageSum(self,event):
        'Sum together image data(?)'
        TextList = []
        DataList = []
        SumList = []
        Names = []
        Inst = []
        SumItemList = []
        Comments = ['Sum equals: \n']
        if self.PatternTree.GetCount():
            item, cookie = self.PatternTree.GetFirstChild(self.root)
            while item:
                name = self.PatternTree.GetItemText(item)
                Names.append(name)
                if 'IMG' in name:
                    TextList.append([0.0,name])
                    DataList.append(self.PatternTree.GetImageLoc(item))        #Size,Image,Tag
                    Data = self.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(self,item,'Image Controls'))
                item, cookie = self.PatternTree.GetNextChild(self.root, cookie)
            if len(TextList) < 2:
                self.ErrorDialog('Not enough data to sum','There must be more than one "IMG" pattern')
                return
            TextList.append('default_sum_name')                
            dlg = self.SumDialog(self,'Sum data','Enter scale for each image in summation','IMG',TextList)
            try:
                if dlg.ShowModal() == wx.ID_OK:
                    imSize = 0
                    result = dlg.GetData()
                    First = True
                    Found = False
                    for i,item in enumerate(result[:-1]):
                        scale,name = item
                        if scale:
                            Found = True                                
                            Comments.append("%10.3f %s" % (scale,' * '+name))
                            Npix,imagefile,imagetag = DataList[i]
                            image = G2IO.GetImageData(self,imagefile,imageOnly=True,ImageTag=imagetag)
                            if First:
                                newImage = np.zeros_like(image)
                                First = False
                            if imSize:
                                if imSize != Npix:
                                    self.ErrorDialog('Image size error','Images to be summed must be same size'+ \
                                        '\nExpected:'+str(imSize)+ \
                                        '\nFound:   '+str(Npix)+'\nfor '+name)
                                    return
                                newImage = newImage+scale*image
                            else:
                                imSize = Npix
                                newImage = newImage+scale*image
                            del(image)
                    if not Found:
                        self.ErrorDialog('Image sum error','No nonzero image multipliers found')
                        return
                        
                    newImage = np.asfarray(newImage,dtype=np.float32)                        
                    outname = 'IMG '+result[-1]
                    Id = 0
                    if outname in Names:
                        dlg2 = wx.MessageDialog(self,'Overwrite data?','Duplicate data name',wx.OK|wx.CANCEL)
                        try:
                            if dlg2.ShowModal() == wx.ID_OK:
                                Id = G2gd.GetPatternTreeItemId(self,self.root,name)
                        finally:
                            dlg2.Destroy()
                    else:
                        Id = self.PatternTree.AppendItem(parent=self.root,text=outname)
                    if Id:
                        dlg = wx.FileDialog(self, 'Choose sum image filename', '.', '', 
                            'G2img files (*.G2img)|*.G2img', 
                            wx.SAVE|wx.FD_OVERWRITE_PROMPT|wx.CHANGE_DIR)
                        if dlg.ShowModal() == wx.ID_OK:
                            newimagefile = dlg.GetPath()
                            newimagefile = G2IO.FileDlgFixExt(dlg,newimagefile)
                            G2IO.PutG2Image(newimagefile,Comments,Data,Npix,newImage)
                            Imax = np.amax(newImage)
                            Imin = np.amin(newImage)
                            newImage = []
                            self.PatternTree.SetItemPyData(Id,[imSize,newimagefile])
                            self.PatternTree.SetItemPyData(self.PatternTree.AppendItem(Id,text='Comments'),Comments)
                        del(newImage)
                        if self.imageDefault:
                            Data = copy.copy(self.imageDefault)
                        Data['showLines'] = True
                        Data['ring'] = []
                        Data['rings'] = []
                        Data['cutoff'] = 10
                        Data['pixLimit'] = 20
                        Data['ellipses'] = []
                        Data['calibrant'] = ''
                        Data['range'] = [(Imin,Imax),[Imin,Imax]]
                        self.PatternTree.SetItemPyData(self.PatternTree.AppendItem(Id,text='Image Controls'),Data)                                            
                        Masks = {'Points':[],'Rings':[],'Arcs':[],'Polygons':[],'Frames':[],'Thresholds':[(Imin,Imax),[Imin,Imax]]}
                        self.PatternTree.SetItemPyData(self.PatternTree.AppendItem(Id,text='Masks'),Masks)
                        self.PatternTree.SetItemPyData(self.PatternTree.AppendItem(Id,text='Stress/Strain'),
                            {'Type':'True','d-zero':[],'Sample phi':0.0,'Sample z':0.0,'Sample load':0.0})
                        self.PatternTree.SelectItem(Id)
                        self.PatternTree.Expand(Id)
                        self.PickId = G2gd.GetPatternTreeItemId(self,self.root,outname)
                        self.Image = self.PickId
            finally:
                dlg.Destroy()
                      
    def OnAddPhase(self,event):
        'Add a new, empty phase to the tree. Called by Data/Add Phase menu'
        if not G2gd.GetPatternTreeItemId(self,self.root,'Phases'):
            sub = self.PatternTree.AppendItem(parent=self.root,text='Phases')
        else:
            sub = G2gd.GetPatternTreeItemId(self,self.root,'Phases')
        PhaseName = ''
        dlg = wx.TextEntryDialog(None,'Enter a name for this phase','Phase Name Entry','New phase',
            style=wx.OK)
        if dlg.ShowModal() == wx.ID_OK:
            PhaseName = dlg.GetValue()
        dlg.Destroy()
        sub = self.PatternTree.AppendItem(parent=sub,text=PhaseName)
        E,SGData = G2spc.SpcGroup('P 1')
        self.PatternTree.SetItemPyData(sub,G2IO.SetNewPhase(Name=PhaseName,SGData=SGData))
        
    def OnDeletePhase(self,event):
        'Delete a phase from the tree. Called by Data/Delete Phase menu'
        #Hmm, also need to delete this phase from Reflection Lists for each PWDR histogram
        if self.dataFrame:
            self.dataFrame.Clear() 
        TextList = []
        DelList = []
        DelItemList = []
        if G2gd.GetPatternTreeItemId(self,self.root,'Phases'):
            sub = G2gd.GetPatternTreeItemId(self,self.root,'Phases')
        else:
            return
        if sub:
            item, cookie = self.PatternTree.GetFirstChild(sub)
            while item:
                TextList.append(self.PatternTree.GetItemText(item))
                item, cookie = self.PatternTree.GetNextChild(sub, cookie)                
            dlg = wx.MultiChoiceDialog(self, 'Which phase to delete?', 'Delete phase', TextList, wx.CHOICEDLG_STYLE)
            try:
                if dlg.ShowModal() == wx.ID_OK:
                    result = dlg.GetSelections()
                    for i in result: DelList.append([i,TextList[i]])
                    item, cookie = self.PatternTree.GetFirstChild(sub)
                    i = 0
                    while item:
                        if [i,self.PatternTree.GetItemText(item)] in DelList: DelItemList.append(item)
                        item, cookie = self.PatternTree.GetNextChild(sub, cookie)
                        i += 1
                    for item in DelItemList:
                        name = self.PatternTree.GetItemText(item)
                        self.PatternTree.Delete(item)
                        self.G2plotNB.Delete(name)
                    item, cookie = self.PatternTree.GetFirstChild(self.root)
                    while item:
                        name = self.PatternTree.GetItemText(item)
                        if 'PWDR' in name:
                            Id = G2gd.GetPatternTreeItemId(self,item, 'Reflection Lists')
                            refList = self.PatternTree.GetItemPyData(Id)
                            if len(refList):
                                for i,item in DelList:
                                    if item in refList:
                                        del(refList[item])
                            self.PatternTree.SetItemPyData(Id,refList)
                        elif 'HKLF' in name:
                            data = self.PatternTree.GetItemPyData(item)
                            data[0] = {}
                            self.PatternTree.SetItemPyData(item,data)
                            
                        item, cookie = self.PatternTree.GetNextChild(self.root, cookie)
            finally:
                dlg.Destroy()
                
    def OnRenameData(self,event):
        'Renames an existing phase. Called by Data/Rename Phase menu'
        name = self.PatternTree.GetItemText(self.PickId)      
        if 'PWDR' in name or 'HKLF' in name or 'IMG' in name:
            dataType = name[:name.index(' ')+1]                 #includes the ' '
            dlg = wx.TextEntryDialog(self,'Data name: '+dataType,'Change data name',
                defaultValue=name[name.index(' ')+1:])
            try:
                if dlg.ShowModal() == wx.ID_OK:
                    self.PatternTree.SetItemText(self.PickId,dataType+dlg.GetValue())
            finally:
                dlg.Destroy()
        
    def GetFileList(self,fileType,skip=None):        #potentially useful?
        'Appears unused. Note routine of same name in GSASIIpwdGUI'
        fileList = []
        Source = ''
        id, cookie = self.PatternTree.GetFirstChild(self.root)
        while id:
            name = self.PatternTree.GetItemText(id)
            if fileType in name:
                if id == skip:
                    Source = name
                else:
                    fileList.append([False,name,id])
            id, cookie = self.PatternTree.GetNextChild(self.root, cookie)
        if skip:
            return fileList,Source
        else:
            return fileList
            
    def OnDataDelete(self, event):
        '''Delete one or more histograms from data tree. Called by the
        Data/DeleteData menu
        '''
        TextList = []
        DelList = []
        DelItemList = []
        nItems = {'PWDR':0,'SASD':0,'IMG':0,'HKLF':0,'PDF':0}
        ifPWDR = False
        ifSASD = False
        ifIMG = False
        ifHKLF = False
        ifPDF = False
        if self.PatternTree.GetCount():
            item, cookie = self.PatternTree.GetFirstChild(self.root)
            while item:
                name = self.PatternTree.GetItemText(item)
                if name not in ['Notebook','Controls','Covariance','Constraints',
                    'Restraints','Phases','Rigid bodies','Sequential results']:
                    if 'PWDR' in name: ifPWDR = True; nItems['PWDR'] += 1
                    if 'SASD' in name: ifSASD = True; nItems['SASD'] += 1
                    if 'IMG' in name: ifIMG = True; nItems['IMG'] += 1
                    if 'HKLF' in name: ifHKLF = True; nItems['HKLF'] += 1
                    if 'PDF' in name: ifPDF = True; nItems['PDF'] += 1
                    TextList.append(name)
                item, cookie = self.PatternTree.GetNextChild(self.root, cookie)
            dlg = G2G.G2MultiChoiceDialog(self, 'Which data to delete?', 'Delete data', TextList, wx.CHOICEDLG_STYLE)
            try:
                if dlg.ShowModal() == wx.ID_OK:
                    result = dlg.GetSelections()
                    for i in result: DelList.append(TextList[i])
                    item, cookie = self.PatternTree.GetFirstChild(self.root)
                    while item:
                        itemName = self.PatternTree.GetItemText(item)
                        if itemName in DelList:
                            if 'PWDR' in itemName: nItems['PWDR'] -= 1
                            elif 'SASD' in itemName: nItems['SASD'] -= 1
                            elif 'IMG' in itemName: nItems['IMG'] -= 1
                            elif 'HKLF' in itemName: nItems['HKLF'] -= 1
                            elif 'PDF' in itemName: nItems['PDF'] -= 1
                            DelItemList.append(item)
                        item, cookie = self.PatternTree.GetNextChild(self.root, cookie)
                    for item in DelItemList:
                        self.PatternTree.Delete(item)
                    self.PickId = 0
                    self.PickIdText = None
                    self.PatternId = 0
                    if nItems['PWDR']:
                        wx.CallAfter(G2plt.PlotPatterns,self,True)
                    else:
                        self.G2plotNB.Delete('Powder Patterns')
                    if not nItems['IMG']:
                        self.G2plotNB.Delete('2D Powder Image')
                    if not nItems['HKLF']:
                        self.G2plotNB.Delete('Structure Factors')
                        if '3D Structure Factors' in self.G2plotNB.plotList:
                            self.G2plotNB.Delete('3D Structure Factors')
            finally:
                dlg.Destroy()

    def OnFileOpen(self, event, filename=None):
        '''Gets a GSAS-II .gpx project file in response to the
        File/Open Project menu button
        '''
        result = wx.ID_OK
        self.EnablePlot = False
        if self.PatternTree.GetChildrenCount(self.root,False):
            if self.dataFrame:
                self.dataFrame.Clear() 
            dlg = wx.MessageDialog(
                self,
                'Do you want to overwrite the current project? '
                'Any unsaved changes will be lost. Press OK to continue.',
                'Overwrite?',  wx.OK | wx.CANCEL)
            try:
                result = dlg.ShowModal()
                if result == wx.ID_OK:
                    self.PatternTree.DeleteChildren(self.root)
                    self.GSASprojectfile = ''
                    self.HKL = []
                    if self.G2plotNB.plotList:
                        self.G2plotNB.clear()
            finally:
                dlg.Destroy()
        if result != wx.ID_OK: return

        if not filename:
            if self.dataDisplay: self.dataDisplay.Destroy()
            dlg = wx.FileDialog(self, 'Choose GSAS-II project file', '.', '', 
                'GSAS-II project file (*.gpx)|*.gpx',wx.OPEN|wx.CHANGE_DIR)
            try:
                if dlg.ShowModal() != wx.ID_OK: return
                self.GSASprojectfile = dlg.GetPath()
                self.GSASprojectfile = G2IO.FileDlgFixExt(dlg,self.GSASprojectfile)
                self.dirname = dlg.GetDirectory()
            finally:
                dlg.Destroy()
        else:
            self.GSASprojectfile = os.path.splitext(filename)[0]+'.gpx'
            self.dirname = os.path.split(filename)[0]

        try:
            self.StartProject()         #open the file if possible
        except:
            print '\nError opening file ',filename
            import traceback
            print traceback.format_exc()
        
    def StartProject(self):
        '''Opens a GSAS-II project file & selects the 1st available data set to 
        display (PWDR, HKLF or SASD)
        '''
        
        Id = 0
        G2IO.ProjFileOpen(self)
        self.PatternTree.SetItemText(self.root,'Loaded Data: '+self.GSASprojectfile)
        self.PatternTree.Expand(self.root)
        self.HKL = []
        item, cookie = self.PatternTree.GetFirstChild(self.root)
        while item and not Id:
            name = self.PatternTree.GetItemText(item)
            if name[:4] in ['PWDR','HKLF','IMG ','PDF ','SASD',]:
                Id = item
            elif name == 'Controls':
                data = self.PatternTree.GetItemPyData(item)
                if data:
                    for item in self.Refine: item.Enable(True)
                    self.EnableSeqRefineMenu()
            item, cookie = self.PatternTree.GetNextChild(self.root, cookie)
        if Id:
            self.EnablePlot = True
            self.PatternTree.SelectItem(Id)
        self.CheckNotebook()
        if self.dirname: os.chdir(self.dirname)           # to get Mac/Linux to change directory!

    def OnFileClose(self, event):
        '''Clears the data tree in response to the
        File/New Project menu button. User is given option to save
        the project.
        '''
        if self.dataFrame:
            self.dataFrame.Clear()
            self.dataFrame.SetLabel('GSAS-II data display') 
        dlg = wx.MessageDialog(self, 'Save current project?', ' ', wx.YES | wx.NO | wx.CANCEL)
        try:
            result = dlg.ShowModal()
            if result == wx.ID_OK:
                self.OnFileSaveMenu(event)
            if result != wx.ID_CANCEL:
                self.GSASprojectfile = ''
                self.PatternTree.SetItemText(self.root,'Loaded Data: ')
                self.PatternTree.DeleteChildren(self.root)
                if self.HKL: self.HKL = []
                if self.G2plotNB.plotList:
                    self.G2plotNB.clear()
        finally:
            dlg.Destroy()

    def OnFileSave(self, event):
        '''Save the current project in response to the
        File/Save Project menu button
        '''
        
        if self.GSASprojectfile: 
            self.PatternTree.SetItemText(self.root,'Loaded Data: '+self.GSASprojectfile)
            self.CheckNotebook()
            G2IO.ProjFileSave(self)
        else:
            self.OnFileSaveas(event)

    def OnFileSaveas(self, event):
        '''Save the current project in response to the
        File/Save as menu button
        '''
        dlg = wx.FileDialog(self, 'Choose GSAS-II project file name', '.', '', 
            'GSAS-II project file (*.gpx)|*.gpx',wx.FD_SAVE|wx.FD_OVERWRITE_PROMPT|wx.CHANGE_DIR)
        try:
            if dlg.ShowModal() == wx.ID_OK:
                self.GSASprojectfile = dlg.GetPath()
                self.GSASprojectfile = G2IO.FileDlgFixExt(dlg,self.GSASprojectfile)
                self.PatternTree.SetItemText(self.root,'Saving project as'+self.GSASprojectfile)
                self.SetTitle("GSAS-II data tree: "+os.path.split(self.GSASprojectfile)[1])
                self.CheckNotebook()
                G2IO.ProjFileSave(self)
                os.chdir(dlg.GetDirectory())           # to get Mac/Linux to change directory!
        finally:
            dlg.Destroy()

    def ExitMain(self, event):
        '''Called if the main window is closed'''
        if self.undofile:
            os.remove(self.undofile)
        sys.exit()
        
    def OnFileExit(self, event):
        '''Called in response to the File/Quit menu button'''
        if self.G2plotNB:
            self.G2plotNB.Destroy()
        if self.dataFrame:
            self.dataFrame.Clear() 
            self.dataFrame.Destroy()
        self.Close()
        
    def OnExportPeakList(self,event):
        dlg = wx.FileDialog(self, 'Choose output peak list file name', '.', '', 
            '(*.*)|*.*',wx.FD_SAVE|wx.FD_OVERWRITE_PROMPT|wx.CHANGE_DIR)
        try:
            if dlg.ShowModal() == wx.ID_OK:
                self.peaklistfile = dlg.GetPath()
                self.peaklistfile = G2IO.FileDlgFixExt(dlg,self.peaklistfile)
                file = open(self.peaklistfile,'w')                
                item, cookie = self.PatternTree.GetFirstChild(self.root)
                while item:
                    name = self.PatternTree.GetItemText(item)
                    if 'PWDR' in name:
                        item2, cookie2 = self.PatternTree.GetFirstChild(item)
                        wave = 0.0
                        while item2:
                            name2 = self.PatternTree.GetItemText(item2)
                            if name2 == 'Instrument Parameters':
                                Inst = self.PatternTree.GetItemPyData(item2)[0]
                                Type = Inst['Type'][0]
                                if 'T' not in Type:
                                    wave = G2mth.getWave(Inst)
                            elif name2 == 'Peak List':
                                peaks = self.PatternTree.GetItemPyData(item2)['peaks']
                            item2, cookie2 = self.PatternTree.GetNextChild(item, cookie2)                            
                        file.write("#%s \n" % (name+' Peak List'))
                        if wave:
                            file.write('#wavelength = %10.6f\n'%(wave))
                        if 'T' in Type:
                            file.write('#%9s %10s %12s %10s %10s %10s %10s %10s\n'%('pos','dsp','int','alp','bet','sig','gam','FWHM'))                                    
                        else:
                            file.write('#%9s %10s %12s %10s %10s %10s\n'%('pos','dsp','int','sig','gam','FWHM'))
                        for peak in peaks:
                            dsp = G2lat.Pos2dsp(Inst,peak[0])
                            if 'T' in Type:  #TOF - more cols
                                FWHM = 2.*G2pwd.getgamFW(peak[10],peak[8])      #to get delta-TOF from Gam(peak)
                                file.write("%10.2f %10.5f %12.2f %10.3f %10.3f %10.3f %10.3f %10.3f\n" % \
                                    (peak[0],dsp,peak[2],np.sqrt(max(0.0001,peak[4])),peak[6],peak[8],peak[10],FWHM))
                            else:               #CW
                                FWHM = 2.*G2pwd.getgamFW(peak[6],peak[4])      #to get delta-2-theta in deg. from Gam(peak)
                                file.write("%10.3f %10.5f %12.2f %10.5f %10.5f %10.5f \n" % \
                                    (peak[0],dsp,peak[2],np.sqrt(max(0.0001,peak[4]))/100.,peak[6]/100.,FWHM/100.)) #convert to deg
                    item, cookie = self.PatternTree.GetNextChild(self.root, cookie)                            
                file.close()
        finally:
            dlg.Destroy()
        
    def OnExportHKL(self,event):
        dlg = wx.FileDialog(self, 'Choose output reflection list file name', '.', '', 
            '(*.*)|*.*',wx.FD_SAVE|wx.FD_OVERWRITE_PROMPT|wx.CHANGE_DIR)
        try:
            if dlg.ShowModal() == wx.ID_OK:
                self.peaklistfile = dlg.GetPath()
                self.peaklistfile = G2IO.FileDlgFixExt(dlg,self.peaklistfile)
                file = open(self.peaklistfile,'w')                
                item, cookie = self.PatternTree.GetFirstChild(self.root)
                while item:
                    name = self.PatternTree.GetItemText(item)
                    if 'PWDR' in name:
                        item2, cookie2 = self.PatternTree.GetFirstChild(item)
                        while item2:
                            name2 = self.PatternTree.GetItemText(item2)
                            if name2 == 'Reflection Lists':
                                data = self.PatternTree.GetItemPyData(item2)
                                phases = data.keys()
                                for phase in phases:
                                    peaks = data[phase]
                                    file.write("%s %s %s \n" % (name,phase,' Reflection List'))
                                    if 'T' in peaks.get('Type','PXC'):
                                        file.write('%s \n'%('   h   k   l   m    d-space     TOF         wid        F**2'))
                                    else:                
                                        file.write('%s \n'%('   h   k   l   m    d-space   2-theta       wid        F**2'))
                                    for peak in peaks['RefList']:
                                        FWHM = 2.*G2pwd.getgamFW(peak[7],peak[6])
                                        if 'T' in peaks.get('Type','PXC'):
                                            file.write(" %3d %3d %3d %3d %10.5f %10.2f %10.5f %10.3f \n" % \
                                                (int(peak[0]),int(peak[1]),int(peak[2]),int(peak[3]),peak[4],peak[5],FWHM,peak[8]))
                                        else:
                                            file.write(" %3d %3d %3d %3d %10.5f %10.5f %10.5f %10.3f \n" % \
                                                (int(peak[0]),int(peak[1]),int(peak[2]),int(peak[3]),peak[4],peak[5],FWHM/100.,peak[8]))
                            item2, cookie2 = self.PatternTree.GetNextChild(item, cookie2)                            
                    item, cookie = self.PatternTree.GetNextChild(self.root, cookie)                            
                file.close()
        finally:
            dlg.Destroy()
        
    def OnExportPDF(self,event):
        #need S(Q) and G(R) to be saved here - probably best from selection?
        names = ['All']
        exports = []
        item, cookie = self.PatternTree.GetFirstChild(self.root)
        while item:
            name = self.PatternTree.GetItemText(item)
            if 'PDF' in name:
                names.append(name)
            item, cookie = self.PatternTree.GetNextChild(self.root, cookie)
        if names:
            dlg = wx.MultiChoiceDialog(self,'Select','PDF patterns to export',names)
            if dlg.ShowModal() == wx.ID_OK:
                sel = dlg.GetSelections()
                if sel[0] == 0:
                    exports = names[1:]
                else:
                    for x in sel:
                        exports.append(names[x])
            dlg.Destroy()
        if exports:
            G2IO.PDFSave(self,exports)
        
    def OnMakePDFs(self,event):
        '''Calculates PDFs
        '''
        sind = lambda x: math.sin(x*math.pi/180.)
        tth2q = lambda t,w:4.0*math.pi*sind(t/2.0)/w
        TextList = ['All PWDR']
        PDFlist = []
        Names = []
        if self.PatternTree.GetCount():
            id, cookie = self.PatternTree.GetFirstChild(self.root)
            while id:
                name = self.PatternTree.GetItemText(id)
                Names.append(name)
                if 'PWDR' in name:
                    TextList.append(name)
                id, cookie = self.PatternTree.GetNextChild(self.root, cookie)
            if len(TextList) == 1:
                self.ErrorDialog('Nothing to make PDFs for','There must be at least one "PWDR" pattern')
                return
            dlg = wx.MultiChoiceDialog(self,'Make PDF controls','Make PDF controls for:',TextList, wx.CHOICEDLG_STYLE)
            try:
                if dlg.ShowModal() == wx.ID_OK:
                    result = dlg.GetSelections()
                    for i in result: PDFlist.append(TextList[i])
                    if 0 in result:
                        PDFlist = [item for item in TextList if item[:4] == 'PWDR']                        
                    for item in PDFlist:
                        PWDRname = item[4:]
                        Id = self.PatternTree.AppendItem(parent=self.root,text='PDF '+PWDRname)
                        Data = {
                            'Sample':{'Name':item,'Mult':1.0,'Add':0.0},
                            'Sample Bkg.':{'Name':'','Mult':-1.0,'Add':0.0},
                            'Container':{'Name':'','Mult':-1.0,'Add':0.0},
                            'Container Bkg.':{'Name':'','Mult':-1.0,'Add':0.0},'ElList':{},
                            'Geometry':'Cylinder','Diam':1.0,'Pack':0.50,'Form Vol':10.0,
                            'DetType':'Image plate','ObliqCoeff':0.2,'Ruland':0.025,'QScaleLim':[0,100],
                            'Lorch':True,}
                        self.PatternTree.SetItemPyData(self.PatternTree.AppendItem(Id,text='PDF Controls'),Data)
                        self.PatternTree.SetItemPyData(self.PatternTree.AppendItem(Id,text='I(Q)'+PWDRname),[])        
                        self.PatternTree.SetItemPyData(self.PatternTree.AppendItem(Id,text='S(Q)'+PWDRname),[])        
                        self.PatternTree.SetItemPyData(self.PatternTree.AppendItem(Id,text='F(Q)'+PWDRname),[])        
                        self.PatternTree.SetItemPyData(self.PatternTree.AppendItem(Id,text='G(R)'+PWDRname),[])        
                for item in self.ExportPDF: item.Enable(True)
            finally:
                dlg.Destroy()
                
    def GetPWDRdatafromTree(self,PWDRname):
        ''' Returns powder data from GSASII tree

        :param str PWDRname: a powder histogram name as obtained from
          :meth:`GSASIIstruct.GetHistogramNames`

        :returns: PWDRdata = powder data dictionary with
          Powder data arrays, Limits, Instrument Parameters,
          Sample Parameters            
        '''
        PWDRdata = {}
        try:
            PWDRdata.update(self.PatternTree.GetItemPyData(PWDRname)[0])            #wtFactor + ?
        except ValueError:
            PWDRdata['wtFactor'] = 1.0
        PWDRdata['Data'] = self.PatternTree.GetItemPyData(PWDRname)[1]          #powder data arrays
        PWDRdata['Limits'] = self.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(self,PWDRname,'Limits'))
        PWDRdata['Background'] = self.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(self,PWDRname,'Background'))
        PWDRdata['Instrument Parameters'] = self.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(self,PWDRname,'Instrument Parameters'))
        PWDRdata['Sample Parameters'] = self.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(self,PWDRname,'Sample Parameters'))
        PWDRdata['Reflection Lists'] = self.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(self,PWDRname,'Reflection Lists'))
        if 'ranId' not in PWDRdata:  # patch, add a random Id
            PWDRdata['ranId'] = ran.randint(0,sys.maxint)
        if 'ranId' not in PWDRdata['Sample Parameters']:  # I hope this becomes obsolete at some point
            PWDRdata['Sample Parameters']['ranId'] = PWDRdata['ranId']
        return PWDRdata

    def GetHKLFdatafromTree(self,HKLFname):
        ''' Returns single crystal data from GSASII tree

        :param str HKLFname: a single crystal histogram name as obtained
          from
          :meth:`GSASIIstruct.GetHistogramNames`

        :returns: HKLFdata = single crystal data list of reflections

        '''
        HKLFdata = {}
        HKLFdata.update(self.PatternTree.GetItemPyData(HKLFname)[0])            #wtFactor + ?
#        try:
#            HKLFdata.update(self.PatternTree.GetItemPyData(HKLFname)[0])            #wtFactor + ?
#        except ValueError:
#            HKLFdata['wtFactor'] = 1.0
        HKLFdata['Data'] = self.PatternTree.GetItemPyData(HKLFname)[1]
        HKLFdata['Instrument Parameters'] = self.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(self,HKLFname,'Instrument Parameters'))
        return HKLFdata
        
    def GetPhaseData(self):
        '''Returns a dict with defined phases. 
        Note routine :func:`GSASIIstrIO.GetPhaseData` also exists to
        get same info from GPX file.
        '''
        phaseData = {}
        if G2gd.GetPatternTreeItemId(self,self.root,'Phases'):
            sub = G2gd.GetPatternTreeItemId(self,self.root,'Phases')
        else:
            print 'no phases found in GetPhaseData'
            sub = None
        if sub:
            item, cookie = self.PatternTree.GetFirstChild(sub)
            while item:
                phaseName = self.PatternTree.GetItemText(item)
                phaseData[phaseName] =  self.PatternTree.GetItemPyData(item)
                if 'ranId' not in phaseData[phaseName]:
                    phaseData[phaseName]['ranId'] = ran.randint(0,sys.maxint)          
                item, cookie = self.PatternTree.GetNextChild(sub, cookie)
        return phaseData

    def GetPhaseInfofromTree(self):
        '''Get the phase names and their rId values,
        also the histograms used in each phase. 

        :returns: (phaseRIdList, usedHistograms) where 

          * phaseRIdList is a list of random Id values for each phase
          * usedHistograms is a dict where the keys are the phase names
            and the values for each key are a list of the histogram names
            used in each phase. 
        '''
        phaseRIdList = []
        usedHistograms = {}
        sub = G2gd.GetPatternTreeItemId(self,self.root,'Phases')
        if sub:
            item, cookie = self.PatternTree.GetFirstChild(sub)
            while item:
                phaseName = self.PatternTree.GetItemText(item)
                ranId = self.PatternTree.GetItemPyData(item).get('ranId')
                if ranId: phaseRIdList.append(ranId)
                data = self.PatternTree.GetItemPyData(item)
                UseList = data['Histograms']
                usedHistograms[phaseName] = UseList.keys()
                item, cookie = self.PatternTree.GetNextChild(sub, cookie)
        return phaseRIdList,usedHistograms

    def GetPhaseNames(self):
        '''Returns a list of defined phases. 
        Note routine :func:`GSASIIstrIO.GetPhaseNames` also exists to
        get same info from GPX file.
        '''
        phaseNames = []
        if G2gd.GetPatternTreeItemId(self,self.root,'Phases'):
            sub = G2gd.GetPatternTreeItemId(self,self.root,'Phases')
        else:
            print 'no phases found in GetPhaseNames'
            sub = None
        if sub:
            item, cookie = self.PatternTree.GetFirstChild(sub)
            while item:
                phase = self.PatternTree.GetItemText(item)
                phaseNames.append(phase)
                item, cookie = self.PatternTree.GetNextChild(sub, cookie)
        return phaseNames
    
    def GetHistogramNames(self,hType):
        """ Returns a list of histogram names found in the GSASII data tree
        Note routine :func:`GSASIIstrIO.GetHistogramNames` also exists to
        get same info from GPX file.
        
        :param str hType: list of histogram types
        :return: list of histogram names
        
        """
        HistogramNames = []
        if self.PatternTree.GetCount():
            item, cookie = self.PatternTree.GetFirstChild(self.root)
            while item:
                name = self.PatternTree.GetItemText(item)
                if name[:4] in hType:
                    HistogramNames.append(name)        
                item, cookie = self.PatternTree.GetNextChild(self.root, cookie)                

        return HistogramNames
                    
    def GetUsedHistogramsAndPhasesfromTree(self):
        ''' Returns all histograms that are found in any phase
        and any phase that uses a histogram.
        This also assigns numbers to used phases and histograms by the
        order they appear in the file. 
        Note routine :func:`GSASIIstrIO.GetUsedHistogramsAndPhasesfromTree` also exists to
        get same info from GPX file.

        :returns: (Histograms,Phases)

            * Histograms = dictionary of histograms as {name:data,...}
            * Phases = dictionary of phases that use histograms
        '''
        Histograms = {}
        Phases = {}
        phaseNames = self.GetPhaseNames()
        phaseData = self.GetPhaseData()
        histoList = self.GetHistogramNames(['PWDR','HKLF'])

        for phase in phaseData:
            Phase = phaseData[phase]
            pId = phaseNames.index(phase)
            Phase['pId'] = pId
            if Phase['Histograms']:
                if phase not in Phases:
                    Phases[phase] = Phase
                for hist in Phase['Histograms']:
                    if 'Use' not in Phase['Histograms'][hist]:      #patch: add Use flag as True
                        Phase['Histograms'][hist]['Use'] = True         
                    if hist not in Histograms and Phase['Histograms'][hist]['Use']:
                        item = G2gd.GetPatternTreeItemId(self,self.root,hist)
                        if item:
                            if 'PWDR' in hist[:4]: 
                                Histograms[hist] = self.GetPWDRdatafromTree(item)
                            elif 'HKLF' in hist[:4]:
                                Histograms[hist] = self.GetHKLFdatafromTree(item)
                            hId = histoList.index(hist)
                            Histograms[hist]['hId'] = hId
                        else: # would happen if a referenced histogram were renamed or deleted
                            print('For phase "'+str(phase)+
                                  '" unresolved reference to histogram "'+str(hist)+'"')
        #G2obj.IndexAllIds(Histograms=Histograms,Phases=Phases)
        G2obj.IndexAllIds(Histograms=Histograms,Phases=phaseData)
        return Histograms,Phases
        
    def MakeLSParmDict(self):
        '''Load all parameters used for computation from the tree into a
        dict of paired values [value, refine flag]. Note that this is
        different than the parmDict used in the refinement, which only has
        values.

        Note that similar things are done in
        :meth:`GSASIIIO.ExportBaseclass.loadParmDict` (from the tree) and 
        :func:`GSASIIstrMain.Refine` and :func:`GSASIIstrMain.SeqRefine` (from
        a GPX file).

        :returns: (parmDict,varyList) where:

         * parmDict is a dict with values and refinement flags
           for each parameter and
         * varyList is a list of variables (refined parameters).
        '''
        parmDict = {}
        Histograms,Phases = self.GetUsedHistogramsAndPhasesfromTree()
        for phase in Phases:
            if 'pId' not in Phases[phase]:
                self.ErrorDialog('View parameter error','You must run least squares at least once')
                raise Exception,'No pId for phase '+str(phase)
        rigidbodyDict = self.PatternTree.GetItemPyData(   
            G2gd.GetPatternTreeItemId(self,self.root,'Rigid bodies'))
        rbVary,rbDict = G2stIO.GetRigidBodyModels(rigidbodyDict,Print=False)
        rbIds = rigidbodyDict.get('RBIds',{'Vector':[],'Residue':[]})
        Natoms,atomIndx,phaseVary,phaseDict,pawleyLookup,FFtable,BLtable,maxSSwave = G2stIO.GetPhaseData(Phases,RestraintDict=None,rbIds=rbIds,Print=False)        
        hapVary,hapDict,controlDict = G2stIO.GetHistogramPhaseData(Phases,Histograms,Print=False)
        histVary,histDict,controlDict = G2stIO.GetHistogramData(Histograms,Print=False)
        varyList = rbVary+phaseVary+hapVary+histVary
        parmDict.update(rbDict)
        parmDict.update(phaseDict)
        parmDict.update(hapDict)
        parmDict.update(histDict)
        for parm in parmDict:
            if parm.split(':')[-1] in ['Azimuth','Gonio. radius','Lam1','Lam2',
                'Omega','Chi','Phi','nDebye','nPeaks']:
                parmDict[parm] = [parmDict[parm],'-']
            elif parm.split(':')[-2] in ['Ax','Ay','Az','SHmodel','SHord']:
                parmDict[parm] = [parmDict[parm],'-']
            elif parm in varyList:
                parmDict[parm] = [parmDict[parm],'T']
            else:
                parmDict[parm] = [parmDict[parm],'F']
        # for i in parmDict: print i,'\t',parmDict[i]
        # fl = open('parmDict.dat','wb')
        # import cPickle
        # cPickle.dump(parmDict,fl,1)
        # fl.close()
        return parmDict,varyList

    def ShowLSParms(self,event):
        '''Displays a window showing all parameters in the refinement.
        Called from the Calculate/View LS Parms menu.
        '''
        parmDict,varyList = self.MakeLSParmDict()
        parmValDict = {}
        for i in parmDict:
            parmValDict[i] = parmDict[i][0]
            
        reqVaryList = tuple(varyList) # save requested variables
        try:
            # process constraints
            sub = G2gd.GetPatternTreeItemId(self,self.root,'Constraints') 
            Constraints = self.PatternTree.GetItemPyData(sub)
            constList = []
            for item in Constraints:
                if item.startswith('_'): continue
                constList += Constraints[item]
            G2mv.InitVars()
            constrDict,fixedList,ignored = G2stIO.ProcessConstraints(constList)
            groups,parmlist = G2mv.GroupConstraints(constrDict)
            G2mv.GenerateConstraints(groups,parmlist,varyList,constrDict,fixedList,parmValDict)
            G2mv.Map2Dict(parmValDict,varyList)
        except:
            pass
        dlg = G2gd.ShowLSParms(self,'Least Squares Parameters',parmValDict,varyList,reqVaryList)
        dlg.ShowModal()
        dlg.Destroy()
        
    def OnRefine(self,event):
        '''Perform a refinement.
        Called from the Calculate/Refine menu.
        '''        
        Id = G2gd.GetPatternTreeItemId(self,self.root,'Sequential results')
        if Id:
            dlg = wx.MessageDialog(
                self,
                'Your last refinement was sequential. Continue with "Refine", removing previous sequential results?',
                'Remove sequential results?',wx.OK|wx.CANCEL)
            if dlg.ShowModal() == wx.ID_OK:
                self.PatternTree.Delete(Id)
                dlg.Destroy()
            else:
                dlg.Destroy()
                return
        self.OnFileSave(event)
        # check that constraints are OK here
        errmsg, warnmsg = G2stIO.ReadCheckConstraints(self.GSASprojectfile)
        if errmsg:
            self.ErrorDialog('Refinement error',errmsg)
            return
        if warnmsg:
            print('Conflict between refinment flag settings and constraints:\n'+
                warnmsg+'\nRefinement not possible')
            self.ErrorDialog('Refinement Flag Error',
                'Conflict between refinement flag settings and constraints:\n'+
                warnmsg+'\nRefinement not possible')
            return
        dlg = wx.ProgressDialog('Residual','All data Rw =',101.0, 
            style = wx.PD_ELAPSED_TIME|wx.PD_AUTO_HIDE|wx.PD_CAN_ABORT,
            parent=self)
        Size = dlg.GetSize()
        if 50 < Size[0] < 500: # sanity check on size, since this fails w/Win & wx3.0
            dlg.SetSize((int(Size[0]*1.2),Size[1])) # increase size a bit along x
        dlg.CenterOnParent()
        Rw = 100.00
        oldId =  self.PatternTree.GetSelection()        #retain current selection
        oldPath = self.GetTreeItemsList(oldId)
        parentName = ''
        oldName = self.PatternTree.GetItemText(oldId)
        parentId = self.PatternTree.GetItemParent(oldId)
        if parentId:
            parentName = self.PatternTree.GetItemText(parentId)     #find the current data tree name
            if 'Phases' in parentName:
                tabId = self.dataDisplay.GetSelection()
        try:
            OK,Msg = G2stMn.Refine(self.GSASprojectfile,dlg)
        finally:
            dlg.Update(101.) # forces the Auto_Hide; needed after move w/Win & wx3.0
            dlg.Destroy()
            wx.Yield()
        if OK:
            Rw = Msg
            dlg2 = wx.MessageDialog(self,'Load new result?','Refinement results, Rw =%.3f'%(Rw),wx.OK|wx.CANCEL)
            try:
                if dlg2.ShowModal() == wx.ID_OK:
                    Id = 0
                    self.PatternTree.DeleteChildren(self.root)
                    self.HKL = []
                    G2IO.ProjFileOpen(self)
                    Id =  self.root
                    txt = None
                    for txt in oldPath:
                        Id = G2gd.GetPatternTreeItemId(self, Id, txt)
                    self.PickIdText = None  #force reload of page
                    self.PickId = Id
                    self.PatternTree.SelectItem(Id)
                    G2gd.MovePatternTreeToGrid(self,Id)
            finally:
                dlg2.Destroy()
        else:
            self.ErrorDialog('Refinement error',Msg)

    def OnSeqRefine(self,event):
        '''Perform a sequential refinement.
        Called from the Calculate/Sequential refine menu.
        '''        
        Id = G2gd.GetPatternTreeItemId(self,self.root,'Sequential results')
        if not Id:
            Id = self.PatternTree.AppendItem(self.root,text='Sequential results')
            self.PatternTree.SetItemPyData(Id,{})            
        Controls = self.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(self,self.root, 'Controls'))
        Controls['ShowCell'] = True
        self.OnFileSave(event)
        # check that constraints are OK here
        errmsg, warnmsg = G2stIO.ReadCheckConstraints(self.GSASprojectfile)
        if errmsg:
            self.ErrorDialog('Refinement error',errmsg)
            return
        if warnmsg:
            print('Conflict between refinment flag settings and constraints:\n'+
                  warnmsg+'\nRefinement not possible')
            self.ErrorDialog('Refinement Flag Error',
                             'Conflict between refinment flag settings and constraints:\n'+
                             warnmsg+'\nRefinement not possible')
            return
        dlg = wx.ProgressDialog('Residual for histogram 0','Powder profile Rwp =',101.0, 
            style = wx.PD_ELAPSED_TIME|wx.PD_AUTO_HIDE|wx.PD_CAN_ABORT,
            parent=self)            
        Size = dlg.GetSize()
        if 50 < Size[0] < 500: # sanity check on size, since this fails w/Win & wx3.0
            dlg.SetSize((int(Size[0]*1.2),Size[1])) # increase size a bit along x
        dlg.CenterOnParent()
        try:
            OK,Msg = G2stMn.SeqRefine(self.GSASprojectfile,dlg)
        finally:
            dlg.Update(101.) # forces the Auto_Hide; needed after move w/Win & wx3.0
            dlg.Destroy()
            wx.Yield()
        if OK:
            dlg = wx.MessageDialog(self,'Load new result?','Refinement results',wx.OK|wx.CANCEL)
            try:
                if dlg.ShowModal() == wx.ID_OK:
                    Id = 0
                    self.PickIdText = None  #force reload of PickId contents
                    self.PatternTree.DeleteChildren(self.root)
                    if len(self.HKL): self.HKL = []
                    if self.G2plotNB.plotList:
                        self.G2plotNB.clear()
                    G2IO.ProjFileOpen(self)
                    Id = G2gd.GetPatternTreeItemId(self,self.root,'Sequential results')
                    self.PatternTree.SelectItem(Id)
    
            finally:
                dlg.Destroy()
        else:
            self.ErrorDialog('Sequential refinement error',Msg)
        
    def ErrorDialog(self,title,message,parent=None, wtype=wx.OK):
        'Display an error message'
        result = None
        if parent is None:
            dlg = wx.MessageDialog(self, message, title,  wtype)
        else:
            dlg = wx.MessageDialog(parent, message, title,  wtype)
            dlg.CenterOnParent() # not working on Mac
        try:
            result = dlg.ShowModal()
        finally:
            dlg.Destroy()
        return result
    
    def OnSaveMultipleImg(self,event):
        '''Select and save multiple image parameter and mask files
        '''
        G2IO.SaveMultipleImg(self)
        
class GSASIImain(wx.App):
    '''Defines a wxApp for GSAS-II

    Creates a wx frame (self.main) which contains the display of the
    data tree.
    '''
    def OnInit(self):
        '''Called automatically when the app is created.'''
        import platform
        if '2.7' not in sys.version[:5]:
            dlg = wx.MessageDialog(None, 
                'GSAS-II requires Python 2.7.x\n Yours is '+sys.version.split()[0],
                'Python version error',  wx.OK)
            try:
                result = dlg.ShowModal()
            finally:
                dlg.Destroy()
            sys.exit()
        self.main = GSASII(None)
        self.main.Show()
        self.SetTopWindow(self.main)
        # save the current package versions
        self.main.PackageVersions = []
        self.main.PackageVersions.append(['Python',sys.version.split()[0]])
        for p in (wx,mpl,np,sp,ogl):
            self.main.PackageVersions.append([p.__name__,p.__version__])
        try:
            self.main.PackageVersions.append([Image.__name__,Image.VERSION])
        except:
            try:
                from PIL import PILLOW_VERSION
                self.main.PackageVersions.append([Image.__name__,PILLOW_VERSION])
            except:
                pass
        self.main.PackageVersions.append([' Platform',sys.platform+' '+platform.architecture()[0]+' '+platform.machine()])
        
#        self.main.PackageVersions = {}
#        self.main.PackageVersions['Python'] = sys.version.split()[0]
#        for p in (wx,mpl,np,sp,ogl):
#            self.main.PackageVersions[p.__name__] = p.__version__
#        try:
#            self.main.PackageVersions[Image.__name__] = Image.VERSION
#        except:
#            try:
#                from PIL import PILLOW_VERSION
#                self.main.PackageVersions[Image.__name__] = PILLOW_VERSION
#            except:
#                pass
#        self.main.PackageVersions[' Platform'] = sys.platform+' '+platform.architecture()[0]+' '+platform.machine()
        # DEBUG: jump to sequential results
        #Id = G2gd.GetPatternTreeItemId(self.main,self.main.root,'Sequential results')
        #self.main.PatternTree.SelectItem(Id)
        # end DEBUG
        return True
    # def MacOpenFile(self, filename):
    #     '''Called on Mac every time a file is dropped on the app when it is running,
    #     treat this like a File/Open project menu action.
    #     Should be ignored on other platforms
    #     '''
    #     # PATCH: Canopy 1.4 script main seems dropped on app; ignore .py files
    #     print 'MacOpen',filename
    #     if os.path.splitext(filename)[1] == '.py': return
    #     # end PATCH
    #     self.main.OnFileOpen(None,filename)
    # removed because this gets triggered when a file is on the command line in canopy 1.4 -- not likely used anyway
       
def main():
    '''Start up the GSAS-II application'''
    #application = GSASIImain() # don't redirect output, someday we
    # may want to do this if we can 
    application = GSASIImain(0)
    if GSASIIpath.GetConfigValue('wxInspector'):
        import wx.lib.inspection as wxeye
        wxeye.InspectionTool().Show()

    #application.main.OnRefine(None)
    application.MainLoop()
    
if __name__ == '__main__':
    # print versions
    print "Python module versions loaded:"
    print "python:     ",sys.version.split()[0]
    print "wxpython:   ",wx.__version__
    print "matplotlib: ",mpl.__version__
    print "numpy:      ",np.__version__
    print "scipy:      ",sp.__version__
    print "OpenGL:     ",ogl.__version__
    try:
        from PIL import Image
        try:
            from PIL import PILLOW_VERSION
            version = PILLOW_VERSION
        except:
            version = Image.VERSION
        print "pillow:     ",version
    except ImportError:
        try:
            import Image
            print "Image (PIL):",Image.VERSION
        except ImportError:
            print "Image module not present; Note that PIL (Python Imaging Library) or pillow is needed for some image operations"
    try:
        import mkl
        print "Max threads ",mkl.get_max_threads()
    except:
        pass
    import platform
    print "Platform info:",sys.platform,platform.architecture()[0],platform.machine()
    #print "wxPython description",wx.PlatformInfo
    print "This is GSAS-II version:     ",__version__,' revision '+str(GSASIIpath.GetVersionNumber())
    GSASIIpath.InvokeDebugOpts()
    main() # start the GUI
