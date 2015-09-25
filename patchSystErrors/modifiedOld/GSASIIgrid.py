# -*- coding: utf-8 -*-
#GSASIIgrid - data display routines
########### SVN repository information ###################
# $Date: 2015-06-04 14:49:58 -0400 (Thu, 04 Jun 2015) $
# $Author: vondreele $
# $Revision: 1878 $
# $URL: https://subversion.xor.aps.anl.gov/pyGSAS/trunk/GSASIIgrid.py $
# $Id: GSASIIgrid.py 1878 2015-06-04 18:49:58Z vondreele $
########### SVN repository information ###################
'''
*GSASIIgrid: Basic GUI routines*
--------------------------------

'''
import wx
import wx.grid as wg
#import wx.wizard as wz
#import wx.aui
import wx.lib.scrolledpanel as wxscroll
# </ Anton Gagin 
import wx.lib.agw.supertooltip as STT
# Anton Gagin /> 
import time
import copy
import cPickle
import sys
import os
import numpy as np
import numpy.ma as ma
import scipy.optimize as so
import GSASIIpath
GSASIIpath.SetVersionNumber("$Revision: 1878 $")
import GSASIImath as G2mth
import GSASIIIO as G2IO
import GSASIIstrIO as G2stIO
import GSASIIlattice as G2lat
import GSASIIplot as G2plt
import GSASIIpwdGUI as G2pdG
import GSASIIimgGUI as G2imG
import GSASIIphsGUI as G2phG
import GSASIIspc as G2spc
import GSASIImapvars as G2mv
import GSASIIconstrGUI as G2cnstG
import GSASIIrestrGUI as G2restG
import GSASIIpy3 as G2py3
import GSASIIobj as G2obj
import GSASIIexprGUI as G2exG
import GSASIIlog as log
import GSASIIctrls as G2G

# trig functions in degrees
sind = lambda x: np.sin(x*np.pi/180.)
tand = lambda x: np.tan(x*np.pi/180.)
cosd = lambda x: np.cos(x*np.pi/180.)

# Define a short name for convenience
WACV = wx.ALIGN_CENTER_VERTICAL

[ wxID_FOURCALC, wxID_FOURSEARCH, wxID_FOURCLEAR, wxID_PEAKSMOVE, wxID_PEAKSCLEAR, 
    wxID_CHARGEFLIP, wxID_PEAKSUNIQUE, wxID_PEAKSDELETE, wxID_PEAKSDA,
    wxID_PEAKSDISTVP, wxID_PEAKSVIEWPT, wxID_FINDEQVPEAKS,wxID_SHOWBONDS,wxID_MULTIMCSA,
    wxID_SINGLEMCSA, wxID_4DMAPCOMPUTE,wxID_4DCHARGEFLIP,
] = [wx.NewId() for item in range(17)]

# </ Anton Gagin 
[ wxID_PWDRADD, wxID_HKLFADD, wxID_PWDANALYSIS, wxID_PWDEXPORTQQ, wxID_PWDCOPY, wxID_PLOTCTRLCOPY, 
    wxID_DATADELETE,wxID_DATACOPY,wxID_DATACOPYFLAGS,wxID_DATASELCOPY,
] = [wx.NewId() for item in range(10)]
# Anton Gagin />  

[ wxID_ATOMSEDITADD, wxID_ATOMSEDITINSERT, wxID_ATOMSEDITDELETE, wxID_ATOMSREFINE, 
    wxID_ATOMSMODIFY, wxID_ATOMSTRANSFORM, wxID_ATOMSVIEWADD, wxID_ATOMVIEWINSERT,
    wxID_RELOADDRAWATOMS,wxID_ATOMSDISAGL,wxID_ATOMMOVE,wxID_MAKEMOLECULE,
    wxID_ASSIGNATMS2RB,wxID_ATOMSPDISAGL, wxID_ISODISP,
] = [wx.NewId() for item in range(15)]

[ wxID_DRAWATOMSTYLE, wxID_DRAWATOMLABEL, wxID_DRAWATOMCOLOR, wxID_DRAWATOMRESETCOLOR, 
    wxID_DRAWVIEWPOINT, wxID_DRAWTRANSFORM, wxID_DRAWDELETE, wxID_DRAWFILLCELL, 
    wxID_DRAWADDEQUIV, wxID_DRAWFILLCOORD, wxID_DRAWDISAGLTOR,  wxID_DRAWPLANE,
    wxID_DRAWDISTVP,
] = [wx.NewId() for item in range(13)]

[ wxID_DRAWRESTRBOND, wxID_DRAWRESTRANGLE, wxID_DRAWRESTRPLANE, wxID_DRAWRESTRCHIRAL,
] = [wx.NewId() for item in range(4)]

[ wxID_ADDMCSAATOM,wxID_ADDMCSARB,wxID_CLEARMCSARB,wxID_MOVEMCSA,wxID_MCSACLEARRESULTS,
] = [wx.NewId() for item in range(5)]

[ wxID_CLEARTEXTURE,wxID_REFINETEXTURE,
] = [wx.NewId() for item in range(2)]

[ wxID_PAWLEYLOAD, wxID_PAWLEYESTIMATE, wxID_PAWLEYUPDATE,
] = [wx.NewId() for item in range(3)]

[ wxID_IMCALIBRATE,wxID_IMRECALIBRATE,wxID_IMINTEGRATE, wxID_IMCLEARCALIB,  
    wxID_IMCOPYCONTROLS, wxID_INTEGRATEALL, wxID_IMSAVECONTROLS, wxID_IMLOADCONTROLS,
] = [wx.NewId() for item in range(8)]

[ wxID_MASKCOPY, wxID_MASKSAVE, wxID_MASKLOAD, wxID_NEWMASKSPOT,wxID_NEWMASKARC,wxID_NEWMASKRING,
    wxID_NEWMASKFRAME, wxID_NEWMASKPOLY,  wxID_MASKLOADNOT,
] = [wx.NewId() for item in range(9)]

[ wxID_STRSTACOPY, wxID_STRSTAFIT, wxID_STRSTASAVE, wxID_STRSTALOAD,wxID_STRSTSAMPLE,
    wxID_APPENDDZERO,wxID_STRSTAALLFIT,wxID_UPDATEDZERO,
] = [wx.NewId() for item in range(8)]

[ wxID_BACKCOPY,wxID_LIMITCOPY, wxID_SAMPLECOPY, wxID_SAMPLECOPYSOME, wxID_BACKFLAGCOPY, wxID_SAMPLEFLAGCOPY,
    wxID_SAMPLESAVE, wxID_SAMPLELOAD,wxID_ADDEXCLREGION,wxID_SETSCALE,wxID_SAMPLE1VAL,wxID_ALLSAMPLELOAD,
    wxID_PEAKSMOVE,
] = [wx.NewId() for item in range(13)]

[ wxID_INSTPRMRESET,wxID_CHANGEWAVETYPE,wxID_INSTCOPY, wxID_INSTFLAGCOPY, wxID_INSTLOAD,
    wxID_INSTSAVE, wxID_INST1VAL, wxID_INSTCALIB,
] = [wx.NewId() for item in range(8)]

[ wxID_UNDO,wxID_LSQPEAKFIT,wxID_LSQONECYCLE,wxID_RESETSIGGAM,wxID_CLEARPEAKS,wxID_AUTOSEARCH,
    wxID_PEAKSCOPY, wxID_SEQPEAKFIT,
] = [wx.NewId() for item in range(8)]

[  wxID_INDXRELOAD, wxID_INDEXPEAKS, wxID_REFINECELL, wxID_COPYCELL, wxID_MAKENEWPHASE,
    wxID_EXPORTCELLS,
] = [wx.NewId() for item in range(6)]

[ wxID_CONSTRAINTADD,wxID_EQUIVADD,wxID_HOLDADD,wxID_FUNCTADD,
  wxID_CONSPHASE, wxID_CONSHIST, wxID_CONSHAP, wxID_CONSGLOBAL,
] = [wx.NewId() for item in range(8)]

[ wxID_RESTRAINTADD, wxID_RESTSELPHASE,wxID_RESTDELETE, wxID_RESRCHANGEVAL, 
    wxID_RESTCHANGEESD,wxID_AARESTRAINTADD,wxID_AARESTRAINTPLOT,
] = [wx.NewId() for item in range(7)]

[ wxID_RIGIDBODYADD,wxID_DRAWDEFINERB,wxID_RIGIDBODYIMPORT,wxID_RESIDUETORSSEQ,
    wxID_AUTOFINDRESRB,wxID_GLOBALRESREFINE,wxID_RBREMOVEALL,wxID_COPYRBPARMS,
    wxID_GLOBALTHERM,wxID_VECTORBODYADD
] = [wx.NewId() for item in range(10)]

[ wxID_RENAMESEQSEL,wxID_SAVESEQSEL,wxID_SAVESEQSELCSV,wxID_SAVESEQCSV,wxID_PLOTSEQSEL,
  wxID_ORGSEQSEL,wxADDSEQVAR,wxDELSEQVAR,wxEDITSEQVAR,wxCOPYPARFIT,wxID_AVESEQSEL,
  wxADDPARFIT,wxDELPARFIT,wxEDITPARFIT,wxDOPARFIT,
] = [wx.NewId() for item in range(15)]

[ wxID_MODELCOPY,wxID_MODELFIT,wxID_MODELADD,wxID_ELEMENTADD,wxID_ELEMENTDELETE,
    wxID_ADDSUBSTANCE,wxID_LOADSUBSTANCE,wxID_DELETESUBSTANCE,wxID_COPYSUBSTANCE,
    wxID_MODELUNDO,wxID_MODELFITALL,wxID_MODELCOPYFLAGS,
] = [wx.NewId() for item in range(12)]

[ wxID_SELECTPHASE,wxID_PWDHKLPLOT,wxID_PWD3DHKLPLOT,wxID_3DALLHKLPLOT,
] = [wx.NewId() for item in range(4)]

[ wxID_PDFCOPYCONTROLS, wxID_PDFSAVECONTROLS, wxID_PDFLOADCONTROLS, 
    wxID_PDFCOMPUTE, wxID_PDFCOMPUTEALL, wxID_PDFADDELEMENT, wxID_PDFDELELEMENT,
] = [wx.NewId() for item in range(7)]

[ wxID_MCRON,wxID_MCRLIST,wxID_MCRSAVE,wxID_MCRPLAY,
] = [wx.NewId() for item in range(4)]

VERY_LIGHT_GREY = wx.Colour(235,235,235)

# Aliases for Classes/Functions moved to GSASIIctrls, all should be tracked down but leaving as a reminder
#SingleFloatDialog = G2G.SingleFloatDialog
#SingleStringDialog = G2G.SingleStringDialog
#MultiStringDialog = G2G.MultiStringDialog
#G2ColumnIDDialog = G2G.G2ColumnIDDialog
#ItemSelector = G2G.ItemSelector
#HorizontalLine = G2G.HorizontalLine
#G2LoggedButton = G2G.G2LoggedButton
#EnumSelector = G2G.EnumSelector
#G2ChoiceButton = G2G.G2ChoiceButton
#GSGrid = G2G.GSGrid
#Table = G2G.Table
#GridFractionEditor = G2G.GridFractionEditor
#GSNoteBook = G2G.GSNoteBook

# Should SGMessageBox, SymOpDialog, DisAglDialog be moved? 

################################################################################
#### GSAS-II class definitions
################################################################################

class SGMessageBox(wx.Dialog):
    ''' Special version of MessageBox that displays space group & super space group text
    in two blocks
    '''
    def __init__(self,parent,title,text,table,):
        wx.Dialog.__init__(self,parent,wx.ID_ANY,title,pos=wx.DefaultPosition,
            style=wx.DEFAULT_DIALOG_STYLE|wx.RESIZE_BORDER)
        self.text=text
        self.table = table
        self.panel = wx.Panel(self)
        mainSizer = wx.BoxSizer(wx.VERTICAL)
        mainSizer.Add((0,10))
        for line in text:
            mainSizer.Add(wx.StaticText(self.panel,label='     %s     '%(line)),0,WACV)
        ncol = self.table[0].count(',')+1
        tableSizer = wx.FlexGridSizer(0,2*ncol+3,0,0)
        for j,item in enumerate(self.table):
            num,flds = item.split(')')
            tableSizer.Add(wx.StaticText(self.panel,label='     %s  '%(num+')')),0,WACV|wx.ALIGN_LEFT)            
            flds = flds.replace(' ','').split(',')
            for i,fld in enumerate(flds):
                if i < ncol-1:
                    tableSizer.Add(wx.StaticText(self.panel,label='%s, '%(fld)),0,WACV|wx.ALIGN_RIGHT)
                else:
                    tableSizer.Add(wx.StaticText(self.panel,label='%s'%(fld)),0,WACV|wx.ALIGN_RIGHT)
            if not j%2:
                tableSizer.Add((20,0))
        mainSizer.Add(tableSizer,0,wx.ALIGN_LEFT)
        btnsizer = wx.StdDialogButtonSizer()
        OKbtn = wx.Button(self.panel, wx.ID_OK)
        OKbtn.SetDefault()
        btnsizer.AddButton(OKbtn)
        btnsizer.Realize()
        mainSizer.Add((0,10))
        mainSizer.Add(btnsizer,0,wx.ALIGN_CENTER)
        self.panel.SetSizer(mainSizer)
        self.panel.Fit()
        self.Fit()
        size = self.GetSize()
        self.SetSize([size[0]+20,size[1]])

    def Show(self):
        '''Use this method after creating the dialog to post it
        '''
        self.ShowModal()
        return

################################################################################
class SymOpDialog(wx.Dialog):
    '''Class to select a symmetry operator
    '''
    def __init__(self,parent,SGData,New=True,ForceUnit=False):
        wx.Dialog.__init__(self,parent,-1,'Select symmetry operator',
            pos=wx.DefaultPosition,style=wx.DEFAULT_DIALOG_STYLE)
        panel = wx.Panel(self)
        self.SGData = SGData
        self.New = New
        self.Force = ForceUnit
        self.OpSelected = [0,0,0,[0,0,0],False,False]
        mainSizer = wx.BoxSizer(wx.VERTICAL)
        if ForceUnit:
            choice = ['No','Yes']
            self.force = wx.RadioBox(panel,-1,'Force to unit cell?',choices=choice)
            self.force.Bind(wx.EVT_RADIOBOX, self.OnOpSelect)
            mainSizer.Add(self.force,0,WACV|wx.TOP,5)
#        if SGData['SGInv']:
        choice = ['No','Yes']
        self.inv = wx.RadioBox(panel,-1,'Choose inversion?',choices=choice)
        self.inv.Bind(wx.EVT_RADIOBOX, self.OnOpSelect)
        mainSizer.Add(self.inv,0,WACV)
        if SGData['SGLatt'] != 'P':
            LattOp = G2spc.Latt2text(SGData['SGLatt']).split(';')
            self.latt = wx.RadioBox(panel,-1,'Choose cell centering?',choices=LattOp)
            self.latt.Bind(wx.EVT_RADIOBOX, self.OnOpSelect)
            mainSizer.Add(self.latt,0,WACV)
        if SGData['SGLaue'] in ['-1','2/m','mmm','4/m','4/mmm']:
            Ncol = 2
        else:
            Ncol = 3
        OpList = []
        for Opr in SGData['SGOps']:
            OpList.append(G2spc.MT2text(Opr))
        self.oprs = wx.RadioBox(panel,-1,'Choose space group operator?',choices=OpList,
            majorDimension=Ncol)
        self.oprs.Bind(wx.EVT_RADIOBOX, self.OnOpSelect)
        mainSizer.Add(self.oprs,0,WACV|wx.BOTTOM,5)
        mainSizer.Add(wx.StaticText(panel,-1,"   Choose unit cell?"),0,WACV)
        cellSizer = wx.BoxSizer(wx.HORIZONTAL)
        cellName = ['X','Y','Z']
        self.cell = []
        for i in range(3):
            self.cell.append(wx.SpinCtrl(panel,-1,cellName[i],size=wx.Size(50,20)))
            self.cell[-1].SetRange(-3,3)
            self.cell[-1].SetValue(0)
            self.cell[-1].Bind(wx.EVT_SPINCTRL, self.OnOpSelect)
            cellSizer.Add(self.cell[-1],0,WACV)
        mainSizer.Add(cellSizer,0,WACV|wx.BOTTOM,5)
        if self.New:
            choice = ['No','Yes']
            self.new = wx.RadioBox(panel,-1,'Generate new positions?',choices=choice)
            self.new.Bind(wx.EVT_RADIOBOX, self.OnOpSelect)
            mainSizer.Add(self.new,0,WACV)

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

    def OnOpSelect(self,event):
#        if self.SGData['SGInv']:
        self.OpSelected[0] = self.inv.GetSelection()
        if self.SGData['SGLatt'] != 'P':
            self.OpSelected[1] = self.latt.GetSelection()
        self.OpSelected[2] = self.oprs.GetSelection()
        for i in range(3):
            self.OpSelected[3][i] = float(self.cell[i].GetValue())
        if self.New:
            self.OpSelected[4] = self.new.GetSelection()
        if self.Force:
            self.OpSelected[5] = self.force.GetSelection()

    def GetSelection(self):
        return self.OpSelected

    def OnOk(self,event):
        parent = self.GetParent()
        parent.Raise()
        self.EndModal(wx.ID_OK)

    def OnCancel(self,event):
        parent = self.GetParent()
        parent.Raise()
        self.EndModal(wx.ID_CANCEL)

class DisAglDialog(wx.Dialog):
    '''Distance/Angle Controls input dialog. After
    :meth:`ShowModal` returns, the results are found in
    dict :attr:`self.data`, which is accessed using :meth:`GetData`.

    :param wx.Frame parent: reference to parent frame (or None)
    :param dict data: a dict containing the current
      search ranges or an empty dict, which causes default values
      to be used.
      Will be used to set element `DisAglCtls` in 
      :ref:`Phase Tree Item <Phase_table>`
    :param dict default:  A dict containing the default
      search ranges for each element.
    '''
    def __init__(self,parent,data,default):
        wx.Dialog.__init__(self,parent,wx.ID_ANY,
                           'Distance Angle Controls', 
            pos=wx.DefaultPosition,style=wx.DEFAULT_DIALOG_STYLE)
        self.default = default
        self.panel = wx.Panel(self)         #just a dummy - gets destroyed in Draw!
        self._default(data,self.default)
        self.Draw(self.data)
                
    def _default(self,data,default):
        '''Set starting values for the search values, either from
        the input array or from defaults, if input is null
        '''
        if data:
            self.data = copy.deepcopy(data) # don't mess with originals
        else:
            self.data = {}
            self.data['Name'] = default['Name']
            self.data['Factors'] = [0.85,0.85]
            self.data['AtomTypes'] = default['AtomTypes']
            self.data['BondRadii'] = default['BondRadii'][:]
            self.data['AngleRadii'] = default['AngleRadii'][:]

    def Draw(self,data):
        '''Creates the contents of the dialog. Normally called
        by :meth:`__init__`.
        '''
        self.panel.Destroy()
        self.panel = wx.Panel(self)
        mainSizer = wx.BoxSizer(wx.VERTICAL)
        mainSizer.Add(wx.StaticText(self.panel,-1,'Controls for phase '+data['Name']),
            0,WACV|wx.LEFT,10)
        mainSizer.Add((10,10),1)
        
        radiiSizer = wx.FlexGridSizer(0,3,5,5)
        radiiSizer.Add(wx.StaticText(self.panel,-1,' Type'),0,WACV)
        radiiSizer.Add(wx.StaticText(self.panel,-1,'Bond radii'),0,WACV)
        radiiSizer.Add(wx.StaticText(self.panel,-1,'Angle radii'),0,WACV)
        self.objList = {}
        for id,item in enumerate(self.data['AtomTypes']):
            radiiSizer.Add(wx.StaticText(self.panel,-1,' '+item),0,WACV)
            bRadii = wx.TextCtrl(self.panel,-1,value='%.3f'%(data['BondRadii'][id]),style=wx.TE_PROCESS_ENTER)
            self.objList[bRadii.GetId()] = ['BondRadii',id]
            bRadii.Bind(wx.EVT_TEXT_ENTER,self.OnRadiiVal)
            bRadii.Bind(wx.EVT_KILL_FOCUS,self.OnRadiiVal)
            radiiSizer.Add(bRadii,0,WACV)
            aRadii = wx.TextCtrl(self.panel,-1,value='%.3f'%(data['AngleRadii'][id]),style=wx.TE_PROCESS_ENTER)
            self.objList[aRadii.GetId()] = ['AngleRadii',id]
            aRadii.Bind(wx.EVT_TEXT_ENTER,self.OnRadiiVal)
            aRadii.Bind(wx.EVT_KILL_FOCUS,self.OnRadiiVal)
            radiiSizer.Add(aRadii,0,WACV)
        mainSizer.Add(radiiSizer,0,wx.EXPAND)
        factorSizer = wx.FlexGridSizer(0,2,5,5)
        Names = ['Bond','Angle']
        for i,name in enumerate(Names):
            factorSizer.Add(wx.StaticText(self.panel,-1,name+' search factor'),0,WACV)
            bondFact = wx.TextCtrl(self.panel,-1,value='%.3f'%(data['Factors'][i]),style=wx.TE_PROCESS_ENTER)
            self.objList[bondFact.GetId()] = ['Factors',i]
            bondFact.Bind(wx.EVT_TEXT_ENTER,self.OnRadiiVal)
            bondFact.Bind(wx.EVT_KILL_FOCUS,self.OnRadiiVal)
            factorSizer.Add(bondFact)
        mainSizer.Add(factorSizer,0,wx.EXPAND)
        
        OkBtn = wx.Button(self.panel,-1,"Ok")
        OkBtn.Bind(wx.EVT_BUTTON, self.OnOk)
        ResetBtn = wx.Button(self.panel,-1,'Reset')
        ResetBtn.Bind(wx.EVT_BUTTON, self.OnReset)
        btnSizer = wx.BoxSizer(wx.HORIZONTAL)
        btnSizer.Add((20,20),1)
        btnSizer.Add(OkBtn)
        btnSizer.Add(ResetBtn)
        btnSizer.Add((20,20),1)
        mainSizer.Add(btnSizer,0,wx.EXPAND|wx.BOTTOM|wx.TOP, 10)
        self.panel.SetSizer(mainSizer)
        self.panel.Fit()
        self.Fit()
    
    def OnRadiiVal(self,event):
        Obj = event.GetEventObject()
        item = self.objList[Obj.GetId()]
        try:
            self.data[item[0]][item[1]] = float(Obj.GetValue())
        except ValueError:
            pass
        Obj.SetValue("%.3f"%(self.data[item[0]][item[1]]))          #reset in case of error
        
    def GetData(self):
        'Returns the values from the dialog'
        return self.data
        
    def OnOk(self,event):
        'Called when the OK button is pressed'
        parent = self.GetParent()
        parent.Raise()
        self.EndModal(wx.ID_OK)              
        
    def OnReset(self,event):
        'Called when the Reset button is pressed'
        data = {}
        self._default(data,self.default)
        self.Draw(self.data)
                
################################################################################
class ShowLSParms(wx.Dialog):
    '''Create frame to show least-squares parameters
    '''
    def __init__(self,parent,title,parmDict,varyList,fullVaryList,
                 size=(300,430)):
        wx.Dialog.__init__(self,parent,wx.ID_ANY,title,size=size,
                           style=wx.DEFAULT_DIALOG_STYLE|wx.RESIZE_BORDER)
        mainSizer = wx.BoxSizer(wx.VERTICAL)

        panel = wxscroll.ScrolledPanel(
            self, wx.ID_ANY,
            #size=size,
            style = wx.TAB_TRAVERSAL|wx.SUNKEN_BORDER)
        num = len(varyList)
        mainSizer.Add(wx.StaticText(self,wx.ID_ANY,'Number of refined variables: '+str(num)))
        if len(varyList) != len(fullVaryList):
            num = len(fullVaryList) - len(varyList)
            mainSizer.Add(wx.StaticText(self,wx.ID_ANY,' + '+str(num)+' parameters are varied via constraints'))
        subSizer = wx.FlexGridSizer(cols=4,hgap=2,vgap=2)
        parmNames = parmDict.keys()
        parmNames.sort()
        subSizer.Add((-1,-1))
        subSizer.Add(wx.StaticText(panel,wx.ID_ANY,'Parameter name  '))
        subSizer.Add(wx.StaticText(panel,wx.ID_ANY,'refine?'))
        subSizer.Add(wx.StaticText(panel,wx.ID_ANY,'value'),0,wx.ALIGN_RIGHT)
        explainRefine = False
        for name in parmNames:
            # skip entries without numerical values
            if isinstance(parmDict[name],basestring): continue
            try:
                value = G2py3.FormatSigFigs(parmDict[name])
            except TypeError:
                value = str(parmDict[name])+' -?' # unexpected
                #continue
            v = G2obj.getVarDescr(name)
            if v is None or v[-1] is None:
                subSizer.Add((-1,-1))
            else:                
                ch = G2G.HelpButton(panel,G2obj.fmtVarDescr(name))
                subSizer.Add(ch,0,wx.LEFT|wx.RIGHT|WACV|wx.ALIGN_CENTER,1)
            subSizer.Add(wx.StaticText(panel,wx.ID_ANY,str(name)))
            if name in varyList:
                subSizer.Add(wx.StaticText(panel,wx.ID_ANY,'R'))
            elif name in fullVaryList:
                subSizer.Add(wx.StaticText(panel,wx.ID_ANY,'C'))
                explainRefine = True
            else:
                subSizer.Add((-1,-1))
            subSizer.Add(wx.StaticText(panel,wx.ID_ANY,value),0,wx.ALIGN_RIGHT)

        # finish up ScrolledPanel
        panel.SetSizer(subSizer)
        panel.SetAutoLayout(1)
        panel.SetupScrolling()
        mainSizer.Add(panel,1, wx.ALL|wx.EXPAND,1)

        if explainRefine:
            mainSizer.Add(
                wx.StaticText(self,wx.ID_ANY,
                          '"R" indicates a refined variable\n'+
                          '"C" indicates generated from a constraint'
                          ),
                0, wx.ALL,0)
        # make OK button 
        btnsizer = wx.BoxSizer(wx.HORIZONTAL)
        btn = wx.Button(self, wx.ID_CLOSE,"Close") 
        btn.Bind(wx.EVT_BUTTON,self._onClose)
        btnsizer.Add(btn)
        mainSizer.Add(btnsizer, 0, wx.ALIGN_CENTER|wx.ALL, 5)
        # Allow window to be enlarged but not made smaller
        self.SetSizer(mainSizer)
        self.SetMinSize(self.GetSize())

    def _onClose(self,event):
        self.EndModal(wx.ID_CANCEL)
 
################################################################################
class DataFrame(wx.Frame):
    '''Create the data item window and all the entries in menus used in
    that window. For Linux and windows, the menu entries are created for the
    current data item window, but in the Mac the menu is accessed from all
    windows. This means that a different menu is posted depending on which
    data item is posted. On the Mac, all the menus contain the data tree menu
    items, but additional menus are added specific to the data item. 

    Note that while the menus are created here, 
    the binding for the menus is done later in various GSASII*GUI modules,
    where the functions to be called are defined.
    '''
    def Bind(self,eventtype,handler,*args,**kwargs):
        '''Override the Bind() function: on the Mac the binding is to
        the main window, so that menus operate with any window on top.
        For other platforms, either wrap calls that will be logged
        or call the default wx.Frame Bind() to bind to the menu item directly.

        Note that bindings can be made to objects by Id or by direct reference to the
        object. As a convention, when bindings are to objects, they are not logged
        but when bindings are by Id, they are logged.
        '''
        if sys.platform == "darwin": # mac
            self.G2frame.Bind(eventtype,handler,*args,**kwargs)
            return
        if eventtype == wx.EVT_MENU and 'id' in kwargs:
            menulabels = log.SaveMenuCommand(kwargs['id'],self.G2frame,handler)
            if menulabels:
                #print 'intercepting bind for',handler,menulabels,kwargs['id']
                wx.Frame.Bind(self,eventtype,self.G2frame.MenuBinding,*args,**kwargs)
                return
            wx.Frame.Bind(self,eventtype,handler,*args,**kwargs)      
        
    def PrefillDataMenu(self,menu,helpType,helpLbl=None,empty=False):
        '''Create the "standard" part of data frame menus. Note that on Linux and
        Windows nothing happens here. On Mac, this menu duplicates the
        tree menu, but adds an extra help command for the data item and a separator. 
        '''
        self.datamenu = menu
        self.G2frame.dataMenuBars.append(menu)
        self.helpType = helpType
        self.helpLbl = helpLbl
        if sys.platform == "darwin": # mac                         
            self.G2frame.FillMainMenu(menu) # add the data tree menu items
            if not empty:
                menu.Append(wx.Menu(title=''),title='|') # add a separator
        
    def PostfillDataMenu(self,empty=False):
        '''Create the "standard" part of data frame menus. Note that on Linux and
        Windows, this is the standard help Menu. On Mac, this menu duplicates the
        tree menu, but adds an extra help command for the data item and a separator. 
        '''
        menu = self.datamenu
        helpType = self.helpType
        helpLbl = self.helpLbl
        if sys.platform == "darwin": # mac
            if not empty:
                menu.Append(wx.Menu(title=''),title='|') # add another separator
            menu.Append(G2G.AddHelp(self.G2frame,helpType=helpType, helpLbl=helpLbl),
                        title='&Help')
        else: # other
            menu.Append(menu=G2G.MyHelp(self,helpType=helpType, helpLbl=helpLbl),
                        title='&Help')

    def _init_menus(self):
        'define all GSAS-II data frame menus'

        # for use where no menu or data frame help is provided
        self.BlankMenu = wx.MenuBar()
        
        # Controls
        self.ControlsMenu = wx.MenuBar()
        self.PrefillDataMenu(self.ControlsMenu,helpType='Controls',empty=True)
        self.PostfillDataMenu(empty=True)
        
        # Notebook
        self.DataNotebookMenu = wx.MenuBar() 
        self.PrefillDataMenu(self.DataNotebookMenu,helpType='Notebook',empty=True)
        self.PostfillDataMenu(empty=True)
        
        # Comments
        self.DataCommentsMenu = wx.MenuBar()
        self.PrefillDataMenu(self.DataCommentsMenu,helpType='Comments',empty=True)
        self.PostfillDataMenu(empty=True)
        
        # Constraints - something amiss here - get weird wx C++ error after refine!
        self.ConstraintMenu = wx.MenuBar()
        self.PrefillDataMenu(self.ConstraintMenu,helpType='Constraints')
        self.ConstraintTab = wx.Menu(title='')
        self.ConstraintMenu.Append(menu=self.ConstraintTab, title='Select tab')
        for id,txt in (
            (wxID_CONSPHASE,'Phase'),
            (wxID_CONSHAP,'Histogram/Phase'),
            (wxID_CONSHIST,'Histogram'),
            (wxID_CONSGLOBAL,'Global')):
            self.ConstraintTab.Append(
                id=id, kind=wx.ITEM_NORMAL,text=txt,
                help='Select '+txt+' constraint editing tab')
        self.ConstraintEdit = wx.Menu(title='')
        self.ConstraintMenu.Append(menu=self.ConstraintEdit, title='Edit')
        self.ConstraintEdit.Append(id=wxID_HOLDADD, kind=wx.ITEM_NORMAL,text='Add hold',
            help='Add hold on a parameter value')
        self.ConstraintEdit.Append(id=wxID_EQUIVADD, kind=wx.ITEM_NORMAL,text='Add equivalence',
            help='Add equivalence between parameter values')
        self.ConstraintEdit.Append(id=wxID_CONSTRAINTADD, kind=wx.ITEM_NORMAL,text='Add constraint',
            help='Add constraint on parameter values')
        self.ConstraintEdit.Append(id=wxID_FUNCTADD, kind=wx.ITEM_NORMAL,text='Add New Var',
            help='Add variable composed of existing parameter')
        self.PostfillDataMenu()

        # item = self.ConstraintEdit.Append(id=wx.ID_ANY,kind=wx.ITEM_NORMAL,text='Update GUI')
        # def UpdateGSASIIconstrGUI(event):
        #     import GSASIIconstrGUI
        #     reload(GSASIIconstrGUI)
        #     import GSASIIobj
        #     reload(GSASIIobj)
        # self.Bind(wx.EVT_MENU,UpdateGSASIIconstrGUI,id=item.GetId())

        # Rigid bodies
        self.RigidBodyMenu = wx.MenuBar()
        self.PrefillDataMenu(self.RigidBodyMenu,helpType='Rigid bodies')
        self.ResidueRBMenu = wx.Menu(title='')
        self.ResidueRBMenu.Append(id=wxID_RIGIDBODYIMPORT, kind=wx.ITEM_NORMAL,text='Import XYZ',
            help='Import rigid body XYZ from file')
        self.ResidueRBMenu.Append(id=wxID_RESIDUETORSSEQ, kind=wx.ITEM_NORMAL,text='Define sequence',
            help='Define torsion sequence')
        self.ResidueRBMenu.Append(id=wxID_RIGIDBODYADD, kind=wx.ITEM_NORMAL,text='Import residues',
            help='Import residue rigid bodies from macro file')
        self.RigidBodyMenu.Append(menu=self.ResidueRBMenu, title='Edit Body')
        self.PostfillDataMenu()

        self.VectorBodyMenu = wx.MenuBar()
        self.PrefillDataMenu(self.VectorBodyMenu,helpType='Vector rigid bodies')
        self.VectorRBEdit = wx.Menu(title='')
        self.VectorRBEdit.Append(id=wxID_VECTORBODYADD, kind=wx.ITEM_NORMAL,text='Add rigid body',
            help='Add vector rigid body')
        self.VectorBodyMenu.Append(menu=self.VectorRBEdit, title='Edit Vector Body')
        self.PostfillDataMenu()

                    
        # Restraints
        self.RestraintTab = wx.Menu(title='')
        self.RestraintEdit = wx.Menu(title='')
        self.RestraintEdit.Append(id=wxID_RESTSELPHASE, kind=wx.ITEM_NORMAL,text='Select phase',
            help='Select phase')
        self.RestraintEdit.Append(id=wxID_RESTRAINTADD, kind=wx.ITEM_NORMAL,text='Add restraints',
            help='Add restraints')
        self.RestraintEdit.Enable(wxID_RESTRAINTADD,True)    #gets disabled if macromolecule phase
        self.RestraintEdit.Append(id=wxID_AARESTRAINTADD, kind=wx.ITEM_NORMAL,text='Add residue restraints',
            help='Add residue based restraints for macromolecules from macro file')
        self.RestraintEdit.Enable(wxID_AARESTRAINTADD,False)    #gets enabled if macromolecule phase
        self.RestraintEdit.Append(id=wxID_AARESTRAINTPLOT, kind=wx.ITEM_NORMAL,text='Plot residue restraints',
            help='Plot selected residue based restraints for macromolecules from macro file')
        self.RestraintEdit.Enable(wxID_AARESTRAINTPLOT,False)    #gets enabled if macromolecule phase
        self.RestraintEdit.Append(id=wxID_RESRCHANGEVAL, kind=wx.ITEM_NORMAL,text='Change value',
            help='Change observed value')
        self.RestraintEdit.Append(id=wxID_RESTCHANGEESD, kind=wx.ITEM_NORMAL,text='Change esd',
            help='Change esd in observed value')
        self.RestraintEdit.Append(id=wxID_RESTDELETE, kind=wx.ITEM_NORMAL,text='Delete restraints',
            help='Delete selected restraints')

        self.RestraintMenu = wx.MenuBar()
        self.PrefillDataMenu(self.RestraintMenu,helpType='Restraints')
        self.RestraintMenu.Append(menu=self.RestraintTab, title='Select tab')
        self.RestraintMenu.Append(menu=self.RestraintEdit, title='Edit')
        self.PostfillDataMenu()
            
        # Sequential results
        self.SequentialMenu = wx.MenuBar()
        self.PrefillDataMenu(self.SequentialMenu,helpType='Sequential',helpLbl='Sequential Refinement')
        self.SequentialFile = wx.Menu(title='')
        self.SequentialMenu.Append(menu=self.SequentialFile, title='Columns')
        self.SequentialFile.Append(id=wxID_RENAMESEQSEL, kind=wx.ITEM_NORMAL,text='Rename selected',
            help='Rename selected sequential refinement columns')
        self.SequentialFile.Append(id=wxID_SAVESEQSEL, kind=wx.ITEM_NORMAL,text='Save selected as text',
            help='Save selected sequential refinement results as a text file')
        self.SequentialFile.Append(id=wxID_SAVESEQCSV, kind=wx.ITEM_NORMAL,text='Save all as CSV',
            help='Save all sequential refinement results as a CSV spreadsheet file')
        self.SequentialFile.Append(id=wxID_SAVESEQSELCSV, kind=wx.ITEM_NORMAL,text='Save selected as CSV',
            help='Save selected sequential refinement results as a CSV spreadsheet file')
        self.SequentialFile.Append(id=wxID_PLOTSEQSEL, kind=wx.ITEM_NORMAL,text='Plot selected',
            help='Plot selected sequential refinement results')
        self.SequentialFile.Append(id=wxID_AVESEQSEL, kind=wx.ITEM_NORMAL,text='Compute average',
            help='Compute average for selected parameter')            
        self.SequentialFile.Append(id=wxID_ORGSEQSEL, kind=wx.ITEM_NORMAL,text='Reorganize',
            help='Reorganize variables where variables change')
        self.SequentialPvars = wx.Menu(title='')
        self.SequentialMenu.Append(menu=self.SequentialPvars, title='Pseudo Vars')
        self.SequentialPvars.Append(
            id=wxADDSEQVAR, kind=wx.ITEM_NORMAL,text='Add',
            help='Add a new pseudo-variable')
        self.SequentialPvars.Append(
            id=wxDELSEQVAR, kind=wx.ITEM_NORMAL,text='Delete',
            help='Delete an existing pseudo-variable')
        self.SequentialPvars.Append(
            id=wxEDITSEQVAR, kind=wx.ITEM_NORMAL,text='Edit',
            help='Edit an existing pseudo-variable')

        self.SequentialPfit = wx.Menu(title='')
        self.SequentialMenu.Append(menu=self.SequentialPfit, title='Parametric Fit')
        self.SequentialPfit.Append(
            id=wxADDPARFIT, kind=wx.ITEM_NORMAL,text='Add equation',
            help='Add a new equation to minimize')
        self.SequentialPfit.Append(
            id=wxCOPYPARFIT, kind=wx.ITEM_NORMAL,text='Copy equation',
            help='Copy an equation to minimize - edit it next')
        self.SequentialPfit.Append(
            id=wxDELPARFIT, kind=wx.ITEM_NORMAL,text='Delete equation',
            help='Delete an equation for parametric minimization')
        self.SequentialPfit.Append(
            id=wxEDITPARFIT, kind=wx.ITEM_NORMAL,text='Edit equation',
            help='Edit an existing parametric minimization equation')
        self.SequentialPfit.Append(
            id=wxDOPARFIT, kind=wx.ITEM_NORMAL,text='Fit to equation(s)',
            help='Perform a parametric minimization')
        self.PostfillDataMenu()
            
        # PWDR & SASD
        self.PWDRMenu = wx.MenuBar()
        self.PrefillDataMenu(self.PWDRMenu,helpType='PWDR Analysis',helpLbl='Powder Fit Error Analysis')
        self.ErrorAnal = wx.Menu(title='')
        self.PWDRMenu.Append(menu=self.ErrorAnal,title='Commands')
        self.ErrorAnal.Append(id=wxID_PWDANALYSIS,kind=wx.ITEM_NORMAL,text='Error Analysis',
            help='Error analysis on powder pattern')
# </ Anton Gagin             
        self.ErrorAnal.Append(id=wxID_PWDEXPORTQQ,kind=wx.ITEM_NORMAL,text='Export QQ-plot as text',
            help='Save data for QQ-plot as a text file')    
# Anton Gagin />  
        self.ErrorAnal.Append(id=wxID_PWDCOPY,kind=wx.ITEM_NORMAL,text='Copy params',
            help='Copy of PWDR parameters')
        self.ErrorAnal.Append(id=wxID_PLOTCTRLCOPY,kind=wx.ITEM_NORMAL,text='Copy plot controls',
            help='Copy of PWDR plot controls')
            
        self.PostfillDataMenu()
            
        # HKLF 
        self.HKLFMenu = wx.MenuBar()
        self.PrefillDataMenu(self.HKLFMenu,helpType='HKLF Analysis',helpLbl='HKLF Fit Error Analysis')
        self.ErrorAnal = wx.Menu(title='')
        self.HKLFMenu.Append(menu=self.ErrorAnal,title='Commands')
        self.ErrorAnal.Append(id=wxID_PWDANALYSIS,kind=wx.ITEM_NORMAL,text='Error Analysis',
            help='Error analysis on single crystal data')
        self.ErrorAnal.Append(id=wxID_PWD3DHKLPLOT,kind=wx.ITEM_NORMAL,text='Plot 3D HKLs',
            help='Plot HKLs from single crystal data in 3D')
        self.ErrorAnal.Append(id=wxID_3DALLHKLPLOT,kind=wx.ITEM_NORMAL,text='Plot all 3D HKLs',
            help='Plot HKLs from all single crystal data in 3D')
        self.ErrorAnal.Append(id=wxID_PWDCOPY,kind=wx.ITEM_NORMAL,text='Copy params',
            help='Copy of HKLF parameters')
        self.PostfillDataMenu()
            
        # PDR / Limits
        self.LimitMenu = wx.MenuBar()
        self.PrefillDataMenu(self.LimitMenu,helpType='Limits')
        self.LimitEdit = wx.Menu(title='')
        self.LimitMenu.Append(menu=self.LimitEdit, title='Edit')
        self.LimitEdit.Append(id=wxID_LIMITCOPY, kind=wx.ITEM_NORMAL,text='Copy',
            help='Copy limits to other histograms')
        self.LimitEdit.Append(id=wxID_ADDEXCLREGION, kind=wx.ITEM_NORMAL,text='Add exclude',
            help='Add excluded region - select a point on plot; drag to adjust')            
        self.PostfillDataMenu()
            
        # PDR / Background
        self.BackMenu = wx.MenuBar()
        self.PrefillDataMenu(self.BackMenu,helpType='Background')
        self.BackEdit = wx.Menu(title='')
        self.BackMenu.Append(menu=self.BackEdit, title='File')
        self.BackEdit.Append(id=wxID_BACKCOPY, kind=wx.ITEM_NORMAL,text='Copy',
            help='Copy background parameters to other histograms')
        self.BackEdit.Append(id=wxID_BACKFLAGCOPY, kind=wx.ITEM_NORMAL,text='Copy flags',
            help='Copy background refinement flags to other histograms')
        self.BackEdit.Append(id=wxID_PEAKSMOVE, kind=wx.ITEM_NORMAL,text='Move peaks',
            help='Move background peaks to Peak List')
            
        self.PostfillDataMenu()
            
        # PDR / Instrument Parameters
        self.InstMenu = wx.MenuBar()
        self.PrefillDataMenu(self.InstMenu,helpType='Instrument Parameters')
        self.InstEdit = wx.Menu(title='')
        self.InstMenu.Append(menu=self.InstEdit, title='Operations')
        self.InstEdit.Append(help='Calibrate from indexed peaks', 
            id=wxID_INSTCALIB, kind=wx.ITEM_NORMAL,text='Calibrate')            
        self.InstEdit.Append(help='Reset instrument profile parameters to default', 
            id=wxID_INSTPRMRESET, kind=wx.ITEM_NORMAL,text='Reset profile')            
        self.InstEdit.Append(help='Load instrument profile parameters from file', 
            id=wxID_INSTLOAD, kind=wx.ITEM_NORMAL,text='Load profile...')            
        self.InstEdit.Append(help='Save instrument profile parameters to file', 
            id=wxID_INSTSAVE, kind=wx.ITEM_NORMAL,text='Save profile...')            
        self.InstEdit.Append(help='Copy instrument profile parameters to other histograms', 
            id=wxID_INSTCOPY, kind=wx.ITEM_NORMAL,text='Copy')
        self.InstEdit.Append(id=wxID_INSTFLAGCOPY, kind=wx.ITEM_NORMAL,text='Copy flags',
            help='Copy instrument parameter refinement flags to other histograms')
#        self.InstEdit.Append(help='Change radiation type (Ka12 - synch)', 
#            id=wxID_CHANGEWAVETYPE, kind=wx.ITEM_NORMAL,text='Change radiation')
        self.InstEdit.Append(id=wxID_INST1VAL, kind=wx.ITEM_NORMAL,text='Set one value',
            help='Set one instrument parameter value across multiple histograms')

        self.PostfillDataMenu()
        
        # PDR / Sample Parameters
        self.SampleMenu = wx.MenuBar()
        self.PrefillDataMenu(self.SampleMenu,helpType='Sample Parameters')
        self.SampleEdit = wx.Menu(title='')
        self.SampleMenu.Append(menu=self.SampleEdit, title='Command')
        self.SetScale = self.SampleEdit.Append(id=wxID_SETSCALE, kind=wx.ITEM_NORMAL,text='Set scale',
            help='Set scale by matching to another histogram')
        self.SampleEdit.Append(id=wxID_SAMPLELOAD, kind=wx.ITEM_NORMAL,text='Load',
            help='Load sample parameters from file')
        self.SampleEdit.Append(id=wxID_SAMPLESAVE, kind=wx.ITEM_NORMAL,text='Save',
            help='Save sample parameters to file')
        self.SampleEdit.Append(id=wxID_SAMPLECOPY, kind=wx.ITEM_NORMAL,text='Copy',
            help='Copy refinable and most other sample parameters to other histograms')
        self.SampleEdit.Append(id=wxID_SAMPLECOPYSOME, kind=wx.ITEM_NORMAL,text='Copy selected...',
            help='Copy selected sample parameters to other histograms')
        self.SampleEdit.Append(id=wxID_SAMPLEFLAGCOPY, kind=wx.ITEM_NORMAL,text='Copy flags',
            help='Copy sample parameter refinement flags to other histograms')
        self.SampleEdit.Append(id=wxID_SAMPLE1VAL, kind=wx.ITEM_NORMAL,text='Set one value',
            help='Set one sample parameter value across multiple histograms')
        self.SampleEdit.Append(id=wxID_ALLSAMPLELOAD, kind=wx.ITEM_NORMAL,text='Load all',
            help='Load sample parmameters over multiple histograms')

        self.PostfillDataMenu()
        self.SetScale.Enable(False)

        # PDR / Peak List
        self.PeakMenu = wx.MenuBar()
        self.PrefillDataMenu(self.PeakMenu,helpType='Peak List')
        self.PeakEdit = wx.Menu(title='')
        self.PeakMenu.Append(menu=self.PeakEdit, title='Peak Fitting')
        self.AutoSearch = self.PeakEdit.Append(help='Automatic peak search', 
            id=wxID_AUTOSEARCH, kind=wx.ITEM_NORMAL,text='Auto search')
        self.UnDo = self.PeakEdit.Append(help='Undo last least squares refinement', 
            id=wxID_UNDO, kind=wx.ITEM_NORMAL,text='UnDo')
        self.PeakFit = self.PeakEdit.Append(id=wxID_LSQPEAKFIT, kind=wx.ITEM_NORMAL,text='Peakfit', 
            help='Peak fitting' )
        self.PFOneCycle = self.PeakEdit.Append(id=wxID_LSQONECYCLE, kind=wx.ITEM_NORMAL,text='Peakfit one cycle', 
            help='One cycle of Peak fitting' )
        self.PeakEdit.Append(id=wxID_RESETSIGGAM, kind=wx.ITEM_NORMAL, 
            text='Reset sig and gam',help='Reset sigma and gamma to global fit' )
        self.PeakCopy = self.PeakEdit.Append(help='Copy peaks to other histograms', 
            id=wxID_PEAKSCOPY, kind=wx.ITEM_NORMAL,text='Peak copy')
        self.SeqPeakFit = self.PeakEdit.Append(id=wxID_SEQPEAKFIT, kind=wx.ITEM_NORMAL,text='Seq PeakFit', 
            help='Sequential Peak fitting for all histograms' )
        self.PeakEdit.Append(id=wxID_CLEARPEAKS, kind=wx.ITEM_NORMAL,text='Clear peaks', 
            help='Clear the peak list' )
        self.PostfillDataMenu()
        self.UnDo.Enable(False)
        self.PeakFit.Enable(False)
        self.PFOneCycle.Enable(False)
        self.AutoSearch.Enable(True)
        
        # PDR / Index Peak List
        self.IndPeaksMenu = wx.MenuBar()
        self.PrefillDataMenu(self.IndPeaksMenu,helpType='Index Peak List')
        self.IndPeaksEdit = wx.Menu(title='')
        self.IndPeaksMenu.Append(menu=self.IndPeaksEdit,title='Operations')
        self.IndPeaksEdit.Append(help='Load/Reload index peaks from peak list',id=wxID_INDXRELOAD, 
            kind=wx.ITEM_NORMAL,text='Load/Reload')
        self.PostfillDataMenu()
        
        # PDR / Unit Cells List
        self.IndexMenu = wx.MenuBar()
        self.PrefillDataMenu(self.IndexMenu,helpType='Unit Cells List')
        self.IndexEdit = wx.Menu(title='')
        self.IndexMenu.Append(menu=self.IndexEdit, title='Cell Index/Refine')
        self.IndexPeaks = self.IndexEdit.Append(help='', id=wxID_INDEXPEAKS, kind=wx.ITEM_NORMAL,
            text='Index Cell')
        self.CopyCell = self.IndexEdit.Append( id=wxID_COPYCELL, kind=wx.ITEM_NORMAL,text='Copy Cell', 
            help='Copy selected unit cell from indexing to cell refinement fields')
        self.RefineCell = self.IndexEdit.Append( id=wxID_REFINECELL, kind=wx.ITEM_NORMAL, 
            text='Refine Cell',help='Refine unit cell parameters from indexed peaks')
        self.MakeNewPhase = self.IndexEdit.Append( id=wxID_MAKENEWPHASE, kind=wx.ITEM_NORMAL,
            text='Make new phase',help='Make new phase from selected unit cell')
        self.ExportCells = self.IndexEdit.Append( id=wxID_EXPORTCELLS, kind=wx.ITEM_NORMAL,
            text='Export cell list',help='Export cell list to csv file')
        self.PostfillDataMenu()
        self.IndexPeaks.Enable(False)
        self.CopyCell.Enable(False)
        self.RefineCell.Enable(False)
        self.MakeNewPhase.Enable(False)
        
        # PDR / Reflection Lists
        self.ReflMenu = wx.MenuBar()
        self.PrefillDataMenu(self.ReflMenu,helpType='Reflection List')
        self.ReflEdit = wx.Menu(title='')
        self.ReflMenu.Append(menu=self.ReflEdit, title='Reflection List')
        self.SelectPhase = self.ReflEdit.Append(help='Select phase for reflection list',id=wxID_SELECTPHASE, 
            kind=wx.ITEM_NORMAL,text='Select phase')
        self.ReflEdit.Append(id=wxID_PWDHKLPLOT,kind=wx.ITEM_NORMAL,text='Plot HKLs',
            help='Plot HKLs from powder pattern')
        self.ReflEdit.Append(id=wxID_PWD3DHKLPLOT,kind=wx.ITEM_NORMAL,text='Plot 3D HKLs',
            help='Plot HKLs from powder pattern in 3D')
        self.PostfillDataMenu()
        
        # SASD / Instrument Parameters
        self.SASDInstMenu = wx.MenuBar()
        self.PrefillDataMenu(self.SASDInstMenu,helpType='Instrument Parameters')
        self.SASDInstEdit = wx.Menu(title='')
        self.SASDInstMenu.Append(menu=self.SASDInstEdit, title='Operations')
        self.InstEdit.Append(help='Reset instrument profile parameters to default', 
            id=wxID_INSTPRMRESET, kind=wx.ITEM_NORMAL,text='Reset profile')
        self.SASDInstEdit.Append(help='Copy instrument profile parameters to other histograms', 
            id=wxID_INSTCOPY, kind=wx.ITEM_NORMAL,text='Copy')
        self.PostfillDataMenu()
        
        #SASD & REFL/ Substance editor
        self.SubstanceMenu = wx.MenuBar()
        self.PrefillDataMenu(self.SubstanceMenu,helpType='Substances')
        self.SubstanceEdit = wx.Menu(title='')
        self.SubstanceMenu.Append(menu=self.SubstanceEdit, title='Edit')
        self.SubstanceEdit.Append(id=wxID_LOADSUBSTANCE, kind=wx.ITEM_NORMAL,text='Load substance',
            help='Load substance from file')
        self.SubstanceEdit.Append(id=wxID_ADDSUBSTANCE, kind=wx.ITEM_NORMAL,text='Add substance',
            help='Add new substance to list')
        self.SubstanceEdit.Append(id=wxID_COPYSUBSTANCE, kind=wx.ITEM_NORMAL,text='Copy substances',
            help='Copy substances')
        self.SubstanceEdit.Append(id=wxID_DELETESUBSTANCE, kind=wx.ITEM_NORMAL,text='Delete substance',
            help='Delete substance from list')            
        self.SubstanceEdit.Append(id=wxID_ELEMENTADD, kind=wx.ITEM_NORMAL,text='Add elements',
            help='Add elements to substance')
        self.SubstanceEdit.Append(id=wxID_ELEMENTDELETE, kind=wx.ITEM_NORMAL,text='Delete elements',
            help='Delete elements from substance')
        self.PostfillDataMenu()
        
        # SASD/ Models
        self.ModelMenu = wx.MenuBar()
        self.PrefillDataMenu(self.ModelMenu,helpType='Models')
        self.ModelEdit = wx.Menu(title='')
        self.ModelMenu.Append(menu=self.ModelEdit, title='Models')
        self.ModelEdit.Append(id=wxID_MODELADD,kind=wx.ITEM_NORMAL,text='Add',
            help='Add new term to model')
        self.ModelEdit.Append(id=wxID_MODELFIT, kind=wx.ITEM_NORMAL,text='Fit',
            help='Fit model parameters to data')
        self.SasdUndo = self.ModelEdit.Append(id=wxID_MODELUNDO, kind=wx.ITEM_NORMAL,text='Undo',
            help='Undo model fit')
        self.SasdUndo.Enable(False)            
        self.ModelEdit.Append(id=wxID_MODELFITALL, kind=wx.ITEM_NORMAL,text='Sequential fit',
            help='Sequential fit of model parameters to all SASD data')
        self.ModelEdit.Append(id=wxID_MODELCOPY, kind=wx.ITEM_NORMAL,text='Copy',
            help='Copy model parameters to other histograms')
        self.ModelEdit.Append(id=wxID_MODELCOPYFLAGS, kind=wx.ITEM_NORMAL,text='Copy flags',
            help='Copy model refinement flags to other histograms')
        self.PostfillDataMenu()
        
        # IMG / Image Controls
        self.ImageMenu = wx.MenuBar()
        self.PrefillDataMenu(self.ImageMenu,helpType='Image Controls')
        self.ImageEdit = wx.Menu(title='')
        self.ImageMenu.Append(menu=self.ImageEdit, title='Operations')
        self.ImageEdit.Append(help='Calibrate detector by fitting to calibrant lines', 
            id=wxID_IMCALIBRATE, kind=wx.ITEM_NORMAL,text='Calibrate')
        self.ImageEdit.Append(help='Recalibrate detector by fitting to calibrant lines', 
            id=wxID_IMRECALIBRATE, kind=wx.ITEM_NORMAL,text='Recalibrate')
        self.ImageEdit.Append(help='Clear calibration data points and rings',id=wxID_IMCLEARCALIB, 
            kind=wx.ITEM_NORMAL,text='Clear calibration')
        self.ImageEdit.Append(help='Integrate selected image',id=wxID_IMINTEGRATE, 
            kind=wx.ITEM_NORMAL,text='Integrate')
        self.ImageEdit.Append(help='Integrate all images selected from list',id=wxID_INTEGRATEALL,
            kind=wx.ITEM_NORMAL,text='Integrate all')
        self.ImageEdit.Append(help='Copy image controls to other images', 
            id=wxID_IMCOPYCONTROLS, kind=wx.ITEM_NORMAL,text='Copy Controls')
        self.ImageEdit.Append(help='Save image controls to file', 
            id=wxID_IMSAVECONTROLS, kind=wx.ITEM_NORMAL,text='Save Controls')
        self.ImageEdit.Append(help='Load image controls from file', 
            id=wxID_IMLOADCONTROLS, kind=wx.ITEM_NORMAL,text='Load Controls')
        self.PostfillDataMenu()
            
        # IMG / Masks
        self.MaskMenu = wx.MenuBar()
        self.PrefillDataMenu(self.MaskMenu,helpType='Image Masks')
        self.MaskEdit = wx.Menu(title='')
        self.MaskMenu.Append(menu=self.MaskEdit, title='Operations')
        submenu = wx.Menu()
        self.MaskEdit.AppendMenu(
            wx.ID_ANY,'Create new', submenu,
            help=''
            )
        self.MaskEdit.Append(help='Copy mask to other images', 
            id=wxID_MASKCOPY, kind=wx.ITEM_NORMAL,text='Copy mask')
        self.MaskEdit.Append(help='Save mask to file', 
            id=wxID_MASKSAVE, kind=wx.ITEM_NORMAL,text='Save mask')
        self.MaskEdit.Append(help='Load mask from file', 
            id=wxID_MASKLOAD, kind=wx.ITEM_NORMAL,text='Load mask')
        self.MaskEdit.Append(help='Load mask from file; ignore threshold', 
            id=wxID_MASKLOADNOT, kind=wx.ITEM_NORMAL,text='Load mask w/o threshold')
        submenu.Append(help='Create an arc mask with mouse input', 
            id=wxID_NEWMASKARC, kind=wx.ITEM_NORMAL,text='Arc mask')
        submenu.Append(help='Create a frame mask with mouse input', 
            id=wxID_NEWMASKFRAME, kind=wx.ITEM_NORMAL,text='Frame mask')
        submenu.Append(help='Create a polygon mask with mouse input', 
            id=wxID_NEWMASKPOLY, kind=wx.ITEM_NORMAL,text='Polygon mask')
        submenu.Append(help='Create a ring mask with mouse input', 
            id=wxID_NEWMASKRING, kind=wx.ITEM_NORMAL,text='Ring mask')
        submenu.Append(help='Create a spot mask with mouse input', 
            id=wxID_NEWMASKSPOT, kind=wx.ITEM_NORMAL,text='Spot mask')
        self.PostfillDataMenu()
            
        # IMG / Stress/Strain
        self.StrStaMenu = wx.MenuBar()
        self.PrefillDataMenu(self.StrStaMenu,helpType='Stress/Strain')
        self.StrStaEdit = wx.Menu(title='')
        self.StrStaMenu.Append(menu=self.StrStaEdit, title='Operations')
        self.StrStaEdit.Append(help='Append d-zero for one ring', 
            id=wxID_APPENDDZERO, kind=wx.ITEM_NORMAL,text='Append d-zero')
        self.StrStaEdit.Append(help='Fit stress/strain data', 
            id=wxID_STRSTAFIT, kind=wx.ITEM_NORMAL,text='Fit stress/strain')
        self.StrStaEdit.Append(help='Update d-zero from ave d-zero',
            id=wxID_UPDATEDZERO, kind=wx.ITEM_NORMAL,text='Update d-zero')        
        self.StrStaEdit.Append(help='Fit stress/strain data for all images', 
            id=wxID_STRSTAALLFIT, kind=wx.ITEM_NORMAL,text='All image fit')
        self.StrStaEdit.Append(help='Copy stress/strain data to other images', 
            id=wxID_STRSTACOPY, kind=wx.ITEM_NORMAL,text='Copy stress/strain')
        self.StrStaEdit.Append(help='Save stress/strain data to file', 
            id=wxID_STRSTASAVE, kind=wx.ITEM_NORMAL,text='Save stress/strain')
        self.StrStaEdit.Append(help='Load stress/strain data from file', 
            id=wxID_STRSTALOAD, kind=wx.ITEM_NORMAL,text='Load stress/strain')
        self.StrStaEdit.Append(help='Load sample data from file', 
            id=wxID_STRSTSAMPLE, kind=wx.ITEM_NORMAL,text='Load sample data')
        self.PostfillDataMenu()
            
        # PDF / PDF Controls
        self.PDFMenu = wx.MenuBar()
        self.PrefillDataMenu(self.PDFMenu,helpType='PDF Controls')
        self.PDFEdit = wx.Menu(title='')
        self.PDFMenu.Append(menu=self.PDFEdit, title='PDF Controls')
        self.PDFEdit.Append(help='Add element to sample composition',id=wxID_PDFADDELEMENT, kind=wx.ITEM_NORMAL,
            text='Add element')
        self.PDFEdit.Append(help='Delete element from sample composition',id=wxID_PDFDELELEMENT, kind=wx.ITEM_NORMAL,
            text='Delete element')
        self.PDFEdit.Append(help='Copy PDF controls', id=wxID_PDFCOPYCONTROLS, kind=wx.ITEM_NORMAL,
            text='Copy controls')
        self.PDFEdit.Append(help='Load PDF controls from file',id=wxID_PDFLOADCONTROLS, kind=wx.ITEM_NORMAL,
            text='Load Controls')
        self.PDFEdit.Append(help='Save PDF controls to file', id=wxID_PDFSAVECONTROLS, kind=wx.ITEM_NORMAL,
            text='Save controls')
        self.PDFEdit.Append(help='Compute PDF', id=wxID_PDFCOMPUTE, kind=wx.ITEM_NORMAL,
            text='Compute PDF')
        self.PDFEdit.Append(help='Compute all PDFs', id=wxID_PDFCOMPUTEALL, kind=wx.ITEM_NORMAL,
            text='Compute all PDFs')
        self.PostfillDataMenu()
        
        # Phase / General tab
        self.DataGeneral = wx.MenuBar()
        self.PrefillDataMenu(self.DataGeneral,helpType='General', helpLbl='Phase/General')
        self.DataGeneral.Append(menu=wx.Menu(title=''),title='Select tab')
        self.GeneralCalc = wx.Menu(title='')
        self.DataGeneral.Append(menu=self.GeneralCalc,title='Compute')
        self.GeneralCalc.Append(help='Compute Fourier map',id=wxID_FOURCALC, kind=wx.ITEM_NORMAL,
            text='Fourier map')
        self.GeneralCalc.Append(help='Search Fourier map',id=wxID_FOURSEARCH, kind=wx.ITEM_NORMAL,
            text='Search map')
        self.GeneralCalc.Append(help='Run charge flipping',id=wxID_CHARGEFLIP, kind=wx.ITEM_NORMAL,
            text='Charge flipping')
        self.GeneralCalc.Append(help='Run 4D charge flipping',id=wxID_4DCHARGEFLIP, kind=wx.ITEM_NORMAL,
            text='4D Charge flipping')
        self.GeneralCalc.Enable(wxID_4DCHARGEFLIP,False)   
        self.GeneralCalc.Append(help='Clear map',id=wxID_FOURCLEAR, kind=wx.ITEM_NORMAL,
            text='Clear map')
        self.GeneralCalc.Append(help='Run Monte Carlo - Simulated Annealing',id=wxID_SINGLEMCSA, kind=wx.ITEM_NORMAL,
            text='MC/SA')
        self.GeneralCalc.Append(help='Run Monte Carlo - Simulated Annealing on multiprocessors',id=wxID_MULTIMCSA, kind=wx.ITEM_NORMAL,
            text='Multi MC/SA')            #currently not useful
        self.PostfillDataMenu()
        
        # Phase / Data tab
        self.DataMenu = wx.MenuBar()
        self.PrefillDataMenu(self.DataMenu,helpType='Data', helpLbl='Phase/Data')
        self.DataMenu.Append(menu=wx.Menu(title=''),title='Select tab')
        self.DataEdit = wx.Menu(title='')
        self.DataMenu.Append(menu=self.DataEdit, title='Edit')
        self.DataEdit.Append(id=wxID_DATACOPY, kind=wx.ITEM_NORMAL,text='Copy data',
            help='Copy phase data to other histograms')
        self.DataEdit.Append(id=wxID_DATACOPYFLAGS, kind=wx.ITEM_NORMAL,text='Copy flags',
            help='Copy phase data flags to other histograms')
        self.DataEdit.Append(id=wxID_DATASELCOPY, kind=wx.ITEM_NORMAL,text='Copy selected data',
            help='Copy selected phase data to other histograms')
        self.DataEdit.Append(id=wxID_PWDRADD, kind=wx.ITEM_NORMAL,text='Add powder histograms',
            help='Select new powder histograms to be used for this phase')
        self.DataEdit.Append(id=wxID_HKLFADD, kind=wx.ITEM_NORMAL,text='Add single crystal histograms',
            help='Select new single crystal histograms to be used for this phase')
        self.DataEdit.Append(id=wxID_DATADELETE, kind=wx.ITEM_NORMAL,text='Remove histograms',
            help='Remove histograms from use for this phase')
        self.PostfillDataMenu()
            
        # Phase / Atoms tab
        self.AtomsMenu = wx.MenuBar()
        self.PrefillDataMenu(self.AtomsMenu,helpType='Atoms')
        self.AtomsMenu.Append(menu=wx.Menu(title=''),title='Select tab')
        self.AtomEdit = wx.Menu(title='')
        self.AtomCompute = wx.Menu(title='')
        self.AtomsMenu.Append(menu=self.AtomEdit, title='Edit')
        self.AtomsMenu.Append(menu=self.AtomCompute, title='Compute')
        self.AtomEdit.Append(id=wxID_ATOMSEDITADD, kind=wx.ITEM_NORMAL,text='Append atom',
            help='Appended as an H atom')
        self.AtomEdit.Append(id=wxID_ATOMSVIEWADD, kind=wx.ITEM_NORMAL,text='Append view point',
            help='Appended as an H atom')
        self.AtomEdit.Append(id=wxID_ATOMSEDITINSERT, kind=wx.ITEM_NORMAL,text='Insert atom',
            help='Select atom row to insert before; inserted as an H atom')
        self.AtomEdit.Append(id=wxID_ATOMVIEWINSERT, kind=wx.ITEM_NORMAL,text='Insert view point',
            help='Select atom row to insert before; inserted as an H atom')
        self.AtomEdit.Append(id=wxID_ATOMMOVE, kind=wx.ITEM_NORMAL,text='Move atom to view point',
            help='Select single atom to move')
        self.AtomEdit.Append(id=wxID_ATOMSEDITDELETE, kind=wx.ITEM_NORMAL,text='Delete atom',
            help='Select atoms to delete first')
        self.AtomEdit.Append(id=wxID_ATOMSREFINE, kind=wx.ITEM_NORMAL,text='Set atom refinement flags',
            help='Select atoms to refine first')
        self.AtomEdit.Append(id=wxID_ATOMSMODIFY, kind=wx.ITEM_NORMAL,text='Modify atom parameters',
            help='Select atoms to modify first')
        self.AtomEdit.Append(id=wxID_ATOMSTRANSFORM, kind=wx.ITEM_NORMAL,text='Transform atoms',
            help='Select atoms to transform first')
        self.AtomEdit.Append(id=wxID_MAKEMOLECULE, kind=wx.ITEM_NORMAL,text='Assemble molecule',
            help='Assemble molecule from scatterd atom positions')
        self.AtomEdit.Append(id=wxID_RELOADDRAWATOMS, kind=wx.ITEM_NORMAL,text='Reload draw atoms',
            help='Reload atom drawing list')
        submenu = wx.Menu()
        self.AtomEdit.AppendMenu(wx.ID_ANY, 'Reimport atoms', submenu, 
            help='Reimport atoms from file; sequence must match')
        # setup a cascade menu for the formats that have been defined
        self.ReImportMenuId = {}  # points to readers for each menu entry
        for reader in self.G2frame.ImportPhaseReaderlist:
            item = submenu.Append(
                wx.ID_ANY,help=reader.longFormatName,
                kind=wx.ITEM_NORMAL,text='reimport coordinates from '+reader.formatName+' file')
            self.ReImportMenuId[item.GetId()] = reader
        item = submenu.Append(
            wx.ID_ANY,
            help='Reimport coordinates, try to determine format from file',
            kind=wx.ITEM_NORMAL,
            text='guess format from file')
        self.ReImportMenuId[item.GetId()] = None # try all readers

        self.AtomCompute.Append(id=wxID_ATOMSDISAGL, kind=wx.ITEM_NORMAL,text='Show Distances && Angles',
            help='Compute distances & angles for selected atoms')
        self.AtomCompute.Append(id=wxID_ATOMSPDISAGL, kind=wx.ITEM_NORMAL,text='Save Distances && Angles',
            help='Compute distances & angles for selected atoms')
        self.AtomCompute.ISOcalc = self.AtomCompute.Append(
            id=wxID_ISODISP, kind=wx.ITEM_NORMAL,
            text='Compute ISODISTORT mode values',
            help='Compute values of ISODISTORT modes from atom parameters')
        self.PostfillDataMenu()
        
        # Phase / Imcommensurate "waves" tab
        self.WavesData = wx.MenuBar()
        self.PrefillDataMenu(self.WavesData,helpType='Wave Data', helpLbl='Imcommensurate wave data')
        self.WavesData.Append(menu=wx.Menu(title=''),title='Select tab')
        self.WavesDataCompute = wx.Menu(title='')
        self.WavesData.Append(menu=self.WavesDataCompute,title='Compute')
        self.WavesDataCompute.Append(id=wxID_4DMAPCOMPUTE, kind=wx.ITEM_NORMAL,text='Compute 4D map',
            help='Compute 4-dimensional map')
        self.PostfillDataMenu()
                 
        # Phase / Draw Options tab
        self.DataDrawOptions = wx.MenuBar()
        self.PrefillDataMenu(self.DataDrawOptions,helpType='Draw Options', helpLbl='Phase/Draw Options')
        self.DataDrawOptions.Append(menu=wx.Menu(title=''),title='Select tab')
        self.PostfillDataMenu()
        
        # Phase / Draw Atoms tab 
        self.DrawAtomsMenu = wx.MenuBar()
        self.PrefillDataMenu(self.DrawAtomsMenu,helpType='Draw Atoms')
        self.DrawAtomsMenu.Append(menu=wx.Menu(title=''),title='Select tab')
        self.DrawAtomEdit = wx.Menu(title='')
        self.DrawAtomCompute = wx.Menu(title='')
        self.DrawAtomRestraint = wx.Menu(title='')
        self.DrawAtomRigidBody = wx.Menu(title='')
        self.DrawAtomsMenu.Append(menu=self.DrawAtomEdit, title='Edit')
        self.DrawAtomsMenu.Append(menu=self.DrawAtomCompute,title='Compute')
        self.DrawAtomsMenu.Append(menu=self.DrawAtomRestraint, title='Restraints')
        self.DrawAtomsMenu.Append(menu=self.DrawAtomRigidBody, title='Rigid body')
        self.DrawAtomEdit.Append(id=wxID_DRAWATOMSTYLE, kind=wx.ITEM_NORMAL,text='Atom style',
            help='Select atoms first')
        self.DrawAtomEdit.Append(id=wxID_DRAWATOMLABEL, kind=wx.ITEM_NORMAL,text='Atom label',
            help='Select atoms first')
        self.DrawAtomEdit.Append(id=wxID_DRAWATOMCOLOR, kind=wx.ITEM_NORMAL,text='Atom color',
            help='Select atoms first')
        self.DrawAtomEdit.Append(id=wxID_DRAWATOMRESETCOLOR, kind=wx.ITEM_NORMAL,text='Reset atom colors',
            help='Resets all atom colors to defaults')
        self.DrawAtomEdit.Append(id=wxID_DRAWVIEWPOINT, kind=wx.ITEM_NORMAL,text='View point',
            help='View point is 1st atom selected')
        self.DrawAtomEdit.Append(id=wxID_DRAWADDEQUIV, kind=wx.ITEM_NORMAL,text='Add atoms',
            help='Add symmetry & cell equivalents to drawing set from selected atoms')
        self.DrawAtomEdit.Append(id=wxID_DRAWTRANSFORM, kind=wx.ITEM_NORMAL,text='Transform draw atoms',
            help='Transform selected atoms by symmetry & cell translations')
        self.DrawAtomEdit.Append(id=wxID_DRAWFILLCOORD, kind=wx.ITEM_NORMAL,text='Fill CN-sphere',
            help='Fill coordination sphere for selected atoms')            
        self.DrawAtomEdit.Append(id=wxID_DRAWFILLCELL, kind=wx.ITEM_NORMAL,text='Fill unit cell',
            help='Fill unit cell with selected atoms')
        self.DrawAtomEdit.Append(id=wxID_DRAWDELETE, kind=wx.ITEM_NORMAL,text='Delete atoms',
            help='Delete atoms from drawing set')
        self.DrawAtomCompute.Append(id=wxID_DRAWDISTVP, kind=wx.ITEM_NORMAL,text='View pt. dist.',
            help='Compute distance of selected atoms from view point')   
        self.DrawAtomCompute.Append(id=wxID_DRAWDISAGLTOR, kind=wx.ITEM_NORMAL,text='Dist. Ang. Tors.',
            help='Compute distance, angle or torsion for 2-4 selected atoms')   
        self.DrawAtomCompute.Append(id=wxID_DRAWPLANE, kind=wx.ITEM_NORMAL,text='Best plane',
            help='Compute best plane for 4+ selected atoms')   
        self.DrawAtomRestraint.Append(id=wxID_DRAWRESTRBOND, kind=wx.ITEM_NORMAL,text='Add bond restraint',
            help='Add bond restraint for selected atoms (2)')
        self.DrawAtomRestraint.Append(id=wxID_DRAWRESTRANGLE, kind=wx.ITEM_NORMAL,text='Add angle restraint',
            help='Add angle restraint for selected atoms (3: one end 1st)')
        self.DrawAtomRestraint.Append(id=wxID_DRAWRESTRPLANE, kind=wx.ITEM_NORMAL,text='Add plane restraint',
            help='Add plane restraint for selected atoms (4+)')
        self.DrawAtomRestraint.Append(id=wxID_DRAWRESTRCHIRAL, kind=wx.ITEM_NORMAL,text='Add chiral restraint',
            help='Add chiral restraint for selected atoms (4: center atom 1st)')
        self.DrawAtomRigidBody.Append(id=wxID_DRAWDEFINERB, kind=wx.ITEM_NORMAL,text='Define rigid body',
            help='Define rigid body with selected atoms')
        self.PostfillDataMenu()

        # Phase / MCSA tab
        self.MCSAMenu = wx.MenuBar()
        self.PrefillDataMenu(self.MCSAMenu,helpType='MC/SA')
        self.MCSAMenu.Append(menu=wx.Menu(title=''),title='Select tab')
        self.MCSAEdit = wx.Menu(title='')
        self.MCSAMenu.Append(menu=self.MCSAEdit, title='MC/SA')
        self.MCSAEdit.Append(id=wxID_ADDMCSAATOM, kind=wx.ITEM_NORMAL,text='Add atom', 
            help='Add single atom to MC/SA model')
        self.MCSAEdit.Append(id=wxID_ADDMCSARB, kind=wx.ITEM_NORMAL,text='Add rigid body', 
            help='Add rigid body to MC/SA model' )
        self.MCSAEdit.Append(id=wxID_CLEARMCSARB, kind=wx.ITEM_NORMAL,text='Clear rigid bodies', 
            help='Clear all atoms & rigid bodies from MC/SA model' )
        self.MCSAEdit.Append(id=wxID_MOVEMCSA, kind=wx.ITEM_NORMAL,text='Move MC/SA solution', 
            help='Move MC/SA solution to atom list' )
        self.MCSAEdit.Append(id=wxID_MCSACLEARRESULTS, kind=wx.ITEM_NORMAL,text='Clear results', 
            help='Clear table of MC/SA results' )
        self.PostfillDataMenu()
            
        # Phase / Texture tab
        self.TextureMenu = wx.MenuBar()
        self.PrefillDataMenu(self.TextureMenu,helpType='Texture')
        self.TextureMenu.Append(menu=wx.Menu(title=''),title='Select tab')
        self.TextureEdit = wx.Menu(title='')
        self.TextureMenu.Append(menu=self.TextureEdit, title='Texture')
        self.TextureEdit.Append(id=wxID_REFINETEXTURE, kind=wx.ITEM_NORMAL,text='Refine texture', 
            help='Refine the texture coefficients from sequential results')
#        self.TextureEdit.Append(id=wxID_CLEARTEXTURE, kind=wx.ITEM_NORMAL,text='Clear texture', 
#            help='Clear the texture coefficients' )
        self.PostfillDataMenu()
            
        # Phase / Pawley tab
        self.PawleyMenu = wx.MenuBar()
        self.PrefillDataMenu(self.PawleyMenu,helpType='Pawley')
        self.PawleyMenu.Append(menu=wx.Menu(title=''),title='Select tab')
        self.PawleyEdit = wx.Menu(title='')
        self.PawleyMenu.Append(menu=self.PawleyEdit,title='Operations')
        self.PawleyEdit.Append(id=wxID_PAWLEYLOAD, kind=wx.ITEM_NORMAL,text='Pawley create',
            help='Initialize Pawley reflection list')
        self.PawleyEdit.Append(id=wxID_PAWLEYESTIMATE, kind=wx.ITEM_NORMAL,text='Pawley estimate',
            help='Estimate initial Pawley intensities')
        self.PawleyEdit.Append(id=wxID_PAWLEYUPDATE, kind=wx.ITEM_NORMAL,text='Pawley update',
            help='Update negative Pawley intensities with -0.5*Fobs and turn off refinemnt')
        self.PostfillDataMenu()
            
        # Phase / Map peaks tab
        self.MapPeaksMenu = wx.MenuBar()
        self.PrefillDataMenu(self.MapPeaksMenu,helpType='Map peaks')
        self.MapPeaksMenu.Append(menu=wx.Menu(title=''),title='Select tab')
        self.MapPeaksEdit = wx.Menu(title='')
        self.MapPeaksMenu.Append(menu=self.MapPeaksEdit, title='Map peaks')
        self.MapPeaksEdit.Append(id=wxID_PEAKSMOVE, kind=wx.ITEM_NORMAL,text='Move peaks', 
            help='Move selected peaks to atom list')
        self.MapPeaksEdit.Append(id=wxID_PEAKSVIEWPT, kind=wx.ITEM_NORMAL,text='View point',
            help='View point is 1st peak selected')
        self.MapPeaksEdit.Append(id=wxID_PEAKSDISTVP, kind=wx.ITEM_NORMAL,text='View pt. dist.',
            help='Compute distance of selected peaks from view point')   
        self.MapPeaksEdit.Append(id=wxID_SHOWBONDS, kind=wx.ITEM_NORMAL,text='Hide bonds',
            help='Hide or show bonds between peak positions')   
        self.MapPeaksEdit.Append(id=wxID_PEAKSDA, kind=wx.ITEM_NORMAL,text='Calc dist/ang', 
            help='Calculate distance or angle for selection')
        self.MapPeaksEdit.Append(id=wxID_FINDEQVPEAKS, kind=wx.ITEM_NORMAL,text='Equivalent peaks', 
            help='Find equivalent peaks')
        self.MapPeaksEdit.Append(id=wxID_PEAKSUNIQUE, kind=wx.ITEM_NORMAL,text='Unique peaks', 
            help='Select unique set')
        self.MapPeaksEdit.Append(id=wxID_PEAKSDELETE, kind=wx.ITEM_NORMAL,text='Delete peaks', 
            help='Delete selected peaks')
        self.MapPeaksEdit.Append(id=wxID_PEAKSCLEAR, kind=wx.ITEM_NORMAL,text='Clear peaks', 
            help='Clear the map peak list')
        self.PostfillDataMenu()

        # Phase / Rigid bodies tab
        self.RigidBodiesMenu = wx.MenuBar()
        self.PrefillDataMenu(self.RigidBodiesMenu,helpType='Rigid bodies')
        self.RigidBodiesMenu.Append(menu=wx.Menu(title=''),title='Select tab')
        self.RigidBodiesEdit = wx.Menu(title='')
        self.RigidBodiesMenu.Append(menu=self.RigidBodiesEdit, title='Edit')
        self.RigidBodiesEdit.Append(id=wxID_ASSIGNATMS2RB, kind=wx.ITEM_NORMAL,text='Assign atoms to rigid body',
            help='Select & position rigid body in structure of existing atoms')
        self.RigidBodiesEdit.Append(id=wxID_AUTOFINDRESRB, kind=wx.ITEM_NORMAL,text='Auto find residues',
            help='Auto find of residue RBs in macromolecule')
        self.RigidBodiesEdit.Append(id=wxID_COPYRBPARMS, kind=wx.ITEM_NORMAL,text='Copy rigid body parms',
            help='Copy rigid body location & TLS parameters')
        self.RigidBodiesEdit.Append(id=wxID_GLOBALTHERM, kind=wx.ITEM_NORMAL,text='Global thermal motion',
            help='Global setting of residue thermal motion models')
        self.RigidBodiesEdit.Append(id=wxID_GLOBALRESREFINE, kind=wx.ITEM_NORMAL,text='Global residue refine',
            help='Global setting of residue RB refinement flags')
        self.RigidBodiesEdit.Append(id=wxID_RBREMOVEALL, kind=wx.ITEM_NORMAL,text='Remove all rigid bodies',
            help='Remove all rigid body assignment for atoms')
        self.PostfillDataMenu()
    # end of GSAS-II menu definitions
        
    def _init_ctrls(self, parent,name=None,size=None,pos=None):
        wx.Frame.__init__(
            self,parent=parent,
            #style=wx.DEFAULT_FRAME_STYLE ^ wx.CLOSE_BOX | wx.FRAME_FLOAT_ON_PARENT ,
            style=wx.DEFAULT_FRAME_STYLE ^ wx.CLOSE_BOX,
            size=size,pos=pos,title='GSAS-II data display')
        self._init_menus()
        if name:
            self.SetLabel(name)
        self.Show()
        
    def __init__(self,parent,frame,data=None,name=None, size=None,pos=None):
        self.G2frame = frame
        self._init_ctrls(parent,name,size,pos)
        self.data = data
        clientSize = wx.ClientDisplayRect()
        Size = self.GetSize()
        xPos = clientSize[2]-Size[0]
        self.SetPosition(wx.Point(xPos,clientSize[1]+250))
        self.AtomGrid = []
        self.selectedRow = 0
        
    def setSizePosLeft(self,Width):
        clientSize = wx.ClientDisplayRect()
        Width[1] = min(Width[1],clientSize[2]-300)
        Width[0] = max(Width[0],300)
        self.SetSize(Width)
#        self.SetPosition(wx.Point(clientSize[2]-Width[0],clientSize[1]+250))
        
    def Clear(self):
        self.ClearBackground()
        self.DestroyChildren()
                   

################################################################################
#####  Notebook Tree Item editor
################################################################################                  
def UpdateNotebook(G2frame,data):
    '''Called when the data tree notebook entry is selected. Allows for
    editing of the text in that tree entry
    '''
    def OnNoteBook(event):
        data = G2frame.dataDisplay.GetValue().split('\n')
        G2frame.PatternTree.SetItemPyData(GetPatternTreeItemId(G2frame,G2frame.root,'Notebook'),data)
        if 'nt' not in os.name:
            G2frame.dataDisplay.AppendText('\n')
                    
    if G2frame.dataDisplay:
        G2frame.dataDisplay.Destroy()
    G2frame.dataFrame.SetLabel('Notebook')
    G2frame.dataDisplay = wx.TextCtrl(parent=G2frame.dataFrame,size=G2frame.dataFrame.GetClientSize(),
        style=wx.TE_MULTILINE|wx.TE_PROCESS_ENTER | wx.TE_DONTWRAP)
    G2frame.dataDisplay.Bind(wx.EVT_TEXT_ENTER,OnNoteBook)
    G2frame.dataDisplay.Bind(wx.EVT_KILL_FOCUS,OnNoteBook)
    for line in data:
        G2frame.dataDisplay.AppendText(line+"\n")
    G2frame.dataDisplay.AppendText('Notebook entry @ '+time.ctime()+"\n")
    G2frame.dataFrame.setSizePosLeft([400,250])
            
################################################################################
#####  Controls Tree Item editor
################################################################################           
def UpdateControls(G2frame,data):
    '''Edit overall GSAS-II controls in main Controls data tree entry
    '''
    #patch
    if 'deriv type' not in data:
        data = {}
        data['deriv type'] = 'analytic Hessian'
        data['min dM/M'] = 0.0001
        data['shift factor'] = 1.
        data['max cyc'] = 3        
        data['F**2'] = False
    if 'shift factor' not in data:
        data['shift factor'] = 1.
    if 'max cyc' not in data:
        data['max cyc'] = 3
    if 'F**2' not in data:
        data['F**2'] = False
    if 'Author' not in data:
        data['Author'] = 'no name'
    if 'FreePrm1' not in data:
        data['FreePrm1'] = 'Sample humidity (%)'
    if 'FreePrm2' not in data:
        data['FreePrm2'] = 'Sample voltage (V)'
    if 'FreePrm3' not in data:
        data['FreePrm3'] = 'Applied load (MN)'
    if 'Copy2Next' not in data:
        data['Copy2Next'] = False
    if 'Reverse Seq' not in data:
        data['Reverse Seq'] = False
    if 'UsrReject' not in data:
        data['UsrReject'] = {'minF/sig':0,'MinExt':0.01,'MaxDF/F':20.,'MaxD':500.,'MinD':0.05}
# </ Anton Gagin     
# additional controls in main Controls data tree entry
# multiplicative controls
    if 'corrParam E_mu' not in data:

        data['corrParam E_mu'] = str(0)  
    if 'corrParam k_mu' not in data:
        data['corrParam k_mu'] = str(0)    
    if 'corrParam k_mu tmp' not in data:
        data['corrParam k_mu tmp'] = str(0)  
# additive controls
    if 'corrParam E_beta' not in data:
        data['corrParam E_beta'] = str(0)  
    if 'corrParam k_beta' not in data:
        data['corrParam k_beta'] = str(0)    
    if 'corrParam k_beta tmp' not in data:
        data['corrParam k_beta tmp'] = str(0)         
# peak-shape controls
    if 'corrParam sigma_delta' not in data:
        data['corrParam sigma_delta'] = str(0)   
    if 'corrParam l_delta' not in data:
        data['corrParam l_delta'] = str(0)  
    if 'corrParam l_delta tmp' not in data:
        data['corrParam l_delta tmp'] = str(0)          
    if 'corrParam num blocks s' not in data:
        data['corrParam num blocks s'] = str(0)
    if 'corrParam FWHMDivN' not in data:
        data['corrParam FWHMDivN'] = "none" 
        # correlation length l_delta can be calculated as mean(FWHM/FWHMDivN)        
    if 'corrParam l_deltaDivN' not in data:
        data['corrParam l_deltaDivN'] = "none"       
        # stdev sigma_delta can be calculated as l_delta/l_deltaDivN      
    if 'EstimateKMu' not in data:
        data['EstimateKMu'] = False  
    if 'EstimateKBeta' not in data:
        data['EstimateKBeta'] = False  
    if 'doIter' not in data:
        data['doIter'] = False   
# MCMC controls       
    if 'doMCMC' not in data:
        data['doMCMC'] = False       
    if 'nwalkers' not in data:
        data['nwalkers'] = str(0)   
    if 'nIterMCMC' not in data:
        data['nIterMCMC'] = str(0)  
        
# Anton Gagin />  
     
    
    #end patch

    def SeqSizer():
        
        def OnSelectData(event):
            choices = GetPatternTreeDataNames(G2frame,['PWDR','HKLF',])
            sel = []
            try:
                if 'Seq Data' in data:
                    for item in data['Seq Data']:
                        sel.append(choices.index(item))
                    sel = [choices.index(item) for item in data['Seq Data']]
            except ValueError:  #data changed somehow - start fresh
                sel = []
            dlg = G2G.G2MultiChoiceDialog(G2frame.dataFrame, 'Sequential refinement',
                'Select dataset to include',choices)
            dlg.SetSelections(sel)
            names = []
            if dlg.ShowModal() == wx.ID_OK:
                for sel in dlg.GetSelections():
                    names.append(choices[sel])
                data['Seq Data'] = names                
                G2frame.EnableSeqRefineMenu()
            dlg.Destroy()
            wx.CallAfter(UpdateControls,G2frame,data)
            
        def OnReverse(event):
            data['Reverse Seq'] = reverseSel.GetValue()
            
        def OnCopySel(event):
            data['Copy2Next'] = copySel.GetValue() 
                    
        seqSizer = wx.BoxSizer(wx.VERTICAL)
        dataSizer = wx.BoxSizer(wx.HORIZONTAL)
        dataSizer.Add(wx.StaticText(G2frame.dataDisplay,label=' Sequential Refinement: '),0,WACV)
        selSeqData = wx.Button(G2frame.dataDisplay,-1,label=' Select data')
        selSeqData.Bind(wx.EVT_BUTTON,OnSelectData)
        dataSizer.Add(selSeqData,0,WACV)
        SeqData = data.get('Seq Data',[])
        if not SeqData:
            lbl = ' (no data selected)'
        else:
            lbl = ' ('+str(len(SeqData))+' dataset(s) selected)'

        dataSizer.Add(wx.StaticText(G2frame.dataDisplay,label=lbl),0,WACV)
        seqSizer.Add(dataSizer,0)
        if SeqData:
            selSizer = wx.BoxSizer(wx.HORIZONTAL)
            reverseSel = wx.CheckBox(G2frame.dataDisplay,-1,label=' Reverse order?')
            reverseSel.Bind(wx.EVT_CHECKBOX,OnReverse)
            reverseSel.SetValue(data['Reverse Seq'])
            selSizer.Add(reverseSel,0,WACV)
            copySel =  wx.CheckBox(G2frame.dataDisplay,-1,label=' Copy results to next histogram?')
            copySel.Bind(wx.EVT_CHECKBOX,OnCopySel)
            copySel.SetValue(data['Copy2Next'])
            selSizer.Add(copySel,0,WACV)
            seqSizer.Add(selSizer,0)
        return seqSizer
        
    def LSSizer():        
        
        def OnDerivType(event):
            data['deriv type'] = derivSel.GetValue()
            derivSel.SetValue(data['deriv type'])
            wx.CallAfter(UpdateControls,G2frame,data)
            
        def OnConvergence(event):
            try:
                value = max(1.e-9,min(1.0,float(Cnvrg.GetValue())))
            except ValueError:
                value = 0.0001
            data['min dM/M'] = value
            Cnvrg.SetValue('%.2g'%(value))
            
        def OnMaxCycles(event):
            data['max cyc'] = int(maxCyc.GetValue())
            maxCyc.SetValue(str(data['max cyc']))
                        
        def OnFactor(event):
            try:
                value = min(max(float(Factr.GetValue()),0.00001),100.)
            except ValueError:
                value = 1.0
            data['shift factor'] = value
            Factr.SetValue('%.5f'%(value))
            
        def OnFsqRef(event):
            data['F**2'] = fsqRef.GetValue()
        
        def OnUsrRej(event):
            Obj = event.GetEventObject()
            item,limits = Indx[Obj]
            try:
                value = min(max(float(Obj.GetValue()),limits[0]),limits[1])
            except ValueError:
                value = data['UsrReject'][item]
            data['UsrReject'][item] = value
            Obj.SetValue('%.2f'%(value))

        LSSizer = wx.FlexGridSizer(cols=4,vgap=5,hgap=5)
        LSSizer.Add(wx.StaticText(G2frame.dataDisplay,label=' Refinement derivatives: '),0,WACV)
        Choice=['analytic Jacobian','numeric','analytic Hessian']
        derivSel = wx.ComboBox(parent=G2frame.dataDisplay,value=data['deriv type'],choices=Choice,
            style=wx.CB_READONLY|wx.CB_DROPDOWN)
        derivSel.SetValue(data['deriv type'])
        derivSel.Bind(wx.EVT_COMBOBOX, OnDerivType)
            
        LSSizer.Add(derivSel,0,WACV)
        LSSizer.Add(wx.StaticText(G2frame.dataDisplay,label=' Min delta-M/M: '),0,WACV)
        Cnvrg = wx.TextCtrl(G2frame.dataDisplay,-1,value='%.2g'%(data['min dM/M']),style=wx.TE_PROCESS_ENTER)
        Cnvrg.Bind(wx.EVT_TEXT_ENTER,OnConvergence)
        Cnvrg.Bind(wx.EVT_KILL_FOCUS,OnConvergence)
        LSSizer.Add(Cnvrg,0,WACV)
        Indx = {}
        if 'Hessian' in data['deriv type']:
            LSSizer.Add(wx.StaticText(G2frame.dataDisplay,label=' Max cycles: '),0,WACV)
            Choice = ['0','1','2','3','5','10','15','20']
            maxCyc = wx.ComboBox(parent=G2frame.dataDisplay,value=str(data['max cyc']),choices=Choice,
                style=wx.CB_READONLY|wx.CB_DROPDOWN)
            maxCyc.SetValue(str(data['max cyc']))
            maxCyc.Bind(wx.EVT_COMBOBOX, OnMaxCycles)
            LSSizer.Add(maxCyc,0,WACV)
        else:
            LSSizer.Add(wx.StaticText(G2frame.dataDisplay,label=' Initial shift factor: '),0,WACV)
            Factr = wx.TextCtrl(G2frame.dataDisplay,-1,value='%.5f'%(data['shift factor']),style=wx.TE_PROCESS_ENTER)
            Factr.Bind(wx.EVT_TEXT_ENTER,OnFactor)
            Factr.Bind(wx.EVT_KILL_FOCUS,OnFactor)
            LSSizer.Add(Factr,0,WACV)
        if G2frame.Sngl:
            userReject = data['UsrReject']
            usrRej = {'minF/sig':[' Min obs/sig (0-5): ',[0,5], ],'MinExt':[' Min extinct. (0-.9): ',[0,.9],],
                'MaxDF/F':[' Max delt-F/sig (3-1000): ',[3.,1000.],],'MaxD':[' Max d-spacing (3-500): ',[3,500],],
                'MinD':[' Min d-spacing (0.1-1.0): ',[0.1,1.0],]}

            fsqRef = wx.CheckBox(G2frame.dataDisplay,-1,label='Refine HKLF as F^2? ')
            fsqRef.SetValue(data['F**2'])
            fsqRef.Bind(wx.EVT_CHECKBOX,OnFsqRef)
            LSSizer.Add(fsqRef,0,WACV)
            LSSizer.Add((1,0),)
            for item in usrRej:
                LSSizer.Add(wx.StaticText(G2frame.dataDisplay,-1,label=usrRej[item][0]),0,WACV)
                usrrej = wx.TextCtrl(G2frame.dataDisplay,-1,value='%.2f'%(userReject[item]),style=wx.TE_PROCESS_ENTER)
                Indx[usrrej] = [item,usrRej[item][1]]
                usrrej.Bind(wx.EVT_TEXT_ENTER,OnUsrRej)
                usrrej.Bind(wx.EVT_KILL_FOCUS,OnUsrRej)
                LSSizer.Add(usrrej,0,WACV)
        return LSSizer
        
# </ Anton Gagin            
    def MargMultSizer():
# additional controls in main Controls data tree entry: multiplicative factor
        def OnKnotsNumC(event):
            data['corrParam E_mu'] = KnotsNumC.GetValue()
            KnotsNumC.SetValue(data['corrParam E_mu'])
#            wx.CallAfter(UpdateControls,G2frame,data)    
        def OnMargKMu(event):
            data['corrParam k_mu'] = MargKMu.GetValue()
            MargKMu.SetValue(data['corrParam k_mu'])
            data['EstimateKMu'] = False
            OptKMu.SetValue(data['EstimateKMu'])
#            wx.CallAfter(UpdateControls,G2frame,data)         
        def OnOptKMu(event):
            data['EstimateKMu'] = OptKMu.GetValue()
            if data['EstimateKMu']:
                data['corrParam k_mu tmp'] = data['corrParam k_mu']                
                MargKMu.SetValue("optimal")
            else:
                data['corrParam k_mu'] = data['corrParam k_mu tmp']
                MargKMu.SetValue(data['corrParam k_mu'])
       
        MargMultSizer = wx.FlexGridSizer(cols=4,vgap=5,hgap=5)      
        MargMultSizer.Add(wx.StaticText(G2frame.dataDisplay,label=' Number of knots E_mu: '),0,WACV)
        KnotsNumC = wx.TextCtrl(G2frame.dataDisplay,-1,value=data['corrParam E_mu'],style=wx.TE_PROCESS_ENTER)
        KnotsNumC.SetValue(data['corrParam E_mu'])
        KnotsNumC.Bind(wx.EVT_TEXT_ENTER, OnKnotsNumC)
        KnotsNumC.Bind(wx.EVT_KILL_FOCUS, OnKnotsNumC)
        MargMultSizer.Add(KnotsNumC,0,WACV)     

        tip1 = STT.SuperToolTip("Set number of knots for the multiplicative correction. In most cases 30 is a reasonable maximum.\
\nIf E_mu=0, the correction is not applied. If several histograms are fitted simultaneously, \
\nseparate the E_mu values by commas, or use a single value for all the histograms.")                                               
        tip1.SetHeader("E_mu")
        tip1.SetTarget(KnotsNumC)
        tip1.SetDrawHeaderLine(True)
        tip1.ApplyStyle("Office 2007 Blue")
        
        MargMultSizer.Add(wx.StaticText(G2frame.dataDisplay,label=' Prior factor k_mu: '),0,WACV)
        MargKMu = wx.TextCtrl(G2frame.dataDisplay,-1,value=data['corrParam k_mu'],style=wx.TE_PROCESS_ENTER)
        MargKMu.SetValue(data['corrParam k_mu'])
        MargKMu.Bind(wx.EVT_TEXT_ENTER, OnMargKMu)
        MargKMu.Bind(wx.EVT_KILL_FOCUS, OnMargKMu)
        MargMultSizer.Add(MargKMu,0,WACV)   
           
        tip2 = STT.SuperToolTip("If you don't know how to choose k_mu, select 'Estimate optimal k_mu?'")
        tip2.SetHeader("Scale parameter k_mu")
        tip2.SetTarget(MargKMu)
        tip2.SetDrawHeaderLine(True)
        tip2.ApplyStyle("Office 2007 Blue")
        
        
        OptKMu = wx.CheckBox(G2frame.dataDisplay,-1,label=' Estimate optimal k_mu?')
        OptKMu.Bind(wx.EVT_CHECKBOX, OnOptKMu)
        OptKMu.SetValue(data['EstimateKMu'])
        MargMultSizer.Add(OptKMu,0,WACV)     
        
        return MargMultSizer
        
    def MargAddSizer():
# additional controls in main Controls data tree entry: additive factor

        def OnKnotsNumB(event):
            data['corrParam E_beta'] = KnotsNumB.GetValue()
            KnotsNumB.SetValue(data['corrParam E_beta'])
#            wx.CallAfter(UpdateControls,G2frame,data)    
        def OnMargKBeta(event):
            data['corrParam k_beta'] = MargKBeta.GetValue()
            MargKBeta.SetValue(data['corrParam k_beta'])
            data['EstimateKBeta'] = False
            OptKBeta.SetValue(data['EstimateKBeta'])            
#            wx.CallAfter(UpdateControls,G2frame,data)         
        def OnOptKBeta(event):
            data['EstimateKBeta'] = OptKBeta.GetValue()
            if data['EstimateKBeta']:
                data['corrParam k_beta tmp'] = data['corrParam k_beta']                
                MargKBeta.SetValue("optimal")
            else:
                data['corrParam k_beta'] = data['corrParam k_beta tmp']
                MargKBeta.SetValue(data['corrParam k_beta'])                
       
       
        MargAddSizer = wx.FlexGridSizer(cols=4,vgap=5,hgap=5)      
        MargAddSizer.Add(wx.StaticText(G2frame.dataDisplay,label=' Number of knots E_beta: '),0,WACV)
        KnotsNumB = wx.TextCtrl(G2frame.dataDisplay,-1,value=data['corrParam E_beta'],style=wx.TE_PROCESS_ENTER)
        KnotsNumB.SetValue(data['corrParam E_beta'])
        KnotsNumB.Bind(wx.EVT_TEXT_ENTER, OnKnotsNumB)
        KnotsNumB.Bind(wx.EVT_KILL_FOCUS, OnKnotsNumB)
        MargAddSizer.Add(KnotsNumB,0,WACV)     

        tip3 = STT.SuperToolTip("Set number of knots for the additive correction. \nIn most cases 30 is a reasonable maximum. \nIf E_beta=0, the additive correction is not applied.")
        tip3.SetHeader("E_beta")
        tip3.SetTarget(KnotsNumB)
        tip3.SetDrawHeaderLine(True)
        tip3.ApplyStyle("Office 2007 Blue")
        
        MargAddSizer.Add(wx.StaticText(G2frame.dataDisplay,label=' Prior factor k_beta: '),0,WACV)
        MargKBeta = wx.TextCtrl(G2frame.dataDisplay,-1,value=data['corrParam k_beta'],style=wx.TE_PROCESS_ENTER)
        MargKBeta.SetValue(data['corrParam k_beta'])
        MargKBeta.Bind(wx.EVT_TEXT_ENTER, OnMargKBeta)
        MargKBeta.Bind(wx.EVT_KILL_FOCUS, OnMargKBeta)
        MargAddSizer.Add(MargKBeta,0,WACV)   

        tip4 = STT.SuperToolTip("If you don't know how to choose k_beta, select 'Estimate optimal k_beta?' ")
        tip4.SetHeader("Scale parameter k_beta")
        tip4.SetTarget(MargKBeta)
        tip4.SetDrawHeaderLine(True)
        tip4.ApplyStyle("Office 2007 Blue")   
   
        OptKBeta = wx.CheckBox(G2frame.dataDisplay,-1,label=' Estimate optimal k_beta?')
        OptKBeta.Bind(wx.EVT_CHECKBOX, OnOptKBeta)
        OptKBeta.SetValue(data['EstimateKBeta'])
        MargAddSizer.Add(OptKBeta,0,WACV)     
        
        return MargAddSizer

    def MargNumBlocksSizer():
# additional controls in main Controls data tree entry: number of blocks s
        def OnBlocksNum(event):
            data['corrParam num blocks s'] = BlocksNum.GetValue()
            BlocksNum.SetValue(data['corrParam num blocks s'])

        MargNumBlocksSizer = wx.FlexGridSizer(cols=4,vgap=5,hgap=5)
        MargNumBlocksSizer.Add(wx.StaticText(G2frame.dataDisplay,label=' Number of blocks s: '),0,WACV)
        BlocksNum = wx.TextCtrl(G2frame.dataDisplay,-1,value=data['corrParam num blocks s'],style=wx.TE_PROCESS_ENTER)
        BlocksNum.SetValue(data['corrParam num blocks s'])
        BlocksNum.Bind(wx.EVT_TEXT_ENTER, OnBlocksNum)
        BlocksNum.Bind(wx.EVT_KILL_FOCUS, OnBlocksNum)
        MargNumBlocksSizer.Add(BlocksNum,0,WACV)

        tip5 = STT.SuperToolTip("To reduce the computational complexity (e.g. one may get \nan out-of-memory error for extremely large histograms) \
\nand speed the calculations use s>1. The fitted x-range \nwill be divided into s independent segments. s=6 normally \
\nworks fine. If s=0 no peak-shape corrections applied.")
        tip5.SetHeader("Number of blocks s")
        tip5.SetTarget(BlocksNum)
        tip5.SetDrawHeaderLine(True)
        tip5.ApplyStyle("Office 2007 Blue")   
        
        return MargNumBlocksSizer   

        
    def MargShapeSizer():
# additional controls in main Controls data tree entry: peak-shape factor
        def OnSigDel(event):
            data['corrParam sigma_delta'] = SigDel.GetValue()
            SigDel.SetValue(data['corrParam sigma_delta'])
            LDelDiv.SetValue("none")            
            data['corrParam l_deltaDivN'] = "none"
            
        def OnLDelDiv(event):
            data['corrParam l_deltaDivN'] = LDelDiv.GetValue()
            LDelDiv.SetValue(data['corrParam l_deltaDivN'])
            SigDel.SetValue("none") 
            data['corrParam sigma_delta'] = "none"
            
#            wx.CallAfter(UpdateControls,G2frame,data)     
        def OnLDel(event):
            data['corrParam l_delta'] = LDel.GetValue()
            LDel.SetValue(data['corrParam l_delta'])
            FWHMDiv.SetValue("none")
            data['corrParam FWHMDivN'] = "none"
#            wx.CallAfter(UpdateControls,G2frame,data)     
        
        def OnFWHMDiv(event):
            data['corrParam FWHMDivN'] = FWHMDiv.GetValue()
            FWHMDiv.SetValue(data['corrParam FWHMDivN'])
            LDel.SetValue("none")            
            data['corrParam l_delta'] = "none"
            
        MargShapeSizer = wx.FlexGridSizer(cols=4,vgap=5,hgap=5)
                
        MargShapeSizer.Add(wx.StaticText(G2frame.dataDisplay,label=' Correlation legth l_delta: '),0,WACV)
        LDel = wx.TextCtrl(G2frame.dataDisplay,-1,value=data['corrParam l_delta'],style=wx.TE_PROCESS_ENTER)
        LDel.SetValue(data['corrParam l_delta'])
        LDel.Bind(wx.EVT_TEXT_ENTER, OnLDel)
        LDel.Bind(wx.EVT_KILL_FOCUS, OnLDel)
        MargShapeSizer.Add(LDel,0,WACV)
        
        MargShapeSizer.Add(wx.StaticText(G2frame.dataDisplay,label=' OR estimate it as FWHM /'),0,WACV)
        FWHMDiv = wx.TextCtrl(G2frame.dataDisplay,-1,value=data['corrParam FWHMDivN'],style=wx.TE_PROCESS_ENTER)
        FWHMDiv.SetValue(data['corrParam FWHMDivN'])
        FWHMDiv.Bind(wx.EVT_TEXT_ENTER, OnFWHMDiv)
        FWHMDiv.Bind(wx.EVT_KILL_FOCUS, OnFWHMDiv)
        MargShapeSizer.Add(FWHMDiv,0,WACV)
        
        tip6 = STT.SuperToolTip("Set n1~1-6. Smaller n1 values result in smoother peak-shape corrections. \nPlease refer to the README.pdf for more tips on how to select l_delta")
        tip6.SetHeader("Estimate l_delta as FWHM/n1")
        tip6.SetTarget(FWHMDiv)
        tip6.SetDrawHeaderLine(True)
        tip6.ApplyStyle("Office 2007 Blue")   
        
        MargShapeSizer.Add(wx.StaticText(G2frame.dataDisplay,label=' Stdev sigma_delta: '),0,WACV)
        SigDel = wx.TextCtrl(G2frame.dataDisplay,-1,value=data['corrParam sigma_delta'],style=wx.TE_PROCESS_ENTER)
        SigDel.SetValue(data['corrParam sigma_delta'])
        SigDel.Bind(wx.EVT_TEXT_ENTER, OnSigDel)
        SigDel.Bind(wx.EVT_KILL_FOCUS, OnSigDel)
        MargShapeSizer.Add(SigDel,0,WACV)

        MargShapeSizer.Add(wx.StaticText(G2frame.dataDisplay,label=' OR estimate it as l_delta /'),0,WACV)
        LDelDiv = wx.TextCtrl(G2frame.dataDisplay,-1,value=data['corrParam l_deltaDivN'],style=wx.TE_PROCESS_ENTER)
        LDelDiv.SetValue(data['corrParam l_deltaDivN'])
        LDelDiv.Bind(wx.EVT_TEXT_ENTER, OnLDelDiv)
        LDelDiv.Bind(wx.EVT_KILL_FOCUS, OnLDelDiv)
        MargShapeSizer.Add(LDelDiv,0,WACV)

        tip7 = STT.SuperToolTip("Set n2~1-3. Smaller n2 values allow for larger peak-shape corrections.")
        tip7.SetHeader("Estimate sigma_delta as l_delta/n2")
        tip7.SetTarget(LDelDiv)
        tip7.SetDrawHeaderLine(True)
        tip7.ApplyStyle("Office 2007 Blue")   
        
        return MargShapeSizer    


    def MargIterSizer():
# additional controls in main Controls data tree entry: iterative procedure
          
        def OnDoIter(event):
            data['doIter'] = DoIter.GetValue()
            DoIter.SetValue(data['doIter'])
            
        MargIterSizer = wx.FlexGridSizer(cols=4,vgap=5,hgap=5)
                  
        DoIter = wx.CheckBox(G2frame.dataDisplay,-1,label=' Iterative (may cause overfitting)')
        DoIter.Bind(wx.EVT_CHECKBOX, OnDoIter)
        DoIter.SetValue(data['doIter'])
        MargIterSizer.Add(DoIter,0,WACV)  
        
        return MargIterSizer     


    def MargMCMCSizer():
# MCMC controls   
        def OnDoMCMC(event):
            data['doMCMC'] = doMCMC.GetValue()
        
        MargMCMCSizer = wx.FlexGridSizer(cols=4,vgap=5,hgap=5)        
        doMCMC = wx.CheckBox(G2frame.dataDisplay,-1,label=' run sampler for MCMC?')
        doMCMC.Bind(wx.EVT_CHECKBOX, OnDoMCMC)
        doMCMC.SetValue(data['doMCMC'])
        MargMCMCSizer.Add(doMCMC,0,WACV)      
                
        return MargMCMCSizer
        
        
    def MargMCMCControlsSizer():
# MCMC controls 2
        def OnNWalkers(event):
            data['nwalkers'] = NWalkers.GetValue()
            NWalkers.SetValue(data['nwalkers'])  
        def OnNIterMCMC(event):
            data['nIterMCMC'] = nIterMCMC.GetValue()
            nIterMCMC.SetValue(data['nIterMCMC'])     
        
        MargMCMCControlsSizer = wx.FlexGridSizer(cols=4,vgap=5,hgap=5)                
        MargMCMCControlsSizer.Add(wx.StaticText(G2frame.dataDisplay,label=' Number of Goodman&Weare walkers: '),0,WACV)
        NWalkers = wx.TextCtrl(G2frame.dataDisplay,-1,value=data['nwalkers'],style=wx.TE_PROCESS_ENTER)
        NWalkers.SetValue(data['nwalkers'])
        NWalkers.Bind(wx.EVT_TEXT_ENTER, OnNWalkers)
        NWalkers.Bind(wx.EVT_KILL_FOCUS, OnNWalkers)
        MargMCMCControlsSizer.Add(NWalkers,0,WACV)     

        tip8 = STT.SuperToolTip("Walkers are members of an ensemble (resemble separate Metropolis-Hastings chains). \nSet nWalkers at least nWalkers>20 and nWalkers>4*nFittedParameters.")
        tip8.SetHeader("Number of walkers")
        tip8.SetTarget(NWalkers)
        tip8.SetDrawHeaderLine(True)
        tip8.ApplyStyle("Office 2007 Blue")  
        
        MargMCMCControlsSizer.Add(wx.StaticText(G2frame.dataDisplay,label=' Number of steps: '),0,WACV)
        nIterMCMC = wx.TextCtrl(G2frame.dataDisplay,-1,value=data['nIterMCMC'],style=wx.TE_PROCESS_ENTER)
        nIterMCMC.SetValue(data['nIterMCMC'])
        nIterMCMC.Bind(wx.EVT_TEXT_ENTER, OnNIterMCMC)
        nIterMCMC.Bind(wx.EVT_KILL_FOCUS, OnNIterMCMC)
        MargMCMCControlsSizer.Add(nIterMCMC,0,WACV)   

        tip9 = STT.SuperToolTip("The number of steps to run. The final sample size will be nWalkers*nIterations. \nMake sure it is large enough to sample your nFittedParameters-dimensional distribution.")
        tip9.SetHeader("Number of iterations")
        tip9.SetTarget(nIterMCMC)
        tip9.SetDrawHeaderLine(True)
        tip9.ApplyStyle("Office 2007 Blue")  
        
        return MargMCMCControlsSizer
        
# Anton Gagin />       
        
    def AuthSizer():

        def OnAuthor(event):
            data['Author'] = auth.GetValue()

        Author = data['Author']
        authSizer = wx.BoxSizer(wx.HORIZONTAL)
        authSizer.Add(wx.StaticText(G2frame.dataDisplay,label=' CIF Author (last, first):'),0,WACV)
        auth = wx.TextCtrl(G2frame.dataDisplay,-1,value=Author,style=wx.TE_PROCESS_ENTER)
        auth.Bind(wx.EVT_TEXT_ENTER,OnAuthor)
        auth.Bind(wx.EVT_KILL_FOCUS,OnAuthor)
        authSizer.Add(auth,0,WACV)
        return authSizer
        
        
    if G2frame.dataDisplay:
        G2frame.dataDisplay.Destroy()
    if not G2frame.dataFrame.GetStatusBar():
        Status = G2frame.dataFrame.CreateStatusBar()
        Status.SetStatusText('')
    G2frame.dataFrame.SetLabel('Controls')
    G2frame.dataDisplay = wx.Panel(G2frame.dataFrame)
    SetDataMenuBar(G2frame,G2frame.dataFrame.ControlsMenu)
    mainSizer = wx.BoxSizer(wx.VERTICAL)
    mainSizer.Add((5,5),0)
    mainSizer.Add(wx.StaticText(G2frame.dataDisplay,label=' Refinement Controls:'),0,WACV)    
    mainSizer.Add(LSSizer())
    mainSizer.Add((5,5),0)
    mainSizer.Add(SeqSizer())
    mainSizer.Add((5,5),0)
    mainSizer.Add(AuthSizer())
    mainSizer.Add((5,5),0)
        

# </ Anton Gagin               
   # additional controls in main Controls data tree entry 
    mainSizer.Add((10,10),0)
    mainSizer.Add(wx.StaticLine(G2frame.dataDisplay, -1, (1, 1), (500,1)))
    mainSizer.Add((10,10),0)    
    mainSizer.Add(wx.StaticText(G2frame.dataDisplay,label=' Systematic Error Marginalization Controls:'),0,WACV)    
    mainSizer.Add((10,10),0)
    mainSizer.Add(wx.StaticText(G2frame.dataDisplay,label=' Multiplicative Correction:'),0,WACV)    
    mainSizer.Add(MargMultSizer())
    mainSizer.Add((10,10),0)
    mainSizer.Add(wx.StaticText(G2frame.dataDisplay,label=' Additive Correction:'),0,WACV)    
    mainSizer.Add(MargAddSizer())
    mainSizer.Add((10,10),0)    
    mainSizer.Add(wx.StaticText(G2frame.dataDisplay,label=' Profile Peak-Shape Correction:'),0,WACV)    
    mainSizer.Add(MargNumBlocksSizer())
    mainSizer.Add((5,5),0)    
    mainSizer.Add(MargShapeSizer())
    mainSizer.Add((5,5),0)    
    mainSizer.Add(MargIterSizer())
    mainSizer.Add((10,10),0)
    mainSizer.Add(wx.StaticText(G2frame.dataDisplay,label=' MCMC ensemble sampler:'),0,WACV)    
    mainSizer.Add(MargMCMCSizer())
    mainSizer.Add((5,5),0) 
    mainSizer.Add(MargMCMCControlsSizer())
    mainSizer.Add((5,5),0) 
    
    
# Anton Gagin />    
    mainSizer.Layout()    
    G2frame.dataDisplay.SetSizer(mainSizer)
    G2frame.dataDisplay.SetSize(mainSizer.Fit(G2frame.dataFrame))
    G2frame.dataFrame.setSizePosLeft(mainSizer.Fit(G2frame.dataFrame))
     
################################################################################
#####  Comments
################################################################################           
       
def UpdateComments(G2frame,data):                   

    if G2frame.dataDisplay:
        G2frame.dataDisplay.Destroy()
    G2frame.dataFrame.SetLabel('Comments')
    G2frame.dataDisplay = wx.TextCtrl(parent=G2frame.dataFrame,size=G2frame.dataFrame.GetClientSize(),
        style=wx.TE_MULTILINE|wx.TE_READONLY|wx.TE_DONTWRAP)
    for line in data:
        G2frame.dataDisplay.AppendText(line+'\n')
    G2frame.dataFrame.setSizePosLeft([400,250])
            
################################################################################
#####  Display of Sequential Results
################################################################################           
       
def UpdateSeqResults(G2frame,data,prevSize=None):
    """
    Called when the Sequential Results data tree entry is selected
    to show results from a sequential refinement.
    
    :param wx.Frame G2frame: main GSAS-II data tree windows

    :param dict data: a dictionary containing the following items:  

            * 'histNames' - list of histogram names in order as processed by Sequential Refinement
            * 'varyList' - list of variables - identical over all refinements in sequence
              note that this is the original list of variables, prior to processing
              constraints.
            * 'variableLabels' -- a dict of labels to be applied to each parameter
              (this is created as an empty dict if not present in data).
            * keyed by histName - dictionaries for all data sets processed, which contains:

              * 'variables'- result[0] from leastsq call
              * 'varyList' - list of variables passed to leastsq call (not same as above)
              * 'sig' - esds for variables
              * 'covMatrix' - covariance matrix from individual refinement
              * 'title' - histogram name; same as dict item name
              * 'newAtomDict' - new atom parameters after shifts applied
              * 'newCellDict' - refined cell parameters after shifts to A0-A5 from Dij terms applied'
    """

    def GetSampleParms():
        '''Make a dictionary of the sample parameters are not the same over the
        refinement series.
        '''
        if 'IMG' in histNames[0]:
            sampleParmDict = {'Sample load':[],}
        else:
            sampleParmDict = {'Temperature':[],'Pressure':[],'Time':[],
                'FreePrm1':[],'FreePrm2':[],'FreePrm3':[],'Omega':[],
                'Chi':[],'Phi':[],'Azimuth':[],}
        Controls = G2frame.PatternTree.GetItemPyData(
            GetPatternTreeItemId(G2frame,G2frame.root, 'Controls'))
        sampleParm = {}
        for name in histNames:
            if 'IMG' in name:
                for item in sampleParmDict:
                    sampleParmDict[item].append(data[name]['parmDict'].get(item,0))
            else:
                Id = GetPatternTreeItemId(G2frame,G2frame.root,name)
                sampleData = G2frame.PatternTree.GetItemPyData(GetPatternTreeItemId(G2frame,Id,'Sample Parameters'))
                for item in sampleParmDict:
                    sampleParmDict[item].append(sampleData.get(item,0))
        for item in sampleParmDict:
            frstValue = sampleParmDict[item][0]
            if np.any(np.array(sampleParmDict[item])-frstValue):
                if item.startswith('FreePrm'):
                    sampleParm[Controls[item]] = sampleParmDict[item]
                else:
                    sampleParm[item] = sampleParmDict[item]
        return sampleParm

    def GetColumnInfo(col):
        '''returns column label, lists of values and errors (or None) for each column in the table
        for plotting. The column label is reformatted from Unicode to MatPlotLib encoding
        '''
        colName = G2frame.SeqTable.GetColLabelValue(col)
        plotName = variableLabels.get(colName,colName)
        plotName = plotSpCharFix(plotName)
        return plotName,colList[col],colSigs[col]
            
    def PlotSelect(event):
        'Plots a row (covariance) or column on double-click'
        cols = G2frame.dataDisplay.GetSelectedCols()
        rows = G2frame.dataDisplay.GetSelectedRows()
        if cols:
            G2plt.PlotSelectedSequence(G2frame,cols,GetColumnInfo,SelectXaxis)
        elif rows:
            name = histNames[rows[0]]       #only does 1st one selected
            G2plt.PlotCovariance(G2frame,data[name])
        else:
            G2frame.ErrorDialog(
                'Select row or columns',
                'Nothing selected in table. Click on column or row label(s) to plot. N.B. Grid selection can be a bit funky.'
                )
            
    def OnPlotSelSeq(event):
        'plot the selected columns or row from menu command'
        cols = sorted(G2frame.dataDisplay.GetSelectedCols()) # ignore selection order
        rows = G2frame.dataDisplay.GetSelectedRows()
        if cols:
            G2plt.PlotSelectedSequence(G2frame,cols,GetColumnInfo,SelectXaxis)
        elif rows:
            name = histNames[rows[0]]       #only does 1st one selected
            G2plt.PlotCovariance(G2frame,data[name])
        else:
            G2frame.ErrorDialog(
                'Select columns',
                'No columns or rows selected in table. Click on row or column labels to select fields for plotting.'
                )

                
    def OnAveSelSeq(event):
        'average the selected columns from menu command'
        cols = sorted(G2frame.dataDisplay.GetSelectedCols()) # ignore selection order
        if cols:
            for col in cols:
                ave = np.mean(GetColumnInfo(col)[1])
                sig = np.std(GetColumnInfo(col)[1])
                print ' Average for '+G2frame.SeqTable.GetColLabelValue(col)+': '+'%.6g'%(ave)+' +/- '+'%.6g'%(sig)
        else:
            G2frame.ErrorDialog(
                'Select columns',
                'No columns selected in table. Click on column labels to select fields for averaging.'
                )
                
    def OnRenameSelSeq(event):
        cols = sorted(G2frame.dataDisplay.GetSelectedCols()) # ignore selection order
        colNames = [G2frame.SeqTable.GetColLabelValue(c) for c in cols]
        newNames = colNames[:]
        for i,name in enumerate(colNames):
            if name in variableLabels:
                newNames[i] = variableLabels[name]
        if not cols:
            G2frame.ErrorDialog('Select columns',
                'No columns selected in table. Click on column labels to select fields for rename.')
            return
        dlg = G2G.MultiStringDialog(G2frame.dataDisplay,'Set column names',colNames,newNames)
        if dlg.Show():
            newNames = dlg.GetValues()            
            variableLabels.update(dict(zip(colNames,newNames)))
        data['variableLabels'] = variableLabels 
        dlg.Destroy()
        UpdateSeqResults(G2frame,data,G2frame.dataDisplay.GetSize()) # redisplay variables
        G2plt.PlotSelectedSequence(G2frame,cols,GetColumnInfo,SelectXaxis)
            
    def OnReOrgSelSeq(event):
        'Reorder the columns'
        G2G.GetItemOrder(G2frame,VaryListChanges,vallookup,posdict)    
        UpdateSeqResults(G2frame,data,G2frame.dataDisplay.GetSize()) # redisplay variables

    def OnSaveSelSeqCSV(event):
        'export the selected columns to a .csv file from menu command'
        OnSaveSelSeq(event,csv=True)
        
    def OnSaveSeqCSV(event):
        'export all columns to a .csv file from menu command'
        OnSaveSelSeq(event,csv=True,allcols=True)
        
    def OnSaveSelSeq(event,csv=False,allcols=False):
        'export the selected columns to a .txt or .csv file from menu command'
        def WriteCSV():
            def WriteList(headerItems):
                line = ''
                for lbl in headerItems:
                    if line: line += ','
                    line += '"'+lbl+'"'
                return line
            head = ['name']
            for col in cols:
                item = G2frame.SeqTable.GetColLabelValue(col)
                # get rid of labels that have Unicode characters
                if not all([ord(c) < 128 and ord(c) != 0 for c in item]): item = '?'
                if col in havesig:
                    head += [item,'esd-'+item]
                else:
                    head += [item]
            SeqFile.write(WriteList(head)+'\n')
            for row,name in enumerate(saveNames):
                line = '"'+saveNames[row]+'"'
                for col in cols:
                    if col in havesig:
                        line += ','+str(saveData[col][row])+','+str(saveSigs[col][row])
                    else:
                        line += ','+str(saveData[col][row])
                SeqFile.write(line+'\n')
        def WriteSeq():
            lenName = len(saveNames[0])
            line = '  %s  '%('name'.center(lenName))
            for col in cols:
                item = G2frame.SeqTable.GetColLabelValue(col)
                if col in havesig:
                    line += ' %12s %12s '%(item.center(12),'esd'.center(12))
                else:
                    line += ' %12s '%(item.center(12))
            SeqFile.write(line+'\n')
            for row,name in enumerate(saveNames):
                line = " '%s' "%(saveNames[row])
                for col in cols:
                    if col in havesig:
                        line += ' %12.6f %12.6f '%(saveData[col][row],saveSigs[col][row])
                    else:
                        line += ' %12.6f '%saveData[col][row]
                SeqFile.write(line+'\n')

        # start of OnSaveSelSeq code
        if allcols:
            cols = range(G2frame.SeqTable.GetNumberCols())
        else:
            cols = sorted(G2frame.dataDisplay.GetSelectedCols()) # ignore selection order
        nrows = G2frame.SeqTable.GetNumberRows()
        if not cols:
            G2frame.ErrorDialog('Select columns',
                             'No columns selected in table. Click on column labels to select fields for output.')
            return
        saveNames = [G2frame.SeqTable.GetRowLabelValue(r) for r in range(nrows)]
        saveData = {}
        saveSigs = {}
        havesig = []
        for col in cols:
            name,vals,sigs = GetColumnInfo(col)
            saveData[col] = vals
            if sigs:
                havesig.append(col)
                saveSigs[col] = sigs
        if csv:
            wild = 'CSV output file (*.csv)|*.csv'
        else:
            wild = 'Text output file (*.txt)|*.txt'
        dlg = wx.FileDialog(
            G2frame,
            'Choose text output file for your selection', '.', '', 
            wild,wx.FD_SAVE|wx.FD_OVERWRITE_PROMPT|wx.CHANGE_DIR)
        try:
            if dlg.ShowModal() == wx.ID_OK:
                SeqTextFile = dlg.GetPath()
                SeqTextFile = G2IO.FileDlgFixExt(dlg,SeqTextFile) 
                SeqFile = open(SeqTextFile,'w')
                if csv:
                    WriteCSV()
                else:
                    WriteSeq()
                SeqFile.close()
        finally:
            dlg.Destroy()
                
    def striphist(var,insChar=''):
        'strip a histogram number from a var name'
        sv = var.split(':')
        if len(sv) <= 1: return var
        if sv[1]:
            sv[1] = insChar
        return ':'.join(sv)
        
    def plotSpCharFix(lbl):
        'Change selected unicode characters to their matplotlib equivalent'
        for u,p in [
            (u'\u03B1',r'$\alpha$'),
            (u'\u03B2',r'$\beta$'),
            (u'\u03B3',r'$\gamma$'),
            (u'\u0394\u03C7',r'$\Delta\chi$'),
            ]:
            lbl = lbl.replace(u,p)
        return lbl
    
    def SelectXaxis():
        'returns a selected column number (or None) as the X-axis selection'
        ncols = G2frame.SeqTable.GetNumberCols()
        colNames = [G2frame.SeqTable.GetColLabelValue(r) for r in range(ncols)]
        dlg = G2G.G2SingleChoiceDialog(
            G2frame.dataDisplay,
            'Select x-axis parameter for plot or Cancel for sequence number',
            'Select X-axis',
            colNames)
        try:
            if dlg.ShowModal() == wx.ID_OK:
                col = dlg.GetSelection()
            else:
                col = None
        finally:
            dlg.Destroy()
        return col
    
    def EnablePseudoVarMenus():
        'Enables or disables the PseudoVar menu items that require existing defs'
        if Controls['SeqPseudoVars']:
            val = True
        else:
            val = False
        G2frame.dataFrame.SequentialPvars.Enable(wxDELSEQVAR,val)
        G2frame.dataFrame.SequentialPvars.Enable(wxEDITSEQVAR,val)

    def DelPseudoVar(event):
        'Ask the user to select a pseudo var expression to delete'
        choices = Controls['SeqPseudoVars'].keys()
        selected = G2G.ItemSelector(
            choices,G2frame.dataFrame,
            multiple=True,
            title='Select expressions to remove',
            header='Delete expression')
        if selected is None: return
        for item in selected:
            del Controls['SeqPseudoVars'][choices[item]]
        if selected:
            UpdateSeqResults(G2frame,data,G2frame.dataDisplay.GetSize()) # redisplay variables

    def EditPseudoVar(event):
        'Edit an existing pseudo var expression'
        choices = Controls['SeqPseudoVars'].keys()
        if len(choices) == 1:
            selected = 0
        else:
            selected = G2G.ItemSelector(
                choices,G2frame.dataFrame,
                multiple=False,
                title='Select an expression to edit',
                header='Edit expression')
        if selected is not None:
            dlg = G2exG.ExpressionDialog(
                G2frame.dataDisplay,PSvarDict,
                Controls['SeqPseudoVars'][choices[selected]],
                header="Edit the PseudoVar expression",
                VarLabel="PseudoVar #"+str(selected+1),
                fit=False)
            newobj = dlg.Show(True)
            if newobj:
                calcobj = G2obj.ExpressionCalcObj(newobj)
                del Controls['SeqPseudoVars'][choices[selected]]
                Controls['SeqPseudoVars'][calcobj.eObj.expression] = newobj
                UpdateSeqResults(G2frame,data,G2frame.dataDisplay.GetSize()) # redisplay variables
        
    def AddNewPseudoVar(event):
        'Create a new pseudo var expression'
        dlg = G2exG.ExpressionDialog(
            G2frame.dataDisplay,PSvarDict,
            header='Enter an expression for a PseudoVar here',
            VarLabel = "New PseudoVar",
            fit=False)
        obj = dlg.Show(True)
        dlg.Destroy()
        if obj:
            calcobj = G2obj.ExpressionCalcObj(obj)
            Controls['SeqPseudoVars'][calcobj.eObj.expression] = obj
            UpdateSeqResults(G2frame,data,G2frame.dataDisplay.GetSize()) # redisplay variables

    def UpdateParmDict(parmDict):
        '''generate the atom positions and the direct & reciprocal cell values,
        because they might be needed to evaluate the pseudovar
        '''
        Ddict = dict(zip(['D11','D22','D33','D12','D13','D23'],
                         ['A'+str(i) for i in range(6)])
                     )
        delList = []
        phaselist = []
        for item in parmDict: 
            if ':' not in item: continue
            key = item.split(':')
            if len(key) < 3: continue
            # remove the dA[xyz] terms, they would only bring confusion
            if key[2].startswith('dA'):
                delList.append(item)
            # compute and update the corrected reciprocal cell terms using the Dij values
            elif key[2] in Ddict:
                if key[0] not in phaselist: phaselist.append(key[0])
                akey = key[0]+'::'+Ddict[key[2]]
                parmDict[akey] -= parmDict[item]
                delList.append(item)
        for item in delList:
            del parmDict[item]                
        for i in phaselist:
            pId = int(i)
            # apply cell symmetry
            A,zeros = G2stIO.cellFill(str(pId)+'::',SGdata[pId],parmDict,zeroDict[pId])
            # convert to direct cell & add the unique terms to the dictionary
            for i,val in enumerate(G2lat.A2cell(A)):
                if i in uniqCellIndx[pId]:
                    lbl = str(pId)+'::'+cellUlbl[i]
                    parmDict[lbl] = val
            lbl = str(pId)+'::'+'vol'
            parmDict[lbl] = G2lat.calc_V(A)
        return parmDict

    def EvalPSvarDeriv(calcobj,parmDict,sampleDict,var,ESD):
        '''Evaluate an expression derivative with respect to a
        GSAS-II variable name.

        Note this likely could be faster if the loop over calcobjs were done
        inside after the Dict was created. 
        '''
        step = ESD/10
        Ddict = dict(zip(['D11','D22','D33','D12','D13','D23'],
                         ['A'+str(i) for i in range(6)])
                     )
        results = []
        phaselist = []
        VparmDict = sampleDict.copy()
        for incr in step,-step:
            VparmDict.update(parmDict.copy())           
            # as saved, the parmDict has updated 'A[xyz]' values, but 'dA[xyz]'
            # values are not zeroed: fix that!
            VparmDict.update({item:0.0 for item in parmDict if 'dA' in item})
            VparmDict[var] += incr
            G2mv.Dict2Map(VparmDict,[]) # apply constraints
            # generate the atom positions and the direct & reciprocal cell values now, because they might
            # needed to evaluate the pseudovar
            for item in VparmDict:
                if item in sampleDict:
                    continue 
                if ':' not in item: continue
                key = item.split(':')
                if len(key) < 3: continue
                # apply any new shifts to atom positions
                if key[2].startswith('dA'):
                    VparmDict[''.join(item.split('d'))] += VparmDict[item]
                    VparmDict[item] = 0.0
                # compute and update the corrected reciprocal cell terms using the Dij values
                if key[2] in Ddict:
                    if key[0] not in phaselist: phaselist.append(key[0])
                    akey = key[0]+'::'+Ddict[key[2]]
                    VparmDict[akey] -= VparmDict[item]
            for i in phaselist:
                pId = int(i)
                # apply cell symmetry
                A,zeros = G2stIO.cellFill(str(pId)+'::',SGdata[pId],VparmDict,zeroDict[pId])
                # convert to direct cell & add the unique terms to the dictionary
                for i,val in enumerate(G2lat.A2cell(A)):
                    if i in uniqCellIndx[pId]:
                        lbl = str(pId)+'::'+cellUlbl[i]
                        VparmDict[lbl] = val
                lbl = str(pId)+'::'+'vol'
                VparmDict[lbl] = G2lat.calc_V(A)
            # dict should be fully updated, use it & calculate
            calcobj.SetupCalc(VparmDict)
            results.append(calcobj.EvalExpression())
        return (results[0] - results[1]) / (2.*step)
        
    def EnableParFitEqMenus():
        'Enables or disables the Parametric Fit menu items that require existing defs'
        if Controls['SeqParFitEqList']:
            val = True
        else:
            val = False
        G2frame.dataFrame.SequentialPfit.Enable(wxDELPARFIT,val)
        G2frame.dataFrame.SequentialPfit.Enable(wxEDITPARFIT,val)
        G2frame.dataFrame.SequentialPfit.Enable(wxDOPARFIT,val)

    def ParEqEval(Values,calcObjList,varyList):
        '''Evaluate the parametric expression(s)
        :param list Values: a list of values for each variable parameter
        :param list calcObjList: a list of :class:`GSASIIobj.ExpressionCalcObj`
          expression objects to evaluate
        :param list varyList: a list of variable names for each value in Values
        '''
        result = []
        for calcobj in calcObjList:
            calcobj.UpdateVars(varyList,Values)
            result.append((calcobj.depVal-calcobj.EvalExpression())/calcobj.depSig)
        return result

    def DoParEqFit(event,eqObj=None):
        'Parametric fit minimizer'
        varyValueDict = {} # dict of variables and their initial values
        calcObjList = [] # expression objects, ready to go for each data point
        if eqObj is not None:
            eqObjList = [eqObj,]
        else:
            eqObjList = Controls['SeqParFitEqList']
        UseFlags = G2frame.SeqTable.GetColValues(0)         
        for obj in eqObjList:
            expr = obj.expression
            # assemble refined vars for this equation
            varyValueDict.update({var:val for var,val in obj.GetVariedVarVal()})
            # lookup dependent var position
            depVar = obj.GetDepVar()
            if depVar in colLabels:
                indx = colLabels.index(depVar)
            else:
                raise Exception('Dependent variable '+depVar+' not found')
            # assemble a list of the independent variables
            indepVars = obj.GetIndependentVars()
            # loop over each datapoint
            for j,row in enumerate(zip(*colList)):
                if not UseFlags[j]: continue
                # assemble equations to fit
                calcobj = G2obj.ExpressionCalcObj(obj)
                # prepare a dict of needed independent vars for this expression
                indepVarDict = {var:row[i] for i,var in enumerate(colLabels) if var in indepVars}
                calcobj.SetupCalc(indepVarDict)                
                # values and sigs for current value of dependent var
                calcobj.depVal = row[indx]
                calcobj.depSig = colSigs[indx][j]
                calcObjList.append(calcobj)
        # varied parameters
        varyList = varyValueDict.keys()
        values = varyValues = [varyValueDict[key] for key in varyList]
        if not varyList:
            print 'no variables to refine!'
            return
        try:
            result = so.leastsq(ParEqEval,varyValues,full_output=True,   #ftol=Ftol,
                                args=(calcObjList,varyList)
                                )
            values = result[0]
            covar = result[1]
            if covar is None:
                raise Exception
            esdDict = {}
            for i,avar in enumerate(varyList):
                esdDict[avar] = np.sqrt(covar[i,i])
        except:
            print('====> Fit failed')
            return
        print('==== Fit Results ====')
        for obj in eqObjList:
            obj.UpdateVariedVars(varyList,values)
            ind = '      '
            print('  '+obj.GetDepVar()+' = '+obj.expression)
            for var in obj.assgnVars:
                print(ind+var+' = '+obj.assgnVars[var])
            for var in obj.freeVars:
                avar = "::"+obj.freeVars[var][0]
                val = obj.freeVars[var][1]
                if obj.freeVars[var][2]:
                    print(ind+var+' = '+avar + " = " + G2mth.ValEsd(val,esdDict[avar]))
                else:
                    print(ind+var+' = '+avar + " =" + G2mth.ValEsd(val,0))
        # create a plot for each parametric variable
        for fitnum,obj in enumerate(eqObjList):
            calcobj = G2obj.ExpressionCalcObj(obj)
            # lookup dependent var position
            indx = colLabels.index(obj.GetDepVar())
            # assemble a list of the independent variables
            indepVars = obj.GetIndependentVars()            
            # loop over each datapoint
            fitvals = []
            for j,row in enumerate(zip(*colList)):
                calcobj.SetupCalc(
                    {var:row[i] for i,var in enumerate(colLabels) if var in indepVars}
                    )
                fitvals.append(calcobj.EvalExpression())
            G2plt.PlotSelectedSequence(
                G2frame,[indx],GetColumnInfo,SelectXaxis,
                fitnum,fitvals)

    def SingleParEqFit(eqObj):
        DoParEqFit(None,eqObj)

    def DelParFitEq(event):
        'Ask the user to select function to delete'
        txtlst = [obj.GetDepVar()+' = '+obj.expression for obj in Controls['SeqParFitEqList']]
        selected = G2G.ItemSelector(
            txtlst,G2frame.dataFrame,
            multiple=True,
            title='Select a parametric equation(s) to remove',
            header='Delete equation')
        if selected is None: return
        Controls['SeqParFitEqList'] = [obj for i,obj in enumerate(Controls['SeqParFitEqList']) if i not in selected]
        EnableParFitEqMenus()
        if Controls['SeqParFitEqList']: DoParEqFit(event)
        
    def EditParFitEq(event):
        'Edit an existing parametric equation'
        txtlst = [obj.GetDepVar()+' = '+obj.expression for obj in Controls['SeqParFitEqList']]
        if len(txtlst) == 1:
            selected = 0
        else:
            selected = G2G.ItemSelector(
                txtlst,G2frame.dataFrame,
                multiple=False,
                title='Select a parametric equation to edit',
                header='Edit equation')
        if selected is not None:
            dlg = G2exG.ExpressionDialog(
                G2frame.dataDisplay,indepVarDict,
                Controls['SeqParFitEqList'][selected],
                depVarDict=depVarDict,
                header="Edit the formula for this minimization function",
                ExtraButton=['Fit',SingleParEqFit])
            newobj = dlg.Show(True)
            if newobj:
                calcobj = G2obj.ExpressionCalcObj(newobj)
                Controls['SeqParFitEqList'][selected] = newobj
                EnableParFitEqMenus()
            if Controls['SeqParFitEqList']: DoParEqFit(event)

    def AddNewParFitEq(event):
        'Create a new parametric equation to be fit to sequential results'

        # compile the variable names used in previous freevars to avoid accidental name collisions
        usedvarlist = []
        for obj in Controls['SeqParFitEqList']:
            for var in obj.freeVars:
                if obj.freeVars[var][0] not in usedvarlist: usedvarlist.append(obj.freeVars[var][0])

        dlg = G2exG.ExpressionDialog(
            G2frame.dataDisplay,indepVarDict,
            depVarDict=depVarDict,
            header='Define an equation to minimize in the parametric fit',
            ExtraButton=['Fit',SingleParEqFit],
            usedVars=usedvarlist)
        obj = dlg.Show(True)
        dlg.Destroy()
        if obj:
            Controls['SeqParFitEqList'].append(obj)
            EnableParFitEqMenus()
            if Controls['SeqParFitEqList']: DoParEqFit(event)
                
    def CopyParFitEq(event):
        'Copy an existing parametric equation to be fit to sequential results'
        # compile the variable names used in previous freevars to avoid accidental name collisions
        usedvarlist = []
        for obj in Controls['SeqParFitEqList']:
            for var in obj.freeVars:
                if obj.freeVars[var][0] not in usedvarlist: usedvarlist.append(obj.freeVars[var][0])
        txtlst = [obj.GetDepVar()+' = '+obj.expression for obj in Controls['SeqParFitEqList']]
        if len(txtlst) == 1:
            selected = 0
        else:
            selected = G2G.ItemSelector(
                txtlst,G2frame.dataFrame,
                multiple=False,
                title='Select a parametric equation to copy',
                header='Copy equation')
        if selected is not None:
            newEqn = copy.deepcopy(Controls['SeqParFitEqList'][selected])
            for var in newEqn.freeVars:
                newEqn.freeVars[var][0] = G2obj.MakeUniqueLabel(newEqn.freeVars[var][0],usedvarlist)
            dlg = G2exG.ExpressionDialog(
                G2frame.dataDisplay,indepVarDict,
                newEqn,
                depVarDict=depVarDict,
                header="Edit the formula for this minimization function",
                ExtraButton=['Fit',SingleParEqFit])
            newobj = dlg.Show(True)
            if newobj:
                calcobj = G2obj.ExpressionCalcObj(newobj)
                Controls['SeqParFitEqList'].append(newobj)
                EnableParFitEqMenus()
            if Controls['SeqParFitEqList']: DoParEqFit(event)
                                            
    def GridSetToolTip(row,col):
        '''Routine to show standard uncertainties for each element in table
        as a tooltip
        '''
        if colSigs[col]:
            return u'\u03c3 = '+str(colSigs[col][row])
        return ''
        
    def GridColLblToolTip(col):
        '''Define a tooltip for a column. This will be the user-entered value
        (from data['variableLabels']) or the default name
        '''
        if col < 0 or col > len(colLabels):
            print 'Illegal column #',col
            return
        var = colLabels[col]
        return variableLabels.get(var,G2obj.fmtVarDescr(var))
        
    def SetLabelString(event):
        '''Define or edit the label for a column in the table, to be used
        as a tooltip and for plotting
        '''
        col = event.GetCol()
        if col < 0 or col > len(colLabels):
            return
        var = colLabels[col]
        lbl = variableLabels.get(var,G2obj.fmtVarDescr(var))
        dlg = G2G.SingleStringDialog(G2frame.dataFrame,'Set variable label',
                                 'Set a new name for variable '+var,lbl,size=(400,-1))
        if dlg.Show():
            variableLabels[var] = dlg.GetValue()
        dlg.Destroy()
        
    #def GridRowLblToolTip(row): return 'Row ='+str(row)
    
    # lookup table for unique cell parameters by symmetry
    cellGUIlist = [
        [['m3','m3m'],(0,)],
        [['3R','3mR'],(0,3)],
        [['3','3m1','31m','6/m','6/mmm','4/m','4/mmm'],(0,2)],
        [['mmm'],(0,1,2)],
        [['2/m'+'a'],(0,1,2,3)],
        [['2/m'+'b'],(0,1,2,4)],
        [['2/m'+'c'],(0,1,2,5)],
        [['-1'],(0,1,2,3,4,5)],
        ]
    # cell labels
    cellUlbl = ('a','b','c',u'\u03B1',u'\u03B2',u'\u03B3') # unicode a,b,c,alpha,beta,gamma

    #======================================================================
    # start processing sequential results here (UpdateSeqResults)
    #======================================================================
    if not data:
        print 'No sequential refinement results'
        return
    variableLabels = data.get('variableLabels',{})
    data['variableLabels'] = variableLabels
    Histograms,Phases = G2frame.GetUsedHistogramsAndPhasesfromTree()
    Controls = G2frame.PatternTree.GetItemPyData(GetPatternTreeItemId(G2frame,G2frame.root,'Controls'))
    # create a place to store Pseudo Vars & Parametric Fit functions, if not present
    if 'SeqPseudoVars' not in Controls: Controls['SeqPseudoVars'] = {}
    if 'SeqParFitEqList' not in Controls: Controls['SeqParFitEqList'] = []
    histNames = data['histNames']
    if G2frame.dataDisplay:
        G2frame.dataDisplay.Destroy()
    if not G2frame.dataFrame.GetStatusBar():
        Status = G2frame.dataFrame.CreateStatusBar()
        Status.SetStatusText("Select column to export; Double click on column to plot data; on row for Covariance")
    sampleParms = GetSampleParms()

    # make dict of varied atom coords keyed by absolute position
    newAtomDict = data[histNames[0]].get('newAtomDict',{}) # dict with atom positions; relative & absolute
    # Possible error: the next might need to be data[histNames[0]]['varyList']
    # error will arise if there constraints on coordinates?
    atomLookup = {newAtomDict[item][0]:item for item in newAtomDict if item in data['varyList']}
    
    # make dict of varied cell parameters equivalents
    ESDlookup = {} # provides the Dij term for each Ak term (where terms are refined)
    Dlookup = {} # provides the Ak term for each Dij term (where terms are refined)
    # N.B. These Dij vars are missing a histogram #
    newCellDict = data[histNames[0]].get('newCellDict',{})
    for item in newCellDict:
        if item in data['varyList']:
            ESDlookup[newCellDict[item][0]] = item
            Dlookup[item] = newCellDict[item][0]
    # add coordinate equivalents to lookup table
    for parm in atomLookup:
        Dlookup[atomLookup[parm]] = parm
        ESDlookup[parm] = atomLookup[parm]

    # get unit cell & symmetry for all phases & initial stuff for later use
    RecpCellTerms = {}
    SGdata = {}
    uniqCellIndx = {}
    initialCell = {}
    RcellLbls = {}
    zeroDict = {}
    Rcelldict = {}
    for phase in Phases:
        phasedict = Phases[phase]
        pId = phasedict['pId']
        pfx = str(pId)+'::' # prefix for A values from phase
        RcellLbls[pId] = [pfx+'A'+str(i) for i in range(6)]
        RecpCellTerms[pId] = G2lat.cell2A(phasedict['General']['Cell'][1:7])
        zeroDict[pId] = dict(zip(RcellLbls[pId],6*[0.,]))
        SGdata[pId] = phasedict['General']['SGData']
        Rcelldict.update({lbl:val for lbl,val in zip(RcellLbls[pId],RecpCellTerms[pId])})
        laue = SGdata[pId]['SGLaue']
        if laue == '2/m':
            laue += SGdata[pId]['SGUniq']
        for symlist,celllist in cellGUIlist:
            if laue in symlist:
                uniqCellIndx[pId] = celllist
                break
        else: # should not happen
            uniqCellIndx[pId] = range(6)
        for i in uniqCellIndx[pId]:
            initialCell[str(pId)+'::A'+str(i)] =  RecpCellTerms[pId][i]

    SetDataMenuBar(G2frame,G2frame.dataFrame.SequentialMenu)
    G2frame.dataFrame.SetLabel('Sequential refinement results')
    if not G2frame.dataFrame.GetStatusBar():
        Status = G2frame.dataFrame.CreateStatusBar()
        Status.SetStatusText('')
    G2frame.dataFrame.Bind(wx.EVT_MENU, OnRenameSelSeq, id=wxID_RENAMESEQSEL)
    G2frame.dataFrame.Bind(wx.EVT_MENU, OnSaveSelSeq, id=wxID_SAVESEQSEL)
    G2frame.dataFrame.Bind(wx.EVT_MENU, OnSaveSelSeqCSV, id=wxID_SAVESEQSELCSV)
    G2frame.dataFrame.Bind(wx.EVT_MENU, OnSaveSeqCSV, id=wxID_SAVESEQCSV)
    G2frame.dataFrame.Bind(wx.EVT_MENU, OnPlotSelSeq, id=wxID_PLOTSEQSEL)
    G2frame.dataFrame.Bind(wx.EVT_MENU, OnAveSelSeq, id=wxID_AVESEQSEL)
    G2frame.dataFrame.Bind(wx.EVT_MENU, OnReOrgSelSeq, id=wxID_ORGSEQSEL)
    G2frame.dataFrame.Bind(wx.EVT_MENU, AddNewPseudoVar, id=wxADDSEQVAR)
    G2frame.dataFrame.Bind(wx.EVT_MENU, DelPseudoVar, id=wxDELSEQVAR)
    G2frame.dataFrame.Bind(wx.EVT_MENU, EditPseudoVar, id=wxEDITSEQVAR)
    G2frame.dataFrame.Bind(wx.EVT_MENU, AddNewParFitEq, id=wxADDPARFIT)
    G2frame.dataFrame.Bind(wx.EVT_MENU, CopyParFitEq, id=wxCOPYPARFIT)
    G2frame.dataFrame.Bind(wx.EVT_MENU, DelParFitEq, id=wxDELPARFIT)
    G2frame.dataFrame.Bind(wx.EVT_MENU, EditParFitEq, id=wxEDITPARFIT)
    G2frame.dataFrame.Bind(wx.EVT_MENU, DoParEqFit, id=wxDOPARFIT)
    EnablePseudoVarMenus()
    EnableParFitEqMenus()

    # scan for locations where the variables change
    VaryListChanges = [] # histograms where there is a change
    combinedVaryList = []
    firstValueDict = {}
    vallookup = {}
    posdict = {}
    prevVaryList = []
    for i,name in enumerate(histNames):
        for var,val,sig in zip(data[name]['varyList'],data[name]['variables'],data[name]['sig']):
            svar = striphist(var,'*') # wild-carded
            if svar not in combinedVaryList:
                # add variables to list as they appear
                combinedVaryList.append(svar)
                firstValueDict[svar] = (val,sig)
        if prevVaryList != data[name]['varyList']: # this refinement has a different refinement list from previous
            prevVaryList = data[name]['varyList']
            vallookup[name] = dict(zip(data[name]['varyList'],data[name]['variables']))
            posdict[name] = {}
            for var in data[name]['varyList']:
                svar = striphist(var,'*')
                posdict[name][combinedVaryList.index(svar)] = svar
            VaryListChanges.append(name)
    if len(VaryListChanges) > 1:
        G2frame.dataFrame.SequentialFile.Enable(wxID_ORGSEQSEL,True)
    else:
        G2frame.dataFrame.SequentialFile.Enable(wxID_ORGSEQSEL,False)
    #-----------------------------------------------------------------------------------
    # build up the data table by columns -----------------------------------------------
    nRows = len(histNames)
    colList = [nRows*[True]]
    colSigs = [None]
    colLabels = ['Use']
    Types = [wg.GRID_VALUE_BOOL]
    # start with Rwp values
    if 'IMG ' not in histNames[0][:4]:
        colList += [[data[name]['Rvals']['Rwp'] for name in histNames]]
        colSigs += [None]
        colLabels += ['Rwp']
        Types += [wg.GRID_VALUE_FLOAT+':10,3',]
    # add % change in Chi^2 in last cycle
    if histNames[0][:4] not in ['SASD','IMG '] and Controls.get('ShowCell'):
        colList += [[100.*data[name]['Rvals'].get('DelChi2',-1) for name in histNames]]
        colSigs += [None]
        colLabels += [u'\u0394\u03C7\u00B2 (%)']
        Types += [wg.GRID_VALUE_FLOAT,]
    deltaChiCol = len(colLabels)-1
    # add changing sample parameters to table
    for key in sampleParms:
        colList += [sampleParms[key]]
        colSigs += [None]
        colLabels += [key]
        Types += [wg.GRID_VALUE_FLOAT,]
    sampleDict = {}
    for i,name in enumerate(histNames):
        sampleDict[name] = dict(zip(sampleParms.keys(),[sampleParms[key][i] for key in sampleParms.keys()])) 
    # add unique cell parameters TODO: review this where the cell symmetry changes (when possible)
    if Controls.get('ShowCell',False):
        for pId in sorted(RecpCellTerms):
            pfx = str(pId)+'::' # prefix for A values from phase
            cells = []
            cellESDs = []
            colLabels += [pfx+cellUlbl[i] for i in uniqCellIndx[pId]]
            colLabels += [pfx+'Vol']
            Types += (1+len(uniqCellIndx[pId]))*[wg.GRID_VALUE_FLOAT,]
            for name in histNames:
                covData = {
                    'varyList': [Dlookup.get(striphist(v),v) for v in data[name]['varyList']],
                    'covMatrix': data[name]['covMatrix']
                    }
                A = RecpCellTerms[pId][:] # make copy of starting A values
                # update with refined values
                for i in range(6):
                    var = str(pId)+'::A'+str(i)
                    if var in ESDlookup:
                        val = data[name]['newCellDict'][ESDlookup[var]][1] # get refined value 
                        A[i] = val # override with updated value
                # apply symmetry
                Albls = [pfx+'A'+str(i) for i in range(6)]
                cellDict = dict(zip(Albls,A))
                A,zeros = G2stIO.cellFill(pfx,SGdata[pId],cellDict,zeroDict[pId])
                # convert to direct cell & add only unique values to table
                c = G2lat.A2cell(A)
                vol = G2lat.calc_V(A)
                cE = G2stIO.getCellEsd(pfx,SGdata[pId],A,covData)
                cells += [[c[i] for i in uniqCellIndx[pId]]+[vol]]
                cellESDs += [[cE[i] for i in uniqCellIndx[pId]]+[cE[-1]]]
            colList += zip(*cells)
            colSigs += zip(*cellESDs)
    # sort out the variables in their selected order
    varcols = 0
    for d in posdict.itervalues():
        varcols = max(varcols,max(d.keys())+1)
    # get labels for each column
    for i in range(varcols):
        lbl = ''
        for h in VaryListChanges:
            if posdict[h].get(i):
                if posdict[h].get(i) in lbl: continue
                if lbl != "": lbl += '/'
                lbl += posdict[h].get(i)
        colLabels.append(lbl)
    Types += varcols*[wg.GRID_VALUE_FLOAT]
    vals = []
    esds = []
    varsellist = None        # will be a list of variable names in the order they are selected to appear
    # tabulate values for each hist, leaving None for blank columns
    for name in histNames:
        if name in posdict:
            varsellist = [posdict[name].get(i) for i in range(varcols)]
            # translate variable names to how they will be used in the headings
            vs = [striphist(v,'*') for v in data[name]['varyList']]
            # determine the index for each column (or None) in the data[]['variables'] and ['sig'] lists
            sellist = [vs.index(v) if v is not None else None for v in varsellist]
            #sellist = [i if striphist(v,'*') in varsellist else None for i,v in enumerate(data[name]['varyList'])]
        if not varsellist: raise Exception()
        vals.append([data[name]['variables'][s] if s is not None else None for s in sellist])
        esds.append([data[name]['sig'][s] if s is not None else None for s in sellist])
        #GSASIIpath.IPyBreak()
    colList += zip(*vals)
    colSigs += zip(*esds)
                
    # tabulate constrained variables, removing histogram numbers if needed
    # from parameter label
    depValDict = {}
    depSigDict = {}
    for name in histNames:
        for var in data[name].get('depParmDict',{}):
            val,sig = data[name]['depParmDict'][var]
            svar = striphist(var,'*')
            if svar not in depValDict:
               depValDict[svar] = [val]
               depSigDict[svar] = [sig]
            else:
               depValDict[svar].append(val)
               depSigDict[svar].append(sig)
    # add the dependent constrained variables to the table
    for var in sorted(depValDict):
        if len(depValDict[var]) != len(histNames): continue
        colLabels.append(var)
        Types += [wg.GRID_VALUE_FLOAT,]
        colSigs += [depSigDict[var]]
        colList += [depValDict[var]]

    # add atom parameters to table
    colLabels += atomLookup.keys()
    Types += len(atomLookup)*[wg.GRID_VALUE_FLOAT]
    for parm in sorted(atomLookup):
        colList += [[data[name]['newAtomDict'][atomLookup[parm]][1] for name in histNames]]
        if atomLookup[parm] in data[histNames[0]]['varyList']:
            col = data[histNames[0]]['varyList'].index(atomLookup[parm])
            colSigs += [[data[name]['sig'][col] for name in histNames]]
        else:
            colSigs += [None] # should not happen
    # evaluate Pseudovars, their ESDs and add them to grid
    for expr in Controls['SeqPseudoVars']:
        obj = Controls['SeqPseudoVars'][expr]
        calcobj = G2obj.ExpressionCalcObj(obj)
        valList = []
        esdList = []
        for seqnum,name in enumerate(histNames):
            sigs = data[name]['sig']
            G2mv.InitVars()
            parmDict = data[name].get('parmDict')
            badVary = data[name].get('badVary',[])
            constraintInfo = data[name].get('constraintInfo',[[],[],{},[],seqnum])
            groups,parmlist,constrDict,fixedList,ihst = constraintInfo
            varyList = data[name]['varyList']
            parmDict = data[name]['parmDict']
            G2mv.GenerateConstraints(groups,parmlist,varyList,constrDict,fixedList,parmDict,SeqHist=ihst)
            derivs = np.array(
                [EvalPSvarDeriv(calcobj,parmDict.copy(),sampleDict[name],var,ESD)
                 for var,ESD in zip(varyList,sigs)]
                )
            esdList.append(np.sqrt(
                np.inner(derivs,np.inner(data[name]['covMatrix'],derivs.T))
                ))
            PSvarDict = parmDict.copy()
            PSvarDict.update(sampleDict[name])
            UpdateParmDict(PSvarDict)
            calcobj.UpdateDict(PSvarDict)
            valList.append(calcobj.EvalExpression())
        if not esdList:
            esdList = None
        colList += [valList]
        colSigs += [esdList]
        colLabels += [expr]
        Types += [wg.GRID_VALUE_FLOAT,]
    #---- table build done -------------------------------------------------------------

    # Make dict needed for creating & editing pseudovars (PSvarDict).
    name = histNames[0]
    parmDict = data[name].get('parmDict')
    PSvarDict = parmDict.copy()
    PSvarDict.update(sampleParms)
    UpdateParmDict(PSvarDict)
    # Also dicts of dependent (depVarDict) & independent vars (indepVarDict)
    # for Parametric fitting from the data table
    parmDict = dict(zip(colLabels,zip(*colList)[0])) # scratch dict w/all values in table
    parmDict.update(
        {var:val for var,val in data[name].get('newCellDict',{}).values()} #  add varied reciprocal cell terms
    )
    name = histNames[0]

    #******************************************************************************
    # create a set of values for example evaluation of pseudovars and 
    # this does not work for refinements that have differing numbers of variables.
    #raise Exception
    indepVarDict = {}     #  values in table w/o ESDs
    depVarDict = {}
    for i,var in enumerate(colLabels):
        if var == 'Use': continue
        if colList[i][0] is None:
            val,sig = firstValueDict.get(var,[None,None])
        elif colSigs[i]:
            val,sig = colList[i][0],colSigs[i][0]
        else:
            val,sig = colList[i][0],None
        if val is None:
            continue
        elif sig is None:
            indepVarDict[var] = val
        elif striphist(var) not in Dlookup:
            depVarDict[var] = val
    # add recip cell coeff. values
    depVarDict.update({var:val for var,val in data[name].get('newCellDict',{}).values()})

    G2frame.dataDisplay = G2G.GSGrid(parent=G2frame.dataFrame)
    G2frame.SeqTable = G2G.Table(
        [list(c) for c in zip(*colList)],     # convert from columns to rows
        colLabels=colLabels,rowLabels=histNames,types=Types)
    G2frame.dataDisplay.SetTable(G2frame.SeqTable, True)
    #G2frame.dataDisplay.EnableEditing(False)
    # make all but first column read-only
    for c in range(1,len(colLabels)):
        for r in range(nRows):
            G2frame.dataDisplay.SetCellReadOnly(r,c)
    G2frame.dataDisplay.Bind(wg.EVT_GRID_LABEL_LEFT_DCLICK, PlotSelect)
    G2frame.dataDisplay.Bind(wg.EVT_GRID_LABEL_RIGHT_CLICK, SetLabelString)
    G2frame.dataDisplay.SetRowLabelSize(8*len(histNames[0]))       #pretty arbitrary 8
    G2frame.dataDisplay.SetMargins(0,0)
    G2frame.dataDisplay.AutoSizeColumns(True)
    if prevSize:
        G2frame.dataDisplay.SetSize(prevSize)
    else:
        G2frame.dataFrame.setSizePosLeft([700,350])
    # highlight unconverged shifts 
    if histNames[0][:4] not in ['SASD','IMG ']:
        for row,name in enumerate(histNames):
            deltaChi = G2frame.SeqTable.GetValue(row,deltaChiCol)
            if deltaChi > 10.:
                G2frame.dataDisplay.SetCellStyle(row,deltaChiCol,color=wx.Colour(255,0,0))
            elif deltaChi > 1.0:
                G2frame.dataDisplay.SetCellStyle(row,deltaChiCol,color=wx.Colour(255,255,0))
    G2frame.dataDisplay.InstallGridToolTip(GridSetToolTip,GridColLblToolTip)
    G2frame.dataDisplay.SendSizeEvent() # resize needed on mac
    G2frame.dataDisplay.Refresh() # shows colored text on mac
    
################################################################################
#####  Main PWDR panel
################################################################################           
       
def UpdatePWHKPlot(G2frame,kind,item):
    '''Called when the histogram main tree entry is called. Displays the
    histogram weight factor, refinement statistics for the histogram
    and the range of data for a simulation.

    Also invokes a plot of the histogram.
    '''
    def onEditSimRange(event):
        'Edit simulation range'
        inp = [
            min(data[1][0]),
            max(data[1][0]),
            None
            ]
        inp[2] = (inp[1] - inp[0])/(len(data[1][0])-1.)
        names = ('start angle', 'end angle', 'step size')
        dictlst = [inp] * len(inp)
        elemlst = range(len(inp))
        dlg = G2G.ScrolledMultiEditor(
            G2frame,[inp] * len(inp), range(len(inp)), names,
            header='Edit simulation range',
            minvals=(0.001,0.001,0.0001),
            maxvals=(180.,180.,.1),
            )
        dlg.CenterOnParent()
        val = dlg.ShowModal()
        dlg.Destroy()
        if val != wx.ID_OK: return
        if inp[0] > inp[1]:
            end,start,step = inp
        else:                
            start,end,step = inp
        step = abs(step)
        N = int((end-start)/step)+1
        newdata = np.linspace(start,end,N,True)
        if len(newdata) < 2: return # too small a range - reject
        data[1] = [newdata,np.zeros_like(newdata),np.ones_like(newdata),
            np.zeros_like(newdata),np.zeros_like(newdata),np.zeros_like(newdata)]
        Tmin = newdata[0]
        Tmax = newdata[-1]
        G2frame.PatternTree.SetItemPyData(GetPatternTreeItemId(G2frame,item,'Limits'),
            [(Tmin,Tmax),[Tmin,Tmax]])
        UpdatePWHKPlot(G2frame,kind,item) # redisplay data screen

    def OnPlot3DHKL(event):
        refList = data[1]['RefList']
        FoMax = np.max(refList.T[8+Super])
        Hmin = np.array([int(np.min(refList.T[0])),int(np.min(refList.T[1])),int(np.min(refList.T[2]))])
        Hmax = np.array([int(np.max(refList.T[0])),int(np.max(refList.T[1])),int(np.max(refList.T[2]))])
        Vpoint = [int(np.mean(refList.T[0])),int(np.mean(refList.T[1])),int(np.mean(refList.T[2]))]
        controls = {'Type' : 'Fosq','Iscale' : False,'HKLmax' : Hmax,'HKLmin' : Hmin,
            'FoMax' : FoMax,'Scale' : 1.0,'Drawing':{'viewPoint':[Vpoint,[]],'default':Vpoint[:],
            'backColor':[0,0,0],'depthFog':False,'Zclip':10.0,'cameraPos':10.,'Zstep':0.05,
            'Scale':1.0,'oldxy':[],'viewDir':[1,0,0]},'Super':Super,'SuperVec':SuperVec}
        G2plt.Plot3DSngl(G2frame,newPlot=True,Data=controls,hklRef=refList,Title=phaseName)
        
    def OnPlotAll3DHKL(event):
        choices = GetPatternTreeDataNames(G2frame,['HKLF',])
        dlg = G2G.G2MultiChoiceDialog(G2frame, 'Select reflection sets to plot',
            'Use data',choices)
        try:
            if dlg.ShowModal() == wx.ID_OK:
                refNames = [choices[i] for i in dlg.GetSelections()]
            else:
                return
        finally:
            dlg.Destroy()
        refList = np.zeros(0)
        for name in refNames:
            Id = GetPatternTreeItemId(G2frame,G2frame.root, name)
            reflData = G2frame.PatternTree.GetItemPyData(Id)[1]
            if len(refList):
                refList = np.concatenate((refList,reflData['RefList']))
            else:
                refList = reflData['RefList']
            
        FoMax = np.max(refList.T[8+Super])
        Hmin = np.array([int(np.min(refList.T[0])),int(np.min(refList.T[1])),int(np.min(refList.T[2]))])
        Hmax = np.array([int(np.max(refList.T[0])),int(np.max(refList.T[1])),int(np.max(refList.T[2]))])
        Vpoint = [int(np.mean(refList.T[0])),int(np.mean(refList.T[1])),int(np.mean(refList.T[2]))]
        controls = {'Type' : 'Fosq','Iscale' : False,'HKLmax' : Hmax,'HKLmin' : Hmin,
            'FoMax' : FoMax,'Scale' : 1.0,'Drawing':{'viewPoint':[Vpoint,[]],'default':Vpoint[:],
            'backColor':[0,0,0],'depthFog':False,'Zclip':10.0,'cameraPos':10.,'Zstep':0.05,
            'Scale':1.0,'oldxy':[],'viewDir':[1,0,0]},'Super':Super,'SuperVec':SuperVec}
        G2plt.Plot3DSngl(G2frame,newPlot=True,Data=controls,hklRef=refList,Title=phaseName)
        
        
    def OnErrorAnalysis(event):
        G2plt.PlotDeltSig(G2frame,kind)

# </ Anton Gagin
    def OnExportQQ(event):
        G2plt.ExportDeltSig(G2frame,kind)
# Anton Gagin />  

    def OnWtFactor(event):
        try:
            val = float(wtval.GetValue())
        except ValueError:
            val = data[0]['wtFactor']
        data[0]['wtFactor'] = val
        wtval.SetValue('%.3f'%(val))
        
    def onCopyPlotCtrls(event):
        '''Respond to menu item to copy multiple sections from a histogram.
        Need this here to pass on the G2frame object. 
        '''
        G2pdG.CopyPlotCtrls(G2frame)

    def onCopySelectedItems(event):
        '''Respond to menu item to copy multiple sections from a histogram.
        Need this here to pass on the G2frame object. 
        '''
        G2pdG.CopySelectedHistItems(G2frame)
           
    data = G2frame.PatternTree.GetItemPyData(item)
#patches
    if 'wtFactor' not in data[0]:
        data[0] = {'wtFactor':1.0}
    #if isinstance(data[1],list) and kind == 'HKLF':
    if 'list' in str(type(data[1])) and kind == 'HKLF':
        RefData = {'RefList':[],'FF':[]}
        for ref in data[1]:
            RefData['RefList'].append(ref[:11]+[ref[13],])
            RefData['FF'].append(ref[14])
        data[1] = RefData
        G2frame.PatternTree.SetItemPyData(item,data)
#end patches
    if G2frame.dataDisplay:
        G2frame.dataDisplay.Destroy()
    if kind in ['PWDR','SASD']:
        SetDataMenuBar(G2frame,G2frame.dataFrame.PWDRMenu)
        G2frame.dataFrame.Bind(wx.EVT_MENU, OnErrorAnalysis, id=wxID_PWDANALYSIS)
# </ Anton Gagin
        G2frame.dataFrame.Bind(wx.EVT_MENU, OnExportQQ, id=wxID_PWDEXPORTQQ)
# Anton Gagin />          
        G2frame.dataFrame.Bind(wx.EVT_MENU, onCopySelectedItems, id=wxID_PWDCOPY)
        G2frame.dataFrame.Bind(wx.EVT_MENU, onCopyPlotCtrls, id=wxID_PLOTCTRLCOPY)
    elif kind in ['HKLF',]:
        SetDataMenuBar(G2frame,G2frame.dataFrame.HKLFMenu)
        G2frame.dataFrame.Bind(wx.EVT_MENU, OnErrorAnalysis, id=wxID_PWDANALYSIS)
        G2frame.dataFrame.Bind(wx.EVT_MENU, OnPlot3DHKL, id=wxID_PWD3DHKLPLOT)
        G2frame.dataFrame.Bind(wx.EVT_MENU, OnPlotAll3DHKL, id=wxID_3DALLHKLPLOT)
#        G2frame.dataFrame.Bind(wx.EVT_MENU, onCopySelectedItems, id=wxID_PWDCOPY)
    G2frame.dataDisplay = wx.Panel(G2frame.dataFrame)
    
    mainSizer = wx.BoxSizer(wx.VERTICAL)
    mainSizer.Add((5,5),)
    wtSizer = wx.BoxSizer(wx.HORIZONTAL)
    wtSizer.Add(wx.StaticText(G2frame.dataDisplay,-1,' Weight factor: '),0,WACV)
    wtval = wx.TextCtrl(G2frame.dataDisplay,-1,'%.3f'%(data[0]['wtFactor']),style=wx.TE_PROCESS_ENTER)
    wtval.Bind(wx.EVT_TEXT_ENTER,OnWtFactor)
    wtval.Bind(wx.EVT_KILL_FOCUS,OnWtFactor)
    wtSizer.Add(wtval,0,WACV)
    mainSizer.Add(wtSizer)
    if data[0].get('Dummy'):
        simSizer = wx.BoxSizer(wx.HORIZONTAL)
        Tmin = min(data[1][0])
        Tmax = max(data[1][0])
        num = len(data[1][0])
        step = (Tmax - Tmin)/(num-1)
        t = u'2\u03b8' # 2theta
        lbl =  u'Simulation range: {:.2f} to {:.2f} {:s}\nwith {:.4f} steps ({:d} points)'
        lbl += u'\n(Edit range resets observed intensities).'
        lbl = lbl.format(Tmin,Tmax,t,step,num)
        simSizer.Add(wx.StaticText(G2frame.dataDisplay,wx.ID_ANY,lbl),
                    0,WACV)
        but = wx.Button(G2frame.dataDisplay,wx.ID_ANY,"Edit range")
        but.Bind(wx.EVT_BUTTON,onEditSimRange)
        simSizer.Add(but,0,WACV)
        mainSizer.Add(simSizer)
    if 'Nobs' in data[0]:
        mainSizer.Add(wx.StaticText(G2frame.dataDisplay,-1,
            ' Data residual wR: %.3f%% on %d observations'%(data[0]['wR'],data[0]['Nobs'])))
        for value in data[0]:
            if 'Nref' in value:
                mainSizer.Add((5,5),)
                pfx = value.split('Nref')[0]
                name = data[0].get(pfx.split(':')[0]+'::Name','?')
                mainSizer.Add(wx.StaticText(G2frame.dataDisplay,-1,' For phase '+name+':'))
                mainSizer.Add(wx.StaticText(G2frame.dataDisplay,-1,
                    u' Unweighted phase residuals RF\u00b2: %.3f%%, RF: %.3f%% on %d reflections  '% \
                    (data[0][pfx+'Rf^2'],data[0][pfx+'Rf'],data[0][value])))
    mainSizer.Add((5,5),)
    mainSizer.Layout()    
    G2frame.dataDisplay.SetSizer(mainSizer)
    Size = mainSizer.Fit(G2frame.dataFrame)
    Size[1] += 10
    G2frame.dataFrame.setSizePosLeft(Size)
    G2frame.PatternTree.SetItemPyData(item,data)
    if kind in ['PWDR','SASD']:
        G2plt.PlotPatterns(G2frame,plotType=kind,newPlot=True)
    elif kind == 'HKLF':
        Name = G2frame.PatternTree.GetItemText(item)
        phaseName = G2pdG.IsHistogramInAnyPhase(G2frame,Name)
        if phaseName:
            pId = GetPatternTreeItemId(G2frame,G2frame.root,'Phases')
            phaseId =  GetPatternTreeItemId(G2frame,pId,phaseName)
            General = G2frame.PatternTree.GetItemPyData(phaseId)['General']
            Super = General.get('Super',0)
            SuperVec = General.get('SuperVec',[])
        else:
            Super = 0
            SuperVec = []       
        refList = data[1]['RefList']
        FoMax = np.max(refList.T[5+data[1].get('Super',0)])
        page = G2frame.G2plotNB.nb.GetSelection()
        tab = ''
        if page >= 0:
            tab = G2frame.G2plotNB.nb.GetPageText(page)
        if '3D' in tab:
            Hmin = np.array([int(np.min(refList.T[0])),int(np.min(refList.T[1])),int(np.min(refList.T[2]))])
            Hmax = np.array([int(np.max(refList.T[0])),int(np.max(refList.T[1])),int(np.max(refList.T[2]))])
            Vpoint = [int(np.mean(refList.T[0])),int(np.mean(refList.T[1])),int(np.mean(refList.T[2]))]
            Page = G2frame.G2plotNB.nb.GetPage(page)
            controls = Page.controls
            G2plt.Plot3DSngl(G2frame,newPlot=False,Data=controls,hklRef=refList,Title=phaseName)
        else:
            controls = {'Type' : 'Fo','ifFc' : True,     
                'HKLmax' : [int(np.max(refList.T[0])),int(np.max(refList.T[1])),int(np.max(refList.T[2]))],
                'HKLmin' : [int(np.min(refList.T[0])),int(np.min(refList.T[1])),int(np.min(refList.T[2]))],
                'FoMax' : FoMax,'Zone' : '001','Layer' : 0,'Scale' : 1.0,'Super':Super,'SuperVec':SuperVec}
            G2plt.PlotSngl(G2frame,newPlot=True,Data=controls,hklRef=refList)
                 
################################################################################
#####  Pattern tree routines
################################################################################           
       
def GetPatternTreeDataNames(G2frame,dataTypes):
    '''Needs a doc string
    '''
    names = []
    item, cookie = G2frame.PatternTree.GetFirstChild(G2frame.root)        
    while item:
        name = G2frame.PatternTree.GetItemText(item)
        if name[:4] in dataTypes:
            names.append(name)
        item, cookie = G2frame.PatternTree.GetNextChild(G2frame.root, cookie)
    return names
                          
def GetPatternTreeItemId(G2frame, parentId, itemText):
    '''Needs a doc string
    '''
    item, cookie = G2frame.PatternTree.GetFirstChild(parentId)
    while item:
        if G2frame.PatternTree.GetItemText(item) == itemText:
            return item
        item, cookie = G2frame.PatternTree.GetNextChild(parentId, cookie)
    return 0                

def MovePatternTreeToGrid(G2frame,item):
    '''Called from GSASII.OnPatternTreeSelChanged when a item is selected on the tree 
    '''
    pickName = G2frame.PatternTree.GetItemText(item)
    if G2frame.PickIdText == pickName:
        return
    
    oldPage = None # will be set later if already on a Phase item
    if G2frame.dataFrame:
        SetDataMenuBar(G2frame)
        if G2frame.dataFrame.GetLabel() == 'Comments':
            try:
                data = [G2frame.dataDisplay.GetValue()]
                G2frame.dataDisplay.Clear() 
                Id = GetPatternTreeItemId(G2frame,G2frame.root, 'Comments')
                if Id: G2frame.PatternTree.SetItemPyData(Id,data)
            except:     #clumsy but avoids dead window problem when opening another project
                pass
        elif G2frame.dataFrame.GetLabel() == 'Notebook':
            try:
                data = [G2frame.dataDisplay.GetValue()]
                G2frame.dataDisplay.Clear() 
                Id = GetPatternTreeItemId(G2frame,G2frame.root, 'Notebook')
                if Id: G2frame.PatternTree.SetItemPyData(Id,data)
            except:     #clumsy but avoids dead window problem when opening another project
                pass
        elif 'Phase Data for' in G2frame.dataFrame.GetLabel():
            if G2frame.dataDisplay: 
                oldPage = G2frame.dataDisplay.GetSelection()
        G2frame.dataFrame.Clear()
        G2frame.dataFrame.SetLabel('')
    else:
        #create the frame for the data item window
        G2frame.dataFrame = DataFrame(parent=G2frame.mainPanel,frame=G2frame)
        G2frame.dataFrame.PhaseUserSize = None
        
    G2frame.dataFrame.Raise()            
    G2frame.PickId = item
    G2frame.PickIdText = None
    parentID = G2frame.root
    #for i in G2frame.ExportPattern: i.Enable(False)
    defWid = [250,150]
    if item != G2frame.root:
        parentID = G2frame.PatternTree.GetItemParent(item)
    if G2frame.PatternTree.GetItemParent(item) == G2frame.root:
        G2frame.PatternId = item
        if G2frame.PatternTree.GetItemText(item) == 'Notebook':
            SetDataMenuBar(G2frame,G2frame.dataFrame.DataNotebookMenu)
            G2frame.PatternId = 0
            #for i in G2frame.ExportPattern: i.Enable(False)
            data = G2frame.PatternTree.GetItemPyData(item)
            UpdateNotebook(G2frame,data)
        elif G2frame.PatternTree.GetItemText(item) == 'Controls':
            G2frame.PatternId = 0
            #for i in G2frame.ExportPattern: i.Enable(False)
            data = G2frame.PatternTree.GetItemPyData(item)
            if not data:           #fill in defaults
                data = copy.copy(G2obj.DefaultControls)    #least squares controls
                G2frame.PatternTree.SetItemPyData(item,data)                             
            for i in G2frame.Refine: i.Enable(True)
            G2frame.EnableSeqRefineMenu()
            UpdateControls(G2frame,data)
        elif G2frame.PatternTree.GetItemText(item) == 'Sequential results':
            data = G2frame.PatternTree.GetItemPyData(item)
            UpdateSeqResults(G2frame,data)
        elif G2frame.PatternTree.GetItemText(item) == 'Covariance':
            data = G2frame.PatternTree.GetItemPyData(item)
            G2frame.dataFrame.setSizePosLeft(defWid)
            text = ''
            if 'Rvals' in data:
                Nvars = len(data['varyList'])
                Rvals = data['Rvals']
                text = '\nFinal residuals: \nwR = %.3f%% \nchi**2 = %.1f \nGOF = %.2f'%(Rvals['Rwp'],Rvals['chisq'],Rvals['GOF'])
                text += '\nNobs = %d \nNvals = %d'%(Rvals['Nobs'],Nvars)
                if 'lamMax' in Rvals:
                    text += '\nlog10 MaxLambda = %.1f'%(np.log10(Rvals['lamMax']))
            wx.TextCtrl(parent=G2frame.dataFrame,size=G2frame.dataFrame.GetClientSize(),
                value='See plot window for covariance display'+text,style=wx.TE_MULTILINE)
            G2plt.PlotCovariance(G2frame,data)
        elif G2frame.PatternTree.GetItemText(item) == 'Constraints':
            data = G2frame.PatternTree.GetItemPyData(item)
            G2cnstG.UpdateConstraints(G2frame,data)
        elif G2frame.PatternTree.GetItemText(item) == 'Rigid bodies':
            data = G2frame.PatternTree.GetItemPyData(item)
            G2cnstG.UpdateRigidBodies(G2frame,data)
        elif G2frame.PatternTree.GetItemText(item) == 'Restraints':
            data = G2frame.PatternTree.GetItemPyData(item)
            Phases = G2frame.GetPhaseData()
            phase = ''
            phaseName = ''
            if Phases:
                phaseName = Phases.keys()[0]
            G2frame.dataFrame.setSizePosLeft(defWid)
            G2restG.UpdateRestraints(G2frame,data,Phases,phaseName)
        elif 'IMG' in G2frame.PatternTree.GetItemText(item):
            G2frame.Image = item
            G2frame.dataFrame.SetTitle('Image Data')
            data = G2frame.PatternTree.GetItemPyData(GetPatternTreeItemId( \
                G2frame,item,'Image Controls'))
            G2imG.UpdateImageData(G2frame,data)
            G2plt.PlotImage(G2frame,newPlot=True)
        elif 'PKS' in G2frame.PatternTree.GetItemText(item):
            G2plt.PlotPowderLines(G2frame)
        elif 'PWDR' in G2frame.PatternTree.GetItemText(item):
            #for i in G2frame.ExportPattern: i.Enable(True)
            if G2frame.EnablePlot:
                UpdatePWHKPlot(G2frame,'PWDR',item)
        elif 'SASD' in G2frame.PatternTree.GetItemText(item):
            #for i in G2frame.ExportPattern: i.Enable(True)
            if G2frame.EnablePlot:
                UpdatePWHKPlot(G2frame,'SASD',item)
        elif 'HKLF' in G2frame.PatternTree.GetItemText(item):
            G2frame.Sngl = True
            UpdatePWHKPlot(G2frame,'HKLF',item)
        elif 'PDF' in G2frame.PatternTree.GetItemText(item):
            G2frame.PatternId = item
            for i in G2frame.ExportPDF: i.Enable(True)
            G2plt.PlotISFG(G2frame,type='S(Q)')
        elif G2frame.PatternTree.GetItemText(item) == 'Phases':
            G2frame.dataFrame.setSizePosLeft(defWid)
            wx.TextCtrl(parent=G2frame.dataFrame,size=G2frame.dataFrame.GetClientSize(),
                value='Select one phase to see its parameters')            
    elif 'I(Q)' in G2frame.PatternTree.GetItemText(item):
        G2frame.PatternId = G2frame.PatternTree.GetItemParent(item)
        data = G2frame.PatternTree.GetItemPyData(GetPatternTreeItemId(G2frame,G2frame.PatternId,'PDF Controls'))
        G2pdG.UpdatePDFGrid(G2frame,data)
        G2plt.PlotISFG(G2frame,type='I(Q)',newPlot=True)
    elif 'S(Q)' in G2frame.PatternTree.GetItemText(item):
        G2frame.PatternId = G2frame.PatternTree.GetItemParent(item)
        data = G2frame.PatternTree.GetItemPyData(GetPatternTreeItemId(G2frame,G2frame.PatternId,'PDF Controls'))
        G2pdG.UpdatePDFGrid(G2frame,data)
        G2plt.PlotISFG(G2frame,type='S(Q)',newPlot=True)
    elif 'F(Q)' in G2frame.PatternTree.GetItemText(item):
        G2frame.PatternId = G2frame.PatternTree.GetItemParent(item)
        data = G2frame.PatternTree.GetItemPyData(GetPatternTreeItemId(G2frame,G2frame.PatternId,'PDF Controls'))
        G2pdG.UpdatePDFGrid(G2frame,data)
        G2plt.PlotISFG(G2frame,type='F(Q)',newPlot=True)
    elif 'G(R)' in G2frame.PatternTree.GetItemText(item):
        G2frame.PatternId = G2frame.PatternTree.GetItemParent(item)
        data = G2frame.PatternTree.GetItemPyData(GetPatternTreeItemId(G2frame,G2frame.PatternId,'PDF Controls'))
        G2pdG.UpdatePDFGrid(G2frame,data)
        G2plt.PlotISFG(G2frame,type='G(R)',newPlot=True)            
    elif G2frame.PatternTree.GetItemText(parentID) == 'Phases':
        data = G2frame.PatternTree.GetItemPyData(item)
        G2phG.UpdatePhaseData(G2frame,item,data,oldPage)
    elif G2frame.PatternTree.GetItemText(item) == 'Comments':
        SetDataMenuBar(G2frame,G2frame.dataFrame.DataCommentsMenu)
        G2frame.PatternId = G2frame.PatternTree.GetItemParent(item)
        data = G2frame.PatternTree.GetItemPyData(item)
        UpdateComments(G2frame,data)
    elif G2frame.PatternTree.GetItemText(item) == 'Image Controls':
        G2frame.dataFrame.SetTitle('Image Controls')
        G2frame.Image = G2frame.PatternTree.GetItemParent(item)
        masks = G2frame.PatternTree.GetItemPyData(
            GetPatternTreeItemId(G2frame,G2frame.Image, 'Masks'))
        data = G2frame.PatternTree.GetItemPyData(item)
        G2imG.UpdateImageControls(G2frame,data,masks)
        G2plt.PlotImage(G2frame)
    elif G2frame.PatternTree.GetItemText(item) == 'Masks':
        G2frame.dataFrame.SetTitle('Masks')
        G2frame.Image = G2frame.PatternTree.GetItemParent(item)
        data = G2frame.PatternTree.GetItemPyData(item)
        G2imG.UpdateMasks(G2frame,data)
        G2plt.PlotImage(G2frame)
    elif G2frame.PatternTree.GetItemText(item) == 'Stress/Strain':
        G2frame.dataFrame.SetTitle('Stress/Strain')
        G2frame.Image = G2frame.PatternTree.GetItemParent(item)
        data = G2frame.PatternTree.GetItemPyData(item)
        G2plt.PlotImage(G2frame)
        G2plt.PlotStrain(G2frame,data,newPlot=True)
        G2imG.UpdateStressStrain(G2frame,data)
    elif G2frame.PatternTree.GetItemText(item) == 'PDF Controls':
        G2frame.PatternId = G2frame.PatternTree.GetItemParent(item)
        for i in G2frame.ExportPDF: i.Enable(True)
        data = G2frame.PatternTree.GetItemPyData(item)
        G2pdG.UpdatePDFGrid(G2frame,data)
        G2plt.PlotISFG(G2frame,type='I(Q)')
        G2plt.PlotISFG(G2frame,type='S(Q)')
        G2plt.PlotISFG(G2frame,type='F(Q)')
        G2plt.PlotISFG(G2frame,type='G(R)')
    elif G2frame.PatternTree.GetItemText(item) == 'Peak List':
        G2frame.PatternId = G2frame.PatternTree.GetItemParent(item)
        for i in G2frame.ExportPeakList: i.Enable(True)
        data = G2frame.PatternTree.GetItemPyData(item)
#patch
        if 'list' in str(type(data)):
            data = {'peaks':data,'sigDict':{}}
            G2frame.PatternTree.SetItemPyData(item,data)
#end patch
        G2pdG.UpdatePeakGrid(G2frame,data)
        G2plt.PlotPatterns(G2frame)
    elif G2frame.PatternTree.GetItemText(item) == 'Background':
        G2frame.PatternId = G2frame.PatternTree.GetItemParent(item)
        data = G2frame.PatternTree.GetItemPyData(item)
        G2pdG.UpdateBackground(G2frame,data)
        G2plt.PlotPatterns(G2frame)
    elif G2frame.PatternTree.GetItemText(item) == 'Limits':
        G2frame.PatternId = G2frame.PatternTree.GetItemParent(item)
        datatype = G2frame.PatternTree.GetItemText(G2frame.PatternId)[:4]
        data = G2frame.PatternTree.GetItemPyData(item)
        G2pdG.UpdateLimitsGrid(G2frame,data,datatype)
        G2plt.PlotPatterns(G2frame,plotType=datatype)
    elif G2frame.PatternTree.GetItemText(item) == 'Instrument Parameters':
        G2frame.PatternId = G2frame.PatternTree.GetItemParent(item)
        data = G2frame.PatternTree.GetItemPyData(item)[0]
        G2pdG.UpdateInstrumentGrid(G2frame,data)
        if 'P' in data['Type'][0]:          #powder data only
            G2plt.PlotPeakWidths(G2frame)
    elif G2frame.PatternTree.GetItemText(item) == 'Models':
        G2frame.PatternId = G2frame.PatternTree.GetItemParent(item)
        data = G2frame.PatternTree.GetItemPyData(item)
        G2pdG.UpdateModelsGrid(G2frame,data)
        G2plt.PlotPatterns(G2frame,plotType='SASD')
        if len(data['Size']['Distribution']):
            G2plt.PlotSASDSizeDist(G2frame)
    elif G2frame.PatternTree.GetItemText(item) == 'Substances':
        G2frame.PatternId = G2frame.PatternTree.GetItemParent(item)
        data = G2frame.PatternTree.GetItemPyData(item)
        G2pdG.UpdateSubstanceGrid(G2frame,data)
    elif G2frame.PatternTree.GetItemText(item) == 'Sample Parameters':
        G2frame.PatternId = G2frame.PatternTree.GetItemParent(item)
        data = G2frame.PatternTree.GetItemPyData(item)
        datatype = G2frame.PatternTree.GetItemPyData(G2frame.PatternId)[2][:4]

        if 'Temperature' not in data:           #temp fix for old gpx files
            data = {'Scale':[1.0,True],'Type':'Debye-Scherrer','Absorption':[0.0,False],'DisplaceX':[0.0,False],
                'DisplaceY':[0.0,False],'Diffuse':[],'Temperature':300.,'Pressure':1.0,
                    'FreePrm1':0.,'FreePrm2':0.,'FreePrm3':0.,
                    'Gonio. radius':200.0}
            G2frame.PatternTree.SetItemPyData(item,data)
    
        G2pdG.UpdateSampleGrid(G2frame,data)
        G2plt.PlotPatterns(G2frame,plotType=datatype)
    elif G2frame.PatternTree.GetItemText(item) == 'Index Peak List':
        G2frame.PatternId = G2frame.PatternTree.GetItemParent(item)
        for i in G2frame.ExportPeakList: i.Enable(True)
        data = G2frame.PatternTree.GetItemPyData(item)
#patch
        if len(data) != 2:
            data = [data,[]]
            G2frame.PatternTree.SetItemPyData(item,data)
#end patch
        G2pdG.UpdateIndexPeaksGrid(G2frame,data)
        if 'PKS' in G2frame.PatternTree.GetItemText(G2frame.PatternId):
            G2plt.PlotPowderLines(G2frame)
        else:
            G2plt.PlotPatterns(G2frame)
    elif G2frame.PatternTree.GetItemText(item) == 'Unit Cells List':
        G2frame.PatternId = G2frame.PatternTree.GetItemParent(item)
        data = G2frame.PatternTree.GetItemPyData(item)
        if not data:
            data.append([0,0.0,4,25.0,0,'P1',1,1,1,90,90,90]) #zero error flag, zero value, max Nc/No, start volume
            data.append([0,0,0,0,0,0,0,0,0,0,0,0,0,0])      #Bravais lattice flags
            data.append([])                                 #empty cell list
            data.append([])                                 #empty dmin
            data.append({})                                 #empty superlattice stuff
            G2frame.PatternTree.SetItemPyData(item,data)                             
#patch
        if len(data) < 5:
            data.append({'Use':False,'ModVec':[0,0,0.1],'maxH':1,'ssSymb':''})                                 #empty superlattice stuff
            G2frame.PatternTree.SetItemPyData(item,data)  
#end patch
        G2pdG.UpdateUnitCellsGrid(G2frame,data)
        if 'PKS' in G2frame.PatternTree.GetItemText(G2frame.PatternId):
            G2plt.PlotPowderLines(G2frame)
        else:
            G2plt.PlotPatterns(G2frame)
    elif G2frame.PatternTree.GetItemText(item) == 'Reflection Lists':   #powder reflections
        G2frame.PatternId = G2frame.PatternTree.GetItemParent(item)
        data = G2frame.PatternTree.GetItemPyData(item)
        G2frame.RefList = ''
        if len(data):
            G2frame.RefList = data.keys()[0]
        G2pdG.UpdateReflectionGrid(G2frame,data)
        G2plt.PlotPatterns(G2frame)
    elif G2frame.PatternTree.GetItemText(item) == 'Reflection List':    #HKLF reflections
        G2frame.PatternId = G2frame.PatternTree.GetItemParent(item)
        name = G2frame.PatternTree.GetItemText(G2frame.PatternId)
        data = G2frame.PatternTree.GetItemPyData(G2frame.PatternId)
        G2pdG.UpdateReflectionGrid(G2frame,data,HKLF=True,Name=name)

    if G2frame.PickId:
        G2frame.PickIdText = G2frame.GetTreeItemsList(G2frame.PickId)
    G2frame.dataFrame.Raise()

def SetDataMenuBar(G2frame,menu=None):
    '''Set the menu for the data frame. On the Mac put this
    menu for the data tree window instead.

    Note that data frame items do not have menus, for these (menu=None)
    display a blank menu or on the Mac display the standard menu for
    the data tree window.
    '''
    if sys.platform == "darwin":
        if menu is None:
            G2frame.SetMenuBar(G2frame.GSASIIMenu)
        else:
            G2frame.SetMenuBar(menu)
    else:
        if menu is None:
            G2frame.dataFrame.SetMenuBar(G2frame.dataFrame.BlankMenu)
        else:
            G2frame.dataFrame.SetMenuBar(menu)

def HowDidIgetHere():
    '''Show a traceback with calls that brought us to the current location.
    Used for debugging.
    '''
    import traceback
    print 70*'*'    
    for i in traceback.format_list(traceback.extract_stack()[:-1]): print(i.strip.rstrip())
    print 70*'*'    
        