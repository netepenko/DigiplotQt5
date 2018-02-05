# -*- coding: utf-8 -*-
"""
Created on Fri Mar 17 17:11:20 2017

@author: Alex
"""
from PyQt5 import QtWidgets, QtCore, QtGui
import FFT
import os
import h5py
import numpy as np
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
from matplotlib.widgets import RectangleSelector
from matplotlib.backends.backend_qt5 import NavigationToolbar2QT
from LT.parameterfile import pfile
import time
import LT.box


#restructarize later
from static_functions import find_peaks, get_window

convert_int=True #change to self detetction

class NavigationToolbar(NavigationToolbar2QT):
    def __init__(self, canvas, parent, coordinates=True):
        NavigationToolbar2QT.__init__(self, canvas, parent, coordinates=True)
        self.axes = canvas.parent().axes
        self.myinitdata = 0.
        self.N = 10000
        
    def back(self, *args):
        NavigationToolbar2QT.back(self, *args)
        self.my_rescale()
        
    def forward(self, *args):
        NavigationToolbar2QT.forward(self, *args)
        self.my_rescale()
    
    def release_pan(self, event):
        NavigationToolbar2QT.release_pan(self, event)
        self.my_rescale()
        
    def home(self, *args):
        NavigationToolbar2QT.home(self, *args)
        self.my_rescale()
        
    def release_zoom(self, event):
        NavigationToolbar2QT.release_zoom(self, event)
        self.my_rescale()     
        
    #rescale and replot  
    def my_rescale(self):
        
        t=self.myinitdata[0]
        V=self.myinitdata[1]
        rest=self.myinitdata[2:]
        ax=self.axes
        (xmin, xmax)=ax.get_xlim()
        rng=xmax-xmin
        #get original data which is in plotting region
        org_points=np.where(((xmin - 2*rng) < t) & (t < (xmax + 2*rng)))
        
        to=t[org_points]
        Vo=V[org_points]
        n=1 #take every nth point from data
        if len(to)/5> self.N:
            n=int(len(to)/5 / self.N)
        tcut=to[::n]
        Vcut=Vo[::n]
        
        #try: peaks=ax.lines[1]
        ax.set_color_cycle(None)
        #ax.cla()
        ax.lines.remove(ax.figure.canvas.parent().all_plot[0])
        ax.autoscale(enable=False, axis='x', tight=True)
        ax.set_autoscaley_on(False)
        ax.figure.canvas.parent().all_plot=ax.plot(tcut,Vcut,*rest, zorder=0)



class NumberDialog(QtWidgets.QDialog):
#     def __init__(self, data, title="Enter Parameters", \
#                  labels=['First Label','Second Label'], \
#                  keys=['first','second'], \
#                  about_txt = 'about text'):
#    
     def __init__(self, data, parent = None, title = 'Enter parameters', labels=['First Label','Second Label'],
                  keys=['first','second'], about_txt = 'about_txt'):
          
          self.parent=parent
          self.keys=keys
          self.data=data
          QtWidgets.QDialog.__init__(self, parent)#, QtCore.Qt.WindowContextHelpButtonHint)
          self.setWindowFlags(self.windowFlags() ^ QtCore.Qt.WindowContextHelpButtonHint)
          self.setWindowTitle(title)
          self.layout=QtWidgets.QGridLayout(self)
          about = QtWidgets.QLabel(about_txt, self)
          ok=QtWidgets.QPushButton('Ok', self)
          ok.setDefault(True)
          cancel=QtWidgets.QPushButton('Cancel', self)
          ok.clicked.connect(self.OnOk)
          cancel.clicked.connect(self.OnCancel)
          sepline=QtWidgets.QFrame()
          sepline.setFrameShape(QtWidgets.QFrame.HLine)
          sepline.setFrameShadow(QtWidgets.QFrame.Sunken)
          
          self.layout.addWidget(about, 0, 0)
          self.layout.addWidget(sepline, 1, 0, 1, 2)
          
          
          # loop over keys to add controls and validators
          nrow = len(keys)+1
          #qle - qlinedit dictionary to retrieve data later
          self.qle={}
          for i, key in enumerate(keys):
               Val_l  = QtWidgets.QLabel(labels[i], self)
               Val_t  = QtWidgets.QLineEdit(self)
               Val_t.setValidator(QtGui.QDoubleValidator(self))
               
               Val_t.textChanged.connect(self.check_state)
               Val_t.setText(data.get(key))
               self.layout.addWidget(Val_l, i+2, 0)
               self.layout.addWidget(Val_t, i+2, 1)
               
               self.qle[key]=Val_t
                       
          self.layout.addWidget(ok, nrow+2, 0)
          self.layout.addWidget(cancel, nrow+2, 1)
          
     def check_state(self, *args, **kwargs):
         
        sender = self.sender()
        validator = sender.validator()
        state = validator.validate(sender.text(), 0)[0]
        if state == QtGui.QValidator.Acceptable:
            color = '#c4df9b' # green
        elif state == QtGui.QValidator.Intermediate:
            color = '#fff79a' # yellow
        else:
            color = '#f6989d' # red
        sender.setStyleSheet('QLineEdit { background-color: %s }' % color)
        
     def OnOk(self):
         for i, key in enumerate(self.keys):
             try:
                 self.data[key]=self.qle[key].text()
             except:
                 mb=QtWidgets.QMessageBox(self)
                 mb.setWindowTitle('Entry error')
                 mb.setText("Please enter the missing parameter")
                 mb.exec_()
                 self.qle[key].setFocus()
                 return
         self.close()
     
     def OnCancel(self):
         self.destroy()
        

class PlotFrame(QtWidgets.QMainWindow):
     def __init__(self, parent):
          QtWidgets.QMainWindow.__init__(self)
         
          self.setWindowTitle('DigiPlot')
          
          # Parameters
          self.par = {}
          self.par["measure"] = True
          self.par["limits"] = False
          self.par["use_limits"] = False
          self.par["auto_histo"] = True
          self.par["filtered"] = False
          self.par["draw_lines"] = False
          self.par["Detector channel"] = 0
          self.par["tmin"] = -1.e30
          self.par["tmax"] = 1.e30
          self.par["xmin"] = -1.e30
          self.par["xmax"] = 1.e30
          self.par["scale_x"] = 1.0
          self.par["scale_y"] = 1.0
          self.par["Vstep"] = 1.0
          self.par["Vthreshold"] = 0.0
          self.par["Vhmin"] = 0.0
          self.par["Vhmax"] = 1.0
          self.par["VNbins"] = 100
          # time slices
          self.par["ts_width"]=1e-3
          self.par["ts_start"] = 0.0
          self.par["ts_stop"] = 1.0
          self.par["ts_Vthreshold"] = 0.0   # for time slice histograms
          self.par["ts_Vhmin"] = 0.0
          self.par["ts_Vhmax"] = 1.0
          self.par["ts_VNbins"] = 100
          self.dt = -1.
          # inital settings
          self.ndata = 0
          self.t = None
          self.V = None
          self.Vinv = None
          self.MINTAB = None
          self.MAXTAB = None
          self.histo = None
          self.histoframe = None
          self.tsplotframe = None
          self.fft = FFT.FFT()
          self.t_slice = None
          self.V_slice = None
          self.ts_tp = None
          self.ts_Vp = None
          self.ts_av = None
          self.ts_counts = None
          self.ts_FFT = None
          self.current_histo = 0
          self.current_slice = 0
          self.histos = []
          self.hwsproc = [] #hws files selection 
          self.hwsnew = []
          self.hwstodo = []
          self.datadir = 'C://Users//Alex//Desktop//'
          #load list of processed files
          self.proclistload() 
          
          # file information data
          self.parent = parent
          self.full_name = None
          self.dir = None
          self.par_dir = None
          self.hist_dir = None
          self.res_dir = None
          self.name = None
          self.ext = None
          self.fileDataOK = False
         
          #Change background color
          self.palette=QtGui.QPalette()
          self.palette.setColor(QtGui.QPalette.Background,QtCore.Qt.white)
          
          #menu bar
          self.mbar = self.menuBar()
          # setup menues
          
          # create status bar
          self.stBar1 = QtWidgets.QLabel('No file loaded')
          self.stBar1.setFrameStyle(2)
          self.stBar1.setFrameShadow(48)
          self.stBar2 = QtWidgets.QLabel('No plot yet')
          self.stBar2.setFrameStyle(2)
          self.stBar2.setFrameShadow(48)
          
          self.statusBar().addWidget(self.stBar1, 1)
          self.statusBar().addWidget(self.stBar2, 1)
          
           # setup figure
          self.figure = Figure(figsize=(10,6))
          self.axes = self.figure.add_subplot(111)
          self.figure_canvas = FigureCanvas(self.figure)
          #  Note that event is a MplEvent
          self.figure_canvas.mpl_connect('motion_notify_event', self.UpdateStatusBar)
          self.figure_canvas.mpl_connect('figure_enter_event', self.ChangeCursor)
          
          # put figure cancas in center of Main Window
          self.setCentralWidget(self.figure_canvas)
          
          # add ToolBar
          self.toolbar = NavigationToolbar(self.figure_canvas, self)
          self.addToolBar(self.toolbar)
          
          self.RS = RectangleSelector(self.axes, self.LineSelectCallback,
                                      drawtype='box', useblit=True,
                                      button=[1,3], # don't use middle button
                                      minspanx=5, minspany=5,
                                      spancoords='pixels')
          self.setPalette(self.palette)
          self.setGeometry(150,150,800,600)
          self.createMenubar()
          self.show()
          
     
     def my_plot(self, *args, **kwargs):
        
        N=kwargs.pop('N', 10000)
        #ax = pl.gca()
        ax=self.axes
        # allow callers to override the hold state by passing hold=True|False
        
        #saves input data into toolbar class object and reuses later
        self.toolbar.myinitdata = args
        self.toolbar.N = N
        
        #dropping points for first plot
        t=args[0]
        V=args[1]
        rest=args[2:]
        n=1 #take every nth point from data
        if len(t)> N:
            n=int(len(t)/N)
        tcut=t[::n]
        Vcut=V[::n]
    
    
        #washold = ax.ishold()
        #hold = kwargs.pop('hold', None)
        #if hold is not None:
        #    ax.hold(hold)
        #try:
        #    ret = ax.plot(tcut, Vcut, *rest, **kwargs)
        #finally:
        #    ax.hold(washold)
        #return ret 
        ax.autoscale(enable=True, axis='x', tight=True)
        #ax.set_autoscaley_on(True)
        return ax.plot(tcut, Vcut, *rest, **kwargs)
        
     def proclistload(self): #load list of processed files
        path_to_watch = self.datadir
        dircont = os.listdir (path_to_watch)
                
        if "ProcFiles.data" in dircont:
            print "\nLoad list of already processed files."
            fopen = open(self.datadir+"//ProcFiles.data")
            self.hwsproc=[x.strip('\n') for x in fopen.readlines()]
            #print self.hwsproc
            #for x in self.hwsproc: print "\n", x
            
        else: 
            self.hwsproc=[]
            print "\nNo files were processed yet in current data directory."
                
     def proclistsave(self): #save list of processed files
        path_to_watch = self.datadir
        dircont = os.listdir (path_to_watch)
                
        if "ProcFiles.data" in dircont:
            print "\nAppend list of processed files."
            fopen = open(self.datadir+"//ProcFiles.data", 'a')
            for item in self.hwsproc:
                fopen.write("%s\n" % item)
            fopen.close()
        else: 
            fopen = open(self.datadir+"//ProcFiles.data", 'w')
            for item in self.hwsproc:
                fopen.write("%s\n" % item)
            fopen.close()
            print "\nCreate list of processed files."
         
         
     def menuData(self): # data for the menu
          return(
               ("&File",                                        # menu label
                ("&Open", "Select data files", self.OnSelectFiles),  # menu items consisting of: label, text, handler
                ("&Reload", "Reload data files", self.OnReload),  
                ("&Load Parameters", "Load parameters", self.OnLoadParameters), 
                ("&Save Parameters", "Save parameters", self.OnSaveParameters), 
                ("&Save Histograms", "Save all histograms", self.OnSaveHistos), 
                ("&Save Event Rate", "Save particle rates from peak detection", self.OnTSsaverate), 
                ("&Quit", "Quit program", self.OnCloseWindow)), # label, status and handler
               #
               ("&Actions",
                ("&Scan data directory", "Scan directory", self.OnScan),
                (None, None, None),  # creates a separator bar in the menu
                ("&Plot", "Plot data", self.OnPlot),
                (None, None, None),  # creates a separator bar in the menu
                ("&FindPeaks", "find peaks in the plotting data", self.OnFindPeaks),
                (None, None, None),
                ("&Histogram", "histogram peak data", self.OnHistogram),
                ("&Delete Histograms", "delete all histograms", self.OnDeleteHistogram),
                (None, None, None),
                ("&Clear Figure", "Clear figure", self.OnClear)),
               #
               ("&Parameters",
                ("&Data directory", "Set data directory", self.OnSetScanDir), 
                (None, None, None), 
		("&Detector Channel", "Set detector channel", self.OnSelectChannel), 
                (None, None, None),			   
                ("&Time Range", "Set time slot", self.OnSelectTimeSlot), 
                (None, None, None),
                ("&Peak Finding", "Set parameters for peak finding", self.OnSetPeakpar), 
                (None, None, None),
                ("&Histogram", "Set histogram parameters", self.OnSetHisto),
                (None, None, None),
                ("FF&T Filter", "Set FFT filter parameters", self.OnSetFFTfilterpar),
                (None, None, None),
                ("Time &Slicing", "Set time slice parameters", self.OnSetTimeSlice)),
               ("&FFT",
                ("&Calculate", "Calculate FFT", self.OnFFTcalc), 
                ("&Freq.Cut", "Apply freq. Cut", self.OnFFTcutfilter),
                ("&R-Cut", "Apply r-value cut", self.OnFFTrfilter),
                (None, None, None),
                ("&Plot PS", "Plot power spectrum", self.OnFFTplotps),
                ("&Plot an", "Plot sin coeff", self.OnFFTplotan),
                ("&Plot bn", "Plot cos coeff", self.OnFFTplotbn),
                (None, None, None),
                ("&Invert", "Calc. inverse FFT", self.OnFFTinvert)),
               ("&Time Slices",
                ("&Slice", "Calculate time slices", self.OnTScalc),
                ("&Plot Slices", "Plot slices", self.OnTSshowslice),
                ("&Find Peaks", "Find peaks", self.OnTSfindpeaks),
                ("&FFT", "Calc FFT", self.OnTSFFTcalc),
                ("&Rate Plot", "Plot rates", self.OnTSplotrate),
                ("&Create Histograms", "Create Histograms",self.OnTShistogram ),
                ("&Delete Histograms", "delete all histograms", self.OnDeleteHistogram),
                ("&Show Histograms", "Plot Histograms", self.OnTSshowhistos)),
               ("&Options",
                # ("&Unfiltered", "Use unfiltered data", self.OnNothing, 'RADIO'),
                (None, None, None),
                ("&SetLimits", "Set limits for data processing", self.OnLimits, 'RADIO'),
                ("&Measure", "Measure differences in t and V", self.OnMeasure, 'RADIO'),
                ("&Use Limits ", "Use filtered data", self.OnTogglelimits, 'CHECK'),
                ("&Plot Lines ", "Draw lines or Points", self.OnToggleLines, 'CHECK'),
                (None, None, None),
                ("&Auto Histogram Limits ", "Do not calculate histogram limits automatically", self.OnToggleHistolimits, 'CHECK'),
                (None, None, None),
                ("&Filtered", "Use filtered data", self.OnUsefiltered, 'CHECK')),
               )

     def createMenubar(self):
         group=QtWidgets.QActionGroup(self) 
         for eachMenuData in self.menuData():
               
              menuLabel = eachMenuData[0]
              menuItems = eachMenuData[1:]
               

              menu = QtWidgets.QMenu(menuLabel, self)
              for eachItem in menuItems:
                   if not eachItem[0]:
                       menu.addSeparator()
                       continue
                   subMenu=QtWidgets.QAction(eachItem[0],self) 
                   if len(eachItem) == 4: #never true
                        if eachItem[3]=='CHECK' or eachItem[3]=='RADIO':
                            subMenu.setCheckable(True)
                        if eachItem[3]=='RADIO':
                            subMenu.setActionGroup(group)
                            
                   if 'Measure' in eachItem[0] and self.par['measure']:
                       subMenu.setChecked(True)
                   if 'SetLimits' in eachItem[0] and self.par['limits']:
                       subMenu.setChecked(True)
                      
                       
                   subMenu.triggered.connect(eachItem[2])
                   menu.addAction(subMenu)
                   
              
              
              self.mbar.addMenu(menu)          
          

 
          

     def set_file_info(self,filename):
          dir, fname = os.path.split(filename)
          name, ext = os.path.splitext(fname)
          self.dir = dir+'//'
          self.name = name
          self.ext = ext
          # that's it

     #----------------------------------------------------------------------
     # File menu routines   
     #----------------------------------------------------------------------
     def OnScan(self): #scan for new files in data directory
        path_to_watch = self.datadir
        before = self.hwsproc
        after = os.listdir (path_to_watch)
        newf = False
        self.hwsnew = []
        for f in after:
            if f.endswith(".hws") and not f in before:
                print "Unprocessed file -- ", f
                newf = True
                self.hwsnew.append(f)
        print "Directory scanned for unprocessed .hws files \n"    
        if newf: 
            process = self.workonfiles()
            if process:
                print 'da na'
                #self.choosefiles=Repository(self, -1, 'Select files to process',self.hwsnew)
    
     #process selected files           
    
     def Process(self): 
        start_time=time.time()
        for f in self.hwstodo:
                #print "\n", f 
                self.OpenFile(self.datadir+"//"+f)
                
                #    self.LoadParameters("C://Users//plasma//Desktop//par.data")
                #    self.OnPlot(None)
                for i in range(0,2):
                    self.RatesForAll(None,f)
                
            
     
     
     
        end_time=time.time()
        work_time=end_time-start_time
        print "Elapsed time ", work_time 
                       
                #self.OnFindPeaks(None)
                #self.OnHistogram(None)
                #os.mkdir(self.dir+self.name+'//')
                #self.histoframe.figure.savefig(self.dir+self.name+'//Histogram.png')
                #self.OnTScalc(None)
                #self.OnTSfindpeaks(None)
                #self.OnClear(None)
                #self.OnTSplotrate(None)
                #self.figure.savefig(self.dir+self.name+'//RatePlot.png')
                #self.hwsproc.append(f)
                
                #for i,h in enumerate(self.histos):
                #    h_file = self.dir + self.name +'//histos_{0}.data'.format(i)
                #    h.save(h_file)
                #print " all histograms saved"
        #self.proclistsave()
        #print self.hwstodo
        self.hwstodo = []                
       
     def RatesForAll(self,parent,ftw):
         self.OnClear(None)
         for i in range(0,2):
             self.LoadParameters(self.datadir+"//par_ch_"+str(i)+".data")
             self.OpenFile(self.datadir+"//"+ftw)
             self.OnTScalc(None)
             self.OnTSfindpeaks(None)
             self.OnTSplotrate(None)
         return None 
     
     def workonfiles(self): # dialog if to work on unrpocessed files
         dlg = QtWidgets.QMessageBox.question(self, 'Message', 'Process the found files?',
                                            QtWidgets.QMessageBox.Yes, QtWidgets.QMessageBox.No)
         
         if dlg == QtWidgets.QMessageBox.Yes:
             return True
         else:
             return False
             
     def OnSetScanDir(self):
         
         print "Set data directory"
         dir_dlg=QtWidgets.QFileDialog.getExistingDirectory(self, 'Select a directory')
         if dir_dlg != '':
             # User has selected something, get the path, set the window's title to the path
             # store relevant file information
             self.datadir=dir_dlg
         else:
             print "so, you changed your mind, I will do nothing"
             return
         print 'Data directory: ', self.datadir
         #self.proclistload()
                    
     def OnSelectFiles(self):
         # Create a file-open dialog in the current directory
         
         filetypes = '*.hws'
         # use a regular file dialog
         if self.dir == None:
             self.dir = os.getcwd()
             
         file_dlg=QtWidgets.QFileDialog.getOpenFileName(self, 'Select a file', self.dir, filetypes )
         
         if file_dlg[0] != '':
             # User has selected something, get the path, set the window's title to the path
             filename=file_dlg[0]
             # store relevant file information
             self.set_file_info(filename)
         else:
             print "so, you changed your mind, I will do nothing"
             filename = None
             return
         
            
         # get channel number
         chan_num = "%0d"%(self.par["Detector channel"])
         self.stBar1.setText('Current file : %s / Channel : %s'%(self.name+self.ext, chan_num))
         # open file
         print "Open file : ",self.dir + self.name + self.ext
         self.f = h5py.File(self.dir + self.name + self.ext, 'r')
         # get the data
         # time information
         print "Get data"
         data_root = 'wfm_group0/traces/trace' + chan_num + '/'
         try:
             self.t0 = self.f[data_root + 'x-axis'].attrs['start']
             self.dt = self.f[data_root + 'x-axis'].attrs['increment']
             # measured data
             # scale dataset
             self.scale = self.f[data_root + 'y-axis/scale_coef'].value
             # get the y dataset
             self.nall = self.f[data_root + 'y-axis/data_vector/data'].shape[0]
         except:
             mb=QtWidgets.QMessageBox(self)
             mb.setText("Problems loading data " + data_root)
             mb.exec_()
             return
         self.ndata = self.nall
         self.ti = self.t0
         self.tf = self.t0 + (self.ndata-1)*self.dt        
         self.par["tmin"] = self.t0
         self.par["tmax"] = self.t0 + (self.ndata-1)*self.dt
         # WB temp solution 8/13/13
         # set the type as 16 bit signed, this is not good the data type should come from
         # the hdf file
         if convert_int:
              self.ydata = self.f[data_root + 'y-axis/data_vector/data'].value.astype('int16')
         else:
              self.ydata = self.f[data_root + 'y-axis/data_vector/data'].value

         # make the time axis
         print "Calculate data"
         self.tall = self.t0 + self.dt*np.arange(self.ydata.shape[0], dtype = float)
         # select the data to be worked on
         self.f.close()
         print "Select data"
         self.select_data()
         print "Done"
         
     def OpenFile(self, fpname):
         
         self.set_file_info(fpname)
         chan_num = "%0d"%(self.par["Detector channel"])
         self.statusBar.SetStatusText('Current file : %s / Channel : %s'%(self.name+self.ext, chan_num), 0)
         # open file
         print "Open file : ",self.dir + self.name + self.ext
         self.f = h5py.File(self.dir + self.name + self.ext, 'r')
         # get the data
         # time information
         print "Get data"
         data_root = 'wfm_group0/traces/trace' + chan_num + '/'
         try:
             self.t0 = self.f[data_root + 'x-axis'].attrs['start']
             self.dt = self.f[data_root + 'x-axis'].attrs['increment']
             # measured data
             # scale dataset
             self.scale = self.f[data_root + 'y-axis/scale_coef'].value
             # get the y dataset
             self.nall = self.f[data_root + 'y-axis/data_vector/data'].shape[0]
         except:
             mb=QtWidgets.QMessageBox(self)
             mb.setText("Problems loading data " + data_root)
             mb.exec_()
             return
         self.ndata = self.nall
         self.ti = self.t0
         self.tf = self.t0 + (self.ndata-1)*self.dt        
         self.par["tmin"] = self.t0
         self.par["tmax"] = self.t0 + (self.ndata-1)*self.dt
         # WB temp solution 8/13/13
         # set the type as 16 bit signed, this is not good the data type should come from
         # the hdf file
         if convert_int:
              self.ydata = self.f[data_root + 'y-axis/data_vector/data'].value.astype('int16')
         else:
              self.ydata = self.f[data_root + 'y-axis/data_vector/data'].value

         # make the time axis
         print "Calculate data"
         self.tall = self.t0 + self.dt*np.arange(self.ydata.shape[0], dtype = float)
         # select the data to be worked on
         self.f.close()
         print "Select data"
         self.select_data()
         print "Done"    
         

     def OnReload(self):
         # reload data file
         # get channel number
         chan_num = "%0d"%(self.par["Detector channel"])
         self.statusBar.SetStatusText('Current file : %s / Channel : %s'%(self.name+self.ext, chan_num), 0)
         # open file
         print "Open file : ",self.dir + self.name + self.ext
         self.f = h5py.File(self.dir + self.name + self.ext, 'r')
         # get the data
         # time information
         print "Get data"
         data_root = 'wfm_group0/traces/trace' + chan_num + '/'
         try:
             self.t0 = self.f[data_root + 'x-axis'].attrs['start']
             self.dt = self.f[data_root + 'x-axis'].attrs['increment']
             # measured data
             # scale dataset
             self.scale = self.f[data_root + 'y-axis/scale_coef'].value
             # get the y dataset
             self.nall = self.f[data_root + 'y-axis/data_vector/data'].shape[0]
         except:
             mb=QtWidgets.QMessageBox(self)
             mb.setText("Problems loading data " + data_root)
             mb.exec_()
             return
         self.ndata = self.nall
         self.ti = self.t0
         self.tf = self.t0 + (self.ndata-1)*self.dt        
         # self.par["tmax"] = self.t0 + (self.ndata-1)*self.dt
         # WB temp solution 8/13/13
         # set the type as 16 bit signed, this is not good the data type should come from
         # the hdf file
         if convert_int:
              self.ydata = self.f[data_root + 'y-axis/data_vector/data'].value.astype('int16')
         else:
              self.ydata = self.f[data_root + 'y-axis/data_vector/data'].value
         # make the time axis
         print "Calculate data"
         self.tall = self.t0 + self.dt*np.arange(self.ydata.shape[0], dtype = float)
         # select the data to be worked on
         self.f.close()
         print "Select data"
         self.select_data()
         print "Done"
         
     def OnLoadParameters(self,event):
          # Create a file-open dialog in the current directory
          filetypes = '*.data'
          # use a regular file dialog
          file_dlg=QtWidgets.QFileDialog.getOpenFileName(self, 'Select a file', self.dir, filetypes )
          if file_dlg[0] != '':
               # User has selected something, get the path, set the window's title to the path
               filename = file_dlg[0]
               # store relevant file information
               print "will open " , filename
          else:
               print "so, you changed your mind, I will do nothing"
               filename = None
               return
          p = pfile(filename)
          for key in self.par:
               self.par[key] = p.get_value(key)
               print key, " ", self.par[key]
          # now set the check marks
          # get the menus
          menus = self.GetMenuBar().GetMenus()
          items = None
          item_dict = {}
          # search for the options menu and get the contents
          for m,l in menus:
               if l == "Options":
                    items = m.GetMenuItems()
                    for i,item in enumerate(items):
                         if item.GetLabel() == "":
                              continue
                         else:
                              key = item.GetLabel().strip()
                              item_dict[key] = item
          # set the values for the options menu
          if item_dict != {}:
               item_dict["SetLimits"].Check(self.par["limits"]); print "SetLimits = ", self.par["limits"]
               item_dict["Measure"].Check(self.par["measure"]); print "Measure = ", self.par["measure"]
               item_dict["Use Limits"].Check(self.par["use_limits"]); print "Use Limits = ",self.par["use_limits"]
               item_dict["Auto Histogram Limits"].Check(self.par["auto_histo"]); print "Auto Histogram Limits = ", self.par["auto_histo"]
               item_dict["Filtered"].Check(self.par["filtered"]); print "Filtered = ", self.par["filtered"]
               item_dict["Plot Lines"].Check(self.par["draw_lines"]); print "draw_lines = ", self.par["draw_lines"]
          self.select_data()
          # all done
          
     def LoadParameters(self,fpname):
          
          p = pfile(fpname)
          for key in self.par:
               self.par[key] = p.get_value(key)
               print key, " ", self.par[key]
          # now set the check marks
          # get the menus
          menus = self.GetMenuBar().GetMenus()
          items = None
          item_dict = {}
          # search for the options menu and get the contents
          for m,l in menus:
               if l == "Options":
                    items = m.GetMenuItems()
                    for i,item in enumerate(items):
                         if item.GetLabel() == "":
                              continue
                         else:
                              key = item.GetLabel().strip()
                              item_dict[key] = item
          # set the values for the options menu
          if item_dict != {}:
               item_dict["SetLimits"].Check(self.par["limits"]); print "SetLimits = ", self.par["limits"]
               item_dict["Measure"].Check(self.par["measure"]); print "Measure = ", self.par["measure"]
               item_dict["Use Limits"].Check(self.par["use_limits"]); print "Use Limits = ",self.par["use_limits"]
               item_dict["Auto Histogram Limits"].Check(self.par["auto_histo"]); print "Auto Histogram Limits = ", self.par["auto_histo"]
               item_dict["Filtered"].Check(self.par["filtered"]); print "Filtered = ", self.par["filtered"]
               item_dict["Plot Lines"].Check(self.par["draw_lines"]); print "draw_lines = ", self.par["draw_lines"]
          self.select_data()
          # all done
          

     def OnSaveParameters(self,event):
          # Create a file-open dialog in the current directory
          filetypes = '*.data'
          # use a regular file dialog
          if self.par_dir == None:
              self.par_dir = os.getcwd()
          
          file_dlg=QtWidgets.QFileDialog.getSaveFileName(self, 'Select parameter file to save', self.dir,
                                                         filetypes)
          if file_dlg[0] != '':
               # User has selected something, get the path, set the window's title to the path
               filename = file_dlg[0]
               # store relevant file information
               print "will save to " , filename
          else:
               print "so, you changed your mind, I will do nothing"
               filename = None
               return
          o = open(filename, 'w')
          # write the parameters
          for key in self.par:
               o.write( "%s  =  %r\n "%(key, self.par[key] ))
          o.close()
          
     def get_time_window(self, tmin, tmax):
          # find range of time values between tmin and tmax
          nmin = int( (tmin - self.t0)/self.dt )
          nmax = min( (int( (tmax - self.t0)/self.dt ) + 1), self.nall -1 )
          return slice(nmin, nmax + 1)

     def select_data(self):
          # select data range to work on
          # debuggin
          if self.ndata == 0:
               return
          if (self.par["tmin"] < self.tall[0]):
               self.par["tmin"] = self.tall[0]
          if (self.par["tmax"] > self.tall[-1:][0]):
               self.par["tmax"] = self.tall[-1:][0]
          print "Get Window"
          sl = self.get_time_window( self.par["tmin"], self.par["tmax"])
          # the final value
          print "Recalculate"
          self.ndata = sl.stop - sl.start
          self.V = self.scale[0] + self.scale[1]*self.ydata[sl]
          self.t = self.tall[sl]
          
          print "finished recalculation"

     def ChangeCursor(self, event):
          self.figure_canvas.setCursor(QtGui.QCursor(QtCore.Qt.CrossCursor))

     def LineSelectCallback(self, evnt_click, evnt_release):
          'evnt_click and evnt_release are the press and release events'
          start_pos = np.array([evnt_click.xdata, evnt_click.ydata])
          end_pos = np.array([evnt_release.xdata, evnt_release.ydata] )
          self.par["xmin"] = start_pos[0]
          self.par["xmax"] = end_pos[0]
          if self.par["xmin"] > self.par["xmax"]:
               xmax = self.par["xmin"]
               self.par["xmin"] = self.par["xmax"]
               self.par["xmax"] = xmax
          self.ymin = start_pos[1]
          self.ymax = end_pos[1]
          if self.ymin > self.ymax:
               ymax = self.ymin
               self.ymin = self.ymax
               self.ymax = ymax
          dpos = end_pos - start_pos
          print "x positions    : %5.2e --> %5.2e" % (start_pos[0],  end_pos[0] )
          print "y positions    : %5.2e --> %5.2e" % (start_pos[1],  end_pos[1] )
          print "displacement :  delta-x = %5.2e   delta-y  = %5.2e " % tuple(dpos)
          if self.par["measure"]:
               # draw
               xpath = [self.par["xmin"], self.par["xmax"], self.par["xmax"], self.par["xmin"]]
               ypath = [self.ymin, self.ymin, self.ymax, self.ymax]
               self.axes.fill(xpath, ypath, 'b', alpha = 0.5, edgecolor = 'k')
               self.axes.text(self.par["xmin"], self.ymin, "%5.2e, %5.2e"%(start_pos[0], start_pos[1]), 
                              ha = 'right', va = 'top')
               self.axes.text(self.par["xmax"], self.ymax, "%5.2e, %5.2e"%(end_pos[0], end_pos[1]) )
          if self.par["limits"]:
               ymin,ymax = self.axes.get_ylim()
               xpath = [self.par["xmin"], self.par["xmax"], self.par["xmax"], self.par["xmin"]]
               ypath = [ymin, ymin, ymax, ymax]
               self.axes.fill(xpath, ypath, 'g', alpha = 0.5, edgecolor = 'g')
               # self.axes.axvspan(self.par["xmin"], self.par["xmax"], color = 'g', alpha = 0.5)
          self.figure_canvas.draw()

            
     def UpdateStatusBar(self, event):
          if event.inaxes:
               x, y = event.xdata, event.ydata
               self.stBar2.setText( "t= " + "%e"%(x) + "  V= " + "%e"%(y))
     def OnCloseWindow(self):
          if self.histoframe != None:
               try:
                    self.histoframe.destroy()
               except:
                    print "cannot destroy histoframe !"
          if self.tsplotframe != None:
               try:
                    self.tsplotframe.destroy()
               except:
                    print "cannot destroy tsplotframe !"
          
          self.destroy()
          # all done
     #----------------------------------------------------------------------
     # Action menu routines   
     #----------------------------------------------------------------------
     def OnPlot(self):
          V = None
          if self.par["filtered"]:
               V = self.Vinv
          else:
               V = self.V
          if (not self.t.any()) or (not V.any()):
               print "Nothing to plot !"
               return
          if self.par["draw_lines"]:
               #self.axes.plot(self.t[:V.shape[0]], V)
               self.all_plot=self.my_plot(self.t[:V.shape[0]], V)
          else:
               self.all_plot=self.my_plot(self.t[:V.shape[0]], V, '.') 
              #self.axes.plot(self.t[:V.shape[0]], V, '.')
          self.axes.set_xlabel('t')
          self.axes.set_ylabel('V')
          self.figure_canvas.draw()

     def OnClear(self):
          self.figure.clf()
          self.axes = self.figure.add_subplot(111)
          self.RS = RectangleSelector(self.axes, self.LineSelectCallback,
                                      drawtype='box', useblit=True,
                                      button=[1,3], # don't use middle button
                                      minspanx=5, minspany=5,
                                      spancoords='pixels')
          self.figure_canvas.draw()

     def OnFindPeaks(self):
          if (self.t == None) or (self.V == None):
               print "No data, nothing to find !"
               return
          if self.par["filtered"]:
               if self.Vinv == None:
                    print "No inverted data, nothing to find !"
                    return
               V = self.Vinv
          else:
               V = self.V
          print '----------> Start peak finding :'
          self.MAXTAB, self.MINTAB = find_peaks( V, self.par["Vstep"], xval = self.t, \
                                                 xmin = self.par["xmin"], xmax = self.par["xmax"], \
                                                 limits = self.par["use_limits"], \
                                                 get_window = get_window, \
                                                 power = 5)
          print '----------> found ', len(self.MAXTAB[0]), ' peaks'
          # get the values above the V threshold
          tp = np.array(self.MAXTAB[0])
          Vp = np.array(self.MAXTAB[1])
          iw = np.where( Vp > self.par["Vthreshold"] )[0]
          print '----------> store ', len(iw), ' peaks above threshold'
          self.t_peak = tp[iw]
          self.V_peak = Vp[iw]
          if self.par["auto_histo"] and not len(iw) == 0:
               self.par["Vhmin"] = self.V_peak.min()
               self.par["Vhmax"] = self.V_peak.max()
          
          self.axes.set_autoscaley_on(True)
          self.axes.plot(self.t_peak,self.V_peak,'r.')
          self.figure_canvas.draw()

     def OnHistogram(self):
          if self.MAXTAB == None:
               print "No peaks yet !"
               return
          histo = LT.box.histo( self.V_peak, \
                                     range = (self.par["Vhmin"], self.par["Vhmax"]), \
                                     bins = self.par["VNbins"])
          if histo == None:
               print "No histogram to plot !"
               return
          histo.title = 'V histogram'
          histo.xlabel = 'Volts'
          histo.ylabel = 'Counts'
          # add the histograms to the list of histos created
          self.histos.append( histo )
          # setup figure for hist plots
          latest_histo = len(self.histos) - 1
          self.show_histos(latest_histo)
          #self.histo.plot(axes = self.histoframe.axes)
          #self.histoframe.figure_canvas.draw()

     def OnDeleteHistogram(self,event):
          self.histos = []
          self.current_histo = 0
          if self.histoframe != None:
               # clear figure
               self.histoframe.figure.clf()
               self.histoframe.axes = self.histoframe.figure.add_subplot(111)

     def OnSaveHistos(self):
          if self.histos == []:
               print "no histograms to save !"
               return
          # Create a file-open dialog in the current directory
          filetypes="*.data"
          # use a regular file dialog
          if self.hist_dir == None:
              self.hist_dir = os.getcwd()
          file_dlg=QtWidgets.QFileDialog.getSaveFileName(self, 'Select histogram file to save', self.hist_dir,
                                                         filetypes)
          if file_dlg[0] != '':
               # User has selected something, get the path, set the window's title to the path
               filename = file_dlg[0]
               # analyze the file name
               dir, fname = os.path.split(filename)
               name, ext = os.path.splitext(fname)
               # store relevant file information
               N = len(self.histos)-1
               print "will save to : "+dir+"/"+name+"_0.data ... " +name+"_{0}.data".format(N) 
          else:
               print "so, you changed your mind, I will do nothing"
               filename = None
               return
          for i,h in enumerate(self.histos):
               h_file = dir + '/'+ name +'_{0}.data'.format(i)
               h.save(h_file)
          print " all histograms saved"

     def show_histos(self, n):
          if self.histos == []:
               print "no histograms to show !"
               return
          # create it's own frame if it doesnot exist already
          self.current_histo = n
          if self.histoframe == None:
               main_rect = self.GetScreenRect()
               x,y,width,height = main_rect.Get()
##               self.histoframe = HistoPlotFrame(self, \
#                                                pos = (x + width + 10 , y), \
#                                                size = (width, height), \
#                                                callback = self)
#               self.histoframe.Show(True)
          else:
               # clear figure
               self.histoframe.figure.clf()
               self.histoframe.axes = self.histoframe.figure.add_subplot(111)
          self.histos[n].plot(axes = self.histoframe.axes)
          self.histoframe.figure_canvas.draw()
          

     def show_slice(self, n):
          if self.t_slice== None:
               print "no slices to show !"
               return
          # create it's own frame if it doesnot exist already
          self.current_slice = n
          if self.tsplotframe == None:
               main_rect = self.GetScreenRect()
               x,y,width,height = main_rect.Get()
##               self.tsplotframe = TSPlotFrame(self, \
#                                                pos = (x + width + 10 , y), \
#                                                size = (width, height), \
#                                                callback = self)
#               self.tsplotframe.Show(True)
          else:
               # clear figure
               self.tsplotframe.figure.clf()
               self.tsplotframe.axes = self.tsplotframe.figure.add_subplot(111)
          self.tsplot(n)
          self.tsplotframe.figure_canvas.draw()
          

     def OnTScalc(self,event):
          if (self.t == None) or (self.V == None):
               print "No data, nothing to find !"
               return
          # make sure the limits are within the data
          if (self.par["ts_start"] > self.par["ts_stop"]):
               stop = self.par["ts_start"]
               self.par["ts_start"] = self.par["ts_stop"]
               self.par["ts_stop"] = stop
          self.par["ts_start"] = min( self.par["ts_start"], self.par["tmax"]) 
          self.par["ts_stop"] = max( self.par["ts_stop"], self.par["tmin"])
          self.par["ts_start"] = max( self.par["ts_start"], self.par["tmin"])
          self.par["ts_stop"] = min( self.par["ts_stop"], self.par["tmax"])
          #
          # get the first slice
          sl = get_window(self.par["ts_start"], self.t, self.par["ts_stop"] )
          ns = int((self.t[sl].max() - self.t[sl].min())/self.par["ts_width"])
          # reshape the data
          #division by zero bug (not clear)
          if ns != 0:
              new_shape = (ns, int(self.t[sl].shape[0]/ns) )
          else:
              new_shape = (ns, int(self.t[sl].shape[0]) )
          # the new data have the shape (slice_nr, data_in_slice)
          # axis-0 slice number
          # axis-1 data within slice 
          print 'resizing time data array, be patient !'
          self.t_slice = np.resize(self.t[sl], new_shape )
          if self.par["filtered"]:
               if self.Vinv == None:
                    print "No inverted data, nothing to find !"
                    return
               V = self.Vinv
          else:
               V = self.V
          print 'resizing Signal data array, be patient !'
          self.V_slice = np.resize(V, new_shape )
          # now the data are ready to be worked on
          self.ts_av = np.average(self.t_slice, axis = 1)
          print "created ", len(self.ts_av), " slices"

     def tsplot(self, n):
          if self.V_slice == None:
               print "No time slices!"
               return
          V = self.V_slice[n]
          t = self.t_slice[n]
          # show the data
          if self.par["draw_lines"]:
               self.tsplotframe.axes.plot(t[:V.shape[0]], V)
          else:
               self.tsplotframe.axes.plot(t[:V.shape[0]], V, '.')
          # show peka positions if there are any
          if self.ts_counts != None:
               if self.ts_counts[n] != 0:
                    self.tsplotframe.axes.plot(self.ts_tp[n], self.ts_Vp[n], 'ro')
          self.tsplotframe.axes.set_xlabel('t')
          self.tsplotframe.axes.set_ylabel('V')
          s_title = 'Slice {0:d}, T_av = {1:10.3e}'.format(n, self.ts_av[n])
          self.tsplotframe.axes.set_title(s_title)
          self.tsplotframe.figure_canvas.draw()

 
     def OnTSfindpeaks(self):
          if self.V_slice == None:
               print "No time slices!"
               return
          ts_tp = []
          ts_Vp = []
          ts_counts = []
          for i,V in enumerate(self.V_slice):
               print '----------> Start peak finding in slice :', i
               MAXTAB, MINTAB = find_peaks( V, self.par["Vstep"], xval = self.t_slice[i], \
                                                      limits = False, \
                                                      get_window = get_window, \
                                                      power = 5)
               print '----------> found ', len(MAXTAB[0]), ' peaks'
               # get the values above the V threshold
               tp = np.array(MAXTAB[0])
               Vp = np.array(MAXTAB[1])
               iw = np.where( Vp > self.par["ts_Vthreshold"] )[0]
               print '----------> store ', len(iw), ' peaks above threshold'
               ts_counts.append(len(iw))
               ts_tp.append(tp[iw])
               ts_Vp.append(Vp[iw])
          self.ts_tp = ts_tp
          self.ts_Vp = ts_Vp
          self.ts_counts = np.array(ts_counts)


     def OnTSplotrate(self):
          if (self.ts_av == None) or (self.ts_counts == None):
               print "Nothing to plot !"
               return
          rate = self.ts_counts/self.par["ts_width"]
          rate_err = np.sqrt(self.ts_counts)/self.par["ts_width"]
          self.axes.errorbar(self.ts_av, rate, yerr = rate_err, marker='o', ls = 'None')
          self.axes.set_xlabel('t')
          self.axes.set_ylabel('Rate')
          self.figure_canvas.draw()
         


     def OnTSsaverate(self):
          if (self.ts_av == None) or (self.ts_counts == None):
              print "Nothing to Save !"
              return
          # calculate rate
          rate = self.ts_counts/self.par["ts_width"]
          rate_err = np.sqrt(self.ts_counts)/self.par["ts_width"]
          # Create a file-open dialog in the current directory
          filetypes="*.data"
          # use a regular file dialog
          if self.res_dir == None:
              self.res_dir = os.getcwd()
          file_dlg=QtWidgets.QFileDialog.getSaveFileName(self, 'Select rate file to save', self.res_dir,
                                                         filetypes)
          if file_dlg[0] !='':
               # User has selected something, get the path, set the window's title to the path
               filename = file_dlg
               # analyze the file name
               dir, fname = os.path.split(filename)
               name, ext = os.path.splitext(fname)
               print "will save to : ", filename
          else:
               print "so, you changed your mind, I will do nothing"
               filename = None
               return
          # write the rates into a data files
          o = open(filename, 'w')
          o.write("# rates calculated from: " + self.dir + self.name + self.ext + "\n")
          o.write("#\n")
          o.write("#! time[0,f]/ rate[f,1]/ drate[f,2]/\n")
          for i,t in enumerate(self.ts_av):
              o.write( "{} {} {} \n".format(t, rate[i], rate_err[i]) )
          o.close()
          # all done    
         
     def OnTShistogram(self):
          if (self.ts_Vp == None):
               print "Nothing to histogram !"
               return
          # clear histos list
          self.histos = []
          for i,Vp in enumerate(self.ts_Vp):
               h = LT.box.histo( Vp, \
                                 range = (self.par["ts_Vhmin"], self.par["ts_Vhmax"]), \
                                 bins = self.par["ts_VNbins"])
               h.title = "V histogram # %d from: %6.3e to: %6.3e (s)"%(i, self.ts_tp[i][0], self.ts_tp[i][-1])
               h.xlabel = 'Volts'
               h.ylabel = 'Counts'
               self.histos.append(h)
          # hisograms created
          print "Time slice histograms completed"

     def OnTSshowhistos(self):
          if self.histos == []:
               print "No TS histos"
               return
          # show ths histograms, start with the first one
          self.show_histos(0)

     def OnShow_next_histo(self):
          n = self.current_histo + 1
          if n >= len(self.histos):
               print "at the last histogram"
               n = len(self.histos)-1
          self.show_histos(n)
          
     def OnShow_prev_histo(self):
          n = self.current_histo - 1
          if n < 0:
               print "at the first histogram"
               n = 0
          self.show_histos(n)

     def OnTSshowslice(self):
          if self.ts_av == None:
               print "No slices"
               return
          # show the histograms, start with the first one
          self.show_slice(0)

     def OnShow_next_slice(self):
          n = self.current_slice + 1
          if n >= len(self.ts_av):
               print "at the last slice"
               n = len(self.ts_av)-1
          self.show_slice(n)
          
     def OnShow_prev_slice(self):
          n = self.current_slice - 1
          if n < 0:
               print "at the first slice"
               n = 0
          self.show_slice(n)
          
     #---------------------------------------------------------------------
     # FFT on time slices     
     #---------------------------------------------------------------------
     def OnTSFFTcalc(self,event):
          if self.ts_av == None:
               print "No slices"
               print "No data, nothing to transform !"
               return
          # calculate FFT and PS
          self.ts_FFT = []
          for i,V in enumerate(self.V_slice):
              t = self.t_slice[i]
 #             V = self.V_slice[i]
              print "Calculate FFT for slice: ", i, " and ", len(t), " data points"
              ffts = FFT.FFT(t)   
              ffts.transform(self.V)
              ffts.get_ps()
              ffts.logp = np.log10(ffts.p)
              # self.ts_FFT.append(ffts)
          print "all TS FFT completed" 

     def OnTSFFTplotps(self,event):
          # plot power spectrum
          if not self.fft.ok:
               print "No power spectrum to plot"
               return
          if self.par["draw_lines"]:
               self.axes.plot(self.fft.f,self.fft.logp)
          else:
               self.axes.plot(self.fft.f,self.fft.logp, '.')
          # self.axes.set_yscale('log')
          self.axes.set_xlabel('f')
          self.axes.set_ylabel('Power')
          self.axes.set_title('log(Power)')
          self.figure_canvas.draw()
          

     #----------------------------------------------------------------------
     # Parameter routines
     #----------------------------------------------------------------------
     def OnSelectChannel(self):
	  # show and change parameters
	  # parameter keys
	  pkeys=['Detector channel']
	  data = {pkeys[0]:"%d"%(self.par["Detector channel"])}
          pdlg = NumberDialog(data, self, title="Detector Channel", \
                              labels = pkeys, \
                              keys = pkeys, \
                              about_txt = "Select the data channel")
          pdlg.exec_()
          # now set the new parameters
          self.par['Detector channel']=int(pdlg.data[pkeys[0]])
          pdlg.destroy()

     def OnSelectTimeSlot(self):
          # show and change the parameters
          # parameter keys
          pkeys=['Number of Data','t start', 't stop', 'delta t']
          # current data to be shown in the dialog
          data = {pkeys[0]:"%d"%(self.ndata), \
                  pkeys[1]:"%5.2e"%(self.par["tmin"]), \
                  pkeys[2]:"%5.2e"%(self.par["tmax"]), \
                  pkeys[3]:"%5.2e"%(self.dt), \
                  }
          pdlg = NumberDialog(data, self, title="Time axis", \
                              labels = pkeys, \
                              keys = pkeys, \
                              about_txt = "Select the time slot")
          pdlg.exec_()
          # now set the new parameters
          self.par["tmin"] = float(pdlg.data[pkeys[1]])
          self.par["tmax"] = float(pdlg.data[pkeys[2]])
          pdlg.destroy()
          
          self.select_data()

     def OnSetPeakpar(self):
          # show and change the parameters
          # parameter keys
          pkeys=['t min', 't max','delta V','threshold']
          # current data to be shown in the dialog
          data = {pkeys[0]:"%5.2e"%(self.par["xmin"]), \
                  pkeys[1]:"%5.2e"%(self.par["xmax"]), \
                  pkeys[2]:"%5.3f"%(self.par["Vstep"]), \
                  pkeys[3]:"%5.3f"%(self.par["Vthreshold"]) \
                  }
          pdlg = NumberDialog(data, self, title="Parameters for Peak Finding", \
                              labels = pkeys, \
                              keys = pkeys, \
                              about_txt = "Peak finding")
          pdlg.exec_()
          # now set the new parameters
          self.par["xmin"] = float(pdlg.data[pkeys[0]])
          self.par["xmax"] = float(pdlg.data[pkeys[1]])
          self.par["Vstep"] = float(pdlg.data[pkeys[2]])
          self.par["Vthreshold"] = float(pdlg.data[pkeys[3]])
          pdlg.destroy()
          
     def OnSetHisto(self):
          # show and change the parameters
          # parameter keys
          pkeys=['V min','V max', 'N bins']
          # current data to be shown in the dialog
          data = {pkeys[0]:"%5.3f"%(self.par["Vhmin"]), \
                  pkeys[1]:"%5.3f"%(self.par["Vhmax"]), \
                  pkeys[2]:"%d"%(self.par["VNbins"])}
          
          pdlg = NumberDialog(data, self, title="Histogram", \
                              labels = pkeys, \
                              keys = pkeys, \
                              about_txt = "Current parameters ")
          pdlg.exec_()
          # now set the new parameters
          self.par["Vhmin"] = float(pdlg.data[pkeys[0]])
          self.par["Vhmax"] = float(pdlg.data[pkeys[1]])
          self.par["VNbins"] = int(pdlg.data[pkeys[2]])
          pdlg.destroy()

     def OnSetTimeSlice(self):
          # show and change the parameters
          # parameter keys
          pkeys=['t start','t stop','t width', 'V min', 'V max', 'N bins', 'V threshold']
          # current data to be shown in the dialog
          data = {pkeys[0]:"%6.2e"%(self.par["ts_start"]), \
                  pkeys[1]:"%6.2e"%(self.par["ts_stop"]), \
                  pkeys[2]:"%6.2e"%(self.par["ts_width"]), \
                  pkeys[3]:"%5.3f"%(self.par["ts_Vhmin"]), \
                  pkeys[4]:"%5.3f"%(self.par["ts_Vhmax"]), \
                  pkeys[5]:"%d"%(self.par["ts_VNbins"]), \
                  pkeys[6]:"%5.3f"%(self.par["ts_Vthreshold"]), \
                  }
          pdlg = NumberDialog(data, self, title="Time Slicing", \
                              labels = pkeys, \
                              keys = pkeys, \
                              about_txt = "Current parameters ")
          pdlg.exec_()
          # now set the new parameters
          self.par["ts_start"] = float(pdlg.data[pkeys[0] ])
          self.par["ts_stop"] = float(pdlg.data[pkeys[1] ])
          self.par["ts_width"] = float(pdlg.data[pkeys[2] ])
          self.par["ts_Vhmin"] = float(pdlg.data[pkeys[3] ])
          self.par["ts_Vhmax"] = float(pdlg.data[pkeys[4] ])
          self.par["ts_VNbins"] = int(pdlg.data[pkeys[5] ])
          self.par["ts_Vthreshold"] = float(pdlg.data[pkeys[6] ])
          pdlg.destroy()

     #----------------------------------------------------------------------
     # FFT routines
     #----------------------------------------------------------------------
     def OnFFTcalc(self):
          if (self.t == None) or (self.V == None):
               print "No data, nothing to transform !"
               return
          # calculate FFT and PS
          print "Calculate FFT for : ", len(self.t), " data points"
          self.fft = FFT.FFT(self.t)
          self.fft.transform(self.V)
          self.fft.get_ps()
          self.fft.logp = np.log10(self.fft.p)

     def OnFFTplotps(self):
          # plot power spectrum
          if not self.fft.ok:
               print "No power spectrum to plot"
               return
          if self.par["draw_lines"]:
               self.axes.plot(self.fft.f,self.fft.logp)
          else:
               self.axes.plot(self.fft.f,self.fft.logp, '.')
          # self.axes.set_yscale('log')
          self.axes.set_xlabel('f')
          self.axes.set_ylabel('Power')
          self.axes.set_title('log(Power)')
          self.figure_canvas.draw()

     def OnFFTplotan(self):
          # plot sin coeff
          if not self.fft.ok:
               print "No power spectrum to plot"
               return
          if self.par["draw_lines"]:
               self.axes.plot(self.fft.f,np.abs(self.fft.an))
          else:
               self.axes.plot(self.fft.f,np.abs(self.fft.an),'.')
          self.axes.set_yscale('log')
          self.axes.set_xlabel('f')
          self.axes.set_ylabel('an')
          self.figure_canvas.draw()

     def OnFFTplotbn(self):
          # plot cos coeff
          if not self.fft.ok:
               print "No power spectrum to plot"
               return
          if self.par["draw_lines"]:
               self.axes.plot(self.fft.f,np.abs(self.fft.bn))
          else:
               self.axes.plot(self.fft.f,np.abs(self.fft.bn),'.')
          self.axes.set_yscale('log')
          self.axes.set_xlabel('f')
          self.axes.set_ylabel('abs(bn)')
          self.figure_canvas.draw()

     def OnSetFFTfilterpar(self):
          # show and change the parameters
          if not self.fft.ok:
               print "No FFT"
               return
          # parameter keys
          pkeys=['f min', 'f max', 'alpha ', 'r - threshold']
          # current data to be shown in the dialog
          data = {pkeys[0]:"%5.2e"%(self.fft.par["xmin"]), \
                  pkeys[1]:"%5.2e"%(self.fft.par["xmax"]), \
                  pkeys[2]:"%5.2e"%(self.fft.par["alpha"]), \
                  pkeys[3]:"%5.2e"%(self.fft.par["rthresh"]), \
                  }
          pdlg = NumberDialog(data, self, title="FFT Filter Settings", \
                              labels = pkeys, \
                              keys = pkeys, \
                              about_txt = "Parameters")
          pdlg.exec_()
          # now set the new parameters
          self.fft.par["xmin"] = float(pdlg.data[pkeys[0]])
          self.fft.par["xmax"] = float(pdlg.data[pkeys[1]])
          self.fft.par["alpha"] = float(pdlg.data[pkeys[2]])
          self.fft.par["rthresh"] = float(pdlg.data[pkeys[3]])
          pdlg.destroy()

     def OnFFTcutfilter(self):
          # Apply freq. cut to FFT
          if not self.fft.ok:
               print "No FFT"
               return
          self.fft.an_cut_freq(replace = True)
          self.fft.bn_cut_freq(replace = True)
          self.fft.store_fft_coeff(self.fft.an, self.fft.bn)
          self.fft.get_ps()
          self.fft.logp = np.log10(self.fft.p)

     def OnFFTrfilter(self):
          # Apply r-cut to FFT
          if not self.fft.ok:
               print "No FFT"
               return
          self.fft.r_cut_freq(replace = True)
          self.fft.store_fft_coeff(self.fft.an, self.fft.bn)
          self.fft.get_ps()
          self.fft.logp = np.log10(self.fft.p)
     
     def OnFFTinvert(self):
          print "Invert FFT"
          if not self.fft.ok:
               print "No FFT"
               return
          self.Vinv = self.fft.inv_transform()
          print "FFT inverted"
     #----------------------------------------------------------------------
     # options routines
     #----------------------------------------------------------------------

     def OnMeasure(self):
          # to measure create a Rectangle Selector
          self.RS.set_active(True)
          self.par["measure"] = True
          self.par["limits"] = False
          print "limits = ", self.par["limits"], " measure = ", self.par["measure"]

     def OnLimits(self):
          # to measure create a Rectangle Selector
          self.RS.set_active(True)
          self.par["measure"] = False
          self.par["limits"] = True
          print "limits = ", self.par["limits"], " measure = ", self.par["measure"]

     def OnTogglelimits(self, event):
          menubar = self.GetMenuBar()
          itemId = event.GetId()
          item = menubar.FindItemById(itemId)
          self.par["use_limits"] = item.IsChecked()
          print "use limits = ", self.par["use_limits"]

     def OnToggleHistolimits(self,event):
          menubar = self.GetMenuBar()
          itemId = event.GetId()
          item = menubar.FindItemById(itemId)
          self.par["auto_histo"] = item.IsChecked()
          print "auto histo = ", self.par["auto_histo"]

     def OnToggleLines(self, event):
          menubar = self.GetMenuBar()
          itemId = event.GetId()
          item = menubar.FindItemById(itemId)
          self.par["draw_lines"] = item.IsChecked()
          print "draw Lines = ", self.par["draw_lines"]

     def OnUsefiltered(self, event):
          menubar = self.GetMenuBar()
          itemId = event.GetId()
          item = menubar.FindItemById(itemId)
          self.par["filtered"] = item.IsChecked()
          print "use filtered = ", self.par["filtered"]
          
     def OnNothing(self):
          print "do nothing"