from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5 import NavigationToolbar2QT
from matplotlib.widgets import RectangleSelector
import matplotlib.pyplot as pl 
from matplotlib.figure import Figure
# for GUI stuff
from PyQt5 import QtWidgets, QtCore, QtGui
# HDF file handling
import h5py
import sys
import os
#import time
import numpy as np
# general LabTools
import LT.box
from LT.parameterfile import pfile
# FFT
import FFT
#Navigation toolbar and number dialog
from NavToolBar_NumberDial import NavigationToolbar, NumberDialog

convert_int = True
colors = ['red', 'green', 'blue', 'magenta', 'cyan', 'orange',
          'lavenderblush', 'maroon', 'plum']
import ffind_peaks as FP

#----------------------------------------------------------------------
# useful functions
#----------------------------------------------------------------------

         
def get_window(xmin, x, xmax):
     # get the slice corresponding to xmin and xmax in x
     # the x-values need to be equally spaced
     dx = x[1] - x[0]
     nmin = max( 0, int( (xmin - x[0])/dx ))
     nmax = min( (int( (xmax - x[0])/dx ) + 1), len(x) -1 )
     return slice(nmin, nmax + 1)
     
     
def find_peaks(yval, ystep, xval = None, \
               xmin = None, xmax = None, limits = None, \
               power = 5,\
               get_window = None ):
     # find peaks in an array of data
     # get_window is a function that returns the slice of data between xmin and xmax
     # peak finding
     nmin = 0
     nmax = 0
     MAXTAB=[]
     MINTAB=[]
     if ((xmin == None) and (xmax == None)) or (not limits) or (get_window == None):
          print("No limits present, analyze all data")
          results = np.zeros((2,), dtype = 'int32')
          pmin = np.zeros((len(yval)//5, ), dtype='int32')
          pmax = np.zeros((len(yval)//5, ), dtype='int32')
          try:
              FP.find_peaks(len(yval), ystep, yval, results, pmin, pmax)
          except:
              print("problem with peak finding")
              return []
          nmin = results[0]
          nmax = results[1]
     else:
          # get the window
          print(("Analyze data between ", xmin, " and ", xmax))
          sl = get_window( xmin, xval, xmax)
          # progress dialog
          results = np.zeros((2,), dtype = 'int32')
          pmin = np.zeros((len(yval[sl])/5, ), dtype='int32')
          pmax = np.zeros((len(yval[sl])/5, ), dtype='int32')
          try:
              FP.find_peaks(len(yval[sl]), ystep, yval[sl], results, pmin, pmax)
          except:
               print("problem with peak finding")
               return []
          nmin = results[0]
          nmax = results[1]
     MAXTAB.append( xval[pmax[:nmax]] )
     MAXTAB.append( yval[pmax[:nmax]] )
     MAXTAB.append(pmax[:nmax])
     MINTAB.append( xval[pmin[:nmin]] )
     MINTAB.append( yval[pmin[:nmin]] )
     MINTAB.append(pmin[:nmin])
     return [MAXTAB,MINTAB]
     






#----------------------------------------------------------------------
# plot frame, for histograms
#----------------------------------------------------------------------


class HistoPlotFrame(QtWidgets.QMainWindow):
     def __init__(self, parent, callback = None, **kwargs):
          QtWidgets.QMainWindow.__init__(self)
          self.setWindowTitle('Histogram Plot')
          self.callback = callback
          self.parent=parent
          #Change background color
          self.palette=QtGui.QPalette()
          self.palette.setColor(QtGui.QPalette.Background,QtCore.Qt.white)
          
          # create status bar
          self.stBar1 = QtWidgets.QLabel('Current histogram')
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
          
          # put figure cancas in center of Main Window
          self.setCentralWidget(self.figure_canvas)
          
          
          # Note that event is a MplEvent
          self.figure_canvas.mpl_connect('motion_notify_event', self.UpdateStatusBar)
          self.figure_canvas.mpl_connect('figure_enter_event', self.ChangeCursor)
          #
          
          # add ToolBar
          self.toolbar = NavigationToolbar2QT(self.figure_canvas, self)
#          self.toolbar.removeAction(self.toolbar.Ncontrol)
#          self.toolbar.removeAction(self.toolbar.Nlabel)
#          self.toolbar.removeAction(self.toolbar.refr)
          self.addToolBar(self.toolbar)
          
          
          self.setPalette(self.palette)
          geo=self.parent.geometry()
          geo.adjust(20,20,20,20)
          self.setGeometry(geo)
          
          # add tools for histogram selection
          self.hist_back=QtWidgets.QAction(QtGui.QIcon("l_arrow.png"),'Previous histogram', self)
          self.toolbar.addAction(self.hist_back)
          self.hist_back.triggered.connect(self.OnPrevhisto)
          
          self.hist_forward=QtWidgets.QAction(QtGui.QIcon("r_arrow.png"),'Next histogram', self)
          self.toolbar.addAction(self.hist_forward)
          self.hist_forward.triggered.connect(self.OnNexthisto)
          
          self.show()

     def UpdateStatusBar(self, event):
          if event.inaxes:
               x, y = event.xdata, event.ydata
               self.stBar2.setText( "x= " + "%e"%(x) + "  Counts= " + "%e"%(y))
    
     def ChangeCursor(self, event):
          self.figure_canvas.setCursor(QtGui.QCursor(QtCore.Qt.CrossCursor))

     def OnPrevhisto(self):
          if self.callback == None:
               print("Prev button clicked !")
          else:
               self.callback.OnShow_prev_histo()
               

     def OnNexthisto(self):
          if self.callback == None:
               print("Next button clicked !")
          else:
               self.callback.OnShow_next_histo()

     def closeEvent(self, Event):
          
          print("Close histo frame")
          self.callback.histoframe=None
          self.destroy()
#----------------------------------------------------------------------
# plot frame, for time slices
#----------------------------------------------------------------------

class TSPlotFrame(QtWidgets.QMainWindow):
     def __init__(self, parent, callback = None, **kwargs):
          QtWidgets.QMainWindow.__init__(self)
          self.setWindowTitle('Time Slice Plot')
          self.callback = callback
          self.parent=parent
          self.par = self.parent.par
          #Change background color
          self.palette=QtGui.QPalette()
          self.palette.setColor(QtGui.QPalette.Background,QtCore.Qt.white)
          
          # create status bar
          self.stBar1 = QtWidgets.QLabel('Current slice')
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
          
          # put figure cancas in center of Main Window
          self.setCentralWidget(self.figure_canvas)
          
          # Note that event is a MplEvent
          self.figure_canvas.mpl_connect('motion_notify_event', self.UpdateStatusBar)
          self.figure_canvas.mpl_connect('figure_enter_event', self.ChangeCursor)
          #
          # add ToolBar
          self.toolbar = NavigationToolbar(self.figure_canvas, self)
          self.addToolBar(self.toolbar)
          
          self.setPalette(self.palette)
          geo=self.parent.geometry()
          geo.adjust(20,20,20,20)
          self.setGeometry(geo)
          
          # add tools for histogram selection
          self.hist_back=QtWidgets.QAction(QtGui.QIcon("l_arrow.png"),'Previous slice', self)
          self.toolbar.addAction(self.hist_back)
          self.hist_back.triggered.connect(self.OnPrevslice)
          
          self.hist_forward=QtWidgets.QAction(QtGui.QIcon("r_arrow.png"),'Next slice', self)
          self.toolbar.addAction(self.hist_forward)
          self.hist_forward.triggered.connect(self.OnNextslice)
          
          self.show() 
    
     def UpdateStatusBar(self, event):
          if event.inaxes:
               x, y = event.xdata, event.ydata
               self.stBar2.setText( "t = " + "%e"%(x) + " V = " + "%e"%(y))
    
     def ChangeCursor(self, event):
          self.figure_canvas.setCursor(QtGui.QCursor(QtCore.Qt.CrossCursor))

     

     def closeEvent(self, Event):
          
          print("Close timeslice frame")
          self.callback.tsplotframe=None
          self.destroy()
     
     
     def OnPrevslice(self):
          if self.callback == None:
               print("Prev button clicked !")
          else:
               self.callback.OnShow_prev_slice()

     def OnNextslice(self):
          if self.callback == None:
               print("Next button clicked !")
          else:
               self.callback.OnShow_next_slice()
     
     #tsplotframe plot
     def my_plot(self, *args, **kwargs):
        
        N=kwargs.pop('N', self.toolbar.N)
        ax=self.axes
        #saves input data into toolbar class object and reuses later
        self.toolbar.myinitdata = args
        
        #dropping points for first plot
        t=args[0]
        V=args[1]
        rest=args[2:]
        n=1 
        #take every nth point from data
        if len(t)> N:
            n=int(len(t)/N)
        tcut=t[::n]
        Vcut=V[::n]
      
        ax.autoscale(enable=True, axis='x', tight=False)
        
     
        return ax.plot(tcut, Vcut, *rest, **kwargs)      
        

#--
## Classes for processing files selection
#class Repository(QtWidgets.QDialog):
#    def __init__(self, parent = None, title = 'Choose files to process', files=None):
#          
#        QtWidgets.QDialog.__init__(self, parent)
#        self.parent=parent
#        self.setWindowFlags(self.windowFlags() ^ QtCore.Qt.WindowContextHelpButtonHint)
#        self.setWindowTitle(title)
#        self.layout=QtWidgets.QGridLayout(self)
#        
#        ok=QtWidgets.QPushButton('Ok', self)
#        ok.setDefault(True)
#        cancel=QtWidgets.QPushButton('Cancel', self)
#        ok.clicked.connect(self.OnOk)
#        cancel.clicked.connect(self.OnCancel)
#        selall=QtWidgets.QPushButton('Select All', self)
#        selall.clicked.connect(self.OnSelectAll)
#        desall=QtWidgets.QPushButton('Deselect All', self)
#        desall.clicked.connect(self.OnDeselectAll)
#          
#        self.layout.addWidget(selall, 0, 0)
#        self.layout.addWidget(desall, 0, 1)                
#        self.layout.addWidget(ok, 3, 0)
#        self.layout.addWidget(cancel, 3, 1)
##          
#        self.files=files
#        self.parent=parent
#        
#        self.model = QtGui.QStandardItemModel()
#        
#        
#        for i in files:
#            item = QtGui.QStandardItem(i)
#            item.setCheckable(True)
#            self.model.appendRow(item)
#        
#        view=QtWidgets.QListView(parent)
#        view.setModel(self.model)
#        self.layout.addWidget(view, 1, 0, 2, 2)
#        self.show()
##
#    def OnSelectAll(self):
#        for index in range(self.model.rowCount()):
#            item = self.model.item(index)
#            if item.isCheckable() and item.checkState() == QtCore.Qt.Unchecked:
#                item.setCheckState(QtCore.Qt.Checked)
#
#    def OnDeselectAll(self, event):
#       for index in range(self.model.rowCount()):
#            item = self.model.item(index)
#            if item.isCheckable() and item.checkState() == QtCore.Qt.Checked:
#                item.setCheckState(QtCore.Qt.Unchecked)
#
#    def OnOk(self):
#        for index in range(self.model.rowCount()):
#            item = self.model.item(index)
#            if item.isCheckable() and item.checkState() == QtCore.Qt.Checked:
#                self.parent.hwstodo.append(self.files[index])
#        self.hide()
#        self.parent.Process()
#        self.destroy()
#        
#    def OnCancel(self):
#         self.destroy()

#----------------------------------------------------------------------
# master frame, contains the plot image
#----------------------------------------------------------------------
class PlotFrame(QtWidgets.QMainWindow):
     def __init__(self, parent):
          super(PlotFrame, self).__init__()
          #QtWidgets.QMainWindow.__init__(self)
         
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
          self.par["Vthreshold"] = 0.3
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
          
          self.datadir = '../Raw_Data/'
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
#          #load list of processed files
#          self.proclistload() 
          
          # file information data
          self.parent = parent
          self.full_name = None
          self.dir = '../Raw_Data/' #None
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
          self.figure_canvas.mpl_connect('motion_notify_event', self.UpdateStatusBar)
          self.figure_canvas.mpl_connect('figure_enter_event', self.ChangeCursor)
          
          self.all_plot=[]
         
          
          #self.setCentralWidget(self.figure_canvas)
          
          
          
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

        
#     def proclistload(self): #load list of processed files
#        path_to_watch = self.datadir
#        dircont = os.listdir (path_to_watch)
#                
#        if "ProcFiles.data" in dircont:
#            print "\nLoad list of already processed files."
#            fopen = open(self.datadir+"//ProcFiles.data")
#            self.hwsproc=[x.strip('\n') for x in fopen.readlines()]
#            
#        else: 
#            self.hwsproc=[]
#            print "\nNo files were processed yet in current data directory."
                
#     def proclistsave(self): #save list of processed files
#        path_to_watch = self.datadir
#        dircont = os.listdir (path_to_watch)
#                
#        if "ProcFiles.data" in dircont:
#            print "\nAppend list of processed files."
#            fopen = open(self.datadir+"//ProcFiles.data", 'a')
#            for item in self.hwsproc:
#                fopen.write("%s\n" % item)
#            fopen.close()
#        else: 
#            fopen = open(self.datadir+"//ProcFiles.data", 'w')
#            for item in self.hwsproc:
#                fopen.write("%s\n" % item)
#            fopen.close()
#            print "\nCreate list of processed files."
         
         
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
#                ("&Scan data directory", "Scan directory", self.OnScan),
#                (None, None, None),  # creates a separator bar in the menu
                ("&Plot", "Plot data", self.OnPlot),
                (None, None, None),  # creates a separator bar in the menu
                ("&FindPeaks", "find peaks in the plotting data", self.OnFindPeaks),
                (None, None, None),
                ("&Histogram", "histogram peak data", self.OnHistogram),
                ("&Delete Histogram", "delete all histograms", self.OnDeleteHistogram),
                (None, None, None),
                ("&Clear Figure", "Clear figure", self.OnClear)),
               #
               ("&Parameters",
#                ("&Data directory", "Set data directory", self.OnSetScanDir), 
#                (None, None, None),
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
                            
                   
                   if 'Measure' in eachItem[0]: subMenu.setChecked(self.par['measure'])
                   if 'SetLimits' in eachItem[0]: subMenu.setChecked(self.par['limits'])
                   if 'Use Limits' in eachItem[0]: subMenu.setChecked(self.par['use_limits'])
                   if 'Plot Lines' in eachItem[0]: subMenu.setChecked(self.par['draw_lines'])
                   if 'Auto Histogram' in eachItem[0]: subMenu.setChecked(self.par['auto_histo'])
                   if 'Filtered' in eachItem[0]: subMenu.setChecked(self.par['filtered'])
                      
                       
                   subMenu.triggered.connect(eachItem[2])
                   menu.addAction(subMenu)
                   
              
              
              self.mbar.addMenu(menu)          
          

 
          

     def set_file_info(self,filename):
          dir, fname = os.path.split(filename)
          name, ext = os.path.splitext(fname)
          self.dir = dir + '//'
          self.name = name
          self.ext = ext
          # that's it

     #----------------------------------------------------------------------
     # File menu routines   
     #----------------------------------------------------------------------
#     def OnScan(self): #scan for new files in data directory
#        path_to_watch = self.datadir
#        before = self.hwsproc
#        after = os.listdir (path_to_watch)
#        newf = False
#        self.hwsnew = []
#        for f in after:
#            if f.endswith(".hws") and not f in before:
#                print "Unprocessed file -- ", f
#                newf = True
#                self.hwsnew.append(f)
#        print "Directory scanned for unprocessed .hws files \n"    
#        if newf: 
#            self.repo=Repository(self, 'Select files to process', self.hwsnew)
#        else: print 'No unprocessed files found'           
#                #print self.choosefiles
                
    
     #process selected files           
    
#     def Process(self): 
#        start_time=time.time()
#        for f in self.hwstodo:
#                
#                #self.LoadParameters("C://Users//plasma//Desktop//par.data")
#                #self.OnPlot(None)
#                for i in range(0,1):
#                    self.OnFindPeaks()
#                    self.OnHistogram()
#                    os.mkdir(self.dir+self.name+'//')
#                    self.histoframe.figure.savefig(self.dir+self.name+'//Histogram_ch'+str(i)+'.png')
#                
#                self.RatesForAllChan(f)
#                self.figure.savefig(self.dir+self.name+'//RatePlot.png')
#                self.hwsproc.append(f)
#                self.OnTShistogram()
##                for i,h in enumerate(self.histos):
##                    h_file = self.dir + self.name +'//histo_{0}.data'.format(i)
##                    h.save(h_file)
##                print " all histograms saved"
#        self.proclistsave()
##        print self.hwstodo
#        self.hwstodo = []                
#        end_time=time.time()
#        work_time=end_time-start_time
#        print "Elapsed time ", work_time 
        
     def RatesForAllChan(self,ftw):
         self.OnClear()
         try:
             print((self.datadir + "/params.data"))
             self.LoadParameters(self.datadir + "params.data")
         except:
             pass
         for i in range(0,1):
             
             self.par['Detector channel'] = i
             self.OpenFile(self.datadir+"//"+ftw)
             self.OnTScalc()
             self.OnTSfindpeaks()
             self.OnTSplotrate()
     
             
#     def OnSetScanDir(self):
#         
#         print "Set data directory"
#         dir_dlg=QtWidgets.QFileDialog.getExistingDirectory(self, 'Select a directory')
#         if dir_dlg != '':
#             # User has selected something, get the path, set the window's title to the path
#             # store relevant file information
#             self.datadir = dir_dlg
#         else:
#             print "so, you changed your mind, I will do nothing"
#             return
#         print 'Data directory: ', self.datadir
#         #self.proclistload()
                    
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
             self.OpenFile(filename)
         else:
             print("so, you changed your mind, I will do nothing")
             filename = None
             return
         
         
     def OpenFile(self, fpname):
         self.set_file_info(fpname)
         self.OnLoad()
         self.select_data()
         print("Done")

     def OnLoad(self):
         
         chan_num = "%0d"%(int(self.par["Detector channel"]))
         self.stBar1.setText('Current file : %s / Channel : %s'%(self.name+self.ext, chan_num))
         # open file
         print(("Open file : ",self.dir + self.name + self.ext))
         self.f = h5py.File(self.dir + self.name + self.ext, 'r')
         # get the data
         # time information
         print("Get data")
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
         
         # WB temp solution 8/13/13
         # set the type as 16 bit signed, this is not good the data type should come from
         # the hdf file
         if convert_int:
              self.ydata = self.f[data_root + 'y-axis/data_vector/data'].value.astype('int16')
         else:
              self.ydata = self.f[data_root + 'y-axis/data_vector/data'].value

         # make the time axis
         print("Calculate data")
         self.tall = self.t0 + self.dt*np.arange(self.ydata.shape[0], dtype = float)
         # select the data to be worked on
         self.f.close()
         print("Select data")
     
     def OnReload(self):
         self.OnLoad()
         self.select_data()
         print("Done")
         
     def select_data(self):
          # select data range to work on
          # debuggin
          if self.ndata == 0:
              print('No data to select') 
              return
          if self.par["tmin"]<self.tall[0]: self.par["tmin"] = self.tall[0]
          if self.par["tmax"]>self.tall[-1]: self.par["tmax"] = self.tall[-1]  
          print("Get Window")
          sl = self.get_time_window( self.par["tmin"], self.par["tmax"])
          # the final value
          print("Recalculate")
          self.ndata = sl.stop - sl.start
          self.V = self.scale[0] + self.scale[1]*self.ydata[sl]
          self.t = self.tall[sl]
          
          print("Finished recalculation")
         
     def OnLoadParameters(self):
          # Create a file-open dialog in the current directory
          filetypes = '*.data'
          # use a regular file dialog
          file_dlg=QtWidgets.QFileDialog.getOpenFileName(self, 'Select a file', self.dir, filetypes )
          if file_dlg[0] != '':
               # User has selected something, get the path, set the window's title to the path
               filename = file_dlg[0]
               # store relevant file information
               print(("will open " , filename))
          else:
               print("so, you changed your mind, I will do nothing")
               filename = None
               return
          p = pfile(filename)
          for key in self.par:
               self.par[key] = p.get_value(key)
               print((key, " ", self.par[key]))
          # now set the check marks
          # get the menus
          
          menus = [i.menu() for i in self.mbar.actions()]
          items = None
          item_dict = {}
          # search for the options menu and get the contents
          for m in menus:
               if m.title() == "&Options":
                    items = m.actions()
                    for item in items:
                         if item.text() == "":
                              continue
                         else:
                              key = item.text().strip()[1:]
                              item_dict[key] = item
          # set the values for the options menu
          if item_dict != {}:
               item_dict["SetLimits"].setChecked(self.par["limits"]); print(("SetLimits = ", self.par["limits"]))
               item_dict["Measure"].setChecked(self.par["measure"]); print(("Measure = ", self.par["measure"]))
               item_dict["Use Limits"].setChecked(self.par["use_limits"]); print(("Use Limits = ",self.par["use_limits"]))
               item_dict["Auto Histogram Limits"].setChecked(self.par["auto_histo"]); print(("Auto Histogram Limits = ", self.par["auto_histo"]))
               item_dict["Filtered"].setChecked(self.par["filtered"]); print(("Filtered = ", self.par["filtered"]))
               item_dict["Plot Lines"].setChecked(self.par["draw_lines"]); print(("draw_lines = ", self.par["draw_lines"]))
          self.select_data()
          # all done
          
     def LoadParameters(self,fpname):
          
          p = pfile(fpname)
          for key in self.par:
               self.par[key] = p.get_value(key)
               print((key, " ", self.par[key]))
          # now set the check marks
          # get the menus
          
          menus = [i.menu() for i in self.mbar.actions()]
          items = None
          item_dict = {}
          # search for the options menu and get the contents
          for m in menus:
               if m.title() == "&Options":
                    items = m.actions()
                    for item in items:
                         if item.text() == "":
                              continue
                         else:
                              key = item.text().strip()[1:]
                              item_dict[key] = item
          # set the values for the options menu
          if item_dict != {}:
               item_dict["SetLimits"].setChecked(self.par["limits"]); print(("SetLimits = ", self.par["limits"]))
               item_dict["Measure"].setChecked(self.par["measure"]); print(("Measure = ", self.par["measure"]))
               item_dict["Use Limits"].setChecked(self.par["use_limits"]); print(("Use Limits = ",self.par["use_limits"]))
               item_dict["Auto Histogram Limits"].setChecked(self.par["auto_histo"]); print(("Auto Histogram Limits = ", self.par["auto_histo"]))
               item_dict["Filtered"].setChecked(self.par["filtered"]); print(("Filtered = ", self.par["filtered"]))
               item_dict["Plot Lines"].setChecked(self.par["draw_lines"]); print(("draw_lines = ", self.par["draw_lines"]))
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
               print(("will save to " , filename))
          else:
               print("so, you changed your mind, I will do nothing")
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
          print(("x positions    : %5.2e --> %5.2e" % (start_pos[0],  end_pos[0] )))
          print(("y positions    : %5.2e --> %5.2e" % (start_pos[1],  end_pos[1] )))
          print(("displacement :  delta-x = %5.2e   delta-y  = %5.2e " % tuple(dpos)))
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
                    print("cannot destroy histoframe !")
          if self.tsplotframe != None:
               try:
                    self.tsplotframe.destroy()
               except:
                    print("cannot destroy tsplotframe !")
          
          self.destroy()
          # all done
     #----------------------------------------------------------------------
     # Action menu routines   
     #----------------------------------------------------------------------
     def OnPlot(self):
          self.setCentralWidget(self.figure_canvas)
          V = None
          if self.par["filtered"]:
               V = self.Vinv
          else:
               V = self.V
          if (self.t is None) or (V is None):
               print("Nothing to plot !")
               return
          self.toolbar.t.append(self.t)
          self.toolbar.V.append(V)
          self.toolbar.ch_n.append(self.par['Detector channel'])
          self.toolbar.fn.append(self.name[-6:])
        
          self.toolbar.thinning()

     def OnClear(self):
          self.figure.clf()
          self.axes=self.figure.add_subplot(111)
          self.RS = RectangleSelector(self.axes, self.LineSelectCallback,
                                      drawtype='box', useblit=True,
                                      button=[1,3], # don't use middle button
                                      minspanx=5, minspany=5,
                                      spancoords='pixels')
          self.figure_canvas.draw()

     def OnFindPeaks(self):
          if (self.t is None) or (self.V is None):
               print("No data, nothing to find !")
               return
          if self.par["filtered"]:
               if self.Vinv is None:
                    print("No inverted data, nothing to find !")
                    return
               V = self.Vinv
          else:
               V = self.V
          print('----------> Start peak finding :')
          self.MAXTAB, self.MINTAB = find_peaks( V, self.par["Vstep"], xval = self.t, \
                                                 xmin = self.par["xmin"], xmax = self.par["xmax"], \
                                                 limits = self.par["use_limits"], \
                                                 get_window = get_window, \
                                                 power = 5)
          print(('----------> found ', len(self.MAXTAB[0]), ' peaks'))
          # get the values above the V threshold
          tp = np.array(self.MAXTAB[0])
          Vp = np.array(self.MAXTAB[1])
          iw = np.where( Vp > self.par["Vthreshold"] )[0]
          print(('----------> store ', len(iw), ' peaks above threshold'))
          self.t_peak = tp[iw]
          self.V_peak = Vp[iw]
          if self.par["auto_histo"] and not len(iw) == 0:
               self.par["Vhmin"] = self.V_peak.min()
               self.par["Vhmax"] = self.V_peak.max()
          
          self.axes.set_autoscaley_on(True)
          self.axes.plot(self.t_peak,self.V_peak, '.', color=colors[self.par['Detector channel']])
          self.figure_canvas.draw()

     def OnHistogram(self):
          if self.MAXTAB == None:
               print("No peaks yet !")
               return
          histo = LT.box.histo( self.V_peak, \
                                     range = (self.par["Vhmin"], self.par["Vhmax"]), \
                                     bins = int(self.par["VNbins"]))
          if histo == None:
               print("No histogram to plot !")
               return
          histo.title = 'V histogram'
          histo.xlabel = 'Volts'
          histo.ylabel = 'Counts'
          # add the histograms to the list of histos created
          self.histos.append( histo )
          # setup figure for hist plots
          latest_histo = len(self.histos) - 1
          self.show_histos(latest_histo)

     def OnDeleteHistogram(self):
          self.histos = []
          self.current_histo = 0
          if self.histoframe != None:
               # clear figure
               self.histoframe.figure.clf()
               self.histoframe.axes = self.histoframe.figure.add_subplot(111)
               self.histoframe.figure_canvas.draw()

     def OnSaveHistos(self):
          if self.histos == []:
               print("no histograms to save !")
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
               print(("will save to : "+dir+"/"+name+"_0.data ... " +name+"_{0}.data".format(N))) 
          else:
               print("so, you changed your mind, I will do nothing")
               filename = None
               return
          for i,h in enumerate(self.histos):
               h_file = dir + '/'+ name +'_{0}.data'.format(i)
               h.save(h_file)
          print(" all histograms saved")

     def show_histos(self, n):
          if self.histos == []:
               print("no histograms to show !")
               return
          # create it's own frame if it doesnot exist already
          self.current_histo = n
          if self.histoframe == None:
               self.histoframe = HistoPlotFrame(self, callback=self)
          else:
               # clear figure
               self.histoframe.figure.clf()
               self.histoframe.axes = self.histoframe.figure.add_subplot(111)
          self.histos[n].plot(axes = self.histoframe.axes)
          self.histoframe.figure_canvas.draw()
          

     def show_slice(self, n):
          if self.t_slice is None:
               print("no slices to show !")
               return
          # create it's own frame if it doesnot exist already
          self.current_slice = n
          if self.tsplotframe == None:
               self.tsplotframe = TSPlotFrame(self, callback = self)
          else:
               # clear figure
               self.tsplotframe.axes.cla()
          self.tsplot(n)
          self.tsplotframe.figure_canvas.draw()
          

     def OnTScalc(self):
          if (self.t is None) or (self.V is None):
               print("No data, nothing to find !")
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
          print('resizing time data array, be patient !')
          self.t_slice = np.resize(self.t[sl], new_shape )
          if self.par["filtered"]:
               if self.Vinv is None:
                    print("No inverted data, nothing to find !")
                    return
               V = self.Vinv
          else:
               V = self.V
          print('resizing Signal data array, be patient !')
          self.V_slice = np.resize(V, new_shape )
          # now the data are ready to be worked on
          self.ts_av = np.average(self.t_slice, axis = 1)
          print(("created ", len(self.ts_av), " slices"))

     def tsplot(self, n):
          if self.V_slice is None:
               print("No time slices!")
               return
          V = self.V_slice[n]
          t = self.t_slice[n]
          # show the data
          if self.par["draw_lines"]:
              self.tsplotframe.all_plot=self.tsplotframe.my_plot(t[:V.shape[0]], V) 
          else:
              self.tsplotframe.all_plot=self.tsplotframe.my_plot(t[:V.shape[0]], V, '.') 
          # show peak positions if there are any
          if self.ts_counts is not None:
               if self.ts_counts[n] != 0:
                   self.tsplotframe.axes.set_autoscaley_on(True) 
                   self.tsplotframe.axes.plot(self.ts_tp[n], self.ts_Vp[n], 'r.')
          self.tsplotframe.axes.set_xlabel('t')
          self.tsplotframe.axes.set_ylabel('V')
          s_title = 'Slice {0:d}, T_av = {1:10.3e}'.format(n, self.ts_av[n])
          self.tsplotframe.axes.set_title(s_title)
          self.tsplotframe.figure_canvas.draw()

 
     def OnTSfindpeaks(self):
          if self.V_slice is None:
               print("No time slices!")
               return
          ts_tp = []
          ts_Vp = []
          ts_counts = []
          for i,V in enumerate(self.V_slice):
               print(('----------> Start peak finding in slice :', i))
               MAXTAB, MINTAB = find_peaks( V, self.par["Vstep"], xval = self.t_slice[i], \
                                                      limits = False, \
                                                      get_window = get_window, \
                                                      power = 5)
               print(('----------> found ', len(MAXTAB[0]), ' peaks'))
               # get the values above the V threshold
               tp = np.array(MAXTAB[0])
               Vp = np.array(MAXTAB[1])
               iw = np.where( Vp > self.par["ts_Vthreshold"] )[0]
               print(('----------> store ', len(iw), ' peaks above threshold'))
               ts_counts.append(len(iw))
               ts_tp.append(tp[iw])
               ts_Vp.append(Vp[iw])
          self.ts_tp = ts_tp
          self.ts_Vp = ts_Vp
          self.ts_counts = np.array(ts_counts)


     def OnTSplotrate(self):
          if (self.ts_av is None) or (self.ts_counts is None):
               print("Nothing to plot !")
               return
          rate = self.ts_counts/self.par["ts_width"]
          rate_err = np.sqrt(self.ts_counts)/self.par["ts_width"]
          
          f1=pl.figure()
          pl.title('Rate')
          f1.axes[0].errorbar(self.ts_av, rate, yerr = rate_err, marker='o', ls = 'None')
          f1.axes[0].set_xlabel('t [s]')
          f1.axes[0].set_ylabel('Rate [Hz]')
          f1.draw()
         


     def OnTSsaverate(self):
          if (self.ts_av is None) or (self.ts_counts is None):
              print("Nothing to Save !")
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
               filename = file_dlg[0]
               # analyze the file name
               dir, fname = os.path.split(filename)
               name, ext = os.path.splitext(fname)
               print(("will save to : ", filename))
          else:
               print("so, you changed your mind, I will do nothing")
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
          if (self.ts_Vp is None):
               print("Nothing to histogram !")
               return
          # clear histos list
          self.histos = []
          for i,Vp in enumerate(self.ts_Vp):
               if Vp.any():
                   h = LT.box.histo( Vp, \
                                     range = (self.par["ts_Vhmin"], self.par["ts_Vhmax"]), \
                                     bins = int(self.par["ts_VNbins"]))
                   h.title = "V histogram # %d from: %6.3e to: %6.3e (s)"%(i, self.ts_tp[i][0], self.ts_tp[i][-1])
                   h.xlabel = 'Volts'
                   h.ylabel = 'Counts'
                   self.histos.append(h)
               else:
                   continue
          # hisograms created
          print(("Time slice histograms completed, %d histograms added" %len(self.histos)))

     def OnTSshowhistos(self):
          if self.histos == []:
               print("No TS histos")
               return
          # show ths histograms, start with the first one
          self.show_histos(0)

     def OnShow_next_histo(self):
          n = self.current_histo + 1
          if n >= len(self.histos):
               print("at the last histogram")
               n = len(self.histos)-1
          self.show_histos(n)
          
     def OnShow_prev_histo(self):
          n = self.current_histo - 1
          if n < 0:
               print("at the first histogram")
               n = 0
          self.show_histos(n)

     def OnTSshowslice(self):
          if self.ts_av is None:
               print("No slices")
               return
          # show the histograms, start with the first one
          self.show_slice(0)

     def OnShow_next_slice(self):
          n = self.current_slice + 1
          if n >= len(self.ts_av):
               print("at the last slice")
               n = len(self.ts_av)-1
          self.show_slice(n)
          
     def OnShow_prev_slice(self):
          n = self.current_slice - 1
          if n < 0:
               print("at the first slice")
               n = 0
          self.show_slice(n)
          
     #---------------------------------------------------------------------
     # FFT on time slices     
     #---------------------------------------------------------------------
     def OnTSFFTcalc(self):
          if self.ts_av is None:
               print("No slices")
               print("No data, nothing to transform !")
               return
          # calculate FFT and PS
          self.ts_FFT = []
          for i,V in enumerate(self.V_slice):
              t = self.t_slice[i]
 #             V = self.V_slice[i]
              print(("Calculate FFT for slice: ", i, " and ", len(t), " data points"))
              ffts = FFT.FFT(t)   
              ffts.transform(self.V)
              ffts.get_ps()
              ffts.logp = np.log10(ffts.p)
              # self.ts_FFT.append(ffts)
          print("all TS FFT completed") 

     def OnTSFFTplotps(self,event):
          # plot power spectrum
          if not self.fft.ok:
               print("No power spectrum to plot")
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
          data = {pkeys[0]:"%d"%(int(self.par["Detector channel"]))}
          pdlg = NumberDialog(data, self, title="Detector channel", \
                              labels = pkeys, \
                                  keys = pkeys, \
                                      about_txt = "Select the data channel")
          pdlg.exec_()
          # now set the new parameters
          self.par['Detector channel']=int(pdlg.data[pkeys[0]])
          pdlg.destroy()
          
          self.OnReload()
          self.OnPlot()

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
                  pkeys[2]:"%d"%(int(self.par["VNbins"]))}
          
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
                  pkeys[5]:"%d"%(int(self.par["ts_VNbins"])), \
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
          if (self.t is None) or (self.V is None):
               print("No data, nothing to transform !")
               return
          # calculate FFT and PS
          print(("Calculate FFT for : ", len(self.t), " data points"))
          self.fft = FFT.FFT(self.t)
          self.fft.transform(self.V)
          self.fft.get_ps()
          self.fft.logp = np.log10(self.fft.p)

     def OnFFTplotps(self):
          # plot power spectrum
          if not self.fft.ok:
               print("No power spectrum to plot")
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
               print("No power spectrum to plot")
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
               print("No power spectrum to plot")
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
               print("No FFT")
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
               print("No FFT")
               return
          self.fft.an_cut_freq(replace = True)
          self.fft.bn_cut_freq(replace = True)
          self.fft.store_fft_coeff(self.fft.an, self.fft.bn)
          self.fft.get_ps()
          self.fft.logp = np.log10(self.fft.p)

     def OnFFTrfilter(self):
          # Apply r-cut to FFT
          if not self.fft.ok:
               print("No FFT")
               return
          self.fft.r_cut_freq(replace = True)
          self.fft.store_fft_coeff(self.fft.an, self.fft.bn)
          self.fft.get_ps()
          self.fft.logp = np.log10(self.fft.p)
     
     def OnFFTinvert(self):
          print("Invert FFT")
          if not self.fft.ok:
               print("No FFT")
               return
          self.Vinv = self.fft.inv_transform()
          print("FFT inverted")
     #----------------------------------------------------------------------
     # options routines
     #----------------------------------------------------------------------

     def OnMeasure(self):
          # to measure create a Rectangle Selector
          self.RS.set_active(True)
          self.par["measure"] = True
          self.par["limits"] = False
          print(("limits = ", self.par["limits"], " measure = ", self.par["measure"]))

     def OnLimits(self):
          # to measure create a Selector
          self.RS.set_active(True)
          self.par["measure"] = False
          self.par["limits"] = True
          print(("limits = ", self.par["limits"], " measure = ", self.par["measure"]))

     def OnTogglelimits(self):
          self.par["use_limits"] = not self.par["use_limits"]
          print(("use limits = ", self.par["use_limits"]))

     def OnToggleHistolimits(self,event):
          self.par["auto_histo"] = not self.par["auto_histo"]
          print(("auto histo = ", self.par["auto_histo"]))

     def OnToggleLines(self, event):
          self.par["draw_lines"] = not self.par["draw_lines"]
          print(("draw Lines = ", self.par["draw_lines"]))
          try:
              self.toolbar.thinning()
          except:
              'Thinning did not work'
     def OnUsefiltered(self, event):
          self.par["filtered"] = not self.par["filtered"]
          print(("use filtered = ", self.par["filtered"]))
          
     def OnNothing(self):
          print("do nothing")

#----------------------------------------------------------------------
# main program
#----------------------------------------------------------------------
def RatesForAlS(app, chnum):
         fobj=PlotFrame(app)
         fobj.LoadParameters("C://Users//Alex//Desktop//par.data")
         fobj.OpenFile("C://Users//Alex//Desktop//DAQ_190813-112521.hws")
         fobj.OnTScalc()
         fobj.OnTSfindpeaks()
         if (fobj.ts_av is None) or (fobj.ts_counts is None):
               print("Nothing to plot !")
               return
         rate = fobj.ts_counts/fobj.par["ts_width"]
         rate_err = np.sqrt(fobj.ts_counts)/fobj.par["ts_width"]
         pl.errorbar(fobj.ts_av, rate, yerr = rate_err, marker='o', ls = 'None')
         fobj.destroy()
         
if __name__ == '__main__':
    #Create App
    app = QtCore.QCoreApplication.instance()
    
    if app is None:
        app = QtWidgets.QApplication(sys.argv)
    
    #Create PlotFrame in App
    frame = PlotFrame(app)
    #frame.LoadParameters('C:/Users/Alex/Desktop/test_data/params.data')
    #RatesForAlS(app, 0)
    sys.exit(app.exec_())
    #Test change
     
     
