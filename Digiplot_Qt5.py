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
import numpy as np
import LT.box as B
from LT.parameterfile import pfile
# FFT
import FFT
#Navigation toolbar and number dialog
from NavToolBar_CustomDialogs import NavigationToolbar, NumberDialog, TextDialog, LimitsListDialog
import ffind_peaks as FP
import itertools

import bkg_subtraction as BS

# matplotlibn colors
import matplotlib.colors as COL
# colors list
color_table = COL.CSS4_COLORS
color_names = list(color_table.keys())


# for argument passing 
import argparse as AG

convert_int = False  
# convert_int = True                       # added 3/3 from True to False based on LabView updates
colors = [color_table['red'],
          color_table['green'], 
          color_table['blue'], 
          color_table['magenta'], 
          color_table['darkcyan'], 
          color_table['orange'],
          color_table['firebrick'], 
          color_table['olive'], 
          color_table['mediumorchid']]

#----------------------------------------------------------------------
# useful functions
#----------------------------------------------------------------------
#%%

def is_even(x):
    return x%2 == 0

def moving_average(x, m):
    # calculate the moving average for M points, based on convolution
    # setup the

    if is_even(m) :
        m +=1
        print (f'you need an odd number of kernel points, new value m = {m}')
    kern = np.ones(m)/m
    xc = np.convolve(x, kern)[:x.shape[0]]
    return np.roll(xc, -m//2 + 1)
         
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
          pmin = np.zeros((len(yval[sl])//5, ), dtype='int32')
          pmax = np.zeros((len(yval[sl])//5, ), dtype='int32')
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


#%% standard peak shape     
def peak(x, a, b):
    # calculate the location of the maximum of the model peak shape
    if 2.*b/a-1.>0:
        x_0 = np.log(2.*b/a-1.)/(2.*b)
    else:
        x_0 = 0
    # value at the maximum location
    y_0 = np.exp(-a*(x_0))*(1.+np.tanh(b*(x_0)))
    # normalized shape centered at 0
    y = 1./y_0*np.exp(-a*(x+x_0))*(1.+np.tanh(b*(x+x_0)))
    # estimate full width
    large = y>=0.5
    try:
        il = np.where(large)[0]
        i1 = il.min()
        i2 = il.max()
        sig = x[i2] - x[i1]
    except:
        sig = 0.
    return sig, y 

def print_array(a, name):
    a_str = ''.join([f'{xx},' for xx in a])[:-1]
    print(f'{name} = [{a_str}]')
    


#----------------------------------------------------------------------
# load data functions, here more versions can be added
#----------------------------------------------------------------------

def load_gage_data(self, chan_num):
    print(("Open file : ",self.dir + self.name + self.ext))
    self.f = h5py.File(self.dir + self.name + self.ext, 'r')
    # load gage digitizer data
    data_root = 'wfm_group0/traces/trace' + chan_num + '/'
    xaxis_access = data_root + 'x-axis'
    self.t0 = self.f[xaxis_access].attrs['start'] + self.par['toffset']
    self.dt = self.f[xaxis_access].attrs['increment']
    # measured data
    # scale dataset
    yaxis_scale_access = data_root + 'y-axis/scale_coef'
    yaxis_data_access = data_root + 'y-axis/data_vector/data'
    # there are 2 values now
    #self.scale = self.f[yaxis_scale_access].value
    self.scale = self.f[yaxis_scale_access][0]
    # get the y dataset
    self.nall = self.f[yaxis_data_access].shape[0]         
    self.ndata = self.nall
    self.ti = self.t0
    self.tf = self.t0 + (self.ndata-1)*self.dt        
    
    # WB temp solution 8/13/13
    # set the type as 16 bit signed, this is not good the data type should come from
    # the hdf file
    if convert_int:
         self.ydata = self.f[data_root + 'y-axis/data_vector/data'].astype('int16')[:]
    else:
         self.ydata = self.f[data_root + 'y-axis/data_vector/data'][:]
         
    #added on 3/3/2021              
    for i, yd in enumerate(self.ydata[:10]):
        print(f'i = {i}, ydata = {yd}')
        
    #added on 12/21/2020:(oscillosope gives different peaks)    
    self.y_max = self.f[data_root + 'measurements/y_maximum'].attrs['value']
    
    # make the time axis
    print("Calculate data")
    self.tall = self.t0 + self.dt*np.arange(self.ydata.shape[0], dtype = float)
    # select the data to be worked on
    self.f.close()
    print("Select data")
    self.select_data()
    if self.par['auto_time_range']:
        self.par["tmin"] = self.ti
        self.par["tmax"] = self.tf       


def load_NI_data(self, chan_num):
    print(("Open file : ",self.dir + self.name + self.ext))
    self.f = h5py.File(self.dir + self.name + self.ext, 'r')
    # get the data
    # ee information
    print("Get data")
    data_root = 'wfm_group0/traces/trace' + chan_num + '/'
    try:
        self.t0 = self.f[data_root + 'x-axis'].attrs['start']
        self.dt = self.f[data_root + 'x-axis'].attrs['increment']
        # measured data
        # scale dataset
        self.scale = self.f[data_root + 'y-axis/scale_coef'][()]
        # get the y dataset
        self.nall = self.f[data_root + 'y-axis/data_vector/data'].shape[0]
    except:
        mb=QtWidgets.QMessageBox(self)
        mb.setText("Problems loading data " + data_root)
        mb.exec_()
        return
    self.ndata = self.nall
    self.ti = self.t0 - self.par['toffset']
    self.tf = self.t0 + (self.ndata-1)*self.dt - self.par['toffset']
    # set range
    #self.par["tmin"] = self.t0
    #self.par["tmax"] = self.t0 + (self.ndata-1)*self.dt
    # WB temp solution 8/13/13
    # set the type as 16 bit signed, this is not good the data type should come from
    # the hdf file
    if self.par["convert_int"]:
         self.ydata = self.f[data_root + 'y-axis/data_vector/data'][()].astype('int16')
    else:
         self.ydata = self.f[data_root + 'y-axis/data_vector/data'][()]
    
    # make the time axis
    print("Calculate data")
    self.tall = self.t0 + self.dt*np.arange(self.ydata.shape[0], dtype = float)
    # select the data to be worked on
    self.f.close()
    print("Select data")
    self.select_data()
    print("Done")
    if self.par['auto_time_range']:
        self.par["tmin"] = self.ti
        self.par["tmax"] = self.tf       
    

def load_npz_data(self, chan_num):
    # load npz data file
    print(("Open npz data file : ",self.dir + self.name + self.ext))
    self.f = np.load(self.dir + self.name + self.ext)
    d = self.f
    print("Get npz data")
    self.t0 = d['time'][0]
    self.dt = np.diff(d['time'])[0]
    self.scale = [0.,1.]
    self.nall = len(d['time'])
    self.ndata = self.nall
    self.ti = self.t0
    self.tf = self.t0 + (self.ndata-1)*self.dt
    self.tall = d['time']
    self.ydata = d['signal']
    print("Select data")
    self.select_data()
    print("Done")
    if self.par['auto_time_range']:
        self.par["tmin"] = self.ti
        self.par["tmax"] = self.tf       
   

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
          self.toolbar.mker = '.'
          self.toolbar.colors = colors
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
        



#----------------------------------------------------------------------
# master frame, contains the plot image
#----------------------------------------------------------------------
class PlotFrame(QtWidgets.QMainWindow):
     def __init__(self, parent, key = 'PF'):
          super(PlotFrame, self).__init__()
          #QtWidgets.QMainWindow.__init__(self)
         
          self.setWindowTitle(f'DigiPlot_{key}')
          
          self.key = key
          # Parameters
          self.par = {}
          self.par["measure"] = True
          self.par["limits"] = False
          self.par["use_limits"] = False
          self.par["auto_histo"] = True
          self.par["auto_time_range"] = True
          self.par["corrected"] = False
          self.par["draw_lines"] = False
          self.par["draw_points"] = True
          self.par["hide_peaks"] = False
          self.par["Detector channel"] = 0
          self.par["Bkg Detector channel"] = 4
          self.par["Bkg Iter"] = 5
          self.par["Bkg Delta T"] = 200e-6
          self.par["toffset"] = 0.
          self.par["tmin"] = -1.e30
          self.par["tmax"] = 1.e30
          self.par["xmin"] = -1.e30
          self.par["xmax"] = 1.e30
          self.par["scale_x"] = 1.0
          self.par["scale_y"] = 1.0
          self.par["Vstep"] = 0.3
          self.par["Vthreshold"] = 0.3
          self.par["Vhmin"] = 0.0
          self.par["Vhmax"] = 1.0
          self.par["VNbins"] = 100
          self.par["Navg"] = 21
          self.par["alpha"] = 0.75
          # time slices
          self.par["ts_width"]=1e-3
          self.par["ts_start"] = 0.0
          self.par["ts_stop"] = 1.0
          self.par["ts_Vthreshold"] = 0.0   # for time slice histograms
          self.par["ts_Vwidth"] = 0.2   # window width for rate calculation
          self.par["ts_Vhmin"] = 0.0
          self.par["ts_Vhmax"] = 1.0
          self.par["ts_VNbins"] = 100
          self.par["use_npz_file"] = False # use npz data
          self.par["use_NI_file"] = False # use NI digitizer data
          self.par["use_GaGe_file"] = True # use Gage digitizer data
          self.par["convert_int"] = False # integer conversion
          self.par["plot_histo_points"] = False # plot histogram data as points
          
          self.datadir = '../Raw_Data/'
          self.dt = -1.
          # inital settings
          self.ndata = 0
          self.t = None
          self.V = None
          self.limits = []
          self.Vcorr = None
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
          self.setCentralWidget(self.figure_canvas)
          self.all_plot=[]
         
          
          #self.setCentralWidget(self.figure_canvas)
          
          
          
          # add ToolBar
          self.toolbar = NavigationToolbar(self.figure_canvas, self)
          self.toolbar.colors=colors
          self.toolbar.markers=itertools.cycle(['.','*'])
          self.toolbar.mker=next(self.toolbar.markers)
          self.toolbar.peak_marker = 'x'
          self.addToolBar(self.toolbar)
          
          #self.RS = RectangleSelector(self.axes, self.LineSelectCallback)
          
          self.RS = RectangleSelector(self.axes, self.LineSelectCallback,
                                      useblit=True,
                                      button=[1,3], # don't use middle button
                                      minspanx=5, minspany=5,
                                      spancoords='pixels')
          
          self.setPalette(self.palette)
          self.setGeometry(150,150,800,600)
          self.createMenubar()
          self.show()

                 
         
     def menuData(self): # data for the menu
          return(
               ("&File",                                        # menu label
                ("&Open", "Select data files", self.OnSelectFile),  # menu items consisting of: label, text, handler
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
                ("S&ubtract Background", "Subtract background", self.OnSubtractBkg),
                ("&Save Limits", "save curent x-axis limits", self.OnSaveLimits),
                ("&Choose  Limits", "use saved x-axis limits", self.OnChooseLimits),
                ("&Clear  Limits", "clear x-axis limits", self.OnClearLimits),
                ("&Clear  Peaks", "clear peka positions", self.OnClearPeaks),
                (None, None, None),  # creates a separator bar in the menu
                ("&FindPeaks", "find peaks in the plotting data", self.OnFindPeaks),
                ("&Moving Average", "smooth data via moving avg.", self.OnMovingAvg),
                (None, None, None),
                ("&Histogram", "histogram peak data", self.OnHistogram),
                ("&Fit Histogram", "Fit histogram peak area", self.OnFitHistogram),
                ("&Sum Histogram", "Sum histogram peak area", self.OnSumHistogram),
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
                ("Time &Slicing", "Set time slice parameters", self.OnSetTimeSlice),
                (None, None, None),
                ("Moving &Average", "Set number points", self.OnSetMovingAvg),
                (None, None, None),
                ("Time offset", "Set time offset", self.OnSetTimeOffset)),               
               ("&Bkg Subtraction Parameters",
#                ("&Data directory", "Set data directory", self.OnSetScanDir), 
#                (None, None, None),
                ("&Background Channel", "Set backgound channel", self.OnSelectBkgChannel), 
                (None, None, None),
                ("&Time Window", "Set time window width", self.OnSelectBkgTimeWindow), 
                (None, None, None),
                ("Moving Average &Iterations", "Set number interations for moving averagew", self.OnSelectBkgMovingAvgIter)),
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
                ("&Create 2D Histogram", "Create a 2D Histogram",self.OnTS2Dhistogram ),
                ("&Create Histograms", "Create Histograms",self.OnTShistogram ),
                ("&Delete Histograms", "delete all histograms", self.OnDeleteHistogram),
                ("&Show Histograms", "Plot Histograms", self.OnTSshowhistos)),
               ("&Options",
                ("&Use npz file", "Use npz data", self.OnUsenpz_file, 'RADIO', 'filetypes'),
                ("&Use NI file", "Use NI data", self.OnUseNI_file, 'RADIO', 'filetypes'),
                ("&Use GaGe file", "Use Gage data", self.OnUseGaGe_file, 'RADIO', 'filetypes'), 
                (None, None, None),
                ("&Convert Integer", "Convert raw data to int16", self.OnToggleConvertInteger, 'CHECK'),
                ("&Auto Time Range", "Set time range automatically", self.OnToggleAutoTimeRange, 'CHECK'),
                (None, None, None),
                ("&ChooseLimits", "Set limits for data processing", self.OnLimits, 'RADIO', 'rect'),
                ("&Measure", "Measure differences in t and V", self.OnMeasure, 'RADIO', 'rect'),
                ("&Use Limits ", "Use Limits", self.OnTogglelimits, 'CHECK'),
                ("&Plot Lines ", "Draw lines", self.OnToggleLines, 'CHECK'),
                ("&Plot Points ", "Draw Points", self.OnTogglePoints, 'CHECK'),
                ("&Hide Peaks ", "Do not show peaks", self.OnToggleHidePeaks, 'CHECK'),
                (None, None, None),
                ("Plot &Histogram Points ", "Plot points instead of bars", self.OnToggleHistoPoints, 'CHECK'),
                ("&Auto Histogram Limits ", "Do not calculate histogram limits automatically", self.OnToggleHistolimits, 'CHECK'),
                (None, None, None),
                ("&Corrected", "Use corrected data", self.OnUseCorrected, 'CHECK')),
               )
               
     def createMenubar(self):
        groups = {}
        for eachMenuData in self.menuData():
              
             menuLabel = eachMenuData[0]
             menuItems = eachMenuData[1:]
              
             menu = QtWidgets.QMenu(menuLabel, self)
             for eachItem in menuItems:
                  if not eachItem[0]:
                       menu.addSeparator()
                       continue
                  subMenu=QtWidgets.QAction(eachItem[0],self) 
                 
                  if len(eachItem) == 4:       
                       if eachItem[3]=='CHECK':
                           subMenu.setCheckable(True)
        
                  if len(eachItem) == 5:   # item has a group assignment
                       groupName = eachItem[-1]
                       if groupName in list(groups.keys()):
                           pass
                       else:
                           groups[groupName] = QtWidgets.QActionGroup(self)
                       if eachItem[3]=='RADIO':
                           subMenu.setCheckable(True)
                           subMenu.setActionGroup(groups[groupName])
                                                    
                  # set initial checks
                  if 'Convert Integer' in eachItem[0]: subMenu.setChecked(self.par['convert_int'])         
                  if 'Measure' in eachItem[0]: subMenu.setChecked(self.par['measure'])           
                  if 'ChooseLimits' in eachItem[0]: subMenu.setChecked(self.par['limits'])
                  if 'Use Limits' in eachItem[0]: subMenu.setChecked(self.par['use_limits'])
                  if 'Plot Lines' in eachItem[0]: subMenu.setChecked(self.par['draw_lines'])
                  if 'Plot Points' in eachItem[0]: subMenu.setChecked(self.par['draw_points'])
                  if 'Hide Peaks' in eachItem[0]: subMenu.setChecked(self.par['hide_peaks'])
                  if 'Auto Histogram' in eachItem[0]: subMenu.setChecked(self.par['auto_histo'])
                  if 'Histogram Points' in eachItem[0]: subMenu.setChecked(self.par['plot_histo_points'])
                  if 'Corrected' in eachItem[0]: subMenu.setChecked(self.par['corrected'])
                  if 'Use npz file' in eachItem[0]: subMenu.setChecked(self.par['use_npz_file'])
                  if 'Use NI file' in eachItem[0]: subMenu.setChecked(self.par['use_NI_file'])
                  if 'Use GaGe file' in eachItem[0]: subMenu.setChecked(self.par['use_GaGe_file'])
                  if 'Auto Time Range' in eachItem[0]: subMenu.setChecked(self.par['auto_time_range'])
                     
                      
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
     
             
                    
     def OnSelectFile(self):
        # Create a file-open dialog in the current directory
        if self.par['use_npz_file']:
            filetypes = '*.npz'
        else:
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
            self.load_data()
            self.select_data()
            print("Done")
        else:
            print("so, you changed your mind, I will do nothing")
            filename = None
            return
        

     def load_data(self):
         # reload data file
         # get channel number  
         if (self.name is None) or (self.name == ''):
             print('No files loaded yet , nothing to bedone !')
             mb=QtWidgets.QMessageBox(self)
             mb.setText('No files loaded yet , nothing to bedone !')
             mb.exec_()

             return
         chan_num = "%0d"%(int(self.par["Detector channel"]))
         self.stBar1.setText('Current file : %s / Channel : %s'%(self.name+self.ext, chan_num))
         # get the data according to the data file type
         print("Get data")
         # clear peaks
         self.OnClearPeaks()
         try:
             if self.par["use_GaGe_file"]:   
                 load_gage_data(self, chan_num)
             elif self.par["use_NI_file"]:
                 load_NI_data(self, chan_num)
             elif self.par["use_npz_file"]:
                 load_npz_data(self, chan_num)
         except Exception as err:
             print(f'Cannot load data from {self.f}: {err}')
             return
         
         
     def OnReload(self):
         self.load_data()
         # self.select_data()
         print("Done")
         
         
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
               try:
                   self.par[key] = p.get_value(key)
               except Exception as err:
                   print(f'Problem reading {key} from parameter file: {err}, skipping!')
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
                              key = item.text().split('&')[-1].strip()   # use only last part
                              item_dict[key] = item
          # set the values for the options menu
          if item_dict != {}:
               item_dict["ChooseLimits"].setChecked(self.par["limits"]); print("ChooseLimits = ", self.par["limits"])
               item_dict["Convert Integer"].setChecked(self.par["convert_int"]); print(("Convert Integer = ", self.par["convert_int"]))
               item_dict["Measure"].setChecked(self.par["measure"]); print("Measure = ", self.par["measure"])
               item_dict["Use Limits"].setChecked(self.par["use_limits"]); print("Use Limits = ",self.par["use_limits"])
               item_dict["Auto Histogram Limits"].setChecked(self.par["auto_histo"]); print("Auto Histogram Limits = ", self.par["auto_histo"])
               item_dict['Histogram Points'].setChecked(self.par['plot_histo_points']);print(("Plot Histogram Points = ", self.par["plot_histo_points"]))
               item_dict["Corrected"].setChecked(self.par["corrected"]); print("corrected = ", self.par["corrected"])
               item_dict["Plot Lines"].setChecked(self.par["draw_lines"]); print("draw_lines = ", self.par["draw_lines"])
               item_dict["Plot Points"].setChecked(self.par["draw_points"]); print("draw_points = ", self.par["draw_points"])
               item_dict["Hide Peaks"].setChecked(self.par["hide_peaks"]); print("hide_peaks = ", self.par["hide_peaks"])
               item_dict['Use npz file'].setChecked(self.par['use_npz_file']); print("use_npz_file = ", self.par['use_npz_file'])
               item_dict['Use NI file'].setChecked(self.par['use_NI_file']); print("use_NI_file = ", self.par['use_NI_file'])
               item_dict['Use GaGe file'].setChecked(self.par['use_GaGe_file']); print("use_GaGe_file = ", self.par['use_GaGe_file'])
               item_dict['Auto Time Range'].setChecked(self.par['auto_time_range']); print("auto_time_range = ", self.par['auto_time_range'])

          self.select_data()
          # all done
          
     def LoadParameters(self,fpname):
          
          p = pfile(fpname)
          for key in self.par:
               self.par[key] = p.get_value(key)
               print((key, " ", self.par[key]))
          # now set the check marks
          # get the menus
          self.update_options_menu()        

 
     def update_options_menu(self):
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
                              key = item.text().split('&')[-1].strip()   # use only last part
                              item_dict[key] = item
          # set the values for the options menu
          print(f'Options item dict = {item_dict.keys()}')
          if item_dict != {}:
               item_dict["ChooseLimits"].setChecked(self.par["limits"]); print(("ChooseLimits = ", self.par["limits"]))
               item_dict["Convert Integer"].setChecked(self.par["convert_int"]); print(("Convert Integer = ", self.par["convert_int"]))
               item_dict["Measure"].setChecked(self.par["measure"]); print(("Measure = ", self.par["measure"]))
               item_dict["Use Limits"].setChecked(self.par["use_limits"]); print(("Use Limits = ",self.par["use_limits"]))
               item_dict["Auto Histogram Limits"].setChecked(self.par["auto_histo"]); print(("Auto Histogram Limits = ", self.par["auto_histo"]))
               item_dict['Histogram Points'].setChecked(self.par['plot_histo_points']);print(("Plot Histogram Points = ", self.par["plot_histo_points"]))
               item_dict["Corrected"].setChecked(self.par["corrected"]); print(("corrected = ", self.par["corrected"]))
               item_dict["Plot Lines"].setChecked(self.par["draw_lines"]); print(("draw_lines = ", self.par["draw_lines"]))
               item_dict["Plot Points"].setChecked(self.par["draw_points"]); print(("draw_points = ", self.par["draw_points"]))
               item_dict["Hide Peaks"].setChecked(self.par["hide_peaks"]); print(("hide_peaks = ", self.par["hide_peaks"]))
               item_dict['Use npz file'].setChecked(self.par['use_npz_file']); print(("use_npz_file = ", self.par['use_npz_file']))
               item_dict['Use NI file'].setChecked(self.par['use_NI_file']); print(("use_NI_file = ", self.par['use_NI_file']))
               item_dict['Use GaGe file'].setChecked(self.par['use_GaGe_file']); print(("use_GaGe_file = ", self.par['use_GaGe_file']))
               item_dict['Auto Time Range'].setChecked(self.par['auto_time_range']); print(("auto_time_range = ", self.par['auto_time_range']))
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
          

     def check_time_limits(self, t):
         # make sure limits are within time range, if not adjust them
         self.par["tmin"] = max(t[0], self.par["tmin"])
         self.par["tmax"] = min(t[-1], self.par["tmax"])


     def select_data(self):
          # select data range to work on
          if not hasattr(self, "tall"):
              print('No data !')
              return
          # debugging
          if self.ndata == 0:
               return
          self.check_time_limits(self.tall)   
          print("Get Window")
          sl = get_window( self.par["tmin"], self.tall,  self.par["tmax"])
          print(f'DigiPlot time slice = {sl}')
          # the final value
          print("Recalculate")
          self.ndata = sl.stop - sl.start
          if self.par['use_NI_file']:
              self.V = self.scale[0] + self.scale[1]*self.ydata[sl]
              self.t = self.tall[sl]
          else:
              self.V = self.ydata[sl]
              self.t = self.tall[sl]
          print("Finished recalculation")

          
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
               delta_t = dpos[0]
               f_measure = 1./delta_t
               mb=QtWidgets.QMessageBox(self)
               mb.setText(f"Time interval = {delta_t:.3e}s / f = {f_measure:.3e} Hz ")
               mb.exec_()

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
          
          print("Closing All")
          self.destroy()
          # self.parent.quit()
          frames.pop(self.key)
          if frames == {}:
              QtCore.QCoreApplication.instance().exit()
          # all done
     def closeEvent(self, event):
         self.OnCloseWindow()

     #----------------------------------------------------------------------
     # Action menu routines   
     #----------------------------------------------------------------------
     def OnPlot(self):
          #if (self.par['Detector channel'] ==  self.toolbar.ch_n) and (self.name[-6:] == self.toolbar.fn):
          #    return
          V = None
          if self.par["corrected"]:
               V = self.Vcorr
          else:
               V = self.V
          if (self.t is None) or (V is None):
               print("Nothing to plot !")
               return
          # store data in toolbar, necessary for thinning
          self.axes.set_xlim(self.t.min(), self.t.max())
          self.toolbar.set_data(self.t, V)
          self.toolbar.xlabel = 't'
          self.toolbar.ylabel = 'V'
          self.toolbar.ch_n = self.par['Detector channel']
          self.toolbar.alpha = self.par['alpha']
          print(f'---> set toolbar alpha to {self.toolbar.alpha}')
          self.toolbar.fn = self.name[-6:]
        
          self.toolbar.thinning()
    
     def OnSubtractBkg(self):
        print('Will perfrom Bkg Subtraction using the follwing parameters:')
        print(f'Bkg data channel =  {self.par["Bkg Detector channel"]}')
        print(f'Time window = {self.par["Bkg Delta T"]}')
        print(f'Number of moving avarage iterations  = {self.par["Bkg Iter"]}')
        print(f'Moving agverage width = {self.par["Bkg Delta T"]}')
        # initialize bkg subtraction
        file_name = self.name + self.ext
        bd = BS.bkg(file_name, self.dir,\
                    data_channel = self.par['Detector channel'],\
                    bkg_channel = self.par['Bkg Detector channel'],\
                    t_min = self.par["tmin"]*BS.us, \
                    t_max = self.par["tmax"]*BS.us, \
                    t_offset = self.par["toffset"]*BS.us, \
                    moving_average = self.par["Navg"],\
                    niter = self.par["Bkg Iter"],\
                    time_window = self.par["Bkg Delta T"]*BS.us)
        bd.load_data()
        bd.moving_average()
        # perform the correction and save the data
        tl, self.Vcorr = bd.correct()

     def OnMovingAvgIter(self):
         self.V = moving_average(self.V, self.par["Navg"])
         print('Moving avg. calculated')
        
    
     def OnSaveLimits(self):
          pkeys=['Enter Comment']
          data = {pkeys[0]:""}
          pdlg = TextDialog(data, self, title="Comment for plot range", \
                              labels = pkeys, \
                                  keys = pkeys, \
                                      about_txt = "Save current plot range")
          pdlg.exec_()
          # now set the new parameters
          comment = pdlg.data[pkeys[0]]
          pdlg.destroy()
          
          self.limits.append({'value':self.axes.get_xlim(), 'comment':comment})
          
     def OnChooseLimits(self):
          if self.limits == []:
              print('Limits not yet set !')
              return
          dlg = LimitsListDialog(self.limits)
          
          dlg.exec_()
          self.axes.set_xlim(self.limits[dlg.selected_index]['value'])
          dlg.destroy()
          # redraw
          self.toolbar.thinning()
          
          
     def OnClearLimits(self):
          self.limits = []

     def OnClear(self):
          self.figure.clf()
          self.axes=self.figure.add_subplot(111)
          self.RS = RectangleSelector(self.axes, self.LineSelectCallback,
                                      useblit=True,
                                      button=[1,3], # don't use middle button
                                      minspanx=5, minspany=5,
                                      spancoords='pixels')
          self.figure_canvas.draw()
          self.toolbar.ch_n=self.par['Detector channel']
          self.toolbar.set_data(np.array([]), np.array([]))
          self.toolbar.fn=''
          self.toolbar.markers=itertools.cycle(['.','*'])
          self.toolbar.mker=next(self.toolbar.markers)
          self.toolbar.peak_marker = 'x'
          
     def OnFindPeaks(self):
          if (self.t is None) or (self.V is None):
               print("No data, nothing to find !")
               return
          if self.par["corrected"]:
               if self.Vcorr is None:
                    print("No corrected data, nothing to find !")
                    return
               V = self.Vcorr
          else:
               V = self.V
          print('----------> Start peak finding :')
          self.MAXTAB, self.MINTAB = find_peaks( V, self.par["Vstep"], xval = self.t, \
                                                 xmin = self.par["xmin"], xmax = self.par["xmax"], \
                                                 limits = self.par["use_limits"], \
                                                 get_window = get_window, \
                                                 power = 5)
          print('----------> found ', len(self.MAXTAB[0]), ' peaks')
          # get the values above the V threshold
          tp = np.array(self.MAXTAB[0])
          Vp = np.array(self.MAXTAB[1])
          iw = np.where( Vp > self.par["Vthreshold"] )[0]
          print('----------> store ', len(iw), ' peaks above threshold')
          self.t_peak = tp[iw]
          self.V_peak = Vp[iw]
          if self.par["auto_histo"] and not len(iw) == 0:
               self.par["Vhmin"] = self.V_peak.min()
               self.par["Vhmax"] = self.V_peak.max()
          
          self.axes.set_autoscaley_on(True)
          self.toolbar.set_peakdata(self.t_peak, self.V_peak)
          self.toolbar.plot_peaks = not self.par["hide_peaks"] 
          if len(self.toolbar.x) == 0:
              # not data for thinning
              self.axes.plot(self.t_peak,self.V_peak, self.toolbar.peak_marker , color='m')
          else:
              self.toolbar.thinning()
          
          self.figure_canvas.draw()
          
     def OnClearPeaks(self):
         self.t_peak = []
         self.V_peak = []
         self.toolbar.set_peakdata(self.t_peak, self.V_peak)
     
     def OnMovingAvg(self):
         self.V = moving_average(self.V, self.par["Navg"])
         print('Moving avg. calculated')
                
         
     def OnHistogram(self):
          if self.MAXTAB == None:
               print("No peaks yet !")
               return
          histo = B.histo( self.V_peak, \
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
          
     def OnFitHistogram(self):
         # fit current histo
         # pl.sca(self.histoframe.axes)
         if self.histos == []:
             print('No histrograms to fit!')
             return         
         h = self.histos[self.current_histo]
         # get current plot limits (useful for setting fit limits)
         cmin, cmax = self.histoframe.axes.get_xlim()
         # fit the hitogram within these limits
         h.fit(cmin, cmax, plot_fit = False)
         h.plot_fit(axes =  self.histoframe.axes)
         # calculate peak parameters for title
         fwhm = h.sigma.value*2.355
         pos = h.mean.value
         res = fwhm/pos*100.  # resolution in %
         print ('FWHM = ', h.sigma.value*2.355)
         self.histoframe.figure_canvas.draw()
         # new histo with finer binning
         try:
             hn = B.histo(self.V_peak, range=(cmin*0.9, cmax*1.1), bins= int(self.par["VNbins"]))
             self.histos.append( hn )
             hn.title = 'FWHM = {:.3e}, Position = {:.3e}, Resolution = {:.2e}%'.format(fwhm, pos, res)
             hn.fit(cmin, cmax, plot_fit = False)
             hn.plot()
             hn.plot_fit()
             latest_histo = len(self.histos) - 1
             self.show_histos(latest_histo)
         except Exception as err:
             print(f'cannot create detail histogram: {err}')
             print('If the error message is : "PlotFrame" object has no attribute "V_peak", run "Find Peaks" for all data first')

     def OnSumHistogram(self):    
         # sum up current histo
         # pl.sca(self.histoframe.axes)
         if self.histos == []:
             print('No histrograms to sum!')
             return
         h = self.histos[self.current_histo]
         # get current plot limits (useful for setting fit limits)
         cmin, cmax = self.histoframe.axes.get_xlim()
         # fit the hitogram within these limits
         S = h.sum(cmin, cmax)
         # present result in a message boxk
         mb=QtWidgets.QMessageBox(self)
         mb.setText(f"Histogram sum between {cmin:.2e} and {cmax:.2e} = {S[0]:.3e} +/- {S[1]:.1e} ")
         mb.exec_()
         
     def OnDeleteHistogram(self):
          self.histos = []
          self.current_histo = 0
          if self.histoframe != None:
               # clear figure
               self.histoframe.figure.clf()
               self.histoframe.axes = self.histoframe.figure.add_subplot(111)
               self.histoframe.figure_canvas.draw()
               print ('All histograms have been deleted')

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
          is_2d = self.histos[n].__class__ == B.histo2d
          if is_2d:
              # for 2d histos
              self.histos[n].plot(axes = self.histoframe.axes)
          else:
              if self.par['plot_histo_points']:
                  self.histos[n].plot_exp(axes = self.histoframe.axes, ignore_zeros = True)
                  self.histos[n].plot_fit(axes = self.histoframe.axes, ignore_zeros = True)
              else:
                  self.histos[n].plot(axes = self.histoframe.axes)
                  self.histos[n].plot_fit(axes = self.histoframe.axes)
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
          if self.par["corrected"]:
               if self.Vcorr is None:
                    print("No corrected data, nothing to find !")
                    return
               V = self.Vcorr
          else:
               V = self.V
          print('resizing Signal data array, be patient !')
          self.V_slice = np.resize(V, new_shape )
          # now the data are ready to be worked on
          self.ts_av = np.average(self.t_slice, axis = 1) + self.par['toffset']
          print(("created ", len(self.ts_av), " slices"))

     def tsplot(self, n):
          if self.V_slice is None:
               print("No time slices!")
               return
          V = self.V_slice[n]
          t = self.t_slice[n]
          self.tsplotframe.toolbar.set_data(t,V)
          self.tsplotframe.toolbar.xlabel = 't'
          self.tsplotframe.toolbar.ylabel = 'V'
          self.tsplotframe.toolbar.ch_n = self.par['Detector channel']
          self.tsplotframe.toolbar.fn = self.name[-6:]
          s_title = 'Slice {0:d}, T_av = {1:10.3e}'.format(n, self.ts_av[n])
          self.tsplotframe.toolbar.thinning(title = s_title)

 
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
               iw = np.where( (Vp > self.par["ts_Vthreshold"]) & 
                              (Vp < self.par["ts_Vthreshold"] +  self.par["ts_Vwidth"]) )[0]
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
          f1.axes[0].ticklabel_format(axis='y', style='sci', scilimits=(0,0))
          f1.axes[0].errorbar(self.ts_av, rate, linestyle='-', marker='.', capsize = 0., color=colors[self.par['Detector channel']])
          f1.axes[0].fill_between(self.ts_av, rate-rate_err, rate+rate_err, alpha=0.3, color=colors[self.par['Detector channel']])
          f1.axes[0].set_xlabel('t [s]')
          f1.axes[0].set_ylabel('Rate [Hz]')
          f1.show()
          
          # pl.locator_params(axis='x', nbins=30)


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
                   h = B.histo( Vp, \
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

     def OnTS2Dhistogram(self):
          if (not hasattr(self, 'V_peak')) or (not hasattr(self, 't_peak')):
              print("No peak data run find peaks first")
              return
          # time slice the data
          dt = self.par['ts_width'] # step width
          # 2d histogram setup
          # y - aaxis
          hy_min = self.par['ts_Vhmin']
          hy_max = self.par['ts_Vhmax']
          hy_bins = self.par['ts_VNbins']        
          # x - axis
          tmin = self.par['ts_start'] 
          tmax = self.par['ts_stop']
          hx_bins = int((tmax - tmin)/dt) + 1           
          h_title = f'channel: {self.par["Detector channel"]}'             
          h2 = B.histo2d(self.t_peak, self.V_peak, range = [[tmin,tmax],[hy_min,hy_max]], bins = [hx_bins, hy_bins],
                              title = h_title, xlabel = r't[$\mu$ s]', ylabel = 'Pulse Height [V]')
          self.histos.append(h2)
     
     def OnShow2dHisto(self):
         pass
          
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
          elif self.par["draw_lines"] and self.par["draw_points"] :
               self.axes.plot(self.fft.f,self.fft.logp, '.-')
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

     def OnSelectBkgChannel(self):
          # show and change parameters
          # parameter keys
          pkeys=['Bkg channel']
          data = {pkeys[0]:"%d"%(int(self.par["Bkg Detector channel"]))}
          pdlg = NumberDialog(data, self, title="Backgound Detector channel", \
                              labels = pkeys, \
                                  keys = pkeys, \
                                      about_txt = "Select the Bkg channel")
          pdlg.exec_()
          # now set the new parameters
          self.par['Bkg Detector channel']=int(pdlg.data[pkeys[0]])
          pdlg.destroy()
          
          self.OnReload()

     def OnSelectBkgTimeWindow(self):
          # show and change the parameters
          # parameter keys
          pkeys=['Backgound time window']
          # current data to be shown in the dialog
          data = {pkeys[0]:"%5.2e"%(self.par["Bkg Delta T"]), \
                  }
          pdlg = NumberDialog(data, self, title="Background Time Window", \
                              labels = pkeys, \
                              keys = pkeys, \
                              about_txt = "Select the time window")
          pdlg.exec_()
          # now set the new parameters
          self.par["Bkg Delta T"] = float(pdlg.data[pkeys[0]])
          
     
     def OnSelectBkgMovingAvgIter(self):
          # show and change parameters
          # parameter keys
          pkeys=['Moving Average Iterations']
          data = {pkeys[0]:"%d"%(int(self.par["Bkg Iter"]))}
          pdlg = NumberDialog(data, self, title="Moving average Iterations", \
                              labels = pkeys, \
                                  keys = pkeys, \
                                      about_txt = "Select number of iterations")
          pdlg.exec_()
          # now set the new parameters
          self.par["Bkg Iter"]=int(pdlg.data[pkeys[0]])
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
          tmin = float(pdlg.data[pkeys[1]])
          tmax = float(pdlg.data[pkeys[2]])
          pdlg.destroy()
          if (tmax == tmin) :
              print("t_min < t_max is necessary !")
              return
          elif(tmax < tmin):
              print("t_min < t_max is necessary, will switch the numbers!")
              self.par["tmin"] = tmax
              self.par["tmax"] = tmin
          else:
              self.par["tmin"] = tmin
              self.par["tmax"] = tmax          
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
          pkeys=['t start','t stop','t width', 'V min', 'V max', 'N bins', 'V threshold', 'V width']
          # current data to be shown in the dialog
          data = {pkeys[0]:"%6.2e"%(self.par["ts_start"]), \
                  pkeys[1]:"%6.2e"%(self.par["ts_stop"]), \
                  pkeys[2]:"%6.2e"%(self.par["ts_width"]), \
                  pkeys[3]:"%5.3f"%(self.par["ts_Vhmin"]), \
                  pkeys[4]:"%5.3f"%(self.par["ts_Vhmax"]), \
                  pkeys[5]:"%d"%(int(self.par["ts_VNbins"])), \
                  pkeys[6]:"%5.3f"%(self.par["ts_Vthreshold"]), \
                  pkeys[7]:"%5.3f"%(self.par["ts_Vwidth"]), \
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
          self.par["ts_Vwidth"] = float(pdlg.data[pkeys[7] ])
          pdlg.destroy()

     def OnSetMovingAvg(self):
          # show and change parameters
          # parameter keys
          pkeys=['N Average']
          data = {pkeys[0]:"%d"%(int(self.par["Navg"]))}
          pdlg = NumberDialog(data, self, title="Points for Moving Average", \
                              labels = pkeys, \
                                  keys = pkeys, \
                                      about_txt = "Set number of points")
          pdlg.exec_()
          # now set the new parameters
          self.par['Navg']=int(pdlg.data[pkeys[0]])
          print(f'Set Navg = {self.par["Navg"]}')

     def OnSetTimeOffset(self):
          # show and change parameters
          # parameter keys
          pkeys=['Time offset']
          data = {pkeys[0]:"%f"%(self.par["toffset"])}
          pdlg = NumberDialog(data, self, title="Set time offset", \
                              labels = pkeys, \
                                  keys = pkeys, \
                                      about_txt = "enter the time offset for time axis")
          pdlg.exec_()
          # now set the new parameters
          self.par['toffset']=float(pdlg.data[pkeys[0]])
          print(f'Set toffset = {self.par["toffset"]}')
          self.load_data()
          self.select_data()
          
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
          if not self.fft.ok:
               print("No power spectrum to plot")
               return
          # plot power spectrum
          # set initiall axis limits:
          self.axes.set_xlim(self.fft.f.min(), self.fft.f.max())
          self.toolbar.set_data(self.fft.f, self.fft.logp)
          self.toolbar.xlabel = 'f (Hz)'
          self.toolbar.ylabel = 'PS'
          self.toolbar.thinning('log(Power)')
          """
          if self.par["draw_lines"]:
               self.axes.plot(self.fft.f[::skip],self.fft.logp[::skip])
          else:
               self.axes.plot(self.fft.f[::skip],self.fft.logp[::skip], '.')
          # self.axes.set_yscale('log')
          self.axes.set_xlabel('f')
          self.axes.set_ylabel('Power')
          self.axes.set_title('log(Power)')
          self.figure_canvas.draw()
          """
          
     def OnFFTplotan(self):
          # plot sin coeff
          if not self.fft.ok:
               print("No power spectrum to plot")
               return
          self.axes.set_xlim(self.fft.f.min(), self.fft.f.max())
          self.toolbar.set_data(self.fft.f, np.log10(np.abs(self.fft.an)))
          self.toolbar.xlabel = 'f (Hz)'
          self.toolbar.ylabel = 'log(abs(an))'
          self.toolbar.thinning('an')
          
     def OnFFTplotbn(self):
          # plot cos coeff
          if not self.fft.ok:
               print("No power spectrum to plot")
               return
          self.axes.set_xlim(self.fft.f.min(), self.fft.f.max())
          self.toolbar.set_data(self.fft.f, np.log10(np.abs(self.fft.bn)))
          self.toolbar.xlabel = 'f (Hz)'
          self.toolbar.ylabel = 'log(abs(bn))'
          self.toolbar.thinning('bn')

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
          self.Vcorr = self.fft.inv_transform()
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

     def OnUsenpz_file(self):
          self.par["use_npz_file"] = True
          self.par["use_NI_file"] = False
          self.par["use_GaGe_file"] = False
          print(("use npz data = ", self.par["use_npz_file"]))
          print(("use NI data = ", self.par["use_NI_file"]))
          print(("use GaGe data = ", self.par["use_GaGe_file"]))
          
     def OnUseNI_file(self):
          self.par["use_npz_file"] = False
          self.par["use_NI_file"] = True
          self.par["use_GaGe_file"] = False
          print(("use npz data = ", self.par["use_npz_file"]))
          print(("use NI data = ", self.par["use_NI_file"]))
          print(("use GaGe data = ", self.par["use_GaGe_file"]))

     def OnUseGaGe_file(self):
          self.par["use_npz_file"] = False
          self.par["use_NI_file"] = False
          self.par["use_GaGe_file"] = True
          print(("use npz data = ", self.par["use_npz_file"]))
          print(("use NI data = ", self.par["use_NI_file"]))
          print(("use GaGe data = ", self.par["use_GaGe_file"]))

     def OnToggleConvertInteger(self):
          self.par["convert_int"] = not self.par["convert_int"]
          print(("convert integers = ", self.par["convert_int"]))

     def OnToggleAutoTimeRange(self):
          self.par["auto_time_range"] = not self.par["auto_time_range"]
          print(("auto_time_range = ", self.par["auto_time_range"]))

     def OnTogglelimits(self):
          self.par["use_limits"] = not self.par["use_limits"]
          print(("use limits = ", self.par["use_limits"]))

     def OnToggleHistolimits(self,event):
          self.par["auto_histo"] = not self.par["auto_histo"]
          print(("auto histo = ", self.par["auto_histo"]))

     def OnToggleHistoPoints(self,event):
          self.par["plot_histo_points"] = not self.par["plot_histo_points"]
          print(("Plot histogram points = ", self.par["plot_histo_points"]))


     def OnToggleLines(self, event):
          self.par["draw_lines"] = not self.par["draw_lines"]
          if not (self.par["draw_lines"] or self.par["draw_points"]):
              self.par["draw_points"] = True
              print('----> force draw_points = True')
              self.update_options_menu()
          try:
              self.toolbar.thinning()
          except:
              'Thinning did not work'

     def OnTogglePoints(self, event):
          self.par["draw_points"] = not self.par["draw_points"]
          if not (self.par["draw_lines"] or self.par["draw_points"]):
              self.par["draw_points"] = True
              print('----> force draw_points = True')
              self.update_options_menu()
          print(("draw points = ", self.par["draw_points"]))
          try:
              self.toolbar.thinning()
          except:
              'Thinning did not work' 

     def OnToggleHidePeaks(self, event):
          self.par["hide_peaks"] = not self.par["hide_peaks"]
          print(("hide peaks = ", self.par["hide_peaks"]))
          if self.par["hide_peaks"]:
              self.toolbar.plot_peaks = False
          else:
              self.toolbar.plot_peaks = True

              
     def OnUseCorrected(self, event):
          self.par["corrected"] = not self.par["corrected"]
          print(("use corrected = ", self.par["corrected"]))
          
     def OnNothing(self):
          print("do nothing")

#----------------------------------------------------------------------
# main program
#----------------------------------------------------------------------

         
if __name__ == '__main__':
    parser = AG.ArgumentParser()
    parser.add_argument('-c', '--channels', help="list of channel names for which a window will be opened. Separate the names with a / e.g. ch1/ch2/ch3 (no spaces)", default = 'one_channel')
    args = parser.parse_args()
    #Create App
    # check if app already exists
    app = QtCore.QCoreApplication.instance()
    
    if app is None:
        # if not create it
        app = QtWidgets.QApplication(sys.argv)
    
    # add PlotFrames to app
    
    frame_names = args.channels.split('/')
    
    plot_frames = [PlotFrame(app, key = k) for k in frame_names]
    
    frames = dict(zip(frame_names, plot_frames))

    sys.exit(app.exec_())
     
     
