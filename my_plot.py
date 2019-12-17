# -*- coding: utf-8 -*-
"""
Created on Tue Oct 25 14:43:47 2016

@author: Alex
"""
import matplotlib
import matplotlib.pyplot as pl
from matplotlib.figure import Figure

from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg
from matplotlib.backends.backend_qt5 import FigureManagerQT, NavigationToolbar2QT

import numpy as np


#modified figuremaneger class which uses modified toolbar class
class myFigureManagerQT(FigureManagerQT):
    def __init__(self, canvas, num):
        FigureManagerQT.__init__(self, canvas, num)
        
    def _get_toolbar(self, canvas, parent):
        # must be inited after the window, drawingArea and figure
        # attrs are set
        if matplotlib.rcParams['toolbar'] == 'toolbar2':
            toolbar = myNavigationToolbar2QT(canvas, parent, False)
        else:
            toolbar = None
        return toolbar

#modified toolbar class with rescaling added on events
class myNavigationToolbar2QT(NavigationToolbar2QT):
    def __init__(self, canvas, parent, coordinates=True):
        NavigationToolbar2QT.__init__(self, canvas, parent, coordinates=True)
        self.myinitdata = 0.
        self.N = 100
        self.tN=0
        self.VN=0
        
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
        self.my_rescale()      #rescale and replot  
#        
#        
    def my_rescale(self):
        
        t=self.myinitdata[0]
        V=self.myinitdata[1]
        rest=self.myinitdata[2:]
        ax=pl.gca()
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
        
        pl.cla()
        ax.set_autoscale_on(False)
        ax.plot(tcut,Vcut,*rest)
        self.tN=t[::n]
        self.VN=V[::n]
        
    
        
        

def my_figure_manager(num, *args, **kwargs):
    """
    Create a new figure manager instance
    """
    FigureClass = kwargs.pop('FigureClass', Figure)
    thisFig = FigureClass(*args, **kwargs)
    return new_figure_manager_given_figure(num, thisFig)


def new_figure_manager_given_figure(num, figure):
    """
    Create a new figure manager instance for the given figure.
    """
    canvas = FigureCanvasQTAgg(figure)
    return myFigureManagerQT(canvas, num)

    
   
# pl plot modification wich plots maximum N points
# N can be passed as a kwarg
def my_plot(*args, **kwargs):
    N=kwargs.pop('N', 10000)
    pl.new_figure_manager = my_figure_manager
    ax = pl.gca()
    # allow callers to override the hold state by passing hold=True|False
    
    #saves input data into toolbar class object and reuses later
    figManager = pl.get_current_fig_manager()
    figManager.toolbar.myinitdata = args
    figManager.toolbar.N = N
    
    #dropping points for first plot
    t=args[0]
    V=args[1]
    rest=args[2:]
    
    n=1 #take every nth point from data
    if len(t)> N:
        n=int(len(t)/N)
    tcut=t[::n]
    Vcut=V[::n]


    try:
        ret = ax.plot(tcut, Vcut, *rest, **kwargs)
        return ret 
    except:
        print("Plotting faild.")
    

    
    
        
#for testing        
my_plot(np.array(list(range(10000000))),np.array(list(range(10000000))),'.', N=10000)
#pl.plot(np.array(range(1000)),np.array(range(1000)),'.')

    