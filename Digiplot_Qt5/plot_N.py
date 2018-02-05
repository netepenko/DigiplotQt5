# -*- coding: utf-8 -*-
"""
Created on Mon Mar 27 19:03:39 2017

@author: Alex
"""

# Plot routine which shows only N or less points on figure by skipping points
# from originally enetered data. Adds or removes points to keep it N or less 
# when zooming in and out and switching forward or backward in views.
# This versiton is specific for digiplot_PD_fast

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
    if len(to)> self.N:
        n=int(len(to)/ self.N)
    tcut=to[::n]
    Vcut=Vo[::n]
    
    pl.cla()
    ax.set_autoscale_on(False)
    ax.plot(tcut,Vcut,*rest)
    self.tN=t[::n]
    self.VN=V[::n]
    
    
def plot_N(*args, **kwargs):
    N=kwargs.pop('N', 10000)
    
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


    washold = ax.ishold()
    hold = kwargs.pop('hold', None)
    if hold is not None:
        ax.hold(hold)
    try:
        ret = ax.plot(tcut, Vcut, *rest, **kwargs)
    finally:
        ax.hold(washold)
    return ret 

    
    
        
