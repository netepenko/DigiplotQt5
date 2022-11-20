# -*- coding: utf-8 -*-
"""
Created on Sat Mar 09 16:14:31 2019



@author: Alex, WB
"""
from matplotlib.backends.backend_qt5 import NavigationToolbar2QT
from PyQt5 import QtWidgets, QtCore, QtGui
import numpy as np


class NavigationToolbar(NavigationToolbar2QT):
    def __init__(self, canvas, parent, coordinates=True):
        NavigationToolbar2QT.__init__(self, canvas, parent, coordinates=True)
        self.parent = parent
        self.axes = parent.axes
        self.N = 1000
        self.t=[]
        self.V=[]
        self.ch_n=0
        self.fn='' # fila name array
        self.Ncontrol= QtWidgets.QLineEdit(self)
        self.Ncontrol.setValidator(QtGui.QIntValidator(self.Ncontrol))
        self.Ncontrol.setFixedWidth(50)
        self.Ncontrol.setText(str(self.N))
        self.Nlabel=QtWidgets.QLabel('Data points on Fig, N=', self)
        
        
        self.refr=QtWidgets.QPushButton(QtGui.QIcon('refresh.png'), None, self)
        self.refr.clicked.connect(self.getN)
          
        self.addWidget(self.Nlabel)
        self.addWidget(self.Ncontrol)
        self.addWidget(self.refr)
    # to remove matplotlib toolbar status line
    def set_message(self, msg):
        pass
    
    def getN(self):
      self.N = int(self.Ncontrol.text())
      try:
          self.axes.set_autoscaley_on(True)
          self.thinning()
      except:
          print("Thinnning didn't work properly.")
          
    def back(self, *args):
        NavigationToolbar2QT.back(self, *args)
        self.thinning()
        
    def forward(self, *args):
        NavigationToolbar2QT.forward(self, *args)
        self.thinning()
        
    def home(self, *args):
        NavigationToolbar2QT.home(self, *args)
        self.thinning()
        
    def release_pan(self, event):
        NavigationToolbar2QT.release_pan(self, event)
        self.thinning()
        
    def release_zoom(self, event):
        NavigationToolbar2QT.release_zoom(self, event)
        self.thinning()    
        
    #rescale and replot
    # this function drops unnecessary ponts in ploting region to reduce
    # plotting delays, but keeps everything above threshold to see peaks
    
    def thinning(self):
        # get current axes
        self.axes = self.parent.axes
        ax = self.axes

        # get current axes limits
        (xmin, xmax)=ax.get_xlim()
        (ymin, ymax)=ax.get_ylim()

        t=self.t
        V=self.V

        #get original data which is in plotting region
        org_points=np.where((xmin < t) & (t < xmax))
        to=t[org_points]
        Vo=V[org_points]
        
        n=1 #take every nth point from data
        if len(to)> self.N: # was: len(to)/k
            n=int(len(to)/self.N) # was: len(to)/k/self.N)
        tcut=to[::n]
        Vcut=Vo[::n]
        ax.set_prop_cycle(None)
        # cleear the data points (daters than removing data points and lines)
        ax.clear()
        # set the new limits careful with units here
        ax.set_xlim(max(tcut[0]-1e-6, xmin), min(tcut[-1] + 1e-6, xmax))
        ax.autoscale(enable=False, axis='x', tight=False)
        
        if self.parent.par['draw_lines']:
            all_plot_new = ax.plot(tcut,Vcut, self.mker + '-',  label='%s Ch %d'%(self.fn,self.ch_n), alpha=0.5, color=self.colors[self.ch_n])
        else:
            all_plot_new = ax.plot(tcut,Vcut, self.mker, label='%s Ch %d'%(self.fn,self.ch_n), alpha=0.5, color=self.colors[self.ch_n])
        Vc=Vcut[np.where((xmin < tcut) & (tcut < xmax))]
        Vcmin=Vc.min()
        Vcmax=Vc.max()
        ax.set_ylim(min((Vcmin - 0.1), ymin), max(Vcmax + 0.1,ymax))
        ax.set_autoscaley_on(False)
        self.parent.all_plot=all_plot_new
        ax.legend()
        ax.set_xlabel('t')
        ax.set_ylabel('V')
        ax.figure.canvas.draw()



class NumberDialog(QtWidgets.QDialog):
     def __init__(self, data, parent = None, title = 'Enter parameters', labels=['First Label','Second Label'],
                  keys=['first','second'], about_txt = 'about_txt'):
          
          self.parent=parent
          self.keys=keys
          self.data=data
          QtWidgets.QDialog.__init__(self, parent)
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
               if title == 'Time axis' and key =='t start':
                   Val_l=QtWidgets.QPushButton(labels[i],self)
                   Val_l.clicked.connect(self.getStart)
               else:
                   if title == 'Time axis' and key=='t stop':
                       Val_l=QtWidgets.QPushButton(labels[i],self)
                       Val_l.clicked.connect(self.getStop)
                   else:
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
     def getStart(self):
         xst=self.parent.axes.get_xlim()[0]
         self.qle['t start'].setText(str(xst))
     def getStop(self):
         xtp=self.parent.axes.get_xlim()[1]
         self.qle['t stop'].setText(str(xtp))
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