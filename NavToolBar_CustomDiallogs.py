# -*- coding: utf-8 -*-
"""
Created on Sat Mar 09 16:14:31 2019



@author: Alex, WB
"""
from matplotlib.backends.backend_qt5 import NavigationToolbar2QT
from PyQt5 import QtWidgets, QtCore, QtGui
from PyQt5.QtWidgets import QDialog

# this is maintained by QtDesigner
import LimitsListWidget as LLW 

import numpy as np


class NavigationToolbar(NavigationToolbar2QT):
    def __init__(self, canvas, parent, coordinates=True):
        NavigationToolbar2QT.__init__(self, canvas, parent, coordinates=True)
        self.parent = parent
        self.axes = parent.axes
        self.N = 1000
        #self.t=[]
        #self.V=[]
        self.ch_n=0
        self.x=[]
        self.y=[]
        self.fn='' # file name array
        self.xlabel = 't'
        self.ylabel = 'V'
        self.xpeak = []
        self.ypeak = []
        self.plot_peaks = False
        # enter the number if data points to show
        self.Ncontrol= QtWidgets.QLineEdit(self)
        self.Ncontrol.setValidator(QtGui.QIntValidator(self.Ncontrol))
        self.Ncontrol.setFixedWidth(50)
        self.Ncontrol.setText(str(self.N))
        self.Nlabel=QtWidgets.QLabel('Data points on Fig, N=', self)
        
        # refresh button
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
        
    def set_data(self, x, y):
        self.x = x
        self.y = y
        
    def set_peakdata(self, x, y):
        self.xpeak = x
        self.ypeak = y
        
    #rescale and replot
    # this function drops unnecessary ponts in ploting region to reduce
    # plotting delays, but keeps everything above threshold to see peaks
    
    def thinning(self, title = ''):
        # get current axes
        if self.x == []:
            print('Sorry not data in toolbar for thinned plotting, issue plot command first !')
            return
        self.axes = self.parent.axes
        ax = self.axes
        print(f'N = {self.N}')
        
        # get current axes limits
        (xmin, xmax)=ax.get_xlim()
        (ymin, ymax)=ax.get_ylim()

        print(f'thinning x - limits : {xmin},{xmax}')


        x=self.x
        y=self.y

        #get original data which is in plotting region
        org_points=(xmin < x) & (x < xmax)
        try:
            if not org_points.max():
                print(f'No points satisfying: ({xmin} < x) & (x < {xmax})!')
                return
        except Exception as err:
            print(70*'-')
            print(70*'-')
            print(f'Error occurred : {err}\n Current values {ax.get_xlim()}, {ax.get_ylim()}, x = {x}\n most likely no data to show')
            print(70*'-')
            print(70*'-')
            return
        xo=x[org_points]
        yo=y[org_points]

        
        n=1 #take every nth point from data
        if len(xo)> self.N: # was: len(to)/k
            n=int(len(xo)/self.N) # was: len(to)/k/self.N)
        xcut=xo[::n]
        ycut=yo[::n]
        if len(xcut) == 0:
            print('no tcut data points left to plot after thinning')
            return
        ax.set_prop_cycle(None)
        # clear the data points (daters than removing data points and lines)
        ax.clear()
        # set the new limits careful with units here
        xtol = 0. #1e-6
        ax.set_xlim(max(xcut[0]-xtol, xmin), min(xcut[-1] + xtol, xmax))
        ax.autoscale(enable=False, axis='x', tight=False)
        l_text = '%s Ch %d'%(self.fn,self.ch_n)
        print('Label data = ', self.fn, self.ch_n)
        if self.parent.par['draw_lines'] and (not self.parent.par['draw_points']):    
            all_plot_new = ax.plot(xcut,ycut, '-',  label=l_text, alpha=0.5, color=self.colors[self.ch_n])
        elif self.parent.par['draw_lines'] and self.parent.par['draw_points']:
            all_plot_new = ax.plot(xcut,ycut, self.mker + '-',  label=l_text, alpha=0.5, color=self.colors[self.ch_n])
        else:
            all_plot_new = ax.plot(xcut,ycut, self.mker, label=l_text, alpha=0.5, color=self.colors[self.ch_n])
        # plot peaks if selected
        if self.plot_peaks:
            peak_plot = ax.plot(self.xpeak, self.ypeak, 'x' , color='m')
        sel_y = (xmin < xcut) & (xcut < xmax)
        if not sel_y.max():
            print(f'No points satisfying: ({xmin} < xcut) & (xcut < {xmax})!') 
            return
        yc=ycut[sel_y]
        ycmin=yc.min()
        ycmax=yc.max()
        ytol = 0. # 0.1
        ax.set_ylim(min((ycmin - ytol), ymin), max(ycmax + ytol,ymax))
        ax.set_autoscaley_on(False)
        self.parent.all_plot=all_plot_new
        ax.legend()
        ax.set_xlabel(self.xlabel)
        ax.set_ylabel(self.ylabel)
        if title != '':
            ax.set_title(title)
        ax.figure.canvas.draw()



class TextDialog(QtWidgets.QDialog):
     def __init__(self, data, parent = None, title = 'Enter Text', labels=['First Label','Second Label'],
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
               Val_l  = QtWidgets.QLabel(labels[i], self)
               Val_t  = QtWidgets.QLineEdit(self)
               Val_t.setText(data.get(key))
               self.layout.addWidget(Val_l, i+2, 0)
               self.layout.addWidget(Val_t, i+2, 1)
               self.qle[key]=Val_t                       
          self.layout.addWidget(ok, nrow+2, 0)
          self.layout.addWidget(cancel, nrow+2, 1)

        
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
         
class LimitsListDialog(QDialog):
    def __init__(self, limits):
        # limits is an array of limit dictionaries
        super().__init__()
        self.ui = LLW.Ui_Dialog()
        self.ui.setupUi(self)
        self.ui.listWidgetLimits.itemClicked.connect(self.getSelectedLimit)
        self.init_list(limits)
        self.selected_index = 0
        self.selected_text = 'default'
        self.show()   
        
    def init_list(self, limits):
        for d in limits:
            xmin, xmax = d['value']
            comment = d['comment']
            list_entry = f'({xmin:.2e}, {xmax:.2e}) / {comment}'
            print(f'adding to list : {list_entry}')
            # add to list
            item = QtWidgets.QListWidgetItem()
            item.setText(list_entry)
            self.ui.listWidgetLimits.addItem(item)

    def getSelectedLimit(self):
        self.ui.limitsChosenLabel.setText("You have selected"+self.ui.listWidgetLimits.currentItem().text())
        self.selected_text = self.ui.listWidgetLimits.currentItem().text()
        self.selected_index = self.ui.listWidgetLimits.currentRow()
        # print (f'{self.selected_text}, {self.selected_index}')

