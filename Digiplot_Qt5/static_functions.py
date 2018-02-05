# -*- coding: utf-8 -*-
"""
Created on Sat Mar 25 20:02:29 2017

@author: Alex
"""
import numpy as np
import ffind_peaks as FP     
def get_window(xmin, x, xmax):
     # get the slice corresponding to xmin and xmax in x
     # the x-values need to equally spaced
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
          print "No limits present, analyze all data"
          results = np.zeros((2,), dtype = 'int32')
          pmin = np.zeros((len(yval)/5, ), dtype='int32')
          pmax = np.zeros((len(yval)/5, ), dtype='int32')
          # try:
          R = FP.find_peaks(len(yval), ystep, yval, results, pmin, pmax)
          #except:
          #print "problem with peak finding"
          #return []
          nmin = results[0]
          nmax = results[1]
     else:
          # get the window
          print "Analyze data between ", xmin, " and ", xmax
          sl = get_window( xmin, xval, xmax)
          # progress dialog
          results = np.zeros((2,), dtype = 'int32')
          pmin = np.zeros((len(yval[sl])/5, ), dtype='int32')
          pmax = np.zeros((len(yval[sl])/5, ), dtype='int32')
          try:
              R = FP.find_peaks(len(yval[sl]), ystep, yval[sl], results, pmin, pmax)
          except:
               print "problem with peak finding"
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