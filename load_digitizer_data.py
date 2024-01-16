#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov  1 13:26:13 2023

@author: boeglinw
"""
import numpy as np
import h5py  # HDF file handling
import os
import sys

# conversion to microseconds constant
us = 1.e6

to_bool = {'True':True, 'False':False, '1':True,'0':False}



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


def get_window_slice(xmin, x, xmax):
        # fast function to 
        # get the slice corresponding to xmin and xmax in x
        # the x-values need to be equally spaced        
        dx = x[1] - x[0]
        nmin = max( 0, int(round( (xmin - x[0])/dx )))
        nmax = min( (int(round ((xmax - x[0])/dx ))), len(x) -1  )
        return slice(nmin, nmax + 1)

class digi_data:
    
    def __init__(self, filename, datatype = 'Gage', 
                 convert = False, 
                 channel = 0,
                 t_offset = 0.,
                 tmin = -5.e6,
                 tmax = 10.e6):
        
        # all times should be in us
        
        self.data_filename = filename
        self.data_type = datatype
        self.channel = channel
        self.t_offset = t_offset
        self.tmin = tmin
        self.tmax = tmax
        self.dtmin = tmin
        self.dtmax = tmax
        self.convert = convert
        # dictionary of load data functions
        self.load_data_dict = {'Gage':self.load_gage_data,
                               'NI':self.load_NI_data,
                          'corrected':self.load_hdf_data,
                          'filtered':self.load_npz_data}
    def show_data_types(self):
        print(self.load_data_dict.keys())
        
    def load_raw_data(self, **kwargs):
        return self.load_data_dict[self.data_type]()

    
    def load_NI_data(self, **kwargs):
        # --------------------------------
        # ######## Load raw data #########
        f = h5py.File(self.data_filename, 'r')
        # setup reading the data
        data_root = 'wfm_group0/traces/trace' + str(self.channel) + '/'

        print("-----------------------Getting NI data------------------------")

        # load time information
        t0 = f[data_root + 'x-axis'].attrs['start']*us + self.t_offset
        dt = f[data_root + 'x-axis'].attrs['increment']*us
        # load scale coeeff and scale dataset
        scale = f[data_root + 'y-axis/scale_coef'][()]

        # get the y dataset length
        nall = f[data_root + 'y-axis/data_vector/data'].shape[0]

        # make time array based on number of points in y data
        tall = t0 + dt*np.arange(nall, dtype=float)
        self.tall = tall

        # get the y dataset (measured data)
        if self.convert:
            ydata = f[data_root + 'y-axis/data_vector/data'][()].astype('int16')
        else:
            ydata = f[data_root + 'y-axis/data_vector/data'][()]

        # calculate voltage for dataset
        V = scale[0] + scale[1]*ydata
        print("-----------------------Data loaded-------------------------")

        # check time limits
        self.check_limits(tall)
        # select time range
        t_sl = get_window_slice(self.tmin, tall, self.tmax)
        """
        i_min = max(0,   int((self.tmin - tall[0])/dt))
        i_max = min(nall,int((self.tmax - tall[0])/dt) + 1)
        print(f'setting i_min =  {i_min} and i_max = {i_max}')
        self.i_min = i_min
        self.i_max = i_max
        """
        print(f'setting time slice =  {t_sl.start} and i_max = {t_sl.stop}')
        # store data
        self.td = tall[t_sl]  # time data (microseconds) in analysis interval
        self.Vps = V[t_sl]  # voltage data
        self.dt = dt  # time step (microseconds)



    def load_gage_data(self, **kwargs):
        # ######## Load raw data #########
        f = h5py.File(self.data_filename, 'r')
        
        print("-----------------------Getting Gage data------------------------")
        
        # load gage digitizer data
        data_root = 'wfm_group0/traces/trace' + str(self.channel)  + '/'
        
        xaxis_access = data_root + 'x-axis'
        t0 = f[xaxis_access].attrs['start']*us + self.t_offset
        dt = f[xaxis_access].attrs['increment']*us
        # measured data
        # scale dataset
        yaxis_scale_access = data_root + 'y-axis/scale_coef'
        yaxis_data_access = data_root + 'y-axis/data_vector/data'
        # there are 2 values now
        #self.scale = self.f[yaxis_scale_access].value
        scale = f[yaxis_scale_access][0]
         
        # WB temp solution 8/13/13
        # set the type as 16 bit signed, this is not good the data type should come from
        # the hdf file
        if self.convert:
             ydata = f[data_root + 'y-axis/data_vector/data'].astype('int16')[:]
        else:
             ydata = f[data_root + 'y-axis/data_vector/data'][:]
        # get the y dataset length
        nall = ydata.shape[0]          
        #added on 3/3/2021 for debugging             
        for i, yd in enumerate(ydata[:10]):
            print(f'i = {i}, ydata = {yd}')
        
        # make the time axis
        print("Calculate t data")
        tall = t0 + dt*np.arange(ydata.shape[0], dtype = float)
        # vself.tall = tall
        print("-----------------------Data loaded-------------------------")
        print(f"t0 = {t0}, tall[0] = {tall[0]}, tall[-1] = {tall[-1]}, t_offset = {self.t_offset}")
        #self.td  = tall
        #self.Vps = ydata
        # check time limits
        self.check_time_limits(tall)
        # select time range
        t_sl = get_window_slice(self.tmin, tall, self.tmax)
        """
        i_min = max(0,   int((self.tmin - tall[0])/dt))
        i_max = min(nall,int((self.tmax - tall[0])/dt) + 1)
        print(f'setting i_min =  {i_min} and i_max = {i_max}')
        print('slices = ',i_min, i_max)
        """
        self.i_min = t_sl.start
        self.i_max = t_sl.stop
        print(f'setting time slice =  {t_sl.start} and i_max = {t_sl.stop}')
        # store data
        td = tall[t_sl]  # time data (microseconds) in analysis interval
        Vps = ydata[t_sl]  # voltage data
        print(f'finished loading data for channel {self.channel}')
        self.td = td
        self.Vps = Vps
        self.dt = dt


    def check_time_limits(self, t):
        # make sure limits are within time range, if not adjust them
        self.tmin = max(t[0], self.tmin)
        self.tmax = min(t[-1], self.tmax)
        
        

        
    def get_data(self):
        return self.td, self.Vps


#   plotting of raw data without overloading the figure
# (skipping some data points according to maximum allowed points on plot)

    def load_npz_data(self):
        # --------------------------------
        # ######## Load raw data #########

        self.f = np.load(self.data_filename)
        d = self.f
        print("--------------------- Get npz data ------------------------")
        self.td = d['time']*us
        self.Vps = d['signal']
        self.dt = d['time'][1] - d['time'][0]
        print("-----------------------Data loaded-------------------------")
        # add pulser to data if add_pulser parameter set to True



    def load_hdf_data(self):
        # --------------------------------
        # ######## Load raw data #########
        print("--------------------- Get hdf data ------------------------")
        f = h5py.File(self.data_filename, 'r')
        # load the data set for corrected data
        Vc = f['V_corr']
        scale = Vc.attrs['V_corr_scale']
        Vd = Vc[:]/scale
        dt = Vc.attrs['dt']
        t0 = Vc.attrs['t0']
        td = t0 + dt*np.arange(Vd.shape[0], dtype=float)
        nall = Vd.shape[0]
        # select time range
        i_min = max(0,int(td[0] - self.tmin)/dt)
        i_max = min(nall,int(td[-1] - self.tmin)/dt)
        # check time limits
        self.tmin = max(td[0], self.tmin)
        self.tmax = min(td[-1], self.tmax)
        # select time range
        i_min = max(0,   int((td[0] - self.tmin)/dt))
        i_max = min(nall,int((td[-1] - self.tmin)/dt))
        #
        self.tmin = max(td[0], self.tmin)
        self.tmax = min(td[-1], self.tmax)
        # select time range
        i_min = max(0,   int((self.tmin - td[0])/dt))
        i_max = min(nall,int((self.tmax - td[0])/dt) + 1)
        self.i_min = i_min
        self.i_max = i_max
        print(f'setting i_min =  {i_min} and i_max = {i_max}')
        print('slices = ',i_min, i_max)
        # store data
        self.td = td[i_min:i_max]  # time data (microseconds) in analysis interval
        self.Vps = Vd[i_min:i_max]  # voltage data
        self.dt = dt  # time step (microseconds)
        print("-----------------------Data loaded-------------------------")
        # add pulser to data if add_pulser parameter set to True


    def get_time_window(self, tmin, tmax):
         # find range of time values between tmin and tmax
         nmin = int( (tmin - self.t0)/self.dt )
         nmax = min( (int( (tmax - self.t0)/self.dt ) + 1), self.nall -1 )
         return slice(nmin, nmax + 1)
      
    def moving_average(self, n = 11):
        self.Vps = moving_average(self.Vps, n)


