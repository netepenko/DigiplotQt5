#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 14 20:45:48 2023

Adapt mesured nosise to data and subtract from signal data

all times are in us

@author: boeglinw
"""

import numpy as np

import load_digitizer_data as LD
import correlation_lags as CL

from scipy import signal


def shift_slice(sl, shift):
    return slice(sl.start+shift, sl.stop+shift) 


h_line = 70*'-'
us = 1e6

class bkg:
    
    def __init__(self, file_name, data_dir, shot = -1, \
                 t_min = 0., t_max = 1e6, t_offset = 0, \
                 data_channel = 0, \
                 bkg_channel = 4,\
                 moving_average = 21, \
                 niter = 5,\
                 time_window = 200.,\
                     ):
        """
        setup background data for subtraction

        Parameters
        ----------
        file_name : str
            file name containing the data.
        data_dir : str
            directory for the data
        t_min : float, optional
            start time of data (in us). The default is 0..
        t_max : fload, optional
            end time of data (in us). The default is 1e6.
        t_offset : float, optional
            time offset in us to be added to the exp. time data. The default is 0.
        channel: int, optional
            channel number for the noise channel. The default is 4.
        moving_average: int
            window size in data points for moving average (odd number), The default is 21
        niter: int
            number of moving average iterations. THe default is 5
        time_window: float
            time window width used to perform the correction. THe default is 200

        Returns
        -------
        None.

        """
        self.file_name = file_name
        self.data_dir = data_dir
        self.tmin = t_min
        self.tmax = t_max
        self.t_offset = t_offset
        self.shot = shot
        self.data_channel = data_channel
        self.bkg_channel = bkg_channel
        self.moving_avg = moving_average
        self.niter = niter
        self.delta_t = time_window
        self.d_data  = None
        self.d_bkg = None
        
    def load_data(self):
        if self.shot == -1:
            # this assumes the standard file name format: DAQ_49182_231218_150221.hws
            self.shot = int(self.file_name.split('_')[1])
        self.d_data = LD.digi_data(self.data_dir + self.file_name, 
                                   channel = self.data_channel, 
                                   tmin = self.tmin, 
                                   tmax = self.tmax,
                                   t_offset = self.t_offset)
        self.d_data.load_raw_data()
        self.dt = self.d_data.dt
        self.d_bkg = LD.digi_data(self.data_dir + self.file_name, 
                                  channel = self.bkg_channel, 
                                  tmin = self.tmin, 
                                  tmax = self.tmax,
                                  t_offset = self.t_offset)
        self.d_bkg.load_raw_data()
        
    def moving_average(self):
        if self.d_data is None:
            print(' no data loaded, nothing to do !')
            return
        # calculate the moving averages for data and bkg
        print(h_line)
        print('Start moving average for data:')
        for i in range(self.niter):
            self.d_data.moving_average(self.moving_avg)
        print('Finished moving average for data:')
        print('Start moving average for bkg:')
        for i in range(self.niter):
            self.d_bkg.moving_average(self.moving_avg)
        print('Finished moving average for bkg:')

    def correct(self, clear_data = False):
        if self.d_data is None:
            print(' no data loaded, nothing to do !')
            return        
        #loop over all data and subtract noise
        dd = self.d_data
        dn = self.d_bkg
        # create an arrau of slices
        n_slice = int(self.delta_t//dd.dt)
        
        i_start = np.arange(0, dd.td.shape[0], n_slice)  # starting index
        i_end = np.roll(i_start, -1)                        # stopping index
        
        # create the slice array
        slices = [slice(ss,ee) for ss, ee in zip(i_start, i_end)][:-1]  # skip the last slice
        
        V_corr = np.zeros_like(dd.Vps)
        
        n_slice = len(slices)
        print(f'Analyzing {n_slice} slices ')
        
        for i,sel in enumerate(slices):
         
            if not i%100:
                print(f'working on slice {i} {(i/n_slice*100):.1f} %')
            Vps_loc = dd.Vps[sel]
            Vn_loc = dn.Vps[sel]
            corr = signal.correlate(Vps_loc, Vn_loc, mode = 'full')
            lags = CL.correlation_lags(Vps_loc.size, Vn_loc.size, mode="full")
            
            lag = lags[np.argmax(corr)]    
            sel_r = slice(sel.start - lag, sel.stop - lag)
            # shift noise to align with data
            #Vnr = np.roll(dn.Vps, lag)
            Vnr = dn.Vps[sel_r]
            # calculate optimal scaling factor
            #a = np.sum(dd.Vps[sel]*Vnr[sel])/np.sum(Vnr[sel]**2)
            a = np.sum(dd.Vps[sel]*Vnr)/np.sum(Vnr**2)
            
            #V_sig[sel] = dd.Vps[sel] - a*Vnr[sel]
            V_corr[sel] = dd.Vps[sel] - a*Vnr
        if clear_data :
            self.d_data = None
            self.d_bkg = None
        print('Correction completed !')
        return self.d_data.td, V_corr
