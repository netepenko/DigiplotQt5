# FFT for digitzer data

import scipy.fftpack as FFTP
import numpy as np
import pdb

def get_window(xmin, x, xmax):
     # get the slice corresponding to xmin and xmax in x
     # the x-values need to equally spaced
     dx = x[1] - x[0]
     nmin = int( (xmin - x[0])/dx )
     nmax = min( (int( (xmax - x[0])/dx ) + 1), len(x) -1 ) + 1
     # make sure that there are an even number of points
     nmax = nmax - (nmax - nmin)%2
     return slice(nmin, nmax)

def z_value(z, mag_val):
    if (z.real == 0):
        return mag_val*(0. + 1.j)
    elif( z.imag == 0.):
        return mag_val*(1. + 0.j)
    else:
        phi = np.arctan(z.imag/z.real)
        return mag_val*np.exp(phi*1j)
    # done
    
set_z_value = np.vectorize(z_value)
minval = 1.e-30
maxval = 1.e30

class FFT:
     def __init__(self, t = None):
          self.par = {}
          self.par["xmin"] = -1.e30
          self.par["xmax"] = 1.e30
          self.par["alpha"] = 1.0
          self.par["rthresh"] = 1e-30
          self.par["limits"] = False
          if t is None:
               self.ok = False
               return
          # frequencies, make sure the numer of data points is even
          self.N = 2* int( len(t)/2.)
          # find the time interval
          self.T = t[self.N-1] - t[0]
          self.f = np.arange(0, self.N/2)/self.T
          self.par["xmin"] = self.f[0]
          self.par["xmax"] = self.f[-1]
          self.ok = False
          
     def transform( self, V):
          # calculate the FFT oeff.
          print("Calculate FFT")
          self.cn = FFTP.fft(V[:self.N])
          self.an0 = self.cn[0]
          # sin coefs
          self.an = -2./self.N*np.imag(self.cn[0:self.N/2])
          # cos coefs
          self.bn = 2./self.N*np.real(self.cn[self.N/2:])
          # reverse bn for future work and to aplly cuts
          self.bnr = self.bn[::-1]
          self.rr = np.sqrt(self.an**2 + self.bnr**2)
          print("FFT completed")
        
     def get_ps(self):
          # get power spectrum
          pt = np.abs( self.cn )/ (self.N / 2.) 
          self.p = (pt[0:self.N/2]**2).clip(minval, maxval)
          self.ok = True
          print("Power spectrum completed")
          
     def cut_freq(self, a, value = 0.):
          # cut freq. between self.xmin and self.xmax and replace them with value
          # in array a
          factl = 0.5*(1. - np.tanh( self.par["alpha"] * (self.f - self.par["xmax"]) ))
          facth = 0.5*(np.tanh( self.par["alpha"] * (self.f - self.par["xmin"])) + 1 )
          # set the approp. parameters to 0
          ftot = (1.-factl*facth).clip(minval, maxval)
          aa = (a*ftot + value)
          return aa

     def low_pass_freq( self, a, value = 0.):
          # pass freq. up to self.xmax the steepness of the cut is controlled by alpha
          # index array of values to cut
          fact = 0.5*(1. - np.tanh( self.par["alpha"] * (self.f - self.par["xmax"]) )) 
          # set the approp. parameters to 0
          return a*(fact.clip(minval, maxval))  + value
     
     def high_pass_freq( self, a, value = 0.):
          # pass freq. from self.xmin,  the steepness of the cut is controlled by alpha
          # index array of values to cut
          fact = 0.5*(np.tanh( self.par["alpha"] * (self.f - self.par["xmin"])) + 1 )  
          # set the approp. parameters to 0
          return a*(fact.clip(minval, maxval)) + value
     
     def an_cut_freq(self, replace = True, **kwargs):
          # pass freq. outside of self.xmin, self.xmax, the steepness is controlled by alpha
          # cut for an
          if replace:
               self.an = self.cut_freq(self.an, **kwargs)
          else:
               return self.cut_freq(self.an, **kwargs)

     def bn_cut_freq(self, replace= True, ** kwargs):
          # pass freq. outside of self.xmin, self.xmax, the steepness is controlled by alpha
          # cut for bn
          # reverse b tp apply smae cut as to a
          bnr = self.cut_freq(self.bnr, **kwargs)
          if replace:
               self.bnr = bnr
               self.bn = bnr[::-1]
          else:
               return bnr[::-1]

     def r_cut_freq(self, replace = True):
          # pass all data there rr is within rcut, replace all data outside with a
          # the magnitude of the r-value and the appropriate phase
          ic = np.where(self.rr >= self.par["rthresh"])[0]
          if replace:
               zz = self.anc[ic] + 1.j * self.bnrc[ic]
               zz = set_z_value(zz, self.par["rthresh"])
               self.an[ic] = zz.real
               self.bnr[ic] = zz.imag
               self.bn = self.bnr[::-1]
          else:
               anc = self.an.copy()
               bnrc = self.bnr.copy()
               zz = self.anc[ic] + 1.j * self.bnrc[ic]
               zz = set_z_value(zz, self.par["rthresh"])
               anc[ic] = zz.real
               bnrc[ic] = zz.imag
               return anc, bnrc[::-1]

     def store_fft_coeff(self, an , bn, const = (0. + 0j), replace = True):
          # return a new complex fft coeff. array from
          # the an and bn arrays
          N = 2.* len(an)
          bnr = np.hstack( (bn[0], bn[1:][::-1]) ) # reverse array
          anr = np.hstack( (an[:1],an[1:][::-1]) ) # reverse array
          # construct new fft coeff
          z1 = (bnr - an*1j)*N/2.
          z2 = (bn + anr*1j)*N/2.
          z1[0] = const
          cn = np.hstack( (z1, z2) )
          if replace:
               self.cn = cn
          else:
               return cn

     def inv_transform(self):
          # perform the inverse FFT
          return FFTP.ifft(self.cn)
     
