import sys
from numpy import NaN, Inf, arange, isscalar, asarray
import pdb

def peakdet(v, delta, x = None, gauge = None, power = 5):
    """

    MAXTAB,MINTAB = peakdet(v, delta, x)

    v is the array of values to analyze
    delta is the min. step around a maximum

    keyword arguments:
    x array of x-values corresponding to v
    gauge GaugeFrame  object for progress display
    power at which power of 10 the analysis progress should be printed

    MAXTAB[0] x-values for the maxima
    MAXTAB[1] v-values for the maxima
    MAXTAB[2] indices into v for the maxima

    MINTAB same for the minima
     
    """
    maxtab_x = []
    maxtab_y = []
    maxtab_i = []
    mintab_x = []
    mintab_y = []
    mintab_i = []

    # step size for progress information
    pstep = 10**power
    ndat = len(v)

    if x is None:
        x = arange(len(v))
    
    v = asarray(v)
    
    if len(v) != len(x):
        sys.exit('Input vectors v and x must have same length')
    
    if not isscalar(delta):
        sys.exit('Input argument delta must be a scalar')
    
    if delta <= 0:
        sys.exit('Input argument delta must be positive')
    
    mn, mx = Inf, -Inf
    mnpos, mxpos = NaN, NaN
    
    lookformax = True
    step = len(v)/10.
    for i in arange(len(v)):
        if i%pstep == 0:
            if gauge == None:
                print('process value number = ', i)
            else:
                gauge.Update(i, "%5.1f"%(float(i)/ndat*100.)+" % Completed" )
        this = v[i]
        if this > mx:
            mx = this
            mxpos = x[i]
            max_i = i
        if this < mn:
            mn = this
            mnpos = x[i]
            min_i = i
        
        if lookformax:
            if this < mx-delta:
                maxtab_x.append(mxpos)
                maxtab_y.append(mx)
                maxtab_i.append(max_i)
                mn = this
                min_i = i
                mnpos = x[i]
                lookformax = False
        else:
            if this > mn+delta:
                mintab_x.append(mnpos)
                mintab_y.append(mn)
                mintab_i.append(min_i)
                mx = this
                mxpos = x[i]
                max_i = i
                lookformax = True
    return [maxtab_x, maxtab_y, maxtab_i], [mintab_x, mintab_y, mintab_i]

if __name__=="__main__":
    series = [0,0,0,2,0,0,0,-2,0,0,0,2,0,0,0,-2,0]
    print(peakdet(series,1))
