#!/usr/bin/env python

import numpy as n
from pylab import *
import pyfits
import nmpfit

f=pyfits.open('O_I_EDCSNJ1216462-1200073.fits')
galaxy = f[3].data
f.close()
nn, bins, patches = hist(galaxy, 100, normed=0)


p1=500.0
p2=0.0
p3=0.01
p = [p1, p2, p3] 

nMaxArg = nn.argmax()

x=bins[nMaxArg-5:nMaxArg+5]
y=nn[nMaxArg-5:nMaxArg+5]
#x=bins
#y=nn
print x.shape, y.shape


def g(p,fjac = None, x = None, y=None, err=None,weights=None):
    if p[2] != 0.0: 
            Z = (x - p[1])
            model = p[0]*n.e ** (-Z**2 / (p[2] * p[2] * 2.0)) 
    else: 
            model = N.zeros(N.size(x))
    status = 0 
#    if n.sqrt(y)!=0.0:
    return [status, (y-model)]#/ (y)]
#    else:
#        return [status, (y-model)]

fa = {'x':x, 'y':y}
m=nmpfit.mpfit(g,p,functkw=fa)
p=m.params
print p
#figure()
print p[2] / abs(bins[1] - bins[0])
yy = p[0] * n.e**(-0.5*(bins-p[1])**2/p[2]**2)
l=plot(bins, yy, 'b--')
setp(l, 'linewidth', 3)
grid(1)
show()

