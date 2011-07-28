import config as c
import os
import numpy as n
import pylab as p
import pyfits
import convolve as conv
import numpy.ma as ma

def QuarterMask(z, zm, xcntr, ycntr, bbya, pa, quarter):
    nxpts, nypts = z.shape
    zmm = n.reshape(n.ones(nxpts * nypts), (nxpts, nypts))
    co = n.cos(pa * n.pi / 180.0)
    si = n.sin(pa * n.pi / 180.0)
    one_minus_eg_sq = (bbya)**2.0
    x = n.reshape(n.arange(nxpts * nypts), (nxpts, nypts)) % nypts
    x = x.astype(n.float32)
    y = n.reshape(n.arange(nxpts * nypts), (nxpts, nypts)) / nypts
    y = y.astype(n.float32)
    xrot = (x - xcntr) * co + (y - ycntr) * si
    xsq = xrot**2.0
    yrot = (xcntr - x) * si + (y - ycntr) * co
    ysq = yrot**2.0 
    r = n.sqrt(xsq + ysq / one_minus_eg_sq)
    if quarter == 0: 
        condition = xrot > -1e5
    if quarter == 1:
        condition = (xrot - 0 >= 0) & (yrot - 0 >= 0)
    if quarter == 2:
        condition = (xrot - 0 < 0) & (yrot - 0 >= 0)
    if quarter == 3:
        condition = (xrot - 0 < 0) & (yrot - 0 < 0)
    if quarter == 4:
        condition = (xrot - 0 >= 0) & (yrot - 0 < 0)
    zmm[condition] = 0
    zmm = zm + zmm
    zmm[n.where(zmm > 0)] = 1
    return n.median(ma.masked_array(z, zmm).compressed())

def FindYetSky(gimg, X0, Y0):
    try:
        SEx_GAIN = c.SEx_GAIN 
    except:
        SEx_GAIN = 1.0
    SEx_PIXEL_SCALE = c.SEx_PIXEL_SCALE
    SEx_FILTER_NAME = c.SEx_FILTER_NAME
    mag_zero = c.mag_zero 
    pymorph_path = c.PYMORPH_PATH
    f = pyfits.open(gimg)
    z = f[0].data
    f.close()
    f_tpl = open(str(c.PYMORPH_PATH) + '/SEx/default_seg.sex','r')
    template = f_tpl.read()
    f_tpl.close()
    f_sex = open('default_seg.sex', 'w')
    f_sex.write(template %vars())
    f_sex.close()
    cmd = c.SEX_PATH + ' ' + gimg + ' -c default_seg.sex > /dev/null'
    os.system(cmd)
    f = pyfits.open('seg.fits')
    zm = f[0].data
    f.close()
    SexSky, SkyYet, SkyMed, SkyMin, SkyQua, SkySig = 9999, 9999, 9999, \
    9999, 9999, 9999
    for l_s in open('SegCat.cat'):
        v_s = l_s.split()
        try:
            SexId = n.float(v_s[0])
            xcntr = n.float(v_s[1])
            ycntr = n.float(v_s[2])
            pa = n.float(v_s[11])
            bbya = 1.0 / n.float(v_s[12])
            a = n.float(v_s[14]) * 3.0
            b = a * bbya
            hr = n.float(v_s[9])
            sky = n.float(v_s[10])
            if n.abs(X0 - xcntr) < 5.0 and n.abs(Y0 - ycntr) < 5.0:
               zm = conv.boxcar(zm, (3, 3), mode='nearest')
               zm[n.where(zm > 0)] = 1
               SkyQua = []
               for ii in n.arange(1, 5): 
                   SkyQua.append(QuarterMask(z, zm, xcntr-1.0, ycntr-1.0, \
                                 bbya, pa, ii))
               SkyQua = n.array(SkyQua)
               SexSky = sky
               SkyYet = n.median(ma.masked_array(z, zm).compressed())
               SkyMed = n.median(SkyQua)
               SkyMin = n.min(SkyQua)
               SkySig = n.std(ma.masked_array(z, zm).compressed())
#               os.system('rm -f SegCat.cat default_seg.sex seg.fits') 
               break
        except:
            pass 
    return SexSky, SkyYet, SkyMed, SkyMin, SkyQua, SkySig
