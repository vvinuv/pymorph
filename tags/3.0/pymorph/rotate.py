import pyfits
import numpy as np 
import os

def ImSec(imdata, xcntr, ycntr, ex_rad):
    """Find the image section for asymmetry"""
    SizeY, SizeX = imdata.shape
    ExceedSize = 0
    #All are floor to make the size even number
    xmin = np.floor(xcntr - ex_rad)
    ymin = np.floor(ycntr - ex_rad)
    xmax = np.floor(xcntr + ex_rad)
    ymax = np.floor(ycntr + ex_rad)
    if xmin < 0:
        xmin = 0
        cut_xcntr = xcntr
        ExceedSize = 1
    else:
        cut_xcntr = ex_rad + np.modf(xcntr)[0]
    if ymin < 0:
        cut_ycntr = ycntr
        ymin = 0
        ExceedSize = 1
    else:
        cut_ycntr = ex_rad + np.modf(ycntr)[0]
    if xmax > SizeX - 1:
        xmax = SizeX
        ExceedSize = 1
    if ymax > SizeY - 1:
        ymax = SizeY
        ExceedSize = 1
    CutImDa = imdata[ymin:ymax, xmin:xmax].copy()
    SizeY, SizeX = CutImDa.shape
    return CutImDa, cut_xcntr, cut_ycntr, SizeY, SizeX, ymin, \
           ymax, xmin, xmax, ExceedSize

def rotate_slow(z, angle, xcntr, ycntr, cval=0.0):
    """Rotate an array. xcntr, ycntr are the iraf kind of values. Values \
       outside the image will be filled with cval. Deprecated"""
    xcntr -= 1.0 
    ycntr -= 1.0
    z = np.asarray(z)
    angle = np.pi / 180 * angle
    co = np.cos(angle)
    si = np.sin(angle)
    SizeY = z.shape[0] # Y size
    SizeX = z.shape[1] # X size
    x = np.reshape(np.arange(SizeX * SizeY), (SizeY, SizeX)) % SizeX
    x = x.astype(np.float32)
    y = np.reshape(np.arange(SizeX * SizeY), (SizeY, SizeX)) / SizeX
    y = y.astype(np.float32)
    rz = np.reshape(np.zeros(SizeX * SizeY), (SizeY, SizeX)) / SizeX
    rz = rz.astype(np.float32)
    rz += cval
    rx = (x - xcntr)* co + (y - ycntr) * si
    ry = (xcntr - x) * si + (y - ycntr) * co
    rx += xcntr
    ry += ycntr
    try:
     os.system('rm -f xnorm.fits xrotate.fits rotate.fits')
    except:
     pass
#    hdu = pyfits.PrimaryHDU(y-xcntr)
#    hdu.writeto('xnorm.fits')
#    hdu = pyfits.PrimaryHDU(ry)
#    hdu.writeto('xrotate.fits')
    # x.shape[0] is the y dimension
    for i in range(x.shape[1]): #runs through x. Also rx[y][x]
        for j in range(x.shape[0]): #runs through y
            x1 = np.floor(rx[j][i]) # jth position in Num. Rece. p. 124
            x2 = np.ceil(rx[j][i])  # j+1
            y1 = np.floor(ry[j][i]) # k
            y2 = np.ceil(ry[j][i])  # k+1
            if x1 >= 0 and x1 < SizeX and x2 >= 0 and x2 < SizeX and \
               y1 >= 0 and y1 < SizeY and y2 >= 0 and y2 < SizeY:
                z1 = z[y1][x1] # Eqn. 3.6.3. Note that here 
                               # x2 = x1 + 1, y2 =y2 +1
                z2 = z[y1][x2]
                z3 = z[y2][x2]
                z4 = z[y2][x1]
                t = rx[j][i] - x1 # Eqn. 3.6.4. The denominators are 
                                  # 1 in our case
                u = ry[j][i] - y1
                rz[j][i] = (1 - t) * (1 - u) * z1 + t * (1 - u) * z2 + \
                           t * u * z3 + (1 - t) * u * z4
    hdu = pyfits.PrimaryHDU(rz)
    hdu.writeto('rotate.fits')
    return rz

def rotate(z, angle, xcntr, ycntr, cval=0.0):
    """Rotate an array. xcntr, ycntr are the iraf kind of values. Values \
       outside the image will be filled with cval"""
    xcntr -= 1.0 
    ycntr -= 1.0
    z = np.asarray(z)
    angle = np.pi / 180 * angle
    co = np.cos(angle)
    si = np.sin(angle)
    SizeY = z.shape[0] # Y size
    SizeX = z.shape[1] # X size
    x = np.reshape(np.arange(SizeX * SizeY), (SizeY, SizeX)) % SizeX
    x = x.astype(np.float32)
    y = np.reshape(np.arange(SizeX * SizeY), (SizeY, SizeX)) / SizeX
    y = y.astype(np.float32)
    mx = np.reshape(np.zeros(SizeX * SizeY), (SizeY, SizeX)) 
    my = np.reshape(np.zeros(SizeX * SizeY), (SizeY, SizeX)) 
    rz = np.reshape(np.zeros(SizeX * SizeY), (SizeY, SizeX)) 
    rz = rz.astype(np.float32)
    rz += cval
    rx = (x - xcntr)* co + (y - ycntr) * si
    ry = (xcntr - x) * si + (y - ycntr) * co
    rx += xcntr
    ry += ycntr
#    try:
#     os.system('rm -f xnorm.fits xrotate.fits rotate.fits')
#    except:
#     pass
#    hdu = pyfits.PrimaryHDU(y-xcntr)
#    hdu.writeto('xnorm.fits')
#    hdu = pyfits.PrimaryHDU(ry)
#    hdu.writeto('xrotate.fits')
    # x.shape[0] is the y dimension
    x1 = np.floor(rx)
    y1 = np.floor(ry)
    x1 = x1.astype(int)
    y1 = y1.astype(int)
    con = (x1 < 0) | (x1 > SizeX-2)
    mx[con] = 1
    x1[con] = 0 
    con = (y1 < 0) | (y1 > SizeY-2)
    my[con] = 1 
    y1[con] = 0
    z1 = z[y1, x1]
    z2 = z[y1, x1 + 1]
    z3 = z[y1 + 1, x1 + 1]
    z4 = z[y1 + 1, x1]
    t = rx - x1
    u = ry - y1
    rz = (1 - t) * (1 - u) * z1 + t * (1 - u) * z2 + \
                         t * u * z3 + (1 - t) * u * z4
    rz[np.where(mx == 1)] = cval
    rz[np.where(my == 1)] = cval
#    hdu = pyfits.PrimaryHDU(rz)
#    hdu.writeto('rotate.fits')
    return rz

#f = pyfits.open('/home/vinu/test-pymorph/rewriting-pymorph-test/00018611_serexp.fits')
#z = f[0].data
#f.close()
#rotate1(z, 180, 95.95, 145.92, cval=100)
#rotate(z, 180, 62.19, 65.23, cval=0) #lower left
#rotate(z, 180, 157.49, 218.48, cval=0) #upper right

#CutImDa, cut_xcntr, cut_ycntr, SizeY, SizeX = ImSec(z, 157.49, 218.48, 100)
#try:
# os.system('rm -f z.fits')
#except:
# pass
#hdu = pyfits.PrimaryHDU(CutImDa)
#hdu.writeto('z.fits')
#rotate(CutImDa, 180, cut_xcntr, cut_ycntr, cval=0)
