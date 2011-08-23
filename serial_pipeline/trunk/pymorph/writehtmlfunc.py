from os.path import exists
import csv
import sys
import config as c
import numpy as n
import fileinput
from cosmocal import *
import datetime
import MySQLdb as mysql
#from utilities import WriteDb
import traceback

class WriteHtmlFunc:
    """The class which will write html and csv output. This class will also 
       check whether the fit is good or bad using the Chisq and Goodness value
       It will also notify the goodness/badness of fit"""
    def __init__(self, cutimage, xcntr, ycntr, distance, alpha_j, delta_j, z, Goodness, C, C_err, A, A_err, S, S_err, G, M, EXPTIME):
        self.cutimage     = cutimage
        self.xcntr        = xcntr
        self.ycntr        = ycntr
        self.distance     = distance
        self.alpha_j      = alpha_j
        self.delta_j      = delta_j
        self.z            = z
        self.Goodness     = Goodness
        self.C            = C
        self.C_err        = C_err
        self.A            = A
        self.A_err        = A_err
        self.S            = S
        self.S_err        = S_err
        self.G            = G
        self.M            = M
        self.EXPTIME      = EXPTIME
        self.write_params = write_params(cutimage, xcntr, ycntr, distance, alpha_j, delta_j, z, Goodness, C, C_err, A, A_err, S, S_err, G, M, EXPTIME)

def write_params(cutimage, xcntr, ycntr, distance, alpha_j, delta_j, z, Goodness, C, C_err, A, A_err, S, S_err, G, M, EXPTIME):

    lanczosG = 7
    lanczos_coef = [0.99999999999980993, 676.5203681218851,\
                    -1259.1392167224028, 771.32342877765313,\
                    -176.61502916214059, 12.507343278686905, \
                    -0.13857109526572012, 9.9843695780195716e-6,\
                    1.5056327351493116e-7]
    def Gamma(z):
        """This is the Lanczos approximation for Gamma function"""
        if z < 0.5:
            return n.pi / (n.sin(n.pi*z)*Gamma(1-z))
        else:
            z -= 1
            x = lanczos_coef[0]
            for i in range(1, lanczosG + 2):
                x += lanczos_coef[i]/(z + i)
            t = z + lanczosG + 0.5
            return n.sqrt(2*n.pi) * t**(z + 0.5) * n.exp(-t) * x
    def bfunc(x):
        """ This function gives value of b_n given the Sersic index"""
        return 0.868242*x -0.142058 # Khosroshahi et al. 2000 approximation
    try:
        ComP = c.components
    except:
        ComP = ['bulge', 'disk']
    if len(ComP) == 0:
        ComP = ['bulge', 'disk']
    Goodness = float(str(round(Goodness, 3))[:5])
    f_tpl = open(str(c.PYMORPH_PATH) + '/html/default.html', 'r')
    template = f_tpl.read()
    f_tpl.close()
    outfile = open('R_' + str(cutimage)[:-4] + 'html','w')
    ra1 = int(float(alpha_j) / 15.0)
    ra2 = int((float(alpha_j) / 15.0 - int(float(alpha_j) / 15.0))*60.0)
    ra3 = (((float(alpha_j) / 15.0 - int(float(alpha_j) / 15.0))*60.0) - ra2) * 60.0
    dec1 = int(float(delta_j))
    dec2 = abs(int((float(delta_j) - dec1) * 60.0))
    dec3 = (abs(((float(delta_j) - dec1) * 60.0)) - dec2) * 60.0

    if ra1 < 10:
        alpha1 = '0' + str(ra1)
    else:
        alpha1 = str(ra1)
    if ra2 < 10:
        alpha2 = '0' + str(ra2)
    else:
        alpha2 = str(ra2)
    if ra3 < 10:
        alpha3 = '0' + str(ra3)[:3]
    else:
        alpha3 = str(ra3)[:4]
    if abs(dec1) < 10:
        delta1 = '0' + str(dec1)
    else:
        delta1 = str(dec1)
    if dec1 > 0:
        delta1 = '+' + str(delta1)
    else:
        pass
    if dec2 < 10:
        delta2 = '0' + str(dec2)
    else:
        delta2 = dec2
    if dec3 < 10:
        delta3 = '0' + str(dec3)[:3]
    else:
        delta3 = str(dec3)[:4]
    if(c.repeat == False or c.repeat):
        for line_i in fileinput.input("index.html",inplace =1):
            line_i = line_i.strip()
            if not '</BODY></HTML>' in line_i:
                print line_i    
        indexfile = open('index.html', 'a+')
        NoImage = 1
        for indexline in indexfile:
            if cutimage[:-5] in indexline:
                NoImage = 0
            else:
                pass
        if NoImage:
            indexfile.writelines(['<a href="R_',\
                                  str(cutimage)[:-5],'.html',\
                                  '"> ', str(cutimage)[:-5],\
                                  ' </a> <br>\n'])
        indexfile.writelines(['</BODY></HTML>\n'])
        indexfile.close()
    object = 1
    object_err = 1
    if exists('fit.log'):
        for line in open('fit.log','r'): 
            values = line.split() 
            try: 
                if(str(values[0]) == 'Input'):
                    alpha_ned = str(alpha_j)[:10]
                    delta_ned = str(delta_j)[:10]
                if(str(values[0]) == 'Init.'):
                    initial_conf = str(values[4])
                if(str(values[0]) == 'Restart'):
                    restart_conf = str(values[3])
                if(str(values[0]) == 'Chi^2/nu'):
                    chi2nu = float(values[2])
                    Distance = str(round(distance, 3))[:5]
                #if(str(values[0]) == 'sersic' and object == 2):
                #    mag_bar = float(values[4])
                #    re_bar = float(values[5])
                #    SersicIndexBar = float(values[6])
                #    SersicEllipticityBar = float(values[7])
		#    SersicBoxyBar = float(values[9])
                #    object += 1
                if(str(values[0]) == 'sersic' and object == 1):
                    mag_b = float(values[4])
                    re = float(values[5])
                    SersicIndex = float(values[6])
                    SersicEllipticity = float(values[7])
		    SersicBoxy = float(values[9])
                    object += 1
                if(str(values[0]) == 'expdisk'):
                    mag_d = float(values[4])
                    rd = float(values[5])
                    DiskEllipticity = float(values[6])
		    DiskBoxy = float(values[8])
                if(str(values[0]) == 'psf'):
                    mag_p = float(values[4])
                if(str(values[0]) == 'sky'):
                    galfit_sky = float(values[4])
                if(str(values[0])[:1] == '('):
                 #   if(str(a) == 'sersic' and object_err == 2):
                 #       mag_bar_err = float(values[2])
                 #       re_bar_err  = float(values[3])
                 #       SersicIndexBarErr = float(values[4])
                 #       SersicEllipticityBarErr = float(values[5])
		 #       SersicBoxyBarErr = float(values[7])
                 #       object_err += 1
                    if(str(a) == 'sersic' and object_err == 1):
                        mag_b_err = float(values[2])
                        re_err  = float(values[3])
                        SersicIndexErr = float(values[4])
                        SersicEllipticityErr = float(values[5])
			SersicBoxyErr = float(values[7])
                        object_err += 1
                    if(str(a) == 'expdisk'):
                        mag_d_err = float(values[2])
                        rd_err = float(values[3])
                        DiskEllipticityErr = float(values[4])
			DiskBoxyErr = float(values[6])
                    if(str(a) == 'psf'):
                        mag_p_err = float(values[2])

                a=values[0]				
            except:
                pass
        if(z == 9999):
            re_kpc = 9999
            re_err_kpc = 9999
            rd_kpc = 9999
            rd_err_kpc = 9999
            DisMoD = 9999
            re_bar_kpc = 9999
            re_bar_err_kpc = 9999
        else:
            phy_parms = cal(z, c.H0, c.WM, c.WV, c.pixelscale)
            DisMoD = phy_parms[2]
            if 'bulge' in ComP:
                re_kpc = phy_parms[3] * re
                re_err_kpc = phy_parms[3] * re_err
            else:
                re_kpc = 9999
                re_err_kpc = 9999
            if 'disk' in ComP:
                rd_kpc = phy_parms[3] * rd
                rd_err_kpc = phy_parms[3] * rd_err
            else:
                rd_kpc = 9999
                rd_err_kpc = 9999
            if 'bar' in ComP:
                re_bar_kpc = phy_parms[3] * re_bar
                re_bar_err_kpc = phy_parms[3] * re_bar_err
            else:
                re_bar_kpc = 9999
                re_bar_err_kpc = 9999

#            if 'point' in ComP:
#                fwhm_kpc = 0.5 * phy_parms[3]
#            else:
#                fwhm_kpc = 9999
        if 'bulge' in ComP and 'disk' in ComP:
            if 'point' in ComP:
                fb = 10**(-0.4 * (mag_b - c.mag_zero))
                fd = 10**(-0.4 * (mag_d - c.mag_zero))
                fp = 10**(-0.4 * (mag_p - c.mag_zero))
                BD = fb / fd 
                BT = fb / (fb + fd + fp)
            elif 'bar' in ComP:
                fb = 10**(-0.4 * (mag_b - c.mag_zero))
                fd = 10**(-0.4 * (mag_d - c.mag_zero))
                fbar = 10**(-0.4 * (mag_bar - c.mag_zero)) 
                BD = fb / fd
                BT = fb / (fb + fd + fbar)
            else:
                BD = 10**(-0.4 * ( mag_b - mag_d))
                BT = 1.0 / (1.0 + 1.0 / BD)
        elif 'bulge' in ComP:
            BD = 'nan'
            BT = 1.0
        elif 'disk' in ComP:
            BD = 0.0
            BT = 0.0
        else:
            BD = 'nan'
            BT = 'nan'
    pngfile = 'P_' + str(cutimage)[:-4] + 'png'
    object = 1
    object_err = 1
    Neighbour_Sersic = ''
    Object_Sersic = ''
    Object_Sersic_err = ''
    Object_Exp = ''
    Object_Exp_err = ''
    Point_Vals = ''
    Point_Vals_err = ''
    if exists('fit.log'):
        for line in open('fit.log','r'): 
            values = line.split() 
            try: 
                if(str(values[0]) == 'sersic'):
                    if object > 1 or 'bulge' not in ComP:
                        Neighbour_Sersic = str(Neighbour_Sersic) + \
                                          '<TR align="center" bgcolor=' + \
                                          '"#99CCFF"><TD>' + str(values[0]) + \
                                          ' </TD> <TD> ' + \
                                          str(values[2])[1:-1] + '</TD> <TD> '\
                                          + str(values[3])[:-1] + \
                                          ' </TD> <TD> ' + str(values[4]) + \
                                          ' </TD> <TD> ' + \
                                          str(values[5]) + ' </TD> <TD> ' + \
                                          ' ' + ' </TD> <TD> ' + \
                                          str(values[6]) + ' </TD> <TD> ' +\
                                          str(values[7]) + ' </TD> <TD> ' +\
                                          str(values[8]) + ' </TD> <TD> ' + \
                                          str(values[9]) + ' </TD> </TR>'
                    if object == 1 and 'bulge' in ComP:
                        bulge_xcntr = float(str(values[2])[1:-1])
                        bulge_ycntr = float(str(values[3])[:-1])
                        Object_Sersic = '<TR align="center" ' +\
                                        'bgcolor="#99CCFF">' +\
                                       '<TD>' + str(values[0]) + '</TD> <TD> '\
                                       + str(values[2])[1:-1] + ' </TD> <TD> '\
                                       + str(values[3])[:-1] + ' </TD> <TD> ' \
                                       + str(values[4]) + ' </TD> <TD> ' \
                                       + str(values[5]) +  ' </TD> <TD> ' + \
                                       str(round(re_kpc, 3))[:5] +'</TD> <TD>'\
                                       + str(values[6]) + ' </TD> <TD> ' \
                                       + str(values[7]) + ' </TD> <TD> ' \
                                       + str(values[8]) + ' </TD> <TD> ' \
                                       + str(values[9]) + ' </TD></TR>'
                        object += 1
                if(str(values[0]) == 'expdisk'):
                    disk_xcntr = float(str(values[2])[1:-1])
                    disk_ycntr = float(str(values[3])[:-1])
                    Object_Exp = '<TR align="center" bgcolor="#99CCFF">' +\
                                '<TD>' + str(values[0]) + ' </TD> <TD> ' + \
                                str(values[2])[1:-1] + ' </TD> <TD> ' + \
                                str(values[3])[:-1] + ' </TD> <TD> ' + \
                                str(values[4]) + ' </TD> <TD> ' + \
                                str(values[5]) +  ' </TD> <TD> ' + \
                                str(round(rd_kpc, 3))[:5] + ' </TD> <TD> ' + \
                                ' ' + ' </TD> <TD> ' + str(values[6]) +  \
                                ' </TD> <TD> ' + str(values[7]) + \
                                ' </TD> <TD> ' + str(values[8]) + ' </TD></TR>'
                if(str(values[0]) == 'psf'):
                    point_xcntr = float(str(values[2])[1:-1])
                    point_ycntr = float(str(values[3])[:-1])
                    Point_Vals = '<TR align="center" bgcolor="#99CCFF">' + \
                                 '<TD>' + str(values[0]) + ' </TD> <TD> ' + \
                               str(values[2])[1:-1] + ' </TD> <TD> ' + \
                               str(values[3])[:-1] + ' </TD> <TD> ' + \
                               str(values[4]) + ' </TD> <TD> ' + \
                               str('9999') +  ' </TD> <TD> ' + \
                               str('9999') + ' </TD> <TD> ' + \
                               ' ' + ' </TD> <TD> ' + str('9999') +  \
                               ' </TD> <TD> ' + str('9999') + \
                               ' </TD> <TD> ' + str('9999') + ' </TD></TR>'
                if(str(values[0])[:1] == '('):
                    if(str(a) == 'sersic' and object_err > 1 or \
                        str(a) == 'sersic'  and 'bulge' not in ComP):
                        Neighbour_Sersic = str(Neighbour_Sersic) + \
                                           '<TR align="center" ' + \
                                           'bgcolor="#CCFFFF"> <TD>' + ' ' + \
                                           ' </TD> <TD> ' + \
                                           str(values[0])[1:-1] + \
                                           ' </TD> <TD> ' + \
                                           str(values[1])[:-1] + \
                                           ' </TD> <TD> ' + str(values[2]) + \
                                           ' </TD> <TD> ' + str(values[3]) + \
                                           ' </TD> <TD> ' + ' ' + \
                                           ' </TD> <TD> ' + str(values[4]) + \
                                           ' </TD> <TD> ' + str(values[5]) + \
                                           ' </TD> <TD> ' + str(values[6]) + \
                                           ' </TD> <TD> ' + str(values[7]) + \
                                           ' </TD> </TR> '
                    if(str(a) == 'sersic' and object_err == 1 and \
                        'bulge' in ComP ):
                        Object_Sersic_err = '<TR align="center" ' + \
                                           'bgcolor="#CCFFFF">' + \
                                           '<TD>' + ' ' + '</TD> <TD>' + \
                                          str(values[0])[1:-1] + '</TD> <TD>'\
                                           + str(values[1])[:-1] + \
                                           ' </TD> <TD> ' + str(values[2]) + \
                                           ' </TD> <TD> ' + str(values[3]) + \
                                           ' </TD> <TD> ' + \
                                           str(round(re_err_kpc, 3))[:5] + \
                                           ' </TD> <TD> ' + str(values[4]) + \
                                           ' </TD> <TD> ' + str(values[5]) + \
                                           ' </TD> <TD> ' + str(values[6]) + \
                                           ' </TD> <TD> ' + str(values[7]) + \
                                           ' </TD></TR>'
                        object_err += 1
                    if(str(a) == 'expdisk'):
                        Object_Exp_err = '<TR align="center" ' + \
                                         'bgcolor="#CCFFFF">' + \
                                         '<TD>' + ' ' + '</TD> <TD>' + \
                                         str(values[0])[1:-1] + '</TD> <TD>'\
                                         + str(values[1])[:-1] + \
                                         ' </TD> <TD> ' + str(values[2]) + \
                                         ' </TD> <TD> ' + str(values[3]) + \
                                         ' </TD> <TD> ' + \
                                         str(round(rd_err_kpc, 3))[:5] + \
                                         ' </TD> <TD> ' + ' ' + \
                                         ' </TD> <TD> ' + str(values[4]) + \
                                         ' </TD> <TD> ' + str(values[5]) + \
                                         ' </TD> <TD> ' + str(values[6]) + \
                                         ' </TD></TR>'
                    if(str(a) == 'psf'):
                        Point_Vals_err = '<TR align="center" ' + \
                                         'bgcolor="#CCFFFF">' + \
                                         '<TD>' + ' ' + '</TD> <TD>' + \
                                         str(values[0])[1:-1] + '</TD> <TD>'\
                                         + str(values[1])[:-1] + \
                                         ' </TD> <TD> ' + str(values[2]) + \
                                         ' </TD> <TD> ' + str('9999') + \
                                         ' </TD> <TD> ' + \
                                         str('9999') + \
                                         ' </TD> <TD> ' + ' ' + \
                                         ' </TD> <TD> ' + str('9999') + \
                                         ' </TD> <TD> ' + str('9999') + \
                                         ' </TD> <TD> ' + str('9999') + \
                                         ' </TD></TR>'
                a=values[0]				
            except:
                pass
    if 'bulge' in ComP:
        if EXPTIME == -9999:
            EXPTIME = 1
        else:
            pass
#        Ftot = EXPTIME * 10 ** ((mag_b - c.magzero) / -2.5)
#        Ie_Avg = Ftot / (2.0 * 3.14 * re * re)
#        MagInsideRe = -2.5  * n.log10(Ie_Avg / EXPTIME) + c.magzero
        try:
              pixelscale = c.pixelscale
        except:
              pixelscale = 1
        AvgMagInsideRe = mag_b + 2.5 * n.log10(2 * 3.14 * pixelscale * \
                      pixelscale * re * re * n.sqrt(1 - SersicEllipticity**2.0))
#        AvgMagInsideReErr = n.sqrt(mag_b_err**2.0 + (2.17 * re_err / re)**2.0)
        AvgMagInsideReErr2 = (1.085 * n.sqrt((2 * re * re_err)**2.0 + \
                            ((SersicEllipticity * SersicEllipticityErr) / \
                             n.sqrt(1 - SersicEllipticity**2.0))**2.0)) / \
                             (n.sqrt(1 - SersicEllipticity**2.0) * 2 * 3.14 * \
                              re * re)
        AvgMagInsideReErr = n.sqrt(mag_b_err**2.0 + AvgMagInsideReErr2**2.0)
    else:
        AvgMagInsideRe = 9999
        AvgMagInsideReErr = 9999
    wC = str(C)[:5]
    wA = str(A)[:5]
    wS = str(S)[:5]
    wG = str(G)[:5]
    wM = str(M)[:5]
    wBD = str(BD)[:5]
    wBT = str(BT)[:5]
    wAvgMagInsideRe = str(AvgMagInsideRe)[:5]
    error_mesg1 = ''
    error_mesg2 = ''
    error_mesg3 = ''
    error_mesg4 = ''
    error_mesg5 = ''
    error_mesg6 = ''
    if c.starthandle:
        error_mesg6 = '<a href="R_' + str(cutimage)[:-5] + \
	              '_1.html"> Crashed </a>' 
    HitLimit = 1
    if 'bulge' in ComP:
        if abs(mag_b - c.UMag) < 0.2 or abs(mag_b - c.LMag) < 0.2 or \
           abs(re - c.URe) < 1.0 or abs(re - c.LRe) < 0.1 or\
           abs(SersicIndex - c.LN) < 0.03 or abs(SersicIndex - c.UN) < 0.5:
            c.Flag += 65536
            HitLimit = 0
        else:
            pass
    else:
        bulge_xcntr = xcntr
        bulge_ycntr = ycntr
    if 'disk' in ComP:
        if abs(mag_d - c.UMag) < 0.2 or abs(mag_d - c.LMag) < 0.2 or \
            abs(rd - c.LRd) < 0.1 or abs(rd - c.URd) < 1.0:
            c.Flag += 131072
            HitLimit = 0
        else:
            pass
    else:
        disk_xcntr = xcntr
        disk_ycntr = ycntr
    if chi2nu <= c.chi2sq and Goodness >= c.Goodness and \
        HitLimit and \
        abs(bulge_xcntr - xcntr) <= c.center_deviation and \
        abs(bulge_ycntr - ycntr) <= c.center_deviation and \
        abs(disk_xcntr - xcntr) <= c.center_deviation and \
        abs(disk_ycntr - ycntr) <= c.center_deviation:
        img_notify = str(c.PYMORPH_PATH) + '/html/goodfit.gif'
        good_fit = 1
    else:
        if chi2nu > c.chi2sq:
            error_mesg1 = str(error_mesg1) + 'Chi2nu is large!'
            if chi2nu != 9999:
                c.Flag += 8192
        if Goodness < c.Goodness:
            error_mesg2 = str(error_mesg2) + 'Goodness is poor!'
            c.Flag += 16384
        if HitLimit == 0:
            error_mesg4 = str(error_mesg4) + 'One of the parameters'
            error_mesg5 = str(error_mesg5) + '          hits limit!'
        if abs(bulge_xcntr - xcntr) > c.center_deviation or \
             abs(bulge_ycntr - ycntr) > c.center_deviation or \
             abs(disk_xcntr - xcntr) > c.center_deviation or \
             abs(disk_ycntr - ycntr) > c.center_deviation:
            error_mesg3 = str(error_mesg3) + 'Fake Center!'
            if bulge_xcntr == 9999 or bulge_ycntr == 9999 or \
               disk_xcntr == 9999 or disk_ycntr == 9999:
                pass
            else:
                c.Flag += 32768
        img_notify = str(c.PYMORPH_PATH) + '/html/badfit.gif'
        good_fit = 0
    outfile.write(template %vars())
    outfile.close()
    run = 1
    to_remove = len(c.rootname) + 2
    if exists('result.csv'):
        for line_res in csv.reader(open('result.csv').readlines()[1:]):    
            if(str(line_res[0]) == cutimage[to_remove:-5]):
                run += 1
    f_res = open("result.csv", "ab")
    writer = csv.writer(f_res)
    galid = str(cutimage)[to_remove:-5]
    ParamToWrite = [galid, alpha_j, delta_j, z, c.SexMagAuto, c.SexMagAutoErr]
    if 'bulge' in ComP:
        for bulgecomp in [mag_b, mag_b_err, re, re_err, re_kpc, re_err_kpc, \
                          SersicIndex, SersicIndexErr, AvgMagInsideRe,\
                          AvgMagInsideReErr, SersicEllipticity, \
                          SersicEllipticityErr, SersicBoxy, \
			  SersicBoxyErr]:
            ParamToWrite.append(bulgecomp)
    else:
        for bulgecomp in [9999, 9999, 9999, 9999, 9999, 9999, 9999, 9999, \
                          9999, 9999, 9999, 9999, 9999, 9999]:
            ParamToWrite.append(bulgecomp)
    if 'disk' in ComP:
        for diskcomp in [mag_d, mag_d_err, rd, rd_err, rd_kpc, rd_err_kpc,\
                         DiskEllipticity, DiskEllipticityErr, \
			 DiskBoxy, DiskBoxyErr]:
            ParamToWrite.append(diskcomp)
    else:
        for diskcomp in [9999, 9999, 9999, 9999, 9999, 9999, 9999, 9999, \
		        9999, 9999]:
            ParamToWrite.append(diskcomp)
    ParamToWrite.append(BD)
    ParamToWrite.append(BT)
    if 'point' in ComP:
        ParamToWrite.append(mag_p)
        ParamToWrite.append(mag_p_err)
        ParamToWrite.append(0.5)
        ParamToWrite.append(9999)
    else:
        ParamToWrite.append(9999)
        ParamToWrite.append(9999)
        ParamToWrite.append(9999)
        ParamToWrite.append(9999)
    try:
        galfit_sky * 1.0
    except:
        galfit_sky = 9999
        print 'GALFIT does not report sky'
    if c.GalSky != 9999:
        galfit_sky = c.GalSky
    for otherparam in [chi2nu, Goodness, run, C, C_err, A, A_err, S, S_err, G,\
                       M, c.SexSky, galfit_sky, DisMoD, \
                       distance, good_fit, c.Flag, c.SexHalfRad]:
        ParamToWrite.append(otherparam)
    if 'bar' in ComP:
        for barcomp in [mag_bar, mag_bar_err, re_bar, re_bar_err, re_bar_kpc, \
                          re_bar_err_kpc, \
                          SersicIndexBar, SersicIndexBarErr, \
                          SersicEllipticityBar, \
                          SersicEllipticityBarErr, SersicBoxyBar]:
            ParamToWrite.append(barcomp)
    else:
        for barcomp in [9999, 9999, 9999, 9999, 9999, 9999, 9999, 9999, \
                          9999, 9999, 9999]:
            ParamToWrite.append(barcomp)

    writer.writerow(ParamToWrite)
    f_res.close()
    #Remove any nan or inf from the parameter
    for p in ParamToWrite:
        if str(p) in ('nan', '-nan', 'inf', '-inf'):
            ParamToWrite[ParamToWrite.index(p)] = 9999
        else:
            pass
    try:
        from utilities import WriteDb
    except:
        print 'DB writing Problem'
    try:
        WriteDb(ParamToWrite)
    except:
        print 'No database can be created!'
        traceback.print_exc()
    outfile1 = open('pymorph.html', 'w')
    outfile1.write('<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01' \
                   ' Transitional//EN">')
    outfile1.write('<HTML><HEAD> \n')
    outfile1.writelines(['<META HTTP-EQUIV="Refresh" CONTENT="10; \
                    URL=pymorph.html"> \n'])
    outfile1.write('</HEAD> \n')
    outfile = 'R_' + str(cutimage)[:-4] + 'html'
    try:
        for line in open(outfile, 'r'):
            if(str(line) != '<html>'):
                outfile1.writelines([str(line)])
    except:
        pass
    outfile1.close()
#write_params('I_EDCSNJ1216453-1201176.fits', 60.0, 60.0, 47.86, 188.6875, -12.7763, 0.67, 0.9887, 9999, 9999, 9999, 9999, 9999, 9999, 9999, 9999, 1024)

