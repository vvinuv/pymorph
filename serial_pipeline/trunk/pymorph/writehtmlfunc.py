from os.path import exists
import csv
import sys
import numpy as n
import fileinput
from cosmocal import cal 
import datetime
import MySQLdb as mysql
import traceback
from flagfunc import GetFlag, isset, Get_FitFlag
from pymorphutils import RaDegToHMS, DecDegToDMS
import config as c
import os

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
        self.WriteParams = WriteParams(cutimage, xcntr, ycntr, distance, alpha_j, delta_j, z, Goodness, C, C_err, A, A_err, S, S_err, G, M, EXPTIME)

def WriteParams(cutimage, xcntr, ycntr, distance, alpha_j, delta_j, z, Goodness, C, C_err, A, A_err, S, S_err, G, M, EXPTIME):
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
    ra1, ra2, ra3 = RaDegToHMS(alpha_j)
    dec1, dec2, dec3 = DecDegToDMS(delta_j)
    # Formatting ra and dec for display purpose
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
    # Writing index file
    if(c.repeat == False or c.repeat):
        for line_i in fileinput.input("index.html",inplace =1):
            line_i = line_i.strip()
            if not '</BODY></HTML>' in line_i:
                print line_i    
        indexfile = open('index.html', 'a+')
        NoImage = 1
        for indexline in indexfile:
            if c.fstring in indexline:
                NoImage = 0
            else:
                pass
        if NoImage:
            indexfile.writelines(['<a href="R_',\
                                  c.fstring,'.html',\
                                  '"> ', c.fstring,\
                                  ' </a> <br>\n'])
        indexfile.writelines(['</BODY></HTML>\n'])
        indexfile.close()
    # Reading fit.log
    if 'bar' in ComP:
        basic_info, fit_info = read_fitlog(filename = 'fit.log', yes_bar = 1)
    else:
        basic_info, fit_info = read_fitlog(filename = 'fit.log', yes_bar = 0)


    if 'Input' in basic_info:
        alpha_ned = str(alpha_j)[:10]
        delta_ned = str(delta_j)[:10]
        
    initial_conf = basic_info['initial_conf']
    restart_conf = basic_info['restart_conf']

    # move the restart file to a reasonably named output file
    new_outname = initial_conf.replace('in','out')
    os.rename(restart_conf, new_outname)
    basic_info['restart_conf'] = new_outname

    
    chi2nu = basic_info['chi2nu']
    Distance = str(round(distance, 3))[:5]
    if 'bulge' in fit_info:
        bulge_xcntr = fit_info['bulge']['xctr'][0]
        bulge_ycntr = fit_info['bulge']['yctr'][0]
        mag_b = fit_info['bulge']['mag'][0]
        re = fit_info['bulge']['rad'][0]
        SersicIndex = fit_info['bulge']['n'][0]
        SersicEllipticity = fit_info['bulge']['ell'][0]
        SersicPA = fit_info['bulge']['pa'][0]
        SersicBoxy = fit_info['bulge']['boxy'][0]

        bulge_xcntr_err = fit_info['bulge']['xctr'][1]
        bulge_ycntr_err = fit_info['bulge']['yctr'][1]
        mag_b_err = fit_info['bulge']['mag'][1]
        re_err = fit_info['bulge']['rad'][1]
        SersicIndexErr = fit_info['bulge']['n'][1]
        SersicEllipticityErr = fit_info['bulge']['ell'][1]
        SersicPAErr = fit_info['bulge']['pa'][1]
        SersicBoxyErr = fit_info['bulge']['boxy'][1]
    else:
        bulge_xcntr = -999.
        bulge_ycntr = -999.
        mag_b = -999.
        re = -999.
        SersicIndex = -999.
        SersicEllipticity = -999.
        SersicPA = -999.
        SersicBoxy = -999.

        bulge_xcntr_err = -999.
        bulge_ycntr_err = -999.
        mag_b_err = -999.
        re_err = -999.
        SersicIndexErr = -999. 
        SersicEllipticityErr = -999.
        SersicPAErr = -999.
        SersicBoxyErr = -999.

    if 'disk' in fit_info:
        disk_xcntr = fit_info['disk']['xctr'][0]
        disk_ycntr = fit_info['disk']['yctr'][0]
        mag_d = fit_info['disk']['mag'][0]
        rd = fit_info['disk']['rad'][0]
        DiskEllipticity = fit_info['disk']['ell'][0]
        DiskPA = fit_info['disk']['pa'][0]
        DiskBoxy = fit_info['disk']['boxy'][0]

        disk_xcntr_err = fit_info['disk']['xctr'][1]
        disk_ycntr_err = fit_info['disk']['yctr'][1]
        mag_d_err = fit_info['disk']['mag'][1]
        rd_err = fit_info['disk']['rad'][1]
        DiskEllipticityErr = fit_info['disk']['ell'][1]
        DiskPAErr = fit_info['disk']['pa'][1]
        DiskBoxyErr = fit_info['disk']['boxy'][1]
    else:
        disk_xcntr = -999.
        disk_ycntr = -999.
        mag_d = -999.
        rd = -999.
        DiskEllipticity = -999.
        DiskPA = -999.
        DiskBoxy = -999.

        disk_xcntr_err = -999.
        disk_ycntr_err = -999.
        mag_d_err = -999.
        rd_err = -999.
        DiskEllipticityErr = -999.
        DiskPAErr = -999.
        DiskBoxyErr = -999.

    if 'psf' in fit_info:
        psf_xcntr = fit_info['disk']['xctr'][0]
        psf_ycntr = fit_info['disk']['yctr'][0]
        mag_p = fit_info['disk']['mag'][0]
        
        psf_xcntr_err = fit_info['disk']['xctr'][1]
        psf_ycntr_err = fit_info['disk']['yctr'][1]
        mag_p_err = fit_info['disk']['mag'][1]
    else:
        psf_xcntr = -999.
        psf_ycntr = -999.
        mag_p = -999.
        
        psf_xcntr_err = -999.
        psf_ycntr_err = -999.
        mag_p_err = -999.

    if 'bar' in fit_info:
        bar_xcntr = fit_info['bulge']['xctr'][0]
        bar_ycntr = fit_info['bulge']['yctr'][0]
        mag_bar = fit_info['bulge']['mag'][0]
        rbar = fit_info['bulge']['rad'][0]
        BarIndex = fit_info['bulge']['n'][0]
        BarEllipticity = fit_info['bulge']['ell'][0]
        BarPA = fit_info['bulge']['pa'][0]
        BarBoxy = fit_info['bulge']['boxy'][0]

        bar_xcntr_err = fit_info['bulge']['xctr'][1]
        bar_ycntr_err = fit_info['bulge']['yctr'][1]
        mag_bar_err = fit_info['bulge']['mag'][1]
        rbar_err = fit_info['bulge']['rad'][1]
        BarIndexErr = fit_info['bulge']['n'][1]
        BarEllipticityErr = fit_info['bulge']['ell'][1]
        BarPAErr = fit_info['bulge']['pa'][1]
        BarBoxyErr = fit_info['bulge']['boxy'][1]
    else:
        bar_xcntr = -999.
        bar_ycntr = -999.
        mag_bar = -999.
        rbar = -999.
        BarIndex = -999.
        BarEllipticity = -999.
        BarPA = -999.
        BarBoxy = -999.

        bar_xcntr_err = -999.
        bar_ycntr_err = -999.
        mag_bar_err = -999.
        rbar_err = -999.
        BarIndexErr = -999. 
        BarEllipticityErr = -999.
        BarPAErr = -999.
        BarBoxyErr = -999.
    if 'sky' in fit_info:
        galfit_sky = fit_info['sky']['mag'][0]
        galfit_sky_err = fit_info['sky']['mag'][1]
    else:
        galfit_sky = -999.
        galfit_sky_err = -999.
                                       
    # Converting fitted params to physical params
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
        if 'point' in ComP:
            fwhm_kpc = 0.5 * phy_parms[3]
        else:
            fwhm_kpc = 9999
    # Finding derived parameters
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
    # Start writing html file. Now the template keywords will get values
    pngfile = 'P_' + c.fstring + '.png'
    Neighbour_Sersic = ''
    Object_Sersic = ''
    Object_Sersic_err = ''
    Object_Exp = ''
    Object_Exp_err = ''
    Point_Vals = ''
    Point_Vals_err = ''
    try:
        for key in fit_info.keys():
            if 'neighbor' in key:
                Neighbour_Sersic = str(Neighbour_Sersic) + \
                                   '<TR align="center" bgcolor=' + \
                                   '"#99CCFF"><TD> neighbor sersic' + \
                                   ' </TD> <TD> ' + \
                                   str(fit_info[key]['xctr'][0]) + '</TD> <TD> '\
                                   + str(fit_info[key]['yctr'][0]) + \
                                   ' </TD> <TD> ' + str(fit_info[key]['mag'][0]) + \
                                   ' </TD> <TD> ' + \
                                   str(fit_info[key]['rad'][0]) + ' </TD> <TD> ' + \
                                   ' ' + ' </TD> <TD> ' + \
                                   str(fit_info[key]['n'][0]) + ' </TD> <TD> ' +\
                                   str(fit_info[key]['ell'][0]) + ' </TD> <TD> ' +\
                                   str(fit_info[key]['pa'][0]) + ' </TD> <TD> ' + \
                                   str(fit_info[key]['boxy'][0]) + ' </TD> </TR>'
                Neighbour_Sersic = str(Neighbour_Sersic) + \
                                   '<TR align="center" ' + \
                                   'bgcolor="#CCFFFF"> <TD>' + ' ' + \
                                   ' </TD> <TD> ' + \
                                   str(fit_info[key]['xctr'][1]) + '</TD> <TD> '\
                                   + str(fit_info[key]['yctr'][1]) + \
                                   ' </TD> <TD> ' + str(fit_info[key]['mag'][1]) + \
                                   ' </TD> <TD> ' + \
                                   str(fit_info[key]['rad'][1]) + ' </TD> <TD> ' + \
                                   ' ' + ' </TD> <TD> ' + \
                                   str(fit_info[key]['n'][1]) + ' </TD> <TD> ' +\
                                   str(fit_info[key]['ell'][1]) + ' </TD> <TD> ' +\
                                   str(fit_info[key]['pa'][1]) + ' </TD> <TD> ' + \
                                   str(fit_info[key]['boxy'][1]) + ' </TD> </TR>'
            if 'bulge' in key:
                Object_Sersic = '<TR align="center" ' +\
                                'bgcolor="#99CCFF">' +\
                                '<TD> sersic bulge</TD> <TD> ' +\
                                str(fit_info[key]['xctr'][0]) + '</TD> <TD> '\
                                + str(fit_info[key]['yctr'][0]) + \
                                ' </TD> <TD> ' + str(fit_info[key]['mag'][0]) + \
                                ' </TD> <TD> ' + \
                                str(fit_info[key]['rad'][0]) + ' </TD> <TD> ' + \
                                str(round(re_kpc, 3))[:5] + ' </TD> <TD> ' + \
                                str(fit_info[key]['n'][0]) + ' </TD> <TD> ' +\
                                str(fit_info[key]['ell'][0]) + ' </TD> <TD> ' +\
                                str(fit_info[key]['pa'][0]) + ' </TD> <TD> ' + \
                                str(fit_info[key]['boxy'][0]) + ' </TD> </TR>'
                Object_Sersic_err = '<TR align="center" ' + \
                                    'bgcolor="#CCFFFF">' + \
                                    '<TD>' + ' ' + '</TD> <TD>' + \
                                    str(fit_info[key]['xctr'][1]) + '</TD> <TD> '\
                                    + str(fit_info[key]['yctr'][1]) + \
                                    ' </TD> <TD> ' + str(fit_info[key]['mag'][1]) + \
                                    ' </TD> <TD> ' + \
                                    str(fit_info[key]['rad'][1]) + ' </TD> <TD> ' + \
                                    str(round(re_err_kpc, 3))[:5] + ' </TD> <TD> ' + \
                                    str(fit_info[key]['n'][1]) + ' </TD> <TD> ' +\
                                    str(fit_info[key]['ell'][1]) + ' </TD> <TD> ' +\
                                    str(fit_info[key]['pa'][1]) + ' </TD> <TD> ' + \
                                    str(fit_info[key]['boxy'][1]) + ' </TD> </TR>'
            if 'disk' in key:
                Object_Exp = '<TR align="center" bgcolor="#99CCFF">' +\
                             '<TD> disk </TD> <TD> ' + \
                             str(fit_info[key]['xctr'][0]) + '</TD> <TD> '\
                             + str(fit_info[key]['yctr'][0]) + \
                             ' </TD> <TD> ' + str(fit_info[key]['mag'][0]) + \
                             ' </TD> <TD> ' + \
                             str(fit_info[key]['rad'][0]) + ' </TD> <TD> ' + \
                             str(round(rd_kpc, 3))[:5] + ' </TD> <TD> ' + \
                             ' </TD> <TD> </TD> <TD> ' +\
                             str(fit_info[key]['ell'][0]) + ' </TD> <TD> ' +\
                             str(fit_info[key]['pa'][0]) + ' </TD> <TD> ' + \
                             str(fit_info[key]['boxy'][0]) + ' </TD> </TR>'

                Object_Exp_err = '<TR align="center" ' + \
                                 'bgcolor="#CCFFFF">' + \
                                 '<TD>' + ' ' + '</TD> <TD>' + \
                                 str(fit_info[key]['xctr'][1]) + '</TD> <TD> '\
                                 + str(fit_info[key]['yctr'][1]) + \
                                 ' </TD> <TD> ' + str(fit_info[key]['mag'][1]) + \
                                 ' </TD> <TD> ' + \
                                 str(fit_info[key]['rad'][1]) + ' </TD> <TD> ' + \
                                 str(round(rd_err_kpc, 3))[:5] + ' </TD> <TD> ' + \
                                 ' </TD> <TD> </TD> <TD> ' +\
                                 str(fit_info[key]['ell'][1]) + ' </TD> <TD> ' +\
                                 str(fit_info[key]['pa'][1]) + ' </TD> <TD> ' + \
                                 str(fit_info[key]['boxy'][1]) + ' </TD> </TR>'


            if(str(values[0]) == 'psf'):
                Point_Vals = '<TR align="center" bgcolor="#99CCFF">' + \
                             '<TD> point </TD> <TD> ' + \
                             str(fit_info[key]['xctr'][0]) + '</TD> <TD> '\
                             + str(fit_info[key]['yctr'][0]) + \
                             ' </TD> <TD> ' + str(fit_info[key]['mag'][0]) +\
                             ' </TD> <TD> ' + str('9999') +  ' </TD> <TD> ' + \
                             str('9999') + ' </TD> <TD> ' + \
                             ' ' + ' </TD> <TD> ' + str('9999') +  \
                             ' </TD> <TD> ' + str('9999') + \
                             ' </TD> <TD> ' + str('9999') + ' </TD></TR>'
                Point_Vals_err = '<TR align="center" ' + \
                                 'bgcolor="#CCFFFF">' + \
                                 '<TD>' + ' ' + '</TD> <TD>' + \
                                 str(fit_info[key]['xctr'][1]) + '</TD> <TD> '\
                                 + str(fit_info[key]['yctr'][1]) + \
                             ' </TD> <TD> ' + str(fit_info[key]['mag'][1]) +\
                             ' </TD> <TD> ' + str('9999') +  ' </TD> <TD> ' + \
                             str('9999') + ' </TD> <TD> ' + \
                             ' ' + ' </TD> <TD> ' + str('9999') +  \
                             ' </TD> <TD> ' + str('9999') + \
                             ' </TD> <TD> ' + str('9999') + ' </TD></TR>'
                            
    except:
        pass
    if 'bulge' in ComP:
        try:
              pixelscale = c.pixelscale
        except:
              pixelscale = 1
        AvgMagInsideRe = mag_b + 2.5 * n.log10(2 * 3.14 * pixelscale * \
                      pixelscale * re * re * n.sqrt(1 - SersicEllipticity**2.0))
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
    error_mesg7 = ''
    if c.starthandle:
        error_mesg6 = '<a href="R_' + c.fstring + \
	              '_1.html"> Crashed </a>' 

    # Now test for fitting problems and set flags for analysis
    FitFlag = 0
    HitLimit = 0

    if not c.detail:
        if 'bulge' in ComP:
            if abs(mag_b - c.UMag) < 0.2 or abs(mag_b - c.LMag) < 0.2):
                FitFlag += 2**Get_FitFlag('IE_AT_LIMIT')
            if abs(re - c.LRe) < 0.1 or abs(re - c.URe) < 1.0:
                FitFlag += 2**Get_FitFlag('RE_AT_LIMIT')
            if abs(SersicIndex - c.LN) < 0.03 or abs(SersicIndex - c.UN) < 0.5:
                FitFlag += 2**Get_FitFlag('N_AT_LIMIT')
            if abs(SersicEllipticity - 0.0) < 0.05 or abs(SersicEllipticity - 0.0) > 0.95:
                FitFlag += 2**Get_FitFlag('EB_AT_LIMIT')
        if 'disk' in ComP:
            if abs(mag_d - c.UMag) < 0.2 or abs(mag_d - c.LMag) < 0.2:
                FitFlag += 2**Get_FitFlag('ID_AT_LIMIT')
            if abs(rd - c.LRd) < 0.1 or abs(rd - c.URd) < 1.0:
                FitFlag += 2**Get_FitFlag('RD_AT_LIMIT')
            if abs(DiskEllipticity - 0.0) < 0.05 or abs(DiskEllipticity - 0.0) > 0.95:
                FitFlag += 2**Get_FitFlag('ED_AT_LIMIT')

    if not FitFlag:
        error_mesg4 = str(error_mesg4) + 'One of the parameters'
        error_mesg5 = str(error_mesg5) + '          hits limit!'
    
    if Goodness < c.Goodness:
        error_mesg2 = str(error_mesg2) + 'Goodness is poor!'
        FitFlag += 2**Get_FitFlag('SMALL_GOODNESS')
    if chi2nu > c.chi2sq:
        error_mesg1 = str(error_mesg1) + 'Chi2nu is large!'
        if chi2nu != 9999:
            FitFlag += 2**Get_FitFlag('LARGE_CHISQ')
    if abs(bulge_xcntr - xcntr) > c.center_deviation or \
           abs(bulge_ycntr - ycntr) > c.center_deviation or \
           abs(disk_xcntr - xcntr) > c.center_deviation or \
           abs(disk_ycntr - ycntr) > c.center_deviation:
        error_mesg3 = str(error_mesg3) + 'Fake Center!'
        if bulge_xcntr == 9999 or bulge_ycntr == 9999 or \
               disk_xcntr == 9999 or disk_ycntr == 9999:
            pass
        else:
            FitFlag += 2**Get_FitFlag('FAKE_CNTR')

    if FitFlag > 0:
        img_notify = str(c.PYMORPH_PATH) + '/html/goodfit.gif'
        good_fit = 1
    else:
        img_notify = str(c.PYMORPH_PATH) + '/html/badfit.gif'
        good_fit = 0
    outfile = open('R_' + c.fstring + '.html','w')
    outfile.write(template %vars())
    outfile.close()
    # Finding number of runs in the csv file 
    run = 1
    if exists('result.csv'):
        for line_res in csv.reader(open('result.csv').readlines()[1:]):    
            if(str(line_res[0]) == c.fstring):
                run += 1
    # Writing csv file 
    f_res = open("result.csv", "ab")
    writer = csv.writer(f_res)
    galid = c.fstring
    ParamToWrite = [galid, alpha_j, delta_j, z, c.SexMagAuto, c.SexMagAutoErr,
                    c.SexTargets]
    if 'bulge' in ComP:
        for bulgecomp in [bulge_xcntr,bulge_xcntr_err,bulge_ycntr,bulge_ycntr_err,
                          mag_b, mag_b_err, re, re_err, re_kpc, re_err_kpc, 
                          SersicIndex, SersicIndexErr, AvgMagInsideRe,\
                          AvgMagInsideReErr, SersicEllipticity, \
                          SersicEllipticityErr, SersicPA, 
			  SersicPAErr, SersicBoxy, SersicBoxyErr]:
            ParamToWrite.append(bulgecomp)
    else:
        for bulgecomp in [9999, 9999, 9999, 9999, 9999, 9999, 9999, 9999, \
                          9999, 9999, 9999, 9999, 9999, 9999, 9999, 9999, \
                          9999, 9999, 9999]:
            ParamToWrite.append(bulgecomp)
    if 'disk' in ComP:
        for diskcomp in [disk_xcntr,disk_xcntr_err,disk_ycntr,disk_ycntr_err,
                         mag_d, mag_d_err, rd, rd_err, rd_kpc, rd_err_kpc, 
                         DiskEllipticity, DiskEllipticityErr, DiskPA, 
                         DiskPAErr, DiskBoxy, DiskBoxyErr]:
            ParamToWrite.append(diskcomp)
    else:
        for diskcomp in [9999, 9999, 9999, 9999, 9999, 9999, 9999, 9999,
                         9999, 9999, 9999, 9999, 9999, 9999, 9999, 9999]:
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
        print "galfit sky was wierd: ", galfit_sky
        galfit_sky = 9999
        galfit_sky_err = 9999
        print 'GALFIT does not report sky'
    if c.GalSky != 9999:
        galfit_sky = c.GalSky
    for otherparam in [chi2nu, Goodness, run, C, C_err, A, A_err, S, S_err, G,\
                       M, c.SexSky, galfit_sky, galfit_sky_err, DisMoD, \
                       distance, good_fit, c.Flag, FitFlag, c.SexHalfRad]:
        ParamToWrite.append(otherparam)
    if 'bar' in ComP:
        for barcomp in [bar_xcntr,bar_xcntr_err,bar_ycntr,bar_ycntr_err,
                        mag_bar, mag_bar_err, re_bar, re_bar_err, re_bar_kpc, \
                        re_bar_err_kpc, \
                        SersicIndexBar, SersicIndexBarErr, \
                        SersicEllipticityBar, \
                        SersicEllipticityBarErr, SersicBoxyBar]:
            ParamToWrite.append(barcomp)
    else:
        for barcomp in [9999, 9999, 9999, 9999, 9999, 9999, 9999, 9999, \
                          9999, 9999, 9999, 9999,9999, 9999, 9999]:
            ParamToWrite.append(barcomp)
    # Remove any nan or inf from the parameter
    for p in ParamToWrite:
        if str(p) in ('nan', '-nan', 'inf', '-inf'):
            ParamToWrite[ParamToWrite.index(p)] = 9999
        else:
            pass
    writer.writerow(ParamToWrite)
    f_res.close()
    # Writing data base
    try:
        from utilities import WriteDb
    except:
        print 'DB writing Problem'
    try:
        WriteDb(ParamToWrite)
    except:
        print 'No database can be created!'
        traceback.print_exc()
    # Writing quick view html pymorph.html
    outfile1 = open('pymorph.html', 'w')
    outfile1.write('<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01' \
                   ' Transitional//EN">')
    outfile1.write('<HTML><HEAD> \n')
    outfile1.writelines(['<META HTTP-EQUIV="Refresh" CONTENT="10; \
                    URL=pymorph.html"> \n'])
    outfile1.write('</HEAD> \n')
    outfile = 'R_' + c.fstring + '.html'
    try:
        for line in open(outfile, 'r'):
            if(str(line) != '<html>'):
                outfile1.writelines([str(line)])
    except:
        pass
    outfile1.close()



def getline(values):
    line = values.pop(0)
    line = line.split()
    return line

def getfit(f):
    values = f.read()
    values = values.replace('\n\n','\n') # remove unnecessary whitespace
    # now replace any unneeded charaters that might affect the I/O process...
    for bad_char in ['(',')','[',']',',','*', '--']:
        values = values.replace(bad_char,' ') 

    return values


def load_component(data_line, err_line):
    """This function will construct and load dictionaries for a object when passed the already split
    data_line containing fit parameters and err_line containing errors on the fit parameters"""
    # construct object dictionary
    # these are all fitted values, but many may be unused for a particular object type
    obj = {'xctr':[-999.,-999.],
           'yctr':[-999.,-999.],
           'mag':[-999.,-999.],
           'rad':[-999.,-999.],
           'n':[-999.,-999.],
           'ell':[-999.,-999.],
           'pa':[-999.,-999.],
           'boxy':[-999.,-999.]
           }
    
    # now load data
    obj['xctr'][0]=float(str(data_line[2]))
    obj['yctr'][0]=float(str(data_line[3]))
    obj['mag'][0]=float(data_line[4])

    if data_line[0] == "sky": 
        obj['mag'][1]=float(err_line[0])

    else:
        obj['xctr'][1]=float(str(err_line[0]))
        obj['yctr'][1]=float(str(err_line[1]))
        obj['mag'][1]=float(err_line[2])

        if data_line[0] in ['sersic','expdisk']:
            obj['rad'][0]=float(data_line[5])
            obj['rad'][1]=float(err_line[3])

            if data_line[0] = 'sersic':
                obj['n'][0]=float(data_line[6])
                obj['n'][1]=float(err_line[4])
                pos = 7
            else:
                pos = 6
            
            
            obj['ell'][0]=float(data_line[pos])
            obj['ell'][1]=float(err_line[pos-2])
            pos +=1
            
            obj['pa'][0]=float(data_line[pos])
            obj['pa'][1]=float(err_line[pos-2])
            pos +=1

            try:
                obj['boxy'][0]=float(data_line[pos])
                obj['boxy'][1]=float(err_line[pos-2])
            except IndexError: #if Galfit version 3.0 or later, then no boxyness is reported
                obj['boxy'][0]=-999.
                obj['boxy'][1]=-999.                

    return obj

def read_fitlog(filename = 'fit.log', yes_bar = 0):
    """ This function will read the fit log and return all the relevant
    information in 2 Dictionaries, 1 with the basic info and one with
    the fit info"""

    neighbor = 0
    basic_info = {}
    fit_info = {}
    if exists(filename):
        f = open(filename,'r')
        values = getfit(f)
        f.close()
        while len(values) > 0:
            line = getline(values)
            try: 
                if(str(line[0]) == 'Input'):
                    basic_info['Input'] = 1
                elif(str(line[0]) == 'Init.'):
                    basic_info['initial_conf'] = str(line[4])
                elif(str(line[0]) == 'Restart'):
                    basic_info['restart_conf'] = str(line[3])
                elif(str(line[0]) == 'Chi^2/nu'):
                    basic_info['chi2nu'] = float(line[2])
                elif(str(line[0]) in ['Output', 'Chi^2']):
                    continue
                # for galaxy bulge
                else: # it must be part of the fit...
                    if str(line[0]) == 'sersic':
                        if 'bulge' not in fit_info:
                            key = 'bulge'
                        elif ('bar' not in fit_info and yes_bar):
                            key = 'bar'
                        else:
                            neighbor +=1
                            key = 'neighbor' + str(neighbor) 
                    elif (str(line[0]) == 'expdisk'):
                        key = 'disk'
                    elif(str(line[0]) == 'psf'):
                        key = 'psf'
                    elif(str(line[0]) == 'sky'):
                        key = 'sky'
                    
                    err_line = getline(values)
                    
                    fit_info[key]= load_component(line, err_line)    
                    
            except: 
                pass
            
    else:
        print "File does not exist!!!!"
    
    return     basic_info, fit_info
