from os.path import exists
import csv
import sys
import numpy as n
import fileinput
from cosmocal import cal 
import datetime
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

def WriteParams(ParamNamesToWrite, cutimage, xcntr, ycntr, distance, alpha_j, delta_j, z, Goodness,
                C, C_err, A, A_err, S, S_err, G, M, EXPTIME):

    # this dictionary will hold any parameters that may be printed
    all_params = dict((ParamNamesToWrite[key][0], -999) for key in ParamNamesToWrite.keys())

    # load some of the passed in parameters
    all_params['Name'] = c.fstring
    all_params['ra_gal'] = alpha_j
    all_params['dec_gal'] = delta_j
    all_params['z'] = z
    all_params['Goodness'] = float(str(round(Goodness, 3))[:5])
    all_params['C'] = C
    all_params['C_err'] = C_err 
    all_params['A'] = A
    all_params['A_err'] = A_err
    all_params['S'] = S
    all_params['S_err'] = S_err
    all_params['G'] = G
    all_params['M'] = M
    all_params['distance'] = float(str(round(distance, 3))[:5])
    all_params['mag_auto'] = c.SexMagAuto
    all_params['magerr_auto'] = c.SexMagAutoErr
    all_params['num_targets'] = c.SexTargets
    all_params['SexSky'] = c.SexSky
    all_params['flag'] = c.Flag
    all_params['SexHalfRad'] = c.SexHalfRad
    all_params['magzp'] = c.mag_zero

    # Now continue
    try:
        ComP = c.components
    except:
        ComP = ['bulge', 'disk']
    if len(ComP) == 0:
        ComP = ['bulge', 'disk']
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
    if abs(ra1) > 360:
        alpha1, alpha2, alpha3 = 9999, '', ''
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
    if abs(dec1) > 90:
        delta1, delta2, delta3 = 9999, '', ''
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
    else:
        alpha_ned = ''
        delta_ned = ''
        
    initial_conf = basic_info['initial_conf']
    restart_conf = basic_info['restart_conf']
    print restart_conf 
    # move the restart file to a reasonably named output file
    new_outname = initial_conf.replace('in','out')
    os.rename(restart_conf, new_outname)
    basic_info['restart_conf'] = new_outname

    
    all_params['chi2nu'] = basic_info['chi2nu']
    if 'bulge' in fit_info:
        all_params['bulge_xctr'] = fit_info['bulge']['xctr'][0]
        all_params['bulge_yctr'] = fit_info['bulge']['yctr'][0]
        all_params['Ie'] = fit_info['bulge']['mag'][0]
        all_params['re_pix'] = fit_info['bulge']['rad'][0]
        all_params['n'] = fit_info['bulge']['n'][0]
        all_params['eb'] = fit_info['bulge']['ell'][0]
        all_params['bpa'] = fit_info['bulge']['pa'][0]
        all_params['bboxy'] = fit_info['bulge']['boxy'][0]

        all_params['bulge_xctr_err'] = fit_info['bulge']['xctr'][1]
        all_params['bulge_yctr_err'] = fit_info['bulge']['yctr'][1]
        all_params['Ie_err'] = fit_info['bulge']['mag'][1]
        all_params['re_pix_err'] = fit_info['bulge']['rad'][1]
        all_params['n_err'] = fit_info['bulge']['n'][1]
        all_params['eb_err'] = fit_info['bulge']['ell'][1]
        all_params['bpa_err'] = fit_info['bulge']['pa'][1]
        all_params['bboxy_err'] = fit_info['bulge']['boxy'][1]
    
    if 'disk' in fit_info:
        all_params['disk_xctr'] = fit_info['disk']['xctr'][0]
        all_params['disk_yctr'] = fit_info['disk']['yctr'][0]
        all_params['Id'] = fit_info['disk']['mag'][0]
        all_params['rd_pix'] = fit_info['disk']['rad'][0]
        all_params['ed'] = fit_info['disk']['ell'][0]
        all_params['dpa'] = fit_info['disk']['pa'][0]
        all_params['dboxy'] = fit_info['disk']['boxy'][0]

        all_params['disk_xctr_err'] = fit_info['disk']['xctr'][1]
        all_params['disk_yctr_err'] = fit_info['disk']['yctr'][1]
        all_params['Id_err'] = fit_info['disk']['mag'][1]
        all_params['rd_pix_err'] = fit_info['disk']['rad'][1]
        all_params['ed_err'] = fit_info['disk']['ell'][1]
        all_params['dpa_err'] = fit_info['disk']['pa'][1]
        all_params['dboxy_err'] = fit_info['disk']['boxy'][1]
        
    if 'point' in fit_info:
        all_params['p_xctr'] = fit_info['point']['xctr'][0]
        all_params['p_yctr'] = fit_info['point']['yctr'][0]
        all_params['Ip'] = fit_info['point']['mag'][0]
        
        all_params['p_xctr_err'] = fit_info['point']['xctr'][1]
        all_params['p_yctr_err'] = fit_info['point']['yctr'][1]
        all_params['Ip_err'] = fit_info['point']['mag'][1]
        
    if 'bar' in fit_info:
        all_params['bar_xctr'] = fit_info['bar']['xctr'][0]
        all_params['bar_yctr'] = fit_info['bar']['yctr'][0]
        all_params['Ibar'] = fit_info['bar']['mag'][0]
        all_params['rbar_pix'] = fit_info['bar']['rad'][0]
        all_params['n_bar'] = fit_info['bar']['n'][0]
        all_params['ebar'] = fit_info['bar']['ell'][0]
        all_params['barpa'] = fit_info['bar']['pa'][0]
        all_params['barboxy'] = fit_info['bar']['boxy'][0]

        all_params['bar_xctr_err'] = fit_info['bar']['xctr'][1]
        all_params['bar_yctr_err'] = fit_info['bar']['yctr'][1]
        all_params['Ibar_err'] = fit_info['bar']['mag'][1]
        all_params['rbar_pix_err'] = fit_info['bar']['rad'][1]
        all_params['nbar_err'] = fit_info['bar']['n'][1]
        all_params['ebar_err'] = fit_info['bar']['ell'][1]
        all_params['barpa_err'] = fit_info['bar']['pa'][1]
        all_params['barboxy_err'] = fit_info['bar']['boxy'][1]

    if 'sky' in fit_info:
        all_params['GalSky'] = fit_info['sky']['mag'][0]
        all_params['GalSky_err'] = fit_info['sky']['mag'][1]
                                           
    # Converting fitted params to physical params
    if(z != 9999 and z > 0):
        phy_parms = cal(z, c.H0, c.WM, c.WV, c.pixelscale)
        all_params['dis_modu'] = phy_parms[2]
        if 'bulge' in ComP:
            all_params['re_kpc'] = phy_parms[3] * all_params['re_pix']
            all_params['re_kpc_err'] = phy_parms[3] * all_params['re_pix_err']
        if 'disk' in ComP:
            all_params['rd_kpc'] = phy_parms[3] * all_params['rd_pix']
            all_params['rd_kpc_err'] = phy_parms[3] * all_params['rd_pix_err']
        if 'bar' in ComP:
            all_params['rbar_kpc'] = phy_parms[3] * all_params['rbar_pix']
            all_params['rbar_kpc_err'] = phy_parms[3] * all_params['rbar_pix_err']
        if 'point' in ComP:
            all_params['Pfwhm_kpc'] = 0.5 * phy_parms[3]
    # Finding derived parameters
    if 'bulge' in ComP and 'disk' in ComP:
        fb = 10**(-0.4 * (all_params['Ie'] - c.mag_zero))
        fd = 10**(-0.4 * (all_params['Id'] - c.mag_zero))
        all_params['BD'] = fb / fd 
        if 'point' in ComP:
            fp = 10**(-0.4 * (all_params['Ip'] - c.mag_zero))
        else:
            fp = 0.0
        if 'bar' in ComP:
            fbar = 10**(-0.4 * (all_params['Ibar'] - c.mag_zero)) 
        else:
            fbar = 0.0
        
        all_params['BT'] = fb / (fb + fd + fp + fbar)
    elif 'bulge' in ComP:
        all_params['BT'] = 1.0
    elif 'disk' in ComP:
        all_params['BD'] = 0.0
        all_params['BT'] = 0.0
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
                                str(round(all_params['re_kpc'], 3))[:5] + ' </TD> <TD> ' + \
                                str(fit_info[key]['n'][0]) + ' </TD> <TD> ' +\
                                str(fit_info[key]['ell'][0]) + ' </TD> <TD> ' +\
                                str(fit_info[key]['pa'][0]) + ' </TD> <TD> ' + \
                                str(0) + ' </TD> </TR>'
                Object_Sersic_err = '<TR align="center" ' + \
                                    'bgcolor="#CCFFFF">' + \
                                    '<TD>' + ' ' + '</TD> <TD>' + \
                                    str(fit_info[key]['xctr'][1]) + '</TD> <TD> '\
                                    + str(fit_info[key]['yctr'][1]) + \
                                    ' </TD> <TD> ' + str(fit_info[key]['mag'][1]) + \
                                    ' </TD> <TD> ' + \
                                    str(fit_info[key]['rad'][1]) + ' </TD> <TD> ' + \
                                    str(round(all_params['re_kpc_err'], 3))[:5] + ' </TD> <TD> ' + \
                                    str(fit_info[key]['n'][1]) + ' </TD> <TD> ' +\
                                    str(fit_info[key]['ell'][1]) + ' </TD> <TD> ' +\
                                    str(fit_info[key]['pa'][1]) + ' </TD> <TD> ' + \
                                    str(0) + ' </TD> </TR>'
            if 'disk' in key:
                Object_Exp = '<TR align="center" bgcolor="#99CCFF">' +\
                             '<TD> disk </TD> <TD> ' + \
                             str(fit_info[key]['xctr'][0]) + '</TD> <TD> '\
                             + str(fit_info[key]['yctr'][0]) + \
                             ' </TD> <TD> ' + str(fit_info[key]['mag'][0]) + \
                             ' </TD> <TD> ' + \
                             str(fit_info[key]['rad'][0]) + ' </TD> <TD> ' + \
                             str(round(all_params['rd_kpc'], 3))[:5] +\
                             ' </TD> <TD> </TD> <TD> ' +\
                             str(fit_info[key]['ell'][0]) + ' </TD> <TD> ' +\
                             str(fit_info[key]['pa'][0]) + ' </TD> <TD> ' + \
                             str(0) + ' </TD> </TR>'

                Object_Exp_err = '<TR align="center" ' + \
                                 'bgcolor="#CCFFFF">' + \
                                 '<TD>' + ' ' + '</TD> <TD>' + \
                                 str(fit_info[key]['xctr'][1]) + '</TD> <TD> '\
                                 + str(fit_info[key]['yctr'][1]) + \
                                 ' </TD> <TD> ' + str(fit_info[key]['mag'][1]) + \
                                 ' </TD> <TD> ' + \
                                 str(fit_info[key]['rad'][1]) + ' </TD> <TD> ' + \
                                 str(round(all_params['rd_kpc_err'], 3))[:5] + \
                                 ' </TD> <TD> </TD> <TD> ' +\
                                 str(fit_info[key]['ell'][1]) + ' </TD> <TD> ' +\
                                 str(fit_info[key]['pa'][1]) + ' </TD> <TD> ' + \
                                 str(0) + ' </TD> </TR>'


            if 'point' in key:
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
    except Exception, inst:
        print type(inst)     # the exception instance
        print inst.args      # arguments stored in\
                                             # .args
        print inst           # __str__ allows args\
                                             # to printed directly
        print "something bad happened writing!!!!\n\n"
        print traceback.print_exc()                   
    if 'bulge' in ComP:
        try:
              pixelscale = c.pixelscale
        except:
              pixelscale = 1
        all_params['AvgIe'] = all_params['Ie'] + 2.5 * n.log10(2 * 3.14 * pixelscale * \
                      pixelscale *  all_params['re_pix'] *  all_params['re_pix'] * n.sqrt(1 -  all_params['eb']**2.0))
        AvgMagInsideReErr2 = (1.085 * n.sqrt((2 *  all_params['re_pix'] *  all_params['re_pix_err'])**2.0 + \
                             (( all_params['eb'] *  all_params['eb_err']) / \
                             n.sqrt(1 -  all_params['eb']**2.0))**2.0)) / \
                             (n.sqrt(1 -  all_params['eb']**2.0) * 2 * 3.14 * \
                              all_params['re_pix'] *  all_params['re_pix'])
        all_params['AvgIe_err'] = n.sqrt( all_params['Ie_err']**2.0 + AvgMagInsideReErr2**2.0)
    wC = str(C)[:5]
    wA = str(A)[:5]
    wS = str(S)[:5]
    wG = str(G)[:5]
    wM = str(M)[:5]
    wBD = str( all_params['BD'])[:5]
    wBT = str( all_params['BT'])[:5]
    wAvgMagInsideRe = str( all_params['AvgIe'])[:5]
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
    all_params['FitFlag'] = 0
    HitLimit = 0

    if not c.detail:
        if 'bulge' in ComP:
            if abs(all_params['Ie'] - c.UMag) < 0.2 or abs(all_params['Ie'] - c.LMag) < 0.2:
                all_params['FitFlag'] += 2**Get_FitFlag('IE_AT_LIMIT')
            if abs(all_params['re_pix'] - c.LRe) < 0.1 or abs(all_params['re_pix'] - c.URe) < 1.0:
                all_params['FitFlag'] += 2**Get_FitFlag('RE_AT_LIMIT')
            if abs(all_params['n'] - c.LN) < 0.03 or abs(all_params['n'] - c.UN) < 0.5:
                all_params['FitFlag'] += 2**Get_FitFlag('N_AT_LIMIT')
            if abs(all_params['eb'] - 0.0) < 0.05 or abs(all_params['eb'] - 0.0) > 0.95:
                all_params['FitFlag'] += 2**Get_FitFlag('EB_AT_LIMIT')
        if 'disk' in ComP:
            if abs(all_params['Id'] - c.UMag) < 0.2 or abs(all_params['Id'] - c.LMag) < 0.2:
                all_params['FitFlag'] += 2**Get_FitFlag('ID_AT_LIMIT')
            if abs(all_params['rd_pix'] - c.LRd) < 0.1 or abs(all_params['rd_pix'] - c.URd) < 1.0:
                all_params['FitFlag'] += 2**Get_FitFlag('RD_AT_LIMIT')
            if abs(all_params['ed'] - 0.0) < 0.05 or abs(all_params['ed'] - 0.0) > 0.95:
                all_params['FitFlag'] += 2**Get_FitFlag('ED_AT_LIMIT')

    if not all_params['FitFlag']:
        error_mesg4 = str(error_mesg4) + 'One of the parameters'
        error_mesg5 = str(error_mesg5) + '          hits limit!'
    
    if all_params['Goodness'] < c.Goodness:
        error_mesg2 = str(error_mesg2) + 'Goodness is poor!'
        all_params['FitFlag'] += 2**Get_FitFlag('SMALL_GOODNESS')
    if all_params['chi2nu'] > c.chi2sq:
        error_mesg1 = str(error_mesg1) + 'Chi2nu is large!'
        if all_params['chi2nu'] != 9999:
            all_params['FitFlag'] += 2**Get_FitFlag('LARGE_CHISQ')
    if abs(all_params['bulge_xctr'] - xcntr) > c.center_deviation or \
           abs(all_params['bulge_yctr'] - ycntr) > c.center_deviation or \
           abs(all_params['disk_xctr'] - xcntr) > c.center_deviation or \
           abs(all_params['disk_yctr'] - ycntr) > c.center_deviation:
        error_mesg3 = str(error_mesg3) + 'Fake Center!'
        if all_params['bulge_xctr'] == 9999 or all_params['bulge_yctr'] == 9999 or \
               all_params['disk_xctr'] == 9999 or all_params['disk_yctr'] == 9999:
            pass
        else:
            all_params['FitFlag'] += 2**Get_FitFlag('FAKE_CNTR')

    if all_params['FitFlag'] > 0:
        img_notify = str(c.PYMORPH_PATH) + '/html/goodfit.gif'
        good_fit = 1
    else:
        img_notify = str(c.PYMORPH_PATH) + '/html/badfit.gif'
        good_fit = 0
    chi2nu = all_params['chi2nu']
    Distance = all_params['distance']
    # Finding number of runs in the csv file 
    all_params['run'] = 1
    if exists('result.csv'):
        for line_res in csv.reader(open('result.csv').readlines()[1:]):    
            if(str(line_res[0]) == c.fstring):
                all_params['run'] += 1
    if c.GalSky != 9999:
        all_params['GalSky'] = c.GalSky
    
    # Writing data base
    try:
        from utilities import WriteDb
    except:
        print 'DB writing Problem'
    try:
        WriteDb(ParamNamesToWrite, all_params)
    except:
        print 'No database can be created!'
        traceback.print_exc()

    # Writing csv file 
    f_res = open("result.csv", "ab")
    writer = csv.writer(f_res)
    writer.writerow([all_params[ParamNamesToWrite[key][0]] for key in ParamNamesToWrite.keys()])
    f_res.close()

    #writing html
    outfile = open('R_' + c.fstring + '.html','w')
    outfile.write(template %vars())
    outfile.close()

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

    values = values.split('\n')
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

            if data_line[0] == 'sersic':
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
                        key = 'point'
                    elif(str(line[0]) == 'sky'):
                        key = 'sky'
                    
                    err_line = getline(values)
                    
                    fit_info[key]= load_component(line, err_line)    
                    
            except: 
                pass
            
    else:
        print "File does not exist!!!!"
    
    return     basic_info, fit_info
