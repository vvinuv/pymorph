import os
import csv
import datetime
import numpy as np
import fileinput
from .cosmocal import CosmoCal 
import traceback
from .flagfunc import *
from .pymorphutils import RaDegToHMS, DecDegToDMS, output_params
try:
    from .writedbfunc import WriteDB
except:
    print('No mysql database or python mysql.connector module')

class WriteHtmlCSV(object):
    """The class which will write html and csv output. This class will also 
       check whether the fit is good or bad using the Chisq and Goodness value
       It will also notify the goodness/badness of fit"""
    def __init__(self, fstring, xcntr, ycntr, alpha_j, delta_j, 
                 SexMagAuto, SexMagAutoErr, SexTargets, SexSky,
                 flag, SexHalfRad, mag_zero,
                 C, C_err, A, A_err, S, S_err, G, M, 
                 components, decompose, repeat, detail, final_result_file,
                 H0, WM, WV, pixelscale, pymorph_config):

        self.chi2sq = 1.0
        self.PYMORPH_PATH = os.path.dirname(__file__)

        self.fstring = fstring
        self.xcntr = xcntr
        self.ycntr = ycntr
        self.alpha_j = alpha_j
        self.delta_j = delta_j

        self.SexMagAuto = SexMagAuto
        self.SexMagAutoErr = SexMagAutoErr
        self.SexTargets = SexTargets
        self.SexSky = SexSky
        self.flag = flag
        self.SexHalfRad = SexHalfRad
        self.mag_zero = mag_zero

        self.C  = C
        self.C_err = C_err
        self.A = A
        self.A_err = A_err
        self.S = S
        self.S_err = S_err
        self.G = G
        self.M = M

        self.components = components 
        self.decompose = decompose
        self.repeat = repeat
        self.detail = detail

        self.final_result_file = final_result_file

        self.H0 = H0
        self.WM = WM
        self.WV = WV
        self.pixelscale = pixelscale
        self.pymorph_config = pymorph_config

    def writeparams(self, params_to_write, distance_psf_gal, z):

        bkg1 = "#E6E6FA"
        bkg2 = "#D3D3D3"

        # this dictionary will hold any parameters that may be printed
        all_params = dict((params_to_write[key][0], -999) for key in params_to_write.keys())
        x = datetime.date.today()
        all_params['Date'] = f'{x.year}-{x.month}-{x.day}' 
        # load some of the passed in parameters
        all_params['Name'] = self.fstring
        all_params['ra_gal'] = self.alpha_j
        all_params['dec_gal'] = self.delta_j
        all_params['z'] = z
        all_params['C'] = self.C
        all_params['C_err'] = self.C_err 
        all_params['A'] = self.A
        all_params['A_err'] = self.A_err
        all_params['S'] = self.S
        all_params['S_err'] = self.S_err
        all_params['G'] = self.G
        all_params['M'] = self.M
        all_params['distance_psf_gal'] = round(distance_psf_gal, 3)
        all_params['mag_auto'] = self.SexMagAuto
        all_params['magerr_auto'] = self.SexMagAutoErr
        all_params['num_targets'] = self.SexTargets
        all_params['SexSky'] = self.SexSky
        all_params['flag'] = self.flag
        all_params['SexHalfRad'] = self.SexHalfRad
        all_params['magzp'] = self.mag_zero

        #print(all_params['C'])
        # Now continue
        comp = self.components

        f_tpl = open(str(self.PYMORPH_PATH) + '/html/default.html', 'r')
        template = f_tpl.read()
        f_tpl.close()

        cutimage = 'I{}.fits'.format(self.fstring)
        ra1, ra2, ra3 = RaDegToHMS(self.alpha_j)
        dec1, dec2, dec3 = DecDegToDMS(self.delta_j)

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
        if(self.repeat == False or self.repeat):
            for line_i in fileinput.input("index.html", inplace=1):
                line_i = line_i.strip()
                if not '</BODY></HTML>' in line_i:
                    print(line_i)
            indexfile = open('index.html', 'a+')
            NoImage = 1
            for indexline in indexfile:
                if self.fstring in indexline:
                    NoImage = 0
                else:
                    pass
            if NoImage:
                indexfile.writelines(['<a href="R_',\
                                      self.fstring,'.html',\
                                      '"> ', self.fstring,\
                                      ' </a> <br>\n'])
            indexfile.writelines(['</BODY></HTML>\n'])
            indexfile.close()

        # Reading fit.log
        if 'bar' in comp:
            basic_info, fit_info, measured_error_bad = read_fitlog(filename='fit.log', 
                                                      yes_bar=1)
        else:
            basic_info, fit_info, measured_error_bad = read_fitlog(filename='fit.log', 
                                                      yes_bar=0)

        #print('basic_info', basic_info)
        if basic_info['dof'] != 9999:
            try:
                from scipy import stats
                #sf=1-cdf which gives the probability of large chi2 given dof
                chisqprob_upper = lambda chisq, df: stats.chi2.sf(chisq, df)
                chi2 = basic_info['chi2nu'] * basic_info['dof']
                goodness_upper = chisqprob_upper(chi2, basic_info['dof']) 
                #cdf which gives the probability of sum upto a chi2 given dof
                goodness_lower = stats.chi2.cdf(chi2, basic_info['dof']) 
                goodness_upper = float(f'{goodness_upper:.2E}')
                goodness_lower = float(f'{goodness_lower:.2E}')
            except:
                goodness_upper = 9999
                goodness_lower
        else:
            goodness_upper = 9999
            goodness_lower = 9999

        all_params['goodness_upper'] = goodness_upper
        all_params['goodness_lower'] = goodness_lower
        #print('goodness_upper', goodness_upper)
        #print('goodness_lower', goodness_lower)

        #print('fit_info', fit_info)
        #print('measured_error_bad', measured_error_bad, self.flag)
        if measured_error_bad:
            try:
                #set the flag
                #print(1, 'all_params', all_params['flag'])
                all_params['flag'] = SetFlag(all_params['flag'], 
                                             GetFlag('ERRORS_FAILED'))
                all_params['flag'] = SetFlag(all_params['flag'], 
                                             GetFlag('GALFIT_FAIL'))

                #print(2, all_params['flag'])
            except:# badflag:
                # the flag is already set
                pass
            
 
        if 'Input' in basic_info:
            alpha_ned = str(self.alpha_j)[:10]
            delta_ned = str(self.delta_j)[:10]
        else:
            alpha_ned = ''
            delta_ned = ''
            
        if self.decompose:
            initial_conf = basic_info['initial_conf']
            restart_conf = basic_info['restart_conf']
            #print(restart_conf)
            # move the restart file to a reasonably named output file
            new_outname = initial_conf.replace('in','out')
            try:
                os.rename(restart_conf, new_outname)
            except:
                print("Failed to find restart file!! Galfit may have crashed!!")
            basic_info['restart_conf'] = new_outname

            all_params['chi2nu'] = basic_info['chi2nu']
            all_params['dof'] = basic_info['dof']
            # first check all err components and replace if nan or inf

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
                phy_parms = CosmoCal.cal(z, self.H0, 
                                         self.WM, self.WV, self.pixelscale)
                all_params['dis_modu'] = phy_parms[2]
                if 'bulge' in comp:
                    if all_params['re_pix'] != 9999:
                        all_params['re_kpc'] = phy_parms[3] * all_params['re_pix']
                        all_params['re_kpc_err'] = phy_parms[3] * all_params['re_pix_err']
                    else:
                        all_params['re_kpc'] = 9999
                        all_params['re_kpc_err'] = 9999
                if 'disk' in comp:
                    if all_params['rd_pix'] != 9999:
                        all_params['rd_kpc'] = phy_parms[3] * all_params['rd_pix']
                        all_params['rd_kpc_err'] = phy_parms[3] * all_params['rd_pix_err']
                    else:
                        all_params['rd_kpc'] = 9999
                        all_params['rd_kpc_err'] = 9999
                if 'bar' in comp:
                    if all_params['rbar_pix'] != 9999:
                        all_params['rbar_kpc'] = phy_parms[3] * all_params['rbar_pix']
                        all_params['rbar_kpc_err'] = phy_parms[3] * all_params['rbar_pix_err']
                    else:
                        all_params['rbar_kpc'] = 9999
                        all_params['rbar_kpc_err'] = 9999
                if 'point' in comp:
                    all_params['Pfwhm_kpc'] = 0.5 * phy_parms[3]
            # Finding derived parameters
            if 'bulge' in comp and 'disk' in comp:
                fb = 10**(-0.4 * (all_params['Ie'] - self.mag_zero))
                fd = 10**(-0.4 * (all_params['Id'] - self.mag_zero))
                if 'point' in comp:
                    fp = 10**(-0.4 * (all_params['Ip'] - self.mag_zero))
                else:
                    fp = 0.0
                if 'bar' in comp:
                    fbar = 10**(-0.4 * (all_params['Ibar'] - self.mag_zero)) 
                else:
                    fbar = 0.0

                try:
                    all_params['BD'] = round(fb / fd, 2)
                    all_params['BT'] = round(fb / (fb + fd + fp + fbar), 2)
                except:
                    all_params['BD'] = 9999
                    all_params['BT'] = 9999

            elif 'bulge' in comp:
                all_params['BT'] = 1.0
            elif 'disk' in comp:
                all_params['BD'] = 0.0
                all_params['BT'] = 0.0
            # Start writing html file. Now the template keywords will get values
            pngfile = 'P_{}.png'.format(self.fstring)
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
                                           '{}><TD> neighbor sersic'.format(bkg2) + \
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
                                           'bgcolor={}> <TD>'.format(bkg1) + ' ' + \
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
                                        'bgcolor={}>'.format(bkg2) +\
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
                                            'bgcolor={}>'.format(bkg1) + \
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
                        Object_Exp = '<TR align="center" bgcolor={}>'.format(bkg2) +\
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
                                         'bgcolor={}>'.format(bkg1) + \
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
                        Point_Vals = '<TR align="center" bgcolor={}>'.format(bkg2) + \
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
                                         'bgcolor={}>'.format(bkg1) + \
                                         '<TD>' + ' ' + '</TD> <TD>' + \
                                         str(fit_info[key]['xctr'][1]) + '</TD> <TD> '\
                                         + str(fit_info[key]['yctr'][1]) + \
                                     ' </TD> <TD> ' + str(fit_info[key]['mag'][1]) +\
                                     ' </TD> <TD> ' + str('9999') +  ' </TD> <TD> ' + \
                                     str('9999') + ' </TD> <TD> ' + \
                                     ' ' + ' </TD> <TD> ' + str('9999') +  \
                                     ' </TD> <TD> ' + str('9999') + \
                                     ' </TD> <TD> ' + str('9999') + ' </TD></TR>'
            except Exception as inst:
                print(type(inst))     # the exception instance
                print(inst.args)      # arguments stored in\
                                                     # .args
                print(inst)           # __str__ allows args\
                                                     # to printed directly
                print("something bad happened (writefunc writing)!!!!\n\n")
                print(traceback.print_exc())                   
            if 'bulge' in comp:
                try:
                      pixelscale = self.pixelscale
                except:
                      pixelscale = 1
                try:
                    all_params['AvgIe'] = all_params['Ie'] + 2.5 * np.log10(2 * 3.14 * pixelscale * \
                                          pixelscale *  all_params['re_pix'] *  all_params['re_pix'] * np.sqrt(1 -  all_params['eb']**2.0))
                    all_params['AvgIe'] = round(all_params['AvgIe'], 2)
                    AvgMagInsideReErr2 = (1.085 * np.sqrt((2 *  all_params['re_pix'] *  all_params['re_pix_err'])**2.0 + \
                                                         (( all_params['eb'] *  all_params['eb_err']) / \
                                                          np.sqrt(1 -  all_params['eb']**2.0))**2.0)) / \
                                                          (np.sqrt(1 -  all_params['eb']**2.0) * 2 * 3.14 * \
                                                           all_params['re_pix'] *  all_params['re_pix'])
                    all_params['AvgIe_err'] = round(np.sqrt( all_params['Ie_err']**2.0 + AvgMagInsideReErr2**2.0), 2)
                except OverflowError:
                    all_params['AvgIe'] = np.inf
                    all_params['AvgIe_err'] = np.inf
                for key in ['AvgIe', 'AvgIe_err']:
                    if np.isnan(all_params[key]) or np.isinf(all_params[key]):
                        if np.isnan(all_params[key]):
                            all_params[key] = -9999.99
                        else:
                            all_params[key] = -6666.66
                        try:
                            #set the AVGIe flag
                            all_params['flag'] = SetFlag(all_params['flag'], GetFlag('AVGIE_FAILED'))
                        except:# badflag:
                            # the flag is already set
                            pass

            wC = str(all_params['C'])[:5]
            wA = str(all_params['A'])[:5]
            wS = str(all_params['S'])[:5]
            wG = str(all_params['G'])[:5]
            wM = str(all_params['M'])[:5]
            wBD = str(all_params['BD'])[:5]
            wBT = str(all_params['BT'])[:5]
            wAvgMagInsideRe = str( all_params['AvgIe'])[:5]
            error_mesg1 = ''
            error_mesg2 = ''
            error_mesg3 = ''
            error_mesg4 = ''
            error_mesg5 = ''
            error_mesg6 = ''
            error_mesg7 = ''
            if 0:#self.starthandle:
                error_mesg6 = '<a href="R_' + self.fstring + \
                              '_1.html"> Crashed </a>' 

            # Now test for fitting problems and set flags for analysis
            all_params['FitFlag'] = 0
            HitLimit = 0

            self.UMag = -10
            self.LMag = 50
            self.LRe = 0.
            self.URe = 200.
            self.LN = 0.2
            self.UN = 20
            self.LRd = 0.
            self.URd = 200.
            self.center_deviation_limit = 2
            chi2nu = all_params['chi2nu']
            if not self.detail:
                if 'bulge' in comp:
                    if abs(all_params['Ie'] - self.UMag) < 0.2 or abs(all_params['Ie'] - self.LMag) < 0.2:
                        all_params['FitFlag'] = SetFlag(all_params['FitFlag'],Get_FitFlag('IE_AT_LIMIT'))
                    if abs(all_params['re_pix'] - self.LRe) < 0.1 or abs(all_params['re_pix'] - self.URe) < 1.0:
                        all_params['FitFlag'] = SetFlag(all_params['FitFlag'],Get_FitFlag('RE_AT_LIMIT'))
                    if abs(all_params['n'] - self.LN) < 0.03 or abs(all_params['n'] - self.UN) < 0.5:
                        all_params['FitFlag'] = SetFlag(all_params['FitFlag'],Get_FitFlag('N_AT_LIMIT'))
                    if abs(all_params['eb'] - 0.0) < 0.05 or abs(all_params['eb'] - 0.0) > 0.95:
                        all_params['FitFlag'] = SetFlag(all_params['FitFlag'],Get_FitFlag('EB_AT_LIMIT'))
                if 'disk' in comp:
                    if abs(all_params['Id'] - self.UMag) < 0.2 or abs(all_params['Id'] - self.LMag) < 0.2:
                        all_params['FitFlag'] = SetFlag(all_params['FitFlag'],Get_FitFlag('ID_AT_LIMIT'))
                    if abs(all_params['rd_pix'] - self.LRd) < 0.1 or abs(all_params['rd_pix'] - self.URd) < 1.0:
                        all_params['FitFlag'] = SetFlag(all_params['FitFlag'],Get_FitFlag('RD_AT_LIMIT'))
                    if abs(all_params['ed'] - 0.0) < 0.05 or abs(all_params['ed'] - 0.0) > 0.95:
                        all_params['FitFlag'] = SetFlag(all_params['FitFlag'],Get_FitFlag('ED_AT_LIMIT'))

            if not all_params['FitFlag']:
                error_mesg4 = str(error_mesg4) + 'One of the parameters'
                error_mesg5 = str(error_mesg5) + '          hits limit!'

            sigma_percentage = {1:0.682689, 2:0.954499, 3:0.997300, 4:0.999936,
                                5:0.999999}
            sigma_percentage = {1: 0.317311, 2: 0.0455, 3: 0.0027, 4: 6.399e-5, 5: 1e-6}

            print('sigma_percentage', sigma_percentage[self.pymorph_config['goodness_limit']])
            #3sigma gaussian has 0.997 probability. Then probabilty 0.003 shows than the model is not agreeing with the data. This means it is significant  
            if all_params['goodness_upper'] < sigma_percentage[self.pymorph_config['goodness_limit']]: 
                error_mesg2 = str(error_mesg2) + 'Upper Goodness is poor!'
                all_params['FitFlag'] = SetFlag(all_params['FitFlag'], Get_FitFlag('SMALL_UPPER_GOODNESS'))
                error_mesg1 = str(error_mesg1) + 'Chi2nu is large!'

            
            if all_params['goodness_lower'] < sigma_percentage[self.pymorph_config['goodness_limit']]: 
                error_mesg2 = str(error_mesg2) + 'Lower Goodness is poor!'
                all_params['FitFlag'] = SetFlag(all_params['FitFlag'], Get_FitFlag('SMALL_LOWER_GOODNESS'))
                error_mesg1 = str(error_mesg1) + 'Chi2nu is low!'

            if all_params['chi2nu'] > self.pymorph_config['chi2nu_limit']:
                all_params['FitFlag'] = SetFlag(all_params['FitFlag'],Get_FitFlag('LARGE_CHISQ'))

            if abs(all_params['bulge_xctr'] - self.xcntr) > self.center_deviation_limit or \
                   abs(all_params['bulge_yctr'] - self.ycntr) > self.center_deviation_limit or \
                   abs(all_params['disk_xctr'] - self.xcntr) > self.center_deviation_limit or \
                   abs(all_params['disk_yctr'] - self.ycntr) > self.center_deviation_limit:
                if all_params['bulge_xctr'] == -999 or all_params['disk_yctr'] == -999:
                    pass
                else:
                    all_params['FitFlag'] = SetFlag(all_params['FitFlag'], Get_FitFlag('FAKE_CNTR'))
                    error_mesg3 = str(error_mesg3) + 'Fake Center!'

            if all_params['FitFlag'] > 0:
                img_notify = str(self.PYMORPH_PATH) + '/html/goodfit.gif'
                good_fit = 1
            else:
                img_notify = str(self.PYMORPH_PATH) + '/html/badfit.gif'
                good_fit = 0

            # This can be a waste of time if the list is wrong...
            # Finding number of runs in the csv file 
            all_params['run'] = 1
            if os.path.exists(self.final_result_file):
                for line_res in csv.reader(open(self.final_result_file).readlines()[1:]):    
                    if(str(line_res[0]) == self.fstring):
                        all_params['run'] += 1
            #if self.GalSky != 9999:
            #    all_params['GalSky'] = self.GalSky

        #print('params_to_write', params_to_write)
        #print('all_params', all_params)

        #for key in params_to_write.keys():
        #    print(key, all_params[params_to_write[key][0]])

        #Final flag and FitFlag
        self.flag = all_params['flag']
        self.fit_flag = all_params['FitFlag']
        # Writing csv file 
        #print('self.final_result_file', self.final_result_file)
        f_res = open(self.final_result_file, "a")
        writer = csv.writer(f_res)
        writer.writerow([all_params[params_to_write[key][0]] for key in params_to_write.keys()])
        f_res.close()
        #sys.exit()
        # Writing data base
        #print('params_to_write', params_to_write)
        #print('all_params', all_params)
        try:
            WDB = WriteDB()
            WDB._first_db(params_to_write, all_params, self.pymorph_config)
        except Exception as e:
            print('No database can be created!')
            print(e)

        
        
        #print(error_mesg1, error_mesg2, error_mesg3, error_mesg4, error_mesg5, error_mesg6, error_mesg7)
        if self.decompose:
            #writing html
            outfile = open('R_{}.html'.format(self.fstring),'w')
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
            outfile = 'R_{}.html'.format(self.fstring)
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
    is_comp = False # for tracking problems with the err params on galfit 
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

    #print('obj', obj)
    # now replace any nan or inf parameters
    for key in obj.keys():
        if np.isnan(obj[key][1]):
            obj[key][1] = -9999.99
            is_comp = True
        elif np.isinf(obj[key][1]):
            obj[key][1] = -6666.66
            is_comp = True
    return obj, is_comp

def read_fitlog(filename = 'fit.log', yes_bar = 0):
    """ This function will read the fit log and return all the relevant
    information in 2 Dictionaries, 1 with the basic info and one with
    the fit info"""

    measured_error_bad = False
    
    neighbor = 0
    basic_info = {}
    fit_info = {}
    if os.path.exists(filename):
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
                elif(str(line[0]) == 'Chi^2'):
                    basic_info['dof'] = int(line[5]) 
                    #print('basic_info', basic_info)
                    #sys.exit()
                elif(str(line[0]) == 'Output'):
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
                    #print('err_line', err_line)
                    #print('line', line) 
                    fit_info[key], is_comp= load_component(line, err_line)    
                    measured_error_bad += is_comp
                    #print('measured_error_bad, is_comp', measured_error_bad, is_comp)
            except: 
                pass
            
    else:
        print("{} File in {} does not exist!!!!".format(filename, os.getcwd()))
    
    return basic_info, fit_info, measured_error_bad

if __name__ == '__main__':
    
    params_to_write = output_params(dbparams=None, decompose=True)

    WF = WriteHtmlCSV('cl1358_1.0', 13.540899999999993, 13.9297, 221.8572889, 8.4680823, 22.9636, 0.0319, 1, 0.656651, 6, 2.546, 25.256, 9999, 9999, 9999, 9999, 9999, 9999, 9999, 9999, ['bulge', 'disk'], True, False, False, 71.0, 0.27, 0.73, 0.045)
    WF.writeparams(params_to_write, 65, 9999, 0.6)
