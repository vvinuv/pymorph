from os.path import exists
import csv
import fileinput
from cosmocal import *
import config as c

class WriteHtmlFunc:
    """The class which will write html and csv output. This class will also 
       check whether the fit is good or bad using the Chisq and Goodness value
       It will also notify the goodness/badness of fit"""
    def __init__(self, cutimage, xcntr, ycntr, distance, alpha_j, delta_j, z, Goodness, C, C_err, A, A_err, S, S_err, G, M, Flag):
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
        self.Flag         = Flag
        self.write_params = write_params(cutimage, xcntr, ycntr, distance, alpha_j, delta_j, z, Goodness, C, C_err, A, A_err, S, S_err, G, M, Flag)

def write_params(cutimage, xcntr, ycntr, distance, alpha_j, delta_j, z, Goodness, C, C_err, A, A_err, S, S_err, G, M, Flag):
    Goodness = float(str(round(Goodness, 3))[:5])
    f_tpl = open(str(c.PYMORPH_PATH) + '/default.html', 'r')
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
    if(c.repeat == False):
        for line_i in fileinput.input("index.html",inplace =1):
            line_i = line_i.strip()
            if not '</BODY></HTML>' in line_i:
                print line_i    
        indexfile = open('index.html', 'a+')
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
                    alpha_ned = str(alpha_j)[:6]
                    delta_ned = str(delta_j)[:6]
                if(str(values[0]) == 'Init.'):
                    initial_conf = str(values[4])
                if(str(values[0]) == 'Restart'):
                    restart_conf = str(values[3])
                if(str(values[0]) == 'Chi^2/nu'):
                    chi2nu = float(values[2])
                    Distance = str(round(distance, 3))[:5]
                if(str(values[0]) == 'sersic' and object == 1):
                    mag_b = float(values[4])
                    re = float(values[5])
                    n = float(values[6])
                    object += 1
                if(str(values[0]) == 'expdisk'):
                    mag_d = float(values[4])
                    rd = float(values[5])
                if(str(values[0])[:1] == '('):
                    if(str(a) == 'sersic' and object_err == 1):
                        mag_b_err = float(values[2])
                        re_err  = float(values[3])
                        n_err = float(values[4])
                        object_err += 1
                    if(str(a) == 'expdisk'):
                        mag_d_err = float(values[2])
                        rd_err = float(values[3])
                a=values[0]				
            except:
                pass
        if(z == 9999):
            re_kpc = 9999
            re_err_kpc = 9999
            rd_kpc = 9999
            rd_err_kpc = 9999
        else:
            phy_parms = cal(z, c.H0, c.WM, c.WV, c.pixelscale)
            re_kpc = phy_parms[3] * re
            re_err_kpc = phy_parms[3] * re_err
            rd_kpc = phy_parms[3] * rd
            rd_err_kpc = phy_parms[3] * rd_err
        BD = 10**(-0.4 * ( mag_b - mag_d))
        BT = 1.0 / (1.0 + 1.0 / BD)
    pngfile = 'P_' + str(cutimage)[:-4] + 'png'
    object = 1
    object_err = 1
    Neighbour_Sersic = ''
    if exists('fit.log'):
        for line in open('fit.log','r'): 
            values = line.split() 
            try: 
                if(str(values[0]) == 'sersic'):
                    if object > 1:
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
                    if object == 1:
                        bulge_xcntr = float(str(values[2])[1:-1])
                        bulge_ycntr = float(str(values[3])[:-1])
                        Object_Sersic = '<TD>' + str(values[0]) + '</TD> <TD> '\
                                        + str(values[2])[1:-1] + ' </TD> <TD> '\
                                        + str(values[3])[:-1] + ' </TD> <TD> ' \
                                        + str(values[4]) + ' </TD> <TD> ' \
                                        + str(values[5]) +  ' </TD> <TD> ' + \
                                        str(round(re_kpc, 3))[:5] +'</TD> <TD>'\
                                        + str(values[6]) + ' </TD> <TD> ' \
                                        + str(values[7]) + ' </TD> <TD> ' \
                                        + str(values[8]) + ' </TD> <TD> ' \
                                        + str(values[9]) + ' </TD>'
                        object += 1
                if(str(values[0]) == 'expdisk'):
                    disk_xcntr = float(str(values[2])[1:-1])
                    disk_ycntr = float(str(values[3])[:-1])
                    Object_Exp = '<TD>' + str(values[0]) + ' </TD> <TD> ' + \
                                 str(values[2])[1:-1] + ' </TD> <TD> ' + \
                                 str(values[3])[:-1] + ' </TD> <TD> ' + \
                                 str(values[4]) + ' </TD> <TD> ' + \
                                 str(values[5]) +  ' </TD> <TD> ' + \
                                 str(round(rd_kpc, 3))[:5] + ' </TD> <TD> ' + \
                                 ' ' + ' </TD> <TD> ' + str(values[6]) +  \
                                 ' </TD> <TD> ' + str(values[7]) + \
                                 ' </TD> <TD> ' + str(values[8]) + ' </TD>'
                if(str(values[0])[:1] == '('):
                    if(str(a) == 'sersic' and object_err > 1):
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
                    if(str(a) == 'sersic' and object_err == 1):
                        Object_Sersic_err = '<TD>' + ' ' + '</TD> <TD>' + \
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
                                            ' </TD>'
                        object_err += 1
                    if(str(a) == 'expdisk'):
                        Object_Exp_err = '<TD>' + ' ' + '</TD> <TD>' + \
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
                                         ' </TD>'
                a=values[0]				
            except:
                pass
    wC = str(C)[:5]
    wA = str(A)[:5]
    wS = str(S)[:5]
    wG = str(G)[:5]
    wM = str(M)[:5]
    wBD = str(BD)[:5]
    wBT = str(BT)[:5]
    error_mesg1 = ''
    error_mesg2 = ''
    error_mesg3 = ''
    if chi2nu <= c.chi2sq and Goodness >= c.Goodness and \
        abs(bulge_xcntr - xcntr) <= c.center_deviation and \
        abs(bulge_ycntr - ycntr) <= c.center_deviation and \
        abs(disk_xcntr - xcntr) <= c.center_deviation and \
        abs(disk_ycntr - ycntr) <= c.center_deviation:
        img_notify = str(c.PYMORPH_PATH) + '/goodfit.gif'
        good_fit = 1
    else:
        if chi2nu > c.chi2sq:
            error_mesg1 = str(error_mesg1) + 'Chi2nu is large!'
            Flag = Flag + 2048
        if Goodness < c.Goodness:
            error_mesg2 = str(error_mesg2) + 'Goodness is poor!'
            Flag = Flag + 4096
        if abs(bulge_xcntr - xcntr) > c.center_deviation or \
             abs(bulge_ycntr - ycntr) > c.center_deviation or \
             abs(disk_xcntr - xcntr) > c.center_deviation or \
             abs(disk_ycntr - ycntr) > c.center_deviation:
            error_mesg3 = str(error_mesg3) + 'Fake Center!'
            Flag = Flag + 8192
        img_notify = str(c.PYMORPH_PATH) + '/badfit.gif'
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
    writer.writerow([galid, alpha_j, delta_j, z, mag_b, mag_b_err, re, re_err, re_kpc, re_err_kpc, n, n_err, mag_d, mag_d_err, rd, rd_err, rd_kpc, rd_err_kpc, BD, BT, chi2nu, Goodness, run, C, C_err, A, A_err, S, S_err, G, M, distance, good_fit, Flag])
    f_res.close()

    outfile1 = open('galfit.html', 'w')
    outfile1.write('<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01' \
                   ' Transitional//EN">')
    outfile1.write('<HTML><HEAD> \n')
    outfile1.writelines(['<META HTTP-EQUIV="Refresh" CONTENT="10; \
                    URL=galfit.html"> \n'])
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

