import numpy as np
import pyfits


def selectpsf(ImG, CaT):
    c.psff = []
    im1 = pyfits.open(ImG)
    image = im1[0].data
    im1.close()
    AreaOfObj = c.AreaOfObj
    def FindPsf(AreaOfObj, CaT):
        for line in open(CaT,'r'):
            values = line.split()
            try:
                BaKgR  = float(values[10])
                if float(values[16]) >= c.StarGalProb and \
                   float(values[13]) > AreaOfObj \
                   and float(values[14]) < 50.0:
                    xcntr = float(values[1]) - 1
                    ycntr = float(values[2]) - 1
                    #size of the psf is 8 times the sigma assuming the star has 
                    #Gaussian profile
                    PsfSize = np.floor(float(values[14])) * c.starsize
                    x1 = int(xcntr) + (PsfSize/2)
                    x2 = int(xcntr) - (PsfSize/2)
                    y1 = int(ycntr) + (PsfSize/2)
                    y2 = int(ycntr) - (PsfSize/2)
                    ra1 = int(float(values[3]) / 15.0)
                    ra2 = int((float(values[3]) / 15.0 - int(float(values[3])\
                          / 15.0))*60.0)
                    ra3 = (((float(values[3]) / 15.0 - int(float(values[3]) / \
                          15.0))*60.0) - ra2) * 60.0
                    dec1 = int(float(values[4]))
                    dec2 = abs(int((float(values[4]) - dec1) * 60.0))
                    dec3 = (abs(((float(values[4]) - dec1) * 60.0)) - dec2) \
                           * 60.0
                    if ra1 < 10:
                        ra11 = '0' + str(ra1)
                    else:
                        ra11 = str(ra1)
                    if ra2 < 10:
                        ra22 = '0' + str(ra2)
                    else:
                        ra22 = str(ra2)
                    if ra3 < 10:
                        ra33 = '0' + (str(np.round(ra3, 1))[:3]).split('.')[0]\
                                      + \
                                     (str(np.round(ra3, 1))[:3]).split('.')[1]
                    else:
                        ra33 = (str(np.round(ra3, 1))[:4]).split('.')[0] + \
                               (str(np.round(ra3, 1))[:4]).split('.')[1]
                    if abs(dec1) < 10:
                        if dec1 < 0.0:
                            dec11 = '-0' + str(abs(dec1))
                        else:
                            dec11 = '+0' + str(abs(dec1))
                    else:
                        if dec1 < 0.0:
                            dec11 = str(dec1)
                        else:
                            dec11 = '+' + str(dec1)
                    if dec2 < 10:
                        dec22 = '0' + str(dec2)
                    else:
                        dec22 = str(dec2)
                    if dec3 < 10:                        dec33 = '0' + (str(np.round(dec3, 1))[:3]).split('.')[0]\                                
                                 + (str(np.round(dec3, 1))[:3]).split('.')[1]
                    else:
                        dec33 = (str(np.round(dec3, 1))[:4]).split('.')[0] +\
                                (str(np.round(dec3, 1))[:4]).split('.')[1]                     psffile = 'psf_' + str(ra11) + str(ra22) + str(ra33) + str(dec11) +str(dec22) + str(dec33) + '.fits'
                    if psffile in c.psff:
                        pass
                    else:
                        c.psff.append(psffile)
                        psf = image[y2:y1, x2:x1]
                        psf = psf - BaKgR
                        if os.access(psffile, os.F_OK):
                            os.remove(psffile)
                        hdu = pyfits.PrimaryHDU(psf.astype(np.float32))
                        hdu.header.update('XCNTR', int(xcntr))
                        hdu.header.update('YCNTR', int(ycntr))
                        hdu.writeto(psffile)
            except:
                pass
    while len(c.psff) < 5 and AreaOfObj > 10:
        FindPsf(AreaOfObj, CaT)
        AreaOfObj -= 5
    if len(c.psff) == 0:
        manualpsf = raw_input("Unfortunately, NO psf is found. Please enter "\
                              "the psf name >>> ")
        c.psff.append(manualpsf)
    if c.Interactive:
        PsfList = []
        TmPLST = []
        for element in c.psff:
            TmPLST.append(element)
        print 'Checking Started. You can just visually check the psf. You can'\
               ' do the thorough checking later'
        for element in c.psff:
            os.system('xpaset -p ds9 frame clear')
            if exists(element):
                os.system('cat ' + str(element) + ' | xpaset ds9 fits')
                os.system('xpaset -p ds9 scale mode zscale')
                os.system('xpaset -p ds9 zoom to fit')
            else:
                print 'The psf you have given does NOT exists!!!'
                pass
            write = raw_input("Do you need this psf? ('y' if yes, 'c'"\
                              " to cancel psf checking) " )
            if write == 'y':
                PsfList.append(element)
                TmPLST.remove(element)
#                c.psff.remove(element)
            elif write == 'c':
                for element1 in TmPLST:
                    try:
                        os.remove(element1)
                    except:
                        pass
                break
            else:
                try:
                    os.remove(element)
                except:
                    pass
        print '\nFinal Checking Started. If you are using psfselect = 2, be '\
              'carefull!!! This is your last chance for selecting psf. '\
              'Select ONLY GOOD psf. Press "y" to accept the previous psf. '\
             'ENTER to go to the next one. This will continue until you press '\
             '"1" when it asked Finished? ALL THE BEST! \n'
        UrPsfChk = raw_input("Do you want to use your own psf? Enter 'y' or " \
                             "'n' >>> ")
        if UrPsfChk == 'y':
            UrPsf = raw_input("Enter your psf >>> ")
            fff = open('psflist.list', 'ab')
            fff.writelines([str(UrPsf), '\n'])
            fff.close()
            if(c.galcut):
                c.ValueS.append(UrPsf)
                fwithpsf = open('CatWithPsf.cat', 'ab')
                for v in c.ValueS:
                    fwithpsf.writelines([str(v), ' '])
                fwithpsf.writelines(['\n'])
                fwithpsf.close()
            finish = 1
            for element in PsfList:
                if os.access(element, os.F_OK):
                    os.remove(element)
        else:
            finish = 0
            TmpPsfList = []
        UpdateCounter = 1
        while finish == 0:
            for element in PsfList:
                if exists(element):
                    os.system('xpaset -p ds9 frame clear')
                    os.system('cat ' + str(element) + ' | xpaset ds9 fits')
                    os.system('xpaset -p ds9 scale mode zscale')
                    os.system('xpaset -p ds9 zoom to fit')
                else:
                    print 'The psf you have given is NOT exists!!!'
                    pass
                write = raw_input("Do you REALLY need this psf? ('y' or," \
                                  "'n' or press any key to continue) ")
                if write == 'y':
                    TmpPsfList.append(element)
                    fff = open('psflist.list', 'ab')
                    fff.writelines([str(element), '\n'])
                    if(c.galcut):
                        if UpdateCounter:
                            c.ValueS.append(element)
                            fwithpsf = open('CatWithPsf.cat', 'ab')
                            for v in c.ValueS:
                                fwithpsf.writelines([str(v), ' '])
                            fwithpsf.writelines(['\n'])
                            fwithpsf.close()
                            UpdateCounter = 0
                    fff.close()
                if write == 'n':
                    TmpPsfList.append(element)
                    try:
                        os.remove(element)
                    except:
                        pass
                else:
                    pass
            for element in TmpPsfList:
                try:
                    PsfList.remove(element)
                except:
                    pass
            TmpPsfList = []
            fi = raw_input("Finished? ('1' to finish, any other key to continue) ")
            if fi == '0' or fi == '1':
                finish = int(fi)
                if finish == 1:
                    for element in PsfList:
                        try:
                            os.remove(element)
                        except:
                            pass
            else:
                finish = 0
    else:
        for element in c.psff:
            fff = open('psflist.list', 'ab')
            fff.writelines([str(element), '\n'])
        fff.close()


def runpsfselect():
    if(c.galcut):   #Given galaxy cutouts
        obj_file = open(os.path.join(c.datadir, c.clus_cata),'r')
        pnames = obj_file.readline().split()
        c.ValueS = []
        for v in pnames:
            c.ValueS.append(v)
        c.ValueS.append('star')
        fwithpsf = open('CatWithPsf.cat', 'ab')
        for v in c.ValueS:
            fwithpsf.writelines([str(v), ' '])
        fwithpsf.writelines(['\n'])
        fwithpsf.close()
        pdb = {}                        #The parameter dictionary
        for line_j in obj_file:
            try:
                values = line_j.split()
                c.ValueS = []
                for v in values:
                    c.ValueS.append(v)
                k = 0
                for pname in pnames:
                    pdb[pname] = values[k]
                    k += 1
                try:
                    gal_id = pdb["gal_id"]
                except:
                    try:
                        gal_id = pdb["gimg"][:-5]
                    except:
                        print "No image or gal_id found in the object" \
                              "catalogue. Exiting"
                        os._exit(0)
                try:
                    gimg = pdb["gimg"]    #Galaxy cutout
                except:
                    if exists(os.path.join(c.datadir, \
                               'I' + str(c.rootname) + '_' \
                               + str(gal_id) + '.fits')):
                        gimg = 'I' + str(c.rootname) + '_' \
                                + str(gal_id) + '.fits'
                    elif exists(os.path.join(c.datadir, \
                         str(gal_id) + '.fits')):
                        gimg = str(gal_id) + '.fits'
                    else:
                        print "No image found. Exiting"
                        os._exit(0)
                try:
                    wimg = pdb["wimg"]   #Weight cut
                except:
                    if exists(os.path.join(c.datadir, \
                              'W' + str(c.rootname) + '_' + \
                              str(gal_id) + '.fits')):
                        wimg = 'W' + str(c.rootname) + '_' + \
                               str(gal_id) + '.fits'
                    else:
                        wimg = 'None'
                GiMg = pyfits.open(os.path.join(c.datadir, gimg))
                headerGiMg = GiMg[0].header
                if (headerGiMg.has_key('GAIN')):
                    c.SEx_GAIN = headerGiMg['GAIN']
                else:
                    c.SEx_GAIN = 1
                GiMg.close()
                if exists(sex_cata):
                    pass
                else:
                    RunSex(os.path.join(c.datadir, gimg), os.path.join(c.datadir, wimg), 'None', 9999, 9999, 0)
                try:
                    selectpsf(os.path.join(c.datadir, gimg), sex_cata)
                except:
                    pass
                if os.access(sex_cata, os.F_OK):
                    os.remove(sex_cata)
            except:
                pass
        obj_file.close()
        AskForUpdate = raw_input("Do you want to update the clus_cata? " \
                        "('y' for yes) ")
        if AskForUpdate == 'y':
            cmd = 'mv CatWithPsf.cat ' +str(c.clus_cata)
            os.system(cmd)
        else:
            pass
    else:
        selectpsf(os.path.join(c.datadir, c.imagefile), sex_cata)


def FindAndFit():
    if c.findandfit == 1:
        if c.psfselect > 2:
            #magmin = raw_input("Enter Minimum Magnitude >>> ")
            magmin = raw_input("What is faintest magnitude galaxy that you want to fit? >>> ")     #modified by abhishek
            try:
                magmin = float(magmin) * 1.0
            except:
                magmin = 9999
            #magmax = raw_input("Enter Maximum Magnitude >>> ")
            magmax = raw_input("What is brightest magnitude galaxy that you want to fit? >>> ")    #modified by abhishek
            try:
                magmax = float(magmax) * 1.0
            except:
                magmax = -9999
            stargal = raw_input("Enter star-galaxy classification "\
                                "(1 for star and 0 is galaxy (0.8 is "\
                                " a good number) >>> ")
            redshift = raw_input("Enter redshift, if you know >>> ")
            try:
                redshift = float(redshift)*1.0
            except:
                redshift = 9999
        elif c.psfselect <= 2:
            try:
                magmin = c.maglim[0]
                magmax = c.maglim[1]
                stargal = c.stargal
                redshift = 9999
            except:
                magmin = 9999
                magmax = -9999
                stargal = 0.8
                redshift = 9999
        NewClusCata = open(c.clus_cata,'w')
        if redshift == 9999:
            NewClusCata.writelines(['gal_id ra1 dec1 mag\n'])
        else:
            NewClusCata.writelines(['gal_id ra1 dec1 mag z\n'])
        for lineSex in open(sex_cata, 'r'):
            ValueSex = lineSex.split()
            try:
                if float(ValueSex[17]) > magmax and \
                   float(ValueSex[17]) < magmin and \
                   float(ValueSex[16]) < float(stargal):
                    if redshift == 9999:
                        NewClusCata.writelines([str(ValueSex[0]), ' ', \
                                   str(float(ValueSex[3]) / 15.0), ' ', \
                                   str(float(ValueSex[4])), ' ',\
                                   str(ValueSex[17]), '\n'])
                    else:
                        NewClusCata.writelines([str(ValueSex[0]), ' ', \
                            str(float(ValueSex[3]) / 15.0), ' ', \
                            str(float(ValueSex[4])), \
                            ' ', str(redshift), ' ', \
                            str(ValueSex[17]), '\n'])
            except:
                pass
        NewClusCata.close()
        c.searchrad = '0.05arc'
    else:
        pass


