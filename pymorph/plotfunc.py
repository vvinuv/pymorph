import os
import fitsio
import numpy as np
from matplotlib.colors import Normalize
import matplotlib.pylab as plt
import numpy.ma as ma


class PlotFunc:

    """
        The class for plotting. It will plot the galaxy image, the model galaxy
        and the residual image. It also plot the histogram of residual image
        and compute the Goodness parameter out of it. It is defind as the ratio
        of number pixels within n times sky sigma around zero value to the total
        number of pixels. The writehtmlfunc.py will check whether the fit is
        good or bad according to some conditions. This will also plot the 1D
        profile of the galaxy and model galaxy from ellipse fitting. The other
        plot it gives is the mask image. After plotting it saves the image
        as png file with name P_string(galid).png

    """

    def __init__(self, oimg, mimg, fstring,     
                 sky, skysig, mag_zero, save_name=None):

        print('Plotting')
        self.oimg = oimg
        self.mimg = mimg
        self.fstring = fstring
        self.sky = sky
        self.skysig = skysig
        self.mag_zero = mag_zero
        self.save_name = save_name
        #print(oimg, mimg, fstring, sky, skysig)
        if self.save_name is None:
            self.save_name = 'P_{}.png'.format(self.fstring)


    def plot_profile(self):

        #print('P0')
        if not os.path.exists(self.oimg):
            print('No output image exists. GALFIT might have been failed')
            print('Exiting plotting')
        else:
            #print('P1')
            # Read the GALFIT output
            fits = fitsio.FITS(self.oimg)
            galaxy = fits[1].read()
            #print('P01')
            model = fits[2].read()
            #print('P02')
            residual = fits[3].read()
            #print('P03')
            fits.close()
            #Color normalization
            anorm = Normalize(self.sky - 2 * self.skysig , self.sky + 12.0 * self.skysig)

            #print('P2')
            # Read mask
            f_mask = fitsio.FITS(self.mimg)
            mask = f_mask[0].read().astype(int)
            f_mask.close()
            #print('P3')
            model = ma.array(model, mask=mask, fill_value=0)#9999)
            model = model.filled()

            #print('P4')
            maskedresidual = ma.masked_array(residual, mask)
            anormRes = Normalize(-2 * self.skysig, 3 * self.skysig)
            residual = ma.filled(maskedresidual, 0)#9999)

            #print('P5')
            #colorbar(shrink=0.90)
            valid_pixels = ma.count(maskedresidual)
            print('No of valid pixels >>> ', valid_pixels)

            pixels_in_skysig = residual[np.where(abs(residual) <= self.skysig)].size
            print('No of pixels within sky sigma >>> ', pixels_in_skysig)

            #goodness is the probability of the model being correct. It is identified using chi2 and dof. The probability is calculated using chisq distribution.
            #print('P6')
            NYPTS, NXPTS = galaxy.shape
            hist_mask = np.zeros((NXPTS, NYPTS))
            hist_mask[np.where(abs(residual) > 12 * self.skysig)] = 1

            #print('P7')
            hist_res = ma.masked_array(residual, mask=hist_mask)

        #print('P8')
        fell = 'E_{}.txt'.format(self.fstring)
        if os.path.exists(fell):
            data = np.genfromtxt(fell, delimiter=',', names=True)
            sma = data['sma']           #sma from ellise fitting
            flux = data['intensity'] #Flux at various sma
            flux_err =data['intensity_err']     #Error in Flux
            mag = data['mag'] + float(self.mag_zero) #Magnitude at various sma
            mag_uerr = data['magu']     #Upper error in magnitude
            mag_lerr = data['magl']     #lower error in Magnitude
            GalEll = 1
        else:
            GalEll = 0

        #print('P9')
        oell = 'OE_{}.txt'.format(self.fstring)
        if os.path.exists(oell):
            data_o = np.genfromtxt(oell, delimiter=',', names=True)
            sma_o = data_o['sma']               #sma from ellise fitting
            flux_o = data_o['intensity']     #Flux at various sma
            flux_err_o =data_o['intensity_err'] #Error in Flux
            mag_o = data_o['mag'] + float(self.mag_zero)#Magnitude at various sma
            mag_uerr_o = data_o['magu'] #Upper error in magnitude
            mag_lerr_o = data_o['magl'] #lower error in Magnitude
            ModelEll = 1
        else:
            ModelEll = 0

        #print('P10')
        if 1:
            try:
                MaxRad = sma.max()
            except:
                MaxRad = sma_o.max()
            #print('P11')
            NoOfPoints = int(30 * MaxRad / 50.)
            SmaCommon = np.logspace(0, np.log10(MaxRad), \
                                    NoOfPoints, endpoint=True)
            #print('P12')
            MagI = np.interp(SmaCommon, sma, mag)
            #print('P13')
            MagI_o = np.interp(SmaCommon, sma_o, mag_o)
            FluxI = np.interp(SmaCommon, sma, flux)
            FluxI_o = np.interp(SmaCommon, sma_o, flux_o)
            FluxErrI = np.interp(SmaCommon, sma, flux_err)
            FluxErrI_o = np.interp(SmaCommon, sma_o, flux_err_o)
            MagDev = MagI - MagI_o
            #print('P14')
            FluxErr = np.sqrt((FluxErrI / FluxI)**2.0 + (FluxErrI_o / FluxI_o)**2.0)
            #print('P15')
            MagLErr = (np.log10(FluxI / FluxI_o) - \
                            np.log10((FluxI / FluxI_o) - FluxErr)) * -2.5
            #print('P16')

            MagUErr = (np.log10((FluxI / FluxI_o) + FluxErr) -\
                           (np.log10(FluxI / FluxI_o))) * -2.5
        #except Exception as e:
        #    print(e)

        #Plotting Starts
        FigSize = [8.0, 4.6]
        MatPlotParams = {'axes.titlesize': 10,
                         'axes.labelsize': 10,
                         'xtick.labelsize': 5,
                         'ytick.labelsize': 5,
                         'figure.figsize': FigSize}
        plt.rcParams.update(MatPlotParams)

        rect1 = [0.125, 0.5, 0.225, 1.5 * 0.225]
        rect2 = [3.25 * 0.125, 0.5, 0.225, 1.5 * 0.225]
        rect3 = [5.5 * 0.125, 0.5, 0.225, 1.5 * 0.225]
        rect4 = [5.5 * 0.125, 0.075, 0.225, 1.5 * 0.225]
        rect5 = [3.25 * 0.125, 0.075, 0.225, 1.5 * 0.225]
    #    rect6 = [0.125, 0.75*0.225 + 0.075, 0.225, 1.5*0.225]
        rect7 = [0.125, 0.075, 0.225, 1.5 * 0.225]

        if 1:
            axUL = plt.axes(rect1)
            im = plt.imshow(galaxy, extent=[0, NXPTS, 0, NYPTS], 
                    cmap='Blues_r')#norm=anorm)
            Dx = abs(axUL.get_xlim()[0]-axUL.get_xlim()[1])
            Dy = abs(axUL.get_ylim()[0]-axUL.get_ylim()[1])
            plt.colorbar(shrink=0.9, format='%.2f')
            axUL.set_aspect(Dx/Dy)
            plt.title('Original Galaxy')

            axUM = plt.axes(rect2)
            im = plt.imshow(model, cmap='Blues_r',
                            extent=[0, NXPTS, 0, NYPTS])# norm=anorm)
            plt.colorbar(shrink=0.9, format='%.2f')
            Dx = abs(axUM.get_xlim()[0]-axUM.get_xlim()[1])
            Dy = abs(axUM.get_ylim()[0]-axUM.get_ylim()[1])
            axUM.set_aspect(Dx/Dy)
            plt.title('Model Galaxy + Mask')

            axUR = plt.axes(rect3)
            im = plt.imshow(residual, cmap='Blues_r',
                            extent=[0, NXPTS, 0, NYPTS])#, norm=anormRes)
            plt.colorbar(shrink=0.9)
            Dx = abs(axUR.get_xlim()[0]-axUR.get_xlim()[1])
            Dy = abs(axUR.get_ylim()[0]-axUR.get_ylim()[1])
            axUR.set_aspect(Dx/Dy)
            plt.title('Residual')

            axLR = plt.axes(rect4)
            hist_res1d = hist_res.compressed()
            nn, bins, patches = plt.hist(hist_res1d, 50, density=False)
            if (nn.argmax() < 16):
                arg_inc = nn.argmax()
            else:
                arg_inc = 16
            # print trapz(bins, nn)
            binmin = bins[nn.argmax() - arg_inc]
            binmax = bins[nn.argmax() + arg_inc]
            plt.axis([binmin, binmax, 0.0, max(nn)])
            plt.setp(patches, 'facecolor', 'g', 'alpha', 0.75)
            plt.xticks((binmin, binmin /2.0, 0.0, binmax/2.0, binmax),\
                    (str(binmin)[:5], str(binmin /2.0)[:5], str(0.0)[:3], \
                     str(binmax/2.0)[:5], str(binmax)[:5]))
            plt.grid(True)
            plt.title('Difference Histogram')
            Dx = abs(axLR.get_xlim()[0]-axLR.get_xlim()[1])
            Dy = abs(axLR.get_ylim()[0]-axLR.get_ylim()[1])
            axLR.set_aspect(Dx/Dy)

        #except Exception as e:
        #    print(e)


        if GalEll and ModelEll:

            axLL = plt.axes(rect7)

            ymin = min(min(mag), min(mag_o))
            ymax = max(max(mag), max(mag_o))
            xmax = max(max(sma), max(sma_o))

            plt.errorbar(sma, mag, [mag_uerr, mag_lerr], fmt='o',
                         ecolor='r', ms=3)
            plt.plot(sma_o, mag_o, color='g',lw=2)

            axLL.set_ylim(ymax, ymin)
            axLL.set_xlim(0, xmax)
            plt.xlabel('Radius')
            plt.ylabel('Surface Brightness')
            plt.title('1-D Profile Comparison')
            plt.grid(True)
            Dx = abs(axLL.get_xlim()[0]-axLL.get_xlim()[1])
            Dy = abs(axLL.get_ylim()[0]-axLL.get_ylim()[1])
            axLL.set_aspect(Dx/Dy)

            axLM = plt.axes(rect5)
            try:
                plt.errorbar(SmaCommon, MagDev, [MagUErr, MagLErr], 
                             fmt='o', ecolor='r', ms=3)
                plt.ylabel('Magnitude Deviation')
                plt.xlabel('Radius')
                plt.grid(True)
                Dx = abs(axLM.get_xlim()[0]-axLM.get_xlim()[1])
                Dy = abs(axLM.get_ylim()[0]-axLM.get_ylim()[1])
                axLM.set_aspect(Dx/Dy)
            except:
                pass
        elif GalEll:
            ymin = min(mag)
            ymax = max(mag)
            xmax = max(sma)
            axLL = plt.axes(rect7)
            plt.errorbar(sma, mag, [mag_uerr,mag_lerr], fmt='o',
                         ecolor='r', ms=3)
            axLL.set_ylim(ymax, ymin)
            axLL.set_xlim(0, xmax)
            plt.xlabel(r'Radius')
            plt.ylabel(r'Surface Brightness')
            plt.title('1-D Profile Comparison')
            plt.grid(True)
            Dx = abs(axLL.get_xlim()[0]-axLL.get_xlim()[1])
            Dy = abs(axLL.get_ylim()[0]-axLL.get_ylim()[1])
            axLL.set_aspect(Dx/Dy)
        elif ModelEll:
            ymin = min(mag_o)
            ymax = max(mag_o)
            xmax = max(sma_o)
            axLL = plt.axes(rect7)
            plt.plot(sma_o, mag_o, color='g',lw=2)
            axLL.set_ylim(ymax, ymin)
            axLL.set_xlim(0, xmax)
            plt.xlabel(r'Radius')
            plt.ylabel(r'Surface Brightness')
            plt.title('1-D Profile Comparison')
            plt.grid(True)
            Dx = abs(axLL.get_xlim()[0]-axLL.get_xlim()[1])
            Dy = abs(axLL.get_ylim()[0]-axLL.get_ylim()[1])
            axLL.set_aspect(Dx/Dy)

        plt.savefig(self.save_name)
#self.fstring = 'test_n5585_lR'
#PlotFunc('O_test_n5585_lR.fits', 'M_test_n5585_lR.fits', 82, 82, 1434., 40.)



