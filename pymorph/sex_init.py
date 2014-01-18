


def SExtractorConf():
    SEx_DETECT_MINAREA = raw_input("DETECT_MINAREA (6) >>> ")
    try:
        c.SEx_DETECT_MINAREA = float(SEx_DETECT_MINAREA)
        c.SEx_DETECT_MINAREA = int(c.SEx_DETECT_MINAREA)
    except:
        c.SEx_DETECT_MINAREA = 6
    SEx_DETECT_THRESH = raw_input('DETECT_THRESH (1.5) >>> ')
    try:
        c.SEx_DETECT_THRESH = float(SEx_DETECT_THRESH)
    except:
        c.SEx_DETECT_THRESH = 1.5
    SEx_ANALYSIS_THRESH = raw_input('ANALYSIS_THRESH (1.5) >>> ')
    try:
        c.SEx_ANALYSIS_THRESH = float(SEx_ANALYSIS_THRESH)
    except:
        c.SEx_ANALYSIS_THRESH = 1.5
    SEx_FILTER = raw_input('FILTER (Y/N) >>> ')
    while SEx_FILTER != 'Y' and SEx_FILTER != 'N' and SEx_FILTER != '':
        SEx_FILTER = raw_input('FILTER (Y/N) >>> ')
    if len(SEx_FILTER) == 0:
        c.SEx_FILTER = 'Y'
    else:
        c.SEx_FILTER = SEx_FILTER
    print 'Available options for convolve filter are gauss_1.5_3x3.conv(1) '\
          'gauss_2.0_3x3.conv(2) gauss_2.0_5x5.conv(3) gauss_2.5_5x5.conv(4) '\
          'gauss_3.0_5x5.conv(5) gauss_3.0_7x7.conv(6) gauss_4.0_7x7.conv(7) '\
          'gauss_5.0_9x9.conv(8) default(0)'
    SEx_FILTER_NAME  = raw_input('FILTER_NAME (default.conv) >>> ')
    if len(SEx_FILTER_NAME) == 0 or SEx_FILTER_NAME == '0':
        c.SEx_FILTER_NAME = 'default.conv'
    elif SEx_FILTER_NAME == '1':
        c.SEx_FILTER_NAME = 'gauss_1.5_3x3.conv'
    elif SEx_FILTER_NAME == '2':
        c.SEx_FILTER_NAME = 'gauss_2.0_3x3.conv'
    elif SEx_FILTER_NAME == '3':
        c.SEx_FILTER_NAME = 'gauss_2.0_5x5.conv'
    elif SEx_FILTER_NAME == '4':
        c.SEx_FILTER_NAME = 'gauss_2.5_5x5.conv'
    elif SEx_FILTER_NAME == '5':
        c.SEx_FILTER_NAME = 'gauss_3.0_5x5.conv'
    elif SEx_FILTER_NAME == '6':
        c.SEx_FILTER_NAME = 'gauss_3.0_7x7.conv'
    elif SEx_FILTER_NAME == '7':
        c.SEx_FILTER_NAME = 'gauss_4.0_7x7.conv'
    elif SEx_FILTER_NAME == '8':
        c.SEx_FILTER_NAME = 'gauss_5.0_9x9.conv'
    SEx_DEBLEND_NTHRESH = raw_input('DEBLEND_NTHRESH (32) >>> ')
    try:
        c.SEx_DEBLEND_NTHRESH = float(SEx_DEBLEND_NTHRESH)
        c.SEx_DEBLEND_NTHRESH = int(c.SEx_DEBLEND_NTHRESH)
    except:
        c.SEx_DEBLEND_NTHRESH = 32
    SEx_DEBLEND_MINCONT = raw_input('DEBLEND_MINCONT (0.005) >>> ')
    try:
        c.SEx_DEBLEND_MINCONT = float(SEx_DEBLEND_MINCONT)
    except:
        c.SEx_DEBLEND_MINCONT = 0.005
    SEx_PHOT_FLUXFRAC = raw_input('PHOT_FLUXFRAC (0.5) >>> ')
    try:
        c.SEx_PHOT_FLUXFRAC = float(SEx_PHOT_FLUXFRAC)
    except:
        c.SEx_PHOT_FLUXFRAC = 0.5
    SEx_pix_scale_disp = 'PIXEL_SCALE (' + str(c.SEx_PIXEL_SCALE) + ') >>> '
    SEx_PIXEL_SCALE = raw_input(SEx_pix_scale_disp)
    try:
        c.SEx_PIXEL_SCALE = float(SEx_PIXEL_SCALE)
    except:
        c.SEx_PIXEL_SCALE = c.pixelscale
    SEx_SEEING_FWHM = raw_input('SEEING_FWHM (0.11) >>> ')
    try:
        c.SEx_SEEING_FWHM = float(SEx_SEEING_FWHM )
    except:
        c.SEx_SEEING_FWHM = c.pixelscale * 3.37
    SEx_BACK_SIZE = raw_input('BACK_SIZE (64) >>> ')
    try:
        c.SEx_BACK_SIZE = float(SEx_BACK_SIZE)
        c.SEx_BACK_SIZE = int(c.SEx_BACK_SIZE)
    except:
        c.SEx_BACK_SIZE = 64
    SEx_BACK_FILTERSIZE = raw_input('BACK_FILTERSIZE (3) >>> ')
    try:
        c.SEx_BACK_FILTERSIZE = float(SEx_BACK_FILTERSIZE)
        c.SEx_BACK_FILTERSIZE = int(c.SEx_BACK_FILTERSIZE)
    except:
        c.SEx_BACK_FILTERSIZE = 3
    SEx_BACKPHOTO_TYPE = raw_input('BACKPHOTO_TYPE (G)LOBAL/(L)OCAL) >>> ')
    while SEx_BACKPHOTO_TYPE != 'G' and SEx_BACKPHOTO_TYPE != 'L' \
          and SEx_BACKPHOTO_TYPE != '':
        SEx_BACKPHOTO_TYPE = raw_input('BACKPHOTO_TYPE (G)LOBAL/(L)OCAL) >>> ')
    if len(SEx_BACKPHOTO_TYPE) == 0:
        c.SEx_BACKPHOTO_TYPE = 'GLOBAL'
    elif SEx_BACKPHOTO_TYPE == 'G':
        c.SEx_BACKPHOTO_TYPE = 'GLOBAL'
    elif SEx_BACKPHOTO_TYPE == 'L':
        c.SEx_BACKPHOTO_TYPE = 'LOCAL'
    if c.SEx_BACKPHOTO_TYPE == 'LOCAL':
        SEx_BACKPHOTO_THICK = raw_input('BACKPHOTO_THICK (24) >>> ')
    try:
        c.SEx_BACKPHOTO_THICK = float(SEx_BACKPHOTO_THICK)
        c.SEx_BACKPHOTO_THICK = int(c.SEx_BACKPHOTO_THICK)
    except:
        c.SEx_BACKPHOTO_THICK = 24
    SEx_WEIGHT_TYPE = raw_input('WEIGHT_TYPE (MAP_RMS) >>> ')
    c.SEx_WEIGHT_TYPE = SEx_WEIGHT_TYPE

