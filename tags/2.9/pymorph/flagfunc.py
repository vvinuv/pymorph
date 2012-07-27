# These functions standardize the flags across all files


def GetFlag(flagname):
        FlagDict = dict([('REPEAT', 0),
                         ('FIT_BULGE_CNTR', 1),
                         ('FIT_DISK_CNTR', 2),
                         ('FIT_SKY', 3),
                         ('EXCEED_SIZE', 4),
                         ('ELLIPSE_FAIL', 5),
                         ('CASGM_FAIL', 6),
                         ('GALFIT_FAIL', 7),
                         ('PLOT_FAIL', 8),
                         ('FIT_BULGE', 9),
                         ('FIT_DISK', 10),
                         ('FIT_POINT', 11),
                         ('NEIGHBOUR_FIT', 12),
                         ('LARGE_CHISQ', 13),
                         ('SMALL_GOODNESS', 14),
                         ('FAKE_CNTR', 15),
                         ('BULGE_AT_LIMIT', 16),
                         ('DISK_AT_LIMIT', 17),
                         ('ASYM_NOT_CONV', 18),
                         ('ASYM_OUT_FRAME', 19),
                         ('BACK_FAILED', 20),
                         ('DETAIL_FAILED', 21),
                         ('RERD_AT_LIMIT', 22),
                         ('BT_AT_LIMIT', 23),
			 ('N_AT_LIMIT', 24),
			 ('RE_AT_LIMIT', 25)])
        return FlagDict[flagname]
def isset(flag, bit):
        """Return True if the specified bit is set in the given bit mask"""
        return (flag & (1 << bit)) != 0
    
