# These functions standardize the flags across all files

class badflag(Exception):
	def __init__(self, value):
		self.value = value
	def __str__(self):
		return repr('Bad flag operation: '+self.value)
    
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
                         ('ASYM_NOT_CONV', 13),
                         ('ASYM_OUT_FRAME', 14),
			 ('FIT_BAR', 15),
			 ('ERRORS_FAILED', 16),
			 ('AVGIE_FAILED', 17),
                         ('NO_TARGET', 18),
                         ('BACK_FAILED', 20),
                         ('DETAIL_FAILED', 21),
                         ])
        return FlagDict[flagname]

def Get_FitFlag(flagname):
        FlagDict = dict([('LARGE_CHISQ', 0),
                         ('SMALL_GOODNESS', 1),
                         ('FAKE_CNTR', 2),
                         ('IE_AT_LIMIT', 3),
                         ('ID_AT_LIMIT', 4),
			 ('BACK_FAILED', 5),
                         ('RERD_AT_LIMIT', 6),
                         ('BT_AT_LIMIT', 7),
			 ('N_AT_LIMIT', 8),
			 ('RE_AT_LIMIT', 9),
			 ('RD_AT_LIMIT', 10),
			 ('EB_AT_LIMIT', 11),
			 ('ED_AT_LIMIT', 12)
			 ])	                 
        return FlagDict[flagname]



def isset(flag, bit):
        """Return True if the specified bit is set in the given bit mask"""
        return (flag & (1 << bit)) != 0

def SetFlag(Flag, bit):
	if isset(Flag, bit):
		raise badflag("Tried to set a flag that was already set")
	Flag += 2**bit
	return Flag

def ClrFlag(Flag, bit):
	if not isset(Flag, bit):
		raise badflag("Tried to clear a flag that was not set")
	Flag -= 2**bit
	return Flag
