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
			 ('FIT_BULGE', 3),
                         ('FIT_DISK', 4),                         
                         ('FIT_SKY', 5),
			 ('FIT_POINT', 6),
			 ('FIT_BAR', 7),
			 ('NEIGHBOUR_FIT', 8),
                         ('EXCEED_SIZE', 9),
                         ('NO_TARGET', 10),
			 ('ASYM_NOT_CONV', 11),
                         ('ASYM_OUT_FRAME', 12),
                         ('ELLIPSE_FAIL', 13),
                         ('CASGM_FAIL', 14),
                         ('GALFIT_FAIL', 15),
                         ('PLOT_FAIL', 16),
			 ('ERRORS_FAILED', 17),
			 ('AVGIE_FAILED', 18),
                         ('BACK_FAILED', 19),
                         ('DETAIL_FAILED', 20)
                         ])
        return FlagDict[flagname]

def Get_FitFlag(flagname):
        FlagDict = dict([('LARGE_CHISQ', 0),
                         ('SMALL_GOODNESS', 1),
                         ('FAKE_CNTR', 2),
                         ('IE_AT_LIMIT', 3),
                         ('ID_AT_LIMIT', 4),
			 ('RERD_AT_LIMIT', 5),
                         ('BT_AT_LIMIT', 6),
			 ('N_AT_LIMIT', 7),
			 ('RE_AT_LIMIT', 8),
			 ('RD_AT_LIMIT', 9),
			 ('EB_AT_LIMIT', 10),
			 ('ED_AT_LIMIT', 11)
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
