# These functions standardize the flags across all files

class badflag(Exception):
	def __init__(self, value):
		self.value = value
	def __str__(self):
		return repr('Bad flag operation: ' + self.value)
    
def GetFlag(flagname):

     FlagDict = dict([
                     ('FIT_BULGE', 3),
                     ('FIT_DISK', 4),                         
                     ('FIT_BAR', 5),
                     ('BULGE_BOXINESS', 6),
                     ('BAR_BOXINESS', 7),
                     ('SKY_FIT', 8),
                     ('DETAILED_FIT', 9),
                     ("CUT_OUT_IMAGE", 10),
                     ("FIND_AND_FIT", 11),
                     ('RA_DEC', 12),
                     ('NEIGHBOUR_FIT', 13),
                     ('CASGM_FAIL', 14),
                     ('PARAM_ERROR', 15),
                     ("LARGE_CHISQ", 16),
                     ("FAKE_CNTR", 17),
                     ("IB_AT_LIMIT", 18),
                     ("N_AT_LIMIT", 19),
                     ("RE_AT_LIMIT", 20),
                     ("ID_AT_LIMIT", 21),
                     ("RD_AT_LIMIT", 22),
                     ("IBAR_AT_LIMIT", 23),
                     ("NBAR_AT_LIMIT", 24),
                     ("RBAR_AT_LIMIT", 25)
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
