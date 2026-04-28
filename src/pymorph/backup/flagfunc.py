# These functions standardize the flags across all files

class badflag(Exception):
	def __init__(self, value):
		self.value = value
	def __str__(self):
		return repr('Bad flag operation: ' + self.value)
    
def GetFlag(flagname):
    FlagDict = dict([('DETAIL', 3),
                    ('FIT_BULGE', 3),
                    ('FIT_DISK', 4),                         
                    ('FIT_BAR', 5),
                    ('FIT_SKY', 6),
                    ('RA_DEC', 7),
                    ('NEIGHBOUR_FIT', 8),
                    ('CASGM_FAIL', 9),
                    #('EXCEED_SIZE', 10),
                    #('NO_TARGET', 11),
                    ('LARGE_CHISQ', 10),
                    ('FAKE_CNTR', 11),
                    ('IB_AT_LIMIT', 12),
                    ('ID_AT_LIMIT', 13),
                    ('IBAR_AT_LIMIT', 14),
                    ('N_AT_LIMIT', 15),
                    ('NBAR_AT_LIMIT', 16),
                    ('RE_AT_LIMIT', 17),
                    ('RD_AT_LIMIT', 18),
                    ('RBAR_AT_LIMIT', 19),
                    ('EB_AT_LIMIT', 20),
                    ('ED_AT_LIMIT', 21)
                    ])
    return FlagDict[flagname]

def Get_FitFlag(flagname):
    FlagDict = dict([('LARGE_CHISQ', 0),
                    ('FAKE_CNTR', 3),
                    ('IE_AT_LIMIT', 4),
                    ('ID_AT_LIMIT', 5),
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
