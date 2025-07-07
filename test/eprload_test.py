#! python3
# This test file tests the functionality of the eprload module
# using data in the eprfiles subdirectory.
from eprsim.EPRload import * 
from pathlib import Path
import numpy as np
import pytest as pt

baseDir = "./eprfiles"

def getFiles(patterns:list,path:Path):
	files = []
	for e in patterns:
		files += list(path.rglob(e))
	return files

### find all filenames
#allExt = ['*.dsc','*.dta','*.DSC','*.DTA', '*.spc','*.par','*.PAR','*.SPC','*.eco', '*.d00','*.d01']
bExt = ['*.dsc','*.dta','*.DSC','*.DTA'] # bruker file extensions
smExt = ['*.d01','.exp'] # specman file extensions
srchPath = Path(baseDir)
bFiles = getFiles(bExt, srchPath)
smFiles = getFiles(smExt,srchPath)

#### Testing Bruker file loading
def test_loadall_bruker():
	objs = []
	for i in bFiles:
		objs.append(eprload(i))

### Testing general file load properties

#def test_filestrings():
#	FilePaths = ['./eprfiles/bes3t/00012107','./eprfiles/bes3t/00011201',
#		"./eprfiles/bes3t/00012107","./eprfiles/bes3t/00012107.dta",
#  		"./eprfiles/bes3t/00011201","./eprfiles/bes3t/00011201.spc"]
#	for i in FilePaths:
#		with pt.raises(IOError):
#			eprload(i)

def test_alias():
	x = eprload(bFiles[0])
	assert(type(x.B) is np.ndarray)
	assert(type(x.S) is np.ndarray)
	assert(type(x.P) is dict)
	with pt.raises(AttributeError):
		assert(type(x.W) is np.ndarray)

def test_parse_dsc():
	simpleFiles = [baseDir + "/bes3t/strong1.dta", baseDir + "/esp/strong1esp.spc"]
	sf1 = eprload(simpleFiles[0])
	assert np.isreal(sf1.Param["XPTS"])
	assert np.isreal(sf1.Param["XMIN"])
	assert np.isreal(sf1.Param["XWID"])
	with pt.raises(NotImplementedError):
		sf2 = eprload(simpleFiles[1]) # have not yet implemented spc files
		assert(np.isreal(sf2.Param["MF"]))
		assert(np.isreal(sf2.Param["GST"]))
		assert(np.isreal(sf2.Param["MP"]))

#def test_fileload():
#	ext = ['*.dsc','*.dta','*.DSC','*.DTA', '*.spc','*.par',
#		'*.PAR','*.SPC','*.eco', '*.d00','*.d01']
#	# search for files
#	# load all files. Foreach:
#	# 	check that data is correct type
#	# 	plot data in default mode??



