#! python3
# This test file tests the functionality of the eprload module
# using data in the eprfiles subdirectory.
from eprload import eprload
from pathlib import Path
import numpy as np
import pytest as pt


baseDir = "./eprfiles"

### load all filenames

#allExt = ['*.dsc','*.dta','*.DSC','*.DTA', '*.spc','*.par','*.PAR','*.SPC','*.eco', '*.d00','*.d01']
implExt = ['*.dsc','*.dta','*.DSC','*.DTA']#, '*.d01','.exp'] # implemented file extensions
srchPath = Path(baseDir)
impFiles = [] 
for ext in implExt:
	impFiles = impFiles + list(srchPath.rglob(ext))

#### Testing Bruker file loading

def test_bruker_object_gen():
	return

### Testing general file load properties

#def test_filestrings():
#	FilePaths = ['./eprfiles/bes3t/00012107','./eprfiles/bes3t/00011201',
#		"./eprfiles/bes3t/00012107","./eprfiles/bes3t/00012107.dta",
#  		"./eprfiles/bes3t/00011201","./eprfiles/bes3t/00011201.spc"]
#	for i in FilePaths:
#		with pt.raises(IOError):
#			eprload(i)

def test_esfmt():
	x = eprload(impFiles[0]).esfmt()
	assert(type(x) is tuple)
	assert(len(x) == 3)
	assert(type(x[0]) is np.ndarray)
	assert(type(x[1]) is np.ndarray)
	assert(type(x[2]) is dict)

def test_parse():
	simpleFiles = [baseDir + "/bes3t/strong1.dta", baseDir + "/esp/strong1esp.spc"]
	sf1 = eprload(simpleFiles[0])
	assert np.isreal(sf1.Param["XPTS"])
	assert np.isreal(sf1.Param["XMIN"])
	assert np.isreal(sf1.Param["XWID"])
	with pt.raises(NotImplementedError):
		sf2 = eprload(simpleFiles[1]) # have not yet implemented spc files
		assert np.isreal(sf2.Param["MF"])
		assert np.isreal(sf2.Param["GST"])
		assert np.isreal(sf2.Param["MP"])
#def test_fileload():
#	ext = ['*.dsc','*.dta','*.DSC','*.DTA', '*.spc','*.par',
#		'*.PAR','*.SPC','*.eco', '*.d00','*.d01']
#	# search for files
#	# load all files. Foreach:
#	# 	check that data is correct type
#	# 	plot data in default mode??



