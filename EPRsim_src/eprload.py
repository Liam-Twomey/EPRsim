# Author: Liam Twomey

from pathlib import Path
from re import search
from warnings import warn
from pprint import pprint

# all logic adapted from Stefan Stoll's EasySpin.

def eprload(filename:str='.', scaling:str=None):
	'''
	arguments: filename: the name of a Bruker BES3T file (.DSC or .DTA), scale: a scaling factor for the data
	outputs: returns tuple of three arguments, the abscissa (typically labelled B) the spectrum (typically S) and the experimental parameters (typically P)
	'''
	filePath = Path(filename)
	stdFileExt = filePath.suffix.lower()
	if scaling == None:
		scale=False
	else:
		scale = True
	match stdFileExt:
		case '.dsc'|'.dta':	
			FileFormat = 'BrukerBES3T'
			if scale:
				if scaling.lower() in "npgtc":
					Data,Abscissa,Parameters = eprload_BrukerBES3T(filePath,filePath.suffix,scaling)
				else:
					warn("Invalid scaling option supplied. Valid options are: N,P,G,T,C.")
			else:
				Data,Abscissa,Parameters = eprload_BrukerBES3T(filePath,filePath.suffix)
		
		case '.par'|'.spc':
			FileFormat = 'BrukerESP'
			warn("File type {FileFormat} not yet implemented.")
			if scale:
				warn("Scaling not supported for this filetype.")
			# [Data,Abscissa,Parameters] = eprload_BrukerESP(FullBaseName,FileExtension,Scaling)
		
		case '.d01':
			FileFormat = 'SpecMan'
			warn("File type {FileFormat} not yet implemented.")
			if scale:
				warn("Scaling not supported for this filetype.")
			# [Data,Abscissa,Parameters] = eprload_specman(FileName)
		
		case '.spe':
			FileFormat = 'MagnettechBinary'
			warn("File type {FileFormat} not yet implemented.")
			if scale:
				warn("Scaling not supported for this filetype.")
			# [Data,Abscissa,Parameters] = eprload_MagnettechBinary(FileName)
		
		case '.xml':
			FileFormat = 'MagnettechXML'
			warn("File type {FileFormat} not yet implemented.")
			if scale:
				warn("Scaling not supported for this filetype.")
			# [Data,Abscissa,Parameters] = eprload_MagnettechXML(FileName)
		
		case '.esr':
			FileFormat = 'ActiveSpectrum'
			warn("File type {FileFormat} not yet implemented.")
			if scale:
				warn("Scaling not supported for this filetype.")
			# [Data,Abscissa,Parameters] = eprload_ActiveSpectrum(FileName)
		
		case '.dat':
			FileFormat = 'AdaniDAT'
			warn("File type {FileFormat} not yet implemented.")
			if scale:
				warn("Scaling not supported for this filetype.")
			# [Data,Abscissa,Parameters] = eprload_AdaniDAT(FileName)

		case '.json':
			FileFormat = 'AdaniJSON'
			warn("File type {FileFormat} not yet implemented.")
			if scale:
				warn("Scaling not supported for this filetype.")
			# [Data,Abscissa,Parameters] = eprload_AdaniJSON(FileName)
		
		case '.eco':
			FileFormat = 'qese/tryscore'
			warn("File type {FileFormat} not yet implemented.")
			if scale:
				warn("Scaling not supported for this filetype.")
			[Data,Abscissa,Parameters] = eprload_qeseETH(FileName)
		
		case '.plt':
			FileFormat = 'MAGRES'
			warn("File type {FileFormat} not yet implemented.")
			if scale:
				warn("Scaling not supported for this filetype.")
			# [Data,Abscissa,Parameters] = eprload_MAGRES(FileName)
		
		case '.spk'|'.ref':
			FileFormat = 'VarianETH'
			warn("File type {FileFormat} not yet implemented.")
			if scale:
				warn("Scaling not supported for this filetype.")
			[Data,Abscissa,Parameters] = eprload_VarianE9ETH(FileName)
		
		case '.d00':
			FileFormat = 'WeizmannETH'
			warn("File type {FileFormat} not yet implemented.")
			if scale:
				warn("Scaling not supported for this filetype.")
			# [Data,Abscissa,Parameters] = eprload_d00WISETH(FileName)
		
		case _:
			# Test for JEOL file
			with open(filename,"rb") as fptr:
				identity = fptr.readline(16) # first 16 bytes = 16 ascii/utf-8 characters.
			iddc = identity.decode('latin1')
			isJeol = search('^spin|^cAcqu|^endor|^pAcqu|^cidep|^sod|^iso|^ani',iddc).group()
			if len(isJeol) >= 0:
				FileFormat = 'JEOL'
				warn("File type {FileFormat} not yet implemented.")
				if scale:
					warn("Scaling not supported for this filetype.")
				# [Data,Abscissa,Parameters] = eprload_jeol(FileName);

			else:
				raise NotImplementedError("Unsupported file extension {stdFileExt}")
	return (Data, Abscissa, Parameters)

def eprload_BrukerBES3T(fileName, fileExt, scaling=None):
	# BES3T file processing
	# (Bruker EPR Standard for Spectrum Storage and Transfer)
	#    .DSC: descriptor file
	#    .DTA: data file
	# used on Bruker ELEXSYS and EMX machines
	# Code based on BES3T version 1.2 (Xepr 2.1)
	Parameters,err = readDSCfile(fileName, fileExt)
	print(Parameters,err)
	return (-1,-1,-1)

def readDSCfile(fileName, extension, scaling=None):
	'''
	@spec
	Read all lines from file into array (\n-delimited)
	if line is terminated by \, append next line
	split each line into key-value pairs
	'''
	Parameters = []
	err = []
	# with open(fileName,'r') as dictFile:
	# 	reader = DictReader(dictFile,fieldnames=fnames)
	# 	for row in reader:
	# 		print('|',row['key'],'|',row['value'],'|')
	# 	print(reader)
	cleanKeys = []
	cleanVals = []
	with open(fileName,'r') as dictFile:
		for line in dictFile:
			lstrip = line.strip()
			if (len(lstrip) != 0) and (lstrip[0] not in '#*.'):
				lineList = line.split()
				# cleanKeys.append(lineList)
				cleanKeys.append(lineList[0])
				if len(lineList) == 2:
					cleanVals.append(lineList[1])
				else:
					cleanVals.append(lineList[1:])
	lineDict = dict(zip(cleanKeys,cleanVals))
	pprint(lineDict)
	# print(cleanKeys)
	return (None, None)

def readDTAfile(DSCfilename):
	None

warn("running class source is not intended use")