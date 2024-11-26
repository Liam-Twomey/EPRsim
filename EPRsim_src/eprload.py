from pathlib import Path
from re import search
from warnings import warn

# all logic adapted from Stefan Stoll's EasySpin.

def eprload(filename:str='.', scaling:str=''):
	'''
	arguments: filename: the name of a Bruker BES3T file (.DSC or .DTA), scale: a scaling factor for the data
	outputs: returns tuple of three arguments, the abscissa (typically labelled B) the spectrum (typically S) and the experimental parameters (typically P)
	'''
	filePath = Path(filename)
	stdFileExt = filePath.suffix.lower()
	if len(scaling) ==0:
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

def eprload_BrukerBES3T():
	# BES3T file processing
	# (Bruker EPR Standard for Spectrum Storage and Transfer)
	#    .DSC: descriptor file
	#    .DTA: data file
	# used on Bruker ELEXSYS and EMX machines
	# Code based on BES3T version 1.2 (Xepr 2.1)
	print ('ok')
	return

def readDSCfile():

def readDTAfile():