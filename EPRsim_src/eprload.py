# Author: Liam Twomey

from pathlib import Path,PurePath
from re import search
from warnings import warn
from pprint import pprint

# all logic adapted from Stefan Stoll's EasySpin.
class eprload:
	def __init__(self, fileName='.', scaling:str=None):
		'''
		arguments: filename: the name of a Bruker BES3T file (.DSC or .DTA), scale: a scaling factor for the data
		outputs: returns tuple of three arguments, the abscissa (typically labelled B) the spectrum (typically S) and the experimental parameters (typically P)
		'''
		self.filePath = Path(fileName) #if already a Path this won't effect anything
		self.fileExt = self.filePath.suffix
		self.stdFileExt = self.fileExt.lower()
		self.scaling = scaling
		self.checkFileType()

	def checkFileType(self):
		if self.scaling == None:
			scale=False
		else:
			scale = True
		match self.stdFileExt:
			case '.dsc'|'.dta':	
				FileFormat = 'BrukerBES3T'
				if scale:
					if scaling.lower() in "npgtc":
						self.loadBrukerBES3T(self.filePath,self.fileExt,self.scaling)
					else:
						warn("Invalid scaling option supplied. Valid options are: N,P,G,T,C.")
				else:
					self.loadBrukerBES3T(self.filePath,self.fileExt)
			
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
	
	def loadBrukerBES3T(self,fileName, fileExt, scaling=None):
		# BES3T file processing
		# (Bruker EPR Standard for Spectrum Storage and Transfer)
		#    .DSC: descriptor file
		#    .DTA: data file
		# used on Bruker ELEXSYS and EMX machines
		# Code based on BES3T version 1.2 (Xepr 2.1)
		self.P = self.readDSCfile(fileName, fileExt)
		# IKKF: Complex-data Flag
		# CPLX indicates complex data, REAL indicates real data.
		if 'IKKF' in self.P.keys():
			isComplex = []
			nValsPerPoint = len(self.P['IKKF'])
			for dim in self.P['IKKF']:
				if dim =='CPLX':
					isComplex.append(1)
				elif dim=='REAL':
					isComplex.append(0)
				else:
					warn("Unknown dimension complexity {dim} in IKKF field of .DSC file!")
		else:
			warn('Keyword IKKF not found in .DSC file! Assuming IKKF=REAL.');
			isComplex = 0;
			nValsPerPoint = 1;
		# XPTS, YPTS, ZPTS specify the number of data points in
		#  x, y and z dimension.
		
		if 'XPTS' in self.P:
			nx = int(self.P['XPTS'][0]) 
		else:
			raise ValueError("No XPTS parameter found in {self.fileName}")
		ny = int(self.P['YPTS'][0]) if 'YPTS' in self.P else 0
		nz = int(self.P['ZPTS'][0]) if "ZPTS" in self.P else 0
		self.dimensions = [nx,ny,nz]

		# BSEQ: Byte Sequence of encoding machine
		# BSEQ describes the byte order of the data. BIG means big-endian,
		# LIT means little-endian. Sun and Motorola-based systems are
		# big-endian (MSB first), Intel-based system little-endian (LSB first).
		if "BSEQ" in self.P:
			if self.P["BSEQ"]=='BIG':
				ByteOrder = 'ieee-be'
			elif self.P["BSEQ"]=="LIT":
				ByteOrder = 'ieee-le'
			else:
				"Unknown byte order specified by BSEQ keyword of {filename}"
		else:
			warn('Keyword BSEQ not found in .DSC file! Assuming BSEQ=LIT.')
			ByteOrder = 'ieee-le'
		# IRFMT: Item Real Format
		# IIFMT: Item Imaginary Format
		# Data format tag of BES3T is IRFMT for the real part and IIFMT
		# for the imaginary part.
		if "IRFMT" in self.P:
			if len(self.P["IRFMT"]) != nValsPerPoint:
				raise ValueError('Problem in BES3T DSC file: inconsistent IKKF and IRFMT fields.')
			for char in self.P["IRFMT"][0].split(','):
				match char.upper():
					case "C": 
						numberFormat = 'int8'
					case "S":
						numberFormat = 'int16'
					case "I":
						numberFormat = 'int32'
					case "F":
						numberFormat = 'float32'
					case "D":
						numberFormat = 'float64'
					case "A":
						raise NotImplementedError("Cannot read BES3T data in ASCII format.")
					case "0"|"N": raise ValueError("No BES3T data found.")
					case _: raise ValueError("Unknown value for IRFMT in .DSC file.")
		else:
			raise ValueError("Keyword IRFMT not found in .DSC file.")
		# We enforce IRFMT and IIFMT to be identical.
		if "IIFMT" in self.P:
			if self.P["IIFMT"] != self.P["IRFMT"]:
				raise ValueError("IRFMT and IIFMT in .DSC file are not equal.")

		# constructing abscissa vectors
		axisNames = ['x','y','z']
		for a in range(len(axisNames)):
			if self.dimensions(a) >=1:
				self.axisType = self.P["{axisNames[a]}TYP"]
				if axisType=="IGD":
					#nonlinear axis, try to load companion file (.xgf, .ygf, .zgf)
					self.abscissa = self.readNonlinearAbscissa(axisType)
				elif axisType=="IDX":
					absc_min[a] = self.P["{axisNames[a]}MIN"]
					absc_width[a] = self.P["{axisNames[a]}WID"]
					if absc_width[a] == 0:
						warn("{axisNames[a]} has a width of 0.")
						absc_min[a] == 1
						absc_width[a] = self.dimensions[a]-1
					self.abscissa[a] = absc_min[a] + np.linspace(0,absc_width[a],num=self.dimensions[a])
				elif axisType=="NTUP":
					raise NotImplementedError('Cannot read data with NTUP axes.')
		if len(self.abscissa)==1:
			self.abscissa = self.abscissa[0][:]
		
		# Here's the hard part. Stefan uses the matlab getmatrix() function, but what I need to do is
		# - get binary data from spec file, specifying byte order
		# - parse the data for use as a complex number
		# - shape data into ndarray with dimensions
		self.data = None

	def readDSCfile(self,fileName, extension, scaling=None):
		'''
		@spec
		Read all lines from file into array (\n-delimited)
		if line is terminated by \\ , append next line - not implemented
		split each line into key-value pairs
		Return datatype: dict with str keys and list values.
		'''
		cleanKeys = []
		cleanVals = []
		with open(fileName,'r') as dictFile:
			for line in dictFile:
				lstrip = line.strip()
				if (len(lstrip) != 0) and (lstrip[0] not in '#*.'):
					lineList = line.split()
					cleanKeys.append(lineList[0])
					# if len(lineList) == 2:
					# 	cleanVals.append(lineList[1])
					# else:
					cleanVals.append(lineList[1:])
		lineDict = dict(zip(cleanKeys,cleanVals))
		pprint(lineDict)
		return lineDict

def readNonlinearAbscissa(self, axisNames):
	### TODO
# 	    companionFileName = [FullBaseName '.' AxisNames{a} 'GF'];
#     % Determine data format form XFMT/YMFT/ZFMT
#     DataFormat = Parameters.([AxisNames{a} 'FMT']);
#     switch DataFormat
#       case 'D', sourceFormat = 'float64';
#       case 'F', sourceFormat = 'float32';
#       case 'I', sourceFormat = 'int32';
#       case 'S', sourceFormat = 'int16';
#       otherwise
#         error('Cannot read data format %s for companion file %s',DataFormat,companionFileName);
#     end
#     % Open and read companion file
#     fg = fopen(companionFileName,'r',ByteOrder);
#     if fg>0
#       Abscissa{a} = fread(fg,Dimensions(a),sourceFormat,ByteOrder);
#       fclose(fg);
#     else
#       warning('Could not read companion file %s for nonlinear axis. Assuming linear axis.',companionFileName);
#       AxisType = 'IDX';
#     end
#   end