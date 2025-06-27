#! python3
# Author: Liam Twomey

from pathlib import Path,PurePath
from re import sub,match,search
import numpy as np
#from sre_constants import CHARSET
from warnings import warn
from pprint import pprint

# all logic adapted from Stefan Stoll's EasySpin.
class eprload:
	def __init__(self, fileName:(str | Path), scaling:str='1',verbose:bool=False): 
		'''
        The parent eprload class acts as a wrapper for instrument-specific `eprload` functions.
        It determines the filetype of the entered file, then calls the appropriate method.

		Args:
            * filename: the name of a EPR data file
            * scaling: a scaling factor for the data
			* verbose: a flag passed for printing of debug statements (specified by self.vprint())
			* esfmt: (WIP) a flag passed to allow return of just the self.Absc, self.Spec, and self 
        returns:
            (absc,spec, par) where...
                * absc: abscissa (typically labelled B)
                * spec: spectrum (typically labelled S) 
                * par:  experimental parameters (typically labelled P)
		'''
		self.filePath = Path(fileName) #if already a Path this won't effect anything
		self.fileExt = self.filePath.suffix
		self.scaling = scaling.upper()
		# self.Param = {}
		self.verbose = verbose
		self.vprint(f'`eprload` initalized on \"{self.filePath}\"')
		if self.fileExt.isupper():
			self.extCase = 1
		elif self.fileExt.islower():
			self.extCase = 0
		else:
			raise IOError("Please don't use mixed-case extensions!")
		self.checkFileType()
	
	def esfmt(self):
		return (self.Absc, self.Spec, self.Param)
		
	def show_params(self):
		'''
		this is my current workaround to not being able to print properly.
		'''
		print('Parameters:')
		pprint(self.Param)

	def vprint(self,output,indentLevel:int=0):
		'''
		A debug printer for easy output of debug/verbose info only when self.verbose=True
		Prints `output` tabbed in by `indentlevel` tabs.
		Does not tab output when type(output) is not a builtin.
		'''
		if indentLevel==0:
			initStr = '>'
		else:
			initStr = '\t'*indentLevel
		if not self.verbose: # get out of here ASAP if self.verbose is false
			return
		elif (output.__class__.__module__=='builtins'): # if obj is builtin class, doesn't need pprinting
			print(initStr,output) # print indentLevel tabs, then the output
			return
		
		else: # We fall here if the object needs pretty printing (i.e. is not builtin)
			pprint(output)
			return

	def checkFileType(self):
		match self.fileExt.upper():
			case '.DSC'|'.DTA':	
				self.vprint('File type is Bruker BES3T')
				# 1 is the default internal value
				if self.scaling.upper() in "NGPTC1":
					self.loadBES3T()
				else:
					raise RuntimeError("Invalid scaling option supplied. Valid options are: N,P,G,T,C.")
			case '.D01':
				self.vprint('File type is SpecMan')
				raise NotImplementedError("File type {FileFormat} not yet implemented.")
				#if self.scaling[0]!='1':
				#	warn("Scaling not supported for this filetype.")
			case ('.PAR'|'.SPC'|'.SPE'|'.XML'|'.ESR'|'.DAT'|'.JSON'|'.ECO'|'.PLT'|'.SPK'|'.REF'|'.D00'):
				#Cover all other file formats 
				raise NotImplementedError("File type {FileFormat} not yet implemented and does not have a roadmap. Contact maintainers if you need it implemented.")
			case _:
				# Test for JEOL file
				with open(self.filePath,"rb") as fptr:
					identity = fptr.readline(16) # first 16 bytes = 16 ascii/utf-8 characters.
				iddc = identity.decode('latin1')
				try:
					isJeol = search('^spin|^cAcqu|^endor|^pAcqu|^cidep|^sod|^iso|^ani',iddc).group()
					if len(isJeol) >= 0:
						self.vprint('File type is JEOL')
						if self.scaling[0]!='1':
							warn("Scaling not supported for this filetype.")
						# [Data,Abscissa,Parameters] = eprload_jeol(FileName);
				except AttributeError:
					raise NotImplementedError(f"File format {self.fileExt} not yet implemented.")

				else:
					raise NotImplementedError(f"Unsupported file extension {self.fileExt}")

# Parse arguments
	def loadBES3T(self) -> None:
		'''
		BES3T file processing
		(Bruker EPR Standard for Spectrum Storage and Transfer)
		* .DSC: descriptor file
		* .DTA: data file
		used on Bruker ELEXSYS and EMX machines
		Code based on BES3T version 1.2 (Xepr 2.1)
		'''
#		BES3TParamLoad()
#		BES3TParamParse()

#	def BES3TParamLoad(self):
		self.Param = {}
		dscExt = ['.DSC' if self.extCase else '.dsc'][0]
		self.vprint(f"Reading parameters from {dscExt} file to self.Param")
		with open(self.filePath.with_suffix(dscExt)) as tmpf:
			lines = tmpf.readlines()
			for line in lines:
				if (not match(r"[*,#]",line)) and (len(line) > 0):
					tmpln = sub(r"\s+",' ',line).split()
					for attr in range(len(tmpln)):
						#cast to ints and floats where possible
						try:
							tmpln[attr]  = int(tmpln[attr])
						except ValueError:
							try:
								tmpln[attr] = float(tmpln[attr])
							except ValueError:
								None
					# this returns all keys as lists of length 1 or longer.
					if len(tmpln[1:]) > 0:
						self.Param[tmpln[0]] = tmpln[1:]
	
#	def BES3TParamParse(self):
		self.vprint('Determining spectrum size and properties')
		# IKKF: Complex-data Flag
		# CPLX indicates complex data, REAL indicates real data.
		isComplex = []
		if 'IKKF' in self.Param.keys():
			nValsPerPoint = len(self.Param['IKKF'])
			for dim in self.Param['IKKF']:
				if dim =='CPLX':
					isComplex.append(1)
				elif dim=='REAL':
					isComplex.append(0)
				else:
					warn(f"Unknown dimension complexity {dim} in IKKF field of .DSC file!")
		else:
			warn('Keyword IKKF not found in .DSC file! Assuming IKKF=REAL.')
			isComplex[0] = 0
			nValsPerPoint = 1 # this is never passed to getBinaryMatrix
		isComplex = np.array(isComplex)
		# XPTS, YPTS, ZPTS specify the number of data points in
		#  x, y and z dimension.
		if 'XPTS' in self.Param:
			nx = int(self.Param['XPTS'][0])
		else:
			raise ValueError("No XPTS parameter found in {self.fileName}")
		ny = int(self.Param['YPTS'][0]) if 'YPTS' in self.Param else 0
		nz = int(self.Param['ZPTS'][0]) if "ZPTS" in self.Param else 0
		self.dimensions = [nx,ny,nz]

		# BSEQ: Byte Sequence of encoding machine
		# BSEQ describes the byte order of the data. BIG means big-endian,
		# LIT means little-endian. Sun and Motorola-based systems are
		# big-endian (MSB first), Intel-based system little-endian (LSB first).
		if "BSEQ" in self.Param:
			if self.Param["BSEQ"][0]=='BIG':
				self.byteOrder = '>'
			elif self.Param["BSEQ"][0]=="LIT":
				self.byteOrder = '<'
			else:
				"Unknown byte order specified by BSEQ keyword of {filename}"
		else:
			warn('Keyword BSEQ not found in .DSC file! Assuming BSEQ=LIT.')
			self.byteOrder = '<'
		# IRFMT: Item Real Format
		# IIFMT: Item Imaginary Format
		# Data format tag of BES3T is IRFMT for the real part and IIFMT
		# for the imaginary part.
		if "IRFMT" in self.Param:
			if len(self.Param["IRFMT"]) != nValsPerPoint:
				raise ValueError('Problem in BES3T DSC file: inconsistent IKKF and IRFMT fields.')
			for char in self.Param["IRFMT"][0].split(','):
				match char.upper():
					# assign datatype of component data. Numpy naming uses byte len not bit.
					case "C": 
						numberFormat = 'i1' # int8
					case "S":
						numberFormat = 'i2' # 'int16'
					case "I":
						numberFormat = 'i4' # 'int32'
					case "F":
						numberFormat =  'f4' # 'float32'
					case "D":
						numberFormat = 'f8' # 'float64'
					case "A":
						raise NotImplementedError("Cannot read BES3T data in ASCII format.")
					case "0"|"N":
						raise ValueError("No BES3T data found.")
					case _:
						raise ValueError("Unknown value for IRFMT in .DSC file.")
		else:
			raise ValueError("Keyword IRFMT not found in .DSC file.")
		# We enforce IRFMT and IIFMT to be identical.
		if "IIFMT" in self.Param:
			if self.Param["IIFMT"] != self.Param["IRFMT"]:
				raise ValueError("IRFMT and IIFMT in .DSC file are not equal.")

#	def BES3TAbscLoad(self):
		self.vprint('Reading abscissas into self.Absc')
		axisNames = ['X','Y','Z']
		# I think I can rewrite this for loop in a way that makes more sense in python @TODO
		self.Absc = np.full([len(axisNames), max(self.dimensions)], np.nan) 
		for a in range(len(axisNames)):
			axName = axisNames[a]
			if self.dimensions[a] >=1:
				axisType = self.Param[f"{axisNames[a]}TYP"][0]
				self.vprint(f'Reading {axName} asbscissa of type {axisType}',1)
				if axisType=="IGD":
					#nonlinear axis, try to load companion file (.xgf, .ygf, .zgf)
					tmp_nla = self.readNonlinearAbscissa(a, axisType) 
					if type(tmp_nla) is np.ndarray:
						self.Absc[a,:] = tmp_nla
					elif type(tmp_nla) is str: 
						axisType = ''
					else:
						raise RuntimeError("Undefined return from readNonLinearAbscissa")
					del tmp_nla
				if axisType=="IDX": # not elif to allow for nonlinearabscissa error handling 2l above
					absc_min = np.zeros(len(axisNames))
					absc_width = np.zeros(len(axisNames))
					absc_min[a] = self.Param[f"{axisNames[a]}MIN"][0]
					absc_width[a] = self.Param[f"{axisNames[a]}WID"][0]
					if absc_width[a] == 0:
						warn(f"{axisNames[a]} has a width of 0.")
						absc_min[a] = 1
						absc_width[a] = self.dimensions[a]-1
					self.Absc[a,:] = np.linspace(absc_min[a],absc_min[a]+absc_width[a],num=self.dimensions[a])
					self.vprint(f'{axName} abscissa size: {self.Absc[a].shape}',1)
				elif axisType=="NTUP":
					raise NotImplementedError('Cannot read data with NTUP axes.')
				else:
					raise FileNotFoundError('AxisType is not defined for axis {axisNames[a]}')
			else:
				self.vprint(f'No abscissa for axis {axName}',1)
			self.Absc = self.Absc.flatten()
	#def BES3TSpecLoad(self):
		# get data from .dta file
		dtaExt = ['.DTA' if self.extCase else '.dta'][0]
		self.vprint(f'Reading data from {dtaExt} file in {self.byteOrder}{numberFormat} format to self.Spec')
		self.Spec = self.readBinaryDataMatrix(dtaExt,numberFormat,isComplex)
		if self.scaling[0] != '1':
			self.Spec = self.scaleData()

	def readNonlinearAbscissa(self, a, axisNames):
		'''
		This function just parses the format of the accessory .XGF, .YGF, and .ZGF files for axisType IGD.
		It determines the encoding method of these 
		'''
		#Nonlinear axis -> Try to read companion file (.XGF, .YGF, .ZGF)
		companionFile = self.filePath.with_suffix(f'.{axisNames[a]}GF')
		# Determine data format form XFMT/YMFT/ZFMT
		DataFormat = self.Param[f'{axisNames[a]}FMT']
		match DataFormat:
			case 'D':
				sourceFormat = 'f8' #'float64'
			case 'F':
				sourceFormat = 'f4' #'float32'
			case 'I':
				sourceFormat = 'i4' #'int32'
			case 'S':
				sourceFormat = 'i2' #'int16'
			case '_':
				raise IOError(f'Cannot read data format {0} for companion file {1}',DataFormat,companionFile);
		# Open and read companion file
		try:
			tmpAbsc = np.fromfile(companionFile,dtype=self.byteOrder+sourceFormat)
			# with open(companionFile,'rb') as ocf:
			# 	self.Abscissa[a] = ocf.read(Dimensions[a], SourceFormat, byteOrder)
			return tmpAbsc
		except:
			warn(f'Could not read companion file {companionFile} for nonlinear axis. Assuming linear axis.')
			axisType = 'IDX'
			return axisType

	def readBinaryDataMatrix(self,fileExt,numberFormat,isComplex) -> np.ndarray:
		'''
		Description of Matlab function, line 147 onward:
		Data = getmatrix([FullBaseName,SpcExtension],Dimensions,numberFormat,byteOrder,isComplex);
		* opens file fullbasename.spcextension with format byteOrder+numberFormat
		* Real and imaginary data are interspersed, so it just reads everything into a matrix
		sized to the total number of real points in all dimensions (nx*ny*nz).
		* This data is then reshaped into a (nRealsPerPoint x (nx*ny*nz)) array, so each row is a different
		real value index, and each column is a new datapoint. 
		* Now we cope with each datapoint having nDataValuesPerPoint values.
			* Note: both Matlab and numpy index [row, column] so format is the same.
			* We run through the isComplex array and compare it against the rows of the data,
			and if a row isComplex, then it and the row after it are combined into a complex number.
			* If the row is not complex, then it is left alone.
		* The resulting matrix is then reshaped to a [nx*ny*nz] array.
		'''
		isComplex = np.array(isComplex)
		self.vprint('Axis complexity is: '+str(isComplex).replace('1','CPLX').replace('0','REAL'),1)
		data = np.fromfile(self.filePath.with_suffix(fileExt),dtype=self.byteOrder+numberFormat)
		# self.vprint('Raw data length: '+str(len(data)),1)
		nDataValuesPerPoint = len(isComplex)
		nRealsPerPoint = sum(isComplex+1)
		#remove non-existent axes since python doesn't start indices at 1 and get prod for total n points
		nPoints =np.prod([dim for dim in self.dimensions if dim >0])
		# find total number of values
		N = nRealsPerPoint*nPoints
		# reshape data into an array of with n rows 
		data = np.reshape(data,(nRealsPerPoint,nPoints)).flatten()
		self.vprint('Loaded {0} values with shape {1}'.format(N,data.shape),1)
		for k in range(nDataValuesPerPoint):
			if isComplex[k]:
				data[k,:] = np.complex64(data[k,:],data[k+1,:]) #@TODO test with actual complex data
				np.delete(data,k+1,axis=1)
		self.vprint(f'Shape after complex condensation: {data.shape}',1)
		return data

	def scaleData(self) -> np.ndarray|None:
		'''
		Scales data, by the following parameters:
			* N: Number of scans
			* G: Reciever gain (dB)
			* C: Conversion time (ms)
			* P: Microwave Power (mW)
			* T: Temperature (K)
		'''
		#check parameters
		fNames = self.Param.keys()
		prescaled = [True if ('SctNorm' in fNames) else False]
		cw = [True if self.Param['EXPT']=='CW' else False]
		self.vprint(f'Scaling data by method: {self.scaling[0]}')
		self.vprint(f'Parameters are isPrescaled: {prescaled}, isCW: {cw}',1)
		if 'N' in self.scaling: # scale by number of scans
			if 'AVGS' not in fNames:
				raise KeyError ('Missing AVGS field in the DSC file.')
			nAverages = self.Param['AVGS'][0]
			if prescaled:
				raise RuntimeError (f'Scaling by number of scans not possible, as DSC/DTA data are are already averaged over {nAverages} scans.')
			else:
				data = self.Spec/nAverages
				return data
		if cw:
			if 'G' in self.scaling: # scale by reciever gain (dB)
				if 'RCAG' not in fNames:
					raise KeyError('Cannot scale by receiver gain, since RCAG is not in the parameter file.')
				else:
					ReceiverGaindB = self.Param['RCAG'][0]
					# Xenon (according to Feb 2011 manual) uses 20*10^(RCAG/20)
					ReceiverGain = 10^(ReceiverGaindB/20);
					data = self.Spec/ReceiverGain;
					return data

			elif 'C' in self.scaling: # scale by conversion time (in seconds)
				if 'SPTP' not in fNames:
					raise KeyError('Cannot scale by conversion time, since SPTP is not in the parameter file.');
				# Xenon (according to Feb 2011 manual) already scaled data by ConvTime if
				# normalization is specified (SctNorm=True). Question: which units are used?
				# Xepr (2.6b.2) scales by conversion time even if data normalization is
				# switched off!
				ConversionTime = self.Param['SPTP'][0] # in seconds
				ConversionTime = ConversionTime*1000; # s -> ms
				data = self.Spec/ConversionTime;
				return data

			elif 'P' in self.scaling: # scale by MW power (watts)
				if 'MWPW' not in fNames:
					raise KeyError('Cannot scale by microwave power, since MWPW is not in the parameter file.')
				mwPower = self.Param['MWPW'][0]*1000 # in milliwatts
				data = self.Spec/np.sqrt(mwPower)
				return data

			elif 'T' in self.scaling: # scale by temp (K)
				if 'STMP' not in fNames:
					raise KeyError('Cannot scale by temperature, since STMP is not in the parameter file.')
				Temperature = self.Param['STMP'][0]
				data = self.Spec*Temperature
				return data
			else:
				raise RuntimeError(f'The method \"{self.scaling}\" only applies to CW data.')
