#! python3
# Author: Liam Twomey

from pathlib import Path,PurePath
from re import sub,match,search
import numpy as np
#from sre_constants import CHARSET
from warnings import warn
from pprint import pprint, pformat

# all logic adapted from Stefan Stoll's EasySpin.
class eprload:
	'''
		The parent eprload class acts as a wrapper for instrument-specific `eprload`
		functions. It determines the filetype of the entered file, then calls the
		appropriate method.
	
		Parameters
		----------
		filename : str or Path
			The name of (one of) the EPR data files to be loaded.
		scaling : {n, p, g, t, c}
			Data attribute to scale by, see :py:meth: `eprload.scaleData` method
		verbose:
			Request printout of debug information for this object. Debug info is implemented
			via :py:meth:`vprint`.

		Returns:
		--------
		`eprload` object with the attributes:
		
		#.. attr:: Absc
		#	Magnetic field abscissa(s) loaded from file, often labelled "B".
		#.. attr:: Spec
		#	Signal component of the data, often labelled "S". 
		#.. attr:: Par
		#	Experimental parameters read from parameter file, often called "P".
	
    '''
	def __init__(self, fileName:(str | Path), scaling:str='1',verbose:bool=False): 
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

	def __repr__(self):
		'''
		Allows for printing of the class object.

		'''
		return 	"Abscissa:\n"+pformat(self.Absc,width=40)+"\nSignal:\n"+pformat(self.Spec,indent=4,width=40)+\
			"\nParameters:\n"+pformat(self.Param,indent=4,width=80)+"\nInternal values:"+"\n\tFilePath: "+str(self.filePath)+\
			"\n\tFile Extension: "+str(self.fileExt)+"\n\tFile Extension Case: "+str(self.extCase).replace('0','Lower').replace('1','Upper')+\
			"\n\tVerbose: "+str(self.verbose)+"\n\tScaling: "+str(self.scaling[0]).replace('1','None')+'\n\tAxis complexity is: '+\
			str(self.isComplex).replace('1','CPLX').replace('0','REAL')

	def __str__(self):
		'''
		Allows for pretty-printing of the object's Absc, Spec, and Param attributes.

		'''
		pfa = pformat(self.Absc)
		pfs = pformat(self.Spec)
		pfp = pformat(self.Param)
		fst = f"eprload object with:\nAbscissa (B):\n {pfa} \nSpectrum (S):\n{pfp}\nParameters (P):\n{pfs}"
		return fst

	def getAbsc(self):
		'''
		Returns the Absc attribute of the object.

		'''
		return self.Absc

	def getSpec(self):
		'''
		Returns the Spec attribute of the object.

		'''
		return self.Spec

	def getParam(self):
		'''
		Returns the Param attribute of the object.

		'''
		return self.Param

	def esfmt(self):
		'''
		(WIP) Returns spectrum in EasySpin format: (B, Spc, P) or (Absc,Spec,Param)

		'''
		return (self.Absc, self.Spec, self.Param)
		
	def vprint(self,output:str,indentLevel:int=0):
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
		'''
		Checks through the list of known EPR filetypes. Calls the appropriate file read method
		if implemented; throws an error if it is not implemented or not a known EPR file.

		Possible calls:
		---------------

		.DTA, .DSC: loadBrukerBES3T()
			Loads data from Bruker ELEXSYS and EMX spectrometers.
		.dO1, .exp: loadSpecMan()
			Loads data collected by SpecMan from a homebuilt instrument.
		.spc, .par: not implemented 
			Load Bruker ESP and WinEPR data.
		.spe, .xml: not implemented
			Load Magnettech data.
		.ESR: not implemented
			Load Active Spectrum data.
		.dat, .json: not implemented
			Load Adani data.
		*no extension*: not implemented
			Load JEOL data.
		.epr: not implemented
			Load CIQTEK data.
		.PLT: not implemented
			Load MAGRES data.
		.eco: not implemented
			Load qese, tryscore data.
		.spk, .ref: not implemented
			Load Varian data.
		.d00, .exp: not implemented
			Load ESE data.

		If your data is in ASCII formats, use pandas.read_csv() to open it.

		'''
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
		File loading for BES3T (Bruker EPR Standard for Spectrum Storage and Transfer)
		format. Interpreter based on BES3T version 1.2 (Xepr 2.1).

		Instuments:
			Bruker ELEXSYS and EMX spectrometers.
		Extensions: (.DSC,.DTA)
			* .DSC: descriptor file
			* .DTA: data file

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
		self.isComplex = np.array(isComplex)
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
						self.numberFormat = 'i1' # int8
					case "S":
						self.numberFormat = 'i2' # 'int16'
					case "I":
						self.numberFormat = 'i4' # 'int32'
					case "F":
						self.numberFormat =  'f4' # 'float32'
					case "D":
						self.numberFormat = 'f8' # 'float64'
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
			# flatten array to minimal dimensionality, and drop all nan values
			tmpabsc = self.Absc.flatten()
			self.Absc = tmpabsc[~np.isnan(tmpabsc)]
			del tmpabsc
	#def BES3TSpecLoad(self):
		# get data from .dta file
		dtaExt = ['.DTA' if self.extCase else '.dta'][0]
		self.vprint(f'Reading data from {dtaExt} file in {self.byteOrder}{self.numberFormat} format to self.Spec')
		self.Spec = self.readBinaryDataMatrix(dtaExt,self.byteOrder,self.numberFormat,self.isComplex)
		if self.scaling[0] != '1':
			self.Spec = self.scaleData()

	def readNonlinearAbscissa(self, a, axisNames) -> None:
		'''
		This function parses the format of the accessory .XGF, .YGF, and .ZGF files
		for axisType IGD.
		
		Parameters:
		-----------

		a: int 
			Index of the axis passed.
		axisNames: list
			List of axis names for the spectrum

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

	def readBinaryDataMatrix(self,fileExt,byteOrder,numberFormat,isComplex) -> np.ndarray:
		'''
		Method called by :meth:`loadBES3T` to actually read in the information from
		a .DTA file.

		Parameters:
		-----------

		fileExt:str
			File extension of the .DTA file to be read (i.e. .DTA or .dta)
		byteOrder: str
			String denoting the endianness of the data to be read, either '>' for big-endian, and '<' for little-endian.
		numberFormat: str
			Numpy-formatted number format string, i.e. 'i4' for int32, 'f8' for float64. 
		isComplex: list
			List of bool, one item per experiment dimension. True if the axis is complex, False if not.

		Returns:
		--------
		
		data: np.ndarray
			A nx*ny*nz array of the data from each 

		Mechanism:
		----------

		#. Open file `self.fileName.fileExt` using `np.readfile` with format byteOrder+numberFormat.
			* Real and imaginary data are interspersed, so it just reads everything into an
			  nx*ny*nz matrix (i.e. total number of points in all dimensions.
		#. This data is then np.reshape()d into a :math:`(nRealsPerPoint :times (nx*ny*nz))`
		   array, so each row is a different real value index, and each column is a new
		   datapoint. 
		#. Now we cope with each datapoint having nDataValuesPerPoint values.
			* Note: both Matlab and numpy index [row, column].
			* We run through the isComplex array and compare it against the rows of the data,
			  and if a row isComplex, then it and the row after it are combined into an
			  np.complex128 complex number.
			* If the row is not complex, then it is left alone.
		#. The resulting matrix is then reshaped to an array of [[Absc X],[Absc Y],[Absc Z]].
		
		'''

		self.vprint('Axis complexity is: '+str(isComplex).replace('1','CPLX').replace('0','REAL'),1)
		data = np.fromfile(self.filePath.with_suffix(fileExt),dtype=byteOrder+numberFormat)
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
			A class method of eprload(), called when a scaling factor is passed
			to the constructor in the *scaling* argument.

			===== =============== =====  =============================
			Value Scale By        Units  Limitations
			===== =============== =====  =============================
			n     number of scans --     non-Bruker only
			g     reciever gain   dB     CW Bruker ESP only 
			c     conversion time ms     CW Bruker ESP only
			p     microwave power mW     CW Bruker only
			t     temperature     K      CW Internal temp control only
			===== =============== =====  =============================

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
