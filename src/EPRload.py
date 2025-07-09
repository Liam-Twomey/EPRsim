#! python3
# Author: Liam Twomey

from pathlib import Path,PurePath
import re
import numpy as np
from warnings import warn
from pprint import pprint, pformat
from types import FunctionType
'''
All logic was adapted from Stefan Stoll's EasySpin package for MATLAB, some was reimplemented
from scratch based on the Bruker BES3T specification.

'''
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

		Attributes
		----------
		Absc: :class:`np.ndarray`
			Magnetic field abscissa(s) loaded from file, often labelled "B".
		Spec: :class:`np.ndarray`
			Signal component of the data, often labelled "S". 
		Par: :class:`dict`
			Experimental parameters read from parameter file, often called "P".
	
	'''
	def __init__(self, fileName:(str | Path), scaling:str='1',debug:bool=False,
			  keepTmp:bool = False, corrFn:FunctionType = None, baseCorr:list=None): 
		self.filePath = Path(fileName) #if already a Path this won't effect anything
		self.fileExt = self.filePath.suffix
		self.scaling = scaling.upper()
		# self.Param = {}
		self.verbose = debug
		self.vprint(f'`eprload` initalized on \"{self.filePath}\"')
		if self.fileExt.isupper():
			self.extCase = 1
		elif self.fileExt.islower():
			self.extCase = 0
		else:
			raise IOError("Please don't use mixed-case extensions!")

		self.fileType = self._checkFileType()
		if isinstance(corrFn,FunctionType) and self.fileType=='bes3t':
			self.AbscCorr = corrFn(self.Absc,self.Param['XWID'][0],self.Param['SweepTime'][0])
			if self.AbscCorr is None:
				raise RuntimeError("No data returned from field correction function.")
		if isinstance(baseCorr, list):
			self.SpecCorr = self.baseCorr(baseCorr[0],baseCorr[1])
		# Delete temporary class items for memory savings; disable with keepTmp = True.
		#if keepTmp is False:
		#	del self.extCase, self.numFormat, self.byteOrder, self.dimensions, self.isComplex

	def __repr__(self):
		'''
		Allows for printing of the eprload class-object.

		'''
		return	"Abscissa:\n"+pformat(self.Absc,width=40)+"\nSignal:\n"+pformat(self.Spec,\
			indent=4,width=40)+"\nParameters:\n"+pformat(self.Param,indent=4,width=80)+\
			"\nInternal values:"+"\n\tFilePath: "+str(self.filePath)+"\n\tFile Extension: "\
			+str(self.fileExt)+"\n\tFile Extension Case: "+str(self.extCase).replace('0','Lower')\
			.replace('1','Upper')+"\n\tVerbose: "+str(self.verbose)+"\n\tScaling: "+str(\
			self.scaling[0]).replace('1','None')+'\n\tAxis complexity is: '+ str(self.isComplex)\
			.replace('1','CPLX').replace('0','REAL')

	def __str__(self):
		'''
		Allows for pretty-printing of the object's Absc, Spec, and Param attributes.

		'''
		pfa = pformat(self.Absc)
		pfs = pformat(self.Spec)
		pfp = pformat(self.Param)
		fst = f"eprload object with:\nAbscissa (B):\n {pfa} \nSpectrum (S):\n{pfp}\n\
		Parameters (P):\n{pfs}"
		return fst

	def __getattr__(self,attr):
		"""
		Define fallbacks for noncanonical aliases to self.Absc, Spec, Param under
		EasySpin format (B,S,P)
		"""
		match attr.upper():
			case 'B':
				return self.Absc
			case 'S':
				return self.Spec
			case 'P':
				return self.Param
			case 'BCORR'|'BC':
				return self.AbscCorr
			case 'SCORR'|'SC':
				return self.SpecCorr
			case _:
				raise AttributeError(f"{attr} is not an attibute of this EPRload object.")
		return self.Absc

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

	def _checkFileType(self):
		'''
		Checks through the list of known EPR filetypes. Calls the appropriate file read method
		if implemented; throws an error if it is not implemented or not a known EPR file.
		
		Parameters
		----------
		self.fileExt: Path
			Extension of the file path passed to eprload

		Returns
		-------
		ftype: str
			nickname of the filetype to be passed back up.

		Notes
		-----
		The following are viable filetypes to the function:

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
		no extension: not implemented
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
				ftyp = 'bes3t'
				self.vprint(f'File type is Bruker {ftyp}')
				# 1 is the default internal value
				if self.scaling.upper() in "NGPTC1":
					self.loadBES3T()
				else:
					raise RuntimeError("Invalid scaling option supplied. Valid options are: N,P,G,T,C.")
			case '.D01':
				ftyp = 'specman'
				self.vprint(f'File type is {ftyp}')
				if self.scaling[0]!='1':
					warn("Scaling not supported for this filetype.")
				self.loadSpecMan()
			case ('.PAR'|'.SPC'|'.SPE'|'.XML'|'.ESR'|'.DAT'|'.JSON'|'.ECO'|'.PLT'|'.SPK'|'.REF'|'.D00'):
				#Cover all other file formats 
				raise NotImplementedError("File type {FileFormat} not yet implemented and does not have a roadmap. Contact maintainers if you need it implemented.")
			case _:
				# Test for JEOL file
				with open(self.filePath,"rb") as fptr:
					identity = fptr.readline(16) # first 16 bytes = 16 ascii/utf-8 characters.
				iddc = identity.decode('latin1')
				try:
					isJeol = re.search('^spin|^cAcqu|^endor|^pAcqu|^cidep|^sod|^iso|^ani',iddc).group()
					if len(isJeol) >= 0:
						self.vprint('File type is JEOL')
						if self.scaling[0]!='1':
							warn("Scaling not supported for this filetype.")
						# [Data,Abscissa,Parameters] = eprload_jeol(FileName);
				except AttributeError:
					raise NotImplementedError(f"File format {self.fileExt} not yet implemented.")

				else:
					raise NotImplementedError(f"Unsupported file extension {self.fileExt}")
		return ftyp

#### LOAD BRUKER BE3T ####
	def loadBES3T(self) -> None:
		'''
		File loading for BES3T (Bruker EPR Standard for Spectrum Storage and Transfer)
		format. Interpreter based on BES3T version 1.2 (Xepr 2.1).
	
		Parameters
		----------
		self.filePath
		self.extCase

		
		Returns
		-------
		None. Sets the self.Absc, self.Spec, and self.Signal attributes.

		Notes
		-----

		Instuments:
			Bruker ELEXSYS and EMX spectrometers.
		Extensions: (.DSC,.DTA)
			* .DSC: descriptor file
			* .DTA: data file

		'''
		self.Param = self._readDSC()
	#	self.Spec = self._load

	
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
				byteOrder = '>'
			elif self.Param["BSEQ"][0]=="LIT":
				byteOrder = '<'
			else:
				byteOrder = ''
				"Unknown byte order specified by BSEQ keyword of {filename}"
		else:
			warn('Keyword BSEQ not found in .DSC file! Assuming BSEQ=LIT.')
			byteOrder = '<'
		self.byteOrder = byteOrder		
		# Data format tag of BES3T is IRFMT for the real part and IIFMT
		# for the imaginary part.
		if "IRFMT" in self.Param:# IRFMT: Item Real Format
			if len(self.Param["IRFMT"]) != nValsPerPoint:
				raise ValueError('Problem in BES3T DSC file: inconsistent IKKF and IRFMT fields.')
			#for char in self.Param["IRFMT"][0].split(','): # this was just assigning
			# based on the last one, might as well cut out the loop.
			match self.Param["IRFMT"][0][0].upper():
				# assign datatype of component data. Numpy naming uses byte len not bit.
				case "C": 
					numberFormat = 'i1' # int8
				case "S":
					numberFormat = 'i2' # 'int16'
				case "I":
					numberFormat = 'i4' # 'int32'
				case "F":
					numberFormat = 'f4' # 'float32'
				case "D":
					numberFormat = 'f8' # 'float64'
				case "A":
					raise NotImplementedError("Cannot read BES3T data in ASCII format.")
				case "0"|"N":
					numberFormat = ''
					raise ValueError("No BES3T data found.")
				case _:
					numberFormat = ''
					raise ValueError("Unknown value for IRFMT in .DSC file.")
		else:
			numberFormat = ''
			raise ValueError("Keyword IRFMT not found in .DSC file, cannot determine number format.")

		self.numFormat = numberFormat
		# We enforce IRFMT and IIFMT to be identical.
		if "IIFMT" in self.Param:# IIFMT: Item Imaginary Format
			if self.Param["IIFMT"] != self.Param["IRFMT"]:
				raise ValueError("IRFMT and IIFMT in .DSC file are not equal.")

#	def BES3TAbscLoad(self):
		self.vprint('Reading abscissas into self.Absc')
		# setup a 
		#axisNames = []
		#for ax in ['X','Y','Z']:
		#	axisNames.append(ax) if self.Param[f'{ax}TYP'][0] != "NODATA" else None
		#print(f"Axes in use: {axisNames}")
		axisNames = ['X','Y','Z']
		# I think I can rewrite this for loop in a way that makes more sense in python @TODO
		self.Absc = np.full((len(axisNames), max(self.dimensions)), np.nan)
		self.vprint(f'Abscissa size:{self.Absc.shape}',1)
		for a in range(len(axisNames)):
			axName = axisNames[a]
			axisType = self.Param[f"{axName}TYP"][0]
			if axisType == "NODATA":
				continue
			else:
				self.vprint(f'Reading {axName} asbscissa of type {axisType}',1)
				if axisType=="IGD":
					#nonlinear axis, try to load companion file (.xgf, .ygf, .zgf)
					tmp_nla = self._readIGD(a, axisType) 
					if type(tmp_nla) is np.ndarray:
						self.Absc[a,:] = tmp_nla
					elif type(tmp_nla) is str: 
						axisType = ''
					else:
						raise RuntimeError("Undefined return from readNonLinearAbscissa")
					del tmp_nla
				if axisType=="IDX": # not elif to allow for nonlinearabscissa error handling 2l above
					absc_min = self.Param[f"{axName}MIN"][0]
					absc_width = self.Param[f"{axName}WID"][0]
					if absc_width == 0:
						warn(f"{axName} has a width of 0.")
						absc_min = 1
						absc_width = self.dimensions[a]-1
					self.Absc[a,:] = np.linspace(start=absc_min,
								  				 stop=absc_min+absc_width,
								  				 num=self.dimensions[a])
					self.vprint(f'{axName} abscissa size: {self.Absc[a].shape}',1)
				elif axisType=="NTUP":
					raise NotImplementedError('Cannot read data with NTUP axes.')
				else:
					raise FileNotFoundError(f'AxisType is not defined for axis {axisNames[a]}')
		# flatten array to minimal dimensionality, and drop all nan values
		tmpabsc = np.squeeze(self.Absc)
		self.Absc = tmpabsc[~np.isnan(tmpabsc)]/10
		del tmpabsc
	#def BES3TSpecLoad(self):
		# get data from .dta file
		dtaExt = ['.DTA' if self.extCase else '.dta'][0]
		self.vprint(f'Reading data from {dtaExt} file in {self.numFormat} format to self.Spec')
		self.Spec, self.auxSpec = self._readDTA(dtaExt,self.byteOrder,self.numFormat,self.isComplex)
		if self.scaling[0] != '1':
			self.Spec = self.scaleData()

	def _readDSC(self) -> dict:
		"""
		Loads the contents of a Bruker DSC file into self.Param

		"""
		param = {}
		dscExt = ['.DSC' if self.extCase else '.dsc'][0]
		self.vprint(f"Reading parameters from {dscExt} file to self.Param")
		with open(self.filePath.with_suffix(dscExt)) as tmpf:
			lines = tmpf.readlines()
			for line in lines:
				# ignore commented lines
				if (not re.match(r"[*,#]",line)) and (len(line) > 0):
					# trim multiple spaces and split into list
					tmpln = re.sub(r"\s+",' ',line)
					tmpln = re.sub(r",",' ',line)
					tmpln = tmpln.split()
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
						param[tmpln[0]] = tmpln[1:]
			#param['IKKF'] = param['IKKF'][0].split(',')
		return param


	def _readIGD(self, a, axisNames) -> np.ndarray| str:
		'''
		This function parses the format of the accessory .XGF, .YGF, and .ZGF files
		for axisType IGD.
		
		Parameters
		----------
		a: int 
			Index of the axis passed.
		axisNames: list
			List of axis names for the spectrum
		self.filePath
		self.byteOrder
		self.numFormat

		Returns
		-------
		absc | np.ndarray (str if fails)
			The nonlinear abscissa of the spectrum.
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
			return tmpAbsc
		except:
			warn(f'Could not read companion file {companionFile} for nonlinear axis. Assuming linear axis.')
			axisType = 'IDX'
			return axisType

	def _readDTA(self,fileExt,byteOrder,numberFormat,isComplex) -> tuple:
		'''
		Private method for loadBES3T. Reads in the information from a .DTA file.
			
		Parameters
		----------
		fileExt:str
			File extension of the .DTA file to be read (i.e. .DTA or .dta)
		byteOrder: str
			String denoting the endianness of the data to be read, either > for big-endian, and <for little-endian.
		numberFormat: str
			Numpy-formatted number format string, i.e. i4 for int32, f8 for float64. 
		isComplex: list
			List of bool, one item per experiment dimension. True if the axis is complex, False if not.

		Returns
		-------
		
		data: np.ndarray
			A :code:`nx*ny*nz` array of the data from each 
		
		Notes
		-----
		#. Open file ``self.fileName.fileExt`` using ``np.readfile`` with
		   format byteOrder+numberFormat. Real and imaginary data are interspersed,
		   so it just reads everything into an ``(nx*ny*nz,)`` matrix (i.e. total
		   number of points in all dimensions into a list).
		#. We determine if the data is saved as result sets (i.e. if it has multiple
		   values for IKKF, and so has auxilliary data.
		#. We reshape the data to ``( nx*ny*nz , sum(isComplex+1) )`` to push each data
		   point into a row, with each column being a different attribute (one column
		   for reals, two for complexes).
		#. If the data doesn't have any aux data recorded with it (len isComplex is 1)
		   we just combine complex numbers into one signal, reshape it to ``(nx, ny, nz)``
		   for as many dimensions are defined, and return (data,None).
		#. If data has auxilliary attributes, we iterate through the columns, and select
		   columns (combining with the following column where necessary to generate
		   complex numbers, and add them to a temporary list. We reshape each member 
		   of the list to (nx,ny,nz), the first member of the list becomes the data,
		   and the rest becomes auxData. We return (data, auxdata)

		
		Data binary is composed of nPoints n tuples of data, each of which is composed
		of one data point for each axis.
		'''
		self.vprint(f"Signal complexity is: {str(isComplex).replace('1','CPLX').replace('0','REAL')}",1)
		self.vprint (f"Axis format is: {"Simple" if len(isComplex)==1 else "Result Set"}",1)
		data = np.fromfile(self.filePath.with_suffix(fileExt),dtype=byteOrder+numberFormat)
		self.vprint(f'Initial data shape: {data.shape}',1)
		dim = np.trim_zeros(self.dimensions)
		totPts = np.prod(dim)
		nCol = isComplex+1
		nValsPerPt = sum(nCol)
		data = np.reshape(data,(totPts,nValsPerPt))
		if len(isComplex) == 1:
			# for normal data not saved as a result set
			data = np.reshape(data,(totPts,isComplex[0]+1))
			if isComplex[0] == 1: # check for complex data
				data = data[:,0] + 1j*data[:,1]
			data = np.reshape(data,dim)
			self.vprint(f'final shape of Spec (simple) {data.shape}',1)
			return (data, None)
		else:
			# deal with data saved as Result Sets
			assert(data.shape[1] == nValsPerPt)
			# split reshaped data into columns
			dholder = []
			for i in range(len(nCol)):
				idx = sum(nCol[:i])
				if nCol[i] ==1:
					dtmp = data[:,idx]
				else:
					dtmp = data[:,idx]+1j*data[:,idx+1]
				dtmp = np.reshape(dtmp,(dim))
				dholder.append(dtmp)
			data = dholder[0]
			auxData = dholder[1:]
			self.vprint(f'Final shape of Spec (Result Set): {data.shape}',1)
			self.vprint(f'Aux data: {len(auxData)} signals, shaped: {auxData[0].shape}',1)
			return (data, auxData)
		
	def loadSpecMan(self) -> None:
		"""
		Load data files generated by SpecMan (.d01/.exp).
		Parameters
		----------
		self.fileName
			A filename to read
		self.fileExt
			The file extension of self.fileName
		Returns
		-------
		None. sets self.Param, self.Spec.

		Notes
		-----
		* calls _readEXP and _readD01 to get file data, then data processing
		* save transient as fullfield

		References
		----------
		Based on original code by Boris Epfel and Alexey Silakov,
		via Stefan Stoll's Easyspin.

		"""
		self.Param = self._readEXP()
		self.Spec = self._readD01()
		finpar = {}
		# post-processing of parameters
		units = {'n':1E-9,'u':1E-6,'m':1E-3,'k':1E3,'M':1E6,'G':1E9}
		finpar['title'] = param['general_name']
		#while isfield() :

	
	def _readEXP(self)->dict:
		"""
		Reads experimental parameters from a SpecMan EXP file.
		
		Parameters
		----------
		self.filePath
			path to .d01 or .exp file

		Returns
		-------
		param: dict
			Experimental parameters read from the .exp file

		Notes
		-----
		SpecMan uses a format very similar to TOML for saving experimental
		parameters, with the exception that it does not quote strings.
		Adapting code from tomllib to deal with this [1]_ .

		Read file as text block.

		#. look for `[program]`; read until next [ as literal string.
		   Ignore other sections.
		#. Split all other lines at ``=``, make first item the key
		#. Split remaining part at ``;``, and first half at ``' '``; this becomes the value.

		References
		----------
		.. [1] https://github.com/python/cpython/blob/3.13/Lib/tomllib/_parser.py

		"""
		self.vprint('Reading .exp file')
		expExt = ['.EXP' if self.extCase else '.exp'][0]
		param = {}
		with open(self.filePath.with_suffix(expExt), 'r') as f:
			raw = f.read() 
			pmatch = r'\[program\]((.|\n+)*?)(?=\[)'
			param['program'] = re.search(pmatch,raw).group(0)
			proc = re.sub(pmatch,'',raw)
			#raw = re.sub(r'\[(.*?)\]','',raw)
			proc = proc.split('\n')
			proc = [x.strip() for x in proc if x.strip()]
			param, prefix = {},None
			for line in proc:
				if line[0].startswith(r'['):
					prefix = re.sub(r'\[|\]','',line)
					print('PRE',prefix)
				else:
					tmpln=line.replace(' = ','=')
					tmpln = tmpln.split('=')
					tmpval = re.split(' ',tmpln[1],maxsplit=1) if len(tmpln) >=2 else '' 
					if tmpval[0].isdigit():
						# convert viable str to int
						tmpval[0] = int(tmpval[0])
					elif re.match(r'[0-9]+[.]?[0-9]+',tmpval[0]):
						# convert viable str to float
						tmpval[0] = float(tmpval[0])
					param[prefix+'_'+tmpln[0]] = tmpval
		print(param) 
		self.vprint('Finsished loading .exp file.')
		return param

	def _readD01(self) -> np.ndarray:
		"""
		Private function for use by `eprload.loadSpecMan`

		Parameters
		----------
		self.filePath
			Path to .d01 or .exp file

		Returns
		-------
		data: np.ndarray
			Array of data read from .d01 file.

		Notes
		-----
		File ``fid`` is little endian. Structure: ``nsig, dformat,
		{nstrm, strmdim,strmtot}`` nsig times, ``tmpdat``

		nSig: <u4
			Number of signals in the file (i.e. real, imaginary)
		dFormat: uint32
			Data format. 1=float32, other=float64
		nStrm: <i4
			Number of data streams per signal.
		strmDim: list of 4 <i4
			Number of points per stream, to a max of 4 streams. 
		strmTot: <i4
			Total number of data streams from this signal 
		tmpDat: np.ndarray
			Remainder of file. sum(strmTot) points in format dformat. Each signal
			is an array of shape strmDim
	
		**Method**

		* Read nSig and dFormat
		* read next 6*nSig bytes as nsig sets of {nstrm, strmdim[1:4],strmtot} as ndarray.
		* Read remainder of file into (1,) ndarray
		* split data into wither a 1d or complex spectrum based on nSig
		* If Nsig != [1|2], then try to figure out if there's a blank dimension
		  at the beginning of the file and try again
		* If all else fails, just read it in as a 1d file.

		"""
		self.vprint('Reading .d01 file')
		dExt = ['.D01' if self.extCase else '.d01'][0]
		dfile = self.filePath.with_suffix(dExt)
		if not hasattr(self,'Param'): # check that parameters have been loaded
			raise RuntimeError("Parameters have not been loaded.")
		dOffset = 0 # Number of bytes read so far
		sigMeta= np.fromfile(dfile,count=2,dtype='<u4')
		dOffset += 2*4 
		nSig  = int(sigMeta[0])
		dformat = '<f4' if (int(sigMeta[1]==1)) else '<f8'

		self.vprint(f'Signal metadata: {sigMeta}')
		self.vprint(f'nSig: {nSig}, dformat: {dformat}',1)
		# After 8-byte header, read ndim1 rows of 6 items as int32 (ndim1*24 bytes)
		strmMeta= np.fromfile(dfile,offset=dOffset,count=nSig*6,
				dtype='<i4').reshape(nSig,6)
		dOffset += nSig*6*4
		#self.vprint(f'Stream metadata: {strmMeta}')
		strms = []
		ntotal = 0
		# the columns of strmMeta are: [nstrm, strmDim1, strmDim2, strmDim3,
		# strmDim4, strmTot] where p denotes points.
		# Now add a column for first, the index of the first data point in the stream.
		colExt	= np.empty(nSig, dtype=int)
		strmMeta = np.c_[strmMeta, colExt]
		for spec in range(nSig):
			# if read correctly, product of nonzero dimensions == total
			nstrm = np.count_nonzero(strmMeta[spec,1:5])
			calcTotal= np.prod(strmMeta[spec,1:1+nstrm])
			if calcTotal != strmMeta[spec,5]:
				raise RuntimeError("Data not loaded correctly, calculated  and\
				specified number of data points differ.")
			strmMeta[spec,6] = ntotal
			ntotal += calcTotal
		# now strmMeta is:
		strMetaLabel = np.array([['nStrm', 'strmD1', 'strmD2', 'strmD3','strmD4',
				  'strmTot', 'startIdx']])
		self.vprint(f'Stream metadata:\n{np.r_[strMetaLabel,strmMeta]}')
		# read spectra into one long ndarray 
		rawData = np.fromfile(dfile, offset = dOffset, dtype=dformat)
		self.vprint('data shape:')
		self.vprint(f'raw: {rawData.shape}',1)
		#specs = []
		#for i in sigMeta:
		#	specs.append(rawData[i[6]:i[6]+i[5]])
		#self.vprint([j.shape for j in specs],1)

		pdim0 = np.trim_zeros(strmMeta[0,1:4])
		match nSig:
			case 0:
				raise RuntimeError("No data present.")
			case 1:
				data = np.reshape(rawData,pdim0)
			case 2:
				# rawData slice from 0th signal first to first+total,
				# then take second stream (same slice) as imaginary
				end0 = strmMeta[0,6]+strmMeta[0,5]
				data = rawData[strmMeta[0,6],end0] + 1j*rawData[end0:end0+strmMeta[1,5]]
				# reshape into transpose of first data stream dimension
				data = np.reshape(data,pdim0)
			case _:
				# if all the dimensions are the same size in all signals, and the
				# number of dimensions in signal 0 is 0, then look at the number of
				# following dimensions. If there were 3 total, read signals 1 and 2
				# as 2d. Else, read all data in as 1d, assuming an error.
				pdim1 = np.trim_zeros(strmMeta[1,1:4])
				dimsSame = strmMeta[:,1:4] == strmMeta[0,1:4]
				if dimsSame and (strmMeta.size[0] == 2):
						# assume nSig is actually 2 and nSig got corrupted
						data = np.reshape(rawData,pdim1)
				elif dimsSame and (strmMeta.size[0] == 3):
						# assume nSig is actually 3 and nSig got corrupted
						end1 = strmMeta[1,6]+strmMeta[1,5]
						data = rawData[strmMeta[1,6],end1] + 1j*rawData[end1:end1+strmMeta[2,5]]

				else:
					warn('Data format cannot be recognized, defaulting to 1d spectrum')
					data = rawData
					
		self.vprint(f'final: {data.shape}',1)
		self.vprint("Finished loading .d01 file.")
		return data

	def scaleData(self) -> np.ndarray|None:
		'''
			A class method of eprload(), called when a scaling factor is passed
			to the constructor in the *scaling* argument.
			
			Notes
			-----

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

	def baseCorr(self, dim:int=0, order:int=1):
		#try:
		#	fit = np.polyfit(self.Absc, self.S[dim,:], deg=order)
		#	bcorr = self.S[dim,:]-fit
		#except IndexError:
		#	print(self.Absc.shape)
		#	print(self.Spec.shape)
		#	fit = np.polyfit(self.Absc, self.S,deg=order)
		#	fitpoly = np.poly1d(fit)
		#	print(fit.shape)
		#	bcorr = self.S-fitpoly
		#return bcorr
		return

