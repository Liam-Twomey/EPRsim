#! python3
# This code aims to reimplement the eprload functionality of EasySpin in Python 3.
#class eprSpectrum:
#    '''
#    The default internal representation for experimental data.
#    '''
#    def __init__(self):
from pathlib import Path
import pandas as pd

def eprload(filename:str, v:bool=False):
    '''
    The parent eprload function; identifies the filetype, then calls the filetype-specific 
    eprload function below.
    Supported filetypes:
        Bruker BES3T:   .DTA, .DSC
        SpecMan:        .d01, .exp
        # more to be added later.
    Args:
        filename: a filename or path to a supported datafile. If it is a multi-file format (i.e.
        BES3T, the partner file must have the same basename and be in the same directory.
    Returns:
        data: an object containing a standard Spectrum object, using the 
    '''
    not_yet_supported = "This file is supported by EasySpin, but not yet supported by this program. Please \
    contact the maintainers if you need it implemented."
    #load filename and split
    datfile = Path(filename)
    # Ensure file has only one extension, extract it.
    if len(datfile.suffixes) > 1:
        raise NotImplementedError('Compressed files or files with multiple extensions are not supported.')
    ext_to_match = datfile.suffixes[0].upper()
    # Set default value for return object
    tmpSpec = None
    # check for known file types. If filetype is in easyspin but not implemeneted here, throw unique error.
    match ext_to_match:
        case '.DSC' | '.DTA':
            tmpSpec = eprload_BrukerBES3T(datfile)
        case '.D01' | '.EXP':
            tmpSpec = eprload_specman(datfile)
        case '.SPC'|'.PAR'|'.SPE'|'.XML'|'.ESR'|'.DAT'|'.JSON'|''|'.EPR'|'.PLT'|'.ECO'|'.SPK'|'.REF'|'.d00'|'.exp':
            raise NotImplementedError(not_yet_supported)
        case '_':
            raise AssertionError('This filetype is not recognized as an EPR data file.')
    if tmpSpec is not None:
        return tmpSpec
    else:
        raise ImportError('File not imported!')

def eprload_BrukerBES3T(datafile:Path):
    '''
    The Bruker BES3T (Bruker EPR Standard for Spectrum Storage and Transfer) file format:
        .DSC: FullBaseName descriptor file
        .DTA: Binary data file:with expression as target:
            pass
    This format is used by Bruker ELEXSYS and EMX machines
    Code adapted from EasySpin, orginally based on BES3T version 1.2 (Xepr 2.1)
    '''
    dscf = datafile.with_suffix('.DSC')
    dtaf = datafile.with_suffix('.DTA')
    # now read descriptor file. Break on one or more space
    dsc = pd.read_csv(dscf,sep=r'\s+',comment='*')
    print(dsc)
        
    return

def eprload_specman(datafile:Path):
    '''
    Loads Specman data. docstring @TODO
    '''
    return

# Parse arguments
