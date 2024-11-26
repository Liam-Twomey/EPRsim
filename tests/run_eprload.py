from pathlib import Path
import sys
import os
sys.path.append(os.path.join(os.path.dirname(__file__), os.pardir)) # allows importing files from parent directory
from EPRsim_src.eprload import *

eprFilePath = Path("tests/eprfiles")
print(os.getcwd())
brukertest = eprload(eprFilePath/'99090211.dta')
brukertest = eprload(eprFilePath/'99090211.dsc')
jeoltest = eprload(eprFilePath/'jeol/C60')