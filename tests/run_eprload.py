#! python3
from pathlib import Path
import sys
import os
# sys.path.append(os.path.join(os.path.dirname(__file__), os.pardir)) # allows importing files from parent directory
# from EPRsim_src.eprload import *
import eprload as el
# from pprint impor`t pprint as pp

eprFilePath = Path("eprfiles")

# print(os.getcwd())
# brukertest = el.eprload(eprFilePath/'99090211.dta')
brukertest = el.eprload(eprFilePath/'bes3t'/'99090211.dsc',verbose=True)
# brukertest.show_params()
