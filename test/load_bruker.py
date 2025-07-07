#! python3
from pathlib import Path
from eprsim.EPRload import *

eprFilePath = Path("eprfiles/bes3t")

# brukertest = el.eprload(eprFilePath/'99090211.dta')
#brukertest = eprload(eprFilePath/'99090211.dsc',debug=True)
#print(type(brukertest.Absc),type(brukertest.Spec),type(brukertest.Param))
# brukertest.show_params()
bt2 = eprload(eprFilePath/'00012107.dsc', debug=True)
