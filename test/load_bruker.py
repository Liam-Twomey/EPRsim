#! python3
from pathlib import Path
from eprsim.EPRload import *

eprFilePath = Path("eprfiles/bes3t")

# brukertest = el.eprload(eprFilePath/'99090211.dta')
#brukertest = eprload(eprFilePath/'99090211.dsc',debug=True)
#print(type(brukertest.Absc),type(brukertest.Spec),type(brukertest.Param))
# brukertest.show_params()
bt2 = eprload(eprFilePath/'00012107.dsc', debug=True)
bt3 = eprload(eprFilePath/'calib_10dB_RT_9370-9440mT_lin_1G.DSC', debug=True)
bt4 = eprload(eprFilePath/'hy3205t132b.DSC', debug=True)


