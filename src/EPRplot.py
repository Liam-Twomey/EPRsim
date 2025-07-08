import seaborn as sns
import matplotlib_inline
from matplotlib.pyplot import rcParams
matplotlib_inline.backend_inline.set_matplotlib_formats('svg')
rcParams["font.family"] = "Arial"
rcParams['font.size'] = 12
rcParams["savefig.bbox"] = 'tight'
rcParams["figure.autolayout"]=True
rcParams["lines.linewidth"]=1

pal2 = sns.color_palette(['#000000','#023eff','#ff7c00','#1ac938','#e8000b','#8b2be2','#9f4800', '#f14cc1','#ffc400','#00d7ff'])
sns.set_palette(pal2)