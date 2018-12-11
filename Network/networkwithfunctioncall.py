import openpnm as op
import random
import numpy as np
import scipy as sp
import sys  
sys.path.append("/home/masood/Pore Network Model/Network")  
from Hierarchical import MacroProsity  
#------------------------------------------------ Functions --------------------------------------

#------------------------------------------------ contorol parameters --------------------------------------
dp=20e-6                   #particle diameter
PtoP=10e-8                 #pore to pore espacing   My idea  =dp/Numberofpores
Numberofpores=int(dp/PtoP)   #number of pores on radius
phi=0.5
PSR=20                     # pore size ratio
averageporeaize=1.6e-8    #average pore size
#------------------------------------------------ circular pore pn generation --------------------------
pn = op.network.Cubic(shape=[Numberofpores, Numberofpores, 1], spacing=PtoP,connectivity=6)
middlepoint=pn.Np/2+Numberofpores/2
Ps = pn.find_nearby_pores(pores=(middlepoint), r=dp/2, flatten=True)
pn['pore.dummy_1']=[False] * pn.Np
pn['pore.dummy_1'][Ps]=True
pn['pore.dummy_1'][int(middlepoint)]=True
Pss=pn.pores('pore.dummy_1',mode='not')
#------------------------------------------------ finding boundary pores --------------------------------------
Pb = pn.find_nearby_pores(pores=(middlepoint), r=(dp/2-PtoP), flatten=True)
pn['pore.surface']=[True] * pn.Np
pn['pore.surface'][Pb]=False
pn['pore.surface'][int(middlepoint)]=False
#------------------------------------------------ triming forcircular pore pn generation ------------------
op.topotools.trim(network=pn, pores=Pss)
Nni=pn.Np
##------------------------------------------------ Hierarchical pn generation ------------------------------
MacroProsity(network=pn, porosity=phi, PoreSizeRatio=PSR,averageporeaize=averageporeaize,Initialnanopores=Nni,trimradius=PtoP)