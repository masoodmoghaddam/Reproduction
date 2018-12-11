import openpnm as op
import random
import numpy as np
import scipy as sp
import sys  
sys.path.append("/home/masood/Pore Network Model/Network")  
from FastMerge import Merge_pores  

def MacroProsity(network, porosity, PoreSizeRatio,averageporeaize,Initialnanopores,trimradius):
    Nnp=int((1-porosity)*Initialnanopores)
    radiuos=PoreSizeRatio*averageporeaize  /1.8 #true value foe radius is /2 
    Ps = network.find_nearby_pores(pores=random.choice(network.pores(labels=['surface'], mode='not')), r=trimradius*2, flatten=True)
    op.topotools.merge_pores(network=network, pores=Ps, labels=['merged'])
    while len(network.pores(labels=['merged'], mode='not'))>Nnp:
        Point=random.choice(network.pores(labels=['surface','merged'], mode='not'))
        pores = network.find_nearby_pores(pores=Point, r=radiuos, flatten=True)
        pores=[i for i in pores if i not in network.pores(labels=['surface','merged'])]
        if len(pores)>0:
            neighbor=network.find_neighbor_pores(pores=pores, mode='union', flatten=True, include_input=False)
            cord=network['pore.coords'][neighbor].mean(axis=0)
            Merge_pores(network=network, pores=pores, Neighbor=neighbor,Cordinate=cord, trimradius=trimradius , labels=['merged'])
            D=(1- len(network.pores(labels=['merged'], mode='not'))/Initialnanopores)/porosity*50
            print('[','='*int(D), '>','.'*(50-int(D)),']')
    # Deleting Pores which placed very close (less than initial pore to pore distance of the grid) to merged porse         
    op.topotools.trim(network=network, pores=network.find_nearby_pores(pores=network.pores(labels=['merged']), r=trimradius, flatten=True))
    # Deleting pores with just one throat
    op.topotools.trim(network=network, pores=[sp.where(network.num_neighbors(network.pores()) == 1)])