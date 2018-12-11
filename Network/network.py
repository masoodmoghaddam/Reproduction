import openpnm as op
import random
import numpy as np
import scipy as sp

#------------------------------------------------ Functions --------------------------------------
def Merge_pores(network, pores, Neighbor,Cordinate, trimradius, labels=['merged']):
    # Assert that `pores` is list of lists
    try:
        len(pores[0])
    except (TypeError, IndexError):
        pores = [pores]        
    N = len(pores)
    NBs, XYZs = [], []
    for Ps in pores:
        NBs.append(Neighbor)
        XYZs.append(Cordinate)
    op.topotools.extend(network, pore_coords=XYZs, labels=labels)
    Pnew = network.Ps[-N::]
    pores_set = [set(items) for items in pores]
    NBs_set = [set(items) for items in NBs]
    ps1, ps2 = [], []    
    from itertools import combinations
    for i, j in combinations(range(N), 2):
        if not NBs_set[i].isdisjoint(pores_set[j]):
            ps1.append([network.Ps[-N+i]])
            ps2.append([network.Ps[-N+j]])
    # Add (possible) connections between the new pores
    op.topotools.connect_pores(network, pores1=ps1, pores2=ps2, labels=labels)
    # Add connections between the new pores and the rest of the network
    op.topotools.connect_pores(network, pores2=sp.split(Pnew, N), pores1=NBs, labels=labels)
    # Trim merged pores from the network
    op.topotools.trim(network=network, pores=sp.concatenate(pores))
    Ps = network.find_nearby_pores(pores=network.Ps[-len(pores)::], r=trimradius, flatten=True)
    op.topotools.trim(network=network, pores=Ps)
def MacroProsity(prosity):
    Nnp=int((1-prosity)*Nni)
    radiuos=PSR*averageporeaize  /1.8 #true value foe radius is /2 
    Ps = pn.find_nearby_pores(pores=random.choice(pn.pores(labels=['surface'], mode='not')), r=PtoP*2, flatten=True)
    op.topotools.merge_pores(network=pn, pores=Ps, labels=['merged'])
    while len(pn.pores(labels=['merged'], mode='not'))>Nnp:
        Point=random.choice(pn.pores(labels=['surface','merged'], mode='not'))
        pores = pn.find_nearby_pores(pores=Point, r=radiuos, flatten=True)
        pores=[i for i in pores if i not in pn.pores(labels=['surface','merged'])]
        if len(pores)>0:
            neighbor=pn.find_neighbor_pores(pores=pores, mode='union', flatten=True, include_input=False)
            cord=pn['pore.coords'][neighbor].mean(axis=0)
            Merge_pores(network=pn, pores=pores, Neighbor=neighbor,Cordinate=cord, trimradius=PtoP , labels=['merged'])
#            print('************ ',1- len(pn.pores(labels=['merged'], mode='not'))/Nni)
            D=(1- len(pn.pores(labels=['merged'], mode='not'))/Nni)/prosity*50
            print('[','='*int(D), '>','.'*(50-int(D)),']')
    # Deleting Pores which placed very close (less than initial pore to pore distance of the grid) to merged porse         
    op.topotools.trim(network=pn, pores=pn.find_nearby_pores(pores=pn.pores(labels=['merged']), r=PtoP, flatten=True))
    # Deleting pores with just one throat
    op.topotools.trim(network=pn, pores=[sp.where(pn.num_neighbors(pn.pores()) == 1)])
#------------------------------------------------ contorol parameters --------------------------------------
dp=40e-6                   #particle diameter
PtoP=10e-8                 #pore to pore espacing   My idea  =dp/Numberofpores
Numberofpores=int(dp/PtoP)   #number of pores on radius
phi=0.5
PSR=100                     # pore size ratio
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
MacroProsity(prosity=phi)

