import openpnm as op
import numpy as np
import scipy as sp

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