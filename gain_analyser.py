import graph_tool as gt
import glob
import pandas as pd
import os
import numpy as np
from scipy import linalg
# Script for calculating gains on networks, also able to append the list of frequencies the gain is found at
networks = ['mhk','smallworld']
nodes = [240]
## Degrees for networks to be analyzed
degrees = np.arange(4,32,2)
d = {'ID':[],'freq':[],'p':[],'CC':[]}

w = np.logspace(-4,1,21)
increase = False
## Lists of frequencies used
increase_w = np.logspace( np.log10(w[8]),np.log10(w[10]),10)[1:-1]
w = np.append(w,increase_w)
increase_w = np.logspace( np.log10(w[23]),np.log10(w[24]),4)[1:-1]
w = np.append(w,increase_w)
increase_w = np.logspace( np.log10(w[23]),np.log10(w[29]),10)[1:-1]

def get_gain(graph,w,N):
    L = gt.laplacian(graph,normalized=False) #the build in normalized gives the symetric normalized laplacian, but we want the random walk normalized laplacian
    
    L = (L/L.diagonal()).T  ## Random walk normalization  D^-1 L = LD^-1 because L is symetric

    h2 = G.new_vertex_property('vector<double>')
    for g in G.vertices():
        ida = np.arange(N) != g
        idb = np.arange(N) == g
        A = L[np.ix_(ida,ida)].astype(complex)
        B = L[np.ix_(ida,idb)]

        H2 = []    
        for f in w:
            np.fill_diagonal(A,1.0+1j*f)
            #A.setdiag(f*1j-1)
            h = linalg.solve(A,-B)
            
            H2.append(linalg.norm(h)**2)

        h2[g] = H2

    return h2

for network in networks:
    for node in nodes:
        for k in degrees:
            print(k)
            database_change = False
            savepath = f'networks/{node}/{k}/{network}/'
            graphs = glob.glob(f'{savepath}*')
            data_file =  f'networks/{node}/{k}/network_gains.csv'
            file_exist = False
            if file_exist:
                network_gains = pd.read_csv(data_file,sep='\t')
                network_gains.set_index(['ID','freq','p','CC'],inplace=True)
                IDs = network_gains.index.values
            else:
                network_gains = pd.DataFrame(data=d)
                network_gains.set_index(['ID','freq','p','CC'],inplace=True)

            for graph in graphs:
                print(graph)
                change = False
                G = gt.load_graph(graph)
                # The following line is to update the gain of the network
                if not G.vertex_properties.get('gains',False):
                    frequencies = G.new_graph_property('vector<double>',val=w)
                    G.graph_properties['frequencies'] = frequencies

                    G.vertex_properties['gains'] = get_gain(G,np.append(w,increase_w),node)

                if increase:
                    if True:
                        if len(G.gp.frequencies) < len(np.append(w,increase_w)):
                            gains = G.vp['gains'].get_2d_array(range(len(w)))
                            frequencies = G.new_graph_property('vector<double>',val=np.append(w,increase_w))
                            G.graph_properties['frequencies'] = frequencies

                            new_gains = np.append(gains,get_gain(G,increase_w,node).get_2d_array(range(len(increase_w))),axis=0)
                            G.vp.gains = G.new_vertex_property('vector<double>',vals=new_gains.T)
                            # G.save(graph)
