import graph_tool as gt
import glob
import pandas as pd
import numpy as np

## Script to extract all data from network files located in savepath.
## Saves all data in appropriate .csv files for use in the Figure2.py script
networks = ['smallworld','mhk']

nodes = [240]
# degrees = [4,6,8,10,12,14,16,18,20,22,24,26,28,30,32]

b = {'ID':[],'CC':[],'T':[],'p':[],'SP':[],'E_glob':[]}
d = {'ID':[],'freq':[],'p':[]}

#Don't make permanent changes
save = False

for network in networks:
    for node in nodes:
        for k in degrees:
            print(k)
            data_change = False
            savepath = f'networks/{node}/{k}/{network}/'
            graphs = glob.glob(f'networks/{savepath}*')
            data_file =  f'networks/{node}/{k}/{network}_corr_gains.csv'
            data_data_file =  f'networks/{node}/{k}/{network}_props.csv'
            network_gains = pd.DataFrame(data=d)
            network_gains.set_index(['ID','freq','p'],inplace=True)
           
            network_props = pd.DataFrame(data=b)
            network_props.set_index(['ID','p'],inplace=True)

            for graph in graphs:
                print(graph)
                change = False
                G = gt.load_graph(graph)
                my_gain_data = pd.DataFrame(data=d)
                my_gain_data.set_index(['ID','freq','p'],inplace=True)

                prop = G.new_graph_property('double',0)
                if network == 'caveman' or network == 'connected_caveman':
                    G.gp.probability = prop

                # if G.vp.get('gains',False) and G.vp.get('local_clustering',False) and G.gp.get('transitivity',False) and G.gp.get('shortest_path',False) and G.gp.get('efficiency',False): 

                if G.vp.get('gains',False) and G.vp.get('local_clustering',False) and G.gp.get('transitivity',False) and G.gp.get('shortest_path',False):
                    data_change = True
                    # my_data.set_index(['ID','freq','p'],inplace=True)
                    ## SPecial for PRL fig 1

                    for f,g in zip(G.graph_properties['frequencies'],np.mean(G.vp['gains'].get_2d_array(range(len(G.graph_properties['frequencies']))),1)): 
                        my_gain_data.loc[(G.graph_properties['ID'],f,G.gp['probability']),'H2'] = g 

                    network_props.loc[(G.gp['ID'],G.gp['probability']),'CC'] = sum((G.vp.local_clustering.get_array()))/len(G.get_vertices())

                    network_props.loc[(G.gp['ID'],G.gp['probability']),'T'] = G.gp.transitivity 

                    network_props.loc[(G.gp['ID'],G.gp['probability']),'SP'] = G.gp.get('shortest_path')  

                    network_props.loc[(G.gp['ID'],G.gp['probability']),'l2'] = np.sort(G.vp.eig_laplacian.a)[1]

                    network_props.loc[(G.gp['ID'],G.gp['probability']),'lmax_l2'] =  np.max(G.vp.eig_laplacian.a) / np.sort(G.vp.eig_laplacian.a)[1]

                    network_props.loc[(G.gp['ID'],G.gp['probability']),'Rg'] =  node*np.sum(1/np.sort(G.vp.eig_laplacian.a)[1:])

                    network_props.loc[(G.gp['ID'],G.gp['probability']),'l_norm'] = np.log(np.sum(np.exp(G.vp.eig_adjacency.a))/G.num_vertices())
                
                    network_gains = network_gains.append(my_gain_data)
                    
                else:
                    print('miss')

            if data_change:
                if save:
                    network_gains.to_csv(data_file,sep='\t',mode='w',header=True)
                    network_props.to_csv(data_data_file,sep='\t',mode='w',header=True)
