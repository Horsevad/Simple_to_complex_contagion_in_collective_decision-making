import numpy as np
import networkx as nx
import scipy as sp
from pathlib import Path
from scipy import sparse
import graph_tool as gt
import graph_tool.clustering
import graph_tool.spectral
import graph_tool.topology

## Remove these later

def ws_network(N,k,p,seed=None):
    """
    Function for creating a Modified Holme-Kim network
    Takes inputs:
       N: int, Number of nodes
       k: integer, The mean degree of nodes
       p: Probability of rewireing of edges on the graph
       seed: Integer for the numpy random seed function
    Returnd:
       G: a graphtool graph

    """ 
    if seed is None:
        seed=np.random.randint(2**63)

    G = gt.Graph(directed=False)
    G.add_edge_list(np.transpose(sp.sparse.tril(nx.adjacency_matrix(nx.connected_watts_strogatz_graph(N,k,p,tries=1000000,seed=seed))).nonzero()))

    return G

def linear_threshold_model(G,threshold,seed_nodes=None,init_spread=True,max_iter=None):
    """
    Runs the linear threshold model on a grap_too
l graph G. S_i(t+1) = {1 if <s_j>i~j > T otherwise 0. If the average state of neighbours of i is more than the threshold T, switch the state from 0 to 1

    Takes inputs:
       G: Graph_tool graph
       threshold: float or list of the thresholds as a fraction of neighbors that need to be active to transmit
       seed_nodes: int or list of nodes to start infection from, if None a single random agent is choosen. If int, a number of random nodes is choosen, if a list of nodes those are the initial seeds.
       init_spread: Bool, Wether to spread the infection to the neighbors of seed nodes
       max_iter = maximum number of iterations before stopping. If all nodes in the graph are active the model will stop
    Returnd:
       infected_step: Graph_tool VertexPropertyMap with an interger denoting the step of activation for a given vertex, None if never activated. Velues are stored as vectors, with each index corresponding to a given threshold.
       seed_nodes: Graph_tool GraphProperty with the initail seed nodes for the model run

    """

    if seed_nodes == None:
        [seed_nodes for x in np.random.choice(G.get_vertices(),1)]

    if not type(seed_nodes) is list:
        seed_nodes = np.random.choice(G.get_vertices(),seed_nodes)

    if max_iter is None:
        max_iter = G.num_vertices()
        
    if not type(threshold) is list:
        [threshold]

    infections = []
    degree_dist = G.get_out_degrees(G.get_vertices())

    T = np.array((graph_tool.spectral.adjacency(G).T / degree_dist).T)  

    for th in threshold:
        # Choose the initial infected nodes
        infected = np.zeros(G.num_vertices(),dtype=int)
        infection_step = np.full(G.num_vertices(),np.inf,dtype=float)
        node_list = np.arange(G.num_vertices(),dtype=int)
        
        #Infect the seed nodes
        infected[seed_nodes] = 1
        #Record seed nodes infected at t=-1
        infection_step[seed_nodes] = -1

        #Initial spread, if choosen
        if init_spread:
            infected[T.dot(infected) > 0] = 1
            infection_step[np.logical_and(infected > 0, np.isinf(infection_step))] = 0
            i = 1
        else:
            i = 0

        while (not all(infected) and (i < max_iter) and i-1 in infection_step):
            infected[T.dot(infected) >= th] = 1
            infection_step[np.logical_and(infected > 0, np.isinf(infection_step))] = i

            i += 1
        
        infected_step = G.new_vp(value_type='int',vals=infection_step)
        infections.append(infected_step)
        
    
    infected_vectormap = gt.group_vector_property(infections)
    # G.vp['infected_step'] = infected_vectormap    
    threshold_vector = G.new_gp(value_type='vector<double>',val=threshold)


    return infected_vectormap,seed_nodes,threshold_vector


def get_laplacian_eigenvalues(G):
    """
    """

    if not G.vertex_properties.get('eig_laplacian',False):
        eig_lap = np.linalg.eigvalsh(gt.spectral.laplacian(G,norm=False).todense())
        G.vp['eig_laplacian'] =  G.new_vertex_property('double',vals=eig_lap)
    return G

def get_kirchhoff_index(G):
    G = get_laplacian_eigenvalues(G)
    G.graph_properties['kirchhoff'] = G.new_graph_property('int64_t',sum(1/np.sort(G.vp.eig_laplacian.get_array())[1:]))
    return G

def get_recursive_graph_paths(root):
    paths = Path(Path(root)).rglob('*.gt')
    return paths

def get_local_clutsering(G):
    if not G.vertex_properties.get('local_clustering',False):
        G.vertex_properties['local_clustering'] = graph_tool.clustering.local_clustering(G)
    return G

def get_transitivity(G):
    if not G.gp.get('transitivity',False):
        trans = G.new_graph_property('double',val=graph_tool.clustering.global_clustering(G)[0])
        G.graph_properties['transitivity'] =  trans
    return G

def get_ave_shortest_path(G):
    if not G.gp.get('shortest_path',False):
        G.gp['shortest_path'] =  G.new_graph_property('double',val=np.sum( graph_tool.topology.shortest_distance(G).get_2d_array(range(G.num_vertices())))/(G.num_vertices()*(G.num_vertices()-1)))
    return G
