#!/usr/bin/env python3
## Script for running the generating Small-world graphs, analyzing their network properties and running the Linear Threshold Model
import os
import numpy as np
import time
import graph_tool as gt
import random

from network_functions import ws_network
from network_functions import get_recursive_graph_paths
from network_functions import linear_threshold_model
from network_functions import get_local_clutsering
from network_functions import get_transitivity
from network_functions import get_laplacian_eigenvalues
from network_functions import get_kirchhoff_index
from network_functions import get_ave_shortest_path
## File used to generate netwokrs, and subsequently do the full analysis
old_graph_list = []
#Parameters for networks
network_type = 'ws'
nodes=[100]
average_degree=[16]
## Thresholds are chosen such that they approximately give the thresholds of [1,...,8]
thresholds = np.linspace(0.01,0.5,16)

## How many realizations to do of each set of parameters
desired_realizations= 1
## How many unique starting points to run the LTM from on a network
unique_network_seeds = 1
## The rewirering probabilities to be used in the watts-strogatz model
## We use from 10**-3 to 1 to get the full range of the Watts-Strogatz networks
probabilities = np.append([0], np.logspace(-3,-0,10))
## Root directory for storing networks
network_root = 'networks/'

## Make a dict with the probabilities as keys, to count the available networks
realization_counter=dict.fromkeys(probabilities,0)
##Cascade sizes of interest
cascades = np.round(np.linspace(0.1,0.9,9),1)
##Loop over entries in nodes and degrees, to create the desired networks
for N in nodes:
    for k in average_degree:
        #Seed path for saving networks, put the
        network_path = f"{network_root}{N}/{k}"
        ## If the network_path does not exist it is created
        if not os.path.exists(network_path):
            os.makedirs(network_path)

        # Count the graphs in network_path and create missing networks according to  parameters (N,k,p)
        for graph_path in get_recursive_graph_paths(network_path):
            g_old = gt.load_graph(str(graph_path))
            realization_counter[g_old.gp.probability] += 1
        #### Create missing networks ####
        for p in realization_counter:
            ## Because only one network exists in the Watts-Strogatz
            ## model for p=0 wer have a special case. Having many instances of the same network
            ## will obscure the correlation results later
            if p == 0:
                if realization_counter[p] > 0:
                    ## Makes the rest of the code not create any further p=0 networks
                    realization_counter[p] = desired_realizations
                else:
                    realization_counter[p] = desired_realizations-1
            # Too many networsk of a given probability has already been created
            # This might obscure the correlation results later
            if realization_counter[p] > desired_realizations:
                print(f'too many networks of p={p}')
            # Create new networks, until the disired number of realizations is achieved
            elif realization_counter[p] < desired_realizations:
                ## Loop until the desired number of networks is created
                while realization_counter[p] < desired_realizations:
                    G = ws_network(N,k,p,seed=None)
                    # Add relevant creation properties with the graph.
                    # Name graphs by their creation time in seconds since the epoc. This should ensure unique filenames
                    G.graph_properties['ID'] = G.new_graph_property('int64_t',val=int(time.time()*1000))
                    G.graph_properties['ntype'] = G.new_graph_property('string',val=network_type)
                    G.graph_properties['probability'] = G.new_graph_property('double',p)

                    G.save(f'{network_path}/{G.gp.ID}.gt')
                    realization_counter[p] += 1

        #### Calculate relevant network metrics and the LTM for the seeds ####
        net_count = 0
        for graph_path in list(get_recursive_graph_paths(network_path))[:]:
            G = gt.load_graph(str(graph_path))

            print(f'analysis of {G.gp.ID}')
            ## Run the network analyzer function on all the networks
            G = get_local_clutsering(G)
            G = get_transitivity(G)
            G = get_kirchhoff_index(G)
            G = get_ave_shortest_path(G)


            ## run the liner thershold model on the model, for an appropriate number of times
            print(f'Running LTM on {G.gp.ID}')
            ## Check if any seeds existes, and compute the LTM such that a sufficient number of unique_seeds was used
            old_seeds = G.gp.get('seeds',False)

            # if they exist place them in a list
            if old_seeds:

                ## Only nessesary if not sure wether seeds is properly formated, and takes a single integer seed and makes it to a list
                seeds = [[s] if type(s) is not list else s for s in old_seeds]
            else:
                seeds = []
            # Find the missing seeds
            new_picks  = unique_network_seeds - len(seeds)
            if new_picks > 0:
                if old_seeds:
                    # If old runs exist make a ocpty to be merged with ne w runs later
                    all_spreads = gt.ungroup_vector_property(G.vp['infected_step'],range(len(seeds)*len(G.gp.thresholds)))
                else:
                    all_spreads = []

                #Calculate new unique seeds different from previous seeds, if desired_seeds exceed nudes, only compute all nodes as seeds
                new_seeds = []
                if (G.num_vertices() - len(seeds)) >= new_picks:
                    new_seeds = random.sample(list(np.delete(G.get_vertices(),seeds)),new_picks)
                    print('New seeds added')
                else:
                    new_seeds = np.delete(G.get_vertices(),seeds)
                    print('Full, preforming last possible seeds')

                #Running of the LTM on the newly selected seeds
                print('#### new runs p={},k={} ###'.format(np.round(G.gp.probability,4),np.max(G.get_out_degrees(new_seeds))))
                for seed in new_seeds:
                    spreads,_,thresholds_map = linear_threshold_model(G,thresholds,seed_nodes=[seed],max_iter=N,init_spread=True)

                    all_spreads += gt.ungroup_vector_property(spreads,range(len(thresholds)))
                    seeds.append(seed)
                    # Merge the spreading vectors
                    new_spreads = gt.group_vector_property(all_spreads)
                # Add new properties ot the graph and save
                G.vp['infected_step'] = new_spreads
                G.gp['thresholds'] = thresholds_map
                G.gp['seeds'] = G.new_gp(value_type='python::object',val=seeds)

            else:
                print(f'No new seeds')

            #### Calculate polarization speeds and save them to seed vertices ###
            ## fetch if exists; othterwise create
            print('calculating polarization speeds')
            if 'polarization_speed' not in G.vp:
                polarization_speeds = np.full((G.num_vertices(),len(cascades)*len(thresholds)),-1.0)
                G.vp['polarization_speed'] = G.new_vp(value_type='vector<double>',vals=polarization_speeds)

            polarization_speeds = G.vp.polarization_speed.get_2d_array(range(len(cascades)*len(thresholds))).T

            for idy,s in enumerate(seeds):
                # Hack to get around index [0] or [0][0]
                first = polarization_speeds[s]
                while isinstance(first,np.ndarray):
                    first = first[0]
                if  first != -20:
                    # print('calc')
                    #initialize lits to aggregate the polarization speeds
                    speeds = []
                    # load simulation data
                    spread = G.vp.infected_step.get_2d_array(np.arange(idy*len(thresholds),(idy+1)*len(thresholds)))
                    for idx,th in enumerate(G.gp.thresholds):
                        cascade_sizes = cascades ## Copy where we can shorten list on iterations in while loop
                        infected = 0             ## Reset on loop start

                        ## Find nodes infected at given step
                        val,counts = np.unique(spread[idx],return_counts=True)
                        ## Refactor to fraction of nodes
                        counts = counts / G.num_vertices()
                        ## For each step in LTM add newly infected to total
                        for i,new in enumerate(counts[val>-2]):
                            infected += new
                            ## Current number of infected are used to find polarization speed if larger than cascade size, and the step is after initialization, and all cascade sizes has not been exceeded yet.
                            while len(cascade_sizes) > 0 and infected > cascade_sizes[0] and val[i] > 0:
                                speeds.append(infected/val[i])
                                cascade_sizes = cascade_sizes[1:]
                        for j in cascade_sizes:
                            speeds.append(-1)
                    polarization_speeds[s]=speeds

            G.vp['polarization_speed'] = G.new_vp(value_type='vector<double>',vals=polarization_speeds)
            G.gp['cascades'] = G.new_gp(value_type='python::object',val=cascades)
            G.save(f'{network_path}/{G.gp.ID}.gt')
            net_count+= 1
            print(net_count)
