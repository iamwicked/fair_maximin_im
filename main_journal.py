##############################################################################
# Copyright (c) 2022, Ruben Becker, Sajjad Ghobadi
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#     * Redistributions of source code must retain the above copyright
#       notice, this list of conditions and the following disclaimer.
#     * Redistributions in binary form must reproduce the above copyright
#       notice, this list of conditions and the following disclaimer in the
#       documentation and/or other materials provided with the distribution.
#     * The name of the author may not be used to endorse or promote products
#       derived from this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE AUTHOR ``AS IS'' AND ANY EXPRESS OR IMPLIED
# WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
# MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO
# EVENT SHALL THE AUTHOR BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
# SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
# PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
# OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
# WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
# OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
# ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
##############################################################################

"""Main module that executes the experiments."""

import os
import sys
import threading
import time

import networkx as nx
import numpy as np
from numpy.random import choice, seed

import influence_max as im
import maximin_fish as mf
import cpp_proxy as cpp
import print_functions as pf
import tsang_maximin as tm
from tsang.algorithms import rounding
import generation as gen


def sample_sets(G, vector, times, type):
    """Sample times many sets from vector depending on type.

    Parameters
    ----------
    G : networkx.DiGraph
        The underlying considered instance.
    vector : depends on type
        If type is deterministic, vector is a set.
        If type is node_based or swap_rounding, vector is a dictionary with
            keys being the nodes and values the probability of the respective
            node.
        If type is set_based, vector is a dictionary with keys being ints
            representing the sets (binary) and values are the probability
            of the respective set.
    times : int
        How many times to sample a set.
    type : str
        The type that specifies how to sample the sets, this depends on whether
        this is a deterministic, node_based, or set_based problem.

    Returns
    -------
    list of lists
        The list of sampled sets.

    """
    if type == 'deterministic':
        sets = [vector for _ in range(times)]
    elif type == 'node_based':
        sets = [[v for v in vector.keys() if choice(
            [0, 1], p=[1 - vector[v], vector[v]])]
            for _ in range(times)]
    elif type == 'swap_rounding':
        x_items_list = sorted(vector.items())
        x = np.array([x_items[1] for x_items in x_items_list])
        rounded_xs = [rounding(x) for _ in range(times)]
        sets = [[v for v in G.nodes() if rounded_xs[i][v]]
                for i in range(times)]
    elif type == 'set_based':
        sets = [im.number_to_set(G, choice(list(vector.keys()),
                                           p=list(vector.values())))
                for _ in range(times)]
    else:
        print("Error: Unknown option.", type)
        assert(False)
    return sets


def comp_ex_post(G, solution, fct_name):
    """Compute ex_post values by sampling one set.

    Parameters
    ----------
    G : networkx.DiGraph
        The underlying considered instance.
    solution : depends on fct_name
        The computed solution.
    fct_name : str
        The name of the function that computed solution.

    Returns
    -------
    
    The ex_post fairness value and spread obtained for one sampled set.

    """
    sets = sample_sets(G, solution, 1,
                       ex_post_sampling_types[fct_name])
    
    sigma_sets_ep = [im.sigma(G, S, weight='weight')
                             for S in sets]
    exp_sigma = np.mean(sigma_sets_ep)
    min_probs_sets = [min(im.sigma_C(G, S).values())
                      for S in sets]
    ex_post_apx_val = np.mean(min_probs_sets)
    
    return ex_post_apx_val, exp_sigma


def comp_ex_ante(G, solution, fct_name):
    """Compute ex_ante values ((0.1, 0.1)-approximately) and expected spread.

    Parameters
    ----------
    G : networkx.DiGraph
        The underlying considered instance.
    solution : depends on fct_name
        The computed solution.
    fct_name : str
        The name of the function that computed solution.

    Returns
    -------
    
    The ex_ante value obtained by (0.1, 0.1)-approximating in case of
    the node_based problem and exact computation in case of the set_based
    problem. 
    The expected influence according to the fct_name and returned solution.

    """
    if fct_name in deterministic_algorithms + ['tim']:
        comm_probs = im.sigma_C(G, solution)
        return min(comm_probs.values()), im.sigma(G, solution, weight='weight')
    else:
        assert(fct_name in probabilistic_algorithms + cpp_algos)
        sets = sample_sets(G, solution, 100,
                                ex_post_sampling_types[fct_name])
        sigma_sets = [im.sigma(G, S, weight='weight')
                                      for S in sets]
        sigma = np.mean(sigma_sets)
        if fct_name == 'set_based':
            comm_probs = im.sigma_C_p(G, solution)
        elif fct_name in ['node_based', 'tsang_maximin', 'uniform_node_based']:
            comm_probs = im.sigma_C_x(G, solution, 0.1, 0.1)
        else:
            print("Error: Unknown option:", fct_name)
            assert(False)
        return min(comm_probs.values()), sigma


def read_graph(graph_file):
    """Read graph.

    Parameters
    ----------
    graph_file : The path of the folder containing graph files

    Returns
    -------
    G : networkx.DiGraph
        The underlying considered instance.
    """
    
    G = nx.read_edgelist(graph_file + '/graph_ic.txt', nodetype = int,
                         data=(('p',float),), create_using = nx.DiGraph())
      
    gen.init_communities(G)
    G.graph['nodes_in_comms'] = set()    
    
    comm_file = open(graph_file + '/community.txt','r')
    lines = comm_file.readlines()[1:]
    for line in lines:
        comm, *nodes = map(str, line.split())
        G.graph['communities'][comm] = [int(v) for v in nodes]
        
        for u in nodes:
            G.nodes[int(u)]['communities'].append(comm)
            G.graph['nodes_in_comms'].add(int(u))        
    
    comm_file.close()
    gen.set_unit_node_weights(G)
    return G


def execute(function, graph_file, k, rep_per_graph_i, num):
    """Execute function on G (using graph_file) and writes results to out_file (global).

    Parameters
    ----------
    function : function that takes a networkx graph and returns a solution.
    graph_file : The path of the folder containing graph files.
    k : int
        Budget (size of the seed set)
    rep_per_graph_i : int
        The number showing repetition per graph.
    num : int
        Number specifying the name of output file (used for cpp implementations).

    """
    fct_name = function.__name__
        
    G = read_graph(graph_file)
    G.graph['T'] = no_le_graphs
    G.graph['k'] = k
    
    start = time.time()
    G.graph['leg'], G.graph['leg_tsang'] = im.gen_leg(G)  
    
    if fct_name in python_algos:        
        solution = function(G)
        ex_time = time.time() - start
        
    elif fct_name in cpp_algos: 
        ex1_time = time.time() - start
        solution, ex2_time = function(G, graph_file, k, num)
        ex_time = ex1_time + ex2_time
    
    elif fct_name == 'tim':
        ex1_time = time.time() - start
        solution, ex2_time = function(graph_file, k, num)
        ex_time = ex1_time + ex2_time
    else:
        print("Error: Unknown option.")
        assert(0)
          
    G.graph['leg'], G.graph['leg_tsang'] = im.gen_leg(G) 

    n, m = len(G), len(G.edges())
    ex_ante_val, sigma = comp_ex_ante(G, solution, fct_name)  
    ex_post_val, exp_sigma = comp_ex_post(G, solution, fct_name)
                                   
    x = graph_file.split("/")
    graph_file = x[3]
    
    res = (graph_file,
           fct_name,
           n,
           m,
           k,
           ex_time,
           ex_ante_val,
           ex_post_val,
           rep_per_graph_i,
           sigma,
           exp_sigma)

    with open(out_file, 'a') as f:
        f.write('\t\t\t'.join('%s' % x for x in res) + '\n')

    print("\n")
    print("function: ".ljust(30, ' '), fct_name)
    print("maximizing solution: ".ljust(30, ' '))
    # if fct_name == 'set_based':
    #     pf.myprint([(im.number_to_set(G, s), solution[s])
    #                 for s in solution.keys()])
    # else:
    #     pf.myprint(solution)
    #print("ex_ante: ".ljust(30, ' '), ex_ante_val)
    print("ex_ante: ".ljust(30, ' '), ex_ante_val)
    print("ex_post: ".ljust(30, ' '), ex_post_val)
    print("sigma: ".ljust(30, ' '), sigma)
    print("running time: ".ljust(30, ' '), ex_time)
    if fct_name not in cpp_algos:
        print("number of live edge graphs: ".ljust(30, ' '), G.graph['T'])
    print("\n")


def generate_executables(
        functions,
        instance_folder,
        k_range,
        rep_per_graph,
        no_le_graphs):
    """Generate executables by generating graphs and combining them.

    Parameters
    ----------
    functions : list
        List of all the functions that are supposed to be executed.
    instance_folder : the folder of the instance to be executed.
    k_range : list
        Range of k-values to be tested.
    rep_per_graph : int
        Number of times the execution should be repeated per graph.
    no_le_graphs : int
        Number of live-edge graphs to be generated and used in all following
        computations.

    Returns
    -------
    list
        List of executables, i.e., (function, graph)-pairs.

    """
    num = 0     
    executables = []
    instance_names = sorted([f for f in os.listdir('./data_set/' + instance_folder)
                            if not f.startswith(".")])
    for graph_folder_name in instance_names:
        for rep_per_graph_i in range(rep_per_graph):
            for k in k_range:
                for function in functions:
                    num += 1
                    executables.append((function,
                                        './data_set/' + instance_folder + '/' + graph_folder_name ,
                                         k, rep_per_graph_i, num)) 
            print('collected', len(executables), 'executables')
    return executables


#############
# main
#############

# forbid python 2 usage
version = sys.version_info[0]
if version == 2:
    sys.exit("This script shouldn't be run by python 2 ")

# do not set seed specifically
s = None
seed(s)

# the following lists are used to specify the way in which sets are sampled
# for the respective algorithms

deterministic_algorithms = [
    'myopic_fish',
    'naive_myopic_fish',
    'greedy_maximin']

probabilistic_algorithms = [
    'uniform_node_based',
    'tsang_maximin']

python_algos = deterministic_algorithms + probabilistic_algorithms
cpp_algos = [
    'node_based',
    'set_based']


sampling_types = [
    'deterministic',
    'node_based',
    'set_based',
    'swap_rounding'
]

ex_ante_sampling_types = {
    'set_based': 'set_based',
    'node_based': 'node_based',
    'uniform_node_based': 'node_based',
    'tsang_maximin': 'node_based'}
for alg in deterministic_algorithms + ['tim']:
    ex_ante_sampling_types[alg] = 'deterministic'

ex_post_sampling_types = ex_ante_sampling_types
ex_post_sampling_types['tsang_maximin'] = 'swap_rounding'

print('++++++++++++++++++++++++++++++++++++++++++++++++++++')
print('++++++ Expecting experiment_type as argument. ++++++')
print('++++++++++++++++++++++++++++++++++++++++++++++++++++')

# read number of desired processes from the shell
experiment_type = sys.argv[1]
if len(sys.argv) == 3:
    number_of_processes = int(sys.argv[2])
else:
    number_of_processes = 1

# default values for experiments     
functions = [cpp.tim,
             mf.myopic_fish,
             cpp.set_based,
             im.uniform_node_based,                
             cpp.node_based,
             tm.tsang_maximin,
             im.greedy_maximin,
             mf.naive_myopic_fish
             ]

no_le_graphs = 100
rep_per_graph = 10                                              

# specify experiment dependent parameters
if experiment_type == 'ba-singletons-0_0.4-incr_n':
    instance_folder = 'ba-singletons-0_0.4-incr_n'
    k_range = [20]

elif experiment_type == 'ba-bfs_cover-k-0_0.4-incr_n':
    instance_folder = 'ba-bfs_cover-k-0_0.4-incr_n'
    k_range = [20]

elif experiment_type == 'ba-rand_imbal-4-0_0.4-incr_n':
    instance_folder = 'ba-rand_imbal-4-0_0.4-incr_n'
    k_range = [20]    
    
elif experiment_type == 'block_stochastic-const_05-incr_p_p_0.1_1_q_0.1':
    instance_folder = 'block_stochastic-const_05-incr_p_p_0.1_1_q_0.1'
    k_range = [8]

elif experiment_type == 'block_stochastic-const_05-incr_q_p_0.1_q_0_p':
    instance_folder = 'block_stochastic-const_05-incr_q_p_0.1_q_0_p'
    k_range = [8] 
            
elif experiment_type == 'tsang-gender-0_0.4':
    instance_folder = 'tsang-gender-0_0.4'
    k_range = [5 * i for i in range(1, 11)]
    
elif experiment_type == 'tsang-gender-region-ethnicity-0-0.4':
    instance_folder = 'tsang-gender-region-ethnicity-0-0.4'
    k_range = [5 * i for i in range(1, 11)] 
            
elif experiment_type == 'arena_0_0.2_bfs_comm_10':           
    instance_folder = 'arena_0_0.2_bfs_comm_10'
    k_range = [5, 10 , 20, 50, 100] 

elif experiment_type == 'arena_0_0.2_bfs_comm_n_10':           
    instance_folder = 'arena_0_0.2_bfs_comm_n_10'
    k_range = [5, 10 , 20, 50, 100] 

elif experiment_type == 'arena_0_0.2_comm_n_2':           
    instance_folder = 'arena_0_0.2_comm_n_2'
    k_range = [5, 10 , 20, 50, 100] 
        
elif experiment_type == 'arena-0_0.2-imbal_16':           
    instance_folder = 'arena-0_0.2-imbal_16'
    k_range = [5, 10 , 20, 50, 100] 
               
elif experiment_type == 'irvine_0_0.2_bfs_comm_10':           
    instance_folder = 'irvine_0_0.2_bfs_comm_10'
    k_range = [5, 10 , 20, 50, 100]

elif experiment_type == 'irvine_0_0.2_bfs_comm_n_10':           
    instance_folder = 'irvine_0_0.2_bfs_comm_n_10'
    k_range = [5, 10 , 20, 50, 100]

elif experiment_type == 'irvine_0_0.2_bfs_comm_n_2':           
    instance_folder = 'irvine_0_0.2_bfs_comm_n_2'
    k_range = [5, 10 , 20, 50, 100]
        
elif experiment_type == 'irvine-0_0.2-imbal_16':           
    instance_folder = 'irvine-0_0.2-imbal_16'
    k_range = [5, 10 , 20, 50, 100]
            
elif experiment_type == 'ca-GrQc_0_0.2_bfs_comm_10':           
    instance_folder = 'ca-GrQc_0_0.2_bfs_comm_10'
    k_range = [5, 10 , 20, 50, 100]

elif experiment_type == 'ca-GrQc_0_0.2_bfs_comm_n_10':           
    instance_folder = 'ca-GrQc_0_0.2_bfs_comm_n_10'
    k_range = [5, 10 , 20, 50, 100]

elif experiment_type == 'ca-GrQc_0_0.2_bfs_comm_n_2':           
    instance_folder = 'ca-GrQc_0_0.2_bfs_comm_n_2'
    k_range = [5, 10 , 20, 50, 100]
    
elif experiment_type == 'ca-GrQc-0_0.2-imbal_16':           
    instance_folder = 'ca-GrQc-0_0.2-imbal_16'
    k_range = [5, 10 , 20, 50, 100]
            
elif experiment_type == 'ca-HepTh_0_0.2_bfs_comm_10':           
    instance_folder = 'ca-HepTh_0_0.2_bfs_comm_10'
    k_range = [5, 10 , 20, 50, 100]

elif experiment_type == 'ca-HepTh_0_0.2_bfs_comm_n_10':           
    instance_folder = 'ca-HepTh_0_0.2_bfs_comm_n_10'
    k_range = [5, 10 , 20, 50, 100]

elif experiment_type == 'ca-HepTh_0_0.2_bfs_comm_n_2':           
    instance_folder = 'ca-HepTh_0_0.2_bfs_comm_n_2'
    k_range = [5, 10 , 20, 50, 100]
        
elif experiment_type == 'ca-HepTh-0_0.2-imbal_16':           
    instance_folder = 'ca-HepTh-0_0.2-imbal_16'
    k_range = [5, 10 , 20, 50, 100]
        
elif experiment_type == 'facebook_combined_0_0.2_bfs_comm_10':           
    instance_folder = 'facebook_combined_0_0.2_bfs_comm_10'
    k_range = [5, 10 , 20, 50, 100]

elif experiment_type == 'facebook_combined_0_0.2_bfs_comm_n_10':           
    instance_folder = 'facebook_combined_0_0.2_bfs_comm_n_10'
    k_range = [5, 10 , 20, 50, 100]

elif experiment_type == 'facebook_combined_0_0.2_bfs_comm_n_2':           
    instance_folder = 'facebook_combined_0_0.2_bfs_comm_n_2'
    k_range = [5, 10 , 20, 50, 100]
    
elif experiment_type == 'facebook_combined-0_0.2-imbal_16':           
    instance_folder = 'facebook_combined-0_0.2-imbal_16'
    k_range = [5, 10 , 20, 50, 100]

elif experiment_type == 'email-Eu-core_0_0.2':           
    instance_folder = 'email-Eu-core_0_0.2'
    k_range = [5, 10 , 20, 50, 100] 
    
elif experiment_type == 'email-Eu-core_0_0.2_bfs_comm_2':           
    instance_folder = 'email-Eu-core_0_0.2_bfs_comm_2'
    k_range = [5, 10 , 20, 50, 100]
    
elif experiment_type == 'email-Eu-core_0_0.2_bfs_comm_10':           
    instance_folder = 'email-Eu-core_0_0.2_bfs_comm_10'
    k_range = [5, 10 , 20, 50, 100]

elif experiment_type == 'email-Eu-core_0_0.2_bfs_comm_n_10':           
    instance_folder = 'email-Eu-core_0_0.2_bfs_comm_n_10'
    k_range = [5, 10 , 20, 50, 100]
    
elif experiment_type == 'email-Eu-core_0_0.2_bfs_comm_n_2':           
    instance_folder = 'email-Eu-core_0_0.2_bfs_comm_n_2'
    k_range = [5, 10 , 20, 50, 100]

elif experiment_type == 'email-Eu--core_0_0.2_imbal_16':           
    instance_folder = 'email-Eu--core_0_0.2_imbal_16'
    k_range = [5, 10 , 20, 50, 100]
    
elif experiment_type == 'youtube_0_0.1_n_3000_no_singl':           
    instance_folder = 'youtube_0_0.1_n_3000_no_singl'
    k_range = [5, 10, 20, 50, 100]
   
else:
    print("Error: Unknown option.")
    assert(0)

# create output file with header
folder = './results_journal/'
out_file = folder + experiment_type + '.txt'
if os.path.exists(out_file):
    out_file = out_file[:-4] + '_' + str(int(time.time())) + '.txt'
print('Output is written to:', out_file, '\n')
header = [
    'graphname',
    'algorithm',
    'n',
    'm',
    'k',
    'running_time',
    'ex_ante',
    'ex_post',
    'rep_per_graph',
    'spread',
    'exp_spread'
    ]
with open(out_file, 'a') as f:
    f.write('\t\t\t'.join('%s' % x for x in header) + '\n')


#generate the various experiments
executables = generate_executables(
    functions,
    instance_folder,
    k_range,
    rep_per_graph,
    no_le_graphs)

# run experiments in parallel (if number_of_processes > 1)
thread_list = []
for executable in executables:
    thread_list.insert(
        0,
        threading.Thread(
            target=execute,
            args=(executable[0], executable[1], executable[2], executable[3], executable[4])))

while thread_list:
    if threading.active_count() < number_of_processes + 1:
        thread_list.pop().start()