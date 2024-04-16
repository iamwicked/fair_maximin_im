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
import os
import pickle
from math import ceil, floor

import networkx as nx
import numpy as np
from numpy.random import choice
import secrets


def set_probabilities(G, prob_type='uar'):
    """Set the probabilities of the edges.

    Parameters
    ----------
    G : networkx.DiGraph
        The underlying considered instance.
    prob_type : str
        Identifies how the probabilities are set.

    """
    if prob_type == 'uar':
        probabilities = {e: {'p': secrets.SystemRandom().random()} for e in G.edges()}
    elif prob_type == 'const_1_3':
        probabilities = {e: {'p': 1.0 / 3.0} for e in G.edges()}
    elif prob_type == 'const':
        probabilities = {e: {'p': 0.1} for e in G.edges()}
    elif prob_type == 'const_05':
        probabilities = {e: {'p': 0.05} for e in G.edges()}
    elif prob_type == 'const_1_2':
        probabilities = {e: {'p': 0.5} for e in G.edges()}
    elif prob_type == '0_0.4':
        probabilities = {e: {'p': secrets.SystemRandom().uniform(0,0.4)} for e in G.edges()}
    elif prob_type == '0_0.2':
        probabilities = {e: {'p': secrets.SystemRandom().uniform(0,0.2)} for e in G.edges()}
    elif prob_type == '0_0.1':
        probabilities = {e: {'p': secrets.SystemRandom().uniform(0,0.1)} for e in G.edges()}
    else:
        print('unknown probability value')
        assert 0
    # probabilities = {e : {'p' : rand_zero_half_one()} for e in G.edges()}
    # probabilities = {e : {'p' : rand_eight_fourth()} for e in G.edges()}
    # probabilities = {e : {'p' : 0.0} for e in G.edges()}
    nx.set_edge_attributes(G, probabilities)


def random_direct(U):
    """Randomly direct an undirected graph.

    Parameters
    ----------
    U : networkx.Graph
        Undirected graph to be directed

    Returns
    -------
    networkx.DiGraph
        Directed version of U.

    """
    G = nx.DiGraph()
    if 'partition' in U.graph.keys():
        G.graph['partition'] = U.graph['partition']
    G.add_nodes_from(U.nodes())
    for u, v in U.edges():
        toss = secrets.SystemRandom().random()
        if toss < 0.5:
            G.add_edge(u, v)
        else:
            G.add_edge(v, u)
    return G


def set_singleton_communities(G):
    """Set singleton community structure.

    Parameters
    ----------
    G : networkx.DiGraph
        The underlying considered instance.

    """
    init_communities(G)
    add_communities(G, [(node, [node]) for node in sorted(G.nodes())])


def grow_component(H, max_size):
    """Grow one community by taking BFS components until max_size is reached.

    Parameters
    ----------
    H : networkx.DiGraph
        The underlying considered instance.
    max_size : int
        The maximum community size.

    Returns
    -------
    list
        The grown community as a list of nodes.

    """
    community = []
    bfs = []
    while len(community) < max_size:
        if len(H.nodes()) == 1:
            for v in H.nodes():
                community.append(v)
            break
        if len(H.nodes()) == 0:
            break
        if not bfs:
            source = choice(H.nodes())
            community.append(source)
            bfs = list(nx.bfs_successors(H, source))
            H.remove_node(source)
        else:
            if bfs[0][1]:
                new_node = bfs[0][1][0]
                bfs[0] = (bfs[0][0], bfs[0][1][1:])
                community.append(new_node)
                H.remove_node(new_node)
            else:
                bfs = bfs[1:]
    return community


def init_communities(G):
    """Initialize the data structure within graph for storing the communities.

    Parameters
    ----------
    G : networkx.DiGraph
        The underlying considered instance.

    """
    for u in sorted(G.nodes()):
        G.nodes[u]['communities'] = []
    G.graph['communities'] = {}


def add_communities(G, communities):
    """Add the communiites to the graph.

    Parameters
    ----------
    G : networkx.DiGraph
        The underlying considered instance.
    communities : list
        List of pairs of community id and community nodes.

    """
    for comm_id, comm_nodes in communities:
        G.graph['communities'][comm_id] = comm_nodes
        for u in comm_nodes:
            G.nodes[u]['communities'].append(comm_id)


def set_bfs_communities(G, num_comms):
    """Grow num_comms communities of equal size iteratively by BFS.

    Parameters
    ----------
    G : networkx.DiGraph
        The underlying considered instance.
    num_comms : int
        Number of equal sized communities to be grown (n is suppsed to be
        divisible by num_comms)

    """

    init_communities(G)
    comm_sizes = []
    for i in range(num_comms):
        comm_sizes.append(floor(len(G) / num_comms))
        
    if sum(comm_sizes) != len(G.nodes()):
        for i in range(len(comm_sizes)):
            comm_sizes[i] += 1
            if sum(comm_sizes) == len(G.nodes()):
                break
            
    H = G.copy()
    for i in range(len(comm_sizes)):
        community = grow_component(H, comm_sizes[i])
        add_communities(G, [(i, community)])


def set_random_communities(G, sizes):
    """Construct a random community structure with communities of sizes sizes.

    Parameters
    ----------
    G : networkx.DiGraph
        The underlying considered instance.
    sizes : list
        List of pairs of community ID and size of the community..

    """
    if sum(sizes) != len(G):
        for i in range(len(sizes)):
            sizes[i] -= 1
            if sum(sizes) == len(G):
                break
    init_communities(G)
    nodes = set(G.nodes())
    for comm_id, size in enumerate(sizes):
        comm = choice(list(nodes), size)
        nodes -= set(comm)
        add_communities(G, [(comm_id, comm)])


def set_communities(G, comm_type):
    """Set the community structure according to the type indicated.

    Parameters
    ----------
    G : networkx.DiGraph
        The underlying considered instance.
    comm_type : str
        Indicating type of community structure

    """
    if comm_type == 'singleton':
        set_singleton_communities(G)
    elif comm_type == 'bfs_k':
        num_comms = G.graph['k']
        set_bfs_communities(G, num_comms)
    elif comm_type == 'bfs_k_2':
        num_comms = int(G.graph['k'] / 2)
        set_bfs_communities(G, num_comms)
    elif comm_type == 'bfs_2':
        num_comms = 2
        set_bfs_communities(G, num_comms)
    elif comm_type == 'bfs_5':
        num_comms = 5
        set_bfs_communities(G, num_comms)
    elif comm_type == 'bfs_10':
        num_comms = 10
        set_bfs_communities(G, num_comms)
    elif comm_type == 'bfs_n_10':
        num_comms = int(len(G.nodes()) / 10)
        set_bfs_communities(G, num_comms)  
    elif comm_type == 'bfs_n_2':
        num_comms =  int(len(G.nodes()) / 2)
        set_bfs_communities(G, num_comms)
    elif comm_type == 'rand_equal_k':
        num_comms = G.graph['k']
        set_random_communities(G, [ceil(len(G) / num_comms)
                                   for _ in range(num_comms)])
    elif comm_type == 'rand_imbal_4':
        n = len(G)
        set_random_communities(G, [ceil(4 * n / 10),
                                   ceil(3 * n / 10),
                                   ceil(2 * n / 10),
                                   ceil(1 * n / 10)])
        
    elif comm_type == 'rand_imbal_16':
        n = len(G)
        set_random_communities(G, [ceil(1 * n / 40), ceil(1 * n / 40),
                                   ceil(1 * n / 40), ceil(2 * n / 40),
                                   ceil(1 * n / 40), ceil(1 * n / 40),
                                   ceil(1 * n / 40), ceil(2 * n / 40),
                                   ceil(2 * n / 40), ceil(4 * n / 40),
                                   ceil(4 * n / 40), ceil(5 * n / 40),
                                   ceil(2 * n / 40), ceil(4 * n / 40),
                                   ceil(4 * n / 40), ceil(5 * n / 40)])
                               
    elif comm_type == 'block_stochastic':
        set_block_stochastic_communities(G)
    elif comm_type == 'tsang_all':
        set_communities_tsang(G)
    elif comm_type == 'tsang_gender':
        set_communities_tsang(G, attributes=['gender'])
    elif comm_type == 'tsang_region':
        set_communities_tsang(G, attributes=['region'])
    elif comm_type == 'tsang_region_ethnicity':
        set_communities_tsang(G, attributes=['region', 'ethnicity'])
    elif comm_type == 'tsang_region_gender':
        set_communities_tsang(G, attributes=['region', 'gender'])
    elif comm_type == 'tsang_gender_region_ethnicity':
        set_communities_tsang(G, attributes=['region', 'gender', 'ethnicity'])
    else:
        print("Error: Unknown option.")
        assert(0)
        
        
def set_communities_email_Eu_core(G, U):
    """Set the community structure according to the communities given in the email_Eu_core network.

    Parameters
    ----------
    G : networkx.DiGraph
        The underlying considered instance.

    """
    init_communities(G)
    comm_file = open('./data_set/fish/email-Eu-core-department-labels.txt','r')
    
    node_index = {}
    j = 0;
    for u in sorted(U.nodes()):
        node_index[u] = j
        j += 1

    lines = comm_file.readlines()
  
    for line in lines:
        node , comm = map(str, line.split())
        
        if int(node) in U.nodes():
            if comm in G.graph['communities']:
                if not isinstance(G.graph['communities'][comm], list):
                    G.graph['communities'][comm] = [G.graph['communities'][comm]]
                print(int(node), node_index[int(node)])
                G.graph['communities'][comm].append(node_index[int(node)])
            else:
                G.graph['communities'][comm] = [node_index[int(node)]]
            
            G.nodes[node_index[int(node)]]['communities'].append(comm)
        
    comm_file.close()
    

def directed_barabasi_albert(parameters):
    """Create a directed graph according to the barabasi albert model.

    Parameters
    ----------
    parameters :
        List of parameters:
            n_range: range of graph sizes,
            m: model parameter (number of nodes new nodes get connected to)
            prob_type: type of edge probabilities to be set
            rep_graph: number of graphs to be generated per n.
            seed: random seed to be used

    Returns
    -------
    list
        list of networkx.DiGraphs

    """
    n_range, m, prob_type, rep_graph, seed = parameters
    graphs = []
    for n in n_range:
        for i in range(rep_graph):
            U = nx.barabasi_albert_graph(n, m=m, seed=seed)
            G = U.to_directed()
            G.graph['graphname'] = '_'.join(
                ['ba', 'n', str(n), 'm', str(m), str(i)])
            set_probabilities(G, prob_type=prob_type)
            set_unit_node_weights(G)
            graphs.append(G)
    return graphs


def set_block_stochastic_communities(G):
    """Set communities as encoded by the partition in the generated graph.

    Parameters
    ----------
    G : networkx.DiGraph
        The underlying considered instance.

    """
    init_communities(G)
    add_communities(G, list(enumerate(G.graph['partition'])))


def block_stochastic(parameters):
    """Construct the block stochastic graphs.

    Parameters
    ----------
    parameters :
        List of parameters:
            n_range: range of graph sizes,
            prob_type: type of edge probabilities to be set
            rep_graph: number of graphs to be generated per n.
            p_range: the probability of an edge within communities.
            p_betw_comm: the probability of an edge between communities.
            seed: random seed to be used

    Returns
    -------
    list
        list of networkx.Digraphs

    """
    n_range, prob_type, rep_graph, p_range, p_betw_comm, seed = parameters
    
    graphs = []
    for n in n_range:
        for i in range(rep_graph):
            sizes = [ceil(1 * n / 40), ceil(1 * n / 40), ceil(1 * n / 40), ceil(2 * n / 40),
                     ceil(1 * n / 40), ceil(1 * n / 40), ceil(1 * n / 40), ceil(2 * n / 40),
                     ceil(2 * n / 40), ceil(4 * n / 40), ceil(4 * n / 40), ceil(5 * n / 40),
                     ceil(2 * n / 40), ceil(4 * n / 40), ceil(4 * n / 40), ceil(5 * n / 40)]
            
            for p in p_range: 
                for q in p_betw_comm:
                    probs = q * (np.ones(16) - np.eye(16)) + p * (np.eye(16))
                    print(probs)
                    print(sizes)
                    U = nx.stochastic_block_model(sizes, probs, seed=seed)
                    G = random_direct(U)
                    G.graph['graphname'] = '_'.join(
                            ['block_stochastic', 'n', str(n), 'q', str(format(q, '.2f')), str(i)])
                    set_probabilities(G, prob_type=prob_type)
                    set_unit_node_weights(G)
                    graphs.append(G)
    return graphs



def youtube(parameters):
    """Read the youtube network from ./data_set/youtube/.

    Parameters
    ----------
    parameters : pair
            prob_type: type of edge probabilities to be set
            rep_graph: number of graphs to be generated per n

    Returns
    -------
    networkx.DiGraph

    """
    (prob_type, rep_graph) = parameters
    graphs = []
    folder = './data_set/youtube/'
    
    graph_file = "com-youtube.ungraph.txt"
    with open(folder + graph_file) as f:
        edge_list = [tuple(map(int, line.split()[:2]))
                      for line in f]
                
    U = nx.Graph(edge_list)
    G = U.to_directed()

    comm_file = open('./data_set/youtube/com-youtube.all.cmty.txt','r')
    
    lines = comm_file.readlines()
    
    nodes_in_comm = []
    for line in lines:
        nodes = map(str, line.split())
        for v in nodes:
            if int(v) not in nodes_in_comm:
                nodes_in_comm.append(int(v))
        
    comm_file.close()

    G.remove_nodes_from([v for v in G.nodes() if v not in nodes_in_comm])    
    rand_node = choice(G.nodes())
    nodes_bfs = list(nx.bfs_tree(G, rand_node).nodes())
    
    new_nodes = []
    for i in range(3000):                
        new_nodes.append(nodes_bfs[i])
    G.remove_nodes_from([v for v in G.nodes() if v not in new_nodes])

    node_index = {}
    j = 0;
    for u in sorted(G.nodes()):
        node_index[u] = j
        j += 1

    M = G.copy()
    G = nx.relabel_nodes(G, node_index)
    init_communities(G)
    
    new_com_nodes = {}   
    new_comm_num = 0

    comm = 0                           
    for line in lines:
        nodes = map(str, line.split())
        nodess = []
        for v in nodes:
            if int(v) in M.nodes():
                nodess.append(node_index[int(v)])
        if len(nodess)>0:
            new_com_nodes[new_comm_num] = nodess
            new_comm_num += 1
            
    comm_singletn = []
    for comms in new_com_nodes.keys():
        if len(new_com_nodes[comms]) <= 2:
            comm_singletn.append(comms)
    
    for comm in comm_singletn:
        new_com_nodes.pop(comm)

    seen = set()
    for key in new_com_nodes.copy():
        value = tuple(new_com_nodes[key])
        if value in seen:
            del new_com_nodes[key]
        else:
            seen.add(value) 
            
    c_num = 0
    for c in new_com_nodes.keys():
        G.graph['communities'][c_num] = new_com_nodes[c]
        c_num += 1

    G.graph['graphname'] = graph_file[:-4]
    set_probabilities(G, prob_type=prob_type)
    set_unit_node_weights(G)
    graphs.append(G)
    
    return graphs


def graph_fish(parameters):
    """Read the networks of Fish et al's study from ./data_set/fish/.

    Parameters
    ----------
    parameters : tuple
        List of parameters:
            prob_type: type of edge probabilities to be set
            rep_graph: number of graphs to be generated per n
            file_list: the name of the file to be read
    Returns
    -------
    list
        list of networkx.DiGraphs

    """
    (prob_type, rep_graph, file_list) = parameters
    graphs = []
    folder = './data_set/fish/'
    directed_graphs = ["irvine.txt",
                       "email-EuAll.txt",
                       "email-Eu-core.txt",
                       "arena.txt",
                       "ca-GrQc.txt",
                       "ca-HepTh.txt"]
    undirected_graphs = ["facebook_combined.txt"]

    print(file_list)
    for graph_file in file_list:
        print('\n', graph_file)
        
        with open(folder + graph_file) as f:
                edge_list = [tuple(map(int, line.split()[:2])) for line in f]
                
        if graph_file in directed_graphs:
            G = nx.DiGraph(edge_list)    
        elif graph_file in undirected_graphs:
            U = nx.Graph(edge_list)
            G = U.to_directed()
        else:
            print("Error: Unknown option.")
            assert(0)
            
        largest_cc = max(nx.weakly_connected_components(G), key=len)
        G.remove_nodes_from([v for v in G.nodes() if v not in largest_cc])
        
        G.graph['graphname'] = graph_file[:-4]
        set_probabilities(G, prob_type=prob_type)
        set_unit_node_weights(G)
        graphs.append(G)
        print(G.graph)
        print('m', len(G.edges()))
        print('n', len(G.nodes()))

    return graphs


def set_unit_node_weights(G):
    """Set all node's weights to 1, weights are used for weighted IM problems.

    Parameters
    ----------
    G : networkx.DiGraph
        The underlying considered instance.

    """
    for u in G.nodes():
        G.nodes[u]['weight'] = 1


def graph_tsang(parameters):
    """Read the networks of Tsang et al's study from ./tsang/networks/.

    Parameters
    ----------
    parameters : pair
        List of parameters:
            prob_type: type of edge probabilities to be set
            rep_graph: number of graphs to be generated per n
    Returns
    -------
    list
        list of networkx.DiGraphs

    """
    (prob_type, rep_graph) = parameters
    graphs = []
    folder = './tsang/networks/'
    graphnames = ['spa_500_{}'.format(graphidx) for graphidx in range(20)]
    for graphname in graphnames[:rep_graph]:
        G = pickle.load(
            open(folder + 'graph_{}.pickle'.format(graphname), 'rb'))

        set_probabilities(G, prob_type=prob_type)
        set_unit_node_weights(G)
        G.graph['graphname'] = 'tsang_' + graphname
        graphs.append(G)
    return graphs


def set_communities_tsang(
        G,
        attributes=['region', 'ethnicity', 'age', 'gender', 'status']):
    """Set communities as stored in the graph.

    Parameters
    ----------
    G : networkx.DiGraph
        The underlying considered instance.
    attributes : list
        list of attributes inducing the desired community structure

    """
    G.graph['communities'] = {}

    for u in sorted(G.nodes()):
        a = []
        for attribute in attributes:
            values = np.unique([G.nodes[v][attribute] for v in sorted(G.nodes())])
            for val in values:
                if G.nodes[u][attribute] == val:
                    a.append(val)
        G.nodes[u]['communities'] = a

    for attribute in attributes:
        values = np.unique([G.nodes[v][attribute] for v in sorted(G.nodes())])
        for val in values:
            G.graph['communities'][val] = [
                v for v in G.nodes() if G.nodes[v][attribute] == val]
        
