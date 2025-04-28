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

"""Module proxy for the cpp implementation."""

import os
import uuid
import influence_max as im
import operator
from decimal import Decimal
import secrets
    

def tim(graph_file, k, num):
    """Call the tim implementation for influence maximization.

    Parameters
    ----------
    graph_file : string
        The path of the folder containing graph files.
    k : int
        budget
    num : int
        Number specifying the name of output file using cpp implementation.

    Returns
    -------
    S : list of nodes
       computed solution
    running_time : float
       The running time of tim to return the computed solution.
    
    """    
    log_file = str(secrets.SystemRandom().random()) + str(uuid.uuid1()) + str(num) + '.txt'
    path_make = 'cd' + ' ' + './TIM' + '&& make'
    process_make = os.popen(path_make)
    output_make = process_make.read()
    process_make.close()
    print(output_make)

    tim_command = 'cd' + ' ' + './TIM' + '&& ./multiplicative_weight -model IC -dataset' + ' ' +'.'+ graph_file + ' ' + '-epsilon 0.1 -k' + ' ' + str(k) + ' ' + '-setting' + ' ' + 'tim' + ' ' + '-file_name' + ' ' + str(log_file)
    process = os.popen(tim_command)
    output = process.read()
    process.close()
    print(output)
    file = './TIM/' + log_file
     
    S = set()
    with open(file) as f:
        line1 = f.readline()
        running_time = float(line1)
        
        line2 = f.readline()
        nodes = map(str, line2.split())
        for v in nodes:
            S.add(int(v))
    os.remove(file)
    
    return S, running_time


def set_based(G, graph_file, k, num):
    """Call the cpp multiplicative weight routine for the set case.

    Parameters
    ----------
    G : networkx.DiGraph
        The underlying considered instance.
    graph_file : string
        The path of the folder containing graph files.
    k : int
        budget
    num : int
        Number specifying the name of output file using cpp implementation.

    Returns
    -------
    p : dict
       dictionary with keys being ints representing the sets (binary) and values are the probability of the respective set.
    running_time : float
       The running time of the multiplicative weight routine to return the computed solution.
    
    """
    path_make = 'cd' + ' ' + './TIM' + '&& make'
    process_make = os.popen(path_make)
    output_make = process_make.read()
    process_make.close()
    print(output_make)
    
    log_file = str(secrets.SystemRandom().random()) + str(uuid.uuid1()) + str(num) + '.txt'
    tim_command = 'cd' + ' ' + './TIM' + '&& ./multiplicative_weight -model IC -dataset' + ' ' +'.' + graph_file + ' ' + '-epsilon 0.1 -k' + ' ' + str(k) + ' ' + '-setting' + ' '+ 'set' + ' ' + '-file_name' + ' ' + str(log_file)
    file = './TIM/'+ log_file
    process = os.popen(tim_command)
    output = process.read()
    process.close()
    print(output)   
        
    p = {}
    with open(file) as f:
        line1 = f.readline()
        running_time = float(line1)
            
    with open(file) as f:
        lines = f.readlines()[1:]
        for line in lines:
            *sets , pro = map(float, line.split())
            if len(sets) == 0.0:
                p[0.0] = pro
            else:
                p[im.set_to_number(G, sets)] = pro
    os.remove(file)
    sum_p = sum(p.values())
    if sum_p != 1:
        p = dict(sorted(p.items(), key=lambda x: x[1]))
        if sum_p > 1:
            x = sum_p - 1
            for S_pro in p:
                if p[S_pro] >= x:
                      p[S_pro] -= x
                      break
        if sum_p < 1:
            x = 1 - sum_p
            max_S = max(p.items(), key=operator.itemgetter(1))[0]
            p[max_S] += x; 
             
    return p, running_time


def node_based(G, graph_file, k, num):
    """Call the cpp multiplicative weight routine for the node case.

    Parameters
    ----------
    G : networkx.DiGraph
        The underlying considered instance.
    graph_file : string
        The path of the folder containing graph files.
    k : int
        budget
    num : int
        Number specifying the name of output file using cpp implementation.

    Returns
    -------
    p : dict
       dictionary with keys being nodes and values the probabilities.
    running_time : float
       The running time of the multiplicative weight routine to return the computed solution.


    """
    
    log_file = str(secrets.SystemRandom().random()) + str(uuid.uuid1()) + str(num) + '.txt'

    path_make = 'cd' + ' ' + './TIM' + '&& make'
    process_make = os.popen(path_make)
    output_make = process_make.read()
    process_make.close()
    print(output_make)
    
    tim_command = 'cd' + ' ' + './TIM' + '&& ./multiplicative_weight -model IC -dataset' + ' ' +'.'+ graph_file + ' ' + '-epsilon 0.1 -k' + ' ' + str(k) + ' ' + '-setting' + ' '+ 'node' + ' ' + '-file_name' + ' ' + str(log_file)
    process = os.popen(tim_command)
    output = process.read()
    process.close()
    print(output)
    
    x = {}
    file = './TIM/'+ log_file
    with open(file) as f:
        first_line = f.readlines()[:1]
        for l in first_line:
            running_time = float(l)
        
    with open(file) as f:
        lines = f.readlines()[1:]
        for line in lines:
            node , pro = map(float, line.split())
            x[node] = pro
    os.remove(file)

    return x, running_time 
   
