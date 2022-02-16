This is the source code used for the experimental part of the paper 
    "Fairness in Influence Maximization through Randomization" by Ruben Becker, Gianlorenzo D'Angelo, Sajjad Ghobadi, and Hugo Gilbert, published at JAIR. A complete version is available on the arxiv under: https://128.84.4.12/pdf/2010.03438.pdf

For execution of the experiments use
    python main_journal.py experiment_type N

where experiment_type is one of the 'experiment_type' used in the paper. For instance:
	experiment_type == 'ba-singletons-0_0.4-incr_n':
    	instance_folder = 'ba-singletons-0_0.4-incr_n'
    	k_range = [20]

and N is the number of experiments that are supposed to be run in parallel.

For creating the instances use
    python main_generate.py experiment_type


- For creating the networks used by Fish et al. [1] download the data sets into a subfolder 'fish' and put it into the folder 'data_set'.  Also download the 'youtube' network into a subfolder 'youtube' and put it into the same folder ('data_set').

- For running tests on the networks used by Tsang et al. [2], download the instances into the subfolder 'networks' and put it into the folder 'tsang'.

- The files have the following content:
    - main_journal.py: contains the infrastructure for the execution of the experiments.
    - generation.py: contains the functions that generate the different instances.
    - influence_max.py: contains the necessary implementations of influence maximization functions: computation of the influence sigma of a set (and related functions sigma_v, sigma_C, see the paper), computation of nodes reachable from a given set, generation of live edge graphs, greedy algorithms for influence maximization and maximin
    - maximin_fish.py: contains the implementations of the myopic and naive myopic routines of Fish et al.
    - main_generate.py: contains the infrastructure for generating the instances used in the paper.
    - print_functions.py: contains some functions for prettier printing of graphs etc.
    - tsang_maximin.py: contains the function used to call the code of Tsang et al. [2]
    - cpp_proxy.py: contains the functions for executing the c++ implementations, the multiplicative weight routine of Young [3] for the set-based and node-based problems and TIM implementation for influence maximization [4].

- Folder 'TIM' contains the c++ implementation of the proposed algorithms in the paper. It uses TIM implementation for influence maximization.
- The execution generates an output file within the folder 'results_journal' with the name being the experiment_type.


- The experiment_types have the following meaning, see the paper for the interpretation of the community structure types:
  Here we only explain some of the experiment_types.

  - 'ba-singletons-0_0.4-incr_n':         		+ Barabasi albert graph (parameter m=2),
                                        		+ singleton community structure
                                       			+ edge weights uniformly at random in [0,0.4]
                                      			+ k = 20
                                       			+ n increasing from 40 to 200 in steps of 20

  - 'ba-bfs_cover-k-0_0.4-incr_n':    			+ Barabasi albert graph (parameter m=2),
                                         		+ BFS community structure with 20 communities
                                       		   	+ edge weights uniformly at random in [0,0.4]
                                       			+ k = 20
                                       			+ n increasing from 40 to 200 in steps of 20

  - 'ba-rand_imbal-4-0_0.4-incr_n':      		+ Barabasi albert graph (parameter m=2),
                                           		+ random imbalanced community structure with 4 communities
                                         		+ edge weights uniformly at random in [0,0.4]
                                       			+ k = 20
                                       			+ n increasing from 40 to 200 in steps of 20

  - 'block_stochastic-const_05-incr_p_p_0.1_1_q_0.1':	+ block stochastic graph,
                                           	   	+ block community structure with q = 0.1 and p increasing 
							  from q to 1 in steps of 0.1
   			                             	+ edge weights 0.05
                                      		  	+ k = 8
                                       			+ n = 200

  - 'block_stochastic-const_05-incr_q_p_0.1_q_0_p':	+ block stochastic graph,
                                           	   	+ block community structure with p = 0.1 and q increasing 
							  from 0 to p in steps of 0.01
   			                              	+ edge weights 0.05
                                       		  	+ k = 8
                                       			+ n = 200

  - 'tsang-gender-0_0.4':              			+ instances of Tsang et al. (2019)
                                            		+ community structure induced by attribute gender
                                            		+ edge weights uniformly at random in [0,0.4]
                                            		+ k increasing from 5 to 50 in steps of 5
                                       			+ n = 500

  - tsang-gender-region-ethnicity-0-0.4:  		+ instances of Tsang et al. (2019)
                                            		+ community structure induced by attributes gender, region 
							  and ethnicity
                                            		+ edge weights uniformly at random in [0,0.4]
                                       			+ k increasing from 5 to 50 in steps of 5
                                       			+ n = 500

  - 'arena_0_0.2_bfs_comm_10':            		+ instances used by Fish et al. (2019)
                                            		+ BFS community structure with 10 communities
                                            		+ edge weights uniformly at random in [0,0.2]
                                       			+ k = [5, 10 , 20, 50, 100] 
                                       			+ n = 1133

  - 'arena_0_0.2_bfs_comm_n_10':          		+ instances used by Fish et al. (2019)
                                            		+ BFS community structure with n/10 communities
                                            		+ edge weights uniformly at random in [0,0.2]
                                       			+ k = [5, 10 , 20, 50, 100] 
                                       			+ n = 1133

  - 'arena_0_0.2_comm_n_2':            			+ instances used by Fish et al. (2019)
                                            		+ BFS community structure with n/2 communities
                                            		+ edge weights uniformly at random in [0,0.2]
                                            		+ k = [5, 10 , 20, 50, 100] 
                                       			+ n = 1133

  - 'arena-0_0.2-imbal_16':               		+ instances used by Fish et al. (2019)
                                            		+ imbalanced community structure with 16 communities
                                           		+ edge weights uniformly at random in [0,0.2]
                                       			+ k = [5, 10 , 20, 50, 100] 
                                       			+ n = 1133

  - 'email-Eu-core_0_0.2':               		+ instances used by Fish et al. (2019)
                                            		+ community structure induced by departments
                                            		+ edge weights uniformly at random in [0,0.2]
                                       			+ k = [5, 10 , 20, 50, 100] 
                                       			+ n = 1005



[1] Fish, Bashardoust, Boyd, Friedler, Scheidegger, Venkatasubramanian. Gaps in Information Access in Social Networks. WWW 2019.
    http://sorelle.friedler.net/papers/access_networks_www19.pdf
[2] Tsang, Wilder, Rice, Tambe, Zick. Group-Fairness in Influence Maximization. IJCAI 2019.
    https://www.ijcai.org/Proceedings/2019/0831.pdf
[3] Young, Neal E. "Randomized Rounding Without Solving the Linear Program." SODA. Vol. 95. 1995.
    https://arxiv.org/pdf/cs/0205036.pdf
[4] Tang, Youze, Xiaokui Xiao, and Yanchen Shi. "Influence maximization: Near-optimal time complexity meets practical efficiency." Proceedings of the 2014 ACM SIGMOD international conference on Management of data. 2014.
