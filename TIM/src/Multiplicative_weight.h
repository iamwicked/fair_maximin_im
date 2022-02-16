//typedef double (*pf)(int,int);
#include <math.h> 
#include "head.h"
#include <map>
#include <random>
#include <numeric>
#include <algorithm>

class Multiplicative_weight: public TimGraph
{
    public:
        int i_multi = 1;
        Multiplicative_weight(string folder, string graph_file):TimGraph(folder, graph_file){
        }

        map<int, double>sigma_v(vector<int> S){
            // Compute sigma_v(S) for all nodes v.
            // Parameters
            // ----------
            // S : list of nodes
            // Returns
            // -------
            // pro_v : dict
            //       Dictionary with keys being nodes and values probabilities of being reached of nodes from S. 

            sort(S.begin(), S.end());
            map<int, double> pro_v;
            for(int index=0; index<nodes.size(); index++){
                pro_v[nodes[index]] = 0.0;
                for(int r:rootR[nodes[index]]){
                    sort(hyperGT[r].begin(), hyperGT[r].end());
                    vector<int> intersect;   
                    set_intersection(hyperGT[r].begin(), hyperGT[r].end(), S.begin(),
                    S.end(), back_inserter(intersect));
                    if(intersect.size()>0){
                       pro_v[nodes[index]] += 1.0;
                    }
                }
                pro_v[nodes[index]] = pro_v[nodes[index]] * (double) n/(hyperGT.size() * weights[index]);
            }
            return pro_v;
       }
       
        map<string, double>sigma_C(vector<int> S){
            // Compute sigma_C(S) for all communities C.
            // Parameters
            // ----------
            // S : list of nodes
            // Returns
            // -------
            // comm_probs : dict
            //       Dictionary with keys being communities and values average probability of 
            //       being reached of nodes in communities from S. 

            map<int, double>node_probs = sigma_v(S);
            map<string, double>comm_probs;
            for(auto C:communities){
                comm_probs[C] = 0.0;
            }
            for (int v:nodes){
                for (auto C:node_communities[v]){
                        comm_probs[C] += (double)(node_probs[v] / communities_node[C].size());
                }
            }
            return comm_probs;
        }

        map<vector<int>, double> set_based(double epsilon, double eps=0.1, string set_or_node = "set"){
            // Call multi_weight for the set case in the paper.
            // Parameters
            // ----------
            // epsilon : float
            //      epsilon of (eps, delta)-approximation (for tim implementation).
            // eps : float
            //      eps of multi_weight.
            // set_or_node: str
            //      string showing that multi_weight should be executed for the set case.

            return multi_weight(epsilon, eps, set_or_node = "set");
        }

        map<vector<int>, double> node_based(double epsilon, double eps=0.1, string set_or_node = "node"){
            // Call multi_weight for the node case in the paper.
            // Parameters
            // ----------
            // epsilon : float
            //      epsilon of (eps, delta)-approximation (for tim implementation).
            // eps : float
            //      eps of multi_weight.
            // set_or_node: str
            //      string showing that multi_weight should be executed for the node case.

            return multi_weight(epsilon, eps, set_or_node = "node");
        }

        map<vector<int>, double> multi_weight(double epsilon, double eps, string set_or_node = "set"){
            // Implement multi_weight for the set and node based problems in the paper.
            // Parameters
            // ----------
            // epsilon : float
            //      epsilon of (eps, delta)-approximation (for tim implementation)     
            // eps : float
            //      eps of multi_weight.
            // set_or_node: str
            //      string showing that multi_weight should be executed for the set or node case.
            // Returns
            // -------
            // p : Dictionary with keys sets and values probabilities.
            //     in the node case the size of the sets is 1.

            map<string, double>z;
            map<string, double>s;

            for(auto C:communities){
                z[C] = 1.0;
                s[C] = 0.0;
            }

            map<vector<int>, double>p;
            i_multi = 1;
            double primal = -INFINITY;
            double dual = INFINITY;

            while(true)
            {
                double sum_z = 0;
                for(auto item: z){
                    sum_z += item.second;
                }

                double sum_weights = 0.0;
                for(int index=0; index<nodes.size(); index++){
                    double temp = 0.0;
                    for(auto c: node_communities[nodes[index]]){
                        temp += (double) (z[c] / communities_node[c].size());
                    }
                    weights[index] = (double) temp;
                    sum_weights += (double)weights[index];
                }
 
                for(int v:nodes){
                    weights[v] *= (double) (n/sum_weights);
                }
                EstimateOPT(epsilon, true);
                vector<int>oracle_solution = seedSet;
                double oracle_value = spread;

                if((double)(oracle_value * (sum_weights / nodes.size()) / sum_z) > 1.0000001){ 
                    cout << "_____________________________"<<endl;
                    for(auto item: z){
                        cout<< item.second<<endl;
                    }
                    cout<< "oracle_value, sum(z.values())" <<endl;
                    cout << oracle_value<< "," << sum_z<< endl;
                    cout << (double)(oracle_value * (sum_weights / nodes.size()) / sum_z) <<endl;
                    cout << "_____________________________"<<endl;
                }
                ASSERT((double)(oracle_value * (sum_weights/nodes.size()) / sum_z) <= 1.0000001); 
 
                if(set_or_node == "set"){
                    vector<int>difference;
                    sort(oracle_solution.begin(), oracle_solution.end());
                    int ii = 0;
                    for(auto item:p){
                        if(oracle_solution == item.first){
                            p[oracle_solution] += 1;
                            ii += 1;
                            break;
                        }
                    }
                    if (ii == 0){
                        p[oracle_solution] = 1;
                    }
                }
                else if(set_or_node == "node"){
                    int jj = 0;
                    for(int v:oracle_solution){
                        jj = 0;
                        for(auto item:p){
                            if(find(item.first.begin(), item.first.end(), v) != item.first.end()){
                                p[item.first] += 1;
                                jj += 1;
                                break;
                            }
                        }
                        if(jj == 0){
                            p[{v}] = 1;
                        }
                    }
                }
                else{
                    cout <<"Error: Unknown option."<< set_or_node<<endl;
                }

                dual = min((double)dual, (double)oracle_value * (sum_weights / nodes.size()) / sum_z);

                map<string, double> community_prob;
                community_prob = sigma_C(oracle_solution);

                for (auto C:communities){
                    z[C] *= (double)(1 - (double)(eps * community_prob[C]));
                }

                for (auto C:communities){
                    s[C] = (double)(i_multi - 1) / i_multi * s[C] + (double)1 / i_multi * community_prob[C];
                }

                auto min_s = *min_element(s.begin(), s.end(),
                    [](decltype(s)::value_type& l, decltype(s)::value_type& r) -> bool { return l.second < r.second; });
                primal = min_s.second;
                if(primal >= (1 - eps) * dual){
                    break;
                }
                i_multi += 1;
            }

            for(auto item:p){
                p[item.first] = (double)item.second/i_multi;
            }
            return p; 
        }
};
