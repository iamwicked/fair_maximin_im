#define HEAD_INFO
#include "sfmt/SFMT.h"
#include "head.h"
#include "memoryusage.h"
#include "graph.h"
#include "Multiplicative_weight.h"
#include <set>
#include <chrono>
#include <iomanip>

void run(Multiplicative_weight & x, string dataset, int k, double epsilon, string model, string setting, string file_name){
    // Infrastructure for the execution of the experiments
    // Parameters
    // ----------
    // dataset : str
    //     path to the dataset to be executed
    // k : int
    //     budget (seed set size)
    // epsilon : float
    //     epsilon of (eps, delta)-approximation (for tim implementation).
    // model : str
    //     underlying diffusion model (IC or LT)
    // setting : str
    //     string showing which function should be executed
    // file_name : 
    //     name of the file the results should be written

    cout << "dataset:" << dataset << " k:" << k << " epsilon:"<< epsilon << " model:" << model << " setting:" << setting << " file_name:" << file_name << endl;
    x.k=k;
    x.setting = setting;

    if(model=="IC")
        x.setInfuModel(InfGraph::IC);
    else if(model=="LT")
        x.setInfuModel(InfGraph::LT);
    else
        ASSERT(false);

    map<vector<int>, double>p_x;
    vector<int> S;

    auto start = chrono::high_resolution_clock::now();
    ios_base::sync_with_stdio(false);
    if(setting == "set"){   // for the set case
        p_x = x.set_based(epsilon); 
    }
    else if (setting == "node"){   // for the node case
        p_x = x.node_based(epsilon);
    }
    else if (setting == "tim"){
        x.EstimateOPT(epsilon, true);
        S = x.seedSet;
    }
    else{
        ASSERT(false);
    }
        
    auto end = chrono::high_resolution_clock::now();
    double time_taken = chrono::duration_cast<chrono::nanoseconds>(end - start).count();
    time_taken *= 1e-9;

    ofstream out_file;
    out_file.open(file_name);
    
    if (setting == "set" or setting == "node"){
        out_file << time_taken << endl;
        for(auto S_pro:p_x){
            for (int s:S_pro.first){
                out_file << s << " ";
            }
            out_file << fixed << setprecision(20) << S_pro.second << endl;
        }
    }
    else if (setting == "tim"){
        out_file << time_taken << endl;
        int i = 0;
        for(int s:S){
            i += 1;
            if(i < S.size()){
                out_file << s << " ";
            }
            else{
                out_file << s << endl;
            }
        }
    }
    else{
        ASSERT(false);
    }  
    Counter::show();
}
void parseArg(int argn, char ** argv)
{
    string file_name=""; 
    string setting=""; 
    string dataset="";
    double epsilon=0;
    string model="";
    int k=0;

    for(int i=0; i<argn; i++)
    {
        if(argv[i]==string("-file_name")) file_name=string(argv[i+1]);
        if(argv[i]==string("-setting")) setting=string(argv[i+1]);
        if(argv[i]==string("-dataset")) dataset=string(argv[i+1])+"/";
        if(argv[i]==string("-epsilon")) epsilon=atof(argv[i+1]);
        if(argv[i]==string("-k")) k=atoi(argv[i+1]);
        if(argv[i]==string("-model")) {
            if(argv[i+1]==string("LT"))
            {
                model=argv[i+1];
            }
            else if(argv[i+1]==string("IC"))
            {
                model=argv[i+1];
            }
            else
                ExitMessage("model should be IC or LT");
        }
    }
    if (dataset=="")
        ExitMessage("argument dataset missing");
    if (k==0)
        ExitMessage("argument k missing");
    if (epsilon==0)
        ExitMessage("argument epsilon missing");
    if (model=="")
        ExitMessage("argument model missing");
    if (setting=="")                             
        ExitMessage("argument setting missing");
    if (file_name=="")                             
        ExitMessage("argument file_name missing");
 
    string graph_file;
    if(model=="IC")
        graph_file=dataset + "graph_ic.txt"; 
    else if(model=="LT")
        graph_file=dataset + "graph_lt.inf";

    Multiplicative_weight x(dataset, graph_file);  
    run(x, dataset, k, epsilon, model, setting, file_name);
}

int main(int argn, char ** argv)
{
    OutputInfo info(argn, argv);
    parseArg( argn, argv );
}
