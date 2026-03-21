#include <cassert>
#include <queue>
#include <iostream>
#include <utility>
#include <vector>
#include <list>
#include <map>
#include <unordered_map>
#include <array>
#include <ctime>
#include <bitset>
#include <chrono>
#include <algorithm>
#include <random>
#include <unordered_set>

#include "gate.h"
#include "QASMReader.h"
#include "CouplingGraph.h"
#include "WeightedGate.h"

// #define DEBUG

using namespace std;
#define RANDOM_INIT

// ablation
// scoring function
// #define ABLATION_BRS

// initial mapping
// #define ABLATION_MODIFY

// tie-breaking rule
// #define ABLATION_STEP

// #define ABLATION_BRIDGE

// statistical info
// #define STAT

// params
// attention: params defined below will be replaced by the input command, 
// the initial value just show the default value
int IN_LIMIT = 5; // beta in paper, 5 is the default value
int MAX_LOOK_SIZE = 4; // varphi-1 in paper, 4 is the default value
int C_MAX = 500; // alpha in paper, 500 is the default value

using Key = pair<vector<int>, vector<gate>>;
map<Key,array<int,2>> flat_map;
map<vector<int>,int> flat_map2;

// Log file
ofstream logout;

// global random seed, defined in input command
int global_algorithm_seed;
mt19937 global_random_device;

string selected_device_enum;
enum Device { TOKYO,Q16,GUADALUPE,ROCHESTER,X19,SYCAMORE,torino };
map<string,Device> device_enum_map={ 
    {"TOKYO",TOKYO},{"GUADALUPE",GUADALUPE},{"ROCHESTER",ROCHESTER},
    {"SYCAMORE",SYCAMORE},{"torino",torino},
};

chrono::system_clock::time_point start_time;
chrono::system_clock::time_point current_time;
// cutoff time with s
const double cutoff_time=3600;

int total_bridge_count;
int best_ans_bridge_count;
int cur_ans_bridge_count;

#ifdef STAT
int total_decision_count;
// operations evaluated
int num_ops;
// meet tie count
int ties;
// bridges evaluated
int bris;
#endif

set<int> exe_st;

unordered_set<long long> doneOps;

auto makeKey = [](int x, int y) {
    return (static_cast<long long>(x) << 32) | (static_cast<unsigned int>(y));
};

// debug function
void debugMapping(const vector<int>& qubits){
    logout<<"qubits_mapping:"<<endl;
    for(int i=0;i<qubits.size();i++){
        logout<<i<<"->"<<qubits[i]<<endl;
    }
    logout<<"-------------------------------"<<endl;
}

// mapping function
void getQubitsSequence(vector<int>&qubits)
{
    #ifdef RANDOM_INIT
        shuffle(qubits.begin(),qubits.end(),global_random_device);
        return ;
    #endif
    logout<<"[ERROR] we should use random init."<<endl;
    exit(-1);
}

// use shuffle to perurbate the mapping
void PerurbationByShuffle(
    vector<int> & qubits, 
    vector<int> & locations,
    int nlocs)
{
    getQubitsSequence(qubits);

    for(int i=0;i<nlocs;i++)
    {
        if(qubits[i] == -1)
            continue;
        locations[qubits[i]] = i;
    }
}

// use randomnization to get initial mapping
void randInitMapping(vector<int> & qubits, int nlocs, int nqubits)
{
    for(int i = 0;i < nqubits; i++)
        qubits[i] = i;
}

/**
* Partition the logic circuit into n gate arrays,
* lines[i] holds all the gate along the ith qubit.
*/
void genLines(vector<list<gate>> & lines, vector<gate> & gateList)
{
    for(gate g1:gateList)
    {
        int c1 = g1.control;
        int t1 = g1.target;
        lines[c1].emplace_back(g1);
        lines[t1].emplace_back(g1);
    }
}

// find the active gates in current logic circuit
void findActivegates(vector<list<gate>> & lines,vector<gate> & activeGates,vector<gate> & phyCir,set<int> & activeLines, vector<int> locations, int ** dist)
{
    while(!activeLines.empty())
    {
        int i = *activeLines.begin();
        if(lines[i].empty())
        {
            activeLines.erase(i);
            continue;
        }
        // in one gate，i could be control, or could betarget
        // j is the other qubit
        gate g1 = lines[i].front();
        int c = g1.control;
        int t = g1.target;
        int j = i^g1.cxt;

        gate g2 = lines[j].front();
        if(g1==g2)
        {
            int locOfC = locations[c];
            int locOfT = locations[t];
            //executable gate
            if(dist[locOfC][locOfT] == 1)
            {
                lines[i].pop_front();
                lines[j].pop_front();
                phyCir.emplace_back(locOfC,locOfT, g1.weight);
                activeLines.emplace(j);
            }
            else
            {
                //need SWAP
                activeGates.emplace_back(g1);
                activeLines.erase(i);
                activeLines.erase(j);
            }
        }
        else activeLines.erase(i);
    }
}

void findActivegates_pre(vector<list<gate>> & lines,vector<gate> & activeGates,set<int> & activeLines, vector<int> locations, int ** dist)
{
    while(!activeLines.empty())
    {
        int i = *activeLines.begin();
        if(lines[i].empty())
        {
            activeLines.erase(i);
            continue;
        }

        gate g1 = lines[i].front();
        int c = g1.control;
        int t = g1.target;
        int j = i^g1.cxt;

        gate g2 = lines[j].front();
        if(g1==g2)
        {
            int locOfC = locations[c];
            int locOfT = locations[t];

            activeGates.emplace_back(g1);
            activeLines.erase(i);
            activeLines.erase(j);
        }
        else activeLines.erase(i);
    }
}
// use additonal variable to recover the state 
void findActivegates_recover(vector<list<gate>> & lines,vector<gate> & activeGates,vector<gate> & phyCir,set<int> & activeLines, vector<int> locations, int ** dist,vector<list<gate>> & re_line)
{
    while(!activeLines.empty())
    {
        int i = *activeLines.begin();
        if(lines[i].empty())
        {
            activeLines.erase(i);
            continue;
        }

        gate g1 = lines[i].front();
        int c = g1.control;
        int t = g1.target;
        int j = i^g1.cxt;

        gate g2 = lines[j].front();
        if(g1==g2)
        {
            int locOfC = locations[c];
            int locOfT = locations[t];
            //executable gate
            if(dist[locOfC][locOfT] == 1)
            {
                lines[i].pop_front();
                lines[j].pop_front();
                phyCir.emplace_back(locOfC,locOfT, g1.weight);

                re_line[i].emplace_back(g1);
                re_line[j].emplace_back(g1);
                activeLines.emplace(j);
            }
            else
            {
                activeGates.emplace_back(g1);
                activeLines.erase(i);
                activeLines.erase(j);
            }
        }
        else activeLines.erase(i);
    }
}

// calculate the number of executed gates and the number of front gates
// and log the executed gates
int numOfExecutedAndFrontGate_3_3_log(vector<list<gate>> & lines, const vector<int> & locations, int ** dist)
{
    int exeNum = 0;
    #ifdef DEBUG
    int gateNum = 0;
    set<pair<int,int>> frontGate;
    #endif
    set<int> activeLines;
    //Save gate deleted from lines for recovery
    vector<gate> executedGates;

    for(int i=0;i<lines.size();i++)
        activeLines.emplace(i);

    while(!activeLines.empty())
    {
        int i = *activeLines.begin();
        if(lines[i].empty())
        {
            activeLines.erase(i);
            continue;
        }
        gate g1 = lines[i].front();
        int c = g1.control;
        int t = g1.target;
        int j = i^g1.cxt;

        gate g2 = lines[j].front();
        if(g1 == g2) //active gate
        {
            int locOfC = locations[c];
            int locOfT = locations[t];
            //executable gate
            if(dist[locOfC][locOfT] == 1)
            {
                // logout<<"[INFO] Executed gate: "<<c<<"->"<<t<<" "<<g1.weight<<endl;
                exe_st.emplace(g1.weight);
                exeNum++;
                lines[i].pop_front();
                lines[j].pop_front();
                activeLines.emplace(j);
                executedGates.emplace_back(g1);
            }
            else
            {
                #ifdef DEBUG
                frontGate.emplace(min(c,t),max(c,t));
                gateNum++;
                #endif
                
                activeLines.erase(i);
                activeLines.erase(j);
            }
        }
        else activeLines.erase(i);
    }
    #ifdef DEBUG
    if(gateNum!=frontGate.size())
    {
        logout<<"[ERROR] gateNum!=frontGate.size()"<<endl;
        exit(-1);
    }
    #endif
    //Recover lines
    for(auto it=executedGates.rbegin();it!=executedGates.rend();++it) 
    {
        const gate& g1 = *it;
        int c = g1.control;
        int t = g1.target;
        lines[c].emplace_front(g1);
        lines[t].emplace_front(g1);
    }
    return exeNum;
}

// remove the executed gates from the logic circuit
void checkActivegates(vector<list<gate>> & lines,vector<gate> & activeGates,vector<gate> & phyCir,set<int> & activeLines, vector<int> locations, int ** dist)
{
    doneOps.clear();
    for(auto it = activeGates.begin(); it!=activeGates.end(); )
    {
        gate g = * it;
        int c = g.control;
        int t = g.target;
        int locOfC = locations[c];
        int locOfT = locations[t];
        //execute gate
        if(dist[locOfC][locOfT] == 1)
        {
            lines[c].pop_front();
            lines[t].pop_front();

            phyCir.emplace_back(locOfC,locOfT,g.weight);
            activeLines.emplace(c);
            activeLines.emplace(t);
            it = activeGates.erase(it);
        }
        else
            it++;
    }
}

// remove the executed gates from the logic circuit
void checkActivegates_recover(vector<list<gate>> & lines,vector<gate> & activeGates,vector<gate> & phyCir,set<int> & activeLines, vector<int> locations, int ** dist,vector<list<gate>> & re_line)
{
    for(auto it=activeGates.begin();it!=activeGates.end();)
    {
        gate g = * it;
        int c = g.control;
        int t = g.target;
        int locOfC = locations[c];
        int locOfT = locations[t];
        //execute gate
        if(dist[locOfC][locOfT] == 1)
        {
            lines[c].pop_front();
            lines[t].pop_front();

            re_line[c].emplace_back(g);
            re_line[t].emplace_back(g);

            phyCir.emplace_back(locOfC,locOfT,g.weight);
            activeLines.emplace(c);
            activeLines.emplace(t);

            it=activeGates.erase(it);
        }
        else it++;
    }
}

// calculate the function first section of SWAP gate
array<int,2> effOfSWAP(int s1,int s2, gate g1, vector<int> locations,  int ** dist)
{
    //two logical qubits of g1
    int i1 = g1.control;
    int j1 = g1.target;

    //two physical qubits of g1 before SWAP gate
    int locOfI1 = locations[i1];
    int locOfJ1 = locations[j1];
    if(locOfI1 > locOfJ1)
        swap(locOfI1, locOfJ1);

    #ifdef DEBUG
    if(locOfI1==-1 || locOfJ1==-1)
    {
        logout<<"[BUG] error in effOfSWAP, active gate must have two physical qubits"<<endl;
        exit(-1);
    }
    #endif // DEBUG

    int oriDist = dist[locOfI1][locOfJ1];
    int sufDist = oriDist;
    
    if(s1 == locOfI1 && s2 == locOfJ1)
        ;
    else if(s1 == locOfI1) // s1 swaped to locOfJ1
        sufDist = dist[locOfJ1][s2];
    else if(s1 == locOfJ1)
        sufDist = dist[locOfI1][s2];
    else if(s2 == locOfI1)
        sufDist = dist[s1][locOfJ1];
    else if(s2 == locOfJ1)
        sufDist = dist[s1][locOfI1];
        
    int bridgeBonus=0;
    if(oriDist!=2&&sufDist==2)
        bridgeBonus=1;
    else if(oriDist==2&&sufDist==3)
        bridgeBonus=-1;

    return {oriDist-sufDist,bridgeBonus};
}

/**
* Generate relevant SWAPs or Cnots which have at least one common qubit with some active gate
*/
void genRelevantGates(const CouplingGraph * p_cg,const vector<gate> & activeGates,set<weightedGate> & candiGates,const  vector<int> & locations,int ** dist )
{
    for(gate g1:activeGates)
    {
        int c = g1.control;
        int t = g1.target;

        int cAndT[2] = {locations[c], locations[t]};
        int oriDst = dist[locations[c]][locations[t]];
        int ori = locations[c]^locations[t];
        //bridge
        #ifndef ABLATION_BRIDGE
        if(oriDst == 2)
        {
            candiGates.emplace("C",cAndT[0],cAndT[1]);
            #ifdef STAT
            bris++;
            #endif
        }
        #endif

        //for all possible SWAPs on c and t
        for(int locOfC:cAndT)
        {
            for(int s1:p_cg->adjList[locOfC])
            {
                if(dist[s1][ori^locOfC]>oriDst)
                    continue;

                if(locOfC < s1)
                    candiGates.emplace("S",locOfC,s1);
                else
                    candiGates.emplace("S",s1,locOfC);
            }
        }
    }
}

void genRelevantGates_pre(const CouplingGraph * p_cg,const vector<gate> & activeGates,set<weightedGate> & candiGates,const  vector<int> & locations,int ** dist )
{
    for(gate g1:activeGates)
    {
        int c = g1.control;
        int t = g1.target;

        int cAndT[2] = {locations[c], locations[t]};
        int oriDst = dist[locations[c]][locations[t]];
        int ori = locations[c]^locations[t];

        //for all possible SWAPs on c and t
        for(int locOfC:cAndT)
        {
            for(int s1:p_cg->adjList[locOfC])
            {
                if(locOfC < s1)
                    candiGates.emplace("S",locOfC,s1);
                else
                    candiGates.emplace("S",s1,locOfC);
            }
        }
    }
}
/**
* Get the first nth layers of logical circuit
*/
vector<vector<gate>> getFrontLayers(vector<list<gate>> & lines, int layerDepth)
{
    vector<vector<gate>> layers;
    int nqubits = lines.size();
    vector<int> layer(nqubits,-1); //depth of each line during iteration
    vector<list<gate>::iterator> itOfLines(nqubits);
    for(int i=0;i<nqubits;i++)
        itOfLines[i] = lines[i].begin();

    bool isAllLineEmpty = false;
    bool isAllOverWindow = false;
    while(!isAllLineEmpty && !isAllOverWindow)
    {
        isAllLineEmpty = true;
        isAllOverWindow = true;
        for(int i=0;i<nqubits;i++)
        {
            if(itOfLines[i]!=lines[i].end() && layer[i]<layerDepth)
            {
                isAllLineEmpty = false;
                gate g1 = * itOfLines[i];
                int c = g1.control;
                int t = g1.target;
                if(*itOfLines[c] == *itOfLines[t])//same gate
                {
                    int layerOfG1 = max(layer[c],layer[t])+1;
                    layer[c] = layerOfG1;
                    layer[t] = layerOfG1;
                    if(layerOfG1 <= layerDepth)
                    {
                        if(layers.size() <= layerOfG1)
                            layers.emplace_back();
                        layers[layerOfG1].emplace_back(g1);
                        itOfLines[c]++;
                        itOfLines[t]++;
                        isAllOverWindow = false;
                    }
                }
            }
        }
    }
    return layers;
}

/**
* Update the mapping between physical and logical qubits with a SWAP gate
*/
void updateMapping(int q1,int q2, vector<int> & qubits, vector<int> & locations )
{
    //Two logical qubits of swapGate
    int lFirst = qubits[q1];
    int lSecond = qubits[q2];

    //Update qubits array
    swap(qubits[q1],qubits[q2]);

    //update locations array
    if(lFirst != -1) // corresponding to a unoccupied physical qubit
        locations[lFirst] = q2;
    if(lSecond != -1)
        locations[lSecond] = q1;
    // debugMapping(qubits);
}

int numOfExecutedAndFrontGate(vector<list<gate>> & lines, const vector<int> & locations, int ** dist)
{
    int exeNum = 0;
    #ifdef DEBUG
    int gateNum = 0;
    set<pair<int,int>> frontGate;
    #endif
    set<int> activeLines;
    //Save gate deleted from lines for recovery
    vector<gate> executedGates;

    for(int i=0;i<lines.size();i++)
        activeLines.emplace(i);

    while(!activeLines.empty())
    {
        int i = *activeLines.begin();
        if(lines[i].empty())
        {
            activeLines.erase(i);
            continue;
        }
        gate g1 = lines[i].front();
        int c = g1.control;
        int t = g1.target;
        int j = i^g1.cxt;

        gate g2 = lines[j].front();
        if(g1 == g2) //active gate
        {
            int locOfC = locations[c];
            int locOfT = locations[t];
            //executable gate
            if(dist[locOfC][locOfT] == 1)
            {
                exeNum++;
                lines[i].pop_front();
                lines[j].pop_front();
                activeLines.emplace(j);
                executedGates.emplace_back(g1);
            }
            else
            {
                #ifdef DEBUG
                frontGate.emplace(min(c,t),max(c,t));
                gateNum++;
                #endif
                
                activeLines.erase(i);
                activeLines.erase(j);
            }
        }
        else activeLines.erase(i);
    }
    #ifdef DEBUG
    if(gateNum!=frontGate.size())
    {
        logout<<"[ERROR] gateNum!=frontGate.size()"<<endl;
        exit(-1);
    }
    #endif
    //Recover lines
    for(int i=executedGates.size()-1;i>=0;i--)
    {
        gate g1 = executedGates[i];
        int c = g1.control;
        int t = g1.target;
        lines[c].emplace_front(g1);
        lines[t].emplace_front(g1);
    }
    return exeNum;
}

array<int,3> costOfInsertedCnot_L0_Nexe_Frontgate(vector<list<gate>> & lines, int q1, int q2, 
                                           const vector<int> & qubits,const vector<int> & locations,  int ** dist)
{
    int c = qubits[q1];
    int t = qubits[q2];
    gate g1 = lines[c].front();
    lines[c].pop_front();
    gate g2 = lines[t].front();
    lines[t].pop_front();
    #ifdef DEBUG
    if(!(g1 == g2) || !(g1.control==c && g1.target==t))
    {
        logout<<"[BUG] error in 980" <<endl;
        logout<<g1.control<<","<<g1.target<<","<<g1.weight<<endl;
        logout<<g2.control<<","<<g2.target<<","<<g2.weight<<endl;
        exit(-1);
    }
    #endif

    int gateCost = numOfExecutedAndFrontGate(lines,locations,dist);
    lines[c].emplace_front(g1);
    lines[t].emplace_front(g1);
    return {1,gateCost+1,0};
}

array<int,3> costOfInsertedSwap_L0_Nexe_Frontgate_3_3_log(vector<list<gate>> & lines, int q1, int q2,const vector<vector<gate>> & layers, 
                                            vector<int> qubits,vector<int> locations,  int ** dist)
{
    bool hasExecutable = false;
    int funcV1 = 0;
    int funcV2 = 0, funcV3 = 0;

    // judge exec
    for(gate g1:layers[0])
    {
        auto v=effOfSWAP(q1, q2, g1, locations, dist);
        int effectOnG1 = v[0];
        // dis from 2 to 1
        if(effectOnG1 == 1)
        {
            int loc1_g1 = locations[g1.control];
            int loc2_g1 = locations[g1.target];
            int dist_g1 = dist[loc1_g1][loc2_g1];
            if(dist_g1 == 2)
                hasExecutable = true;
        }
        
        funcV1 += effectOnG1;
        funcV3 += v[1];
    }

    // cal funcv1
    for(int i=1;i<layers.size();i++)
    {
        for(gate g1:layers[i])
        {
            auto v=effOfSWAP(q1, q2, g1, locations, dist);
            funcV1 += v[0];
            funcV3 += v[1];
        }
    }

    if(hasExecutable)
    {
        updateMapping(q1,q2,qubits,locations);
        funcV2 = numOfExecutedAndFrontGate_3_3_log(lines,locations,dist);

        for(auto layer:layers)
        {
            for(auto g:layer)
            {
                if(exe_st.find(g.weight)!=exe_st.end())
                {
                    auto v= effOfSWAP(q1, q2, g, locations, dist);
                    // no change or 2 to 1, then sub value 1 to 2 
                    funcV1+=v[0];
                }
            }
        }
        exe_st.clear();
    }

    return {funcV1,funcV2,funcV3};
}

array<int,2> costOfInsertedSwap_L0_Nexe_Frontgate_pre(int q1, int q2,const vector<vector<gate>> & layers, 
                                            vector<int> locations,  int ** dist)
{
    bool hasExecutable = false;
    int funcV1 = 0;
    int funcV3 = 0;

    // cal funcv1
    for(int i=0;i<layers.size();i++)
    {
        for(gate g1:layers[i])
        {
            auto v=effOfSWAP(q1, q2, g1, locations, dist);
            funcV1 += v[0];
            funcV3 += v[1];
        }
    }

    return {funcV1,funcV3};
}
/**
*Insert a SWAP gate into the resultant circuit
*/
void insertSwap(int p1,int p2, vector<gate> & finalList)
{
    int num = -1;
    finalList.emplace_back(p1,p2,num);
    finalList.emplace_back(p2,p1,num);
    finalList.emplace_back(p1,p2,num);
}

/**
*Insert a bridge cnot gate into the resultant circuit
*/
void insertBridgeCnot(int c, int t, int via, vector<gate> & finalList,int number)
{
    int trans_number=-number-2;//0->-2,1->-3,2->-4...
    // or 1 2|0 1|1 2|0 1
    finalList.emplace_back(c,via,trans_number);
    finalList.emplace_back(via,t,trans_number);
    finalList.emplace_back(c,via,trans_number);
    finalList.emplace_back(via,t,trans_number);
}

// get the valued gate list
auto getValuedGateList_L0_Nexe_FrontGate_3_3_log(const set<weightedGate> & candiSwap,
                            vector<list<gate>> & lines,
                            int maxLookSize,  const vector<int> & qubits,
                            const vector<int> & locations,  int ** dist)
{
    vector<vector<gate>> layers = getFrontLayers(lines, maxLookSize);

    using value_gate = pair<weightedGate,array<int,3>>;
    vector<value_gate> candidate_gate;
    for(auto gate1:candiSwap)
    {
        array<int,3> effOfGate1;
        if(gate1.type == "S")
            effOfGate1 = costOfInsertedSwap_L0_Nexe_Frontgate_3_3_log(lines,gate1.q1,gate1.q2,layers,qubits,locations,dist);
        else //gate1.type == "CNOT"
        {
            #ifdef ABLATION_BRIDGE

            #else
                effOfGate1 = costOfInsertedCnot_L0_Nexe_Frontgate(lines,gate1.q1,gate1.q2,qubits,locations,dist);
            #endif
        }
        candidate_gate.emplace_back(gate1,effOfGate1);
    }
    return candidate_gate;
}

auto getValuedGateList_L0_Nexe_FrontGate_pre(vector<vector<gate>> & layers,
                            const set<weightedGate> & candiSwap,
                            const vector<int> & locations,  int ** dist)
{

    using value_gate = pair<weightedGate,array<int,2>>;
    vector<value_gate> candidate_gate;
    for(auto gate1:candiSwap)
    {
        array<int,2> effOfGate1;
        effOfGate1 = costOfInsertedSwap_L0_Nexe_Frontgate_pre(gate1.q1,gate1.q2,layers,locations,dist);

        candidate_gate.emplace_back(gate1,effOfGate1);
    }
    return candidate_gate;
}
/**
* The mapping procedure by using swap and bridge cnot gate
* return the number of ancillary cnot gates
*/

void insert_SWAP_to_phy(weightedGate tar,vector<gate> & phyCir,int &nCnots,vector<int> & qubits,vector<int> & locations)
{
    insertSwap(tar.q1,tar.q2, phyCir);
    nCnots += 3;
    updateMapping(tar.q1,tar.q2, qubits, locations);
}

void insert_CNOT_to_phy(weightedGate tar,vector<gate> & phyCir,int &nCnots,vector<int> & qubits,vector<int> & locations,
                        CouplingGraph * p_cg,vector<list<gate>>  &lines,vector<gate> &activeGates,set<int> & activeLines)
{
    vector<int> linkedPath = p_cg->getPath(tar.q1,tar.q2);
    
    nCnots +=3;
    int c = qubits[tar.q1];
    int t = qubits[tar.q2];
    int w = lines[c].front().weight;

    #ifdef DEBUG
    gate g1 = lines[c].front();
    gate g2 = lines[t].front();
    #endif
    lines[c].pop_front();
    lines[t].pop_front();
    #ifdef DEBUG
    if(!(g1 == g2)){
        logout<<"[BUG] error in insert_CNOT_to_phy"<<endl;
        exit(-1);
    }
    #endif
    for(vector<gate>::iterator it = activeGates.begin(); it!=activeGates.end();it++ )
    {
        gate g2 = *it;
        if(c == g2.control && t == g2.target)
        {
            activeGates.erase(it);
            break;
        }
    }
    activeLines.emplace(c);
    activeLines.emplace(t);
    insertBridgeCnot(tar.q1,tar.q2,linkedPath[1],phyCir,w);
}

// insert CNOT gate, but also record the line situation in re_line, used for later recovery
void insert_CNOT_to_phy_recover(weightedGate tar,vector<gate> & phyCir,int & nCnots,vector<int> & qubits,vector<int> & locations,
                        CouplingGraph * p_cg,vector<list<gate>> & lines,vector<gate> & activeGates,set<int> & activeLines,
                        vector<list<gate>> & re_line)
{
    vector<int> linkedPath = p_cg->getPath(tar.q1,tar.q2);
    
    nCnots +=3;
    int c = qubits[tar.q1];
    int t = qubits[tar.q2];
    int w = lines[c].front().weight;
    
    #ifdef DEBUG
    gate g1 = lines[c].front();
    gate g2 = lines[t].front();
    #endif
    
    re_line[c].emplace_back(lines[c].front());
    lines[c].pop_front();

    re_line[t].emplace_back(lines[t].front());
    lines[t].pop_front();
    
    #ifdef DEBUG
    if(!(g1 == g2)){
        logout<<"[BUG] error in insert_CNOT_to_phy"<<endl;
        exit(-1);
    }
    #endif
    for(auto it = activeGates.begin(); it!=activeGates.end();it++ )
    {
        gate g2 = *it;
        if(c == g2.control && t == g2.target)
        {
            activeGates.erase(it);
            break;
        }
    }
    activeLines.emplace(c);
    activeLines.emplace(t);
    insertBridgeCnot(tar.q1,tar.q2,linkedPath[1],phyCir,w);
}

// recover the state after do first insert operation
void doRecoverWork(int &nCnots,int &beforeCnots,
                   set<int> &activeLines,
                   vector<gate> &activeGates,vector<gate> &tmpActiveGates,
                   vector<list<gate>> &lines,vector<list<gate>> &re_line,
                   int beforePhyCirSize,vector<gate> & phyCir,
                   weightedGate &gateInfo,
                   vector<int> &qubits,vector<int> &locations)
{
    // recover nCnots
    nCnots=beforeCnots;
    // recover activeLines(set to empty)
    activeLines.clear();
    // recover activeGates
    activeGates=tmpActiveGates;

    for(int i=0;i<lines.size();i++)
    {
        if(re_line[i].size())
        {
            if(re_line[i].size()==1)
                lines[i].emplace_front(*re_line[i].begin());
            else
                lines[i].splice(lines[i].begin(),re_line[i]);
        }
    }

    // recover phyCircle
    phyCir.resize(beforePhyCirSize);
    
    // if SWAP gate is inserted, recover qubits and locations
    if(gateInfo.type=="S")
    {
        updateMapping(gateInfo.q1,gateInfo.q2,qubits,locations);
    }
}

//final version
bool doTwoInsertOperate_new_func_fixed_ablation(vector<list<gate>> &lines,
                        vector<gate> &activeGates,set<int> & activeLines,
                        CouplingGraph * p_cg,vector<int> & qubits,vector<int> & locations,
                        int ** dist,vector<gate> & phyCir,int &nCnots,int seq){

    findActivegates(lines,activeGates,phyCir,activeLines,locations,dist);
    if(activeGates.empty()) // all gate have been mapped
        return true;

    // log nCnots status
    auto tryFind=flat_map.find({locations,activeGates});
    if(tryFind!=flat_map.end())
    {
        if(tryFind->second[seq]>nCnots)
            tryFind->second[seq]=nCnots;
        else
        {
            nCnots=-1;
            return true;
        }
    }
    else
    {
        array<int,2> log;
        log[seq]=nCnots,log[1-seq]=numeric_limits<int>::max();
        flat_map.emplace(Key(locations,activeGates),log);
    }

    //generate candidate swap and bridge gates
    set<weightedGate> candiGates;

    genRelevantGates(p_cg, activeGates, candiGates, locations, dist);

    using value_gate = pair<weightedGate,array<int,3>>;

    vector<value_gate> gateList = getValuedGateList_L0_Nexe_FrontGate_3_3_log(candiGates,lines,MAX_LOOK_SIZE,qubits,locations,dist);

    int bestCost = numeric_limits<int>::min();

    #ifndef ABLATION_STEP
    vector<pair<int,weightedGate>> bestGate;
    #else
    int bestGate;
    int bestCnt=0;
    #endif

    #ifdef STAT
    num_ops=0;
    #endif
    for(int i=0;i<gateList.size();i++)
    {
        auto curGate=gateList[i];

        auto gateInfo=curGate.first;
        auto gateValue=curGate.second;

        int weight=gateValue[0]+gateValue[1];

        if(weight<0)
            continue;
        
        if(gateInfo.type=="S")
        {
            if(doneOps.find(makeKey(gateInfo.q1,gateInfo.q2))!=doneOps.end())
                continue;
        }

        int beforeCnots=nCnots;
        int beforePhyCirSize=phyCir.size();
        vector<list<gate>> re_line(lines.size());

        auto tmpActiveGates=activeGates;

        if(gateInfo.type == "S")
        {
            insert_SWAP_to_phy(gateInfo,phyCir,nCnots,qubits,locations);
            if(gateValue[1])
                checkActivegates_recover(lines,activeGates,phyCir,activeLines,locations,dist,re_line);
        }
        else
        {
            insert_CNOT_to_phy_recover(gateInfo,phyCir,nCnots,qubits,locations,p_cg,lines,activeGates,activeLines,re_line);
            checkActivegates_recover(lines,activeGates,phyCir,activeLines,locations,dist,re_line);
        }
        
        findActivegates_recover(lines,activeGates,phyCir,activeLines,locations,dist,re_line);

        if(activeGates.empty())
        {
            bestCost=2*weight+gateValue[2];
            #ifdef ABLATION_BRS
            bestCost=weight;
            #endif

            #ifndef ABLATION_STEP
            bestGate.clear();
            bestGate.emplace_back(i,weightedGate("N"));
            #else
            bestGate=i;
            #endif

            doRecoverWork(nCnots,beforeCnots,activeLines,activeGates,tmpActiveGates,lines,re_line,beforePhyCirSize,phyCir,gateInfo,qubits,locations);
            break;
        }
        set<weightedGate> nextCandiGates;
        genRelevantGates(p_cg,activeGates,nextCandiGates,locations,dist);
        
        auto nextGateList=getValuedGateList_L0_Nexe_FrontGate_3_3_log(nextCandiGates,lines,MAX_LOOK_SIZE,qubits,locations,dist);
        for(auto nextGate:nextGateList)
        {
            auto nextGateInfo=nextGate.first;
            auto nextGateValue=nextGate.second;

            int weight2=nextGateValue[0]+nextGateValue[1];
            
            if(weight2<=0)
                continue;

            #ifdef STAT
            num_ops++;
            #endif
            int weight3=gateValue[2]+nextGateValue[2];
            int curCost=2*(weight+weight2)+weight3;

            #ifdef ABLATION_BRS
            curCost=weight+weight2;
            #endif

            #ifndef ABLATION_STEP
            if(curCost>bestCost)
            {
                bestCost=curCost;
                bestGate.clear();
                bestGate.emplace_back(i,nextGateInfo);
            }
            else if(curCost==bestCost)
            {
                bestGate.emplace_back(i,nextGateInfo);
            }
            #else
            if(curCost>bestCost)
            {
                bestCost=curCost;
                bestGate=i;
                bestCnt=1;
            }
            else if(curCost==bestCost)
            {
                bestCnt++;
                if(global_random_device()%bestCnt==0)
                    bestGate=i;
            }
            #endif
        }
        
        doRecoverWork(nCnots,beforeCnots,activeLines,activeGates,tmpActiveGates,lines,re_line,beforePhyCirSize,phyCir,gateInfo,qubits,locations);
    }

    #ifdef STAT
    total_decision_count++;
    #endif

    if(bestCost>0)
    {
        #ifndef ABLATION_STEP
        int bestInd  = 0;
        if(bestGate.size()>1)
        {
            int bestCnt  = 0;
            int bestCost3= numeric_limits<int>::min();
            // alt-3 trial
            for(int i=0;i<bestGate.size();i++)
            {
                auto gatePair=bestGate[i];

                auto curGate= gateList[gatePair.first];
                auto gateInfo=curGate.first;
                auto gateValue=curGate.second;

                // recover prepare
                int beforeCnots2=nCnots;
                int beforePhyCirSize2=phyCir.size();
                vector<list<gate>> re_line2(lines.size());
                auto tmpActiveGates2=activeGates;

                if(gateInfo.type == "S")
                {
                    insert_SWAP_to_phy(gateInfo,phyCir,nCnots,qubits,locations);
                    if(gateValue[1])
                        checkActivegates_recover(lines,activeGates,phyCir,activeLines,locations,dist,re_line2);
                }
                else
                {
                    insert_CNOT_to_phy_recover(gateInfo,phyCir,nCnots,qubits,locations,p_cg,lines,activeGates,activeLines,re_line2);
                    checkActivegates_recover(lines,activeGates,phyCir,activeLines,locations,dist,re_line2);
                }
                
                findActivegates_recover(lines,activeGates,phyCir,activeLines,locations,dist,re_line2);

                auto gateInfo2=gatePair.second;
                if(gateInfo2.type == "S")
                {
                    insert_SWAP_to_phy(gateInfo2,phyCir,nCnots,qubits,locations);
                    checkActivegates_recover(lines,activeGates,phyCir,activeLines,locations,dist,re_line2);
                }
                else if(gateInfo2.type == "C")
                {
                    insert_CNOT_to_phy_recover(gateInfo2,phyCir,nCnots,qubits,locations,p_cg,lines,activeGates,activeLines,re_line2);
                    checkActivegates_recover(lines,activeGates,phyCir,activeLines,locations,dist,re_line2);
                }
                
                findActivegates_recover(lines,activeGates,phyCir,activeLines,locations,dist,re_line2);
 
                if(activeGates.empty())
                {
                    bestInd=i;
                    bestCnt=1;
                    doRecoverWork(nCnots,beforeCnots2,activeLines,activeGates,tmpActiveGates2,lines,re_line2,beforePhyCirSize2,phyCir,gateInfo2,qubits,locations);
                    if(gateInfo.type=="S")
                    {
                        updateMapping(gateInfo.q1,gateInfo.q2,qubits,locations);
                    }

                    break;
                }

                // step3
                set<weightedGate> nextNextCandiGates;
                genRelevantGates(p_cg,activeGates,nextNextCandiGates,locations,dist);
                
                auto nextNextGateList=getValuedGateList_L0_Nexe_FrontGate_3_3_log(nextNextCandiGates,lines,MAX_LOOK_SIZE,qubits,locations,dist);
                int cost3=0;
                for(auto nextNextGate:nextNextGateList)
                {
                    auto nextNextGateInfo =nextNextGate.first;
                    auto nextNextGateValue=nextNextGate.second;

                    int weight3=nextNextGateValue[0]+nextNextGateValue[1];
                    
                    if(weight3<=0)
                        continue;

                    int weight4=nextNextGateValue[2];
                    cost3=max(cost3,2*weight3+weight4);
                    #ifdef ABLATION_BRS
                    cost3=max(cost3,weight3);
                    #endif
                }

                if(cost3>bestCost3)
                {
                    bestInd=i;
                    bestCnt=1;
                    bestCost3=cost3;
                }
                else if(cost3==bestCost3)
                {
                    // 1/count probability replace
                    bestCnt++;
                    if(global_random_device()%bestCnt==0)
                        bestInd=i;
                }

                doRecoverWork(nCnots,beforeCnots2,activeLines,activeGates,tmpActiveGates2,lines,re_line2,beforePhyCirSize2,phyCir,gateInfo2,qubits,locations);
                if(gateInfo.type=="S")
                {
                    updateMapping(gateInfo.q1,gateInfo.q2,qubits,locations);
                }
            }

            #ifdef STAT
            if(bestCnt>1)
                ties++;
            #endif
        }
        auto exeGate1=gateList[bestGate[bestInd].first].first;
        #else
        auto exeGate1=gateList[bestGate].first;
        #endif
        if(exeGate1.type == "S")
        {
            //log swap gate
            doneOps.insert(makeKey(exeGate1.q1,exeGate1.q2));

            insert_SWAP_to_phy(exeGate1,phyCir,nCnots,qubits,locations);
            #ifndef ABLATION_STEP
            int value=gateList[bestGate[bestInd].first].second[1];
            #else
            int value=gateList[bestGate].second[1];
            #endif

            if(value)
                checkActivegates(lines,activeGates,phyCir,activeLines,locations,dist);            
        }
        else
        {
            insert_CNOT_to_phy(exeGate1,phyCir,nCnots,qubits,locations,p_cg,lines,activeGates,activeLines);
            // to do bridge
            total_bridge_count++;
            cur_ans_bridge_count++;
            checkActivegates(lines,activeGates,phyCir,activeLines,locations,dist);
        }
    }
    else
    {
        int minLenOfFrontCnot = numeric_limits<int>::max();
        int controlOfMinLenGate;
        int targetOfMinlenGate;
        // shortest dist
        for(gate frontGate1:activeGates)
        {
            int locOfC = locations[frontGate1.control];
            int locOfT = locations[frontGate1.target];
            if(dist[locOfC][locOfT] < minLenOfFrontCnot)
            {
                minLenOfFrontCnot = dist[locOfC][locOfT];
                controlOfMinLenGate = locOfC;
                targetOfMinlenGate = locOfT;
            }
        }

        vector<int> shortestPath = p_cg->getPath(controlOfMinLenGate,targetOfMinlenGate);
        // path swap value

        set<weightedGate> pathGates;
        for(int i=0;i<shortestPath.size()-1;i++)
            pathGates.emplace("S",shortestPath[i],shortestPath[i+1]);

        vector<value_gate> gateList=getValuedGateList_L0_Nexe_FrontGate_3_3_log(pathGates,lines,MAX_LOOK_SIZE,qubits,locations,dist);
        int minWeight=numeric_limits<int>::max();
        int minSeq=0;
        for(int i=0;i<gateList.size();i++)
        {
            auto v=gateList[i].second;
            int w=2*(v[0]+v[1])+v[2];
            #ifdef ABLATION_BRS
            w=v[0]+v[1];
            #endif

            if(w<minWeight)
            {
                minWeight=w;
                minSeq=i;
            }
        }
        //insert swap except i

        //from sides to middle
        nCnots += 3*(gateList.size()-1);

        for(int i=0;i<minSeq;i++)
        {
            auto v=gateList[i].first;
            insertSwap(v.q1,v.q2, phyCir);
            updateMapping(v.q1,v.q2, qubits, locations);
        }

        for(int i=gateList.size()-1;i>minSeq;i--)
        {
            auto v=gateList[i].first;
            insertSwap(v.q1,v.q2, phyCir);
            updateMapping(v.q1,v.q2, qubits, locations);
        }
        
        checkActivegates(lines,activeGates,phyCir,activeLines,locations,dist);
    }
    
    return false;
}

int swapAndCnotBasedMapping_advance(
    int nqubits, vector<gate> & gateList, 
    CouplingGraph * p_cg, vector<gate> & phyCir, 
    int ** dist, vector<int> & qubits, 
    vector<int> & locations, bool rev_use)
{
    int nCnots = 0;
    set<int> activeLines;
    for(int i=0;i<nqubits;i++)
        activeLines.emplace(i);

    vector<list<gate>> lines(nqubits);
    genLines(lines,gateList);
    vector<gate> activeGates;

    int seq=0;
    if(rev_use) seq=1;

    while(!doTwoInsertOperate_new_func_fixed_ablation(lines,activeGates,activeLines,
        p_cg,qubits,locations,dist,phyCir,nCnots,seq));
    
    return nCnots;
}

bool doTwoInsertOperate_pre(vector<vector<gate>>& layers,vector<gate> &activeGates,set<int> & activeLines,
                        CouplingGraph * p_cg,vector<int> & qubits,vector<int> & locations,
                        int ** dist)
{
    auto tryFind=flat_map2.find(locations);
    if(tryFind!=flat_map2.end())
    {
        if(tryFind->second>0)
            tryFind->second=0;
        else
            return false;
    }
    else
        flat_map2.emplace(locations,0);

    set<weightedGate> candiGates;
    genRelevantGates_pre(p_cg, activeGates, candiGates, locations, dist);

    using value_gate = pair<weightedGate,array<int,2>>;

    vector<value_gate> gateList = getValuedGateList_L0_Nexe_FrontGate_pre(layers,candiGates,locations,dist);

    int bestCost = numeric_limits<int>::min();
    int bestCnt  = 0; // record the occurrence count of the maximum value
    int bestGate; //index of best

    for(int i=0;i<gateList.size();i++)
    {
        auto gateValue=gateList[i].second;

        int weight=gateValue[0];

        if(weight<=0)
            continue;
        
        int curCost=2*weight+gateValue[1];

        #ifdef ABLATION_BRS
        curCost=weight;
        #endif

        if(curCost>bestCost)
        {
            bestCost=curCost;
            bestGate=i;
            bestCnt=1;
        }
        else if(curCost==bestCost)
        {
            bestCnt++;
            if(global_random_device()%bestCnt==0)
                bestGate=i;
        }
    }

    if(bestCost>0)
    {
        auto exeGate1=gateList[bestGate].first;
        updateMapping(exeGate1.q1,exeGate1.q2, qubits, locations);
    }
    else
        return false;
    
    return true;
}

void preTreat(int nqubits, vector<gate> & gateList, CouplingGraph * p_cg,
              int ** dist, vector<int> & qubits, vector<int> & locations)
{
    set<int> activeLines;
    for(int i=0;i<nqubits;i++)
        activeLines.emplace(i);

    vector<list<gate>> lines(nqubits);
    genLines(lines,gateList);

    vector<gate> activeGates;
    findActivegates_pre(lines,activeGates,activeLines,locations,dist);

    vector<vector<gate>> layers=getFrontLayers(lines,2*MAX_LOOK_SIZE);
    while(doTwoInsertOperate_pre(layers,activeGates,activeLines,p_cg,qubits,locations,dist));
}

// check whether the physical circuit is valid
void answer_checker(int minCnotCount,vector<gate> & gateList,
                    vector<gate> & phyCirOfBest,const CouplingGraph & cg,
                    int nqubits,vector<int> & qubits,vector<int> & locations)
{
    for(int i=0;i<qubits.size();i++)
    {
        if(qubits[i] == -1)
            continue;
        locations[qubits[i]] = i;
    }

    if(gateList[0].weight>gateList.back().weight)
        reverse(gateList.begin(),gateList.end());
    
    if(minCnotCount + gateList.size() != phyCirOfBest.size())
    {
        logout<<"[ERROR] minCnotCount + gateList.size() != phyCirOfBest.size()"<<endl;
        logout<<"minCnotCount: "<<minCnotCount<<endl;
        logout<<"gateList.size(): "<<gateList.size()<<endl;
        logout<<"phyCirOfBest.size(): "<<phyCirOfBest.size()<<endl;
        logout<<"qubitsArrOfBest: "<<endl;
        for(auto i : phyCirOfBest)
            logout << i.control << " " << i.target << endl;
        
        exit(-1);
    }

    for(auto i:phyCirOfBest)
    {
        // logout<<i.weight<<' '<<i.control<<' '<<i.target<<endl;
        if(i.control == i.target)
        {
            logout<<"[ERROR] control == target"<<endl;
            exit(-1);
        }
        if(i.control<0||i.target<0)
        {
            logout<<"[ERROR] control or target < 0"<<endl;
            exit(-1);
        }
        if(i.control>=cg.getPositions()||i.target>=cg.getPositions())
        {
            logout<<"[ERROR] control or target >= nlocations"<<endl;
            exit(-1);
        }
        if(cg.dist[i.control][i.target] != 1)
        {
            logout<<"[ERROR] dist != 1"<<endl;
            exit(-1);
        }
    }

    vector<pair<int,int>> graph(gateList.size(),{-1,-1});
    vector<int> last_gate(nqubits,-1);
    unordered_map<int,bool> is_checked;
    is_checked[-1] = true;
    for(auto i:gateList)
    {
        int target = i.target;
        int control = i.control;
        int number = i.weight;
        graph[number].first = last_gate[target];
        graph[number].second = last_gate[control];
        last_gate[target] = number;
        last_gate[control] = number;
        // logout<<"[INFO] "<<number<<"<-"<<graph[number].first<<"<-"<<graph[number].second<<endl;
    }

    for(int i=0;i<phyCirOfBest.size();i++)
    {
        auto phy_tar = phyCirOfBest[i].target;
        auto phy_con = phyCirOfBest[i].control;
        auto number = phyCirOfBest[i].weight;
        // logout<<"[INFO] "<<number<<" "<<phy_tar<<" "<<phy_con<<endl;
        if(number==-1)
        {
            if(i+2>=phyCirOfBest.size())
            {
                logout<<"[ERROR] i+2>=phyCirOfBest.size()"<<endl;
                exit(-1);
            }
            swap(qubits[phy_con],qubits[phy_tar]);
            for(int i=0;i<qubits.size();i++)
            {
                if(qubits[i] == -1)
                    continue;
                locations[qubits[i]]=i;
            }
            i+=2;
            continue;
        }
        bool is_bridge = false;
        if(number < 0)
        {
            is_bridge=true;
            if(i+3>=phyCirOfBest.size())
            {
                logout<<"[ERROR] i+3>=phyCirOfBest.size()"<<endl;
                exit(-1);
            }
            number += 2;
            number = -number;
            if(cg.dist[locations[gateList[number].control]][locations[gateList[number].target]] != 2)
            {
                logout<<"[ERROR] dist != 2"<<endl;
                logout<<"number: "<<number<<endl;
                logout<<"control: "<<gateList[number].control<<" target: "<<gateList[number].target<<endl;
                logout<<"location of control: "<<locations[gateList[number].control]<<" location of target: "<<locations[gateList[number].target]<<endl;
                logout<<"cur mapping: "<<endl;
                debugMapping(qubits);
                logout<<"------------"<<endl;

                exit(-1);
            }
            i+=3;
            // continue;
        }
        if(is_checked[number])
        {
            logout<<"[ERROR] repeated gate in number: "<<number<<endl;
            logout<<"control: "<<phy_con<<" target: "<<phy_tar<<endl;
            exit(-1);
        }
        if(!is_checked[graph[number].first])
        {
            logout<<"[ERROR] cross error in number: "<<number;
            logout<<"control: "<<phy_con<<" target: "<<phy_tar<<endl;
            logout<<" its constrians is "<<graph[number].first<<" and "<<graph[number].first<<endl;
            logout<<"not checked: "<<graph[number].first<<endl;
            exit(-1);
        }
        if(!is_checked[graph[number].second])
        {
            logout<<"[ERROR] cross error in number: "<<number;
            logout<<"control: "<<phy_con<<" target: "<<phy_tar<<endl;
            logout<<" its constrians is "<<graph[number].first<<" and "<<graph[number].first<<endl;
            logout<<"not checked: "<<graph[number].first<<endl;
            exit(-1);
        }
        is_checked[number] = true;
        if(is_bridge)
            continue;

        auto origin_gate = gateList[number];
        if(qubits[phy_con] != origin_gate.control)
        {
            logout<<"[ERROR] control error in number: "<<number;
            logout<<"phy control: "<<phy_con<<" phy target: "<<phy_tar<<endl;
            logout<<"logi control: "<<origin_gate.control<<" logi target: "<<origin_gate.target<<endl;
            logout<<"qubits[phy_con]: "<<qubits[phy_con]<<endl;
            exit(-1);
        }
        if(qubits[phy_tar] != origin_gate.target)
        {
            logout<<"[ERROR] target error in number: "<<number;
            logout<<"phy control: "<<phy_con<<" phy target: "<<phy_tar<<endl;
            logout<<"logi control: "<<origin_gate.control<<" logi target: "<<origin_gate.target<<endl;
            logout<<"qubits[phy_tar]: "<<qubits[phy_tar]<<endl;
            exit(-1);
        }
        
        // logout<<"[INFO] over"<<endl;
    }
    logout<<"[INFO] check over. Circult has no fault."<<endl;
}

// Begin of algorithm
void EffectiveQMFramework(string circuitFileName, int inLimit=IN_LIMIT)
{
    // logout<<"[INFO] greedy version"<<endl;
    logout<<"[INFO] parameter: IN_LIMIT: "<<IN_LIMIT<<" MAX_LOOK_SIZE: "<<MAX_LOOK_SIZE<<" C_MAX: "<<C_MAX<<endl;
    vector<gate> gateList;  //logical circuit
    vector<gate> phyCir;   //physical circuit

    CouplingGraph cg;
    int ** dist;

    // Read qasm file into gate list
    // circuitFileName = "./"+circuitFileName;
    QASMReader reader(std::move(circuitFileName));
    reader.parse(gateList);
    vector<gate> rev_gateList(gateList.rbegin(),gateList.rend());
    int nqubits = reader.getVerNum();//number of logical
    
    logout << "Seed: "<<global_algorithm_seed<<endl;
    logout << "The number of logic qubits: " << nqubits << endl;
    logout << "The number of CNOT gates: " << gateList.size() << endl;
    
    //generate coupling map , distance matrix and interaction matrix
    // choose a device
    switch(device_enum_map[selected_device_enum])
    {
        case TOKYO:  cg.buildQX20();break;
        case GUADALUPE: cg.build_IBMQ_GUADALUPE();break;
        case ROCHESTER: cg.build_IBM_Rochester();break;
        case SYCAMORE:  cg.build_Google_Sycamore();break;
        case torino: cg.build_torino();break;
        default: logout<<"[ERROR] no such device"<<endl;exit(-1);
    }

    if(cg.positions<nqubits){
        logout<<"[INFO] no enough qubits"<<endl;
        exit(0);
    }
    
    int nlocations = cg.getPositions();
    // cout<<"[INFO] start computeDist"<<endl;
    cg.computeDist();
    // cout<<"[INFO] start getDist"<<endl;
    dist = cg.getDist();
    // cout<<"[INFO] start initInteractionMatrix"<<endl;

    //Initialize qubits[]   -device  qubit
    //and        locations[]-logical qubit
    vector<int> locations(nqubits,-1);
    vector<int> qubits(nlocations,-1);
    int minCnotCount = numeric_limits<int>::max();

    vector<int> qubitsArrOfBest;
    int dirOfBest = 0;
    vector<gate> phyCirOfBest;
    int noBetterCounter = 0; 
    // cout<<"[INFO] start randInitMapping"<<endl;
    randInitMapping(qubits,nlocations,nqubits);
    // cout<<"[INFO] end randInitMapping"<<endl;
    int iter_count = 0;
    bool rev_use;

    while(true)
    {
        PerurbationByShuffle(qubits,locations,nlocations);

        noBetterCounter++;
        bool isNeedNoGate = false;
        int inner_count = 0;
        int epoch_min_cnot_count=numeric_limits<int>::max();

        // weight为order
        rev_use=gateList[0].weight>gateList.back().weight;

        #ifndef ABLATION_MODIFY
        preTreat(nqubits, rev_use? rev_gateList:gateList, &cg,dist,qubits,locations);
        #endif

        while(true)
        {
            iter_count++;
            inner_count++;
            vector<int> tmpQubitsArr = qubits;
            int curCnotCount = swapAndCnotBasedMapping_advance(nqubits, rev_use? rev_gateList:gateList , &cg, phyCir, dist, qubits, locations,rev_use);
            if(curCnotCount==-1)
            {
                cur_ans_bridge_count = 0;
                phyCir.clear();
                break;
            }
            if(curCnotCount < minCnotCount)
            {
                noBetterCounter = 0;
                minCnotCount = curCnotCount;
                phyCirOfBest = phyCir;
                current_time = chrono::system_clock::now();
                auto time_cost=((double)chrono::duration_cast<chrono::milliseconds>(current_time - start_time).count())/1000;
                logout<<"[INFO] better ans: "<<minCnotCount<<" in time: "<<time_cost<<"s. in iter "<<iter_count<<endl;
                
                if(rev_use) //reverse
                {
                    dirOfBest = 1;
                    qubitsArrOfBest = qubits;
                }
                else //forward traversal
                {
                    dirOfBest = 0;
                    qubitsArrOfBest = tmpQubitsArr;
                }
                // to do bridge
                best_ans_bridge_count = cur_ans_bridge_count;
                cur_ans_bridge_count = 0;
            }
            // to do bridge
            cur_ans_bridge_count = 0;

            rev_use=!rev_use;

            phyCir.clear();
            if(minCnotCount == 0)
            {
                isNeedNoGate = true;
                break;
            }
            if(epoch_min_cnot_count>curCnotCount)
            {
                epoch_min_cnot_count=curCnotCount;
                inner_count = 0;
            }
            if(inner_count >= inLimit)
                break;
        }

        current_time = chrono::system_clock::now();
        auto time_cost=((double)chrono::duration_cast<chrono::milliseconds>(current_time - start_time).count())/1000;
        if(time_cost > cutoff_time || noBetterCounter == C_MAX || isNeedNoGate)
            break;
    }

    if(dirOfBest == 1)
        reverse(phyCirOfBest.begin(),phyCirOfBest.end());
    
    answer_checker(minCnotCount,gateList,phyCirOfBest,cg,nqubits,qubitsArrOfBest,locations);

    logout << "Cnot gate count: " << minCnotCount << endl;
    logout << "Size of the best physical circuit:" << phyCirOfBest.size() << endl;
    logout << "Direction of best mapping: " << dirOfBest << endl;//最优方向
    logout << "[INFO] restart count: "<<iter_count<<endl;
    logout << "[INFO] total bridge decision count: "<<total_bridge_count<<endl;
    #ifdef STAT
    logout << "[INFO] total decision count: "<<total_decision_count<<endl;
    logout << "[INFO] tie: "<<ties<<endl;
    logout << "[INFO] bridges num: "<<bris<<endl;
    #endif
    logout << "[INFO] best decision bridge count: "<<best_ans_bridge_count<<endl;
}

int doExperiments(vector<string>files,int seed,string log_path,int epoch=1)
{
    global_algorithm_seed=seed;

    cout<<"#File:"<<files.size()<<endl;
    auto filenameSplit=[](string &file)->pair<string,string>{
        int index;
        // -1 .qasm -> -5
        for(int i=file.size()-5;i;i--){
            if(file[i]=='/'||file[i]=='\\'){
                index=i;
                break;
            }
        }
        string path=file.substr(0,index);
        string name=file.substr(index+1);

        return {path,name};
    };

    using namespace chrono;
    for(auto &file:files){
        
        global_random_device.seed(seed);
        srand(seed);

        auto p=filenameSplit(file);
        cout<<"Benchmark circuit " << p.second << " is in mapping...."<<endl;
        // string logfile = "./" +log_path;
        string logfile = log_path + '/' + p.second + "_log.txt";
        logout.open(logfile,ios::app);
        for(int i=0;i<epoch;i++){
            start_time = system_clock::now();

            EffectiveQMFramework(file);
            current_time = system_clock::now();

            auto epoch_time=((double)duration_cast<milliseconds>(current_time - start_time).count())/1000;
            logout <<"Time cost :" << epoch_time << " seconds." << endl;
        }
        logout<<endl;
        logout<<endl<<endl;
        
        logout.close();
        cout<<"Benchmark circuit " << p.second << " has been mapped."<<endl;
    }
    return 0;
}

// commmand params:
// 1. seed
// 2. log path
// 3. device enum
// 4. IN_LIMIT
// 5. MAX_LOOK_SIZE
// 6. C_MAX
// 7 ... n. file path
int main(int argc, char* argv[])
{
    if(argc<7)
    {
        cerr<<"no enough args."<<endl;
        exit(-1);
    }

    vector<string> files;
    int seed=atoi(argv[1]);
    string log_path=argv[2];
    selected_device_enum=argv[3];

    IN_LIMIT=atoi(argv[4]);
    MAX_LOOK_SIZE=atoi(argv[5])-1;
    C_MAX=atoi(argv[6]);

    for(int i=7;i<argc;i++)
        files.emplace_back(argv[i]);
    
    doExperiments(files,seed,log_path);
    
    return 0;
}