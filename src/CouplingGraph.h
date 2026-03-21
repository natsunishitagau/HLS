#ifndef COUPLINGGRAPH_H_INCLUDED
#define COUPLINGGRAPH_H_INCLUDED

#include <iostream>
#include <fstream>
#include <vector>
#include <stdlib.h>
#include <numeric>
#include <algorithm>

using namespace std;

#define MAX_DIST 1000

class CouplingGraph{
    private:
        void initDist();

    public:
        vector<vector<int>> adjList;
        vector<int> deg;
        int positions;  //number of physical qubits
        int ** dist; //distance matrix
        int ** path; //path matrix
    public:
        CouplingGraph() {
            positions = 0;
        }
        CouplingGraph(int pos) {
            positions = pos;
        }
        ~CouplingGraph();

        int ** getDist() {return dist;}

        int getPositions() const {return positions;}

        vector<int> getDeg() {return deg;}

        //build coupling graph for IBM Q Tokyo
        void buildQX20();

        void build_IBMQ_GUADALUPE();

        void build_IBM_Rochester();

        void build_Google_Sycamore();

        void build_torino();

        //get the distance matrix with Floyd-Warshall algorithm
        void computeDist();

        vector<int> getPath(int i,int j);

        void showInfo();

};

#endif
