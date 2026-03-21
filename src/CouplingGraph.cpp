#include "CouplingGraph.h"

CouplingGraph::~CouplingGraph(){
    for(int i=0; i<positions; i++)
        delete [] dist[i];
    delete [] dist;
}

//init dist and path
void CouplingGraph::initDist(){
    dist = new int*[positions];
    path = new int*[positions];
    for(int i=0; i<positions; i++)
    {
        dist[i] = new int[positions];
        path[i] = new int[positions];
    }

    for(int i=0; i<positions; i++)
        for(int j=0; j<positions; j++)
        {
            if(i==j)
                dist[i][j]=0;
            else
                dist[i][j]=MAX_DIST;
            path[i][j] = j;
        }


    //Initialize distance matrix to a adjacency matrix
    for(int i=0;i<positions;i++){
        for(unsigned int j=0;j<adjList[i].size();j++){
            int k=adjList[i][j];
            dist[i][k]=1;
        }
    }

    computeDist();
}

void CouplingGraph::buildQX20(){
    fstream fin;
    fin.open("./device/QX20",ios::in);
    int n,m;
    fin>>n>>m;
    positions=n;
    adjList.resize(positions);
    deg.resize(positions);
    while(m--){
        int u,v;
        fin>>u>>v;
        adjList[u].push_back(v);
        adjList[v].push_back(u);
        deg[u]+=1;
        deg[v]+=1;
    }
    initDist();
    fin.close();
    cout<<"[INFO] build tokyo Over"<<endl;
}

void CouplingGraph::build_IBMQ_GUADALUPE()
{
    fstream fin;
    fin.open("./device/IBMQ_GUADALUPE",ios::in);
    int n,m;
    fin>>n>>m;
    positions=n;
    adjList.resize(positions);
    deg.resize(positions);
    while(m--){
        int u,v;
        fin>>u>>v;
        adjList[u].push_back(v);
        adjList[v].push_back(u);
        deg[u]+=1;
        deg[v]+=1;
    }
    initDist();
    fin.close();
    
    cout<<"[INFO] build IBMQ GUADALUPE Over"<<endl;
}

void CouplingGraph::build_Google_Sycamore(){
    fstream fin;
    fin.open("./device/Google_Sycamore",ios::in);
    int n,m;
    fin>>n>>m;
    positions=n;
    adjList.resize(positions);
    deg.resize(positions);
    while(m--){
        int u,v;
        fin>>u>>v;
        adjList[u].push_back(v);
        adjList[v].push_back(u);
        deg[u]+=1;
        deg[v]+=1;
    }
    initDist();
    fin.close();
    
    cout<<"[INFO] build Google Sycamore Over"<<endl;
}
void CouplingGraph::build_IBM_Rochester()
{
    fstream fin;
    fin.open("./device/IBM_Rochester",ios::in);
    int n,m;
    fin>>n>>m;
    positions=n;
    adjList.resize(positions);
    deg.resize(positions);
    while(m--){
        int u,v;
        fin>>u>>v;
        adjList[u].push_back(v);
        adjList[v].push_back(u);
        deg[u]+=1;
        deg[v]+=1;
    }
    initDist();
    fin.close();
    cout<<"[INFO] build IBM Rochester Over"<<endl;
}

void CouplingGraph::build_torino(){
    fstream fin;
    fin.open("./device/torino",ios::in);
    int n,m;
    fin>>n>>m;
    positions=n;
    adjList.resize(positions);
    deg.resize(positions);
    while(m--){
        int u,v;
        fin>>u>>v;
        adjList[u].push_back(v);
        adjList[v].push_back(u);
        deg[u]+=1;
        deg[v]+=1;
    }
    initDist();
    fin.close();
    cout<<"[INFO] build IBM Torino Over"<<endl;
}

//Floyd-Warshall Algorithm
void CouplingGraph::computeDist(){
    for(int k = 0 ; k < positions ; k ++)
    {
        for(int i = 0 ; i < positions ; i ++)
        {
            for(int j = 0 ; j < positions ; j ++)
            {
                if(dist[i][j] > dist[i][k] + dist[k][j])
                {
                    dist[i][j] = dist[i][k] + dist[k][j];
                    path[i][j] = path[i][k];
                }
            }
        }
    }
 }

vector<int> CouplingGraph::getPath(int i, int j){
    vector<int> pathBetween;
    pathBetween.push_back(i);
    int k = path[i][j];
    while(k!=j){
        pathBetween.push_back(k);
        k = path[k][j];
    }
    pathBetween.push_back(j);
    return pathBetween;
 }

//print
void CouplingGraph::showInfo(){
    for(unsigned int i=0; i<adjList.size(); i++)
    {
        cout<<i<<":";
        for(unsigned int j=0;j<adjList[i].size();j++)
            cout<<adjList[i][j]<<",";
        cout<<std::endl;
    }
 }
