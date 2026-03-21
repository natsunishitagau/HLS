#ifndef GATE_H_INCLUDED
#define GATE_H_INCLUDED

//CNOT gate
struct gate{
    int control;
    int target;
    int weight; //order in original circuit,
                //0~inf: original gate number,
                //-1:SWAP gate, 
                //-2~-inf: Bridge gate,-x means this Bridge replace x-2th original gate
    int cxt;

    gate(int c=-1,int t=-1,int w=-1){
        control = c;
        target = t;
        weight = w;
        cxt=c^t;
    }

    bool operator ==(const gate & otherGate) const{

        return weight == otherGate.weight ;
    }
    bool operator <(const gate & otherGate) const{

        return weight < otherGate.weight;
    }
};


#endif // GATE_H_INCLUDED
