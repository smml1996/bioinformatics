//
// Created by Stefanie Muroya lei on 3/31/20.
//

#include "SequenceAligner.h"
#include <iostream>
#include <fstream>
#include <string>
#include <vector>

using namespace std;

AlignmentAnswer SequenceAligner::align(const string &s, const string &t, bool justOne) {

    this->as = AlignmentAnswer(getScore(s,t),getAlignment(s,t, s.size(), t.size(), justOne) );
    return as;
}

void SequenceAligner::printDpMatrix(const string &name) {
    ofstream myfile;
    myfile.open ("../results/" + name);
    for(auto & i : matrix){
        for(auto &j : i){
            myfile << j << " ";
        }
        myfile << endl;
    }
    myfile.close();
}

void SequenceAligner::printAlignments(const string &name) {
    ofstream myfile;
    myfile.open ("../results/" + name);
    for(const auto & alignment : this->as.alignment){
        myfile << alignment.first << " " << alignment.second << endl;
    }

    myfile.close();

}

vector< vector<int> > SequenceAligner::getMatrix() {
    return this->matrix;
}

void SequenceAligner::setMatrix(vector< vector<int> > &m){
    this->matrix = m;
}


