//
// Created by Stefanie Muroya lei on 5/11/20.
//

#ifndef BIOINFORMATICS_DISSOCIATIVECLUSTERING_H
#define BIOINFORMATICS_DISSOCIATIVECLUSTERING_H
#include "ProgressiveMSA.h"
#include <vector>

class DissociativeClustering {

    double getLess(matrix distances, int t);
    void printClusters(vector< vector<int> > &clusters, bool isAppend = false);
public:
    void run(matrix distances,const string& filePath = "../results/distances.txt", bool appendDistances = false);
    void run(const vector<string> &seqs,const string& filePath = "../results/distances.txt");

};


#endif //BIOINFORMATICS_DISSOCIATIVECLUSTERING_H
