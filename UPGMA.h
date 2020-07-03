//
// Created by Stefanie Muroya lei on 6/1/20.
//

#ifndef BIOINFORMATICS_UPGMA_H
#define BIOINFORMATICS_UPGMA_H

#include <vector>
#include <string>
using namespace  std;

typedef vector< vector<double> > matrix;

class UPGMA {
    static void refreshParent(int index, vector<int > &ancestors, vector<vector<int>> groups, int count);
    static bool isInGroup(int index, vector<int> g);
    static void printGraph(vector< vector<double> > unions, int nodes);
    static vector<int> getGroup(int index, vector<vector<int> > groups);
    static double getGroupsSum(vector< vector<int> > groups, matrix d, int index1, int index2, int index3);
    static matrix getNewD(matrix d, int size, vector< vector<int> > groups, matrix initialD, int tempIndexI, int tempIndexJ);
    static pair<int, int> getIndices(matrix d, vector< vector<int> > groups, bool isAppend);
public:
    static void run(matrix d,const string& filePath = "../results/distances.txt", bool appendDistances = false);
    static void run(const vector<string> &seqs, const string& filePath = "../results/distances.txt");
};


#endif //BIOINFORMATICS_UPGMA_H
