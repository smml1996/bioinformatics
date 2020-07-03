    //
// Created by Stefanie Muroya lei on 6/13/20.
//

#ifndef BIOINFORMATICS_ADDITIVETREE_H
#define BIOINFORMATICS_ADDITIVETREE_H

#include <string>


using namespace std;

typedef vector< vector<double> > matrix;

class AdditiveTree {
    static void writeTree(matrix d, int countElements);
    static vector<int> updateMapping(vector<int> mapping, const int &index);
    static matrix getNewD(matrix d, const int &index);
    static array<int, 3> getTransitives(matrix d);
    static int reduce(matrix &d);
    static int getAlpha(const matrix &d);
    static bool getTree(matrix d, matrix &adj, vector<int> mapping, const int &totalElements, int &currentInternal);
public:
    static void run(matrix d,const string& filePath = "../results/distances.txt", bool appendDistances = false);
    static void run(const vector<string> &seqs, const string& filePath = "../results/distances.txt");

};


#endif //BIOINFORMATICS_ADDITIVETREE_H
