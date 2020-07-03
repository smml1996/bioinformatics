//
// Created by Stefanie Muroya lei on 4/15/20.
//

#ifndef BIOINFORMATICS_MSA_H
#define BIOINFORMATICS_MSA_H

#include "NeedlemanWunsch.h"
#include <vector>
#include <string>

using namespace std;



class MSA {
private:
    NeedlemanWunsch nw;
    vector< vector<int> >m;
    vector<int> scoresSum;
    int center;
    void saveMatrix();
    void saveAlignments(const vector<string> &ans);
    void findCenter();
    void calculateMatrix(const vector<string> &seqs);
    static void fillGaps(vector<string> &ans, const int& maxSize);
public:
    vector<string>  getAlignments(const vector<string> &seqs);

};


#endif //BIOINFORMATICS_MSA_H
