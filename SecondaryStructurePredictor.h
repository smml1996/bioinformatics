//
// Created by Stefanie Muroya lei on 6/20/20.
//

#ifndef BIOINFORMATICS_SECONDARYSTRUCTUREPREDICTOR_H
#define BIOINFORMATICS_SECONDARYSTRUCTUREPREDICTOR_H

#include <string>
#include <vector>

using namespace std;

#define WRITE_MATRIX_PATH "../results/matrix.txt"

class SecondaryStructurePredictor {

    static void writeComplements(const string &sequence, const vector< pair<int, int> > &pairs);
    static void initMatrix(vector< vector<int> > & matrix);
    static void writeMatrix(const string &sequence, const vector< vector<int> > & result);
    static bool isComplement(const char &first, const char &second);
    static int alpha(const string &sequence,const int &i,const int &j);
    static int best(string sequence, int i, int j);
    static void backtrack(const vector< vector<int> > &dp, vector< pair<int, int> > &aligned, int i, int j, const string &sequence );
    static void getDotBracket(const vector< pair<int, int> > &aligned, const string &sequence);

    static bool isValidIndex(const vector< vector<int> > &matrix, const int &i , const int &j);

public:

    static void predictSequence(const string &sequence);

};


#endif //BIOINFORMATICS_SECONDARYSTRUCTUREPREDICTOR_H
