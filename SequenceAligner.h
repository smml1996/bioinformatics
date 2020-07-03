//
// Created by Stefanie Muroya lei on 3/31/20.
//

#ifndef BIOINFORMATICS_SEQUENCEALIGNER_H
#define BIOINFORMATICS_SEQUENCEALIGNER_H

#include <string>
#include <utility>
#include <vector>

using namespace std;

struct AlignmentAnswer{
    int score;
    vector< pair<string, string> > alignment;
    AlignmentAnswer(const int &score, vector< pair<string, string> > alignment): score(score), alignment(std::move(alignment)){};
    AlignmentAnswer(){
        score = 0;
    };
};

class SequenceAligner {
protected:
    vector< vector<int> > matrix;
    AlignmentAnswer as;

public:
    virtual int getScore(const string &s, const string &t) = 0;
    virtual vector< pair<string, string> > getAlignment(const string &s, const string &t, const int &i, const int &j, bool justOne = false) = 0;
    AlignmentAnswer align(const string &s, const string &t, bool justOne = false);
    void printDpMatrix(const string &name);
    void printAlignments(const string &name);
    vector< vector<int> > getMatrix();
    void setMatrix( vector< vector<int> > &m);
};


#endif //BIOINFORMATICS_SEQUENCEALIGNER_H
