//
// Created by Stefanie Muroya lei on 3/31/20.
//

#ifndef BIOINFORMATICS_NEEDLEMANWUNSCH_H
#define BIOINFORMATICS_NEEDLEMANWUNSCH_H


#include "SequenceAligner.h"

class NeedlemanWunsch : public SequenceAligner{
    // global alignment dp solution
private:
    static int get_score(const string &s);
    static bool compare_alignments(const pair<string, string> &s1, const pair<string, string> &s2 );

public:
    int getScore(const string &s, const string &t) override;
    vector< pair<string, string> > getAlignment(const string &s, const string &t, const int &i, const int &j, bool justOne = false) override;

};


#endif //BIOINFORMATICS_NEEDLEMANWUNSCH_H
