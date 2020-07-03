//
// Created by Stefanie Muroya lei on 4/1/20.
//

#ifndef BIOINFORMATICS_SMITHWATERMAN_H
#define BIOINFORMATICS_SMITHWATERMAN_H

#include <string>
#include "SequenceAligner.h"

using namespace std;

class SmithWaterman : public SequenceAligner {
    //local alignment
private:
    vector< pair<string, string> > alignFromIndex(const string &s, const string &t, const int &i, const int &j);

public:
    int getScore(const string &s, const string &t) override;
    vector< pair<string, string> > getAlignment(const string &s, const string &t, const int &i, const int &j, bool justOne = false) override;
};


#endif //BIOINFORMATICS_SMITHWATERMAN_H
