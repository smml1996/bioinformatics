//
// Created by Stefanie Muroya lei on 4/23/20.
//

#ifndef BIOINFORMATICS_PROGRESSIVEMSA_H
#define BIOINFORMATICS_PROGRESSIVEMSA_H
#include <vector>
#include<string>

using namespace std;

typedef vector< vector<double> > matrix;

class ProgressiveMSA {
    public:
        static matrix getNewD(const matrix &d, const int &newSize, int tempIndexI, int tempIndexJ);
        static void printD(const matrix &d, const string &filePath, const bool &isApp);
        static void normalizeD(matrix &d);
        static void printGroups(const vector< vector<int> > &groups, const bool &isApp, const string& filePath = "../results/groups.txt" );
        static void printQ(const matrix &Q, const string &filePath, const bool &isApp);
        static vector<double> getRowSums(const matrix &d);
        static double getQij(const matrix &d, const int &n, const int &i, const int&j, const vector<double> &rowSums);
        matrix calculateQ(const matrix &d, const int &size, const vector<double> &rowSums, double &minimum, int &tempIndexI, int &tempIndexJ);
        static vector< vector <int> > initializeGroups(const int&size);
        static vector< vector<int> > joinIndices(vector< vector<int> > &groups, int tempIndexI, int tempIndexJ);
        void run(matrix d, const string& filePath = "../results/q.txt");
        void run(const vector<string> &seqs, const string& filePath = "../results/q.txt");
};


#endif //BIOINFORMATICS_PROGRESSIVEMSA_H
