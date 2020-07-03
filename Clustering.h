//
// Created by Stefanie Muroya lei on 5/4/20.
//

#ifndef BIOINFORMATICS_CLUSTERING_H
#define BIOINFORMATICS_CLUSTERING_H

#include <string>
#include <vector>

#define MINIMO_MINIMOS 0
#define MINIMO_MAXIMOS 1
#define MINIMO_PROMEDIO 2

using namespace  std;

typedef vector< vector<double> > matrix;
class Clustering {
    static double dotProduct(const matrix &m1, const matrix &m2);
    static double getElementsSum(const matrix &m, bool isSquareElements = false);
    static double pearsonCorr(const matrix &m1, const matrix &m2);
    static pair<int, int> getIndices(matrix d, int tipo, vector< vector<int> >groups,bool isApp);
    static matrix getNewD(matrix d, int size, int tipo, int tempIndexI, int tempIndexJ);
    static void updateCofeneticMatrix(const pair<int, int> &indices, const vector< vector<int> > &groups, matrix &cofeneticMatrix, const double &val);
public:
    void run(matrix, int tipo = MINIMO_MINIMOS,const string& filePath = "../results/distances.txt", bool appendDistances = false);
    void run(const vector<string> &seqs, int tipo = MINIMO_MINIMOS,const string& filePath = "../results/distances.txt");

};


#endif //BIOINFORMATICS_CLUSTERING_H
