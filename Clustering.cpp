//
// Created by Stefanie Muroya lei on 5/4/20.
//

#include "Clustering.h"
#include "NeedlemanWunsch.h"
#include "ProgressiveMSA.h"
#include <fstream>
#include <cmath>
#include <iostream>

using namespace std;


void Clustering::run(const vector<string> &seqs, int tipo, const string& filePath) {


    matrix d = vector< vector<double> > (seqs.size(), vector<double>(seqs.size(), 0));

    pair<string, string>  alignment;

    double up, down;
    ofstream myfile = ofstream(filePath);
    for(int i = 0; i < seqs.size(); i++){
        for(int j = i + 1; j < seqs.size(); j++){

            NeedlemanWunsch nw;
            alignment = nw.align(seqs[i], seqs[j], true).alignment[0];


            up = 0;
            down = 0;



            for(int k = 0; k < alignment.first.size(); k++){
                if(alignment.first[k] == '-'){
                    if(alignment.second[k] != '-'){
                        up++;
                    }
                }else if(alignment.second[k] == '-'){
                    up++;
                }else{
                    down++;
                }
            }


            d[i][j] = trunc(up/down * 10000)/100;
            myfile << d[i][j] << "\t";
        }

        myfile << endl;
    }

    myfile.close();
    return this->run(d, tipo,filePath, true);
}




pair<int, int> Clustering::getIndices(matrix d, int tipo, vector< vector<int> >groups, bool isApp) {


    int indexI = d[0][0], indexJ = d[0][1];
    for(int i = 0; i < d.size(); i++){
        for(int j = i+1; j < d[0].size(); j++){

                if(d[i][j] < d[indexI][indexJ] || (i == 0 && j == 1)){
                    indexI = i;
                    indexJ = j;
                }

        }
    }
    ofstream myfile;
    if(!isApp){
        myfile = ofstream("../results/groups_trace.txt");
    }else{
        myfile = ofstream("../results/groups_trace.txt",ios::out | ios::app);
    }


    myfile << "(";
    for(int i = 0; i < groups[indexI].size(); i++){
        myfile << (char)(groups[indexI][i] + int('A')) << "\t";
    }

    myfile << ") \t\t y \t\t (";

    for(int i = 0; i < groups[indexJ].size(); i++){
        myfile << (char)(groups[indexJ][i] + int('A'))  << "\t";
    }

    myfile << ")" << endl;

    myfile.close();

    return make_pair(indexI, indexJ);
}

void Clustering::run(matrix d , int tipo, const string &filePath, bool appendDistances){


        matrix cofeneticMatrix = vector< vector<double> >(d.size(), vector<double>(d.size(),-1));
        ProgressiveMSA::normalizeD(d);

        const matrix initialMatrix = d;

        int size = d.size();

        bool isAppend = false;

        vector< vector<int> > groups = ProgressiveMSA::initializeGroups(size);

        pair<int, int>indices;


        while(size > 2){

            indices = getIndices(d, tipo, groups, isAppend);
            updateCofeneticMatrix(indices, groups, cofeneticMatrix, d[indices.first][indices.second]);
            groups = ProgressiveMSA::joinIndices(groups, indices.first, indices.second);

            ProgressiveMSA::printGroups(groups, isAppend);

            size--;
            d = getNewD(d, size,tipo,  indices.first, indices.second);

            ProgressiveMSA::printD(d, "../results/distances.txt", appendDistances);

            isAppend = true;
            appendDistances = true;
        }

    ofstream myfile;

        myfile = ofstream("../results/groups_trace.txt",ios::out | ios::app);



    myfile << "(";
    for(int i = 0; i < groups[0].size(); i++){
        myfile << (char)(groups[0][i] + int('A')) << "\t";
    }

    myfile << ") \t\t y \t\t (";

    for(int i = 0; i < groups[1].size(); i++){
        myfile << (char)(groups[1][i] + int('A'))  << "\t";
    }

    myfile << ")" << endl;

    myfile.close();

    updateCofeneticMatrix(make_pair(0,1), groups, cofeneticMatrix, d[0][1]);

    ProgressiveMSA::printD(cofeneticMatrix, "../results/matrizCofenetica.txt", false);

    cout << "CCC: " << pearsonCorr(initialMatrix, cofeneticMatrix) << endl;
}

matrix Clustering::getNewD(matrix d, int size, int tipo, int tempIndexI, int tempIndexJ) {
    matrix ans = vector<vector<double>>(size,vector<double>(size));

    int indexIA = 0;
    int indexJA;
    for(int i = 0; i < d.size(); i++){

        if(tempIndexJ == i) continue;
        indexJA = indexIA + 1;
        for(int j = i+1; j < d.size(); j++){
            if(tempIndexJ == j) continue;
            if(tempIndexI == i){
                if(tipo == MINIMO_MINIMOS){
                    ans[indexIA][indexJA] = min(d[tempIndexI][j] , d[tempIndexJ][j]);
                }else if(tipo == MINIMO_MAXIMOS){
                    ans[indexIA][indexJA] = max(d[tempIndexI][j] , d[tempIndexJ][j]);
                }else{
                    ans[indexIA][indexJA] = trunc((d[tempIndexI][j] + d[tempIndexJ][j])/2.0 * 100)/100.0;
                }

            }else if(tempIndexI == j){
                if(tipo == MINIMO_MINIMOS){
                    ans[indexIA][indexJA] = min(d[tempIndexI][i] , d[tempIndexJ][i]);
                }else if(tipo == MINIMO_MAXIMOS){
                    ans[indexIA][indexJA] = max(d[tempIndexI][i] , d[tempIndexJ][i]);
                }else{
                    ans[indexIA][indexJA] = trunc((d[tempIndexI][i] + d[tempIndexJ][i])/2.0 * 100)/100.0;
                }
            }else{
                ans[indexIA][indexJA] = d[i][j];
            }

            indexJA++;

        }
        indexIA++;
    }
    ProgressiveMSA::normalizeD(ans);
    return ans;
}

void Clustering::updateCofeneticMatrix(const pair<int,int> &indices, const vector<vector<int>> &groups,
                                       matrix &cofeneticMatrix, const double &val) {


    for(auto const & v1 : groups[indices.first]){
        for(auto const &v2 : groups[indices.second]){
            if(cofeneticMatrix[v1][v2] == -1){
                cofeneticMatrix[v1][v2] = val;
                cofeneticMatrix[v2][v1] = cofeneticMatrix[v1][v2];
            }
        }
    }

}

double Clustering::pearsonCorr(const matrix &m1, const matrix &m2) {
    double N = ((m1.size() -1.0) * (m1.size()))/2.0;

    double elementsSum1 = getElementsSum(m1);
    double elementsSum2 = getElementsSum(m2);
    double elementsSumSquared1 = getElementsSum(m1, true);
    double elementsSumSquared2 = getElementsSum(m2, true);
    double matricesProduct = dotProduct(m1, m2);

    double numerador = N*matricesProduct - (elementsSum1*elementsSum2);

    double denominator = sqrt(N*elementsSumSquared1 - elementsSum1*elementsSum1) *
            sqrt(N*elementsSumSquared2 - elementsSum2*elementsSum2);


    return numerador/denominator;



}

double Clustering::getElementsSum(const matrix &m, bool isSquareElements) {
    double ans = 0.0;

    for(int i = 0; i < m.size(); i++){

        for(int j = i+1; j < m[0].size(); j++){

            ans+= isSquareElements ? m[i][j] * m[i][j] : m[i][j];
        }

    }

    return ans;
}

double Clustering::dotProduct(const matrix &m1, const matrix &m2) {
    double ans = 0;

    for(int i = 0; i < m1.size(); i++){
        for(int j = i+1; j < m1[0].size(); j++){
            ans+= (m1[i][j] * m2[i][j]);
        }
    }

    return ans;
}




