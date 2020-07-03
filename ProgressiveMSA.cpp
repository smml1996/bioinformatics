//
// Created by Stefanie Muroya lei on 4/23/20.
//

#include "ProgressiveMSA.h"
#include "NeedlemanWunsch.h"
#include <fstream>
#include <cmath>
#include <iostream>

using namespace  std;

void ProgressiveMSA::run(matrix d, const string& filePath) {

    normalizeD(d);


    int size = d.size();


    vector<double> rowSums;
    matrix Q;

    bool isAppend = false;

    vector< vector<int> > groups = initializeGroups(size);

    double minimum;
    int tempIndexI, tempIndexJ;
    while(size > 2){

        rowSums = getRowSums(d);
        Q = calculateQ(d, size, rowSums,minimum, tempIndexI, tempIndexJ);
        printQ(Q, filePath, isAppend);

        groups = joinIndices(groups, tempIndexI, tempIndexJ);

        printGroups(groups, isAppend);

        size--;
        d = getNewD(d, size, tempIndexI, tempIndexJ);

        printD(d, "../results/distances.txt", isAppend);
        isAppend = true;
    }

}

vector<double> ProgressiveMSA::getRowSums(const matrix &d) {

    vector<double> ans (d.size(), 0);

    for(int i = 0; i < d.size(); i++){
        for(int j = 0; j < d[0].size(); j++){
            if(i != j)
                ans[i]+= d[i][j];
        }
    }

    return ans;
}

double ProgressiveMSA::getQij(const matrix &d, const int &n, const int &i, const int &j, const vector<double> &rowSums) {
    return trunc((d[i][j] - (rowSums[i] + rowSums[j])/(double(n-2))) * 100)/100.0;
}

matrix ProgressiveMSA::calculateQ(const matrix &d, const int &size, const vector<double> &rowSums, double &minimum, int &tempIndexI, int &tempIndexJ) {

    matrix Q = vector< vector<double> > (size, vector<double>(size));

    for(int i = 0; i < size; i++){
        for(int j = i+1; j < size; j++){

            Q[i][j] = getQij(d, size, i,j, rowSums);
            if(i == 0 && j == i+1){
                minimum = Q[i][j];
                tempIndexI = i;
                tempIndexJ = j;
            }else if(minimum > Q[i][j]){
                minimum = Q[i][j];
                tempIndexI = i;
                tempIndexJ = j;
            }

        }
    }

    return Q;

}

void ProgressiveMSA::printQ(const matrix &Q, const string &filePath, const bool &isApp) {

    ofstream myfile;
    if(!isApp){
        myfile = ofstream(filePath);
    }else{
        myfile = ofstream(filePath,ios::out | ios::app);
    }

    myfile << "Q: " << endl;
    for(int i = 0; i < Q.size(); i++){
        for(int j = 0; j < Q[0].size(); j++){
            if( j <= i){
                myfile << "-\t";
            }else{
                myfile << Q[i][j] << "\t";
            }
        }

        myfile << endl;
    }
    myfile << "**&*&*&*&*&*&*&*&*&*&**&*&*&*&*&*&*&*&*&&*&*&*&**&*&&**&" << endl << endl;

    myfile.close();
}

vector<vector<int> > ProgressiveMSA::initializeGroups(const int &size) {
    vector< vector<int> > ans(size, vector<int>());

    for(int i = 0; i < size; i++){
        ans[i].push_back(i);
    }

    return ans;
}

vector<vector<int> > ProgressiveMSA::joinIndices(vector<vector<int>> &groups, int tempIndexI, int tempIndexJ) {

    vector< vector<int> > ans(groups.size() -1);
    if(tempIndexI > tempIndexJ){
        swap(tempIndexJ, tempIndexI);
    }

    int localIndex = 0;
    for(int i =0; i < groups.size(); i++){
        if(i == tempIndexJ) continue;
        for(int j = 0; j < groups[i].size(); j++){
            ans[localIndex].push_back(groups[i][j]);
        }
        localIndex++;
    }

    for(auto n : groups[tempIndexJ]){
        ans[tempIndexI].push_back(n);
    }

    return ans;
}

void ProgressiveMSA::printGroups(const vector<vector<int>> &groups, const bool &isApp, const string& filePath) {
        ofstream myfile;

        if(!isApp){
            myfile = ofstream(filePath);
        }else{
            myfile = ofstream(filePath,ios::out | ios::app);
        }

        for(int i = 0; i < groups.size(); i++){
            myfile << i <<": ";

            for(auto n : groups[i]){
                myfile << n << " ";
            }

            myfile << endl;
        }

        myfile << "*&*&*&**&**&*&*&**&*&**&*&*&*&&*&***&*&*&*&**&*&*&&*&*&*&" << endl;

        myfile << endl;

}

matrix ProgressiveMSA::getNewD(const matrix &d, const int &newSize, int tempIndexI, int tempIndexJ) {
    matrix ans = vector< vector<double> >(newSize, vector<double>(newSize));

    double dfg = d[tempIndexJ][tempIndexI];

    int indexIA = 0;
    int indexJA;
    for(int i = 0; i < d.size(); i++){

        if(tempIndexJ == i) continue;
        indexJA = indexIA + 1;
        for(int j = i+1; j < d.size(); j++){
            if(tempIndexJ == j) continue;
            if(tempIndexI == i){

                //cout << "ponis: " <<endl;
                //cout << indexIA << " " << indexJA << endl;
                //cout << d[tempIndexI][j] << " " << d[tempIndexJ][j] << " " << dfg << endl;

                ans[indexIA][indexJA] = (d[tempIndexI][j] + d[tempIndexJ][j] - dfg)/2.0;

            }else if(tempIndexI == j){
                ans[indexIA][indexJA] = (d[tempIndexI][i] + d[tempIndexJ][i] - dfg)/2.0;
            }else{
                ans[indexIA][indexJA] = d[i][j];
            }

            indexJA++;

        }
        indexIA++;
    }
    normalizeD(ans);
    return ans;
}

void ProgressiveMSA::normalizeD(matrix &d) {

    for(int i = 0; i < d.size(); i++){
        for(int j = i+1; j < d[0].size(); j++){
            d[i][j] = trunc(d[i][j] * 100)/100.0;
            d[j][i] = d[i][j];
        }
    }

}

void ProgressiveMSA::printD(const matrix &d, const string &filePath, const bool &isApp) {
    ofstream myfile;
    if(!isApp){
        myfile = ofstream(filePath);
    }else{
        myfile = ofstream(filePath,ios::out | ios::app);
    }

    myfile << "D: " << endl;
    for(int i = 0; i < d.size(); i++){
        for(int j = 0; j < d[0].size(); j++){
            if( j <= i){
                myfile << "-\t";
            }else{
                myfile << d[i][j] << "\t";
            }
        }

        myfile << endl;
    }
    myfile << "**&*&*&*&*&*&*&*&*&*&**&*&*&*&*&*&*&*&*&&*&*&*&**&*&&**&" << endl << endl;

    myfile.close();
}

void ProgressiveMSA::run(const vector<string> &seqs, const string& filePath) {

    matrix d = vector< vector<double> > (seqs.size(), vector<double>(seqs.size(), 0));

    pair<string, string>  alignment;

    double up, down;
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
            cout << d[i][j] << "\t";
        }

        cout << endl;
    }

    return this->run(d, filePath);
}




