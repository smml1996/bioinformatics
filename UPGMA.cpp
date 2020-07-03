//
// Created by Stefanie Muroya lei on 6/1/20.
//

#include "UPGMA.h"

#include <fstream>
#include "NeedlemanWunsch.h"
#include "ProgressiveMSA.h"
#include <cmath>
#include <stack>
#include <iostream>

using namespace std;


void UPGMA::run(const vector<string> &seqs, const string &filePath) {

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
    return run(d, filePath, true);

}

void UPGMA::run(matrix d, const string &filePath, bool appendDistances) {

    ProgressiveMSA::normalizeD(d);

    stack< int > nodes;

    const auto initialMatrix = d;

    int size = d.size();

    bool isAppend = false;

    vector< vector<int> > groups = ProgressiveMSA::initializeGroups(size);

    pair<int, int>indices;


    if(!appendDistances){
        ProgressiveMSA::printD(d, "../results/distances.txt", appendDistances);
        appendDistances = true;
    }


    vector<int> ancestors(d.size());

    for(int i = 0; i < ancestors.size(); i++){
        ancestors[i] = i;
    }



    vector< vector<double> > unions;

    int countAnc = size;

    while(size > 2){

        indices = getIndices(d, groups, isAppend);



        size--;

        unions.emplace_back(5);
        unions[unions.size()-1][0] = ancestors[groups[indices.first][0]];
        unions[unions.size()-1][1] = ancestors[groups[indices.second][0]];
        unions[unions.size()-1][2] = d[indices.first][indices.second]/2.0;
        if(unions.size() == 1){
            unions[unions.size()-1][3] = d[indices.first][indices.second]/2.0;
            unions[unions.size()-1][4] = d[indices.first][indices.second]/2.0;
        }else{
            unions[unions.size()-1][3] = d[indices.first][indices.second]/2.0;

            if(ancestors[groups[indices.first][0]] >= initialMatrix.size()){

                unions[unions.size()-1][3]-= unions[ancestors[groups[indices.first][0]] -  initialMatrix.size()][2];
            }
            unions[unions.size()-1][4] = d[indices.first][indices.second]/2.0;

            if(ancestors[groups[indices.second][0]] >= initialMatrix.size()){
                unions[unions.size()-1][4]-= unions[ancestors[groups[indices.second][0]] -  initialMatrix.size()][2];
            }
        }


        refreshParent(indices.first, ancestors,groups, countAnc);
        refreshParent(indices.second, ancestors,groups, countAnc);




        countAnc++;



        d = getNewD(d, size, groups, initialMatrix, indices.first, indices.second);


        groups = ProgressiveMSA::joinIndices(groups, indices.first, indices.second);

        ProgressiveMSA::printGroups(groups, isAppend);
        ProgressiveMSA::printD(d, "../results/distances.txt", appendDistances);

        isAppend = true;
        appendDistances = true;
    }


    indices.first = 0;
    indices.second = 1;
    unions.emplace_back(5);
    unions[unions.size()-1][0] = ancestors[groups[indices.first][0]];
    unions[unions.size()-1][1] = ancestors[groups[indices.second][0]];
    unions[unions.size()-1][2] = d[indices.first][indices.second]/2.0;
    if(unions.size() == 1){
        unions[unions.size()-1][3] = d[indices.first][indices.second]/2.0;
        unions[unions.size()-1][4] = d[indices.first][indices.second]/2.0;
    }else{
        unions[unions.size()-1][3] = d[indices.first][indices.second]/2.0;

        if(ancestors[groups[indices.first][0]] >= initialMatrix.size()){

            unions[unions.size()-1][3]-= unions[ancestors[groups[indices.first][0]] -  initialMatrix.size()][2];
        }
        unions[unions.size()-1][4] = d[indices.first][indices.second]/2.0;

        if(ancestors[groups[indices.second][0]] >= initialMatrix.size()){
            unions[unions.size()-1][4]-= unions[ancestors[groups[indices.second][0]] -  initialMatrix.size()][2];
        }
    }

    printGraph(unions, initialMatrix.size());
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

}

pair<int, int> UPGMA::getIndices(matrix d, vector< vector<int> > groups, bool isAppend) {
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
    if(!isAppend){
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

matrix UPGMA::getNewD(matrix d, int size, vector< vector<int> > groups, matrix initialD, int tempIndexI, int tempIndexJ) {
    matrix ans = vector< vector<double> >(size, vector<double>(size));


    int indexIA = 0;
    int indexJA;

    double total;
    for(int i = 0; i < d.size(); i++){

        if(tempIndexJ == i) continue;
        indexJA = indexIA + 1;
        for(int j = i+1; j < d.size(); j++){
            if(tempIndexJ == j) continue;
            if(tempIndexI == i){
                total = (groups[tempIndexI].size()+groups[tempIndexJ].size()) * groups[j].size();

                ans[indexIA][indexJA] = getGroupsSum(groups, initialD, tempIndexJ,tempIndexI, j)/total;
            }else if(tempIndexI == j){

                total = (groups[tempIndexI].size()+groups[tempIndexJ].size()) * groups[i].size();
                ans[indexIA][indexJA] = getGroupsSum(groups, initialD, tempIndexJ,tempIndexI, i)/total;
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

double UPGMA::getGroupsSum(vector<vector<int>> groups, matrix d, int index1, int index2, int index3) {

    double ans = 0;


    for(int i = 0; i < groups[index1].size(); i++){
        for(int j = 0; j < groups[index3].size(); j++){
            ans+= d[groups[index1][i]][groups[index3][j]];

        }
    }

    for(int i = 0; i < groups[index2].size(); i++){
        for(int j = 0; j < groups[index3].size(); j++){
            ans+= d[groups[index2][i]][groups[index3][j]];

        }
    }


    return ans;

}

bool UPGMA::isInGroup(int index, vector<int> g) {

    for(auto a : g){
        if(a == index)return true;
    }

    return false;

}

vector<int> UPGMA::getGroup(int index, vector<vector<int>> groups) {
    for(int g = 0; g < groups.size(); g++){
        if(isInGroup(index, groups[g])){
            return groups[g];
        }
    }

    return vector<int>();
}

void UPGMA::refreshParent(int index, vector<int > &ancestors, vector< vector<int> > groups, int count) {


    for(int i = 0; i < groups[index].size(); i++){
        ancestors[groups[index][i]] = count;
    }

}

void UPGMA::printGraph(vector< vector<double> >unions, int nodes) {

    ofstream file("../results/graph.txt");

    file << nodes << endl;

    for(int i = 0; i < unions.size(); i++){
        file << "N" << i << " ";

        if(unions[i][0] < nodes){
            file << (char)(unions[i][0]  + 'A') << " ";
        }else{
            file << "N" << (unions[i][0] - nodes)<< " ";
        }

        file << trunc(unions[i][3]*100)/100.0 << endl;

        file << "N" << i << " ";

        if(unions[i][1] < nodes){
            file << (char)(unions[i][1]  + 'A') << " ";
        }else{
            file << "N" << (unions[i][1] - nodes)<< " ";
        }

        file << trunc(unions[i][4]*100)/100.0 << endl;
    }


    file.close();

    system("python3.6 ../generate_graph.py");
}
