//
// Created by Stefanie Muroya lei on 6/13/20.
//

#include "AdditiveTree.h"
#include "NeedlemanWunsch.h"
#include "ProgressiveMSA.h"
#include <vector>
#include <fstream>
#include <cmath>
#include <array>
#include <iostream>

using namespace  std;

void AdditiveTree::run(const vector<string> &seqs, const string &filePath) {

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

void AdditiveTree::run(matrix d, const string &filePath, bool appendDistances) {
    ProgressiveMSA::normalizeD(d);
    vector< vector<double> > adj(d.size()*2 -2, vector<double>(d.size()*2-2, -1));

    vector<int> mapping(d.size());

    for(int i = 0; i < d.size(); i++){
        mapping[i] = i;
    }

    int currentInternal = d.size();
    ProgressiveMSA::printD(d, "../results/distances.txt", false);
    if(getTree(d, adj,mapping, d.size(), currentInternal)){

        writeTree(adj, d.size());
        return;
    }

    cout << "Matriz no es transitiva" << endl;


}

int AdditiveTree::reduce(matrix &d) {

    int alpha = getAlpha(d);

    for(int i = 0; i < d.size(); i++){
        for(int j = 0; j < d.size(); j++){
            if(i!=j){
                d[i][j] = d[i][j] - alpha*2;
            }
        }
    }

    return alpha;

}

int AdditiveTree::getAlpha(const matrix &d) {

    int alpha;

    double temp = INT_MAX;
    for(int i = 0; i < d.size(); i++){
        for(int j = i+1; j < d.size(); j++){

            temp = min(d[i][j], temp);
        }
    }

    alpha = temp/2;

    if(temp - alpha*2 <= 0){
        alpha--;
    }

    return alpha;

}

bool AdditiveTree::getTree(matrix d, matrix &adj, vector<int> mapping, const int &totalElements, int &currentInternal) {

    if(d.size() < 2) return false;
    if(d.size() == 2){

        adj[mapping[0]][mapping[1]] = d[0][1];
        adj[mapping[1]][mapping[0]] = d[0][1];
        return true;
    }

    int alpha = reduce(d);
    //cout << alpha << endl;
    ProgressiveMSA::printD(d, "../results/distances.txt", true);
    array<int, 3> transitiveIndices = getTransitives(d);

    cout << transitiveIndices[0] << " " << transitiveIndices[1] << " " << transitiveIndices[2] << endl;
    if(transitiveIndices[0] == -1) return false; // 0: i, 1: j, 2:k

    vector<int> newMapping = updateMapping(mapping, transitiveIndices[1]);
    matrix newD = getNewD(d, transitiveIndices[1]);

    ProgressiveMSA::printD(newD, "../results/distances.txt", true);
    bool isPossible = getTree(newD, adj, newMapping, totalElements, currentInternal);

    if(isPossible){

        adj[mapping[transitiveIndices[1]]][currentInternal] = 0;
        adj[currentInternal][mapping[transitiveIndices[1]]] = 0;

        int closestNode = d[transitiveIndices[0]][transitiveIndices[1]] < d[transitiveIndices[1]][transitiveIndices[2]]?
                transitiveIndices[0] : transitiveIndices[2];

        int endPoint = -1;

        for(int i = 0; i < adj.size(); i++ ){
            if( adj[mapping[closestNode]][i]!=-1){
                if(endPoint == -1 || adj[mapping[closestNode]][i] < adj[mapping[closestNode]][endPoint]){
                    endPoint = i;
                }
            }
        }

        int temp = adj[mapping[closestNode]][endPoint]- d[closestNode][transitiveIndices[1]];


        adj[mapping[closestNode]][endPoint] = -1;
        adj[endPoint][mapping[closestNode]] = -1;

        adj[mapping[closestNode]][currentInternal] = d[closestNode][transitiveIndices[1]];
        adj[currentInternal][mapping[closestNode]] = d[closestNode][transitiveIndices[1]];

        adj[endPoint][currentInternal] = temp;
        adj[currentInternal][endPoint] = temp;


        for(int i = 0; i < totalElements; i++){
            for(int j = 0; j < adj.size(); j++){
                if(i!=j && adj[i][j]!=-1){

                    adj[i][j]+=alpha;
                    adj[j][i]+=alpha;
                }
            }
        }


        currentInternal++;

        return true;
    }
    return false;
}

array<int, 3> AdditiveTree::getTransitives(matrix d) {

    array<int, 3> indices;
    indices.fill(-1);
    for(int i = 0; i < d.size(); i++){
        for(int j = 0; j < d.size(); j++){
            cout << d[i][j] << " ";
            if(i!=j){
                for(int k = 0; k < d.size(); k++){
                    if(k!=j && k!=i){
                        if(d[i][j] + d[j][k] == d[i][k]){
                            indices[0] = i;
                            indices[1] = j;
                            indices[2] = k;
                            return indices;
                        }
                    }
                }
            }
        }
        cout << endl;
    }
    cout << endl;

    return indices;
}

vector<int> AdditiveTree::updateMapping(vector<int> mapping, const int &index) {

    vector<int> ans(mapping.size()-1);
    for(int i = 0, indexMapping = 0; i < mapping.size(); i++){
        if(i != index){
            ans[indexMapping] = mapping[i];
            indexMapping++;
        }
    }

    return ans;

}

matrix AdditiveTree::getNewD(matrix d, const int &index) {

    matrix newD(d.size()-1, vector<double>(d.size()-1));

    int indexI = 0, indexJ;


    for(int i = 0; i < d.size(); i++){
        if(i!=index) {
            indexJ = 0;
            for (int j = 0; j < d.size(); j++) {
                if(j!=index){

                    newD[indexI][indexJ] = d[i][j];
                    indexJ++;
                }
            }
            indexI++;
        }
    }





    return newD;
}

void AdditiveTree::writeTree(matrix d, int countElements) {

    ofstream file("../results/tree_matrix.txt");

    file << countElements << endl;

    for(int i = 0; i < d.size(); i++){
        for(int j = 0; j < d.size(); j++){
            file << d[i][j] <<" ";
        }
        file << endl;
    }

    file.close();

    system("python3.6 ../draw_unrooted.py");

}
