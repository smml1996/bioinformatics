//
// Created by Stefanie Muroya lei on 4/1/20.
//

#include "SmithWaterman.h"
#include <algorithm>
#include <iostream>

using namespace std;

int SmithWaterman::getScore(const string &s, const string &t) {
    matrix = vector< vector<int> >(s.size() +1, vector<int>(t.size()+1));

    int score = -1;
    for(int i = 0; i < matrix.size(); i++){
        for(int j = 0; j < matrix[0].size(); j++){
            if( i == 0 || j == 0){
                matrix[i][j] = 0;
            }else{
                matrix[i][j] = max(0, max(matrix[i][j-1], matrix[i-1][j]) - 2 );
                matrix[i][j] = max(matrix[i][j], matrix[i-1][j-1] + (s[i-1] == t[j-1] ? 1 : -1));
            }
            score = max(matrix[i][j], score);
        }
    }

    return score;
}

vector< pair<string, string> > SmithWaterman::getAlignment(const string &s, const string &t, const int &i,
        const int &j, bool justOne) {

    vector< pair<string, string> > ans;

    if(matrix.size() == 0 || matrix[0].size() == 0) return ans;

    int indexI = 0, indexJ = 0;

    for(int k = 0; k < matrix.size(); k++){
        for(int h = 0; h < matrix[0].size(); h++){
            if(matrix[indexI][indexJ] <= matrix[k][h]){
                indexI = k;
                indexJ = h;
            }
        }
    }
    if(matrix[indexI][indexJ] != 0){
        return alignFromIndex(s,t,indexI, indexJ);
    }

    return ans;
}

vector<pair<string, string> >
SmithWaterman::alignFromIndex(const string &s, const string &t, const int &sIndex, const int &tIndex) {


    vector< pair<string, string> > ans;
    vector< pair<string, string> > temp;
    if(sIndex < 1  || tIndex < 1 ) return vector< pair<string, string> >(1);
    if(matrix[sIndex][tIndex] < 1) return vector< pair<string, string> >(1);;

    if(
        (s[sIndex-1] == t[tIndex-1] && matrix[sIndex-1][tIndex-1] + 1 == matrix[sIndex][tIndex]) ){


        temp = alignFromIndex(s,t, sIndex-1, tIndex-1);
        for(const auto &alignment : temp){
            ans.emplace_back(alignment.first + s[sIndex-1] , alignment.second + t[tIndex-1]);
        }

    }else{
        return  vector< pair<string, string> >(1);
    }
    return ans;
}
