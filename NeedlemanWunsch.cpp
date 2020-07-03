//
// Created by Stefanie Muroya lei on 3/31/20.
//

#include "NeedlemanWunsch.h"
#include <string>
#include <algorithm>
using namespace std;

int NeedlemanWunsch::getScore(const string &s, const string &t) {

    matrix = vector< vector<int> > (s.size() +1, vector<int>(t.size()+1));

    for(int i = 0; i < s.size() +1; i++){
        for(int j = 0; j < t.size()+1; j++){
            if(i == 0){
                if(j == 0){
                    matrix[i][j] = 0;
                }else{
                    matrix[i][j] = matrix[i][j-1] - 2;
                }
            }else if(j == 0){
                matrix[i][j] = matrix[i-1][j] - 2;
            }else{
                matrix[i][j] =  max(matrix[i-1][j], matrix[i][j-1]) -2;
                matrix[i][j] = max(matrix[i][j], matrix[i-1][j-1]  + (s[i-1]==t[j-1] ? 1 : -1));
            }
        }
    }
    return matrix[s.size()][t.size()];
}

vector< pair<string, string> >NeedlemanWunsch::getAlignment(const string &s, const string &t, const int &i, const int &j, bool justOne) {
    vector< pair<string, string> > temp;


    if(i == 0){
        if(j == 0){
            return vector<pair<string, string>>(1);
        }
        temp = getAlignment(s,t, i, j-1, justOne);

        for(auto &alignment : temp){
            alignment.first += "-";
            alignment.second += t[j-1];
        }
        return temp;
    }

    if(j == 0){
        temp = getAlignment(s,t, i-1, j, justOne);
        for(auto &alignment : temp){
            alignment.first+=s[i-1];
            alignment.second+="-";
        }

        return temp;
    }

    vector< pair<string, string> > ans;

    if((s[i-1] == t[j-1] && matrix[i-1][j-1] + 1 == matrix[i][j] )
                            || (s[i-1] != t[j-1] && matrix[i-1][j-1] - 1 == matrix[i][j])
                                                                                        ){
        temp = getAlignment(s, t, i-1, j-1, justOne);
        for(const auto &alignment : temp){
            ans.emplace_back(alignment.first + s[i-1],alignment.second + t[j-1]);
        }

        if(justOne)
            return ans;

    }

    if(matrix[i-1][j] -2 == matrix[i][j]){
        temp = getAlignment(s,t, i-1, j, justOne);

        for(const auto &alignment : temp){
            ans.emplace_back(alignment.first + s[i-1],alignment.second + "-");
        }

        if(justOne)
            return ans;

    }

    if(matrix[i][j-1] - 2 == matrix[i][j]){

        temp = getAlignment(s,t, i, j-1, justOne);

        for(const auto &alignment : temp){
            ans.emplace_back(alignment.first + "-",alignment.second + t[j-1]);
        }
        if(justOne){
            return ans;
        }


    }

    if(i == matrix.size()-1 && j == matrix[0].size()-1){
        sort(ans.begin(), ans.end(), NeedlemanWunsch::compare_alignments);
    }


    return ans;

}

int NeedlemanWunsch::get_score(const string &s) {
    int p = 0;
    int q = 0;
    for(int i = 0; i < s.size(); i++){
        if(s[i] == '-'){
            if(i == 0 || s[i-1] != '-'){
                p++;
            }else{
                q++;
            }
        }
    }
    return s.size() - p*10 - q*5;
}

bool NeedlemanWunsch::compare_alignments(const pair<string, string> &s1, const pair<string, string> &s2) {

    int score1 = min(get_score(s1.first), get_score(s1.second));
    int score2 = min(get_score(s2.first), get_score(s2.second));

    return score1 > score2;

}
