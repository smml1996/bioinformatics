//
// Created by Stefanie Muroya lei on 6/20/20.
//

#include "SecondaryStructurePredictor.h"

#include <fstream>
#include <climits>
#include <iostream>
using namespace std;

void SecondaryStructurePredictor::initMatrix(vector<vector<int>> &matrix) {

    int n = matrix.size();

    for(int i = 0; i < n; i++){
        matrix[i][i] = 0;
        matrix[i][i-1] = 0;
    }

}

void SecondaryStructurePredictor::writeMatrix(const string &sequence,const vector<vector<int>> &result) {

    ofstream file(WRITE_MATRIX_PATH);
    file << sequence << endl;
    for(const auto &r : result){
        for(const auto &c : r){
            file << c << "\t";
        }
        file << endl;
    }

}

void SecondaryStructurePredictor::predictSequence(const string &sequence) {

    vector< vector<int> > matrix(sequence.size(), vector<int>(sequence.size(), 0));

    initMatrix(matrix);

    const int l = sequence.size();

    int j;

    int temp1, temp2;
    for(int d = 0; d < l ; d++){
        for(int i = 0; i + d < l ; i++){
            j = i + d ;

            //cout << i << " " << j << ": " << alpha(sequence, i, j) << endl;

            /*if(isValidIndex(matrix, i+1, j)){
                matrix[i][j] = max(matrix[i][j], matrix[i+1][j]);
            }*/

            if(isValidIndex(matrix, i, j-1)){
                matrix[i][j] = max(matrix[i][j], matrix[i][j-1]);
            }

            /*if(isValidIndex(matrix, i+1, j-1)){
                matrix[i][j] = max(matrix[i+1][j-1] + alpha(sequence, i, j), matrix[i][j]);
            }*/

            for(int k = i; k < j; k++){

                    if(isComplement(sequence[k], sequence[j])) {
                        if(!isValidIndex(matrix, i, k-1)){
                            temp1 = 0;
                        }else{
                            temp1 = matrix[i][k-1];
                        }

                        if(!isValidIndex(matrix, k+1, j-1)){
                            temp2 = 0;
                        }else{
                            temp2 = matrix[k+1][j-1];
                        }

                        matrix[i][j] = max(matrix[i][j], temp1 + temp2 + 1);
                    }

            }

        }
    }

    writeMatrix(sequence, matrix);

    vector< pair<int, int> > aligned;
    backtrack(matrix, aligned, 0,l-1 ,sequence);
    getDotBracket(aligned, sequence);
    writeComplements(sequence, aligned);
}

bool SecondaryStructurePredictor::isValidIndex(const vector<vector<int>> &matrix, const int &i, const int &j) {
    if(i > -1 && j > -1){
        if(i < matrix.size() && j < matrix[0].size()){
            return true;
        }
    }

    return false;
}

int SecondaryStructurePredictor::alpha(const string &sequence, const int &i, const int &j) {

    if(isComplement(sequence[i], sequence[j])){
        return 1;
    }

    return 0;
}

bool SecondaryStructurePredictor::isComplement(const char &first, const char &second) {

    if(first == 'a' || first == 'A'){
        return second == 'u' || second == 'U';
    }

    if(first == 'u' || first == 'U'){
        return second == 'a' || second == 'A';
    }

    if(first == 'C' || first == 'c'){
        return second == 'G' || second == 'g';
    }

    return second == 'C' || second == 'c';
}

void
SecondaryStructurePredictor::backtrack(const vector<vector<int>> &dp, vector<pair<int, int>> &aligned, int i, int j, const string &sequence) {
    if( j <= i){
        return;
    }

    /*if( dp[i][j] == dp[i+1][j]){
        return backtrack(dp, aligned, i+1, j, sequence);
    }*/
    if( dp[i][j] == dp[i][j-1]){
        return backtrack(dp, aligned, i, j-1, sequence);
    }
    /*if(isValidIndex(dp, i+1, j-1) && dp[i+1][j-1] == dp[i][j] + alpha(sequence, i,j)){

        if(alpha(sequence, i, j) > 0){
            aligned.push_back(make_pair(i, j));
        }
        backtrack(dp, aligned, i+1, j-1, sequence);
        return;
    }*/

    for(int k = i; k < j; k++){
        if(isComplement(sequence[k], sequence[j])) {
            if (k - 1 < 0) {
                if (k+1 < sequence.size() && j-1 > 0 && dp[i][j] == dp[k + 1][j - 1] + 1) {
                    aligned.emplace_back(k, j);
                    backtrack(dp, aligned, k + 1, j - 1, sequence);
                    return;
                }else{
                    return;
                }
            } else if (dp[i][j] == dp[i][k - 1] + dp[k + 1][j - 1] + 1) {
                aligned.emplace_back(k, j);
                backtrack(dp, aligned, i, k - 1, sequence);
                backtrack(dp, aligned, k + 1, j - 1, sequence);
                return;
            }
        }
    }

    /*
    if(dp[i][j] == dp[i+1][j]){
        return backtrack(dp, aligned, i+1, j, sequence);
    }*/

    /*if(dp[i][j] == dp[i+1][j-1] + alpha(sequence, i, j)){
        aligned.emplace_back(i, j);
        return backtrack(dp, aligned, i+1, j-1, sequence);
    }*/



}



void SecondaryStructurePredictor::getDotBracket(const vector<pair<int,int>> &aligned, const string &sequence){

    string ans;

    for(const auto &s : sequence){
        ans+=".";
    }

   for( const auto &p : aligned){

       ans[min(p.first, p.second)] = '(';
       ans[max(p.first, p.second)] = ')';

   }
   cout << sequence << endl;
   cout << ans << endl;
}

int SecondaryStructurePredictor::best(string sequence, int i, int j) {
    if(j <= i){
        return 0;
    }

    int answer = 0;


    vector<vector<int>> dummyMatrix(sequence.size(), vector<int>(sequence.size()));
    if(isValidIndex(dummyMatrix, i, j))
        answer = max(best(sequence, i, j-1), answer);
    for(int k = i; k < j; k++){
        if(isComplement(sequence[k], sequence[j])){
            if(isValidIndex(dummyMatrix, i, k-1) && isValidIndex(dummyMatrix, k+1, j-1)){
                answer = max(1 + best(sequence, i, k-1) + best(sequence, k+1, j-1), answer);
            }
        }
    }

    return answer;
}

void SecondaryStructurePredictor::writeComplements(const string &sequence, const vector< pair<int,int> > & pairs) {

    ofstream file("../results/complements.txt");

    file << sequence << endl;

    for(auto p : pairs){

        file << p.first << " " << p.second << endl;
    }

    file.close();

    system("python3.6 ../draw_rna.py");
}
