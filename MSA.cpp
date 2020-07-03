//
// Created by Stefanie Muroya lei on 4/15/20.
//

#include "MSA.h"
#include <iostream>
#include <algorithm>
#include <fstream>
#define MSA_MATRIX_PATH "results/msa1_matrix.txt"
#define MSA_AL_PATH "results/msa1_als.txt"

using namespace std;


int getColumnScore(vector<char> v){
    int score = 0;
    for(int i = 0; i < v.size(); i++){
        for(int j = i+1; j < v.size(); j++){
            if(v[i] == '-'){
                if(v[j] != '-'){
                    score = score -2;
                }
            }else if(v[j] == '-'){
                score-=2;
            }else{
                if(v[j] == v[i]){
                    score+=1;
                }else{
                    score-=1;
                }
            }
        }
    }

    return score;
}

int getScore(vector<string> v){
    int score = 0;

    vector<char>temp;
    for(int c = 0; c < v[0].size(); c++){

        temp.clear();
        for(int r = 0; r < v.size(); r++){
            temp.push_back(v[r][c]);
        }

        score+=getColumnScore(temp);
    }

    return score;
}
bool sortVs(vector<string> & v1 , vector<string> v2){

    int score1 =  getScore(v1), score2 = getScore(v2);

    return score1 > score2;

}

vector<string>  MSA::getAlignments(const vector<string> &seqs) {

    if(seqs.size() < 2){
        return seqs;

    }

    scoresSum = vector<int>(seqs.size(), 0);
    m = vector< vector<int> >( seqs.size(), vector<int>(seqs.size(), 0));
    center = 0;
    int currIndex = 0;
    int temp;
    AlignmentAnswer as;
    vector<string>  ans;

    this->calculateMatrix(seqs);
    this->saveMatrix();
    this->findCenter();

    cout << "centro: " << (center + 1) << endl;

    if(currIndex == center) currIndex++;

    as = nw.align(seqs[center], seqs[currIndex], true);
    cout << "center score: " << as.score << endl;

    int maxSize;

    maxSize = max(as.alignment[0].first.size(), as.alignment[0].second.size());

    ans.push_back(as.alignment[0].first);
    ans.push_back(as.alignment[0].second);


    currIndex++;

    int top;
    vector<string>temporal;

    vector< vector<string> > t2;


    while(currIndex < seqs.size()){

        if(currIndex != center){
            as = nw.align(seqs[center], seqs[currIndex], true);
            ans.push_back(as.alignment[0].second);
            maxSize = max(maxSize, (int)as.alignment[0].second.size());
        }
        currIndex++;
    }


    this->fillGaps(ans, maxSize);

    this->saveAlignments(ans);

    return ans;
}

void MSA::saveMatrix() {

    ofstream myfile;
    myfile.open (MSA_MATRIX_PATH);
    for(int i = 0; i < m.size(); i++){
        for(int j = 0; j < m[i].size(); j++){
            myfile << m[i][j] << " ";
        }

        myfile << scoresSum[i];

        myfile << endl;
    }

    for(int i=0; i < scoresSum.size(); i++){
        myfile << scoresSum[i] << " ";
    }

    cout << endl;

    myfile<< endl << endl;
    myfile.close();

}

void MSA::calculateMatrix(const vector<string> &seqs) {
    int score;
    for(int i = 0; i < seqs.size(); i++){

        for(int j = i+1; j < seqs.size(); j++){
            score = nw.getScore(seqs[i], seqs[j]);
            scoresSum[i]+= score;
            scoresSum[j]+= score;
            m[i][j] += score;
            m[j][i]+=score;
        }
    }
}

void MSA::findCenter() {
    for(int i = 1; i < scoresSum.size(); i++){
        if(scoresSum[i] > scoresSum[center]){
            center = i;
        }
    }
}

void MSA::fillGaps(vector<string>& ans, const int& maxSize) {

    for(auto & an : ans){
        while(an.size() < maxSize){
            an+="-";
        }
    }
}

void MSA::saveAlignments(const vector<string> &ans) {
    ofstream myfile;
    myfile.open (MSA_AL_PATH);

    for(int j = 0; j < ans.size(); j++){
        myfile << ans[j] << endl;
    }

    myfile.close();
}
