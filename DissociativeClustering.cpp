//
// Created by Stefanie Muroya lei on 5/11/20.
//

#include "DissociativeClustering.h"
#include <iostream>
#include <fstream>
#include <cmath>
#include <climits>
#include <queue>

#include "NeedlemanWunsch.h"

using namespace std;

void DissociativeClustering::run(const vector<string> &seqs, const string &filePath) {
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
    return this->run(d, filePath, true);

}

void DissociativeClustering::run(matrix distances, const string &filePath, bool appendDistances) {

    ProgressiveMSA::normalizeD(distances);
    priority_queue< pair<double, int> >pq;
    ProgressiveMSA::printD(distances, "../results/distances.txt", false);
    vector<bool> isUsed(distances.size(), false);
    vector<double> minDistances(distances.size());
    for(int i = 0; i < distances.size(); i++){
        minDistances[i] = getLess(distances,i);
        pq.push(make_pair(minDistances[i], i));
    }

    for(auto r : distances){
        for(auto c: r){
            cout << c <<  "\t";
        }
        cout << endl;
    }

    ofstream file("../results/track.txt");

    for(int i = 0; i < minDistances.size(); i++){
        file << minDistances[i] << "\t";
    }
    file << endl << endl;

    vector< vector<int> > clusters;
    printClusters(clusters);
    vector<int> temp;
    while(!pq.empty()){
        temp.clear();

        if(!isUsed[pq.top().second]){
            file << "Escogido: " << (char)(pq.top().second + 'A') << endl;
            isUsed[pq.top().second] = true;
            temp.push_back(pq.top().second);

            for(int i = 0; i < isUsed.size(); i++){
                if(!isUsed[i]){
                    if(minDistances[i] - distances[i][pq.top().second] > 0){
                        isUsed[i] = true;
                        temp.push_back(i);
                    }
                }
            }
            clusters.push_back(temp);
            for(int i = 0; i < minDistances.size(); i++){
                if(!isUsed[i]){
                    file << (char) (i + 'A') <<": "<< minDistances[i] << "\t\t";
                }

            }
        }

        file << endl << endl;
        printClusters(clusters, true);
        pq.pop();

    }

    file.close();

}

void DissociativeClustering::printClusters(vector<vector<int> > &clusters, bool isAppend) {

    ofstream file;

    if(isAppend){
        file =  ofstream("../results/clusters.txt",ios::out | ios::app);
    }else{
        file = ofstream ("../results/clusters.txt");
    }

    file << "NUEVO_PONIS" << endl;
    for(const auto & c : clusters){
        file << "(";
        for(const auto &n : c){
            file << (char)(n + 'A') << "\t";
        }
        file << ")";
    }
    file << endl <<"END_PONIS" << endl;
    file.close();
}

double DissociativeClustering::getLess(matrix distances, int t) {

    double ans = INT_MAX;
    for(int i = 0; i < distances.size(); i++){
        if(i != t)
            ans = min(distances[t][i], ans);
    }

    return ans;
}


