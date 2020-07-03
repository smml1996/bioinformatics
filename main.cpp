#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include "AdditiveTree.h"
#include "SecondaryStructurePredictor.h"

using namespace std;

string s1,s2;

vector<string>seqs;

void loadSequences(const string &path){
    seqs.clear();
    string line;

    ifstream myfile (path);

    if(myfile.is_open()){
        while(getline(myfile, line)){
            seqs.push_back(line);
        }

        myfile.close();
    }else{
        cout << "cant open" << endl;
    }
}


matrix readMatrix(){
    matrix v;
    ifstream myfile ("../matrix.txt");
    string line;
    if(myfile.is_open()){
        while(getline(myfile, line)){
            v.push_back(vector<double>());
            string temp;
            for(char i : line){

                if(i ==' '){
                    if(!temp.empty()){
                        if(temp != "-"){
                            v[v.size()-1].push_back(stod(temp));
                        }else{
                            v[v.size()-1].push_back(0);
                        }

                        temp.clear();
                    }
                }else{

                    temp+=i;
                }
            }

            if(!temp.empty()){
                if(temp != "-"){
                    v[v.size()-1].push_back(stod(temp));
                }else{
                    v[v.size()-1].push_back(0);
                }
            }

        }

        myfile.close();
    }else{
        cout << "cant open" << endl;
    }

    return v;

}

int main() {

    //matrix input = readMatrix();
    //loadSequences("../sequences.txt");
    ///AdditiveTree c;

    SecondaryStructurePredictor ssp;
/*
    int n;
    cout << "Enter num. seqs. to align: " ;
    cin >> n;
    int index;
    vector<string> seqsToAlign;
    while(n--){
        cout << "seq. index (1-indexed): ";

        cin >> index;

        index--;

        seqsToAlign.push_back(seqs[index]);
    }
*/
    ssp.predictSequence("CUGUUCUUGACA");
/*
    ProgressiveMSA pMSA{};

    pMSA.run(input);
*/
/*
    loadSequences("../sequences.txt");
    int n;
    cout << "Enter num. seqs. to align: " ;
    cin >> n;
    int index;
    vector<string> seqsToAlign;
    while(n--){
        cout << "seq. index (1-indexed): ";

        cin >> index;

        index--;

        seqsToAlign.push_back(seqs[index]);
    }
    ProgressiveMSA pMSA{};

    pMSA.run(seqsToAlign);*/
/*

    SmithWaterman nw;

    AlignmentAnswer as = nw.align("AATGCATGATCATGGATGGCAAGTCTCTGATTGCACAGCAGAAGGACTAAAGGTTGCACTCCTACTG", "GGGATTTCTTCATCATGTGGGAGAGCGTGTTCTGAACACTTGGCCATTTTCAATGCTAAGACAGAAGGCAATAGAAGTTGCTATTAATCATGTACGTTACG");
    cout <<"score: " << as.score << endl;
    cout <<"num. alineamientos: " << as.alignment.size() << endl;
    nw.printDpMatrix("lm23-36.txt"); // ver archivo results/matrix.txt
    nw.printAlignments("la23-36.txt"); // ver archivo results/alignments.txt
*/

    return 0;
}
