#include <iostream>
#include <cstdlib>
#include <cmath>

#include "TFile.h"
#include "TTree.h"
#include "TRandom3.h"
#include "TVectorT.h"
#include "TVector3.h"

using namespace std;

int main(int argc, char** argv){

	if (argc < 3)
    {
      cerr << "\nWrong number of arguments! Instead use:\n"
	   << "\tcombineFiles /path/to/combined/root/ /path/to/allInput/root/ \n";
      exit(-1);
    }

    int numFiles = argc - 2;
    cout << "Number of Files to Combine: " << numFiles << "\n";
    for( int i = 0 ; i < numFiles ; ++i){
    	cout << "\tWorking on file: " << argv[i+2] << "\n";
    	
    	TFile * inFile = new TFile(argv[i+2]);
   		TTree * inTree = (TTree*)inFile->Get("ResTree");
    }

    /*int numFiles = len(argv[1]);
    for(int i; i < numFiles; ++i){
    	cout << argv[1][i] << "\n";

    }
	*/
    
}