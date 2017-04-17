#include <iostream>
#include <prop_tree.h>
#include <TFile.h>
#include <TTree.h>


using namespace std;

int main(int argc, char** argv){
	if (argc != 3)
    {
      cerr << "Wrong number of arguments! Instead use:\n"
	   << "\tdigitize /path/to/inputGeant/root/file /path/to/output/file\n\n";
      exit(-1);
    }


	// Get the files and trees
	TFile * inFile = new TFile(argv[1]);
	TFile * outFile = new TFile(argv[2],"RECREATE");
	TTree * inTree = (TTree*)inFile->Get("PropTree");
	TTree * outTree = new TTree("ReconTree","Reconstructed Values");

	BAND_Event * event = NULL;

	inTree-> SetBranchAddress("band",&event);

	for(int i =0; i<inTree->GetEntries(); i++){
		inTree->GetEntry(i);
		cout << "Working on event " << i << endl;

		for(int j=0; j<event->hits.size(); j++){

			cout << event->hits[j].E_dep << endl;

		}

	}


	cout << "a " << endl;
	return 0;
}