#include <iostream>
#include <prop_tree.h>
#include <TFile.h>
#include <TTree.h>


using namespace std;

int main(int argc, char** argv){
	if (argc != 4)
    {
      cerr << "Wrong number of arguments! Instead use:\n"
	   << "\tdigitize /path/to/generator/root/file /path/to/geant/root/file /path/to/output/file\n\n";
      exit(-1);
    }

	// Get the files and trees

	// generated neutron 
	TFile * inFile1 = new TFile(argv[1]);
	TTree * inTree1 = (TTree*)inFile1->Get("MCout");

	// output from geant based on generated neutron
	TFile * inFile2 = new TFile(argv[2]);
	TTree * inTree2 = (TTree*)inFile2->Get("PropTree");

	
	TTree * outTree = new TTree("ReconTree","Reconstructed Values");

	BAND_Event * event = NULL;

	inTree1-> SetBranchAddress("event",&event);

	/*for(int i =0; i<inTree1->GetEntries(); i++){
		inTree1->GetEntry(i);
		cout << "Working on event " << i << endl;

		for(int j=0; j<event->particles.size(); j++){

			cout << event->particles[j].E_dep << endl;

		}

	}*/


	return 0;
}