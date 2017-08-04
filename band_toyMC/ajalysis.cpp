#include <iostream>
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TVector3.h"
#include "constants.h"
using namespace std;


int main(int argc, char **argv)
{
    if (argc != 3)
    {
        cerr << "Wrong number of arguments, try:\n\t"
             << "ajalysis path/to/inputfile /path/to/outputfile\n";
        return -1;
    }
        
    TFile *infile = new TFile(argv[1]);
    TFile *outfile = new TFile(argv[2], "RECREATE");
    TTree *inTree = (TTree*)infile -> Get("ReconTree");
    
    double pe_true[3], pn_true[3], pe_recon[3], pn_recon[3];
    inTree->SetBranchAddress("truePe_x", &(pe_true[0]));
    inTree->SetBranchAddress("truePe_y", &(pe_true[1]));
    inTree->SetBranchAddress("truePe_z", &(pe_true[2]));
    inTree->SetBranchAddress("truePr_x", &(pn_true[0]));
    inTree->SetBranchAddress("truePr_y", &(pn_true[1]));
    inTree->SetBranchAddress("truePr_z", &(pn_true[2]));
    inTree->SetBranchAddress("reconPe_x", &(pe_recon[0]));
    inTree->SetBranchAddress("reconPe_y", &(pe_recon[1]));
    inTree->SetBranchAddress("reconPe_z", &(pe_recon[2]));
    inTree->SetBranchAddress("reconPr_x", &(pn_recon[0]));
    inTree->SetBranchAddress("reconPr_y", &(pn_recon[1]));
    inTree->SetBranchAddress("reconPr_z", &(pn_recon[2])); 
    
    
    TH1D *outhist1 = new TH1D("analysis", "Analsis; x-Axis; y-Axis", 100, 0, 1);
    TH2D *outhist2 = new TH2D("truevrecon_pn", "Comparison of True to Recon Neutron Momentum in x-Axis;Reconstructed Neutron Momentum; True Neutron Momentum; Events/bin", 100, -0.3, 0.3, 100, -0.3, 0.3);
    for(int iii = 0; iii < inTree->GetEntries(); ++iii)
    {
        inTree -> GetEvent(iii);
        TVector3 pn(pn_true[0], pn_true[1], pn_true[2]);
        TVector3 pe(pe_true[0], pe_true[1], pe_true[2]);
        TVector3 q = TVector3(0, 0, E1) - pe;
//        pe.Theta();
//        pe.Angle(q);
        outhist1 -> Fill(pn.Mag());
        outhist2 -> Fill(pn_recon[0], pn_true[0]);
    }
    
    infile -> Close();
    outhist1 -> Write();
    outhist2 -> Write();
    outfile -> Close();

    return 0;
}
