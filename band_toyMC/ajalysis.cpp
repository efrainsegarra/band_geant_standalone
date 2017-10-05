#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TMath.h"
#include "TVectorT.h"
#include "TVector3.h"
#include "constants.h"  // contains constants E1
using namespace std;

const double luminosity = 1e35*1e-33;  // convert luminosity from cm^-2*s-1 -> nb-1*s-1
const double runtime = 40*24*60*60; // convert 40 day runtime to s
const double azim_CLAS12 = 0.5;  // 50% azimuthal electron coverage for CLAS12
const double azim_BAND = 0.6;    // 60% neutron coverage for BAND

int main(int argc, char **argv)
{
    // CLI arguments handling
    if (argc != 4)
    {
        cerr << "Wrong number of arguments, try:\n\t"
             << "ajalysis  path/to/inputfile  /path/to/outputfile  "
             << "Threshold_eeEDep" << endl;
        return -1;
    }
        
    // Unpack user input for input file, output file and band detection efficiency
    TFile *infile = new TFile(argv[1]);
    TFile *outfile = new TFile(argv[2], "RECREATE");
    double thresh_eeEDep = atof(argv[3]);
    
    // Get values from ReconTree
    TTree *inTree = (TTree*)infile -> Get("ReconTree");

    int num_events_sim = inTree->GetEntries();  // number of events generated in simulation
    double truePe[3], truePn[3], reconPe[3], reconPn[3];
    double reconWp, reconAs, reconQSq, reconXp, recon_eeEDep;
    inTree->SetBranchAddress("truePe_x", &(truePe[0]));
    inTree->SetBranchAddress("truePe_y", &(truePe[1]));
    inTree->SetBranchAddress("truePe_z", &(truePe[2]));
    inTree->SetBranchAddress("truePr_x", &(truePn[0]));
    inTree->SetBranchAddress("truePr_y", &(truePn[1]));
    inTree->SetBranchAddress("truePr_z", &(truePn[2]));
    inTree->SetBranchAddress("reconPe_x", &(reconPe[0]));
    inTree->SetBranchAddress("reconPe_y", &(reconPe[1]));
    inTree->SetBranchAddress("reconPe_z", &(reconPe[2]));
    inTree->SetBranchAddress("reconPr_x", &(reconPn[0]));
    inTree->SetBranchAddress("reconPr_y", &(reconPn[1]));
    inTree->SetBranchAddress("reconPr_z", &(reconPn[2]));
    inTree->SetBranchAddress("reconQSq", &(reconQSq));
    inTree->SetBranchAddress("reconXp", &(reconXp)); 
    inTree->SetBranchAddress("reconWp", &(reconWp));
    inTree->SetBranchAddress("reconAs", &(reconAs));
    inTree->SetBranchAddress("recon_eeEDep", &(recon_eeEDep));

    // Calculate the efficiency of BAND
    int num_src = 0;
    int num_src_hits = 0;

    for(int iii = 0; iii < num_events_sim; ++iii)
    {
        inTree -> GetEvent(iii);

        TVector3 rePn(reconPn[0], reconPn[1], reconPn[2]);
        TVector3 rePe(reconPe[0], reconPe[1], reconPe[2]);
        TVector3 q = TVector3(0, 0, E1) - rePe;
        double theta_nq = rePn.Angle(q) * 180/TMath::Pi();

        // DIS Selection Cuts
        if (reconQSq < 2) continue;
        if (reconWp < 1.8) continue;        
        if (theta_nq < 110) continue;
        ++num_src;

        // Band Selection Cut
        if (recon_eeEDep < thresh_eeEDep) continue;
        ++num_src_hits;
    }

    double band_eff = static_cast<double>(num_src_hits) / num_src;

    // Attempt to pull CS data from input root file
    TVectorT<double> *CSVec = NULL;
    CSVec = (TVectorT<double>*)infile->Get("totalCS");
    
    // Test if input root file from signal-generated run
    double tot_num_events;
    if (CSVec)
    {
        double CS = (*CSVec)[0];    // cross section in units of nb  
        tot_num_events = CS * luminosity * runtime;
    }
    // or from background run
    else
    {
         CSVec = (TVectorT<double>*)infile->Get("totalCSSq");
         
         // Error handling if no CS found in file
         if (!CSVec)
         {
            cerr << "No cross section was found in the file";
            return -2;
         }

         double CSSqDt = (*CSVec)[0] * 1e-9;   // CSSqDt is product of the cross-sections and coincidence time window
         tot_num_events = CSSqDt * pow(luminosity, 2) * runtime;
    }

    // Calculate weighting ratio
    double acceptance = azim_CLAS12 * azim_BAND * band_eff;
    double num_events_detected = tot_num_events * acceptance;
    double weighting = num_events_detected / num_events_sim;

    // Hist output 
    TH1D *XpHist = new TH1D("XpHist",
    "After DIS, Backward Recoil Neutron & eeEDep_{n,r|thresh} Selection Cuts;#bf{#it{x_{B}}};Counts", 14, 0.1, 0.8);
    TH2D *outhist1 = new TH2D("h1", 
    "Before Cuts;#bf{#it{#theta_{e,r}}} #left[#it{deg.}#right];#left|#bf{#it{#vec{p}_{e,r}}}#right| #left[#frac{GeV}{c}#right]", 
    100, 4, 35, 100, 1.5, 8);
    TH2D *outhist2 = new TH2D("h2", 
    "After DIS, Backward Recoil Neutron & eeEDep_{n,r|thresh};#bf{#it{#theta_{#vec{p}_{n,r}#vec{q}}}} #left[#it{deg.}#right];#left|#bf{#it{#vec{p}_{n,r}}}#right| #left[#frac{GeV}{c}#right]",
    100, 135, 180, 100, 0.26, 0.65);
    TH2D *outhist3 = new TH2D("h3", 
    "After DIS, Backward Recoil Neutron & eeEDep_{n,r|thresh};#bf{#it{x_{B}}};#bf{#it{Q^{2}}} #left[#left(#frac{GeV}{c}#right)^{2} #right]", 
    100, 0.1, 0.8, 100, 2, 8);

    // Data (text file) output
    ofstream datfile;
    string root_filename(argv[2]);
    size_t raw_index = root_filename.find_last_of(string("."));
    string filename = root_filename.substr(0, raw_index);
    string dat_filename = filename + "_dat.txt";
    datfile.open(dat_filename.c_str());

    for(int iii = 0; iii < num_events_sim; ++iii)
    {
        inTree -> GetEvent(iii);

        TVector3 rePn(reconPn[0], reconPn[1], reconPn[2]);
        TVector3 rePe(reconPe[0], reconPe[1], reconPe[2]);
        TVector3 q = TVector3(0, 0, E1) - rePe;
        double theta_nq = rePn.Angle(q) * 180/TMath::Pi();

        // DIS Selection Cuts
        if (reconQSq < 2) continue;
        if (reconWp < 1.8) continue;        
        if (theta_nq < 110) continue;

        // Band Selection Cut
        if (recon_eeEDep < thresh_eeEDep) continue;
        
        // Fill histograms
        XpHist->Fill(reconXp, weighting);
        
        double theta_e = rePe.Theta() * 180/TMath::Pi();
        double rePe_mag = rePe.Mag();
        outhist1 -> Fill(theta_e, rePe_mag, weighting);

        double pn_mag = rePn.Mag();
        outhist2 -> Fill(theta_nq, pn_mag, weighting);

        outhist3 -> Fill(reconXp, reconQSq, weighting);
    }

    // File I/O
    infile -> Close();
    CSVec -> Write("CSdata");
    XpHist -> Write();
    outhist1 -> Write();
    outhist2 -> Write();
    outhist3 -> Write();
    outfile -> Close();

    datfile << thresh_eeEDep << "\t" << band_eff << endl;
    datfile.close();

    return 0;
}
