#include <iostream>
#include <cmath>
#include <cstdlib>

#include "TFile.h"
#include "TTree.h"
#include "TH2.h"
#include "TH1.h"
#include "TVector3.h"
#include "TF1.h"
#include "TVectorT.h"

#include "constants.h"


const double luminosity = 1e35*1e-33; 
const double runtime = 40*24*60*60; 
const double azim_CLAS12 = 0.5;

using namespace std;

int main(int argc, char ** argv){
  if (argc != 4){
    cerr << "Wrong number of arguments. Instead use\n"
    << "\t./finalHists /path/to/kinematic/file /path/to/generator/file /path/to/output/file\n";
    exit(-1);
  }

  TFile * outFile = new TFile(argv[3],"RECREATE");
  TH1D * Xp = new TH1D("Xp","Recon Q^2 > 2., W_prime > 1.8, Recon Pn > 0.25;xp;xp;Counts",14,0.1,0.8);
  
  TH1D * QSq_highX = new TH1D("QSq_highX","Recon Q^2 > 2., W_prime > 1.8, Recon Pn > 0.25;QSq;QSq;Counts",42,2,8);
  TH1D * Wp_highX = new TH1D("Wp_highX","Recon Q^2 > 2., W_prime > 1.8, Recon Pn > 0.25;Wp;Wp;Counts",13,1.8,3.1);
  TH1D * As_highX = new TH1D("As_highX","Recon Q^2 > 2., W_prime > 1.8, Recon Pn > 0.25;As;As;Counts",5,1.3,1.55);
  TH1D * Pn_highX = new TH1D("Pn_highX","Recon Q^2 > 2., W_prime > 1.8, Recon Pn > 0.25;Pn;Pn;Counts",6,0.28,0.42);
  TH1D * Theta_qn_highX = new TH1D("Theta_qn_highX","Recon Q^2 > 2., W_prime > 1.8, Recon Pn > 0.25;Theta_qn;Theta_qn;Counts",40,140,180);
  
  TH1D * QSq_lowX = new TH1D("QSq_lowX","Recon Q^2 > 2., W_prime > 1.8, Recon Pn > 0.25;QSq;QSq;Counts",42,2,8);
  TH1D * Wp_lowX = new TH1D("Wp_lowX","Recon Q^2 > 2., W_prime > 1.8, Recon Pn > 0.25;Wp;Wp;Counts",13,1.8,3.1);
  TH1D * As_lowX = new TH1D("As_lowX","Recon Q^2 > 2., W_prime > 1.8, Recon Pn > 0.25;As;As;Counts",5,1.3,1.55);
  TH1D * Pn_lowX = new TH1D("Pn_lowX","Recon Q^2 > 2., W_prime > 1.8, Recon Pn > 0.25;Pn;Pn;Counts",6,0.28,0.42);
  TH1D * Theta_qn_lowX = new TH1D("Theta_qn_lowX","Recon Q^2 > 2., W_prime > 1.8, Recon Pn > 0.25;Theta_qn;Theta_qn;Counts",40,140,180);

  // Set up the input file for all the kinematic variables for histograming
  TFile * inFile = new TFile(argv[1]);
  if (!(inFile->GetListOfKeys()->Contains("ResTree"))) cerr << "File has no entries\n";
  TTree * inTree = (TTree*)inFile->Get("ResTree");

  // Memory for the output tree
  double trueWp, trueXp, trueAs;
  double reconWp, reconXp, reconAs;
  double truePe, truePn, reconPe, reconPn;
  double truereconMom_e[3],trueMom_n[3];
  double reconreconMom_e[3],reconMom_n[3];

  double trueQSq, trueNu, trueX, true_q, truePhi_q, trueTheta_q;
  double reconQSq, reconNu, reconX, recon_q, reconPhi_q, reconTheta_q;
  double truePhi_e, truePhi_n, truePhi_en;
  double reconPhi_e, reconPhi_n, reconPhi_en;
  double trueTheta_n, reconTheta_n;
  double trueTheta_e, reconTheta_e;
  double trueEn, reconEn;
  double trueCosTheta_qn, reconCosTheta_qn;

  inTree->SetBranchAddress("trueWp",&trueWp);
  inTree->SetBranchAddress("trueXp",&trueXp);
  inTree->SetBranchAddress("trueAs",&trueAs);

  inTree->SetBranchAddress("reconWp",&reconWp);
  inTree->SetBranchAddress("reconXp",&reconXp);
  inTree->SetBranchAddress("reconAs",&reconAs);

  inTree->SetBranchAddress("truePe",&truePe);
  inTree->SetBranchAddress("truePn",&truePn);
  inTree->SetBranchAddress("reconPe",&reconPe);
  inTree->SetBranchAddress("reconPn",&reconPn);

  inTree->SetBranchAddress("truePe_x",&(truereconMom_e[0]));
  inTree->SetBranchAddress("truePe_y",&(truereconMom_e[1]));
  inTree->SetBranchAddress("truePe_z",&(truereconMom_e[2]));

  inTree->SetBranchAddress("truePn_x",&(trueMom_n[0]));
  inTree->SetBranchAddress("truePn_y",&(trueMom_n[1]));
  inTree->SetBranchAddress("truePn_z",&(trueMom_n[2]));

  inTree->SetBranchAddress("reconPe_x",&(reconreconMom_e[0]));
  inTree->SetBranchAddress("reconPe_y",&(reconreconMom_e[1]));
  inTree->SetBranchAddress("reconPe_z",&(reconreconMom_e[2]));

  inTree->SetBranchAddress("reconPn_x",&(reconMom_n[0]));
  inTree->SetBranchAddress("reconPn_y",&(reconMom_n[1]));
  inTree->SetBranchAddress("reconPn_z",&(reconMom_n[2]));

  inTree->SetBranchAddress("trueQSq",&(trueQSq));
  inTree->SetBranchAddress("trueNu",&(trueNu));
  inTree->SetBranchAddress("trueX",&(trueX));
  inTree->SetBranchAddress("true_q",&(true_q));
  inTree->SetBranchAddress("truePhi_q",&(truePhi_q));
  inTree->SetBranchAddress("trueTheta_q",&(trueTheta_q));

  inTree->SetBranchAddress("reconQSq",&(reconQSq));
  inTree->SetBranchAddress("reconNu",&(reconNu));
  inTree->SetBranchAddress("reconX",&(reconX));
  inTree->SetBranchAddress("recon_q",&(recon_q));
  inTree->SetBranchAddress("reconPhi_q",&(reconPhi_q));
  inTree->SetBranchAddress("reconTheta_q",&(reconTheta_q));

  inTree->SetBranchAddress("truePhi_e",&(truePhi_e));
  inTree->SetBranchAddress("truePhi_n",&(truePhi_n));
  inTree->SetBranchAddress("truePhi_en",&(truePhi_en));

  inTree->SetBranchAddress("reconPhi_e",&(reconPhi_e));
  inTree->SetBranchAddress("reconPhi_n",&(reconPhi_n));
  inTree->SetBranchAddress("reconPhi_en",&(reconPhi_en));

  inTree->SetBranchAddress("trueTheta_n",&(trueTheta_n));
  inTree->SetBranchAddress("reconTheta_n",&(reconTheta_n));
  inTree->SetBranchAddress("trueTheta_e",&(trueTheta_e));
  inTree->SetBranchAddress("reconTheta_e",&(reconTheta_e));

  inTree->SetBranchAddress("trueEn",&(trueEn));
  inTree->SetBranchAddress("reconEn",&(reconEn));

  inTree->SetBranchAddress("trueCosTheta_qn",&(trueCosTheta_qn));
  inTree->SetBranchAddress("reconCosTheta_qn",&(reconCosTheta_qn));



  // Create the weight unique to this file
  TFile * inGen = new TFile(argv[2]);
  TTree * inTreeGen = (TTree*)inGen->Get("MCout");
  double simulated_nEvents = inTreeGen->GetEntries(); 
  
  double total_num_events;
  TVectorT<double> *CSVec_dis = NULL;
  TVectorT<double> *CSVec_rand = NULL;
  CSVec_dis  = (TVectorT<double>*)inGen->Get("totalCS");
  CSVec_rand = (TVectorT<double>*)inGen->Get("totalCSSq");

  if (CSVec_dis){
    double CS = (*CSVec_dis)[0];
      total_num_events = CS * luminosity * runtime;
  }
  else{
    double CSSqDt = (*CSVec_rand)[0] * 1e-9;   // CSSqDt is product of the cross-sections and coincidence time window
    total_num_events = CSSqDt * pow(luminosity, 2) * runtime;
  }
    // Calculate weighting ratio
  double acceptance = azim_CLAS12 ;
  total_num_events *= acceptance;
  double weighting = total_num_events / simulated_nEvents;


  const int nEvents = inTree->GetEntries();
  //cout << "Filling 2D histograms..." << endl;
  for (int i=0 ; i < nEvents ; i++){
    //if (i % 10000 == 0) cout << "Event " << i << "\n";
    inTree->GetEvent(i);

    // Fill histograms
    Xp->Fill(reconXp,weighting);

    if (reconXp > 0.5){
      QSq_highX->Fill(reconQSq,weighting);
      Wp_highX->Fill(reconWp,weighting);
      As_highX->Fill(reconAs,weighting);
      Pn_highX->Fill(reconPn,weighting);
      Theta_qn_highX->Fill(acos(reconCosTheta_qn)*180./M_PI,weighting);
    }
    if (abs(reconXp-0.3) < 0.05){
      QSq_lowX->Fill(reconQSq,weighting);
      Wp_lowX->Fill(reconWp,weighting);
      As_lowX->Fill(reconAs,weighting);
      Pn_lowX->Fill(reconPn,weighting);
      Theta_qn_lowX->Fill(acos(reconCosTheta_qn)*180./M_PI,weighting);
    }

  }

  inFile->Close();
  inGen->Close();
  outFile->cd();
  Xp->Write();

  QSq_highX->Write();
  Wp_highX->Write();
  As_highX->Write();
  Pn_highX->Write();
  Theta_qn_highX->Write();

  QSq_lowX->Write();
  Wp_lowX->Write();
  As_lowX->Write();
  Pn_lowX->Write();
  Theta_qn_lowX->Write();

  outFile->Close();

}
