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
#include "TLatex.h"

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

  //
  TH2D * pE_thetaE = new TH2D("pE_thetaE","Recon Q^2 > 2., W_prime > 1.8, 0.6 > Recon Pn > 0.25, Recon Xp < 1.",150,5,35,50,2,7);
  pE_thetaE->GetXaxis()->SetTitle("Electron Scattering Angle [^{#bullet}]");
  pE_thetaE->GetYaxis()->SetTitle("Electron Momentum [GeV/c]");

  TH2D * pN_thetaNQ = new TH2D("pN_thetaNQ","Recon Q^2 > 2., W_prime > 1.8, 0.6 > Recon Pn > 0.25, Recon Xp < 1.",100,130,180,90,0.2,0.65);
  pN_thetaNQ->GetXaxis()->SetTitle("Angle Between Neutron and Q Vector [^{#bullet}]");
  pN_thetaNQ->GetYaxis()->SetTitle("Neutron Momentum [GeV/c]");
  
  TH2D * pE_Xp = new TH2D("pE_Xp","Recon Q^2 > 2., W_prime > 1.8, 0.6 > Recon Pn > 0.25, Recon Xp < 1.",70,0.1,0.8,60,2,8);
  pE_Xp->GetXaxis()->SetTitle("x_{p}");
  pE_Xp->GetYaxis()->SetTitle("Electron Momentum [GeV/c]");
  
  TH2D * thetaE_Xp = new TH2D("thetaE_Xp","Recon Q^2 > 2., W_prime > 1.8, 0.6 > Recon Pn > 0.25, Recon Xp < 1.",70,0.1,0.8,150,5,35);
  thetaE_Xp->GetXaxis()->SetTitle("x_{p}");
  thetaE_Xp->GetYaxis()->SetTitle("Electron Scattering Angle [^{#bullet}]");

  TH2D * QSq_Xp = new TH2D("QSq_Xp","Recon Q^2 > 2., W_prime > 1.8, 0.6 > Recon Pn > 0.25, Recon Xp < 1.",70,0.1,0.8,60,2,8);
  QSq_Xp->GetXaxis()->SetTitle("x_{p}");
  QSq_Xp->GetYaxis()->SetTitle("Q^{2} [GeV/c]^{2}");

  TH2D * As_Xp = new TH2D("As_Xp","Recon Q^2 > 2., W_prime > 1.8, 0.6 > Recon Pn > 0.25, Recon Xp < 1.",70,0.1,0.8,50,1.3,1.55);
  As_Xp->GetXaxis()->SetTitle("x_{p}");
  As_Xp->GetYaxis()->SetTitle("#alpha_{S}");

  Float_t xp_bins[ ] = {0.25,0.35,0.5,0.9};
  Int_t binnumXp = sizeof(xp_bins)/sizeof(Float_t) - 1;
  Float_t as_bins[ ] = {1.3,1.35,1.4,1.45,1.5,1.55};
  Int_t binnumAs = sizeof(as_bins)/sizeof(Float_t) - 1;
  TH2D * falsePos = new TH2D("falsePos","#frac{(True Fail) & (Smeared Pass)}{Smeared Pass};Recon x';Recon a_{s};Counts",binnumXp,xp_bins,binnumAs,as_bins);
  TH2D * falseNeg = new TH2D("falseNeg","#frac{(True Pass) & (Smeared Fail)}{Smeared Pass};Recon x';Recon a_{s};Counts",binnumXp,xp_bins,binnumAs,as_bins);  
  TH2D * changeYield = new TH2D("changeYield","False Negative - False Positive [%];Recon x';Recon a_{s};Counts",binnumXp,xp_bins,binnumAs,as_bins);  

  double denom[binnumXp][binnumAs];
  for(int n = 0; n < binnumAs; n++){
    for(int m = 0; m < binnumXp; m++){
      denom[m][n] = 0;
    }
  }


  TH1D * Xp = new TH1D("Xp","Recon Q^2 > 2., W_prime > 1.8, 0.6 > Recon Pn > 0.25, Recon Xp < 1. ;xp;xp;Counts",14,0.1,0.8);
  
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

    // Determine if there is bin migration
    for(int n = 0; n < binnumAs; n++){
      for(int m = 0; m < binnumXp; m++){

        if( m == 1 ) continue;

        double L = xp_bins[m];
        double R = xp_bins[m+1];

        double B = as_bins[n];
        double U = as_bins[n+1];

        // If smeared is in the bin
        if( (reconXp >= L) && (reconXp < R) ){
          if( (reconAs >= B) && (reconAs < U) ){

            denom[m][n]++;
            //cout << "In bin " << L << " " << R << " " << B << " " << U << "\n";

            // If true ISN'T in the bin
            if( !( (trueXp >= L) && (trueXp < R) ) || !( (trueAs >= B) && (trueAs < U) ) ){
                falsePos->Fill(reconXp,reconAs);
                //cout << "\tfalse pos: " << trueXp << " " << trueAs << "\n";

            }
              //else cout << "\tgood: " << trueXp << " " << trueAs << "\n";

          }
        }

        // If true is in the bin
        if( (trueXp >= L) && (trueXp < R) ){
          if( (trueAs >= B) && (trueAs < U) ){

            // If smeared ISN'T in the bin
            if( !( (reconXp >= L) && (reconXp < R) ) || !( (reconAs >= B) && (reconAs < U) ) ){
                falseNeg->Fill(trueXp,trueAs);
            }

          }
        }


      }
    }

    // Fill histograms
    pE_thetaE->Fill(reconTheta_e*180./M_PI,reconPe,weighting);
    pE_Xp->Fill(reconXp,reconPe,weighting);
    thetaE_Xp->Fill(reconXp,reconTheta_e*180./M_PI,weighting);
    QSq_Xp->Fill(reconXp,reconQSq,weighting);
    As_Xp->Fill(reconXp,reconAs,weighting);
    pN_thetaNQ->Fill( acos(reconCosTheta_qn)*180./M_PI , reconPn );


    Xp->Fill(reconXp,weighting);

    if ((reconXp > 0.5) && (reconXp < 0.9)){
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

  pE_thetaE->Write();
  pN_thetaNQ->Write();
  pE_Xp->Write();
  thetaE_Xp->Write();
  QSq_Xp->Write();
  As_Xp->Write();



  cout << "# AsBinCenter\tXp_High[0.7]_Low[0.3]\tFalseNeg\tFalsePos\tChangeYield\n";
  for(int n = 0; n < binnumAs; n++){
    for(int m = 0; m < binnumXp; m++){
      if( m == 1 ) continue;

      falsePos->SetBinContent( m+1,n+1, falsePos->GetBinContent(m+1,n+1) * 100. / (denom[m][n])  );
      falseNeg->SetBinContent( m+1,n+1, falseNeg->GetBinContent(m+1,n+1) * 100. / (denom[m][n])  );

      changeYield->SetBinContent( m+1,n+1, falseNeg->GetBinContent(m+1,n+1) - falsePos->GetBinContent(m+1,n+1) );

      cout << (as_bins[n]+as_bins[n+1])/2. << "\t" << (xp_bins[m]+xp_bins[m+1])/2. << "\t" <<  falseNeg->GetBinContent(m+1,n+1) << "\t" << falsePos->GetBinContent(m+1,n+1) << "\t" << changeYield->GetBinContent(m+1,n+1) << "\n";

      }
  } 

 
  falsePos->Write();
  falseNeg->Write();
  changeYield->Write();
  
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
