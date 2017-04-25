#include <iostream>
#include <cstdlib>
#include <cmath>

#include "TFile.h"
#include "TTree.h"
#include "TRandom3.h"
#include "TVectorT.h"

#include "constants.h"
#include "lad_helpers.h"

using namespace std;

const double tResPMT = 0.5; // 500 ps
const double cScint = 15.;

int main(int argc, char ** argv)
{
  TRandom3 * myRand = new TRandom3(0);

  if (argc != 3)
    {
      cerr << "Wrong number of arguments. Instead use\n\tdet_sim /path/to/input/file /path/to/output/file\n";
      exit(-1);
    }

  // Get the files
  TFile * inFile = new TFile(argv[1]);
  TFile * outFile = new TFile(argv[2],"RECREATE");

  // Copy the cross section
  TVectorT<double> * csVec;
  if (inFile->GetListOfKeys()->Contains("totalCS"))
    {
      csVec = (TVectorT<double>*)inFile->Get("totalCS");
      gFile=outFile;
      csVec->Write("totalCS");
    }
  else if (inFile->GetListOfKeys()->Contains("totalCSSq"))
    {
      csVec = (TVectorT<double>*)inFile->Get("totalCSSq");
      gFile=outFile;
      csVec->Write("totalCSSq");
    }

  // Get the input tree
  TTree * inTree = (TTree*) inFile->Get("PropTree");
  double gem1_x, gem1_y, gem1_z, gem2_x, gem2_y, gem2_z, lad_x, lad_y, lad_z, path_length, flight_time;
  int lad_plane;
  inTree->SetBranchAddress("gem1_x",&gem1_x);
  inTree->SetBranchAddress("gem1_y",&gem1_y);
  inTree->SetBranchAddress("gem1_z",&gem1_z);
  inTree->SetBranchAddress("gem2_x",&gem2_x);
  inTree->SetBranchAddress("gem2_y",&gem2_y);
  inTree->SetBranchAddress("gem2_z",&gem2_z);
  inTree->SetBranchAddress("lad_x",&lad_x);
  inTree->SetBranchAddress("lad_y",&lad_y);
  inTree->SetBranchAddress("lad_z",&lad_z);
  inTree->SetBranchAddress("path_length",&path_length);
  inTree->SetBranchAddress("flight_time",&flight_time);
  inTree->SetBranchAddress("lad_plane",&lad_plane);

  // Set up the output tree
  outFile->cd();
  TTree * outTree = new TTree("ReconTree","Reconstructed values for the recoil protons");
  double mom_recon, theta_recon, phi_recon, z_recon;
  outTree->Branch("mom_recon",&mom_recon,"mom_recon/D");
  outTree->Branch("theta_recon",&theta_recon,"theta_recon/D");
  outTree->Branch("phi_recon",&phi_recon,"phi_recon/D");
  outTree->Branch("z_recon",&z_recon,"z_recon/D");
  
  int nEvents = inTree->GetEntries();
  for (int event=0 ; event<nEvents ; event++)
    {
      if (event%100000==0)
	cerr << "Working on event " << event << "\n";

      inTree->GetEvent(event);

      if (lad_plane >= 0)
	{      
	  // Smear GEMs
	  double gem1_local_x = (lad_gem_radii[0]*sin(lad_angles[1]) - gem1_x)/cos(lad_angles[1]);
	  double gem1_local_x_smeared = gem1_local_x + myRand->Gaus()*0.05; // 0.05 cm smearing
	  double gem1_local_y_smeared = gem1_y + myRand->Gaus()*0.05;
	  double gem1_x_smeared = lad_gem_radii[0]*sin(lad_angles[1]) - gem1_local_x_smeared*cos(lad_angles[1]);
	  double gem1_y_smeared = gem1_local_y_smeared;
	  double gem1_z_smeared = lad_gem_radii[0]*cos(lad_angles[1]) + gem1_local_x_smeared*sin(lad_angles[1]);
  
	  double gem2_local_x = (lad_gem_radii[1]*sin(lad_angles[1]) - gem2_x)/cos(lad_angles[1]);
	  double gem2_local_x_smeared = gem2_local_x + myRand->Gaus()*0.05;
	  double gem2_local_y_smeared = gem2_y + myRand->Gaus()*0.05;
	  double gem2_x_smeared, gem2_y_smeared, gem2_z_smeared;
	  gem2_x_smeared = lad_gem_radii[1]*sin(lad_angles[1]) - gem2_local_x_smeared*cos(lad_angles[1]);
	  gem2_y_smeared = gem2_local_y_smeared;
	  gem2_z_smeared = lad_gem_radii[1]*cos(lad_angles[1]) + gem2_local_x_smeared*sin(lad_angles[1]);
	  
	  // Reconstruct the vertex
	  z_recon = (gem1_z_smeared*gem2_x_smeared - gem2_z_smeared*gem1_x_smeared)/(gem2_x_smeared - gem1_x_smeared);
	    //gem1_z_smeared - (gem2_z_smeared-gem1_z_smeared)/(gem2_x_smeared-gem1_x_smeared)*gem1_x_smeared;
	  
	  // Angle reconstruction
	  theta_recon = acos((gem2_z_smeared-gem1_z_smeared)/
			     sqrt(sq(gem2_x_smeared-gem1_x_smeared)+
				  sq(gem2_y_smeared-gem1_y_smeared)+
				  sq(gem2_z_smeared-gem1_z_smeared)));
	  phi_recon = atan((gem2_y_smeared-gem1_y_smeared)/
			   sqrt(sq(gem2_z_smeared-gem1_z_smeared) + sq(gem2_x_smeared-gem1_x_smeared)));
	  
	  // Reconstruct at LAD
	  double lad_local_x = distInPlane(TVector3(lad_x,lad_y,lad_z), lad_radii[lad_plane], lad_angles[lad_plane])
	    * ((lad_z > lad_radii[lad_plane]*cos(lad_angles[lad_plane]))? 1. : -1.);
	  double lad_local_x_smeared = 22. * floor((lad_local_x+11.)/22.);
	  double lad_local_y_smeared = lad_y + myRand->Gaus() * cScint * tResPMT/sqrt(2.);
	  double lad_x_smeared = lad_radii[lad_plane]*sin(lad_angles[lad_plane]) - lad_local_x_smeared*cos(lad_angles[lad_plane]);
	  double lad_y_smeared = lad_local_y_smeared;
	  double lad_z_smeared = lad_radii[lad_plane]*cos(lad_angles[lad_plane]) + lad_local_x_smeared*sin(lad_angles[lad_plane]);
	  double recon_path_length = sqrt(sq(lad_x_smeared) + sq(lad_y_smeared) + sq(lad_z_smeared - z_recon));
	  
	  // Smear LAD Timing
	  double flight_time_recon = flight_time + myRand->Gaus()*tResPMT/sqrt(2.);
	  double beta_recon = recon_path_length/(cAir * flight_time_recon);
	  if (beta_recon > 1.) beta_recon=0.99999999;
	  mom_recon = mP/sqrt(1./sq(beta_recon) - 1.);
	 
	} 
      // Write the tree
      outTree->Fill();
    }

  inFile->Close();

  outFile->cd();
  outTree->Write();
  outFile->Close();

  delete myRand;

  return 0;
}
