#include <iostream>
#include <cstdlib>
#include <cmath>

#include "TFile.h"
#include "TTree.h"
#include "TRandom3.h"

using namespace std;

// Define some constants for BAND (positions in cm, time in ns, mass in GeV)
const double bandZ = -240.;
const double barWidth = 7.;
const double tResPMT = 0.200;
const double cAir = 29.9792458;
const double cScint = cAir / 1.58; // A more accurate number will eventually need to be supplied
const double mP = 0.93827208;
const double mN = 0.93956541;
const double mD = 1.875612928;
const double E1 = 10.9;

inline double sq(double x){ return x*x;};

int main(int argc, char ** argv)
{
  if (argc != 3)
    {
      cerr << "Wrong number of arguments. Instead use\n\tdet_sim /path/to/input/file /path/to/output/file\n";
      exit(-1);
    }

  // Random number generator
  TRandom3 * myRand = new TRandom3(0);

  // Get the files and trees
  TFile * inFile = new TFile(argv[1]);
  TFile * outFile = new TFile(argv[2],"RECREATE");
  TTree * inTree = (TTree*)inFile->Get("MCout");
  TTree * outTree = new TTree("ReconTree","Reconstructed Values");

  // Memory for the input branches
  double mom_e[3];
  double mom_r[3];
  inTree->SetBranchAddress("x_e",&(mom_e[0]));
  inTree->SetBranchAddress("y_e",&(mom_e[1]));
  inTree->SetBranchAddress("z_e",&(mom_e[2]));
  inTree->SetBranchAddress("x_r",&(mom_r[0]));
  inTree->SetBranchAddress("y_r",&(mom_r[1]));
  inTree->SetBranchAddress("z_r",&(mom_r[2]));

  // Memory for the output tree
  double trueWp, trueXp, trueAs, truePr;
  double reconWp, reconXp, reconAs, reconPr;
  outTree->Branch("trueWp",&trueWp,"trueWp/D");
  outTree->Branch("trueXp",&trueXp,"trueXp/D");
  outTree->Branch("trueAs",&trueAs,"trueAs/D");
  outTree->Branch("truePr",&truePr,"truePr/D");
  outTree->Branch("reconWp",&reconWp,"reconWp/D");
  outTree->Branch("reconXp",&reconXp,"reconXp/D");
  outTree->Branch("reconAs",&reconAs,"reconAs/D");
  outTree->Branch("reconPr",&reconPr,"reconPr/D");

  // Loop over the events
  for (int i=0 ; i<inTree->GetEntries() ; i++)
    {
      if (i%100000==0)
	cerr << "Working on event " << i << "\n";

      // Load the event
      inTree->GetEvent(i);

      // Figure out where the recoiling neutron hits BAND
      double trueZ = bandZ;
      double trueX = mom_r[0]*trueZ/mom_r[2];
      double trueY = mom_r[1]*trueZ/mom_r[2];
      truePr = sqrt(sq(mom_r[0]) + sq(mom_r[1]) + sq(mom_r[2]));
      double trueEr = sqrt(sq(truePr) + sq(mN));
      double trueT = sqrt(sq(trueX)+sq(trueY)+sq(trueZ)) / (cAir * truePr / trueEr); // d / beta

      // Flight time smearing is based on PMT resolution for mean time
      double reconT = myRand->Gaus(trueT,tResPMT/sqrt(2.));

      // X reconstruction is smeared based on PMT resolution for time difference
      double reconX = myRand->Gaus(trueX,cScint * tResPMT/sqrt(2.));

      // Y reconstruction is based on center of the bar
      int bar = floor(trueY / barWidth);
      double reconY = barWidth* (((float)bar) + 0.5);

      // Z is not smeared
      double reconZ = trueZ;

      // Momentum
      double reconPath = sqrt(sq(reconX)+sq(reconY)+sq(reconZ));
      double reconBeta = reconPath/(cAir*reconT);
      reconPr = mN / sqrt( sq(1./reconBeta) - 1.);
      double reconEr = sqrt(sq(reconPr)+sq(mN));

      // Lepton arm is not smeared
      double p_e =sqrt(sq(mom_e[0]) + sq(mom_e[1]) + sq(mom_e[2]));
      double cosTheta_e = mom_e[2]/p_e;
      double theta_e = acos(cosTheta_e);
      double phi_e = atan2(mom_e[1],mom_e[0]);
      double QSq = 2.*E1*p_e*(1.-cosTheta_e);
      double nu = E1 - p_e;
      double x = QSq / (2.*mP*nu);
      double q = sqrt(QSq + sq(nu));
      double phi_q = phi_e+M_PI;
      if (phi_q > M_PI) phi_q -= 2.*M_PI;
      double theta_q = acos((E1 - p_e*cos(theta_e))/q);

      // Proton arm is smeared
      double trueCosTheta_r = mom_r[2]/truePr;
      double trueTheta_r = acos(trueCosTheta_r);
      double truePhi_r = atan2(mom_r[1],mom_r[0]);
      double reconCosTheta_r = reconZ/reconPath;
      double reconTheta_r = acos(trueCosTheta_r);
      double reconPhi_r = atan2(trueY,trueX);
      
      // Let's take care of those pesky relative phi angles
      double truePhi_er = truePhi_r - phi_e;
      if (truePhi_er < -M_PI) truePhi_er += 2.*M_PI;
      if (truePhi_er > M_PI) truePhi_er -= 2.*M_PI;
      double reconPhi_er = reconPhi_r - phi_e;
      if (reconPhi_er < -M_PI) reconPhi_er += 2.*M_PI;
      if (reconPhi_er > M_PI) reconPhi_er -= 2.*M_PI;

      // Theta qr
      double trueCosTheta_qr = cos(trueTheta_r)*cos(theta_q) - sin(trueTheta_r)*sin(theta_q)*cos(truePhi_er);
      double reconCosTheta_qr = cos(reconTheta_r)*cos(theta_q) - sin(reconTheta_r)*sin(theta_q)*cos(reconPhi_er);
      
      // x prime
      trueXp = QSq/(2.*( nu*(mD-trueEr) + truePr*q*trueCosTheta_qr));
      reconXp = QSq/(2.*( nu*(mD-reconEr) + reconPr*q*reconCosTheta_qr));

      // W prime
      trueWp = sqrt(sq(mD) - QSq + sq(mP) + 2.*mD*(nu-trueEr) -2.* nu * trueEr + 2.*q*truePr*trueCosTheta_qr);
      reconWp = sqrt(sq(mD) - QSq + sq(mP) + 2.*mD*(nu-reconEr) -2.* nu * reconEr + 2.*q*reconPr*reconCosTheta_qr);

      // alpha_s
      trueAs = (trueEr - truePr*trueCosTheta_qr)/mN;
      reconAs = (reconEr - reconPr*reconCosTheta_qr)/mN;

      // Fill the tree
      outTree->Fill();
    }

  // Write the tree
  gFile=outFile;
  outTree->Write();

  // Clean up
  inFile->Close();
  outFile->Close();
  delete myRand;

  return 0;
}
