#include <iostream>
#include <cstdlib>
#include <cmath>

#include "TFile.h"
#include "TTree.h"
#include "TRandom3.h"
#include "TVectorT.h"

#include "constants.h"
#include "gen_tree.h"

using namespace std;

// Define some constants for BAND (positions in cm, time in ns, mass in GeV)
const double bandZ = -250.; // The z position of the face of BAND that is closest to the target (in cm)
const double bandZWidth=40.;
const double cScint = 15.;
TRandom3 *myRand;

inline double sq(double x){ return x*x;};
bool didItHit(double x, double y);
double getCLAS12_PRes(double theta, double p);
double eeEDep(double En, double bar_thick);

int main(int argc, char ** argv)
{
  if (argc != 6)
    {
      cerr << "Wrong number of arguments. Instead use\n\tdet_sim /path/to/input/file /path/to/output/file [bar Width (cm)]\n"
	   << "\t\t[PMT time res (ps)] [CLAS-12 res. factor (1=normal)]\n";
      exit(-1);
    }

  const double barWidth = atof(argv[3]);
  const double tResPMT = 1.E-3 * atof(argv[4]); // convert time res from ps to ns
  const double clasFac=atof(argv[5]); // Multiplies all CLAS-12 resolutions

  // Random number generator
  myRand = new TRandom3(0);

  // Get the files and trees
  TFile * inFile = new TFile(argv[1]);
  TFile * outFile = new TFile(argv[2],"RECREATE");
  TTree * inTree = (TTree*)inFile->Get("MCout");
  TTree * outTree = new TTree("ReconTree","Reconstructed Values");

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

  // Memory for the input branches
  Gen_Event * this_event = NULL;
  inTree->SetBranchAddress("event",&this_event);

  // Memory for the output tree
  double trueWp, trueXp, trueAs, trueQSq;
  double truePr, trueEr, trueZHit;
  double mom_e[3];
  double mom_r[3];

  double reconWp, reconXp, reconAs, reconQSq;
  double reconPr, reconEr, reconZHit;
  double reconMom_e[3];
  double reconMom_r[3];

  // Significant data
  outTree->Branch("trueWp",&trueWp,"trueWp/D");
  outTree->Branch("trueXp",&trueXp,"trueXp/D");
  outTree->Branch("trueAs",&trueAs,"trueAs/D");
  outTree->Branch("trueQSq",&trueQSq,"trueQSq/D");

  outTree->Branch("reconWp",&reconWp,"reconWp/D");
  outTree->Branch("reconXp",&reconXp,"reconXp/D");
  outTree->Branch("reconAs",&reconAs,"reconAs/D");
  outTree->Branch("reconQSq",&reconQSq,"reconQSq/D");

  // Modulus of momentum vectors and energy
  outTree->Branch("truePr",&truePr,"truePr/D");
  outTree->Branch("trueEr",&trueEr,"trueEr/D");
  outTree->Branch("reconPr",&reconPr,"reconPr/D");
  outTree->Branch("reconEr",&reconEr,"reconEr/D");

  // Momentum vector components
  outTree->Branch("truePe_x",&(mom_e[0]),"truePe_x/D");
  outTree->Branch("truePe_y",&(mom_e[1]),"truePe_y/D");
  outTree->Branch("truePe_z",&(mom_e[2]),"truePe_z/D");
  outTree->Branch("truePr_x",&(mom_r[0]),"truePr_x/D");
  outTree->Branch("truePr_y",&(mom_r[1]),"truePr_y/D");
  outTree->Branch("truePr_z",&(mom_r[2]),"truePr_z/D");

  outTree->Branch("reconPe_x",&(reconMom_e[0]),"reconPe_x/D");
  outTree->Branch("reconPe_y",&(reconMom_e[1]),"reconPe_y/D");
  outTree->Branch("reconPe_z",&(reconMom_e[2]),"reconPe_z/D");
  outTree->Branch("reconPr_x",&(reconMom_r[0]),"reconPr_x/D");
  outTree->Branch("reconPr_y",&(reconMom_r[1]),"reconPr_y/D");
  outTree->Branch("reconPr_z",&(reconMom_r[2]),"reconPr_z/D");

  // Hit position
  outTree->Branch("trueZHit",&trueZHit,"trueZHit/D");
  outTree->Branch("reconZHit",&reconZHit,"reconZHit/D");

  // Loop over the events
    for (int i=0 ; i<inTree->GetEntries() ; i++)
    {
      if (i%100000==0)
	cerr << "Working on event " << i << "\n";

      // Load the event
      inTree->GetEvent(i);

      // Copy over the true momentum vectors
      mom_e[0] = this_event->particles[0].momentum.X();
      mom_e[1] = this_event->particles[0].momentum.Y();
      mom_e[2] = this_event->particles[0].momentum.Z();

      mom_r[0] = this_event->particles[1].momentum.X();
      mom_r[1] = this_event->particles[1].momentum.Y();
      mom_r[2] = this_event->particles[1].momentum.Z();

      // Figure out where the recoiling neutron hits BAND
      trueZHit = bandZ - bandZWidth*myRand->Rndm();
      double trueXHit = mom_r[0]*trueZHit/mom_r[2];
      double trueYHit = mom_r[1]*trueZHit/mom_r[2];

      // Test if it hit the detector
      //if (!didItHit(trueXHit,trueYHit))
      //continue;

      // Energy and timing
      truePr = sqrt(sq(mom_r[0]) + sq(mom_r[1]) + sq(mom_r[2]));
      trueEr = sqrt(sq(truePr) + sq(mN));
      double trueT = sqrt(sq(trueXHit)+sq(trueYHit)+sq(trueZHit)) / (cAir * truePr / trueEr); // d / (c * beta)
      double trueEDepee = eeEDep(trueEr, barWidth);

      // Test if we were efficient for the neutron hit
      // This could in the future depend on neutron momentum
      //if (myRand->Rndm() > 0.3)
      //continue;

      // Flight time smearing is based on PMT resolution for mean time
      double reconT = myRand->Gaus(trueT,tResPMT/sqrt(2.));

      // X reconstruction is smeared based on PMT resolution for time difference
      double reconXHit = myRand->Gaus(trueXHit,cScint * tResPMT/sqrt(2.));

      // Y reconstruction assumes hit occurred at the center of the bar
      int bar = floor(trueYHit / barWidth);
      double reconYHit = barWidth* (((float)bar) + 0.5);

      // Z is also reconstructed assuming a hit at the center of the bar
      bar = floor((bandZ-trueZHit)/barWidth);
      reconZHit = bandZ - barWidth*(((float)bar) + 0.5);

      // Momentum
      double reconPath = sqrt(sq(reconXHit)+sq(reconYHit)+sq(reconZHit));
      double reconBeta = reconPath/(cAir*reconT);
      reconPr = mN / sqrt( sq(1./reconBeta) - 1.);
      reconEr = sqrt(sq(reconPr)+sq(mN));

      // True lepton quantities
      double truePe =sqrt(sq(mom_e[0]) + sq(mom_e[1]) + sq(mom_e[2]));
      double trueCosTheta_e = mom_e[2]/truePe;
      double trueTheta_e = acos(trueCosTheta_e);
      double truePhi_e = atan2(mom_e[1],mom_e[0]);
      trueQSq = 2.*E1*truePe*(1.-trueCosTheta_e);
      double trueNu = E1 - truePe;
      double trueX = trueQSq / (2.*mP*trueNu);
      double true_q = sqrt(trueQSq + sq(trueNu));
      double truePhi_q = truePhi_e+M_PI;
      if (truePhi_q > M_PI) truePhi_q -= 2.*M_PI;
      double trueTheta_q = acos((E1 - truePe*trueCosTheta_e)/true_q);

      // Lepton needs to be smeared
      double reconPe = myRand->Gaus(truePe,clasFac*getCLAS12_PRes(trueTheta_e,truePe));
      double reconTheta_e = myRand->Gaus(trueTheta_e,clasFac*0.001); // 1mrad smearing
      double reconPhi_e = myRand->Gaus(truePhi_e,clasFac*0.001/sin(trueTheta_e)); // 1mrad/sin(theta)
      if (reconPhi_e > M_PI) reconPhi_e -=2.*M_PI;
      if (reconPhi_e < -M_PI) reconPhi_e +=2.*M_PI;
      reconQSq = 2.*E1*reconPe*(1.-cos(reconTheta_e));
      double reconNu = E1 - reconPe;
      double reconX = reconQSq / (2.*mP*reconNu);
      double recon_q = sqrt(reconQSq + sq(reconNu));
      double reconPhi_q = reconPhi_e + M_PI;
      if (reconPhi_q > M_PI) reconPhi_q -= 2.*M_PI;
      double reconTheta_q = acos((E1 - reconPe*cos(reconTheta_e))/recon_q);
      reconMom_e[0] = reconPe * sin(reconTheta_e) * cos(reconPhi_e);
      reconMom_e[1] = reconPe * sin(reconTheta_e) * sin(reconPhi_e);
      reconMom_e[2] = reconPe * cos(reconTheta_e);

      // Recoil neutron needs to be smeared
      double trueCosTheta_r = mom_r[2]/truePr;
      double trueTheta_r = acos(trueCosTheta_r);
      double truePhi_r = atan2(mom_r[1],mom_r[0]);
      double reconCosTheta_r = reconZHit/reconPath;
      double reconTheta_r = acos(reconCosTheta_r);
      double reconPhi_r = atan2(reconYHit,reconXHit);
      reconMom_r[0] = reconPr * sin(reconTheta_r) * cos(reconPhi_r);
      reconMom_r[1] = reconPr * sin(reconTheta_r) * sin(reconPhi_r);
      reconMom_r[2] = reconPr * reconCosTheta_r;      
      
      // Let's take care of those pesky relative phi angles
      double truePhi_er = truePhi_r - truePhi_e;
      if (truePhi_er < -M_PI) truePhi_er += 2.*M_PI;
      if (truePhi_er > M_PI) truePhi_er -= 2.*M_PI;
      double reconPhi_er = reconPhi_r - reconPhi_e;
      if (reconPhi_er < -M_PI) reconPhi_er += 2.*M_PI;
      if (reconPhi_er > M_PI) reconPhi_er -= 2.*M_PI;

      // Theta qr
      double trueCosTheta_qr = cos(trueTheta_r)*cos(trueTheta_q) - sin(trueTheta_r)*sin(trueTheta_q)*cos(truePhi_er);
      double reconCosTheta_qr = cos(reconTheta_r)*cos(reconTheta_q) - sin(reconTheta_r)*sin(trueTheta_q)*cos(reconPhi_er);
      
      // x prime
      trueXp = trueQSq/(2.*( trueNu*(mD-trueEr) + truePr*true_q*trueCosTheta_qr));
      reconXp = reconQSq/(2.*( reconNu*(mD-reconEr) + reconPr*recon_q*reconCosTheta_qr));

      // W prime
      trueWp = sqrt(sq(mD) - trueQSq + sq(mP) + 2.*mD*(trueNu-trueEr) -2.* trueNu * trueEr + 2.*true_q*truePr*trueCosTheta_qr);
      reconWp = sqrt(sq(mD) - reconQSq + sq(mP) + 2.*mD*(reconNu-reconEr) -2.* reconNu * reconEr + 2.*recon_q*reconPr*reconCosTheta_qr);

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

bool didItHit(double x, double y)
{
  // This function will do the geometric test to see if the hit position is in the array.

  // Test if below bottom bar
  if (y<-50.)
    return false;

  // Test if outside radius
  double r=sqrt(sq(x)+sq(y));
  if (r > 140.)
    return false;

  // Test if in the beam pipe gap
  if ((fabs(y) < 20.) and (fabs(x)<50.))
    return false;

  return true;
}

double getCLAS12_PRes(double theta, double p)
{
  const double a=0.0549141;
  const double b=-0.00332126;
  const double c=-0.052163;

  return p*(a+b*sin(theta)+c*cos(theta));  
}

double eeEDep(double En, double bar_thick)
{
    // test if neutron deposits energy
    const double a = 0.576428;
    const double b = 0.177634;
    const double c = 0.164202;
    const double f = a*exp(-b*(En - mN)) + c;
    
    if (myRand->Rndm() > f*(bar_thick/5))
        return 0;

    // pick trueEDep
    double trueEDep = myRand->Rndm() * (En - mN);

    // translate to e-equivalent E
    double trueE_MeVee = 0.95*trueEDep-8.0*(1-exp(0.10*(pow(trueEDep,0.90))));

    return trueE_MeVee;
}
