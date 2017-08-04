#ifndef __RECON_TREE_H__
#define __RECON_TREE_H__

#include <string>
#include <vector>
#include "TObject.h"
#include "TVector3.h"

class Recon_Particle : public TObject
{
 public:
  Recon_Particle();
  ~Recon_Particle();
  std::string type;
  TVector3 momTrue;
  TVector3 momRecon;

  TVector3 pos; // in mm
  double time; // in ns
  double E_dep; // in MeV
  int barNo; // indicating what bar

  ClassDef(Recon_Particle,3);
};

class Recon_Event : public TObject
{
 public:
  Recon_Event();
  ~Recon_Event();
  std::vector<Recon_Particle> particles;

  ClassDef(Recon_Event,3);
};

#endif
