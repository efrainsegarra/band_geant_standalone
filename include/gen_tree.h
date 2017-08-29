#ifndef __GEN_TREE_H__
#define __GEN_TREE_H__

#include <string>
#include <vector>
#include "TObject.h"
#include "TVector3.h"

class Gen_Particle : public TObject
{
 public:
  Gen_Particle();
  ~Gen_Particle();
  std::string type;
  TVector3 momentum;
  double t0;
  ClassDef(Gen_Particle,2);
};

class Gen_Event : public TObject
{
 public:
  Gen_Event();
  ~Gen_Event();
  std::vector<Gen_Particle> particles;
  ClassDef(Gen_Event,1);
};

#endif
