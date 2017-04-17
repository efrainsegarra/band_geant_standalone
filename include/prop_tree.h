#ifndef __PROP_TREE_H__
#define __PROP_TREE_H__

#include <vector>
#include "TObject.h"
#include "TVector3.h"

class BAND_Hit : public TObject
{
 public:
  BAND_Hit();
  ~BAND_Hit();
  TVector3 pos; // in mm
  double time; // in ns
  double E_dep; // in MeV
  int barNo; // indicating what bar
  
  ClassDef(BAND_Hit,5);
};

class BAND_Event : public TObject
{
 public:
  BAND_Event();
  ~BAND_Event();
  std::vector<BAND_Hit> hits;

  ClassDef(BAND_Event,1);
};

#endif
