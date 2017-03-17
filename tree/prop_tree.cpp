#include "prop_tree.h"

BAND_Hit::BAND_Hit()
{}

BAND_Hit::~BAND_Hit()
{}

ClassImp(BAND_Hit);

BAND_Event::BAND_Event()
{
  hits.clear();
}

BAND_Event::~BAND_Event()
{}

ClassImp(BAND_Event);
