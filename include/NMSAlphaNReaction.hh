#ifndef NMSAlphaNReaction_h
#define NMSAlphaNReaction_h 1

#include "G4ThreeVector.hh"

struct NMSAlphaNReaction {
  G4ThreeVector position;
  G4ThreeVector alphaDirection;
  G4double energy;
  G4double time;
};

#endif
