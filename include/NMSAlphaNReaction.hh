#ifndef NMSAlphaNReaction_h
#define NMSAlphaNReaction_h 1

#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"

class NMSAlphaNReaction {
public:
  NMSAlphaNReaction();
  NMSAlphaNReaction(G4ThreeVector, G4ThreeVector, G4double);
  ~NMSAlphaNReaction();  

  inline G4ThreeVector GetPosition() {return position; };
  inline G4ThreeVector GetAlphaDirection() {return alphaDirection; };
  inline G4double GetEnergy() {return energy; };

  inline void SetPosition(G4ThreeVector pos) { position = pos; };
  inline void SetAlphaDirection(G4ThreeVector adir) { alphaDirection = adir; };
  inline void SetEnergy(G4double e) { energy = e; };
  // fix operators = etc.

private:
  G4ThreeVector position;
  G4ThreeVector alphaDirection;
  G4double energy;

};

#endif
