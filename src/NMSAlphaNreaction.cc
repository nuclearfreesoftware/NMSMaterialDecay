#include "NMSAlphaNReaction.hh"

NMSAlphaNReaction::NMSAlphaNReaction() {
  position = G4ThreeVector();
  alphaDirection = G4ThreeVector(1, 0, 0);
  energy = 1 * MeV;
}

NMSAlphaNReaction::NMSAlphaNReaction(G4ThreeVector pos, G4ThreeVector adir, G4double e) {
  position = pos;
  alphaDirection = adir;
  energy = e;
}

NMSAlphaNReaction::~NMSAlphaNReaction() {
  
}  
