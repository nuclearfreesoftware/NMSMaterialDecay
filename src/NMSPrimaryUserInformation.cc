#include "NMSPrimaryUserInformation.hh"

NMSPrimaryUserInformation::NMSPrimaryUserInformation() : G4VUserPrimaryParticleInformation() {
  
}

NMSPrimaryUserInformation::~NMSPrimaryUserInformation() {
  
}

void NMSPrimaryUserInformation::Print() const {
  
}

NMSOrigin NMSPrimaryUserInformation::GetOrigin() {
  return ori;
}

void NMSPrimaryUserInformation::SetOrigin(NMSOrigin newor) {
  ori = newor;
}

