#ifndef NMSPrimaryUserInformation_h
#define NMSPrimaryUserInformation_h 1

#include "G4VUserPrimaryParticleInformation.hh"

enum NMSOrigin { ORIGIN_ALPHA_N, ORIGIN_SF, ORIGIN_ALPHA };

class NMSPrimaryUserInformation : public G4VUserPrimaryParticleInformation {
public:
  NMSPrimaryUserInformation();
  ~NMSPrimaryUserInformation();

  void Print() const;

  NMSOrigin GetOrigin();
  void SetOrigin(NMSOrigin newor);

private:
  NMSOrigin ori;

};

#endif /* NMSPrimaryUserInformation_h */
