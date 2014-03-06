#ifndef NMSSingleDecaySource_h
#define NMSSingleDecaySource_h 1

#include "G4PrimaryVertex.hh"
#include "G4Event.hh"
#include "G4DynamicParticle.hh"
#include "G4Alpha.hh"
#include "G4Neutron.hh"
#include "G4Gamma.hh"

#include "G4DataVector.hh"
#include "Randomize.hh"

#include "G4IonTable.hh"
#include "G4RadioactiveDecay.hh"
#include "G4VDecayChannel.hh"
#include "G4DecayProducts.hh"
#include "G4AlphaDecayChannel.hh"

#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"

#include "G4SPSPosDistribution.hh"
#include "G4SPSAngDistribution.hh"

#include "G4RNGWrapper.hh"
#include "Fission.hh"

#include "G4SingleParticleSource.hh"

#define NMSDECAY_SF 1
#define NMSDECAY_ALPHA 2
#define NMSDECAY_BETA 3
#define NMSDECAY_SF_N 11
#define NMSDECAY_SF_GAMMA 12
#define NMSDECAY_ALPHA_N 21

class NMSSingleDecaySource : public virtual G4SingleParticleSource
{
public:
  NMSSingleDecaySource();
  NMSSingleDecaySource(G4int iso, G4int type);
  ~NMSSingleDecaySource();

  void setIsotope(G4int iso);
  void setDecayType(G4int type);

  void GenerateAlphaPrimaryVertex(G4Event* anEvent, G4PrimaryVertex* vertex);
  void GenerateSFPrimaryVertex(G4Event *anEvent,  G4PrimaryVertex* vertex);
  void GenerateBetaPrimaryVertex(G4Event *anEvent,  G4PrimaryVertex* vertex);
  void GeneratePrimaryVertex(G4Event* anEvent);

  void SetVerboseLevel(G4int i);

private:
  void LoadDecay();

private:
  G4int DecayIsotope;
  G4int DecayType;

  //G4SPSPosDistribution* posGenerator;
  //  G4SPSAngDistribution* angGenerator;

  G4int verboseLevel;

  std::vector<G4double> energyList;
  std::vector<G4double> branchingList;

  G4bool decayloaded;
};

#endif
