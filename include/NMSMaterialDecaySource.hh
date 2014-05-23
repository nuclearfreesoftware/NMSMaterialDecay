//******************************************************************************
// NMSMaterialDecaySource.hh
//
// Wrapper Class for G4GeneralParticleSource
//******************************************************************************

#ifndef NMSMaterialDecaySource_h
#define NMSMaterialDecaySource_h 1

#include "G4GeneralParticleSource.hh"
#include "G4SPSPosDistribution.hh"
#include "G4Material.hh"
#include "Randomize.hh"
#include "G4SystemOfUnits.hh"

#include "G4AlphaDecayChannel.hh"
#include "G4BetaMinusDecayChannel.hh"

#include "G4RunManager.hh"

#include "NMSSingleDecaySource.hh"
#include "NMSMultipleDecaySource.hh"
#include "NMSPosDistribution.hh"
#include "NMSAlphaNSet.hh"

class NMSMaterialDecaySource
{
public:
  NMSMaterialDecaySource();
  NMSMaterialDecaySource(G4Material* sourceMat);
  ~NMSMaterialDecaySource();

  G4Material* GetSourceMaterial() {return currentSourceMaterial; };
  void SetSourceMaterial(G4Material* sourceMat);

  void SetActivity(G4double act);
  void SetActiveVolume(G4double vol);

  NMSPosDistribution* GetPosDist() { return posGenerator; };

  // Set / Unset particles produced
  G4bool GetSpontaneousFissionNeutron() {return spontaneousFissionNeutron; };
  G4bool GetSpontaneousFissionGamma() {return spontaneousFissionGamma;}
  G4bool GetAlphaDecay() {return alphaDecay; };
  G4bool GetBetaDecay() {return betaDecay; };
  G4bool GetAlphaN() {return alphaN; };

  void SetSpontaneousFission(G4bool neutron, G4bool gamma);
  void SetAlphaDecay(G4bool status);
  void SetBetaDecay(G4bool status);

  void UnsetAlphaN();

  void SetAlphaNSet(NMSAlphaNSet* startpoints);

  void GeneratePrimaryVertex(G4Event*);

  void SetEventTimeLimits(G4double start, G4double end);
  G4double GetStartTime() {return startSourceTimeDistribution; };
  G4double GetEndTime() {return endSourceTimeDistribution; };
  G4double GetRuntime() {return endSourceTimeDistribution - startSourceTimeDistribution; };

  void SetVerboseLevel(G4int i);

  void LoadSources();

private:
  G4bool spontaneousFissionNeutron;
  G4bool spontaneousFissionGamma;
  G4bool alphaDecay;
  G4bool betaDecay;
  G4bool alphaN;

  G4bool sourcesloaded;

  G4Material* currentSourceMaterial;

  G4double activity;
  G4double activeVolume;

  G4double materialIntensity; // Intensity per ???

  NMSMultipleDecaySource* sourceGenerator;
  NMSPosDistribution* posGenerator;

  NMSAlphaNSet* alphaNPoints;

  G4double startSourceTimeDistribution;
  G4double endSourceTimeDistribution;

  G4int verboseLevel;


private:
  // Timing issues
  G4double GetNextEventTime();
  void SetTime(G4double);

  G4double GetSFBranching(G4int a, G4int z);
  G4double GetAlphaBranching(G4int a, G4int z);
  G4double GetBetaBranching(G4int a, G4int z);
};


#endif
