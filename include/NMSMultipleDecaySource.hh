#ifndef NMSMultipleDecaySource_h
#define NMSMultipleDecaySource_h 1

#include "G4Event.hh"
#include "Randomize.hh"

#include "NMSSingleDecaySource.hh"
#include "NMSPosDistribution.hh"

class NMSMultipleDecaySource : public G4VPrimaryGenerator
{
public:
  NMSMultipleDecaySource();
  NMSMultipleDecaySource(NMSSingleDecaySource* src, G4double strength);

  ~NMSMultipleDecaySource();

  G4int GetNumberofSource() { return G4int(sourceVector.size()); };

  void IntensityNormalization();
  void GeneratePrimaryVertex(G4Event* anEvent);

  void AddaSource (G4double strength);
  void DeleteaSource(G4int);
  void ClearAll();

  void SetCurrentSourceto(G4int) ;
  void SetCurrentSourceIntensity(G4double);

  NMSSingleDecaySource* GetCurrentSource() {return currentSource;};
  G4int GetCurrentSourceIndex() { return currentSourceIdx; };
  G4double GetCurrentSourceIntensity() { return sourceIntensity[currentSourceIdx]; };

  void SetParticleTime(G4double time);

  // Wrapper Class for all Position Distributions - allows for easier adjustment
  inline NMSPosDistribution* GetAllPosDist() {return posGenerator; };

  void SetVerboseLevel(G4int i);

private:
  NMSSingleDecaySource* currentSource;
  G4int currentSourceIdx;

  std::vector <NMSSingleDecaySource*> sourceVector;
  std::vector <G4double> sourceIntensity;
  std::vector <G4double> sourceProbability;

  G4bool normalised;

  NMSPosDistribution* posGenerator;

  G4int verboseLevel;
};
#endif
