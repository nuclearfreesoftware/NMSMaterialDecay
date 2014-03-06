#include "NMSMultipleDecaySource.hh"

NMSMultipleDecaySource::NMSMultipleDecaySource() {
  verboseLevel = 0;
  posGenerator = new NMSPosDistribution();
  sourceVector.clear();
  sourceIntensity.clear();
  sourceProbability.clear();

  currentSourceIdx = -1;
  currentSource = 0;
  normalised = false;
}

NMSMultipleDecaySource::NMSMultipleDecaySource(NMSSingleDecaySource* src, G4double strength) {
  verboseLevel = 0;

  posGenerator = new NMSPosDistribution();
  sourceVector.clear();
  sourceIntensity.clear();
  sourceProbability.clear();

  currentSourceIdx = -1;
  currentSource = 0;

  sourceVector.push_back(src);
  sourceIntensity.push_back(strength);
  currentSource = src;
  currentSourceIdx = G4int(sourceVector.size() - 1);

  normalised = false;
  IntensityNormalization();
}

NMSMultipleDecaySource::~NMSMultipleDecaySource(){
  delete posGenerator;
}

void NMSMultipleDecaySource::IntensityNormalization() {
  G4double total  = 0.;
  size_t i;

  if(verboseLevel >= 2) {
    G4cout << "================================================================" << G4endl;
    G4cout << "Intensity normalization" << G4endl;
  }

  if(verboseLevel >= 2) {
    G4cout << " Number of sources: " << sourceIntensity.size() << G4endl;
  }

  for (i = 0; i < sourceIntensity.size(); i++) {
    total += sourceIntensity[i] ;
  }
  sourceProbability.clear();
  sourceProbability.push_back(sourceIntensity[0] / total);
  for (i = 1; i < sourceIntensity.size(); i++) {
    sourceProbability.push_back(sourceIntensity[i] / total + sourceProbability[i-1]);
  }

  if(verboseLevel >= 2) {
    G4cout << " Source Probability Table" << G4endl;
    for (i = 0; i < sourceProbability.size(); i++) {
      G4cout << sourceProbability[i] << G4endl;
    }
  }

  if(verboseLevel >= 2) {
    G4cout << "================================================================" << G4endl;
  }

  normalised = true;
}

void NMSMultipleDecaySource::GeneratePrimaryVertex(G4Event* anEvent) {
  if(verboseLevel >= 2) {
    G4cout << "================================================================" << G4endl;
    G4cout << "*** NMSMultipleDecaySource GeneratePrimaryVertex" << G4endl;
  }

  size_t i = 0;

  if(!normalised) {
    IntensityNormalization();
  }

  if(verboseLevel >= 2) {
    G4cout << " Source Probability Table" << G4endl;
    for (i = 0; i < sourceProbability.size(); i++) {
      G4cout <<  "   " << sourceProbability[i] << G4endl;
    }
  }

  if(verboseLevel >= 2) {
    G4cout << " Source Intensity Table" << G4endl;
    for (i = 0; i < sourceIntensity.size(); i++) {
      G4cout << "   " << sourceIntensity[i] << G4endl;
    }
  }

  G4double ran = G4UniformRand();

  i=0;
  while(ran > sourceProbability[i]) {
    i++;
  }

  if(verboseLevel >= 2) {
    G4cout << " Selected Source " << i << G4endl;
    G4cout << "================================================================" << G4endl;
  }

  SetCurrentSourceto(i);

  currentSource->GeneratePrimaryVertex(anEvent);
}

void NMSMultipleDecaySource::AddaSource(G4double strength){
  currentSource = new NMSSingleDecaySource();
  sourceVector.push_back(currentSource);
  sourceIntensity.push_back(strength);
  currentSourceIdx = G4int(sourceVector.size() - 1);
  currentSource->SetVerboseLevel(verboseLevel);
  posGenerator->AddaPosDist(currentSource->GetPosDist());

  normalised = false;
  //  IntensityNormalization();
}

void NMSMultipleDecaySource::DeleteaSource(G4int idx) {
  size_t id = size_t (idx);

  if( id <= sourceIntensity.size() ) {
    sourceVector.erase(sourceVector.begin() + idx);
    sourceIntensity.erase(sourceIntensity.begin() + idx);

    if ( currentSourceIdx == idx ) {
      if (sourceVector.size() > 0) {
	currentSource = sourceVector[0];
	currentSourceIdx = 1;
      }
      else {
	currentSource = 0;
	currentSourceIdx = -1;
      }
    }
  }
  else {
    G4cout << " source index is invalid " << G4endl;
    G4cout << "    it shall be <= " << sourceIntensity.size() << G4endl;
  }

  posGenerator->DeleteaPosDist(idx);

  normalised = false;
  //  IntensityNormalization();
}

void NMSMultipleDecaySource::ClearAll() {
  currentSourceIdx = 0;
  currentSource = 0;
  sourceVector.clear();
  sourceIntensity.clear();
  sourceProbability.clear();

  posGenerator->ClearAll();

  normalised = false;
}

void NMSMultipleDecaySource::SetCurrentSourceto(G4int idx) {
  size_t id = size_t (idx);

  if ( id < sourceVector.size() ) {
    currentSourceIdx = idx;
    currentSource = sourceVector[id];
  }
  else {
    G4cout << " source index is invalid " << G4endl;
    G4cout << "    it shall be < " << sourceIntensity.size() << G4endl;
  }
}

void NMSMultipleDecaySource::SetCurrentSourceIntensity(G4double strength) {
  sourceIntensity[currentSourceIdx] = strength;

  normalised = false;
  //  IntensityNormalization();
}

void NMSMultipleDecaySource::SetVerboseLevel(G4int verbosity) {
  verboseLevel = verbosity;

  for(size_t i = 0; i < sourceVector.size(); i++) {
    sourceVector[i]->SetVerboseLevel(verbosity);
  }
}

void NMSMultipleDecaySource::SetParticleTime(G4double time) {
  for(size_t i = 0; i < sourceVector.size(); i++) {
    sourceVector[i]->SetParticleTime(time);
  }
}
