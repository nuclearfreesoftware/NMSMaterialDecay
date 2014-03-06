#include "NMSPosDistribution.hh"

NMSPosDistribution::NMSPosDistribution(){

}

NMSPosDistribution::~NMSPosDistribution(){
  posVector.clear();
}

void NMSPosDistribution::ClearAll(){
  posVector.clear();
}

void NMSPosDistribution::AddaPosDist(G4SPSPosDistribution* posDist){
  posVector.push_back(posDist);

  posDist->SetPosDisType(SourcePosType);
  posDist->SetPosDisShape(Shape);
  posDist->SetCentreCoords(CentreCoords);
  posDist->SetRadius(Radius);
  posDist->SetRadius0(Radius0);
  posDist->SetHalfX(halfx);
  posDist->SetHalfY(halfy);
  posDist->SetHalfZ(halfz);

  G4cout << "Volume " << VolName << G4endl;

  if(Confine) {

    G4cout << "Confine to Volume " << VolName << G4endl;

    posDist->ConfineSourceToVolume(VolName);
  }

  // Set all Settings for posDist
}

void NMSPosDistribution::DeleteaPosDist(G4int idx){
  size_t sizeIdx = size_t (idx);
  if(sizeIdx <= posVector.size()){
    posVector.erase(posVector.begin() + idx);
  }
  else {
    G4cout << " source index is invalid " << G4endl;
    G4cout << "    it shall be <= " << posVector.size() << G4endl;
  }
}

void NMSPosDistribution::SetPosDisType(G4String type){
  SourcePosType = type;
  for(size_t i = 0; i < posVector.size(); i++) {
    posVector[i]->SetPosDisType(type);
  }
}

void NMSPosDistribution::SetPosDisShape(G4String shape){
  Shape = shape;
  for(size_t i = 0; i < posVector.size(); i++) {
    posVector[i]->SetPosDisShape(shape);
  }
}


void NMSPosDistribution::ConfineSourceToVolume(G4String volume){
  Confine = true;
  VolName = volume;
  for(size_t i = 0; i < posVector.size(); i++) {
    posVector[i]->ConfineSourceToVolume(volume);
  }
}

void NMSPosDistribution::SetCentreCoords(G4ThreeVector ccoords){
  CentreCoords = ccoords;
  for(size_t i = 0; i < posVector.size(); i++) {
    posVector[i]->SetCentreCoords(ccoords);
  }

}

// FIX Add PosRot

// FIX implement
void NMSPosDistribution::SetPosRot2(G4ThreeVector posrot){

}

void NMSPosDistribution::SetHalfX(G4double hx){
  halfx = hx;
  for(size_t i = 0; i < posVector.size(); i++) {
    posVector[i]->SetHalfX(hx);
  }
}

void NMSPosDistribution::SetHalfY(G4double hy){
  halfy = hy;
  for(size_t i = 0; i < posVector.size(); i++) {
    posVector[i]->SetHalfY(hy);
  }

}

void NMSPosDistribution::SetHalfZ(G4double hz){
  halfz = hz;
  for(size_t i = 0; i < posVector.size(); i++) {
    posVector[i]->SetHalfZ(hz);
  }
}

void NMSPosDistribution::SetRadius(G4double r){
  Radius = r;
  for(size_t i = 0; i < posVector.size(); i++) {
    posVector[i]->SetRadius(r);
  }

}

void NMSPosDistribution::SetRadius0(G4double r0){
  Radius0 = r0;
  for(size_t i = 0; i < posVector.size(); i++) {
    posVector[i]->SetRadius0(r0);
  }
}
