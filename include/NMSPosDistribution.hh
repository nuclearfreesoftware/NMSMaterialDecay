//******************************************************************************
// NMSPosDistribution.hh
//
// Wrapper Class for all Position Distributions - allows for easier adjustment
//
// This file has been inspired by G4SPSPosDistribution.hh
//******************************************************************************

#ifndef NMSPosDistribution_h
#define NMSPosDistribution_h 1

#include "G4SPSPosDistribution.hh"
#include <vector>

class NMSPosDistribution
{
public:
  NMSPosDistribution();
  ~NMSPosDistribution();

  void ClearAll();
  void AddaPosDist(G4SPSPosDistribution* posDist);
  void DeleteaPosDist(G4int);

  void SetPosDisType(G4String); // Point, Plane, Surface, Volume
  inline G4String GetPosDisType() { return SourcePosType; };
  void SetPosDisShape(G4String);
  inline G4String GetPosDisShape() { return Shape; };

  void ConfineSourceToVolume(G4String);
  void SetCentreCoords(G4ThreeVector);
  inline G4ThreeVector GetCentreCoords() { return CentreCoords; } ;
  void SetPosRot2(G4ThreeVector);
  void SetHalfX(G4double);
  inline G4double GetHalfX() { return halfx; } ;
  void SetHalfY(G4double);
  inline G4double GetHalfY()  { return halfy; } ;
  void SetHalfZ(G4double);
  inline G4double GetHalfZ()  { return halfz; } ;
  void SetRadius(G4double);
  inline G4double GetRadius()  { return Radius; };
  void SetRadius0(G4double);
  inline G4double GetRadius0() { return Radius0; };

private:
  std::vector<G4SPSPosDistribution*> posVector;

  G4String SourcePosType; //Point,Plane,Surface,Volume
  G4String Shape; //Circle,Square,Rectangle etc..

  G4double halfx, halfy, halfz; //half lengths
  G4double Radius; //Radius for circles or spheres
  G4double Radius0; // The inner radius of an annulus

  G4ThreeVector CentreCoords; // Coords of centre of input shape

  G4bool Confine; //If true confines source distribution to VolName
  G4String VolName;

};


#endif
