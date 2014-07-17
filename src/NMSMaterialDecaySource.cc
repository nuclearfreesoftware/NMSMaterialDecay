/* Copyright (C) 2014, Moritz Kütt
 * 
 * This file is part of NMSMaterialDecay.
 * 
 * NMSMaterialDecay is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * NMSMaterialDecay is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with NMSMaterialDecay.  If not, see <http://www.gnu.org/licenses/>.
 * 
 */


#include "NMSMaterialDecaySource.hh"

NMSMaterialDecaySource::NMSMaterialDecaySource() {
  startSourceTimeDistribution = 0;
  endSourceTimeDistribution = 0;

  spontaneousFissionNeutron = true;
  spontaneousFissionGamma = true;
  alphaDecay = false;
  betaDecay = false;
  alphaN = false;
  sourceGenerator = new NMSMultipleDecaySource();
  posGenerator = sourceGenerator->GetAllPosDist();
  currentSourceMaterial = 0;

  activity = 1;
  materialIntensity = 1;
  activeVolume = 1 * cm * cm * cm;

  sourcesloaded = false;
}

NMSMaterialDecaySource::NMSMaterialDecaySource(G4Material* sourceMat) {
  sourcesloaded = false;

  startSourceTimeDistribution = 0;
  endSourceTimeDistribution = 0;

  spontaneousFissionNeutron = true;
  spontaneousFissionGamma = true;
  alphaDecay = true;
  betaDecay = true;
  alphaN = false;
  sourceGenerator = new NMSMultipleDecaySource();
  posGenerator = sourceGenerator->GetAllPosDist();
  posGenerator->SetPosDisType("Point");
  posGenerator->SetCentreCoords(G4ThreeVector(0,0,0));

  if(verboseLevel >= 2) {
    G4cout << "NMSMaterialDecaySource: Set source material to " << sourceMat << G4endl;
  }
  currentSourceMaterial = sourceMat;

  activity = -1;
  materialIntensity = -1;
  activeVolume = 1 * cm * cm * cm;

  sourcesloaded = false;

}

NMSMaterialDecaySource::~NMSMaterialDecaySource() {
  G4cout << "Activity: " << activity * second << G4endl;
  sourceGenerator->ClearAll();
  delete sourceGenerator;
}

void NMSMaterialDecaySource::SetSourceMaterial(G4Material* sourceMat) {
  currentSourceMaterial = sourceMat;
  sourcesloaded = false;
  // define isotopes as in G4Material
}


void NMSMaterialDecaySource::SetActiveVolume(G4double vol){
  activeVolume = vol;
  if(materialIntensity != -1.) {
    activity = materialIntensity * activeVolume;
  }
}

void NMSMaterialDecaySource::SetActivity(G4double act){
  activity = act;
}

void NMSMaterialDecaySource::SetSpontaneousFission(G4bool neutron = true, G4bool gamma = true) {
  if ( spontaneousFissionNeutron != neutron ) {
    spontaneousFissionNeutron = neutron;
    sourcesloaded = false;
  }
  if ( spontaneousFissionGamma != gamma ) {
    spontaneousFissionGamma = gamma;
    sourcesloaded = false;
  }
}

void NMSMaterialDecaySource::SetAlphaDecay(G4bool status) {
  if ( alphaDecay != status ) {
    alphaDecay = status;
    sourcesloaded = false;
  }
}

void NMSMaterialDecaySource::SetBetaDecay(G4bool status) {
  if ( betaDecay != status ) {
    betaDecay = status;
    sourcesloaded = false;
  }
}

void NMSMaterialDecaySource::SetAlphaNSource(G4bool status) {
  if ( alphaN != status) {
    alphaN = status;
    sourcesloaded = false;
  }
}

void NMSMaterialDecaySource::SetAlphaNFile(G4String filename){
  alphaNFilename = filename;
  sourcesloaded = false;
}

void NMSMaterialDecaySource::GeneratePrimaryVertex(G4Event* anEvent) {
  if(currentSourceMaterial == 0 && !alphaN) {
    G4cout << "NMSMaterialDecaySource - ERROR: No source material has been selected. Please use messenger commands to select a source material. If you haven't defined an appropriate source material in your detector construction, please do so as well" << G4endl;
    G4cout << "Aborting run!" << G4endl;
    G4RunManager::GetRunManager()-> AbortRun();
    return;
  }
  if(verboseLevel >= 2) {
    G4cout << "NMSMaterialDecaySource GeneratePrimaryVertex" << G4endl; // VERBOSEF
  }
  if(!sourcesloaded) {
    if(verboseLevel >= 2) {
      G4cout << "NMSMaterialDecaySource: Load Sources" << G4endl; // VERBOSEF
    }
    LoadSources();
  }
  G4double time = GetNextEventTime();
  //  G4AnalysisManager* anaman = G4AnalysisManager::Instance();
  //  anaman->FillH1(1, time * second);

  SetTime(time);
  sourceGenerator->GeneratePrimaryVertex(anEvent);
}

void NMSMaterialDecaySource::SetEventTimeLimits(G4double start, G4double end) {
  startSourceTimeDistribution = start;
  endSourceTimeDistribution = end;

  //fix check ranges
}

G4double NMSMaterialDecaySource::GetNextEventTime(){
  G4double time;
  G4double ran = G4UniformRand();

  if((startSourceTimeDistribution == 0) && (endSourceTimeDistribution == 0)) {
    time = - 1 / activity * std::log(ran) * second;
  }
  else {
    if(activity > pow(10., -15)) {  // workaround to avoid numerical difficulties
      G4double alpha = std::exp(- startSourceTimeDistribution / second * activity);
      G4double beta = std::exp(- endSourceTimeDistribution / second * activity);

      time = - 1 / activity * std::log(beta + ran * ( alpha - beta)) * second;
    }
    else {
      time = ran * (endSourceTimeDistribution / second- startSourceTimeDistribution / second) * second;
    }
  }
  if(verboseLevel >= 2) {
    G4cout << " time of next source event: " << time / second << G4endl;
  }

  return time;
}

void NMSMaterialDecaySource::SetTime(G4double time){
  sourceGenerator->SetParticleTime(time);
}

void NMSMaterialDecaySource::SetVerboseLevel(G4int verbose) {
  verboseLevel = verbose;
  sourceGenerator->SetVerboseLevel(verbose);
  posGenerator->SetVerboseLevel(verbose);
}

void NMSMaterialDecaySource::LoadSources() {
  // FIX: Check if possible
  // FIX: Load AlphaN!!

  G4int verboseLevel = 2;

  G4ParticleDefinition* ion = 0;
  G4IsotopeVector* ivec;
  G4double* relabvec;
  G4double elementatoms;
  G4int a; G4int z;
  G4double lifetime = -1; G4double halflife;
  G4double sfbranching = -1;  G4double alphabranching = -1;  G4double betabranching = -1;
  G4double tmpactivity; G4double totalactivity = 0;

  sourceGenerator->ClearAll();
  materialIntensity = 0;

  if(verboseLevel >= 1) {
    G4cout << "================================================================" << G4endl;
    G4cout << "LoadSources" << G4endl;
    if(verboseLevel >= 2) {
      G4cout << "================================================================" << G4endl;
      G4cout << "Source Material" << G4endl;
      if(currentSourceMaterial == 0) {
	G4cout << "Not defined!" << G4endl;
      }
      else {
	G4cout << currentSourceMaterial << G4endl;
      }
      G4cout << "================================================================" << G4endl;
    }
  }

  if(currentSourceMaterial != 0) {
    const G4ElementVector* evec = currentSourceMaterial->GetElementVector();
    // cycle through elments
    for(size_t i=0; i < evec->size(); i++) {
      elementatoms = currentSourceMaterial->GetVecNbOfAtomsPerVolume()[i];
      ivec = (*evec)[i]->GetIsotopeVector();
      // Geant Documenation: a vector of relative abundances referring to such isotopes (where relative abundance means the number of atoms per volume)
      relabvec = (*evec)[i]->GetRelativeAbundanceVector();

      // cycle through isotopes
      for(size_t j=0; j < ivec->size(); j++) {
	a = (*ivec)[j]->GetN();
	z = (*ivec)[j]->GetZ();
	ion = G4IonTable::GetIonTable()->GetIon(z,a,0);
	lifetime = ion->GetPDGLifeTime();
	halflife = lifetime*log(2);
	if(lifetime != -1) {
	  sfbranching = GetSFBranching(a, z);
	  alphabranching = GetAlphaBranching(a, z);
	  betabranching = GetBetaBranching(a, z);
	  //adjust alpha and betabranching - geants own libraries do not include sf values up to now
	  if(sfbranching != 0.) {
	    alphabranching *= (1 - sfbranching);
	    betabranching *= (1 - sfbranching);
	  }
	  totalactivity += relabvec[j] * elementatoms / lifetime;
	}
	if(verboseLevel >= 1) {
	  G4cout << z * 10000 + a * 10 << " " << relabvec[j] << " " << relabvec[j] * elementatoms * (cm3)<< " (per cm^3) tao=" << lifetime / second << " s t_1/2=" << halflife/second << " s" << G4endl;
	  if(lifetime != -1) {
	    G4cout.width(8);
	    G4cout << "";
	    G4cout.width(8);
	    G4cout << sfbranching;
	    G4cout.width(8);
	    G4cout << alphabranching;
	    G4cout.width(8);
	    G4cout << betabranching << G4endl;
	  }
	  else {
	    G4cout << "Stable isotope, will not be used as source." << G4endl;
	  }
	}
	// ??? what happens if no decay - fix
	if(lifetime != -1) {
	  if(spontaneousFissionNeutron || spontaneousFissionGamma) {
	    // N * branching * lambda = N * branching * tao
	    tmpactivity = relabvec[j] * elementatoms * sfbranching / lifetime;
	    materialIntensity += tmpactivity;
	    sourceGenerator->AddaSource(tmpactivity);
	    sourceGenerator->GetCurrentSource()->setIsotope(z * 10000 + a * 10);
	    if(spontaneousFissionNeutron && spontaneousFissionGamma) {
	      sourceGenerator->GetCurrentSource()->setDecayType(NMSDECAY_SF);
	    }
	    else if (spontaneousFissionNeutron) {
	      sourceGenerator->GetCurrentSource()->setDecayType(NMSDECAY_SF_N);
	    }
	    else if (spontaneousFissionGamma) {
	      sourceGenerator->GetCurrentSource()->setDecayType(NMSDECAY_SF_GAMMA);
	    }
	  }
	  if(alphaDecay) {
	    // N * branching * lambda = N * branching * tao
	    tmpactivity = relabvec[j] * elementatoms * alphabranching / lifetime;
	    materialIntensity += tmpactivity;
	    sourceGenerator->AddaSource(tmpactivity);
	    sourceGenerator->GetCurrentSource()->setIsotope(z * 10000 + a * 10);
	    sourceGenerator->GetCurrentSource()->setDecayType(NMSDECAY_ALPHA);
	  }
	  if(betaDecay) {
	  }
	}
      } // end isotopes
    } // end elements

  }
  else {
    if(alphaN) {
      sourceGenerator->AddaSource(1);
      sourceGenerator->GetCurrentSource()->setDecayType(NMSDECAY_ALPHA_N);
      sourceGenerator->GetCurrentSource()->setAlphaNFile(alphaNFilename);
    }
    else {
      G4cout << "Error: No Source Material defined - Can not load sources." << G4endl;
      return;
    }
  }
  sourcesloaded = true;
  activity = materialIntensity * activeVolume;
  if(verboseLevel >= 1) {
    G4cout << "The MaterialDecaySource represents an activity of: " << activity * second << "/s"<< G4endl;
    //FIX Level
    if(verboseLevel >= 1) {
      G4cout << "based on a material intensity of:                  " << materialIntensity * second * ( cm * cm * cm ) << "/(s * cm^3)"<< G4endl;
      G4cout << "and an active volume of:                           " << activeVolume / ( cm * cm * cm ) << " cm^3 " << G4endl;
      G4cout << "Out of a total activity of:                        " << activeVolume * totalactivity *second << "/s" << G4endl;
    }
    G4cout << "================================================================" << G4endl;
  }
}

G4double NMSMaterialDecaySource::GetSFBranching(G4int a, G4int z) {
  G4int zaid = z * 1000 + a;
  G4double br = 0;

  switch(zaid) {
  case 90232:
    br = 1.1e-11; // nndc (03.03.2014)
      break;
  case 92232:
    br = 3e-14; // nndc (03.03.2014)
    break;
  case 92233:
    br = 6e-13; // nndc (< 6e-11, also mg28 and ne24 decay, 03.03.2014);
    // br = 1.3e-12; // shultis / faw
    // janis lists no branching ratio
    break;
  case 92234:
    br = 1.6e-9; // nndc (03.03.2014)
    break;
  case 92235:
    br = 7.2e-11; // janis (jeff 3.1.1 and endf/b-vii.1)
    // br = 7.0e-11; // nndc (03.03.2014)
    // br = 2.0e-9; // shultis / faw
    break;
  case 92236:
    br = 9.4e-10; // nndc (03.03.2014)
    break;
  case 92238:
    br = 5.46e-7; // janis (jeff 3.1.1 and endf/b-vii.1)
    // br = 5.5e-7; // nndc (03.03.2014)
    // br = 5.4e-7; // shultis / faw
    break;
  case 93237:
    br = 2e-12; // nndc (<= 2e-12, 03.03.2014)
    // br = 2.1e-14; // shultis / faw
    // janis lists no branching ratio
    break;
  case 94236: // ++
    br = 8.2e-10; // janis (jeff 3.1.1 and endf/b-vii.1)
    // br = 1.9e-9;  // nndc (03.03.2014)
    // br = 8.1e-10; // shultis / faw
    break;
  case 94238:
    br = 1.86e-9; // janis (jeff 3.1.1 and endf/b-vii.1)
    // br = 1.9e-9; // nndc (03.03.2014)
    // br = 1.8e-9; // shultis / faw
    break;
  case 94239:
    br = 3.1e-12; // janis (jeff 3.1.1 and endf/b-vii.1)
    // 3.1 is also referenced in nuclear data sheets vol 98, 2002
    // br = 3e-12;  // nndc (03.03.2014)
    // br = 4.4e-12; // shultis / faw
    break;
  case 94240:
    br = 5.7e-8; // janis (jeff 3.1.1 and endf/b-vii.1)
    // 5.7 is also reference in nuclear data sheets 2008 (a = 240)
    // br = 5.6552e-8; // panda 1991
    // br = 5.7e-8; // nndc (03.03.2014)
    // br = 5e-8; // shultis / faw
    break;
  case 94241:
    br = 2e-16; //nndc (03.03.2014)
    // < 2.4e-16 referenced in nuclear data sheets 206 (2005), a = 241
    // br = 5.7e-15; // shultis / faw
    // br = 5.74e-15; // panda 1991/2007
    break;
  case 94242:
    br = 5.5e-6; // nndc (03.03.2014)
    // br = 5.5e-6; // shultis/faw
    // br = 5.5e-6; // janis (jeff 3.1.1 and endf/b-vii.1)
    // br = 5.4971e-6; // panda 1991/2007
    break;
  case 95241:
    br = 4e-12; // nndc (03.03.2014)
    // br = 4.1e-12; // shultis / faw
    break;
  case 96242: // ++
    br = 6.2e-8; // nndc (03.03.2014)
    break;
  case 96244: // ++
    br = 1.4e-6; // nndc (03.03.2014)
    break;
  case 96246: // ++
    br = 3e-4; // nndc (03.03.2014)
    break;
  case 96248: // ++
    br = 8.39e-2; // nndc (03.03.2014)
    break;
  case 97249:
    br = 4.7e-10; // nndc (03.03.2014)
    break;
  case 98250: // ++
    br = 8e-4; // nndc (03.03.2014)
    break;
  case 98246: // ++
    br = 2.4e-6; // nndc (03.03.2014)
    break;
  case 98252: // ++
    br = 3.092e-2; // janis (jeff 3.1.1 and endf/b-vii.1)
    // br = 3.092e-2; // referenced in Nuclear Data Sheets 106 (2005) 813–834 (A = 252)
    // br = 3.09e-2; // nndc (03.03.2014)
    // br = 3.09e-2; // shultis / faw
    break;
  case 98254: // ++
    br = 9.969e-1; // nndc (03.03.2014)
    break;
  case 98256: // ++
    br = 1.; // nndc (03.03.2014)
    break;
  case 100257: // ++
    br = 2.1e-3; // nndc (03.03.2014)
    break;
  case 102252: //++
    br = 2.93e-1; //nndc (03.03.2014)
    break;
  default:
    if(verboseLevel >= 1) {
      G4cout << "No spontaneous fission data available for isotope " << zaid << G4endl;
    }
    br = 0;
  }
  return br;
}

G4double NMSMaterialDecaySource::GetAlphaBranching(G4int a, G4int z) {
  G4ParticleDefinition* ion =  G4IonTable::GetIonTable()->GetIon(z,a,0);
  G4RadioactiveDecay* decay = new G4RadioactiveDecay();
  G4DecayTable* dtable = decay->GetDecayTable(ion);
  G4VDecayChannel* dc;

  G4double br = 0;
  for(G4int i = 0; i < dtable->entries(); i++) {
    dc = dtable->GetDecayChannel(i);
    if( G4AlphaDecayChannel* testalpha = dynamic_cast< G4AlphaDecayChannel* >( dc ) ) {
      br += dc->GetBR();
    }
  }

  return br;
}

G4double NMSMaterialDecaySource::GetBetaBranching(G4int a, G4int z) {
  G4ParticleDefinition* ion =  G4IonTable::GetIonTable()->GetIon(z,a,0);
  G4RadioactiveDecay* decay = new G4RadioactiveDecay();
  G4DecayTable* dtable = decay->GetDecayTable(ion);
  G4VDecayChannel* dc;

  G4double br = 0;
  for(G4int i = 0; i < dtable->entries(); i++) {
    dc = dtable->GetDecayChannel(i);
    if( G4BetaMinusDecayChannel* testchannel = dynamic_cast< G4BetaMinusDecayChannel* >( dc ) ) {
      br += dc->GetBR();
    }
  }

  return br;

}
