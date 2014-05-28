#include "NMSSingleDecaySource.hh"

NMSSingleDecaySource::NMSSingleDecaySource() : G4SingleParticleSource(){
  /*  posGenerator = new G4SPSPosDistribution();
  posGenerator->SetBiasRndm(biasRndm);
  angGenerator = new G4SPSAngDistribution();
  angGenerator->SetPosDistribution(posGenerator);
  angGenerator->SetBiasRndm(biasRndm);
  */

  verboseLevel = 0;
  alphanset = 0;

  DecayType = NMSDECAY_ALPHA; // Spontaenous Fission n + gamma

  setIsotope(982520); // Cf-255
  //  setIsotope(942410); // Cf-241
  decayloaded = false;

  // set random number generator
  G4RNGWrapper<CLHEP::HepRandomEngine>::set(CLHEP::HepRandom::getTheEngine(),&CLHEP::HepRandomEngine::flat);
  setrngd_(G4RNGWrapper<CLHEP::HepRandomEngine>::rng);
}

NMSSingleDecaySource::NMSSingleDecaySource(G4int iso, G4int type) : G4SingleParticleSource() {
  verboseLevel = 0;
  alphanset = 0;

  //check if we can do isotope + type
  DecayType = type;
  setIsotope(iso);
  decayloaded = false;

  // set random number generator
  G4RNGWrapper<CLHEP::HepRandomEngine>::set(CLHEP::HepRandom::getTheEngine(),&CLHEP::HepRandomEngine::flat);
  setrngd_(G4RNGWrapper<CLHEP::HepRandomEngine>::rng);
}

NMSSingleDecaySource::~NMSSingleDecaySource(){

}

void NMSSingleDecaySource::setIsotope(G4int iso){
  if(verboseLevel >= 3) {
    G4cout << "NMSSingleDecaySource::setIsotope" << G4endl;
  }

  //check CheckIsotope(iso, DecayType);


  //check if we can do isotope + type
  DecayIsotope = iso;
  decayloaded = false;
}

void NMSSingleDecaySource::LoadDecay() {
  G4int A;
  G4int Z;
  G4double energy;
  G4ParticleDefinition* ion = 0;
  G4RadioactiveDecay* decay = new G4RadioactiveDecay();
  G4DecayTable* dtable;
  G4VDecayChannel* dc;
  G4DecayProducts* products;

  energyList.clear();
  branchingList.clear();

  Z = DecayIsotope / 10000;
  A = (DecayIsotope - Z * 10000) / 10;
  if (verboseLevel >= 1) {
    G4cout << "================================================================" << G4endl;
  }

  if(DecayType == NMSDECAY_ALPHA) {
    ion =  G4IonTable::GetIonTable()->GetIon(Z,A,0);
    if(verboseLevel >= 1) {
      G4cout << "Loading Alpha Decay for " << ion->GetParticleName() << G4endl;
    }
    dtable = decay->GetDecayTable(ion);
    for(G4int i = 0; i < dtable->entries(); i++) {
      dc = dtable->GetDecayChannel(i);
      if( G4AlphaDecayChannel* testalpha = dynamic_cast< G4AlphaDecayChannel* >( dc ) ) {
	products = dc->DecayIt(0);
	for(G4int j = 0; j < products->entries(); j++) {
	  if((*products)[j]->GetParticleDefinition()->GetParticleName() == "alpha") {
	    energy = (*products)[j]->GetKineticEnergy();
	    energyList.push_back(energy);
	    branchingList.push_back(dc->GetBR());
	    if(verboseLevel >= 1) {
	      G4cout << " Decay energy: " << energy / MeV << " MeV (Branching Ratio: " << dc->GetBR() << ")" << G4endl;
	    }
	  }
	  else {
	    //FIX AlphaChannel without alpha?
	  }
	}

      }
    }
  }

  //normalization + convert to probability
  G4double total = 0.;
  std::vector<G4double> tempList(branchingList.size());
  for (size_t i = 0; i < branchingList.size(); i++) {
    total += branchingList[i];
  }

  branchingList[0] = branchingList[0] / total;

  if (verboseLevel >= 1) {
    G4cout << " Branching Ratios are normalized internally: " << G4endl;
  }

  for (size_t i = 1 ;  i < branchingList.size(); i++) {
    if (verboseLevel >= 1) {
      G4cout << " " << branchingList[i] / total << ", summed: ";
    }
    branchingList[i] = branchingList[i] / total + branchingList[i-1];
    if (verboseLevel >= 1) {
      G4cout << " " << branchingList[i] << G4endl;
    }
  }

  if (verboseLevel >= 1) {
    G4cout << "================================================================" << G4endl;
  }

  decayloaded = true;
}

void NMSSingleDecaySource::setDecayType(G4int type){
  //check if we can do isotope + type
  DecayType = type;
}

void NMSSingleDecaySource::setAlphaNFile(G4String filename) {
  if(alphanset != 0) {
    delete alphanset;
  }
  alphanset = new NMSAlphaNSet;
  alphanset->loadFromFile(filename);
}

void NMSSingleDecaySource::GenerateSFPrimaryVertex(G4Event *anEvent,  G4PrimaryVertex* vertex) {

  G4ParticleDefinition* neutron_definition = G4Neutron::Neutron();
  G4ParticleDefinition* photon_definition = G4Gamma::Gamma();
  G4double mom, momx, momy, momz, eng;
  G4int nPrompt, gPrompt;

  G4int isotope = DecayIsotope / 10;
  G4double time = 0.;

  G4DynamicParticle* it;

  // zero polarisation for particles (as in example of llnl fission)
  G4ThreeVector polar;

  if(verboseLevel >= 1) {
    G4cout << "Spontaneous Fission!" << G4endl;
  }

  genspfissevt_(&isotope, &time);

  // this can be carried out without check for decay type
  nPrompt = getnnu_();
  gPrompt = getpnu_();

  if(nPrompt == -1) {
    G4cout << "================================================================" << G4endl;
    G4cout << "Error: NMSSingleDecaySource" << G4endl;
    G4cout << "Isotope is not available in Spontaneous Fission library: " << isotope << G4endl;
    exit(0);
  }

  if(verboseLevel >= 2)
    G4cout << "Creating primaries and assigning to vertex" << G4endl;

  // neutrons
  if(DecayType == NMSDECAY_SF || DecayType == NMSDECAY_SF_N) {
    for(G4int i=0; i<nPrompt; i++)
      {
	it = new G4DynamicParticle();
	it->SetDefinition(neutron_definition);
	eng = getneng_(&i)*MeV;
	it->SetKineticEnergy(eng);
	mom = it->GetTotalMomentum();

	momx = mom*getndircosu_(&i);
	momy = mom*getndircosv_(&i);
	momz = mom*getndircosw_(&i);

	G4PrimaryParticle* particle = new G4PrimaryParticle(
							    neutron_definition,
							    momx, momy, momz,
							    eng);
	particle->SetMass(neutron_definition->GetPDGMass());
	particle->SetCharge(neutron_definition->GetPDGCharge());
	particle->SetPolarization(polar.x(), polar.y(), polar.z());

	if(verboseLevel >= 2){
	  G4cout << "Particle name: "<<particle->GetG4code()->GetParticleName() << G4endl;
	  G4cout << "     Momentum: "<<particle->GetMomentum() << G4endl;
	  G4cout << "     Position: "<<vertex->GetPosition() << G4endl;
	}
	vertex->SetPrimary(particle);
      }

  }

  // photons
  if(DecayType == NMSDECAY_SF || DecayType == NMSDECAY_SF_GAMMA) {
    for(G4int i=0; i<gPrompt; i++) {
      it = new G4DynamicParticle();
      it->SetDefinition(photon_definition);
      eng = getpeng_(&i)*MeV;
      it->SetKineticEnergy(eng);
      mom = it->GetTotalMomentum();

      momx = mom*getpdircosu_(&i);
      momy = mom*getpdircosv_(&i);
      momz = mom*getpdircosw_(&i);

      G4PrimaryParticle* particle = new G4PrimaryParticle(
							  photon_definition,
							  momx, momy, momz,
							  eng);
      particle->SetMass(photon_definition->GetPDGMass());
      particle->SetCharge(photon_definition->GetPDGCharge());
      particle->SetPolarization(polar.x(), polar.y(), polar.z());

      if(verboseLevel >= 2){
	G4cout << "Particle name: "<<particle->GetG4code()->GetParticleName() << G4endl;
	G4cout << "     Momentum: "<<particle->GetMomentum() << G4endl;
	G4cout << "     Position: "<<vertex->GetPosition() << G4endl;
      }
      vertex->SetPrimary(particle);

    }
  }

  anEvent->AddPrimaryVertex( vertex );

}


void NMSSingleDecaySource::GenerateAlphaPrimaryVertex(G4Event *anEvent, G4PrimaryVertex* vertex) {
  G4ParticleDefinition* alpha_definition = G4Alpha::AlphaDefinition();
  G4DynamicParticle* alpha;

  if(!decayloaded) {
    LoadDecay();
  }

  G4double ran = G4UniformRand();
  G4double decayenergy = 0;
  G4int i=0;
  while(ran > branchingList[i])
    {
      i++;
    }
  decayenergy = energyList[i];

  //Calculate Momentum
  alpha = new G4DynamicParticle();
  alpha->SetDefinition(alpha_definition);
  alpha->SetKineticEnergy(decayenergy);
  G4double mom = alpha->GetTotalMomentum();

  //Random, Isotropic direction
  G4double cosTheta = 2*G4UniformRand() - 1., phi = twopi*G4UniformRand();
  G4double sinTheta = std::sqrt(1. - cosTheta*cosTheta);
  G4double ux = sinTheta*std::cos(phi),
    uy = sinTheta*std::sin(phi),
    uz = cosTheta;

  G4double momx = mom*ux;
  G4double momy = mom*uy;
  G4double momz = mom*uz;


  G4PrimaryParticle* particle = new G4PrimaryParticle(alpha_definition,
                                                      momx, momy, momz,
                                                      decayenergy);

  particle->SetMass(alpha_definition->GetPDGMass());
  particle->SetCharge(alpha_definition->GetPDGCharge());
  // Fix (where should polarization come from)
  // particle->SetPolarization(polar.x(), polar.y(), polar.z());

  if(verboseLevel >= 2) {
    G4cout << " Alpha Decay of " << DecayIsotope << G4endl;
    G4cout << "     Energy:   " << decayenergy << G4endl;
    G4cout << "     Momentum: "<<particle->GetMomentum() << G4endl;
    G4cout << "     Position: "<<vertex->GetPosition() << G4endl;
  }

  vertex->SetPrimary(particle);
  anEvent->AddPrimaryVertex(vertex);
}

void NMSSingleDecaySource::GenerateBetaPrimaryVertex(G4Event *anEvent,  G4PrimaryVertex* vertex) {

  // FIX
}

void NMSSingleDecaySource::GenerateAlphaNPrimaryVertex(G4Event *anEvent) {
  G4ParticleDefinition* neutron_definition = G4Neutron::Neutron();

  G4int posindex = floor(G4UniformRand() * alphanset->size());
  NMSAlphaNReaction reaction = alphanset->operator[](posindex);
  G4PrimaryVertex* vertex = new G4PrimaryVertex(reaction.position, GetParticleTime());
  G4DynamicParticle* neutron;
  neutron = new G4DynamicParticle();
  neutron->SetDefinition(neutron_definition);
  neutron->SetMomentum(reaction.alphaDirection);
  
  NMSPrimaryUserInformation* info = new NMSPrimaryUserInformation();
  info->SetOrigin(ORIGIN_ALPHA_N);

  G4PrimaryParticle* particle = new G4PrimaryParticle(neutron_definition,
                                                      reaction.alphaDirection.x(), 
						      reaction.alphaDirection.y(),
						      reaction.alphaDirection.z(),
                                                      reaction.energy);
  particle->SetUserInformation(info);
  particle->SetMass(neutron_definition->GetPDGMass());
  particle->SetCharge(neutron_definition->GetPDGCharge());
  vertex->SetPrimary(particle);
  anEvent->AddPrimaryVertex(vertex);
}

void NMSSingleDecaySource::GeneratePrimaryVertex(G4Event *anEvent) {

  if(DecayType == NMSDECAY_ALPHA_N) {
    if(alphanset == 0) {
      G4cout << "================================================================" << G4endl;
      G4cout << "Error: NMSSingleDecaySource" << G4endl;
      G4cout << "No data set of (alpha,n) reaction data has been loaded." << G4endl;
      exit(0);
    }
    GenerateAlphaNPrimaryVertex(anEvent);
  }
  else {
    G4ThreeVector source_position = GetPosDist()->GenerateOne();

    if(verboseLevel >= 2) {
      G4cout << "New Source Event" << G4endl;
      G4cout << "Time:     " << GetParticleTime() / second << G4endl;
      G4cout << "Position: " << source_position / cm << G4endl;
      G4cout << "Position Distribution Type: " << GetPosDist()->GetPosDisType() << G4endl;
    }

    G4PrimaryVertex* vertex = new G4PrimaryVertex(source_position, GetParticleTime());

    // depending on decay type generate particles
    if(DecayType == NMSDECAY_SF or DecayType == NMSDECAY_SF_N or DecayType == NMSDECAY_SF_GAMMA) {
      // if (verboseLevel >= 1) {
      //   G4cout << "Spontaneous Fission" << G4endl;
      // }
      // output in special function
      GenerateSFPrimaryVertex(anEvent, vertex);
    }
    if(DecayType == NMSDECAY_ALPHA) {
      GenerateAlphaPrimaryVertex(anEvent, vertex);
    }
    if(DecayType == NMSDECAY_BETA) {
      // if (verboseLevel >= 1) {
      //   G4cout << "Beta Decay" << G4endl;
      // }
      // output ín special function
      GenerateBetaPrimaryVertex(anEvent, vertex);
    }
  }

  if (verboseLevel > 1)
    G4cout << "Primary Vertex generated (NMS Decay) !" << G4endl;
}


void NMSSingleDecaySource::SetVerboseLevel(G4int verbose) {
  verboseLevel = verbose;
}
