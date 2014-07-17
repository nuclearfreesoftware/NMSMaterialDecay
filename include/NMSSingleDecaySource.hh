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
 *
 * This file incorporates work covered by the following copyright and  
 * permission notice:
 * 
 *     Copyright (c) 2006-2010 Lawrence Livermore National Security, LLC.
 *     Produced at the Lawrence Livermore National Laboratory 
 *     UCRL-CODE-224807.
 *     
 *     All rights reserved. Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
 *     
 *     o   Redistributions of source code must retain the above copyright notice, this list of conditions and the disclaimer below.
 *     
 *     o  Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the disclaimer (as noted below) in the documentation and/or other materials provided with the distribution.
 *     
 *     o  Neither the name of the LLNS/LLNL nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
 *     
 *     THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL LAWRENCE LIVERMORE NATIONAL SECURITY, LLC, THE U.S. DEPARTMENT OF ENERGY OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *     Additional BSD Notice
 *     
 *     1. This notice is required to be provided under our contract with the U.S. Department of Energy (DOE). This work was produced at Lawrence Livermore National Laboratory under Contract No. DE-AC52-07NA27344 with the DOE. 
 *     
 *     2. Neither the United States Government nor Lawrence Livermore National Security, LLC nor any of their employees, makes any warranty, express or implied, or assumes any liability or responsibility for the accuracy, completeness, or usefulness of any information, apparatus, product, or process disclosed, or represents that its use would not infringe privately-owned rights. 
 *     
 *     3. Also, reference herein to any specific commercial products, process, or services by trade name, trademark, manufacturer or otherwise does not necessarily constitute or imply its endorsement, recommendation, or favoring by the United States Government or Lawrence Livermore National Security, LLC. The views and opinions of authors expressed herein do not necessarily state or reflect those of the United States Government or Lawrence Livermore National Security, LLC, and shall not be used for advertising or product endorsement purposes.
 *
 * 
 */

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

#include "NMSAlphaNSet.hh"
#include "NMSPrimaryUserInformation.hh"

class NMSPrimaryUserInformation;

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
  void setAlphaNFile(G4String filename);

  void GenerateAlphaPrimaryVertex(G4Event* anEvent, G4PrimaryVertex* vertex);
  void GenerateSFPrimaryVertex(G4Event *anEvent,  G4PrimaryVertex* vertex);
  void GenerateBetaPrimaryVertex(G4Event *anEvent,  G4PrimaryVertex* vertex);
  void GenerateAlphaNPrimaryVertex(G4Event* anEvent);
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

  NMSAlphaNSet* alphanset;

  G4bool decayloaded;
};

#endif
