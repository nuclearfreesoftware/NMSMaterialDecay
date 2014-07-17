<<<<<<< HEAD
=======
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

>>>>>>> 151a70ae311b91d9941de27eb1dec38929162e02
#include "NMSPrimaryUserInformation.hh"

NMSPrimaryUserInformation::NMSPrimaryUserInformation() : G4VUserPrimaryParticleInformation() {
  
}

NMSPrimaryUserInformation::~NMSPrimaryUserInformation() {
  
}

void NMSPrimaryUserInformation::Print() const {
  
}

NMSOrigin NMSPrimaryUserInformation::GetOrigin() {
  return ori;
}

void NMSPrimaryUserInformation::SetOrigin(NMSOrigin newor) {
  ori = newor;
}

