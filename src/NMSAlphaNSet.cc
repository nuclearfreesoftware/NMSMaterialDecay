#include "NMSAlphaNSet.hh"

NMSAlphaNSet::NMSAlphaNSet() : std::vector<NMSAlphaNReaction>() {
  
}

NMSAlphaNSet::~NMSAlphaNSet() {
  
}

void NMSAlphaNSet::saveToFile(G4String filename) {
  std::ofstream outfile;
  outfile.open(filename.c_str());
  for(int i = 0; i < this->size(); i++) {
    // add two threevectors position and alphaDirection
    outfile << std::setprecision(15) << this->std::vector<NMSAlphaNReaction>::operator[](i).energy << "," 
	    << this->std::vector<NMSAlphaNReaction>::operator[](i).time << std::endl;
  }
  outfile.close();

}

void NMSAlphaNSet::loadFromFile(G4String filename) {
  
}
