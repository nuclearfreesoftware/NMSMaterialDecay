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
    outfile << std::setprecision(15) 
            << this->std::vector<NMSAlphaNReaction>::operator[](i).position.x() << ","
            << this->std::vector<NMSAlphaNReaction>::operator[](i).position.y() << ","
            << this->std::vector<NMSAlphaNReaction>::operator[](i).position.z() << ","
            << this->std::vector<NMSAlphaNReaction>::operator[](i).alphaDirection.x() << ","
            << this->std::vector<NMSAlphaNReaction>::operator[](i).alphaDirection.y() << ","
            << this->std::vector<NMSAlphaNReaction>::operator[](i).alphaDirection.z() << ","
	    << this->std::vector<NMSAlphaNReaction>::operator[](i).energy << "," 
	    << this->std::vector<NMSAlphaNReaction>::operator[](i).time << std::endl;
  }
  outfile.close();

}

void NMSAlphaNSet::loadFromFile(G4String filename) {
  std::ifstream input;
  G4String test;
  G4String fieldstr;
  G4double value;
  NMSAlphaNReaction tempreaction;

  if(file_exists(filename)) {
    this->clear();

    input.open(filename.c_str());
    while ( input.good() ) {
      getline(input, test);
      //      input >> test;
      if(test.length() > 0) {
	std::istringstream iss( test );
	getline(iss, fieldstr, ',');
	std::istringstream( fieldstr ) >> value;
	tempreaction.position.setX(value * cm);
	getline(iss, fieldstr, ',');
	std::istringstream( fieldstr ) >> value;
	tempreaction.position.setY(value * cm);
	getline(iss, fieldstr, ',');
	std::istringstream( fieldstr ) >> value;
	tempreaction.position.setZ(value * cm);
	getline(iss, fieldstr, ',');
	std::istringstream( fieldstr ) >> value;
	tempreaction.alphaDirection.setX(value * cm);
	getline(iss, fieldstr, ',');
	std::istringstream( fieldstr ) >> value;
	tempreaction.alphaDirection.setY(value * cm);
	getline(iss, fieldstr, ',');
	std::istringstream( fieldstr ) >> value;
	tempreaction.alphaDirection.setZ(value * cm);
	getline(iss, fieldstr, ',');
	std::istringstream( fieldstr ) >> value;
	tempreaction.energy = value * MeV;
	getline(iss, fieldstr, ',');
	std::istringstream( fieldstr ) >> value;
	tempreaction.time = value * microsecond;
	push_back(tempreaction);
      }
    }
    input.close();
  }
  else {
    G4cout << "ERROR: could not find file " << filename << G4endl;
  }


}
