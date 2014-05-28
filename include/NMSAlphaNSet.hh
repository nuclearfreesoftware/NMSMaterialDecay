#ifndef NMSAlphaNSet_h
#define NMSAlphaNSet_h 1

#include <vector>
#include <fstream>
#include <iomanip>

#include "NMSAlphaNReaction.hh"

class NMSAlphaNSet : public std::vector<NMSAlphaNReaction> {
public:  
  NMSAlphaNSet();
  virtual ~NMSAlphaNSet();

  void saveToFile(G4String filename);
  void loadFromFile(G4String filename);
  
};

#endif
