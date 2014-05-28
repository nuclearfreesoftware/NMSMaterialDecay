#ifndef NMSAlphaNSet_h
#define NMSAlphaNSet_h 1

#include <vector>
#include <fstream>
#include <iomanip>
#include <sys/stat.h>

#include "G4SystemOfUnits.hh"

#include "NMSAlphaNReaction.hh"

class NMSAlphaNSet : public std::vector<NMSAlphaNReaction> {
public:  
  NMSAlphaNSet();
  virtual ~NMSAlphaNSet();

  void saveToFile(G4String filename);
  void loadFromFile(G4String filename);

private:  
  inline bool file_exists(G4String filename) {
    struct stat buffer;
    return (stat (filename.c_str(), &buffer) == 0);
  }

};

#endif
