#ifndef CYGNODetectorProperty_h
#define CYGNODetectorProperty_h 1

#include "globals.hh"
#include <vector>
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"

class CYGNODetectorProperty
{
public:

  CYGNODetectorProperty();
  ~CYGNODetectorProperty();
  static CYGNODetectorProperty * GetInstance();
  void Refresh();
  void AddVolumeNameMass(G4String, G4double);
  void AddVolumeNameDensity(G4String, G4double);
  void AddVolumeNameMassAndDensity(G4LogicalVolume*);
  void AddPhysVolumeNameMassAndDensity(G4VPhysicalVolume*);
  std::vector<std::pair <G4String,G4double> > *GetVolumeNameMass() {return vol_name_mass;}
  std::vector<std::pair <G4String,G4double> > *GetVolumeNameDensity() {return vol_name_dens;}
  static CYGNODetectorProperty* fCYGNODetectorProperty;
  const G4double tolerance;//This is subtracted from the size of each solid, instead of translating each one apart since that is more complicated.

private:

  std::vector<std::pair <G4String,G4double> > *vol_name_mass;
  std::vector<std::pair <G4String,G4double> > *vol_name_dens;
  G4int debugmode;
};

#endif 
