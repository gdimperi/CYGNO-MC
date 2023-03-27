#ifndef CYGNOPrimaryGeneratorAction_h
#define CYGNOPrimaryGeneratorAction_h 1

//#include "CYGNOGeneratorPositionSampling.hh"
#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4ThreeVector.hh"
#include <fstream>
#include <stdio.h>
#include <vector>
#include "globals.hh"

#include "G4GenericMessenger.hh"
#include "G4ParticleGun.hh"

class CYGNOPrimaryGeneratorActionMessenger;
class CYGNODetectorConstruction;
class G4GeneralParticleSource;
class G4Event;
class G4ParticleGun;

class CYGNOPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
  public:
    CYGNOPrimaryGeneratorAction(CYGNODetectorConstruction*);    
   ~CYGNOPrimaryGeneratorAction();

  public:
  void GeneratePrimaries(G4Event*);  
  void SetEnergy(G4double ene) {particle_energy = ene;}
  G4double GetEnergy() {return particle_energy;}
  void ChangeFileName(G4String newFile);
  void SetGenerator(const G4String& name);

  private:
  G4GeneralParticleSource* particleGun;
  G4ParticleGun*      MuonParticleGun;
  G4ParticleGun*      ER_Gun;
  CYGNODetectorConstruction* myDetector;
  G4int n_particle;
  G4double particle_energy;

  G4GenericMessenger* fMessenger;
  G4String           fGenerator;
  std::ifstream      fInputFile;
  G4String           fFileName;
  void DefineCommands();
 
  G4ThreeVector position;
//  CYGNOPrimaryGeneratorActionMessenger* fMessenger;
};

#endif


