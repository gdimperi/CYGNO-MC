#ifndef CYGNOTritiumPhysics_h
#define CYGNOTritiumPhysics_h 1

//#include "G4VUserPhysicsList.hh"
#include "G4VPhysicsConstructor.hh"
#include "globals.hh"

//class CYGNOTritiumPhysics: public G4VUserPhysicsList
class CYGNOTritiumPhysics: public G4VPhysicsConstructor
{
public:
  /// constructor
  CYGNOTritiumPhysics();
  /// destructor
  virtual ~CYGNOTritiumPhysics();
  virtual void ConstructParticle();
  virtual void ConstructProcess();    
  virtual void ConstructDecay();
  /// Set user cuts
  //virtual void SetCuts();
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

