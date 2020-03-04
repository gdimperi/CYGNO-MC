#ifndef CYGNOPhysicsList_h
#define CYGNOPhysicsList_h 1

#include "G4VModularPhysicsList.hh"
#include "globals.hh"

//class CYGNOPhysicsListMessenger;

class G4Cerenkov;
class G4Scintillation;
class G4OpAbsorption;
class G4OpMieHG;
class G4OpRayleigh;
class G4OpWLS;
class G4OpBoundaryProcess;
//class MyOpBoundaryProcess;
class CYGNOStepMax;

class CYGNOPhysicsList: public G4VModularPhysicsList
{
public:
  /// constructor
  CYGNOPhysicsList(G4int verbose = 1 , G4String low_energy_neutron_model = "HP", G4String HadrPhysVariant = "");
  /// destructor
  virtual ~CYGNOPhysicsList();
  void ConstructParticle();
  void ConstructProcess();    
  //void AddDecay();
  //void AddStepMax();       
  // Set user cuts
  virtual void SetCuts();
  void ConstructOp();
  void AddStepMax();


private:
  G4double cutForGamma;
  G4double cutForElectron;
  G4double cutForPositron;
  //CYGNOPhysicsListMessenger* pMessenger;
  static G4ThreadLocal G4Cerenkov* fCerenkovProcess;
  static G4ThreadLocal G4Scintillation* fScintillationProcess;
  static G4ThreadLocal G4OpAbsorption* fAbsorptionProcess;
  static G4ThreadLocal G4OpMieHG* fMieHGScatteringProcess;
  static G4ThreadLocal G4OpRayleigh* fRayleighScatteringProcess;
  static G4ThreadLocal G4OpWLS* fWLSProcess;
  static G4ThreadLocal G4OpBoundaryProcess* fBoundaryProcess;
//  static G4ThreadLocal MyOpBoundaryProcess* fMyBoundaryProcess;
  static G4ThreadLocal CYGNOStepMax* fStepMaxProcess;

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

