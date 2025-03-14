#ifndef CYGNOStepMaxMessenger_h
#define CYGNOStepMaxMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class CYGNOStepMax;
class G4UIcmdWithADoubleAndUnit;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class CYGNOStepMaxMessenger: public G4UImessenger
{
public:

  CYGNOStepMaxMessenger(CYGNOStepMax*);
  ~CYGNOStepMaxMessenger();
    
  virtual void SetNewValue(G4UIcommand*, G4String);
    
private:
  CYGNOStepMax* fCYGNOStepMax;
  G4UIcmdWithADoubleAndUnit* fCYGNOStepMaxCmd;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
