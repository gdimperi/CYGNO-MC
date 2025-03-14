
#include "CYGNOStepMaxMessenger.hh"

#include "CYGNOStepMax.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "globals.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

CYGNOStepMaxMessenger::CYGNOStepMaxMessenger(CYGNOStepMax* stepM)
  :G4UImessenger(),fCYGNOStepMax(stepM),fCYGNOStepMaxCmd(0)
{ 
  fCYGNOStepMaxCmd = new G4UIcmdWithADoubleAndUnit("/CYGNO/stepmax",this);
  fCYGNOStepMaxCmd->SetGuidance("Set max allowed step length");
  fCYGNOStepMaxCmd->SetParameterName("mxStep",false);
  fCYGNOStepMaxCmd->SetRange("mxStep>0.");
  fCYGNOStepMaxCmd->SetUnitCategory("Length");
  fCYGNOStepMaxCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

CYGNOStepMaxMessenger::~CYGNOStepMaxMessenger()
{
  delete fCYGNOStepMaxCmd;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void CYGNOStepMaxMessenger::SetNewValue(G4UIcommand* command, G4String newValue)
{ 
  if (command == fCYGNOStepMaxCmd)
    { fCYGNOStepMax->SetMaxStep(fCYGNOStepMaxCmd->GetNewDoubleValue(newValue));}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
