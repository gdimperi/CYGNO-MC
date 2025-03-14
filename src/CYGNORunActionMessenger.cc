#include "globals.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIdirectory.hh"

#include "CYGNORunAction.hh"
#include "CYGNORunActionMessenger.hh" 

//---------------------------------------------------------------------------//

CYGNORunActionMessenger::CYGNORunActionMessenger(CYGNORunAction *runAct):fRunAction(runAct)
{
  fDirectory = new G4UIdirectory("/CYGNO/");
  
  fOutFileCmd = new G4UIcmdWithAString("/CYGNO/outfile",this);
  fOutFileCmd->SetGuidance("Set the output file name without extension");
  fOutFileCmd->SetGuidance("(default: out)");
  fOutFileCmd->SetParameterName("choice",false);
  
}

CYGNORunActionMessenger::~CYGNORunActionMessenger()
{
  delete fOutFileCmd;
}

//---------------------------------------------------------------------------//

void CYGNORunActionMessenger::SetNewValue(G4UIcommand *command,
					      G4String newValue)
{
  if (command == fOutFileCmd ) 
    {
      fRunAction->SetOutFile(newValue);
    }
}
