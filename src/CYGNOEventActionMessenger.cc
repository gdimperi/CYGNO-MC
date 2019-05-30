#include "globals.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIdirectory.hh"
#include "G4PhysicalVolumeStore.hh"

#include "CYGNOEventAction.hh"
#include "CYGNOEventActionMessenger.hh" 

//---------------------------------------------------------------------------//

CYGNOEventActionMessenger::CYGNOEventActionMessenger
(CYGNOEventAction *evact):fEvAct(evact)
{
  fDirectory = new G4UIdirectory("/CYGNO/");
  
  fRepFreqCmd = new G4UIcmdWithAnInteger("/CYGNO/reportingfrequency",this);
  fRepFreqCmd->SetGuidance("Set the reporting frequency");
  fRepFreqCmd->SetGuidance("(default: 1000)");
  fRepFreqCmd->SetParameterName("choice",false);
  fRepFreqCmd->SetRange("choice>0");
}

CYGNOEventActionMessenger::~CYGNOEventActionMessenger()
{
  delete fRepFreqCmd;
}

//---------------------------------------------------------------------------//

void CYGNOEventActionMessenger::SetNewValue(G4UIcommand *command,
					      G4String newValue)
{
  if (command == fRepFreqCmd ) 
    {
      fEvAct->SetRepFreq(fRepFreqCmd->GetNewIntValue(newValue));
    }
}
