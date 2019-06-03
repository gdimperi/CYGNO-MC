#include "globals.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWith3VectorAndUnit.hh"
#include "G4UIdirectory.hh"
#include "G4PhysicalVolumeStore.hh"

#include "CYGNODetectorLNGS.hh"
#include "CYGNODetectorLNGSMessenger.hh"

//---------------------------------------------------------------------------//

CYGNODetectorLNGSMessenger::CYGNODetectorLNGSMessenger
(CYGNODetectorLNGS *LNGS):fLNGS(LNGS)
{
  
  fexternalrockthickCmd = new G4UIcmdWithADoubleAndUnit("/CYGNO/lab/example",this);
  fexternalrockthickCmd->SetGuidance("This is just an example");
  fexternalrockthickCmd->SetDefaultUnit("m");
  fexternalrockthickCmd->SetUnitCandidates("mm cm m");
  fexternalrockthickCmd->SetParameterName("choice",false);
  fexternalrockthickCmd->AvailableForStates(G4State_Init,G4State_Idle);
  
}

CYGNODetectorLNGSMessenger::~CYGNODetectorLNGSMessenger()
{
  delete fexternalrockthickCmd;
}

//---------------------------------------------------------------------------//

void CYGNODetectorLNGSMessenger::SetNewValue(G4UIcommand *command,
                                                     G4String newValue)
{
  if (command == fexternalrockthickCmd )
    {
	  fLNGS->SetExternalRockThickness(fexternalrockthickCmd->GetNewDoubleValue(newValue));
    }
}
