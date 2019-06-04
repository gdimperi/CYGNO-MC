#include "globals.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithADouble.hh"

#include "CYGNODetectorMaterial.hh"
#include "CYGNODetectorMaterialMessenger.hh"

//---------------------------------------------------------------------------//

CYGNODetectorMaterialMessenger::CYGNODetectorMaterialMessenger(CYGNODetectorMaterial *material):
fDetectorMaterials(material)
{
    fGasDirectory = new G4UIdirectory("/CYGNO/gas/");
    fGasDirectory->SetGuidance("Control commands for gas:");
    
    fGasHeCmd = new G4UIcmdWithADouble("/CYGNO/gas/frac_He",this);
    fGasHeCmd->SetGuidance("Set fraction of He gas (<=1)");
    fGasHeCmd->SetParameterName("choice",false);
    fGasHeCmd->AvailableForStates(G4State_Init,G4State_Idle);
    
    fGasCF4Cmd = new G4UIcmdWithADouble("/CYGNO/gas/frac_CF4",this);
    fGasCF4Cmd->SetGuidance("Set fraction of CF4 gas (<=1)");
    fGasCF4Cmd->SetParameterName("choice",false);
    fGasCF4Cmd->AvailableForStates(G4State_Init,G4State_Idle);

}

CYGNODetectorMaterialMessenger::~CYGNODetectorMaterialMessenger()
{
    delete fGasHeCmd;
    delete fGasCF4Cmd;
}

//---------------------------------------------------------------------------//

void CYGNODetectorMaterialMessenger::SetNewValue(G4UIcommand *command,
                                                     G4String newValue)
{
  if (command == fGasHeCmd ){ 
        fDetectorMaterials->SetGasHeFrac( fGasHeCmd->GetNewDoubleValue(newValue));
  }
  else if (command == fGasCF4Cmd ){ 
        fDetectorMaterials->SetGasCF4Frac( fGasCF4Cmd->GetNewDoubleValue(newValue));
  }
}
