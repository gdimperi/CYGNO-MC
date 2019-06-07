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

  fOutFileCutCmd = new G4UIcmdWithAnInteger("/CYGNO/cutoutfile",this);
  fOutFileCutCmd->SetGuidance("  Choice : 0 1 2");
  fOutFileCutCmd->SetGuidance("If 0 saves everything; if 1 saves only events with energy in sensitive volume");
  fOutFileCutCmd->SetParameterName("cutoutfile",true);
  fOutFileCutCmd->SetDefaultValue(1);

  fRegisterOnCmd = new G4UIcmdWithAnInteger("/CYGNO/registeron",this);
  fRegisterOnCmd->SetGuidance("  Choice : 0 1");
  fRegisterOnCmd->SetGuidance("If 1 saves flux, neutron and isotope data; if 0 does not save these data");
  fRegisterOnCmd->SetParameterName("registeron",true);
  fRegisterOnCmd->SetDefaultValue(0);

  fOutFileSaveHitCmd = new G4UIcmdWithAnInteger("/CYGNO/save_hits_branches",this);
  fOutFileSaveHitCmd->SetGuidance("Choice : 0 1");
  fOutFileSaveHitCmd->SetGuidance("If 0 the tree branches containing the hits information for each event are not saved");
  fOutFileSaveHitCmd->SetParameterName("save_hits_branches",true);
  fOutFileSaveHitCmd->SetDefaultValue(0);
  
}

CYGNORunActionMessenger::~CYGNORunActionMessenger()
{
  delete fOutFileCmd;
  delete fOutFileCutCmd;
  delete fRegisterOnCmd;
  delete fOutFileSaveHitCmd;
}

//---------------------------------------------------------------------------//

void CYGNORunActionMessenger::SetNewValue(G4UIcommand *command,
					      G4String newValue)
{
  if (command == fOutFileCmd ) 
    {
      fRunAction->SetOutFile(newValue);
    }
  else if (command == fOutFileCutCmd )
    {
      G4int cutoutfile = fOutFileCutCmd->GetNewIntValue(newValue);
      if(cutoutfile>=0 && cutoutfile<=2)
	{
	  fRunAction->SetOutFileCut(cutoutfile);
	}
      else
	{
	  G4cout << "WARNING! The setting for /CYGNO/cutoutfile is not in the range allowed. Using the default cut (/CYGNO/cutoutfile 1)" << G4endl;
	}
      //fRunAction->SetOutFileCut(fOutFileCutCmd->GetNewIntValue(newValue));
    }
  else if (command == fRegisterOnCmd )
    {
      fRunAction->SetRegisterOn(fRegisterOnCmd->GetNewIntValue(newValue));
    }
  else if (command == fOutFileSaveHitCmd )
    {
        fRunAction->SetHitsInfo(fOutFileSaveHitCmd->GetNewIntValue(newValue));
    }
}
