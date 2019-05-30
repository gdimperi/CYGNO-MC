#ifndef CYGNORunActionMessenger_h
#define CYGNORunActionMessenger_h 1

//---------------------------------------------------------------------------//

#include "globals.hh"
#include "G4UImessenger.hh"

//---------------------------------------------------------------------------//

class G4UIdirectory;
class G4UIcommand;
class G4UIcmdWithAString;
class G4UIcmdWithABool;
class G4UIcmdWithAnInteger;
class CYGNORunAction;

//---------------------------------------------------------------------------//

class CYGNORunActionMessenger : public G4UImessenger
{
public:

  //default constructor
  CYGNORunActionMessenger(CYGNORunAction *runAct);

  //destructor
  ~CYGNORunActionMessenger();

  //public interface
  void SetNewValue(G4UIcommand *command, G4String newValues);

  //protected members
protected:

  //private  members
private:
  CYGNORunAction       *fRunAction;
  G4UIdirectory        *fDirectory;
  G4UIcmdWithAString   *fOutFileCmd;
  G4UIcmdWithAnInteger *fOutFileCutCmd;
  G4UIcmdWithAnInteger *fRegisterOnCmd;
  G4UIcmdWithAnInteger *fTotTCmd;
  G4UIcmdWithAnInteger  *fOutFileSaveHitCmd;
};
#endif
