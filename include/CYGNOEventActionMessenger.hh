#ifndef CYGNOEventActionMessenger_h
#define CYGNOEventActionMessenger_h 1

//---------------------------------------------------------------------------//

#include "globals.hh"
#include "G4UImessenger.hh"

//---------------------------------------------------------------------------//

class G4UIdirectory;
class G4UIcommand;
class G4UIcmdWithAnInteger;
class CYGNOEventAction;

//---------------------------------------------------------------------------//

class CYGNOEventActionMessenger : public G4UImessenger
{
public:

  //default constructor
  CYGNOEventActionMessenger(CYGNOEventAction *evact);

  //copy constructor
  //CYGNOEventActionMessenger();

  //destructor
  ~CYGNOEventActionMessenger();

  //public interface
  void SetNewValue(G4UIcommand *command, G4String newValues);

  //protected members
protected:

  //private  members
private:
  CYGNOEventAction        *fEvAct;
  G4UIdirectory        *fDirectory;
  G4UIcmdWithAnInteger *fRepFreqCmd;
};
#endif
