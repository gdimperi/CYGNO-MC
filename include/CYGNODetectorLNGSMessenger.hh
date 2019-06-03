#ifndef CYGNODetectorLNGSMessenger_h
#define CYGNODetectorLNGSMessenger_h 1

//---------------------------------------------------------------------------//

#include "globals.hh"
#include "G4UImessenger.hh"

//---------------------------------------------------------------------------//

class G4UIdirectory;
class G4UIcommand;
class G4UIcmdWithAString;
class G4UIcmdWithAnInteger;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWith3VectorAndUnit;
class G4UIcmdWithoutParameter;
class CYGNODetectorLNGS;

//---------------------------------------------------------------------------//

class CYGNODetectorLNGSMessenger : public G4UImessenger
{
public:
    
    //default constructor
    CYGNODetectorLNGSMessenger(CYGNODetectorLNGS *LNGS);
    
    //destructor
    ~CYGNODetectorLNGSMessenger();
    
    //public interface
    void SetNewValue(G4UIcommand *command, G4String newValues);
    
    //protected members
protected:
    
    //private  members
private:
    CYGNODetectorLNGS  *fLNGS;
    
    G4UIcmdWithADoubleAndUnit   *fexternalrockthickCmd;
 
};
#endif
