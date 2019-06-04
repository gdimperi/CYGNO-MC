#ifndef CYGNODetectorMaterialMessenger_h
#define CYGNODetectorMaterialMessenger_h 1

//---------------------------------------------------------------------------//

#include "globals.hh"
#include "G4UImessenger.hh"

//---------------------------------------------------------------------------//

class G4UIdirectory;
class G4UIcommand;
class G4UIcmdWithAString;
class G4UIcmdWithABool;
class CYGNODetectorMaterial;
class G4UIcmdWithADouble;

//---------------------------------------------------------------------------//

class CYGNODetectorMaterialMessenger : public G4UImessenger
{
public:
    
    //default constructor
    CYGNODetectorMaterialMessenger(CYGNODetectorMaterial *materials);
    
    //destructor
    ~CYGNODetectorMaterialMessenger();
    
    //public interface
    void SetNewValue(G4UIcommand *command, G4String newValues);
    
    //protected members
protected:
    
    //private  members
private:
    CYGNODetectorMaterial       *fDetectorMaterials;
    
    G4UIdirectory               *fGasDirectory;
    G4UIcmdWithADouble          *fGasHeCmd;
    G4UIcmdWithADouble          *fGasCF4Cmd;
   
};
#endif
