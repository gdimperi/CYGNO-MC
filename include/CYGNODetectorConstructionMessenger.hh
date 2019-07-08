#ifndef CYGNODetectorConstructionMessenger_h
#define CYGNODetectorConstructionMessenger_h 1

//---------------------------------------------------------------------------//

#include "globals.hh"
#include "G4UImessenger.hh"

//---------------------------------------------------------------------------//

class G4UIdirectory;
class G4UIcommand;
class G4UIcmdWithAString;
class G4UIcmdWithAnInteger;
class G4UIcmdWith3Vector;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWith3VectorAndUnit;
class G4UIcmdWithADouble;
class G4UIcmdWithoutParameter;
class G4UIcmdWithABool;
class CYGNODetectorConstruction;

//---------------------------------------------------------------------------//

class CYGNODetectorConstructionMessenger : public G4UImessenger
{
public:
    
    //default constructor
    CYGNODetectorConstructionMessenger(CYGNODetectorConstruction *detector);
    
    //copy constructor
    //CYGNODetectorConstructionMessenger();
    
    //destructor
    ~CYGNODetectorConstructionMessenger();
    
    //public interface
    void SetNewValue(G4UIcommand *command, G4String newValues);
    
    //protected members
protected:
    
    //private  members
private:
    CYGNODetectorConstruction   *fDetectorPrimary;

    G4UIdirectory               *fGeomPathDirectory;
    G4UIcmdWithAString          *fCYGNOCADPathCmd;

    G4UIdirectory               *fLabDirectory;
    G4UIcmdWithAString          *fCYGNOLabCmd;
    G4UIcmdWithADoubleAndUnit   *fexternalrockthickCmd;
    G4UIcmdWithADoubleAndUnit   *fproductionrockthickCmd;
    G4UIcmdWithADoubleAndUnit   *finternalrockthickCmd;

    G4UIdirectory               *fShieldDirectory;
    G4UIcmdWithAString          *fCYGNOShieldingCmd;
    G4UIcmdWithADoubleAndUnit   *fthick0Cmd;
    G4UIcmdWithADoubleAndUnit   *fthick1Cmd;
    G4UIcmdWithADoubleAndUnit   *fthick2Cmd;
    G4UIcmdWithADoubleAndUnit   *fthick3Cmd;
    G4UIcmdWithAString          *fMat0Cmd;
    G4UIcmdWithAString          *fMat1Cmd;
    G4UIcmdWithAString          *fMat2Cmd;
    G4UIcmdWithAString          *fMat3Cmd;
    
    G4UIcmdWithADoubleAndUnit   *fInsideVolumeRadiusCmd;
    G4UIcmdWithADoubleAndUnit   *fInsideVolumeHeightCmd;

    G4UIcmdWithoutParameter     *fupdateCmd;
};
#endif
