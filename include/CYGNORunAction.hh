#ifndef CYGNORunAction_h
#define CYGNORunAction_h 1

#include "G4UserRunAction.hh"
#include "globals.hh"

#include <vector>

class G4Run;
class CYGNORunActionMessenger;
class CYGNODetectorConstruction;
class CYGNOEventAction;

class CYGNORunAction : public G4UserRunAction
{
  public:
  CYGNORunAction(CYGNOEventAction*, CYGNODetectorConstruction*);
  virtual ~CYGNORunAction();
  
  virtual void BeginOfRunAction(const G4Run*);
  virtual void EndOfRunAction(const G4Run*);

  void SetOutFile(G4String fname) {FileName = fname;};

  private:

  void Book(const G4Run*);
  
  CYGNORunActionMessenger* fMessenger;  
  CYGNODetectorConstruction* fDetector;
  CYGNOEventAction* fEventAction;
  G4String FileName;
};

#endif
