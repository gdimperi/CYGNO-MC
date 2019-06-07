#ifndef CYGNOSteppingAction_h
#define CYGNOSteppingAction_h 1

#include "G4UserSteppingAction.hh"
#include "CYGNODetectorConstruction.hh"
//#include "CYGNOEventAction.hh"
#include "CYGNOAnalysis.hh"

class CYGNOSteppingAction : public G4UserSteppingAction
{
  public:
	  CYGNOSteppingAction(CYGNODetectorConstruction*); //, CYGNOEventAction*);
	  ~CYGNOSteppingAction(){};

	  void UserSteppingAction(const G4Step*);
  private:

	  CYGNODetectorConstruction* fDetector;
	  //CYGNOEventAction* fEventAction;
};


#endif
