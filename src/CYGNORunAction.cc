#include "CYGNORunAction.hh"
#include "CYGNOEventAction.hh"
#include "CYGNORunActionMessenger.hh"
#include "CYGNOAnalysis.hh"
#include "G4AnalysisManager.hh"
#include "CYGNODetectorConstruction.hh"
#include "G4Run.hh"


CYGNORunAction::CYGNORunAction(CYGNOEventAction* myEA, CYGNODetectorConstruction* myDC)
  : G4UserRunAction(), fDetector(myDC), fEventAction(myEA)
{
  fMessenger = new CYGNORunActionMessenger(this);
}

CYGNORunAction::~CYGNORunAction()
{
  delete fMessenger;
}

void CYGNORunAction::BeginOfRunAction(const G4Run* aRun)
{
		//Master mode or sequential
  if (IsMaster())    
    G4cout << "### Run " << aRun->GetRunID() << " starts (master)." << G4endl;
  else
    G4cout << "### Run " << aRun->GetRunID() << " starts (worker)." << G4endl;
  
  // Book histograms and ntuples
    CYGNOAnalysis::Instance()->BookNtuples();
    CYGNOAnalysis::Instance()->OpenFile(FileName);
}

void CYGNORunAction::EndOfRunAction(const G4Run* aRun)
{ 
    CYGNOAnalysis::Instance()->CloseFile();

}


