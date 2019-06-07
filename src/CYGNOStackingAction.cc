#include "CYGNOStackingAction.hh"
#include "CYGNOAnalysis.hh"

#include "CYGNOUserEventInformation.hh"
#include "CYGNOSteppingAction.hh"

#include "G4ios.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTypes.hh"
#include "G4RunManager.hh"
#include "G4Event.hh"
#include "G4EventManager.hh"


#include "G4Track.hh"
#include "G4NeutrinoE.hh"
#include "G4AntiNeutrinoE.hh"
#include "G4TransportationManager.hh"
#include "G4VProcess.hh"
#include "G4Ions.hh"
#include "G4SystemOfUnits.hh"
//#include <CLHEP/Units/SystemOfUnits.h>

//using namespace CLHEP;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

CYGNOStackingAction::CYGNOStackingAction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

CYGNOStackingAction::~CYGNOStackingAction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ClassificationOfNewTrack 
CYGNOStackingAction::ClassifyNewTrack(const G4Track* track)
{

  CYGNOUserEventInformation* eventInformation=
    (CYGNOUserEventInformation*)G4EventManager::GetEventManager()
    ->GetConstCurrentEvent()->GetUserInformation();

  
  //kill secondary neutrino
  if (track->GetDefinition() == G4NeutrinoE::NeutrinoE() && track->GetParentID()>0) return fKill;
  if (track->GetDefinition() == G4AntiNeutrinoE::AntiNeutrinoE() && track->GetParentID()>0) return fKill;

  //save secondary radionuclides info
  G4String particleType = track->GetDefinition()->GetParticleType();
  if (particleType == "nucleus" && track->GetParentID()>0)
    {
      SaveIsotope(track);
    }

  return fUrgent;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void CYGNOStackingAction::SaveIsotope(const G4Track* track)
{
  // save secondary radionuclides info on a text file
            
  // get analysis instance                   
  CYGNOAnalysis* analysis = CYGNOAnalysis::getInstance();
  
  G4Ions* ion = (G4Ions*) track->GetDefinition();
  G4double lifetime = ion->GetPDGLifeTime();
  G4double trackID = track->GetTrackID();
  G4int pdg = track->GetDefinition()->GetPDGEncoding();
  G4int volNo = analysis->GetVolNo(track);
  G4int copyNo = analysis->GetCopyNo(track);
  G4int Z = ion->GetAtomicNumber();
  G4double excitationEnergy = ion->GetExcitationEnergy();
  if (lifetime > 1e18*second) lifetime = -1; // do not consider very long-lived isotopes > 3e10 yrs
  if (lifetime > 0 || excitationEnergy > 0)
    {
      G4int A = ion->GetAtomicMass();
      G4ThreeVector position = track->GetPosition(); 
      G4double kinE = track->GetKineticEnergy();

      analysis->RegisterIsotope(A, Z, pdg, kinE, position, volNo, copyNo, trackID);
    }
}

