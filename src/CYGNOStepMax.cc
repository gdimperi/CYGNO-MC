#include "CYGNOStepMax.hh"
#include "CYGNOStepMaxMessenger.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

CYGNOStepMax::CYGNOStepMax(const G4String& processName)
 : G4VDiscreteProcess(processName),fMaxChargedStep(DBL_MAX),fMess(0)
{
  fMess = new CYGNOStepMaxMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

CYGNOStepMax::~CYGNOStepMax() { delete fMess; }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool CYGNOStepMax::IsApplicable(const G4ParticleDefinition& particle)
{
  return (particle.GetPDGCharge() != 0.);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void CYGNOStepMax::SetMaxStep(G4double step) {fMaxChargedStep = step;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double CYGNOStepMax::PostStepGetPhysicalInteractionLength( 
                                                 const G4Track& aTrack,
                                                       G4double,
                                                       G4ForceCondition* condition )
{
  // condition is set to "Not Forced"
  *condition = NotForced;
  G4double ProposedStep = fMaxChargedStep;

  return ProposedStep;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VParticleChange* CYGNOStepMax::PostStepDoIt(const G4Track& aTrack, const G4Step&)
{
   // do nothing
   aParticleChange.Initialize(aTrack);
   return &aParticleChange;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


