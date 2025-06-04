#include "CYGNOSensitiveDetector.hh"
#include "G4HCofThisEvent.hh"
#include "G4Step.hh"
#include "G4VProcess.hh"
#include "G4ThreeVector.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"
#include "G4EmCalculator.hh"
#include "G4ParticleTable.hh"
#include "G4Material.hh"

CYGNOSensitiveDetector::CYGNOSensitiveDetector(G4String name)
  :G4VSensitiveDetector(name), fHCID(-1)
{
  G4String HCname;
  collectionName.insert(HCname="CYGNOCollection");
}

CYGNOSensitiveDetector::~CYGNOSensitiveDetector(){ }

void CYGNOSensitiveDetector::Initialize(G4HCofThisEvent* HCE)
{
  CYGNOCollection = new CYGNOHitsCollection(SensitiveDetectorName,collectionName[0]); 
  if(fHCID<0)
  { fHCID = G4SDManager::GetSDMpointer()->GetCollectionID(CYGNOCollection); }
  HCE->AddHitsCollection( fHCID, CYGNOCollection ); 
}

G4bool CYGNOSensitiveDetector::ProcessHits(G4Step* aStep,G4TouchableHistory*)
{
  G4double edep = aStep->GetTotalEnergyDeposit();
  if ( aStep->GetTrack()->GetDefinition()->GetParticleName() == "opticalphoton" ) {
  //    G4cout<<"opt photon"<<G4endl;
  // optical photon do not count as hit
  return false;
  }
  
  if(edep==0.) return true;

  CYGNOHit* newHit = new CYGNOHit();
  // parent ID
  newHit->SetParentID  (aStep->GetTrack()->GetParentID());
  // particle ID
  newHit->SetParticleID  (aStep->GetTrack()->GetDefinition()->GetPDGEncoding());
  // track ID
  newHit->SetTrackID  (aStep->GetTrack()->GetTrackID());
  // global time in sec
  newHit->SetGlobalTime  (aStep->GetTrack()->GetGlobalTime()/s);
  // kinetic energy in keV
  newHit->SetKineticEne  (aStep->GetPreStepPoint()->GetKineticEnergy()/keV);
    G4double energy = aStep->GetTrack()->GetKineticEnergy();
    G4String particleName = aStep->GetTrack()->GetParticleDefinition()->GetParticleName();
    G4String materialName = aStep->GetTrack()->GetMaterial()->GetName();

    static G4EmCalculator emCal;
    G4double range = emCal.GetRange(energy, particleName, materialName);
      
    // You can use this to compare with step length
    G4double stepLength = aStep->GetStepLength();  
//    if (aStep->GetTrack()->GetKineticEnergy()/keV <= 1.07 ) { //choose the threshold because of the step loss function 
//  // step length in mm
//  newHit->SetLength  (range/mm);
//   G4cout << "Sensitive det class: remaining range = " << newHit->GetLength() << " mm" << G4endl;
//    }
//    else{
  newHit->SetLength  (aStep->GetStepLength()/mm);
    
//    }
  G4double eDepTotal = aStep->GetTotalEnergyDeposit();
  G4double eDepNonIon = aStep->GetNonIonizingEnergyDeposit();
  G4double eDepIon = eDepTotal - eDepNonIon;
  
  // Save eDepIon in your CYGNOHit object (in keV)
  newHit->SetIonizingEnergy(eDepIon/keV);

   
  if (aStep->GetPreStepPoint()->GetProcessDefinedStep() != NULL) {
    newHit->SetProcessIni  (aStep->GetPreStepPoint()->GetProcessDefinedStep()->GetProcessName());
  }
  else  newHit->SetProcessIni  ("none");

  newHit->SetProcessFin  (aStep->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName());

  G4String volumeName = aStep->GetPreStepPoint()->GetTouchableHandle()->GetVolume()->GetName();
  newHit->SetEdep     (edep/keV);
  newHit->SetPos      (aStep->GetPostStepPoint()->GetPosition());
  CYGNOCollection->insert( newHit );
  return true;


}

void CYGNOSensitiveDetector::EndOfEvent(G4HCofThisEvent*)
{
  verboseLevel=0;

  if (verboseLevel>0) { 
    G4int NbHits = CYGNOCollection->entries();
    if(NbHits>0)
      {
	G4cout << "\n-------->Hits Collection: in this event there are " << NbHits 
	       << " hits in CYGNO gas: " << G4endl;
	for (G4int i=0;i<NbHits;i++) (*CYGNOCollection)[i]->Print();   
      } 
  }
}

