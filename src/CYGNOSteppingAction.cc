#include "CYGNOSteppingAction.hh"
#include "CYGNOVolumes.hh"
#include "G4SteppingManager.hh"
#include "G4Track.hh"
#include "globals.hh"
#include "G4Ions.hh"
#include "G4Neutron.hh"
#include "G4Proton.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4MuonMinus.hh"
#include "G4MuonPlus.hh"
#include "G4Gamma.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "G4UnitsTable.hh"

CYGNOSteppingAction::CYGNOSteppingAction(CYGNODetectorConstruction* det):
fDetector(det) {
  totalEnergyDepositAccum = 0.0;
  IonizingEnergyDepositAccum = 0.0;
}

void CYGNOSteppingAction::UserSteppingAction(const G4Step* fStep)
{ 
  G4Track* track = fStep->GetTrack();
    
  // Get the particle definition and the particle name
  G4ParticleDefinition* particleDefinition = track->GetDefinition();
  G4String particleName = particleDefinition->GetParticleName();
  G4int particleID = particleDefinition->GetPDGEncoding();

  // Get the volume the particle is entering
  G4StepPoint* preStepPoint = fStep->GetPreStepPoint();
  G4LogicalVolume* preVolume = preStepPoint->GetTouchableHandle()->GetVolume()->GetLogicalVolume();
  std::string volumeName = preVolume->GetName();

  if (volumeName == "CYGNO_log" && particleID > 1000000000) { 
      // Retrieve total energy deposited in the step
      G4double totalEnergyDeposit = fStep->GetTotalEnergyDeposit();
      
      // Retrieve non-ionizing energy deposited in the step
      G4double nonIonizingEnergyDeposit = fStep->GetNonIonizingEnergyDeposit();
      
      // Calculate ionizing energy deposit
      G4double ionizingEnergyDeposit = totalEnergyDeposit - nonIonizingEnergyDeposit;

      totalEnergyDepositAccum += totalEnergyDeposit;
      IonizingEnergyDepositAccum += ionizingEnergyDeposit;

  }

  if (track->GetTrackStatus() == fStopAndKill) {
      G4cout << "End of Track - Accumulated Energy Deposits for Particle: " 
              << track->GetDefinition()->GetParticleName() << " (PDG ID: " 
              << track->GetDefinition()->GetPDGEncoding() << ")" << G4endl;
      G4cout << "Total Energy Deposited during track: " 
              << totalEnergyDepositAccum / CLHEP::keV << " keV" << G4endl;
      G4cout << "Total Ionizing Energy Deposited during track: " 
              << IonizingEnergyDepositAccum / CLHEP::keV << " keV" << G4endl;

      // Reset accumulators for the next track
      totalEnergyDepositAccum = 0.0;
      IonizingEnergyDepositAccum = 0.0;
  }

}

