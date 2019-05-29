#ifndef CYGNOSensitiveDetector_h
#define CYGNOSensitiveDetector_h 1

#include "G4VSensitiveDetector.hh"
#include "CYGNODetectorConstruction.hh"
#include "CYGNOHit.hh"

class G4Step;
class G4HCofThisEvent;
class G4TouchableHistory;

class CYGNOSensitiveDetector : public G4VSensitiveDetector
{
  public:
      CYGNOSensitiveDetector(G4String name);
      virtual ~CYGNOSensitiveDetector();

      virtual void Initialize(G4HCofThisEvent*);
      virtual G4bool ProcessHits(G4Step*, G4TouchableHistory*);
      virtual void EndOfEvent(G4HCofThisEvent*);

  private:
      CYGNOHitsCollection* CYGNOCollection;
      CYGNODetectorConstruction* fDetector;
      G4int fHCID;
};

#endif

