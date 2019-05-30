#ifndef CYGNOActionInitialization_h
#define CYGNOActionInitialization_h 1

#include "G4VUserActionInitialization.hh"
#include "CYGNODetectorConstruction.hh"

/// Action initialization class.

class CYGNOActionInitialization : public G4VUserActionInitialization
{
  public:
    CYGNOActionInitialization(CYGNODetectorConstruction*);
    virtual ~CYGNOActionInitialization();

    virtual void BuildForMaster() const;
    virtual void Build() const;

  private:
      CYGNODetectorConstruction* fDetector;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

    
