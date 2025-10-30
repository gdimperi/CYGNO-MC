#ifndef CYGNOEventAction_h
#define CYGNOEventAction_h 1

#include "G4UserEventAction.hh"
#include "CYGNORunAction.hh"

class G4Event;
class CYGNODetectorConstruction;
class G4Track;

class CYGNOEventAction : public G4UserEventAction
{
  public:
  CYGNOEventAction(CYGNODetectorConstruction*);
   ~CYGNOEventAction();

  public:
  void BeginOfEventAction(const G4Event* evt);
  void EndOfEventAction(const G4Event* evt);

  std::vector<G4int>& Get_pdgID_hits() { return v_pdgID_hits; };
  std::vector<G4double>& Get_tracklen_hits() { return v_tracklen_hits; };
  std::vector<G4double>& Get_px_particle() { return v_px_particle; };
  std::vector<G4double>& Get_py_particle() { return v_py_particle; };
  std::vector<G4double>& Get_pz_particle() { return v_pz_particle; };
  std::vector<G4double>& Get_energyDep_hits() { return v_energyDep_hits; };
  std::vector<G4double>& Get_energyDep_QF_hits() { return v_energyDep_hits_QF; };
  std::vector<G4double>& Get_energyDep_QF_geant_hits() { return v_energyDep_hits_QF_geant; };
  std::vector<G4double>& Get_energyDep_NR_hits() { return v_energyDep_hits_NR; };
  std::vector<G4double>& Get_energyDep_NR_QF_hits() { return v_energyDep_hits_NRQF; };
  std::vector<G4double>& Get_energyDep_NR_QF_geant_hits() { return v_energyDep_hits_NRQF_geant; };
  std::vector<G4double>& Get_x_hits() { return v_x_hits; };
  std::vector<G4double>& Get_y_hits() { return v_y_hits; };
  std::vector<G4double>& Get_z_hits() { return v_z_hits; };

  private: 
  CYGNODetectorConstruction* fDetector;

  // hits collections
  G4int CYGNOID;
  G4int CYGNO_hits;

  std::vector<G4int>    v_pdgID_hits;
  std::vector<G4double> v_tracklen_hits;
  std::vector<G4double> v_px_particle;
  std::vector<G4double> v_py_particle;
  std::vector<G4double> v_pz_particle;
  std::vector<G4double> v_energyDep_hits;
  std::vector<G4double>   v_energyDep_hits_QF;
  std::vector<G4double>   v_energyDep_hits_QF_geant;
  std::vector<G4double>   v_energyDep_hits_NR;
  std::vector<G4double>   v_energyDep_hits_NRQF;
  std::vector<G4double>   v_energyDep_hits_NRQF_geant;
  std::vector<G4double> v_x_hits;
  std::vector<G4double> v_y_hits;
  std::vector<G4double> v_z_hits;
};

#endif


