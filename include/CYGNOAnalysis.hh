#ifndef CYGNOAnalysis_h
#define CYGNOAnalysis_h 1

#include "globals.hh"
#include "g4root.hh"
#include <vector>
#include <string>
#include <unordered_map>
#include <typeinfo>

#include "TFile.h"
#include "TH1F.h"
#include "TTree.h"

#include "G4Track.hh"
#include "G4ThreeVector.hh"
#include "G4LorentzVector.hh"

typedef std::unordered_map<std::string,int> reverse_list;

class G4Run;
class G4Event;
class G4DynamicParticle;
class CYGNODetectorConstruction;

class CYGNOAnalysis
{

public:

  static CYGNOAnalysis* getInstance();
  
private:

  CYGNOAnalysis();


public: 

  ~CYGNOAnalysis();

  void InitRun(G4String FileName, CYGNODetectorConstruction*);
  void EndOfRun();

  void BeginOfEvent(const G4Event *event, CYGNODetectorConstruction*);
  void EndOfEvent(const G4Event *event);

  void SetVerbose(G4int val) {verbose = val;};
  G4int GetVerbose() const {return verbose;};

  void SetOutFileCut(G4int cut) {fOutFileCut = cut;};
  G4int GetOutFileCut() const {return fOutFileCut;};

  void SetRegisterOn(G4int regOn) {fRegisterOn = regOn;};
  G4int GetRegisterOn() const {return fRegisterOn;};

  //void SetTotT(G4int cut) {fTotT = cut;};
  //G4int GetTotT() const {return fTotT;};
    
  void SetHitsInfo(G4int hitsOn) {fHitsInfo = hitsOn;};
  G4int GetHitsInfo() const {return fHitsInfo;};
  
  G4String GetFileName() {return filename;}

  
  void SetNeutronFlag(G4int neuflag) { neutronflag = neuflag; }
  void SetProtonFlag(G4int pflag) { protonflag = pflag; }
  void SetMuonFlag(G4int muflag) { muonflag = muflag; }
  void SetInelasticFlag(G4int inelastic) { inelasticflag = inelastic; }
  G4int GetCopyNo(const G4Track*);
  G4int GetVolNo(const G4Track*);


  void RegisterIsotope(G4int A, G4int Z, G4int PDG, G4double kinE, G4ThreeVector Position, G4int volNo, G4int copyNo, G4int trackID);
  void RegisterParticle(G4int trackID, G4int nextVolNo, G4int nextCopyNo, G4int PDG, G4ThreeVector preStepPt,  G4ThreeVector postStepPt, G4LorentzVector QuadriMomentum);
  void RegisterNeutron(G4int TrackId, G4int ParentId, G4ThreeVector postStepPt, G4LorentzVector QuadriMomentum);
 
  std::vector<std::pair <G4String,G4double> > *vol_name_mass;  
  std::vector<std::pair <G4String,G4double> > *vol_name_dens;  

private:

  // MEMBERS
  static CYGNOAnalysis* fManager;

  G4int fCYGNOID;

  G4int verbose;
  
  G4int fOutFileCut;
  G4int fRegisterOn;
  //G4int fTotT;
  G4int fHitsInfo;
    
  //For Ntuple
  std::vector<std::pair <std::string,std::string> > var_list;//includes all the variables to be added to the output tree and their type 
  std::vector<std::pair <std::string,G4int*> > i_list;
  std::vector<std::pair <std::string,G4double*> > d_list;
  std::vector<std::pair <std::string,G4float*> > f_list;
  std::vector<std::pair <std::string,G4String*> > s_list;
  reverse_list idx_list;//associate the name of the variable to the index of the tree
  //For Histos
  std::vector<std::pair <std::string,G4int*> > hi_list;
  std::vector<std::pair <std::string,G4double*> > hd_list;
  std::vector<std::pair <std::string,G4double> > hrun_list;

  unsigned int NAlwaysFilledHistI;
  unsigned int NAlwaysFilledHistD;
  G4int NTot;


  //  Int_t     eventnumber;
  G4int eventnumber;

  G4int     numvertex;
  std::vector<double>  v_xpos_vertex;
  std::vector<double>  v_ypos_vertex;
  std::vector<double>  v_zpos_vertex;
  std::vector<double>  v_time_vertex;
  std::vector<int>     v_numparticle_vertex;  

  G4int     numparticles;
  std::vector<int>     v_pdgid_particle;
  std::vector<double>  v_ivertex_particle;  // which primary vertex
  std::vector<double>  v_px_particle;
  std::vector<double>  v_py_particle;
  std::vector<double>  v_pz_particle;
  std::vector<double>  v_ekin_particle;
  std::vector<double>  v_etot_particle;

  std::vector<double> v_impact_parameter;
  std::vector<double> v_direc_angle;

  // hits in CYGNO detector
  G4int    numhits;
  
  std::vector<int>    v_pdgID_hits;
  std::vector<G4String>    v_processIni_hits;
  std::vector<G4String>    v_processFin_hits;
  std::vector<int>     v_parentID_hits;
  std::vector<int>     v_trackID_hits;
  std::vector<double>  v_kinEne_hits;
  std::vector<double>  v_time_hits;
  std::vector<double>  v_len_hits;
  std::vector<double>  v_x_hits;
  std::vector<double>  v_y_hits;
  std::vector<double>  v_z_hits;
  std::vector<double>  v_x_vertex_hits;
  std::vector<double>  v_y_vertex_hits;
  std::vector<double>  v_z_vertex_hits;

  std::vector<double>  v_energyDep_hits;
  G4double  energyDep;

  // secondary radionuclides info (former ".iso" file)
  std::vector<G4int>         v_trackid_iso;
  std::vector<G4int>         v_A_iso;
  std::vector<G4int>         v_Z_iso;
  std::vector<G4int>         v_pdg_iso;
  std::vector<G4int>         v_volNo_iso;
  std::vector<G4int>         v_copyNo_iso;
  std::vector<G4double>      v_kinEne_iso;
  std::vector<G4double>      v_x_iso;
  std::vector<G4double>      v_y_iso;
  std::vector<G4double>      v_z_iso;

  // particle flux info (former ".flu" file)
  std::vector<G4int>         v_trackid_flu;
  std::vector<G4int>         v_volNo_flu;
  std::vector<G4int>         v_copyNo_flu;
  std::vector<G4int>         v_pdg_flu;
  std::vector<G4double>      v_prestepX_flu;
  std::vector<G4double>      v_prestepY_flu;
  std::vector<G4double>      v_prestepZ_flu;
  std::vector<G4double>      v_poststepX_flu;
  std::vector<G4double>      v_poststepY_flu;
  std::vector<G4double>      v_poststepZ_flu;
  std::vector<G4double>      v_px_flu;
  std::vector<G4double>      v_py_flu;
  std::vector<G4double>      v_pz_flu;
  std::vector<G4double>      v_E_flu;
  std::vector<G4double>      v_kinE_flu;
  std::vector<G4double>      v_m_flu;
  
  // neutron info (former ".neu" file)
  std::vector<G4int>         v_trackid_neu;
  std::vector<G4int>         v_parentid_neu;
  std::vector<G4double>      v_poststepX_neu;
  std::vector<G4double>      v_poststepY_neu;
  std::vector<G4double>      v_poststepZ_neu;
  std::vector<G4double>      v_px_neu;
  std::vector<G4double>      v_py_neu;
  std::vector<G4double>      v_pz_neu;
  std::vector<G4double>      v_E_neu;
  std::vector<G4double>      v_kinE_neu;
  
  // neutrons
  G4int neutronflag;
  G4int muonflag;
  G4int protonflag;
  G4int inelasticflag;


  G4String filename;
};

  

#endif
