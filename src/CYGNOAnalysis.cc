#include "CYGNOAnalysis.hh"
#include "CYGNOPhysicsList.hh"
#include "CYGNOHit.hh"
#include "CYGNODetectorConstruction.hh"
#include "CYGNODetectorProperty.hh"
#include "globals.hh"
#include "G4UnitsTable.hh"
#include "G4ProcessManager.hh"
#include "G4ProcessVector.hh"
#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4Event.hh"
#include "G4DynamicParticle.hh"
#include "G4PrimaryParticle.hh"
#include "G4PrimaryVertex.hh"
#include "G4Step.hh"
#include "G4ThreeVector.hh"
#include "G4Track.hh"
#include "G4TrackStatus.hh"
#include "G4VProcess.hh"
#include "G4Ions.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleWithCuts.hh"
#include "G4ParticleTypes.hh"
#include "G4ParticleTable.hh"
#include "G4ios.hh"
#include "G4SDManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4VHitsCollection.hh"
#include "G4UnitsTable.hh"
#include <TFile.h>
#include <TH1D.h>
#include <TTree.h>
#include "CYGNOVolumes.hh"
#include <CLHEP/Units/SystemOfUnits.h>
//#include "CYGNOUserEventInformation.hh"

//using namespace CLHEP;

CYGNOAnalysis* CYGNOAnalysis::fManager = 0;

bool comp_nums (int i,int j) { return (i<j); }

//int comp_nums(const void *, const void *);

CYGNOAnalysis* CYGNOAnalysis::getInstance()
{
    if(!fManager) {
        fManager = new CYGNOAnalysis();
    }
    return fManager;
}

    
CYGNOAnalysis::CYGNOAnalysis():
       	fCYGNOID(-1), 
	fHitsInfo(1), 
	fOutFileCut(1), 
	fRegisterOn(1)
	//fTotT(0)
{
    //fMessenger = new CYGNOAnalysisMessenger(this);
}

CYGNOAnalysis::~CYGNOAnalysis()
{
    //  delete G4AnalysisManager::Instance();
    //delete fMessenger;
  delete vol_name_mass;
  delete vol_name_dens;
}

void CYGNOAnalysis::EndOfRun()
{
    
}


void CYGNOAnalysis::InitRun(G4String FileName="out", CYGNODetectorConstruction* Detector=0){
    
    filename = FileName;
    
    G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
    G4cout << "Using " << analysisManager->GetType() << G4endl;
    
    // Default settings
    analysisManager->SetVerboseLevel(1);
    // set root file name
    G4String RootFileName = FileName + ".root";
    analysisManager->SetFileName(RootFileName);
    
    // Create directories inside the output file
    analysisManager->SetHistoDirectoryName("histo");
    
    CYGNODetectorProperty* CYGNOProperties = CYGNODetectorProperty::GetInstance();

    //-------------------------------------------------------------------------------------------
    // Book histograms, ntuple
    //-------------------------------------------------------------------------------------------
    
    
    // Creating histos
    
    
    //These are the histograms always filled, no matter filter applied
    
    //TH1I histograms
    NAlwaysFilledHistI=0;
    hi_list.clear();//list of TH1I
    //G4cout << "hi_list.size() = "<< hi_list.size() << G4endl;
    //G4cout << "NAlwaysFilledHistI = "<< NAlwaysFilledHistI << G4endl;
    analysisManager->CreateH1("NTot","", 1, 0, 1);
    hi_list.push_back(std::make_pair("NTot",&NTot));
    NAlwaysFilledHistI++;
    //G4cout << "hi_list.size() = "<< hi_list.size() << G4endl;
    //G4cout << "NAlwaysFilledHistI = "<< NAlwaysFilledHistI << G4endl;
    
    if(NAlwaysFilledHistI>hi_list.size())
        G4cout << "Error! The number of TH1I histograms to be filled initially ("<< NAlwaysFilledHistI<< ") > of the number of total histograms ("<< hi_list.size() <<")" << G4endl;
    
    //TH1D histograms
    NAlwaysFilledHistD=0;
    hd_list.clear();//list of TH1D
    
    if(NAlwaysFilledHistD>hd_list.size())
        G4cout << "Error! The number of TH1D histograms to be filled initially ("<< NAlwaysFilledHistD <<") > of the number of total histograms ("<< hd_list.size()<<")" << G4endl;
    
    //-------------------------------------------------------------------------------------------
    
    //Create histograms for variables of all events, regardless of filters applied
    
    
    //None so far
    
    //-------------------------------------------------------------------------------------------
    
    //Creating histograms containing run properties: detector settings, run information
    
    hrun_list.clear();
    //Adding boolean histograms. These histograms are created with the name equal to the configuration. If for instance the vessel is built with csg, then an histogram with name csg and 1 bin with content 1 is created. No histrogram with name equal to another vessel design like gdml is created.
    std::string bool_conf;
    G4double bool_bin=1;//DO NOT CHANGE
    bool_conf= Detector->GetCYGNOLab();
    analysisManager->CreateH1(bool_conf,"", 1, 0, 1);
    hrun_list.push_back(std::make_pair(bool_conf,bool_bin));
    
    //Physics volume masses
    vol_name_mass = CYGNOProperties->GetVolumeNameMass();
    for(unsigned int h=0; h<vol_name_mass->size(); h++)
      {
       	G4cout<<"Filling histogram of masses "<<vol_name_mass->at(h).first<<" "<<vol_name_mass->at(h).second/kg<< " kg" << G4endl;
       	analysisManager->CreateH1(vol_name_mass->at(h).first,"", 1, 0, 1);
       	hrun_list.push_back(std::make_pair(vol_name_mass->at(h).first,vol_name_mass->at(h).second/kg));//Mass in kg
        //NAlwaysFilledHistD++;
      }


    //Physics volume densities
    vol_name_dens = CYGNOProperties->GetVolumeNameDensity();
    for(unsigned int h=0; h<vol_name_dens->size(); h++)
      {
    	G4cout<<"Filling histogram of densities "<<vol_name_dens->at(h).first<<" "<<vol_name_dens->at(h).second*m*m*m/kg << " kg/m^3" << G4endl;
    	analysisManager->CreateH1(vol_name_dens->at(h).first,"", 1, 0, 1);
        hrun_list.push_back(std::make_pair(vol_name_dens->at(h).first,vol_name_dens->at(h).second*m*m*m/kg));//Density in kg/m^3
        //NAlwaysFilledHistD++;
      }


    //Filling the run info histograms
    for(unsigned int h=0; h<hrun_list.size(); h++)
        analysisManager->FillH1(h+hi_list.size()+hd_list.size(), 0., hrun_list.at(h).second);
    
    
    ////-------------------------------------------------------------------------------------------
    //
    ////Ntuple always filled, no matter filter applied
    //analysisManager->CreateNtuple("TotT", "treeD");
    //
    ////These branches are all vectors. See what done with the tree nTuple if you want to add sigle variables
    //analysisManager->CreateNtupleDColumn("EKinTot",v_ekin_particle);
    //analysisManager->CreateNtupleDColumn("ImParaTot",v_impact_parameter);
    //analysisManager->CreateNtupleDColumn("DiAngleTot",v_direc_angle);
    //analysisManager->CreateNtupleDColumn("XPosTot",v_xpos_vertex);
    //analysisManager->CreateNtupleDColumn("YPosTot",v_ypos_vertex);
    //analysisManager->CreateNtupleDColumn("ZPosTot",v_zpos_vertex);
    //analysisManager->CreateNtupleDColumn("PxTot",v_px_particle);
    //analysisManager->CreateNtupleDColumn("PyTot",v_py_particle);
    //analysisManager->CreateNtupleDColumn("PzTot",v_pz_particle);
    //analysisManager->FinishNtuple();
    
    
    //-------------------------------------------------------------------------------------------
    
    //Ntuple filled when the filter selection is passed
    //
    analysisManager->CreateNtuple("nTuple", "tree");
    
    //Here define the variables that go in the tree.
    
    //This is how you add the single parameters. The vector branches are initialized separately, see below
    i_list.clear();
    i_list.push_back(std::make_pair("eventnumber",&eventnumber));
    i_list.push_back(std::make_pair("numvertex",&numvertex));
    i_list.push_back(std::make_pair("numparticles",&numparticles));
    i_list.push_back(std::make_pair("numhits",&numhits));
    i_list.push_back(std::make_pair("neutronflag",&neutronflag));
    i_list.push_back(std::make_pair("muonflag",&muonflag));
    i_list.push_back(std::make_pair("inelasticflag",&inelasticflag));
    d_list.clear();
    d_list.push_back(std::make_pair("energyDep",&energyDep));


    f_list.clear();
    s_list.clear();
    
    for(unsigned int i=0; i<i_list.size(); i++)
        analysisManager->CreateNtupleIColumn(i_list.at(i).first);
    for(unsigned int i=0; i<d_list.size(); i++)
        analysisManager->CreateNtupleDColumn(d_list.at(i).first);
    for(unsigned int i=0; i<f_list.size(); i++)
        analysisManager->CreateNtupleFColumn(f_list.at(i).first);
    for(unsigned int i=0; i<s_list.size(); i++)
        analysisManager->CreateNtupleSColumn(s_list.at(i).first);
    
    
    //-------------------------------------------------------------------------------------------
    
    //Add vector branches to the ntuple

    
    
    // MC true primary particle information
    analysisManager->CreateNtupleDColumn("xpos_vertex",v_xpos_vertex);
    analysisManager->CreateNtupleDColumn("ypos_vertex",v_ypos_vertex);
    analysisManager->CreateNtupleDColumn("zpos_vertex",v_zpos_vertex);
    analysisManager->CreateNtupleDColumn("time_vertex",v_time_vertex);
    analysisManager->CreateNtupleIColumn("numparticle_vertex",v_numparticle_vertex);
    analysisManager->CreateNtupleIColumn("pdgid_particle",v_pdgid_particle);
    analysisManager->CreateNtupleDColumn("ivertex_particle",v_ivertex_particle);
    analysisManager->CreateNtupleDColumn("px_particle",v_px_particle);
    analysisManager->CreateNtupleDColumn("py_particle",v_py_particle);
    analysisManager->CreateNtupleDColumn("pz_particle",v_pz_particle);
    analysisManager->CreateNtupleDColumn("ekin_particle",v_ekin_particle);
    analysisManager->CreateNtupleDColumn("etot_particle",v_etot_particle);
    analysisManager->CreateNtupleDColumn("impact_param_particle",v_impact_parameter);
    analysisManager->CreateNtupleDColumn("direc_angle_particle",v_direc_angle);
    // hits in detector
    //if(fHitsInfo)
      
      analysisManager->CreateNtupleIColumn("parentID_hits",v_parentID_hits);//It is the generator parent track ID
      analysisManager->CreateNtupleIColumn("pdgID_hits",v_pdgID_hits);
      analysisManager->CreateNtupleIColumn("trackID_hits",v_trackID_hits);//It is the generator track ID
      analysisManager->CreateNtupleDColumn("time_hits",v_time_hits);
      analysisManager->CreateNtupleDColumn("kinEne_hits",v_kinEne_hits);
      //analysisManager->CreateNtupleSColumn("processIni_hits",v_processIni_hits);
      //analysisManager->CreateNtupleSColumn("processFin_hits",v_processFin_hits);
      analysisManager->CreateNtupleDColumn("x_hits",v_x_hits);
      analysisManager->CreateNtupleDColumn("y_hits",v_y_hits);
      analysisManager->CreateNtupleDColumn("z_hits",v_z_hits);
      analysisManager->CreateNtupleDColumn("x_vertex_hits",v_x_vertex_hits);
      analysisManager->CreateNtupleDColumn("y_vertex_hits",v_y_vertex_hits);
      analysisManager->CreateNtupleDColumn("z_vertex_hits",v_z_vertex_hits);
      analysisManager->CreateNtupleDColumn("tracklen_hits",v_len_hits);
      analysisManager->CreateNtupleDColumn("energyDep_hits",v_energyDep_hits);
      
      
    if(fRegisterOn){
      // secondary radionuclides (former ".iso" file)
      analysisManager->CreateNtupleIColumn("trackid_iso",v_trackid_iso);
      analysisManager->CreateNtupleIColumn("A_iso",v_A_iso);// Mass number of the isotope
      analysisManager->CreateNtupleIColumn("Z_iso",v_Z_iso);// Atomic number of the isotope
      analysisManager->CreateNtupleIColumn("pdg_iso",v_pdg_iso);
      analysisManager->CreateNtupleIColumn("volNo_iso",v_volNo_iso); // ID number of the volume of generation of the radionuclide --> check definitions in CYGNOVolumes.hh
      analysisManager->CreateNtupleIColumn("copyNo_iso",v_copyNo_iso); // copy number of the detector of generation of the radionuclide
      analysisManager->CreateNtupleDColumn("kinEne_iso",v_kinEne_iso);
      analysisManager->CreateNtupleDColumn("posx_iso",v_x_iso);
      analysisManager->CreateNtupleDColumn("posy_iso",v_y_iso);
      analysisManager->CreateNtupleDColumn("posz_iso",v_z_iso);
      
      // particle flux info (former ".flu" file)
      analysisManager->CreateNtupleIColumn("trackid_flu",v_trackid_flu);
      analysisManager->CreateNtupleIColumn("volNo_flu",v_volNo_flu); // number of the volume in which the particle is entering
      analysisManager->CreateNtupleIColumn("copyNo_flu",v_copyNo_flu); // copy number of the volume in which the particle is entering
      analysisManager->CreateNtupleIColumn("pdg_flu",v_pdg_flu);
      analysisManager->CreateNtupleDColumn("prestepx_flu",v_prestepX_flu);
      analysisManager->CreateNtupleDColumn("prestepy_flu",v_prestepY_flu);
      analysisManager->CreateNtupleDColumn("prestepz_flu",v_prestepZ_flu);
      analysisManager->CreateNtupleDColumn("poststepx_flu",v_poststepX_flu);
      analysisManager->CreateNtupleDColumn("poststepy_flu",v_poststepY_flu);
      analysisManager->CreateNtupleDColumn("poststepz_flu",v_poststepZ_flu);
      analysisManager->CreateNtupleDColumn("px_flu",v_px_flu);
      analysisManager->CreateNtupleDColumn("py_flu",v_py_flu);
      analysisManager->CreateNtupleDColumn("pz_flu",v_pz_flu);
      analysisManager->CreateNtupleDColumn("E_flu",v_E_flu);
      analysisManager->CreateNtupleDColumn("kinE_flu",v_kinE_flu);
      analysisManager->CreateNtupleDColumn("m_flu",v_m_flu);
      
      // neutron info (former ".neu" file)
      analysisManager->CreateNtupleIColumn("trackid_neu",v_trackid_neu);
      analysisManager->CreateNtupleIColumn("parentid_neu",v_parentid_neu);
      analysisManager->CreateNtupleDColumn("poststepx_neu",v_poststepX_neu);
      analysisManager->CreateNtupleDColumn("poststepy_neu",v_poststepY_neu);
      analysisManager->CreateNtupleDColumn("poststepz_neu",v_poststepZ_neu);
      analysisManager->CreateNtupleDColumn("px_neu",v_px_neu);
      analysisManager->CreateNtupleDColumn("py_neu",v_py_neu);
      analysisManager->CreateNtupleDColumn("pz_neu",v_pz_neu);
      analysisManager->CreateNtupleDColumn("E_neu",v_E_neu);
      analysisManager->CreateNtupleDColumn("kinE_neu",v_kinE_neu);
    }
    analysisManager->FinishNtuple();
}

void CYGNOAnalysis::BeginOfEvent(const G4Event *event, CYGNODetectorConstruction* Detector)
{
    
    G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
    //New event, add the user information object
    //G4EventManager::GetEventManager()->SetUserInformation(new CYGNOUserEventInformation);
    
    //eventnumber = G4ToRoot(event->GetEventID());
    eventnumber = event->GetEventID();
    NTot = 0;
    //hitograms always filled, no matter filter applied
    for(unsigned int h = 0; h < NAlwaysFilledHistI; h++)
        analysisManager->FillH1(h, *(hi_list.at(h).second));
    for(unsigned int h = 0; h < NAlwaysFilledHistD; h++)
        analysisManager->FillH1(h+hi_list.size(), *(hd_list.at(h).second));
    
    G4PrimaryParticle *primaryParticle;
    G4PrimaryVertex   *primaryVertex;
    
    primaryVertex = event->GetPrimaryVertex(0);
    primaryParticle = primaryVertex->GetPrimary(0);
    
    numvertex = event->GetNumberOfPrimaryVertex();
    numparticles=0;
    neutronflag = 0;
    muonflag = 0;
    inelasticflag = 0;
    

    //cleaning vectors
    vol_name_mass->clear();
    vol_name_dens->clear();
    v_xpos_vertex.clear();
    v_ypos_vertex.clear();
    v_zpos_vertex.clear();
    v_time_vertex.clear();
    v_numparticle_vertex.clear();
    v_pdgid_particle.clear();
    v_ivertex_particle.clear();
    v_px_particle.clear();
    v_py_particle.clear();
    v_pz_particle.clear();
    v_ekin_particle.clear();
    v_etot_particle.clear();
    v_impact_parameter.clear();
    v_direc_angle.clear();
    
    v_processIni_hits.clear();
    v_processFin_hits.clear();
    v_parentID_hits.clear();
    v_trackID_hits.clear();
    v_kinEne_hits.clear();
    v_time_hits.clear();
    v_x_hits.clear();
    v_y_hits.clear();
    v_z_hits.clear();
    v_x_vertex_hits.clear();
    v_y_vertex_hits.clear();
    v_z_vertex_hits.clear();
    v_energyDep_hits.clear();
    v_len_hits.clear();    


    v_A_iso.clear();
    v_Z_iso.clear();
    v_pdg_iso.clear();
    v_kinEne_iso.clear();
    v_x_iso.clear();
    v_y_iso.clear();
    v_z_iso.clear();
    v_volNo_iso.clear();
    v_copyNo_iso.clear();
    v_trackid_iso.clear();
    
    v_volNo_flu.clear();
    v_copyNo_flu.clear();
    v_pdg_flu.clear();
    v_prestepX_flu.clear();
    v_prestepY_flu.clear();
    v_prestepZ_flu.clear();
    v_poststepX_flu.clear();
    v_poststepY_flu.clear();
    v_poststepZ_flu.clear();
    v_px_flu.clear();
    v_py_flu.clear();
    v_pz_flu.clear();
    v_E_flu.clear();
    v_kinE_flu.clear();
    v_m_flu.clear();
    v_trackid_flu.clear();
    
    v_trackid_neu.clear();
    v_parentid_neu.clear();
    v_poststepX_neu.clear();
    v_poststepY_neu.clear();
    v_poststepZ_neu.clear();
    v_px_neu.clear();
    v_py_neu.clear();
    v_pz_neu.clear();
    v_E_neu.clear();
    v_kinE_neu.clear();
   
     
    
    //Position of detector's centre for imapct parameter calculation
    G4ThreeVector x0;
    //x0 = Detector->GetTranslation();
    
    //Position of event vertex (x, y, z)
    G4ThreeVector x1;
    //Another vector in the direction of the event particle's travel: (x, y, z) + (px, py, pz)
    G4ThreeVector x2;
    
    double x_vertex;
    double y_vertex;
    double z_vertex;
    
    double px_particle;
    double py_particle;
    double pz_particle;
    
    for(G4int i = 0; i < numvertex; i++) {
        primaryVertex = event->GetPrimaryVertex(i);
        
        x_vertex = primaryVertex->GetX0();
        y_vertex = primaryVertex->GetY0();
        z_vertex = primaryVertex->GetZ0();
        
        v_xpos_vertex.push_back(x_vertex);
        v_ypos_vertex.push_back(y_vertex);
        v_zpos_vertex.push_back(z_vertex);
        v_time_vertex.push_back(primaryVertex->GetT0());
        v_numparticle_vertex.push_back(primaryVertex->GetNumberOfParticle());
        
        for (G4int j=0; j<v_numparticle_vertex.at(i); j++) {
            primaryParticle = primaryVertex->GetPrimary(j);
            
            px_particle = primaryParticle->GetPx();
            py_particle = primaryParticle->GetPy();
            pz_particle = primaryParticle->GetPz();
            
            v_pdgid_particle.push_back(primaryParticle->GetPDGcode());
            v_px_particle.push_back(px_particle);
            v_py_particle.push_back(py_particle);
            v_pz_particle.push_back(pz_particle);
            
            v_ekin_particle.push_back(primaryParticle->GetKineticEnergy());
            v_etot_particle.push_back(primaryParticle->GetTotalEnergy());
            v_ivertex_particle.push_back(i);
            
            x1 = G4ThreeVector(x_vertex, y_vertex, z_vertex);
            x2 = G4ThreeVector(x_vertex+px_particle, y_vertex+py_particle, z_vertex+pz_particle);
            
            v_impact_parameter.push_back((((x0-x1).cross(x0-x2)).mag())/((x2-x1).mag()));
            v_direc_angle.push_back((x2-x1).angle(x0-x1));
            
            numparticles++;
            
        }
    }
    
}

void CYGNOAnalysis::RegisterIsotope(G4int A, G4int Z, G4int PDG, G4double kinE, G4ThreeVector Position, G4int volNo, G4int copyNo, G4int trackID)
{
  if(fRegisterOn)
    {
      v_trackid_iso.push_back(trackID);
      v_volNo_iso.push_back(volNo);
      v_copyNo_iso.push_back(copyNo);
      v_A_iso.push_back(A);
      v_Z_iso.push_back(Z);
      v_pdg_iso.push_back(PDG);
      v_kinEne_iso.push_back(kinE/keV);
      v_x_iso.push_back(Position.x()/mm);
      v_y_iso.push_back(Position.y()/mm);
      v_z_iso.push_back(Position.z()/mm);
    }
}

void CYGNOAnalysis::RegisterParticle(G4int trackID, G4int volNo, G4int copyNo, G4int PDG, G4ThreeVector preStepPt,  G4ThreeVector postStepPt, G4LorentzVector QuadriMomentum)
{
  if(fRegisterOn) 
    {
      v_trackid_flu.push_back(trackID);
      v_volNo_flu.push_back(volNo);
      v_copyNo_flu.push_back(copyNo);
      v_pdg_flu.push_back(PDG);
      v_prestepX_flu.push_back(preStepPt.x()/mm);
      v_prestepY_flu.push_back(preStepPt.y()/mm);
      v_prestepZ_flu.push_back(preStepPt.z()/mm);
      v_poststepX_flu.push_back(postStepPt.x()/mm);
      v_poststepY_flu.push_back(postStepPt.y()/mm);
      v_poststepZ_flu.push_back(postStepPt.z()/mm);
      v_px_flu.push_back(QuadriMomentum.px()/keV);
      v_py_flu.push_back(QuadriMomentum.py()/keV);
      v_pz_flu.push_back(QuadriMomentum.pz()/keV);
      v_E_flu.push_back(QuadriMomentum.e()/keV);
      G4double kinE = (QuadriMomentum.e()/keV - sqrt(QuadriMomentum.m2())/keV);
      v_kinE_flu.push_back(kinE);
      v_m_flu.push_back(sqrt(QuadriMomentum.m2())/keV);
    }
}

void CYGNOAnalysis::RegisterNeutron(G4int TrackId, G4int ParentId, G4ThreeVector postStepPt, G4LorentzVector QuadriMomentum)
{

  if(fRegisterOn)
    {
      v_trackid_neu.push_back(TrackId);
      v_parentid_neu.push_back(ParentId);
      v_poststepX_neu.push_back(postStepPt.x()/mm);
      v_poststepY_neu.push_back(postStepPt.y()/mm);
      v_poststepZ_neu.push_back(postStepPt.z()/mm);
      v_px_neu.push_back(QuadriMomentum.px()/keV);
      v_py_neu.push_back(QuadriMomentum.py()/keV);
      v_pz_neu.push_back(QuadriMomentum.pz()/keV);
      v_E_neu.push_back(QuadriMomentum.e()/keV);
      G4double kinE = (QuadriMomentum.e()/keV - sqrt(QuadriMomentum.m2())/keV);
      v_kinE_neu.push_back(kinE);
    }
}

void CYGNOAnalysis::EndOfEvent(const G4Event *event)
{
    
    G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
    //CYGNOUserEventInformation* eventInformation
//	      =(CYGNOUserEventInformation*)event->GetUserInformation();
    
    if (fCYGNOID==-1) {
        G4SDManager * SDman = G4SDManager::GetSDMpointer();
        fCYGNOID = SDman->GetCollectionID("CYGNOCollection");
 
    }
       
    G4HCofThisEvent* HCE = event->GetHCofThisEvent();
    CYGNOHitsCollection* CYGNOHC = 0;
    
    CYGNOAnalysis* analysis = CYGNOAnalysis::getInstance() ;

    
    if (HCE) {
        if (fCYGNOID != -1) CYGNOHC = static_cast<CYGNOHitsCollection*>(HCE->GetHC(fCYGNOID));
    }
    
    energyDep=0.;
    G4ThreeVector tempvec;


    if (CYGNOHC)
    {
         numhits = CYGNOHC->entries();  // number of hits
        if(numhits >0)
        {
            // fill hit_array with the id of the detectors hit
            for (G4int i=0; i<numhits; i++)
            {
            if(fHitsInfo){
                    v_parentID_hits.push_back((int)(*CYGNOHC)[i]->GetParentID());//It is the generator parent track ID
                    v_pdgID_hits.push_back((int)(*CYGNOHC)[i]->GetParticleID());
                    v_trackID_hits.push_back((int)(*CYGNOHC)[i]->GetTrackID());//It is the generator track ID
                    v_time_hits.push_back((*CYGNOHC)[i]->GetGlobalTime());
                    v_kinEne_hits.push_back((*CYGNOHC)[i]->GetKineticEne());
                    v_processIni_hits.push_back((*CYGNOHC)[i]->GetProcessIni());
                    v_processFin_hits.push_back((*CYGNOHC)[i]->GetProcessFin());
                    v_x_vertex_hits.push_back(v_xpos_vertex[0]);
                    v_y_vertex_hits.push_back(v_ypos_vertex[0]);
                    v_z_vertex_hits.push_back(v_zpos_vertex[0]);
                    tempvec = (*CYGNOHC)[i]->GetPos();
                    v_x_hits.push_back(tempvec.getX());
                    v_y_hits.push_back(tempvec.getY());
                    v_z_hits.push_back(tempvec.getZ());
                    v_len_hits.push_back((*CYGNOHC)[i]->GetLength());
	    }

	    v_energyDep_hits.push_back((*CYGNOHC)[i]->GetEdep());
	    
	    // sum total energy deposited in hits
	    energyDep += (*CYGNOHC)[i]->GetEdep();
	    }
    	}
    } 
  
    //Now filling histograms for any event
    for(unsigned int h=NAlwaysFilledHistI; h<hi_list.size(); h++)
        analysisManager->FillH1(h, *(hi_list.at(h).second));
    for(unsigned int h=NAlwaysFilledHistD; h<hd_list.size(); h++)
        analysisManager->FillH1(h+hi_list.size(), *(hd_list.at(h).second));
    
    //Fill all event (filter blind) histograms in TotT tree
    //if(fTotT)
    //    analysisManager->AddNtupleRow(0);
    
    //Selection cuts to be applied. The tree will be filled only if the event passes the selection
    if(fOutFileCut == 0 || (fOutFileCut == 1 && energyDep>0.) )
    {
        //Now filling the ntuple
        for(unsigned int i=0; i<i_list.size(); i++)
            analysisManager->FillNtupleIColumn(0, i,*(i_list.at(i).second));
        for(unsigned int i=0; i<d_list.size(); i++)
            analysisManager->FillNtupleDColumn(0, i+i_list.size(),*(d_list.at(i).second));
        for(unsigned int i=0; i<f_list.size(); i++)
            analysisManager->FillNtupleFColumn(0, i+i_list.size()+d_list.size(),*(f_list.at(i).second));
        for(unsigned int i=0; i<s_list.size(); i++)
            analysisManager->FillNtupleSColumn(0, i+i_list.size()+d_list.size()+f_list.size(),*(s_list.at(i).second));
        
        analysisManager->AddNtupleRow(0);
    }  
}


G4int CYGNOAnalysis::GetCopyNo(const G4Track* track)
{
    G4int copyNo = UNKNOWN;

    if(track->GetNextVolume())
        copyNo = track->GetNextVolume()->GetCopyNo();

    return copyNo;
}

G4int CYGNOAnalysis::GetVolNo(const G4Track* track)
{
    G4int volNo = UNKNOWN;
    G4String PVname;
    if(track->GetNextVolume())
        PVname = track->GetNextVolume()->GetName();

    if(PVname=="WorldVolume") volNo = WORLD;
    else if(PVname=="externalRock") volNo = ROCK;
    else if(PVname=="expHall") volNo = HALLB;
    else if(PVname=="Shield0") volNo = SHIELD0;
    else if(PVname=="Shield1") volNo = SHIELD1;
    else if(PVname=="Shield2") volNo = SHIELD2;
    else if(PVname=="Shield3") volNo = SHIELD3;
    else if(PVname=="AirBox") volNo = AIRBOX;
    else if(PVname=="cad_shell_physical") volNo = ACRYLICSHELL;
    else if(PVname=="CYGNO_gas") volNo = CYGNOGAS;

    return volNo;
}

