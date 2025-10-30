#include "CYGNOAnalysis.hh"
#include "G4AnalysisManager.hh"
#include "G4Threading.hh"

    // Singleton accessor
    CYGNOAnalysis* CYGNOAnalysis::Instance() {
        static CYGNOAnalysis instance;
        return &instance;
    }

    void CYGNOAnalysis::BookNtuples() {
        //auto man = G4AnalysisManager::Instance();
        G4AnalysisManager* man = G4AnalysisManager::Instance();
        man->SetDefaultFileType("root");
         
         man->SetNtupleMerging(true);

         man->SetFirstNtupleId(1);


         // ---- primary ntuple ------
         // id==1
         man->CreateNtuple("tree1", "Particle Source Info");
         man->CreateNtupleDColumn("Energy");
         man->CreateNtupleDColumn("xpos_vertex");
         man->CreateNtupleDColumn("ypos_vertex");
         man->CreateNtupleDColumn("zpos_vertex");
         man->FinishNtuple();

         // ---- secondary ntuple ------   
         //id==2
         man->CreateNtuple("nTuple", "Hits Info");
         man->CreateNtupleIColumn("eventnumber");
         man->CreateNtupleIColumn("numhits");
         man->CreateNtupleDColumn("ekin_particle");
         man->CreateNtupleIColumn("particle_type");
         man->CreateNtupleDColumn("energyDep");
         man->CreateNtupleDColumn("energyDep_QF");
         man->CreateNtupleDColumn("energyDep_NR");
         man->CreateNtupleDColumn("energyDep_NRQF");
         man->CreateNtupleDColumn("energyDep_NRQF_geant");
         man->CreateNtupleIColumn("pdgID_hits");
         man->CreateNtupleDColumn("tracklen_hits");
         man->CreateNtupleDColumn("px_particle");
         man->CreateNtupleDColumn("py_particle");
         man->CreateNtupleDColumn("pz_particle");
         man->CreateNtupleDColumn("energyDep_hits");
         man->CreateNtupleDColumn("energyDep_QF_hits");
         man->CreateNtupleDColumn("energyDep_QF_geant_hits");
         man->CreateNtupleDColumn("energyDep_NR_hits");
         man->CreateNtupleDColumn("energyDep_NRQF_hits");
         man->CreateNtupleDColumn("energyDep_NRQF_geant_hits");
         man->CreateNtupleDColumn("x_hits");
         man->CreateNtupleDColumn("y_hits");
         man->CreateNtupleDColumn("z_hits");
         man->FinishNtuple();
    }

    void CYGNOAnalysis::OpenFile(const G4String& filename) {
        // Get/create analysis manager
        auto man = G4AnalysisManager::Instance();
        if (G4Threading::IsMasterThread()) {
         // Open an output file
         //FIXME -->  take the out name from run action messenger
	 //man->OpenFile(FileName);
	 G4cout << "####### Opening output file" << G4endl;
         man->OpenFile("output.root");
        }
    }
    void CYGNOAnalysis::CloseFile() {
        auto man = G4AnalysisManager::Instance();
        if (G4Threading::IsMasterThread()) {
            man->Write();
	    G4cout << "####### Writing output file" << G4endl;
            man->CloseFile();
	    G4cout << "######### Closing output file" << G4endl;
        }
    }

