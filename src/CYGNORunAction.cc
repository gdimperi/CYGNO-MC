#include "CYGNORunAction.hh"
#include "CYGNORunActionMessenger.hh"
#include "CYGNOAnalysis.hh"
#include "CYGNODetectorConstruction.hh"

#include "CYGNORun.hh"
#include "G4RunManager.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4IonConstructor.hh"
#include "G4BosonConstructor.hh"
#include "G4LeptonConstructor.hh"
#include "G4MesonConstructor.hh"
#include "G4BaryonConstructor.hh"
#include "G4ShortLivedConstructor.hh"
#include "G4ProcessManager.hh"
#include "G4ProcessVector.hh"
#include "G4ionIonisation.hh"
#include "G4SystemOfUnits.hh"
#include "G4LossTableManager.hh"
#include "G4EmCalculator.hh"
#include "G4SystemOfUnits.hh"
#include "G4Alpha.hh"
#include "G4Material.hh"
#include "G4MaterialTable.hh"
#include "CYGNODetectorMaterial.hh"

CYGNORunAction::CYGNORunAction(CYGNODetectorConstruction* myDC)
  : G4UserRunAction(), fDetector(myDC) 
{

  fMessenger = new CYGNORunActionMessenger(this);
  fNumberOfMuonsDetector = 0;
  FileName = "out";
  fOutFileCut = 1;
  fRegisterOn = 1;
  fHitsInfo = 1;


}

CYGNORunAction::~CYGNORunAction()
{
  delete G4AnalysisManager::Instance();
  CYGNOAnalysis* analysis = CYGNOAnalysis::getInstance();
  delete analysis;//??? Is it the place where to do it???
  //delete fMessenger;
}

void CYGNORunAction::BeginOfRunAction(const G4Run* aRun)
{
	G4cout << "### Run " << aRun->GetRunID() << " start." << G4endl;

//  FIXME: modify step function for ion ionization low energy
    //    //Low energy ion
    //    auto particleTable = G4ParticleTable::GetParticleTable();
    //    G4ParticleTable::G4PTblDicIterator* it = particleTable->GetIterator();
    //    it->reset();  
    //    while ((*it)()) {
    //        G4ParticleDefinition* particle = it->value();
    //        G4ProcessManager* pmanager = particle->GetProcessManager();
    //        G4String particleName = particle->GetParticleName();
    //        G4String particleType = particle->GetParticleType();
    //        //G4cout << "Particle: " << particleName 
    //        //   << ", Type: " << particleType 
    //        //   << ", ProcessManager: " << (pmanager ? "OK" : "NULL") 
    //        //   << G4endl;
    //  
    //        // For ions (GenericIon, alpha, etc.)
    //        if (particleType == "nucleus")  {
    //            // Modify existing ionIonisation process
    //  	        if (!pmanager) {
    //  	            G4cerr << "WARNING: No process manager for particle " << particleName << G4endl;
    //  	            continue; // skip this particle
    //  	        }
    //            else{
    //                G4ProcessVector* processList = pmanager->GetProcessList();
    //                for (size_t i = 0; i < processList->size(); ++i) {
    //                    auto* proc = (*processList)[i];
    //                    auto* ionIoni = dynamic_cast<G4ionIonisation*>(proc);
    //                    if (ionIoni) {
    //                        G4cout << " Customizing ionIoni for " << particleName << G4endl;
    //                        ionIoni->SetStepFunction(1e-3, 1.*um);
    //                        ionIoni->SetLinearLossLimit(1e-3);
    //    		    ionIoni->SetMinKinEnergy(0.01 * keV);  // Energy threshold  
    //                    }
    //                }
    //    
    //  	  	}
    //        }
    //    }
 // print dEdx table 
    G4EmCalculator emCal;

    CYGNODetectorMaterial* matBuilder = new CYGNODetectorMaterial();
    matBuilder->ConstructMaterials();

    //G4ParticleDefinition* particle = G4Alpha::AlphaDefinition();
    //G4Material* material = G4Material::GetMaterial("CYGNO_gas");
    G4String particleName = "alpha";
    G4String materialName = "CYGNO_gas";  // Must match your material name
    const G4ParticleDefinition* alpha = G4Alpha::Alpha();
    const G4Material* mat = G4Material::GetMaterial("CYGNO_gas");
    G4String processName = "ionIoni";  // or "msc", "eIoni", etc.

    if (!mat) {
        G4cerr << "Material CYGNO_gas not found!" << G4endl;
        return;
    }

    G4cout << "\nEnergy (keV)\tStopping Power (MeV/mm)\t CSDA range (mm)" << G4endl;
    for (double E_keV = 0.1; E_keV <= 10.0; E_keV += 0.1)
    {
        G4double energy = E_keV * keV;
        G4double dedx = emCal.ComputeDEDX(energy, alpha, processName, mat);
        G4double range = emCal.GetCSDARange(energy , alpha, mat);
        G4cout << E_keV << "\t\t" << dedx / (MeV / mm) <<"\t\t"<< range / mm <<  G4endl;

    }


	fNumberOfMuonsDetector = 0;

	//inform the runManager to save random number seed
	//G4RunManager::GetRunManager()->SetRandomNumberStore(true);
	
	// Get analysis manager
	G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();


	CYGNOAnalysis* analysis = CYGNOAnalysis::getInstance();
	analysis->SetOutFileCut(fOutFileCut);
	analysis->SetRegisterOn(fRegisterOn);
	analysis->SetHitsInfo(fHitsInfo);
	analysis->InitRun(FileName,fDetector);
	// Open an output file
	// The default file name is set in RunAction::RunAction(),
	// it can be overwritten in a macro
	analysisManager->OpenFile();
}

void CYGNORunAction::EndOfRunAction(const G4Run*)
{ 
  
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  analysisManager->Write();
  analysisManager->CloseFile();

}

G4Run* CYGNORunAction::GenerateRun() {
    return new CYGNORun;
}

