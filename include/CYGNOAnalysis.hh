#ifndef CYGNOANALYSIS_HH
#define CYGNOANALYSIS_HH

#include "G4AnalysisManager.hh"
#include "G4Threading.hh"
#include "G4ThreeVector.hh"
#include "G4String.hh"

class CYGNOAnalysis {
public:
    static CYGNOAnalysis* Instance();

    void BookNtuples();
    void OpenFile(const G4String& filename);
    void CloseFile();
    void SetOutFile(G4String fname) {FileName = fname;};

private:
    CYGNOAnalysis() = default;
    ~CYGNOAnalysis() = default;

    CYGNOAnalysis(const CYGNOAnalysis&) = delete;
    CYGNOAnalysis& operator=(const CYGNOAnalysis&) = delete;
    G4String FileName;
};

#endif

