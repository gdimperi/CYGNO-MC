#ifndef CYGNODetectorMaterial_h
#define CYGNODetectorMaterial_h 1

#include "G4Material.hh"
#include "G4Element.hh"
#include "G4NistManager.hh"
#include "G4VisAttributes.hh"
#include "globals.hh"

using namespace std;

class G4OpticalSurface;
//class CYGNODetectorMaterialMessenger;

class CYGNODetectorMaterial
{
public:

  CYGNODetectorMaterial();
  ~CYGNODetectorMaterial();
  static CYGNODetectorMaterial* GetInstance();
  G4Material* Material(G4String);
  G4VisAttributes* VisAttributes(G4String);
  G4Material* FindOrBuildMaterial(const G4String& name, G4bool isotopes=true, G4bool warning=false);
  static CYGNODetectorMaterial* fCYGNODetectorMaterial;
  void Refresh();
  void ConstructMaterials();
  void UseOpticalVeto(bool flag);
  void UseReflectionsVeto(bool flag);
  void UseOpticalCrystal(bool flag);
  void UseReflectionsCrystal(bool flag);
  void UsePC(bool flag);
  bool IsUsingOpticalVeto();
  bool IsUsingReflectionsVeto();
  bool IsUsingOpticalCrystal();
  bool IsUsingReflectionsCrystal();
  bool IsUsingPC();
  G4OpticalSurface* GetSurf(G4String);
  void PrintSurfs();
  void SetPCyield(double PCLY);
  void SetNaIyield(double NaILY);
  void SetLAByield(double LABLY);


private:

//**********************************************************************
//   NAME OF MATERIALS
//**********************************************************************

G4NistManager* man;

G4Material* Air;

G4Material* O;
G4Material* Na;
G4Material* K;
G4Material* B;
G4Material* Al;
G4Material* Cu;
G4Material* Fe;
G4Material* Pb;
G4Material* Co;
G4Material* Ni;
G4Material* Si;
G4Material* In;

G4Material* Teflon;
G4Material* PyrexGlass;
G4Material* BSglass;
G4Material* VetoPMTglass;
G4Material* Quartz;
G4Material* Ceramic;
G4Material* Kovar;
G4Material* lngsRock;
G4Material* Vacuum;
G4Material* Water;
G4Material* Steel;
G4Material* PE;
G4Material* Concrete;
G4Material* Camera;
G4Material* Perspex;
G4Material* CYGNO_gas;
G4Material* SF6_gas;
G4Material* He_gas;
G4Material* CF4_gas; 
G4Material* Kapton; 


G4VisAttributes* PEVis;
G4VisAttributes* PbVis;
G4VisAttributes* WaterVis;
G4VisAttributes* AirVis;
G4VisAttributes* VacuumVis;
G4VisAttributes* CopperVis;
G4VisAttributes* PerspexVis;
G4VisAttributes* CameraVis;
G4VisAttributes* ConcreteVis;
G4VisAttributes* LNGSRockVis;
G4VisAttributes* CYGNOGasVis;

//CYGNODetectorMaterialMessenger* fMessenger;


};

#endif   /* CYGNODetectorMaterial.hh */

