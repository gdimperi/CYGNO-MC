#include "CYGNODetectorProperty.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "G4Material.hh"

CYGNODetectorProperty* CYGNODetectorProperty::fCYGNODetectorProperty = NULL;

CYGNODetectorProperty::CYGNODetectorProperty() :
  tolerance(0.1*mm),
  debugmode(0)
{
  //mass and name of the physical volumes
  vol_name_mass = new std::vector<std::pair <G4String,G4double> >();
  vol_name_mass->clear();
  //density and name of the physical volumes
  vol_name_dens = new std::vector<std::pair <G4String,G4double> >();
  vol_name_dens->clear();
}

CYGNODetectorProperty::~CYGNODetectorProperty()
{
  //  delete vol_name_mass;
  //  delete vol_name_dens;
}

CYGNODetectorProperty* CYGNODetectorProperty::GetInstance()
{
  if (fCYGNODetectorProperty == NULL) {
    fCYGNODetectorProperty = new CYGNODetectorProperty();
  }
  return fCYGNODetectorProperty;
}

void CYGNODetectorProperty::Refresh()
{
  vol_name_mass->clear();
  vol_name_dens->clear();
}

void CYGNODetectorProperty::AddVolumeNameMass(G4String name, G4double mass)
{
  vol_name_mass->push_back(std::make_pair(name,mass));
}
void CYGNODetectorProperty::AddVolumeNameDensity(G4String name, G4double dens)
{
  vol_name_dens->push_back(std::make_pair(name,dens));
}

void CYGNODetectorProperty::AddVolumeNameMassAndDensity(G4LogicalVolume* log)
{
  //G4cout << "Saving volume mass and density"<<G4endl;
  //G4cout << "Log Volume address: "<<log<<G4endl;
  G4cout << "Volume name: "<<log->GetName()<<G4endl;
  G4cout << "Volume mass (kg): "<<log->GetMass(true,false)/kg<<G4endl;
  G4String name=log->GetName();
  //if(name.size()>=4 && name.substr(name.size()-4,4)=="_log")
//	name.remove(name.size()-4,4);
  AddVolumeNameMass(name+"_mass",(log->GetMass(true,false))/kg);
  AddVolumeNameDensity(name+"_dens",((log->GetMaterial())->GetDensity())/(g/cm3));
}

void CYGNODetectorProperty::AddPhysVolumeNameMassAndDensity(G4VPhysicalVolume* phys)
{
  //G4cout << "Saving volume mass and density"<<G4endl;
  //G4cout << "Phys Volume address: "<<phys<<G4endl;
  G4cout << "Volume name: "<<phys->GetName()<<G4endl;
  G4cout << "Log Volume address: "<<phys->GetLogicalVolume()<<G4endl;
  G4cout << "Volume mass (kg): "<<(phys->GetLogicalVolume())->GetMass(true,false)/kg<<G4endl;
    
  G4String name=phys->GetName();
  //if(name.size()>=5 && name.substr(name.size()-5,5)=="_phys")
//	name.remove(name.size()-5,5);
  AddVolumeNameMass(name+"_mass",(phys->GetLogicalVolume())->GetMass(true,false)/kg);
  AddVolumeNameDensity(name+"_dens",(((phys->GetLogicalVolume())->GetMaterial())->GetDensity())/(g/cm3));
}
