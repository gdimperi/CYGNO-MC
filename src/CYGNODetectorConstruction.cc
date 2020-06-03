#include <iostream>
#include <string.h>
#include <fstream> 

// GEANT4 //
#include "globals.hh"
#include "G4ThreeVector.hh"
#include "G4Transform3D.hh"
//#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"

#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Polyhedra.hh"
#include "G4Sphere.hh"
#include "G4EllipticalTube.hh"
#include "G4Ellipsoid.hh"
#include "G4IntersectionSolid.hh"
#include "G4UnionSolid.hh"
#include "G4SubtractionSolid.hh"

#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"

#include "G4NistManager.hh"
#include "G4Material.hh"
#include "G4VisAttributes.hh"
#include "G4SDManager.hh"

#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"
#include "G4GeometryManager.hh"
#include "G4RunManager.hh"

#include "G4RotationMatrix.hh"

// USER //
#include "CYGNODetectorConstruction.hh"
#include "CYGNODetectorConstructionMessenger.hh"
#include "CYGNODetectorLNGS.hh"
#include "CYGNODetectorMaterial.hh"
#include "CYGNODetectorProperty.hh"
#include "CYGNOSensitiveDetector.hh"
#include "CYGNOVolumes.hh"

CYGNODetectorConstruction::CYGNODetectorConstruction() :
   CYGNOGeomPath("../geometry/lime/"),
   rockThicknessOuter(-999*m),
   rockThicknessInner(-999*m),
   CYGNOLab("NoCave"),
   //CYGNOLab("LNGS"),
   CYGNOShielding("FullShield"),
   //CYGNOShielding("NoShield"),
   thick0(0.90*m), thick1(0.40*m), thick2(0.20*m), thick3(0.05*m), 
   Mat0("Water"), Mat1("PE"), Mat2("Pb"), Mat3("Cu")
{
	
     fMessenger = new CYGNODetectorConstructionMessenger(this);
}

CYGNODetectorConstruction::~CYGNODetectorConstruction()
{
	delete fMessenger;
}

G4VPhysicalVolume* CYGNODetectorConstruction::Construct()
{
   
    //UpdateGeometryPath(CYGNOGeomPath); 
    //register the SD
    G4SDManager* SDmanager=G4SDManager::GetSDMpointer();
    
    G4String CYGNOSDname = "CYGNO/CYGNOSD";
    CYGNOSD = new CYGNOSensitiveDetector( CYGNOSDname );
    SDmanager->AddNewDetector( CYGNOSD );
	
	
    G4NistManager * nist_manager = G4NistManager::Instance();
   
    //-----------------------------
    // construction of materials
    //-----------------------------
    G4cout << "Constructing materials..." << G4endl;

    //G4Material * air = nist_manager->FindOrBuildMaterial("G4_AIR");
    //G4Material * water = nist_manager->FindOrBuildMaterial("G4_WATER");
    CYGNODetectorMaterial* CYGNOMaterials = CYGNODetectorMaterial::GetInstance();
    G4cout << "... done" << G4endl;

    ////-----------------------------
    //// construction of general properties
    ////-----------------------------
    //G4cout << "Constructing general properties...";
    //CYGNODetectorProperty* CYGNOProperties = CYGNODetectorProperty::GetInstance();
    //G4cout << "... done" << G4endl;




    //**********************************************************************
    //   DEFINITION OF THE GEOMETRY
    //**********************************************************************
      
    //INITIALIZING TRANSLATION VECTORS TO 0:
    Rock_log = 0;
    size_Rock = G4ThreeVector();
    tr_Rock = G4ThreeVector();
    rot_Rock = G4RotationMatrix();
    absrot_Rock = G4RotationMatrix();
    Laboratory_log = 0;
    size_Laboratory = G4ThreeVector();
    tr_Laboratory = G4ThreeVector();
    rot_Laboratory = G4RotationMatrix();
    Shielding_log = 0;
    size_Shielding = G4ThreeVector();
    tr_Shielding = G4ThreeVector();
    rot_Shielding = G4RotationMatrix();
    absrot_Shielding = G4RotationMatrix();
    InsideVolume_log = 0;
    size_InsideVolume = G4ThreeVector();
    tr_InsideVolume = G4ThreeVector();
    rot_InsideVolume = G4RotationMatrix();

    //Name of the volumes
    G4String name_solid="";
    G4String name_log="";
    G4String name_phys="";

    G4ThreeVector tr;
    G4RotationMatrix rot;

    //**********************************************************************
    // WORLD ***************************************
    //**********************************************************************
    G4double world_x = 100.0*m;
    G4double world_y = 100.0*m;
    G4double world_z = 400.0*m;
      
    name_phys="WorldVolume";
    name_log=name_phys+"_log";
    name_solid=name_phys+"_solid";
    G4Box* WorldVolume_box = new G4Box(name_solid,0.5*world_x,0.5*world_y,0.5*world_z);
    WorldVolume_log = new G4LogicalVolume(WorldVolume_box,CYGNOMaterials->Material("Vacuum"),name_log,0,0,0);
    //THE WORLD CANNOT BE TRANSLATED, therefore the first volume inside it (Rock_log) is the volume that must be translated in order to adjust the origin. Rock_log cannot be as large as the world to avoid that this volume is moved out of the world
    WorldVolume_phys = new G4PVPlacement(0,G4ThreeVector(),WorldVolume_log,name_phys,0,false,0,true);//The world volume cannot be translated
          
    //These variables are used to set the thin tube for depth studies in the externalRock_log
    G4double rockdist_z;
    G4double rockdepth_z;
        
    //**********************************************************************
    // LABORATORY ***************************************
    //**********************************************************************
    G4cout << "Constructing laboratory..." << G4endl;
    G4bool isThinTubeCompatible=false;
    // ---------------------------------- LNGS
    if (CYGNOLab == "LNGS"){
    
        /////////////////////////////////////////////////////////
	// In this configuration the translations are wrong, WIP
	// /////////////////////////////////////////////////////
	G4cout << "========== WARNING =============" << G4endl;
	G4cout << "Translation vectors in this configuration ar wrong, WIP, use 'NoCave'" << G4endl;
	G4cout << "================================" << G4endl;
	CYGNODetectorLNGS* LNGS = CYGNODetectorLNGS::GetInstance();
        if (rockThicknessOuter != -999*m)
              LNGS->SetExternalRockThickness(rockThicknessOuter);
        if (productionLayerThickness != -999*m)
              LNGS->SetProductionRockThickness(productionLayerThickness);
        if (rockThicknessInner != -999*m)
              LNGS->SetInternalRockThickness(rockThicknessInner);	  
        
	LNGS->SetDetectorMaterial(CYGNOMaterials);
	LNGS->ConstructRock();
        Rock_log=LNGS->GetRock();
        size_Rock=LNGS->GetRockSizeXYZ();
        absrot_Rock=LNGS->GetRockAbsRotation();
        Laboratory_log=LNGS->GetLaboratory();
        size_Laboratory=LNGS->GetLaboratorySizeXYZ();
        tr_Laboratory=LNGS->GetLaboratoryTranslation();
        rot_Laboratory=LNGS->GetLaboratoryRotation();
        
        ////for thin tube
        //if(LNGS->GetProductionRockThickness()==0.*cm && LNGS->GetInternalRockThickness()==0.*cm)
        //  isThinTubeCompatible=true;
        //rockdist_z=LNGS->rockdist_z;
        //rockdepth_z=LNGS->rockdepth_z;
    }
    // ---------------------------------- NoCave
    else if (CYGNOLab == "NoCave") 
	{
	  //**********************************************************************
	  // Double Air sphere surrounding the whole detector and shielding
	  //**********************************************************************        
	  G4double airInnerRadius = 10.0*m;
	  G4double airThickness = 5.0*m;
	  if (productionLayerThickness != -999*m)
		airThickness = productionLayerThickness;
	  if (rockThicknessInner != -999*m)
		airInnerRadius=rockThicknessInner;

	  //Air permeates the are around the detector
	  name_phys="OuterAirSphere";
	  name_log=name_phys+"_log";
	  name_solid=name_phys+"_solid";
	  G4Sphere* OuterAirSphere = new G4Sphere(name_solid, 0., airInnerRadius+airThickness, 0*degree, 360*degree, 0*degree,180*degree);
	  G4LogicalVolume* OuterAirSphere_log = new G4LogicalVolume(OuterAirSphere,CYGNOMaterials->Material("Air"),name_log);
	  Rock_log=OuterAirSphere_log;
	  size_Rock=G4ThreeVector(airInnerRadius+airThickness,airInnerRadius+airThickness,airInnerRadius+airThickness);
	  absrot_Rock = G4RotationMatrix();
		
	  name_phys="InnerAirSphere";
	  name_log=name_phys+"_log";
	  name_solid=name_phys+"_solid";
	  G4Sphere* InnerAirSphere = new G4Sphere(name_solid,  0., airInnerRadius,  0*degree, 360*degree, 0*degree,180*degree);
	  G4LogicalVolume* InnerAirSphere_log = new G4LogicalVolume(InnerAirSphere, CYGNOMaterials->Material("Air"), name_log);
	  Laboratory_log=InnerAirSphere_log;
	  size_Laboratory=G4ThreeVector(airInnerRadius,airInnerRadius,airInnerRadius);
	  tr = G4ThreeVector(0.,0.,0.);//translation in mother frame
	  tr_Laboratory+=(rot_Laboratory*tr);
	  rot = G4RotationMatrix();// rotation of daughter volume
	  rot_Laboratory*=rot; //equivalent to rot_Laboratory=rot_Laboratory*rot
	  InnerAirSphere_phys = new G4PVPlacement(G4Transform3D(rot,tr),InnerAirSphere_log,name_phys,OuterAirSphere_log,false,0,true);        
    }
    else
    {
      G4cout << "ERROR: Something went wrong with the definition of the variable CYGNOLab" << G4endl;
      //throw std::exception();
      exit(1);		
    }


    G4cout << "Laboratory done." << G4endl;
 
    //**********************************************************************
    // SHIELDING **********************************************************
          
    if (CYGNOShielding == "FullShield") 
    {
        G4cout<<"Shielding Construction Started"<<G4endl;
        
        if (thick3==-999*m || thick2==-999*m || thick1==-999*m || thick0==-999*m || Mat0=="" || Mat1=="" || Mat2=="" || Mat3=="") 
        {
            G4cout << "ERROR: Please specify the thickennses and the materials of the shielding layers if you want to use the shielding design: FullShield" << G4endl;
            //throw std::exception();
            exit(1);
        }

        // ----------------------------------- Inner room dimensions
        G4double AirBox_x;
        G4double AirBox_y;
        G4double AirBox_z;
        G4Box* AirBox;
        AirBox_x = 2.0*m; //2.65*m;
        AirBox_y = 0.8*m; //1.45*m;
        AirBox_z = 0.8*m; //1.45*m;        
        tr_InsideVolume = G4ThreeVector(0.,0.,0.);
        rot_InsideVolume = G4RotationMatrix();		
        size_InsideVolume = G4ThreeVector(AirBox_x/2.,
              								AirBox_y/2.,
              								AirBox_z/2.);		
        size_Shielding = G4ThreeVector(AirBox_x/2. + thick3 + thick2 + thick1 + thick0,
              							 AirBox_y/2. + thick3 + thick2 + thick1 + thick0,
              							 AirBox_z/2. + thick3 + thick2 + thick1 + thick0);
        absrot_Shielding = G4RotationMatrix();

        // ----------------------------------- Shield 0
        G4double Shield0_x = AirBox_x + 2.*thick3 + 2.*thick2 + 2.*thick1 + 2.*thick0 ;
        G4double Shield0_y = AirBox_y + 2.*thick3 + 2.*thick2 + 2.*thick1 + 2.*thick0 ;
        G4double Shield0_z = AirBox_z + 2.*thick3 + 2.*thick2 + 2.*thick1 + 2.*thick0 ;        
        G4Material* Shield0Mat = CYGNOMaterials->Material(Mat0);
        name_phys="Shield0";
        name_log=name_phys+"_log";
        name_solid=name_phys+"_solid";
        G4Box* Shield0 = new G4Box(name_solid,0.5*Shield0_x,0.5*Shield0_y,0.5*Shield0_z);
        Shield0_log = new G4LogicalVolume(Shield0,Shield0Mat,name_log);
        Shielding_log = Shield0_log;
        Shield0_log->SetVisAttributes(CYGNOMaterials->VisAttributes(Mat0));

        // ----------------------------------- Shield 1
        G4double Shield1_x = AirBox_x + 2.*thick3 + 2.*thick2 + 2.*thick1 ;
        G4double Shield1_y = AirBox_y + 2.*thick3 + 2.*thick2 + 2.*thick1 ;
        G4double Shield1_z = AirBox_z + 2.*thick3 + 2.*thick2 + 2.*thick1 ;
        G4Material* Shield1Mat = CYGNOMaterials->Material(Mat1);
        name_phys="Shield1";
        name_log=name_phys+"_log";
        name_solid=name_phys+"_solid";
        G4Box* Shield1 = new G4Box(name_solid,0.5*Shield1_x,0.5*Shield1_y,0.5*Shield1_z);
        Shield1_log = new G4LogicalVolume(Shield1,Shield1Mat,name_log);
        Shield1_log->SetVisAttributes(CYGNOMaterials->VisAttributes(Mat1));
        
        // ----------------------------------- Shield 2        
        G4double Shield2_x = AirBox_x + 2.*thick3 + 2.*thick2 ;
        G4double Shield2_y = AirBox_y + 2.*thick3 + 2.*thick2 ;
        G4double Shield2_z = AirBox_z + 2.*thick3 + 2.*thick2 ;
        G4Material* Shield2Mat = CYGNOMaterials->Material(Mat2);
        name_phys="Shield2";
        name_log=name_phys+"_log";
        name_solid=name_phys+"_solid";
        G4Box* Shield2 = new G4Box(name_solid,0.5*Shield2_x,0.5*Shield2_y,0.5*Shield2_z);
        Shield2_log = new G4LogicalVolume(Shield2,Shield2Mat,name_log);
        Shield2_log->SetVisAttributes(CYGNOMaterials->VisAttributes(Mat2));
        
        // ----------------------------------- Shield 3        
        G4double Shield3_x = AirBox_x + 2.*thick3 ;
        G4double Shield3_y = AirBox_y + 2.*thick3 ;
        G4double Shield3_z = AirBox_z + 2.*thick3 ;
        G4Material* Shield3Mat = CYGNOMaterials->Material(Mat3);
        name_phys="Shield3";
        name_log=name_phys+"_log";
        name_solid=name_phys+"_solid";
        G4Box* Shield3 = new G4Box(name_solid,0.5*Shield3_x,0.5*Shield3_y,0.5*Shield3_z);
        Shield3_log = new G4LogicalVolume(Shield3,Shield3Mat,name_log);
        Shield3_log->SetVisAttributes(CYGNOMaterials->VisAttributes(Mat3));
        
        // ----------------------------------- Airbox
        name_phys="AirBox";
        name_log=name_phys+"_log";
        name_solid=name_phys+"_solid";
        AirBox = new G4Box(name_solid,0.5*AirBox_x,0.5*AirBox_y,0.5*AirBox_z);
        AirBox_log = new G4LogicalVolume(AirBox,CYGNOMaterials->Material("Air"),name_log);
        AirBox_log->SetVisAttributes(CYGNOMaterials->VisAttributes("Air"));
	InsideVolume_log = AirBox_log;
    
    }
    // ---------------------------------- No shielding
    else if (CYGNOShielding == "NoShield") 
    {
	  G4double AirBox_x;
	  G4double AirBox_y;
	  G4double AirBox_z;
	  G4Box* AirBox;
          AirBox_x = 5.*m;
          AirBox_y = 3.*m;
          AirBox_z = 3.*m;        
	  name_phys="AirBox";
	  name_log=name_phys+"_log";
	  name_solid=name_phys+"_solid";
	  AirBox = new G4Box(name_solid,0.5*AirBox_x,0.5*AirBox_y,0.5*AirBox_z);
	  AirBox_log = new G4LogicalVolume(AirBox,CYGNOMaterials->Material("Air"),name_log,0,0,0);

	  Shielding_log=AirBox_log;
	  size_Shielding=G4ThreeVector(AirBox_x/2.,AirBox_y/2.,AirBox_z/2.);
	  absrot_Shielding = G4RotationMatrix();
	  size_InsideVolume=G4ThreeVector(AirBox_x/2.,AirBox_y/2.,AirBox_z/2.);
	  tr_InsideVolume=G4ThreeVector(0.,0.,0.);
	  rot_InsideVolume=G4RotationMatrix();
    }
    else
    {
      G4cout << "ERROR: Something went wrong with the definition of the variable CYGNOShielding" << G4endl;
      //throw std::exception();
      exit(1);		
    }
     
    /// Cameras
    G4double x_camera_body = 125*mm;
    G4double y_camera_body = 85*mm;
    G4double z_camera_body = 85*mm;
    G4Box* camera_body = new G4Box("camera_body_solid",0.5*x_camera_body,0.5*y_camera_body,0.5*z_camera_body);

    G4double z_camera_lens = 44.9*mm;
    G4double diam_camera_lens = 48*mm;
    G4Tubs* camera_lens = new G4Tubs("camera_lens_solid",0.,0.5*diam_camera_lens,0.5*z_camera_lens,0.*deg,360.*deg);
    
  //keep camera and lens separated
    camera_log = new G4LogicalVolume(camera_body,CYGNOMaterials->Material("Camera"),"camera_log",0,0,0) ;
    camera_lens_log = new G4LogicalVolume(camera_lens,CYGNOMaterials->Material("Camera"),"camera_lens_log",0,0,0) ;

    G4double ztr_cam = 785.*mm + z_camera_lens+0.5*x_camera_body + 1*mm;
    G4ThreeVector trcam0(ztr_cam,-15.*mm,0.);
    G4RotationMatrix* rotcam0 = new G4RotationMatrix();
    
    G4ThreeVector trlens0(ztr_cam-0.5*x_camera_body-0.5*z_camera_lens,-15.*mm,0.);
    G4RotationMatrix* rotlens0 = new G4RotationMatrix();
    rotlens0->rotateY(90*deg);

    camera_phys = new G4PVPlacement(rotcam0,trcam0,camera_log,"camera",AirBox_log, false, 0, true);
    camera_lens_phys = new G4PVPlacement(rotlens0,trlens0,camera_lens_log,"camera_lens",AirBox_log, false, 0, true);
    G4double tolerance = 1*mm;

    //**********************************************************************
    // ********* CYGNO volumes form CADMesh *****************************
    //**********************************************************************
    
    char namestl[50];
    sprintf(namestl,"%s/LIMEDetectorBody.stl",CYGNOGeomPath.c_str());
    G4cout << namestl << G4endl;
    //CADMesh * mesh_LIMEDetectorBody = new CADMesh("../geometry/v2/LIMEDetectorBody.stl");    
    ifstream infile(CYGNOGeomPath.c_str());
    if (infile.good())
    	mesh_LIMEDetectorBody = new CADMesh(namestl);    
    //namestl = "LIMEinternalStructure.stl";
    sprintf(namestl,"%s/LIMEinternalStructure.stl",CYGNOGeomPath.c_str());
    G4cout << namestl << G4endl;
    if (infile.good())
    	mesh_LIMEinternalStructure = new CADMesh(namestl);
    //namestl = "cameras.stl";    
//    sprintf(namestl,"%s/cameras.stl",CYGNOGeomPath.c_str());
//    G4cout << namestl << G4endl;
//    if (infile.good())
//      mesh_camera = new CADMesh(namestl);    
    //CADMesh * mesh_window = new CADMesh("../geometry/v1/glass_windows.stl");  
    //namestl = "LIMEendPMT.stl";  
    sprintf(namestl,"%s/LIMEendPMT.stl",CYGNOGeomPath.c_str());
    G4cout << namestl << G4endl;
    if (infile.good())
      mesh_LIMEendPMT = new CADMesh(namestl);
    //namestl = "turns_support.stl";    
    //sprintf(namestl,"%s/turns_support.stl",CYGNOGeomPath.c_str());
    //G4cout << namestl << G4endl;
    //if (infile.good())
    // mesh_turns_support = new CADMesh(namestl);
    //namestl = "FieldRings.stl";    
    sprintf(namestl,"%s/FieldRings.stl",CYGNOGeomPath.c_str());
    G4cout << namestl << G4endl;
    if (infile.good())
      mesh_FieldRings = new CADMesh(namestl);
    //namestl = "gem_frame.stl";    
    //sprintf(namestl,"%s/gem_frame.stl",CYGNOGeomPath.c_str());
    //G4cout << namestl << G4endl;
    sprintf(namestl,"%s/GEMstretchers.stl",CYGNOGeomPath.c_str());
    G4cout << namestl << G4endl;
    if (infile.good())
      mesh_GEMstretchers= new CADMesh(namestl);  
    sprintf(namestl,"%s/GEMsupportStructure.stl",CYGNOGeomPath.c_str());
    G4cout << namestl << G4endl;     
    if (infile.good())
      mesh_GEMsupportStructure = new CADMesh(namestl); 
    //namestl = "GEMsupportStructure.stl";

    //namestl = "GEMfoils.stl";
    sprintf(namestl,"%s/GEMfoils.stl",CYGNOGeomPath.c_str());
    G4cout << namestl << G4endl;
    if (infile.good())
      mesh_GEMfoils = new CADMesh(namestl);

//    //namstl = "SupportBenchLime.stl";    
//    sprintf(namestl,"%s/SupportBenchLime.stl",CYGNOGeomPath.c_str());
//    G4cout << namestl << G4endl;
//    if (infile.good())
//      mesh_SupportBenchLime = new CADMesh(namestl);

    //namestl = "Cathode.stl";   
    sprintf(namestl,"%s/Cathode.stl",CYGNOGeomPath.c_str());
    G4cout << namestl << G4endl;
    if (infile.good())
      mesh_Cathode = new CADMesh(namestl);   

    if (infile.good()){
      mesh_LIMEDetectorBody->SetScale(mm);
      mesh_LIMEinternalStructure->SetScale(mm);
      //mesh_camera->SetScale(mm);
      //mesh_window->SetScale(mm);
      mesh_LIMEendPMT->SetScale(mm);
      mesh_FieldRings->SetScale(mm);
      mesh_GEMstretchers->SetScale(mm);
      mesh_GEMsupportStructure->SetScale(mm);
      mesh_GEMfoils->SetScale(mm);
//      mesh_SupportBenchLime->SetScale(mm);
      mesh_Cathode->SetScale(mm);
    

      //LIME Detector Body
      cad_LIMEDetectorBody_solid = mesh_LIMEDetectorBody->TessellatedMesh();
      cad_LIMEDetectorBody_logical = new G4LogicalVolume(cad_LIMEDetectorBody_solid, CYGNOMaterials->Material("Perspex"), "cad_LIMEDetectorBody_logical");
      cad_LIMEDetectorBody_logical->SetVisAttributes(CYGNOMaterials->VisAttributes("Perspex"));


      //LIME Internal Structure
      cad_LIMEinternalStructure_solid = mesh_LIMEinternalStructure->TessellatedMesh();
      cad_LIMEinternalStructure_logical = new G4LogicalVolume(cad_LIMEinternalStructure_solid, CYGNOMaterials->Material("Perspex"), "cad_LIMEinternalStructure_logical", 0, 0, 0);
      //cad_LIMEinternalStructure_logical->SetVisAttributes(CYGNOMaterials->VisAttributes(cad_LIMEinternalStructure_logical->GetMaterial()->GetName()));
      cad_LIMEinternalStructure_logical->SetVisAttributes(CYGNOMaterials->VisAttributes("Perspex"));
    }
  
    //TPC gas
    G4double TPC_x = 640.*mm;//640.*mm;
    G4double TPC_y = 500.*mm; //500
    G4double TPC_z = 460.*mm; //470.*mm;
      
    name_phys="TPC";
    name_log=name_phys+"_log";
    name_solid=name_phys+"_solid";
    G4Box* TPC_box = new G4Box(name_solid,0.5*TPC_x,0.5*TPC_y,0.5*TPC_z);
    TPC_log = new G4LogicalVolume(TPC_box,CYGNOMaterials->Material("CYGNO_gas"),name_log,0,0,0);
    
    //CYGNO fiducial gas
    G4double CYGNO_x = 500.*mm;
    G4double CYGNO_y = 300.*mm;
    G4double CYGNO_z = 300.*mm;//300.*mm;
      
    name_phys="CYGNO";
    name_log=name_phys+"_log";
    name_solid=name_phys+"_solid";
    G4Box* CYGNO_box = new G4Box(name_solid,0.5*CYGNO_x,0.5*CYGNO_y,0.5*CYGNO_z);
    CYGNO_log = new G4LogicalVolume(CYGNO_box,CYGNOMaterials->Material("CYGNO_gas"),name_log,0,0,0);
    //CYGNO_log->SetVisAttributes(CYGNOMaterials->VisAttributes(CYGNO_log->GetMaterial()->GetName()));
    CYGNO_log->SetVisAttributes(CYGNOMaterials->VisAttributes("CYGNO_gas"));

    if (infile.good()){
      //LIMEendPMT
      cad_LIMEendPMT_solid = mesh_LIMEendPMT->TessellatedMesh();
      cad_LIMEendPMT_logical = new G4LogicalVolume(cad_LIMEendPMT_solid, CYGNOMaterials->Material("Perspex"), "cad_LIMEendPMT_logical", 0, 0, 0);
      //cad_LIMEendPMT_logical->SetVisAttributes(CYGNOMaterials->VisAttributes(cad_LIMEendPMT_logical->GetMaterial()->GetName()));
      cad_LIMEendPMT_logical->SetVisAttributes(CYGNOMaterials->VisAttributes("Perspex"));
 
      //turns support
//      cad_turns_support_solid = mesh_turns_support->TessellatedMesh();
//      cad_turns_support_logical = new G4LogicalVolume(cad_turns_support_solid, CYGNOMaterials->Material("Perspex"), "cad_turns_support_logical", 0, 0, 0);
      //cad_turns_support_logical->SetVisAttributes(CYGNOMaterials->VisAttributes(cad_turns_support_logical->GetMaterial()->GetName()));
//      cad_turns_support_logical->SetVisAttributes(CYGNOMaterials->VisAttributes("Perspex"));

      //field cage  
      cad_FieldRings_solid = mesh_FieldRings->TessellatedMesh();
      cad_FieldRings_logical = new G4LogicalVolume(cad_FieldRings_solid, CYGNOMaterials->Material("Cu"), "cad_FieldRings_logical", 0, 0, 0);
      //cad_FieldRings_logical->SetVisAttributes(CYGNOMaterials->VisAttributes(cad_FieldRings_logical->GetMaterial()->GetName()));
      cad_FieldRings_logical->SetVisAttributes(CYGNOMaterials->VisAttributes("Cu"));

      //GEM stretchers
      cad_GEMstretchers_solid = mesh_GEMstretchers->TessellatedMesh();
      cad_GEMstretchers_logical = new G4LogicalVolume(cad_GEMstretchers_solid, CYGNOMaterials->Material("Perspex"), "cad_GEMstretchers_logical", 0, 0, 0);
      //cad_GEMstretchers_logical->SetVisAttributes(CYGNOMaterials->VisAttributes(cad_GEMstretchers_logical->GetMaterial()->GetName()));
      cad_GEMstretchers_logical->SetVisAttributes(CYGNOMaterials->VisAttributes("Perspex"));
      
      //GEM structure support
      cad_GEMsupportStructure_solid = mesh_GEMsupportStructure->TessellatedMesh();
      cad_GEMsupportStructure_logical = new G4LogicalVolume(cad_GEMsupportStructure_solid, CYGNOMaterials->Material("Perspex"), "cad_GEMsupportStructure_logical", 0, 0, 0);
      cad_GEMsupportStructure_logical->SetVisAttributes(CYGNOMaterials->VisAttributes("Perspex"));


      //GEMfoils
      cad_GEMfoils_solid = mesh_GEMfoils->TessellatedMesh();
      cad_GEMfoils_logical = new G4LogicalVolume(cad_GEMfoils_solid, CYGNOMaterials->Material("GEM"), "cad_GEMfoils_logical", 0, 0, 0); //GEM material is an effective material of kapton + copper 
      //cad_GEMfoils_logical->SetVisAttributes(CYGNOMaterials->VisAttributes(cad_GEMfoils_logical->GetMaterial()->GetName()));
      
//      //Support Bench LIME
//      cad_SupportBenchLime_solid = mesh_SupportBenchLime->TessellatedMesh();
//      cad_SupportBenchLime_logical = new G4LogicalVolume(cad_SupportBenchLime_solid, CYGNOMaterials->Material("Cu"), "cad_SupportBenchLime_logical", 0, 0, 0);
//      //cad_SupportBenchLime_logical->SetVisAttributes(CYGNOMaterials->VisAttributes(cad_SupportBenchLime_logical->GetMaterial()->GetName()));
//      cad_SupportBenchLime_logical->SetVisAttributes(CYGNOMaterials->VisAttributes("Cu"));
      
      //cathode
      cad_Cathode_solid = mesh_Cathode->TessellatedMesh();
      cad_Cathode_logical = new G4LogicalVolume(cad_Cathode_solid, CYGNOMaterials->Material("Cu"), "cad_Cathode_logical", 0, 0, 0);
      cad_Cathode_logical->SetVisAttributes(CYGNOMaterials->VisAttributes("Cu"));
    }


    if (CYGNOLab == "LNGS"){
	tr+=G4ThreeVector(0.,-1*size_Laboratory.y()+size_Shielding.y(),size_Laboratory.z()-10*m);
	tr_cad+=G4ThreeVector(0.,1.0*m-1*size_Laboratory.y()+size_Shielding.y(),size_Laboratory.z()-10*m);	  
    
	rot = G4RotationMatrix();// rotation of daughter volume
	tr_Shielding+=(rot_Shielding*tr);
    }
    else if (CYGNOLab == "NoCave") {
	tr=G4ThreeVector(0.,0.,0.);
	tr_cad+=G4ThreeVector(0.,0.,0.);
	rot = G4RotationMatrix();
	tr_Shielding+=(rot_Shielding*tr);

    }

    Shield0_phys = new G4PVPlacement(G4Transform3D(rot,tr),Shielding_log,"Shield0",Laboratory_log,false,0,true);
    if (CYGNOShielding == "NoShield")
    {
        tr = G4ThreeVector(0.,0.,0.);//translation in mother frame
        tr_InsideVolume+=(rot_InsideVolume*tr);
        rot = G4RotationMatrix();// rotation of daughter volume
        rot_InsideVolume*=rot; //equivalent to rot_InsideVolume=rot_InsideVolume*rot
        AirBox_log = Shielding_log;
    }
    
    else if (CYGNOShielding == "FullShield") 
    {   
        // ----------------------------------- Volume placements

        tr = G4ThreeVector(0.,0.,0.);//translation in mother frame
        tr_InsideVolume+=(rot_InsideVolume*tr);
        rot = G4RotationMatrix();// rotation of daughter volume
        rot_InsideVolume*=rot; //equivalent to rot_InsideVolume=rot_InsideVolume*rot
        Shield1_phys = new G4PVPlacement(G4Transform3D(rot,tr),Shield1_log,"Shield1",Shielding_log,false,0,true);

        tr = G4ThreeVector(0.,0.,0.);//translation in mother frame
        tr_InsideVolume+=(rot_InsideVolume*tr);
        rot = G4RotationMatrix();// rotation of daughter volume
        rot_InsideVolume*=rot; //equivalent to rot_InsideVolume=rot_InsideVolume*rot
        Shield2_phys = new G4PVPlacement(G4Transform3D(rot,tr),Shield2_log,"Shield2",Shield1_log,false,0,true);

        tr = G4ThreeVector(0.,0.,0.);//translation in mother frame
        tr_InsideVolume+=(rot_InsideVolume*tr);
        rot = G4RotationMatrix();// rotation of daughter volume
        rot_InsideVolume*=rot; //equivalent to rot_InsideVolume=rot_InsideVolume*rot
        Shield3_phys = new G4PVPlacement(G4Transform3D(rot,tr),Shield3_log,"Shield3",Shield2_log,false,0,true);

        tr = G4ThreeVector(0.,0.,0.);//translation in mother frame
        tr_InsideVolume+=(rot_InsideVolume*tr);
        rot = G4RotationMatrix();// rotation of daughter volume
        rot_InsideVolume*=rot; //equivalent to rot_InsideVolume=rot_InsideVolume*rot
        AirBox_phys = new G4PVPlacement(G4Transform3D(rot,tr), AirBox_log, "AirBox", Shield3_log, false, 0,true); 
    
    }
    G4ThreeVector  size;

    tr_cad=G4ThreeVector(-3847*mm+295*mm,-3847*mm,230.*mm);
    if (infile.good()){
      cad_LIMEDetectorBody_physical = new G4PVPlacement(G4Transform3D(rot,tr_cad), 
  		    cad_LIMEDetectorBody_logical,"cad_LIMEDetectorBody_physical", AirBox_log, false, 0, true);
      cad_LIMEendPMT_physical = new G4PVPlacement(G4Transform3D(rot,tr_cad), 
        	    cad_LIMEendPMT_logical,"cad_LIMEendPMT_physical", AirBox_log, false, 0, true);
//      cad_SupportBenchLime_physical = new G4PVPlacement(G4Transform3D(rot,tr), 
//  		    cad_SupportBenchLime_logical,"cad_SupportBenchLime_physical", AirBox_log, false, 0, true);
    }
    tr=G4ThreeVector(0.,0.,0.);
    TPC_phys = new G4PVPlacement(G4Transform3D(rot,tr),
      	    TPC_log,"TPC_gas", AirBox_log, false, 0, true);
    tr_CYGNO_gas_1=G4ThreeVector(TPC_x/2.-CYGNO_x/2.-50.*mm,-20*mm,0.);
    CYGNO_phys = new G4PVPlacement(G4Transform3D(rot,tr_CYGNO_gas_1),
      	    CYGNO_log,"CYGNO_gas", TPC_log, false, 0, true);
          
    tr=G4ThreeVector(0.,0.,0.);
    rot = G4RotationMatrix();
    if (infile.good()){
      cad_LIMEinternalStructure_physical = new G4PVPlacement(G4Transform3D(rot,tr_cad), 
        	    cad_LIMEinternalStructure_logical,"cad_LIMEinternalStructure_physical", TPC_log, false, 0, true);
      cad_FieldRings_physical = new G4PVPlacement(G4Transform3D(rot,tr_cad), 
        	    cad_FieldRings_logical,"cad_FieldRings_physical", TPC_log, false, 0, true);
      cad_GEMstretchers_physical = new G4PVPlacement(G4Transform3D(rot,tr_cad), 
        	    cad_GEMstretchers_logical,"cad_GEMstretchers_physical", TPC_log, false, 0, true);
      cad_GEMsupportStructure_physical = new G4PVPlacement(G4Transform3D(rot,tr_cad), 
        	    cad_GEMsupportStructure_logical,"cad_GEMsupportStructure_physical", TPC_log, false, 0, true);
      cad_GEMfoils_physical = new G4PVPlacement(G4Transform3D(rot,tr_cad), 
  		    cad_GEMfoils_logical,"cad_GEMfoils_physical", TPC_log, false, 0, true);
      cad_Cathode_physical = new G4PVPlacement(G4Transform3D(rot,tr_cad), 
  		    cad_Cathode_logical,"cad_Cathode_physical", TPC_log, false, 0, true);
    }  
   
    //
    //**********************************************************************
    // GLOBAL TRANSLATIONS ***************************************
    G4cout<<"Placement of Laboratory in the World started"<<G4endl;
    
    tr_Rock=-1*(tr_Laboratory+rot_Laboratory*(tr+rot_cad*(tr_InsideVolume+rot_InsideVolume*(tr_Shielding))));//The shift of Rock_log in the world volume to make the origin be the center of the detector

    G4RotationMatrix rot_check = absrot_Rock*(rot_Laboratory*(rot_cad*(rot_InsideVolume)));  
    name_phys="externalRock";
    name_log=name_phys+"_log";
    name_solid=name_phys+"_solid";
    externalRock_phys = new G4PVPlacement(G4Transform3D(absrot_Rock,tr_Rock),Rock_log,name_phys,WorldVolume_log,false,0,true);

    G4cout << "The main CYGNO volume is translated w.r.t the center of the rock volume of:\t x="<< -1*tr_Rock.x()/cm << " cm\t y=" << -1*tr_Rock.y()/cm << " cm\t z=" << -1*tr_Rock.z()/cm << " cm"<< G4endl;

    G4cout << "The Rock volume has been translated to put the main CYGNO volume in the center of the coordinate system"<< G4endl;
    G4cout<<"Placement of Laboratory in the World ended"<<G4endl;
    
    //======= Sensitive detector ========
    CYGNO_log->SetSensitiveDetector(CYGNOSD);
 
    //======= Save volumes mass and density ======
    
    G4cout<<"Saving masses and densities of the volumes"<<G4endl;
    SaveMassAndDensity();


    //===========
    return WorldVolume_phys;


}

void CYGNODetectorConstruction::SaveMassAndDensity()
{
  CYGNODetectorProperty* CYGNOProperties = CYGNODetectorProperty::GetInstance();

  G4cout << "Saving masses and densities of the volumes of the CYGNODetectorConstruction class"<< G4endl;
  //CYGNOProperties->AddVolumeNameMassAndDensity(Rock_log);
  //CYGNOProperties->AddVolumeNameMassAndDensity(Laboratory_log);
  if (CYGNOShielding=="FullShield"){
      CYGNOProperties->AddVolumeNameMassAndDensity(Shielding_log);
      CYGNOProperties->AddVolumeNameMassAndDensity(Shield1_log);
      CYGNOProperties->AddVolumeNameMassAndDensity(Shield2_log);
      CYGNOProperties->AddVolumeNameMassAndDensity(Shield3_log);
  }
  CYGNOProperties->AddVolumeNameMassAndDensity(AirBox_log);
  ifstream infile(CYGNOGeomPath.c_str());
  if (infile.good()) {
    CYGNOProperties->AddVolumeNameMassAndDensity(cad_LIMEDetectorBody_logical);
    CYGNOProperties->AddVolumeNameMassAndDensity(cad_LIMEinternalStructure_logical);
    CYGNOProperties->AddVolumeNameMassAndDensity(camera_log);
    CYGNOProperties->AddVolumeNameMassAndDensity(camera_lens_log);
    CYGNOProperties->AddVolumeNameMassAndDensity(TPC_log);
    CYGNOProperties->AddVolumeNameMassAndDensity(CYGNO_log);
    CYGNOProperties->AddVolumeNameMassAndDensity(cad_LIMEendPMT_logical);
    CYGNOProperties->AddVolumeNameMassAndDensity(cad_GEMstretchers_logical);
    CYGNOProperties->AddVolumeNameMassAndDensity(cad_GEMsupportStructure_logical);
    CYGNOProperties->AddVolumeNameMassAndDensity(cad_GEMfoils_logical);
    CYGNOProperties->AddVolumeNameMassAndDensity(cad_Cathode_logical);
    CYGNOProperties->AddVolumeNameMassAndDensity(cad_FieldRings_logical);
  }

//  if ( productionRockThinTube_phys )
//	{
//	  CYGNOProperties->AddPhysVolumeNameMassAndDensity(productionRockThinTube_phys);
//	}
  G4cout << "All volume masses and densities saved"<< G4endl;
}

void CYGNODetectorConstruction::UpdateGeometry()
{

  G4cout << "Updating the Geometry"<< G4endl;
  CYGNODetectorProperty* CYGNOProperties = CYGNODetectorProperty::GetInstance();
  CYGNOProperties->Refresh();

  //Removing sensitive detectors
  CYGNO_log->SetSensitiveDetector(0);
  CYGNO_log=0;

  //Deleting all the solids, logical and physical objects
  //G4RunManager::GetRunManager()->ReinitializeGeometry(true);
  //Equivalent to
  G4GeometryManager::GetInstance()->OpenGeometry();
  G4PhysicalVolumeStore* PhysicalVolumeStore = G4PhysicalVolumeStore::GetInstance();
  PhysicalVolumeStore->Clean();
  G4LogicalVolumeStore* LogicalVolumeStore = G4LogicalVolumeStore::GetInstance();
  LogicalVolumeStore->Clean();
  G4SolidStore* SolidStore = G4SolidStore::GetInstance();
  SolidStore->Clean();
  
  CYGNODetectorMaterial* CYGNOMaterials = CYGNODetectorMaterial::GetInstance();
  CYGNOMaterials->Refresh();


  //The memory for these pointers has been freed by the above Clean() methods
  CYGNODetectorLNGS* CYGNOLNGS = CYGNODetectorLNGS::GetInstance();
  CYGNOLNGS->Refresh();

  Shield0_log=0; 
  Shield1_log=0; 
  Shield2_log=0; 
  Shield3_log=0; 
  AirBox_log=0;

  cad_LIMEDetectorBody_logical=0;
  cad_LIMEinternalStructure_logical=0;
 // cad_cameras_all_logical=0;
  TPC_log=0;
  CYGNO_log=0;
  cad_LIMEendPMT_logical=0;
  cad_GEMsupportStructure_logical=0;
  cad_GEMstretchers_logical=0;
  cad_GEMfoils_logical=0;
//  cad_SupportBenchLime_logical=0;
  cad_Cathode_logical=0;
  cad_FieldRings_logical=0;
 // camera_log=0;
 // camera_lens_log=0;
 // camera_shield_log=0;


  InsideVolume_log=0;
  Shielding_log=0;
  Laboratory_log=0;
  Rock_log=0;
  WorldVolume_log=0;
  productionRockThinTube_phys=0;
  //log->ClearDaughters();

  // Delete all the geometry you had defined and build everything from scratch
  G4RunManager::GetRunManager()->DefineWorldVolume(Construct());
  G4RunManager::GetRunManager()->GeometryHasBeenModified();
}
void CYGNODetectorConstruction::UpdateGeometryPath(G4String newpath)
{
  G4cout << "Updating the Geometry path to "<< newpath << G4endl;
  SetGeomPath(newpath); 

}
