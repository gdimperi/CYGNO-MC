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
   CYGNOGeomPath("../geometry/v2/"),
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
        AirBox_x = 2.55*m;
        AirBox_y = 1.45*m;
        AirBox_z = 1.45*m;        
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
          AirBox_x = 2.55*m;
          AirBox_y = 1.45*m;
          AirBox_z = 1.45*m;        
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
    G4RotationMatrix* rotlens = new G4RotationMatrix();
    rotlens->rotateY(90.*deg);
    G4ThreeVector translens(0.5*z_camera_lens+0.5*x_camera_body,0.,0.);

    //keep camera and lens separated
    //G4UnionSolid* camera = new G4UnionSolid("camera_solid", camera_body, camera_lens,  rotlens, translens);
    camera_log = new G4LogicalVolume(camera_body,CYGNOMaterials->Material("Camera"),"camera_log",0,0,0) ;
    camera_lens_log = new G4LogicalVolume(camera_lens,CYGNOMaterials->Material("Camera"),"camera_lens_log",0,0,0) ;


    ////Camera shielding
    //box
    G4double x_camera_shield_box = 44.9*mm;
    G4double y_camera_shield_box = 1020*mm;
    G4double z_camera_shield_box = 1020*mm;
    G4Box* camera_shield_box = new G4Box("camera_shield_box_solid",0.5*x_camera_shield_box,0.5*y_camera_shield_box,0.5*z_camera_shield_box);
    G4double tolerance = 1*mm;
    //holes
    G4Tubs* camera_shield_hole = new G4Tubs("camera_shield_hole_solid",0.,0.5*diam_camera_lens+tolerance,0.5*z_camera_lens,0.*deg,360.*deg); 
    G4RotationMatrix* rotholes = rotlens;
    G4ThreeVector trhole0(0.,0.,0.);
    G4ThreeVector trhole1(0.,0.,-354*mm);
    G4ThreeVector trhole2(0.,0.,354*mm);
    G4ThreeVector trhole3(0.,-354*mm,0.);
    G4ThreeVector trhole4(0.,-354*mm,-354*mm);
    G4ThreeVector trhole5(0.,-354*mm,354*mm);
    G4ThreeVector trhole6(0.,354*mm,0.);
    G4ThreeVector trhole7(0.,354*mm,-354*mm);
    G4ThreeVector trhole8(0.,354*mm,354*mm);
    //subtraction of holes
    G4SubtractionSolid* camera_shield = new G4SubtractionSolid("camera_shield_solid", camera_shield_box, camera_shield_hole, rotholes, trhole0);
    camera_shield = new G4SubtractionSolid("camera_shield_solid", camera_shield, camera_shield_hole, rotholes, trhole1);
    camera_shield = new G4SubtractionSolid("camera_shield_solid", camera_shield, camera_shield_hole, rotholes, trhole2);
    camera_shield = new G4SubtractionSolid("camera_shield_solid", camera_shield, camera_shield_hole, rotholes, trhole3);
    camera_shield = new G4SubtractionSolid("camera_shield_solid", camera_shield, camera_shield_hole, rotholes, trhole4);
    camera_shield = new G4SubtractionSolid("camera_shield_solid", camera_shield, camera_shield_hole, rotholes, trhole5);
    camera_shield = new G4SubtractionSolid("camera_shield_solid", camera_shield, camera_shield_hole, rotholes, trhole6);
    camera_shield = new G4SubtractionSolid("camera_shield_solid", camera_shield, camera_shield_hole, rotholes, trhole7);
    camera_shield = new G4SubtractionSolid("camera_shield_solid", camera_shield, camera_shield_hole, rotholes, trhole8);
    //FIXME set to Air for external simulations for consistency with old simulation
    camera_shield_log = new G4LogicalVolume(camera_shield,CYGNOMaterials->Material("Cu"),"camera_shield_log",0,0,0) ;
    //camera_shield_log = new G4LogicalVolume(camera_shield,CYGNOMaterials->Material("Air"),"camera_shield_log",0,0,0) ;


    //**********************************************************************
    // ********* CYGNO volumes form CADMesh *****************************
    //**********************************************************************
    
    char namestl[50];
    sprintf(namestl,"%s/shell.stl",CYGNOGeomPath.c_str());
    G4cout << namestl << G4endl;
    //CADMesh * mesh_shell = new CADMesh("../geometry/v2/shell.stl");    
    ifstream infile(CYGNOGeomPath.c_str());
    if (infile.good())
    	mesh_shell = new CADMesh(namestl);    
    //namestl = "camera_carter.stl";
    sprintf(namestl,"%s/camera_carter.stl",CYGNOGeomPath.c_str());
    G4cout << namestl << G4endl;
    if (infile.good())
    	mesh_camera_carter = new CADMesh(namestl);
    //namestl = "cameras.stl";    
    sprintf(namestl,"%s/cameras.stl",CYGNOGeomPath.c_str());
    G4cout << namestl << G4endl;
    if (infile.good())
      mesh_camera = new CADMesh(namestl);    
    //CADMesh * mesh_window = new CADMesh("../geometry/v1/glass_windows.stl");  
    //namestl = "FC_support.stl";  
    sprintf(namestl,"%s/FC_support.stl",CYGNOGeomPath.c_str());
    G4cout << namestl << G4endl;
    if (infile.good())
      mesh_fc_support = new CADMesh(namestl);
    //namestl = "turns_support.stl";    
    sprintf(namestl,"%s/turns_support.stl",CYGNOGeomPath.c_str());
    G4cout << namestl << G4endl;
    if (infile.good())
      mesh_turns_support = new CADMesh(namestl);
    //namestl = "field_cage.stl";    
    sprintf(namestl,"%s/field_cage.stl",CYGNOGeomPath.c_str());
    G4cout << namestl << G4endl;
    if (infile.good())
      mesh_field_cage = new CADMesh(namestl);
    //namestl = "gem_frame.stl";    
    sprintf(namestl,"%s/gem_frame.stl",CYGNOGeomPath.c_str());
    G4cout << namestl << G4endl;
    if (infile.good())
      mesh_gem_support = new CADMesh(namestl);    
    sprintf(namestl,"%s/gem_structure.stl",CYGNOGeomPath.c_str());
    G4cout << namestl << G4endl;
    if (infile.good())
      mesh_gem_structure = new CADMesh(namestl);    
    //namestl = "gem.stl";
    sprintf(namestl,"%s/gem.stl",CYGNOGeomPath.c_str());
    G4cout << namestl << G4endl;
    if (infile.good())
      mesh_gem = new CADMesh(namestl);
    //namstl = "cathode_frame.stl";    
    sprintf(namestl,"%s/cathode_frame.stl",CYGNOGeomPath.c_str());
    G4cout << namestl << G4endl;
    if (infile.good())
      mesh_cathode_frame = new CADMesh(namestl);
    //namestl = "cathode.stl";   
    sprintf(namestl,"%s/cathode.stl",CYGNOGeomPath.c_str());
    G4cout << namestl << G4endl;
    if (infile.good())
      mesh_cathode = new CADMesh(namestl);   

    if (infile.good()){
      mesh_shell->SetScale(mm);
      mesh_camera_carter->SetScale(mm);
      mesh_camera->SetScale(mm);
      //mesh_window->SetScale(mm);
      mesh_fc_support->SetScale(mm);
      mesh_field_cage->SetScale(mm);
      mesh_gem_support->SetScale(mm);
      mesh_gem_structure->SetScale(mm);
      mesh_gem->SetScale(mm);
      mesh_cathode_frame->SetScale(mm);
      mesh_cathode->SetScale(mm);
    

      //shell
      cad_shell_solid = mesh_shell->TessellatedMesh();
      cad_shell_logical = new G4LogicalVolume(cad_shell_solid, CYGNOMaterials->Material("Perspex"), "cad_shell_logical");
      cad_shell_logical->SetVisAttributes(CYGNOMaterials->VisAttributes("Perspex"));


      //camera carter
      cad_camera_carter_solid = mesh_camera_carter->TessellatedMesh();
      cad_camera_carter_logical = new G4LogicalVolume(cad_camera_carter_solid, CYGNOMaterials->Material("Perspex"), "cad_camera_carter_logical", 0, 0, 0);
      //cad_camera_carter_logical->SetVisAttributes(CYGNOMaterials->VisAttributes(cad_camera_carter_logical->GetMaterial()->GetName()));
      cad_camera_carter_logical->SetVisAttributes(CYGNOMaterials->VisAttributes("Perspex"));

      //------- NOT PLACED, used cameras with simple geant volumes ---------
      //cameras
      cad_cameras_all_solid = mesh_camera->TessellatedMesh();
      cad_cameras_all_logical = new G4LogicalVolume(cad_cameras_all_solid, CYGNOMaterials->Material("Camera"), "cad_cameras_all_logical", 0, 0, 0);
      //cad_cameras_all_logical->SetVisAttributes(CYGNOMaterials->VisAttributes(cad_cameras_all_logical->GetMaterial()->GetName()));
      cad_cameras_all_logical->SetVisAttributes(CYGNOMaterials->VisAttributes("Camera"));
      //----------------------------------------------------------------------

      //window
      //cad_window_solid = mesh_window->TessellatedMesh();
      //cad_window_logical = new G4LogicalVolume(cad_window_solid, CYGNOMaterials->Material("Silica"), "cad_window_logical", 0, 0, 0);
    }
  
    //TPC gas
    G4double TPC_x = 1248.*mm;
    G4double TPC_y = 1258.*mm;
    G4double TPC_z = 1310.*mm;
      
    name_phys="TPC";
    name_log=name_phys+"_log";
    name_solid=name_phys+"_solid";
    G4Box* TPC_box = new G4Box(name_solid,0.5*TPC_x,0.5*TPC_y,0.5*TPC_z);
    TPC_log = new G4LogicalVolume(TPC_box,CYGNOMaterials->Material("CYGNO_gas"),name_log,0,0,0);
    
    //CYGNO fiducial gas
    G4double CYGNO_x = 510.*mm;
    G4double CYGNO_y = 1020.*mm;
    G4double CYGNO_z = 1020.*mm;
      
    name_phys="CYGNO";
    name_log=name_phys+"_log";
    name_solid=name_phys+"_solid";
    G4Box* CYGNO_box = new G4Box(name_solid,0.5*CYGNO_x,0.5*CYGNO_y,0.5*CYGNO_z);
    CYGNO_log = new G4LogicalVolume(CYGNO_box,CYGNOMaterials->Material("CYGNO_gas"),name_log,0,0,0);
    //CYGNO_log->SetVisAttributes(CYGNOMaterials->VisAttributes(CYGNO_log->GetMaterial()->GetName()));
    CYGNO_log->SetVisAttributes(CYGNOMaterials->VisAttributes("CYGNO_gas"));

    if (infile.good()){
      //fc support
      cad_fc_support_solid = mesh_fc_support->TessellatedMesh();
      cad_fc_support_logical = new G4LogicalVolume(cad_fc_support_solid, CYGNOMaterials->Material("Perspex"), "cad_fc_support_logical", 0, 0, 0);
      //cad_fc_support_logical->SetVisAttributes(CYGNOMaterials->VisAttributes(cad_fc_support_logical->GetMaterial()->GetName()));
      cad_fc_support_logical->SetVisAttributes(CYGNOMaterials->VisAttributes("Perspex"));
 
      //turns support
      cad_turns_support_solid = mesh_turns_support->TessellatedMesh();
      cad_turns_support_logical = new G4LogicalVolume(cad_turns_support_solid, CYGNOMaterials->Material("Perspex"), "cad_turns_support_logical", 0, 0, 0);
      //cad_turns_support_logical->SetVisAttributes(CYGNOMaterials->VisAttributes(cad_turns_support_logical->GetMaterial()->GetName()));
      cad_turns_support_logical->SetVisAttributes(CYGNOMaterials->VisAttributes("Perspex"));

      //field cage  
      cad_field_cage_solid = mesh_field_cage->TessellatedMesh();
      cad_field_cage_logical = new G4LogicalVolume(cad_field_cage_solid, CYGNOMaterials->Material("Cu"), "cad_field_cage_logical", 0, 0, 0);
      //cad_field_cage_logical->SetVisAttributes(CYGNOMaterials->VisAttributes(cad_field_cage_logical->GetMaterial()->GetName()));
      cad_field_cage_logical->SetVisAttributes(CYGNOMaterials->VisAttributes("Cu"));

      //GEM support
      cad_gem_support_solid = mesh_gem_support->TessellatedMesh();
      cad_gem_support_logical = new G4LogicalVolume(cad_gem_support_solid, CYGNOMaterials->Material("Perspex"), "cad_gem_support_logical", 0, 0, 0);
      //cad_gem_support_logical->SetVisAttributes(CYGNOMaterials->VisAttributes(cad_gem_support_logical->GetMaterial()->GetName()));
      cad_gem_support_logical->SetVisAttributes(CYGNOMaterials->VisAttributes("Perspex"));
      
      //GEM structure support
      cad_gem_structure_solid = mesh_gem_structure->TessellatedMesh();
      cad_gem_structure_logical = new G4LogicalVolume(cad_gem_structure_solid, CYGNOMaterials->Material("Perspex"), "cad_gem_structure_logical", 0, 0, 0);
      cad_gem_structure_logical->SetVisAttributes(CYGNOMaterials->VisAttributes("Perspex"));


      //GEM
      cad_gem_solid = mesh_gem->TessellatedMesh();
      cad_gem_logical = new G4LogicalVolume(cad_gem_solid, CYGNOMaterials->Material("Kapton"), "cad_gem_logical", 0, 0, 0);
      //cad_gem_logical->SetVisAttributes(CYGNOMaterials->VisAttributes(cad_gem_logical->GetMaterial()->GetName()));
      cad_gem_logical->SetVisAttributes(CYGNOMaterials->VisAttributes("Kapton"));
      
      //cathode frame
      cad_cathode_frame_solid = mesh_cathode_frame->TessellatedMesh();
      cad_cathode_frame_logical = new G4LogicalVolume(cad_cathode_frame_solid, CYGNOMaterials->Material("Cu"), "cad_cathode_frame_logical", 0, 0, 0);
      //cad_cathode_frame_logical->SetVisAttributes(CYGNOMaterials->VisAttributes(cad_cathode_frame_logical->GetMaterial()->GetName()));
      cad_cathode_frame_logical->SetVisAttributes(CYGNOMaterials->VisAttributes("Cu"));
      
      //cathode frame
      cad_cathode_solid = mesh_cathode->TessellatedMesh();
      cad_cathode_logical = new G4LogicalVolume(cad_cathode_solid, CYGNOMaterials->Material("Cu"), "cad_cathode_logical", 0, 0, 0);
      cad_cathode_logical->SetVisAttributes(CYGNOMaterials->VisAttributes("Cu"));
    }

    tr_CYGNO_gas_1=G4ThreeVector(0.5*CYGNO_x+1.*mm,0.,0.);	  
    tr_CYGNO_gas_2=G4ThreeVector(-0.5*CYGNO_x-1.*mm,0.,0.);	  

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
    
    TPC_phys = new G4PVPlacement(G4Transform3D(rot,tr),
      	    TPC_log,"TPC_gas", AirBox_log, false, 0, true);
    CYGNO_phys = new G4PVPlacement(G4Transform3D(rot,tr_CYGNO_gas_1),
      	    CYGNO_log,"CYGNO_gas", TPC_log, false, 0, true);
    CYGNO_phys = new G4PVPlacement(G4Transform3D(rot,tr_CYGNO_gas_2),
  		    CYGNO_log,"CYGNO_gas", TPC_log, false, 1, true);
    
    if (infile.good()){
      cad_shell_physical = new G4PVPlacement(G4Transform3D(rot,tr), 
  		    cad_shell_logical,"cad_shell_physical", AirBox_log, false, 0, true);
      cad_camera_carter_physical = new G4PVPlacement(G4Transform3D(rot,tr), 
  		    cad_camera_carter_logical,"cad_camera_carter_physical", AirBox_log, false, 0, true);
     // cad_cameras_all_physical = new G4PVPlacement(G4Transform3D(rot,tr), 
    //	    cad_cameras_all_logical,"cad_cameras_all_physical", AirBox_log, false, 0, true);
  //    //cad_window_physical = new G4PVPlacement(G4Transform3D(rot_cad,tr_cad), 
  ////		    cad_window_logical,"cad_window_physical", Laboratory_log, false, 0, true);
         
      tr=G4ThreeVector(0.,0.,0.);
      tr_cad=G4ThreeVector(0.,0.,0.);
      rot = G4RotationMatrix();
      cad_fc_support_physical = new G4PVPlacement(G4Transform3D(rot,tr), 
  		    cad_fc_support_logical,"cad_fc_support_physical", TPC_log, false, 0, true);
      cad_turns_support_physical = new G4PVPlacement(G4Transform3D(rot,tr), 
  		    cad_turns_support_logical,"cad_turns_support_physical", TPC_log, false, 0, true);
      cad_field_cage_physical = new G4PVPlacement(G4Transform3D(rot,tr), 
  		    cad_field_cage_logical,"cad_field_cage_physical", TPC_log, false, 0, true);
      cad_gem_support_physical = new G4PVPlacement(G4Transform3D(rot,tr), 
  		    cad_gem_support_logical,"cad_gem_support_physical", TPC_log, false, 0, true);
      cad_gem_structure_physical = new G4PVPlacement(G4Transform3D(rot,tr), 
  		    cad_gem_structure_logical,"cad_gem_structure_physical", TPC_log, false, 0, true);
      cad_gem_physical = new G4PVPlacement(G4Transform3D(rot,tr), 
  		    cad_gem_logical,"cad_gem_physical", TPC_log, false, 0, true);
      cad_cathode_frame_physical = new G4PVPlacement(G4Transform3D(rot,tr), 
  		    cad_cathode_frame_logical,"cad_cathode_frame_physical", TPC_log, false, 0, true);
      cad_cathode_physical = new G4PVPlacement(G4Transform3D(rot,tr), 
  		    cad_cathode_logical,"cad_cathode_physical", TPC_log, false, 0, true);
    }  
   

    G4double ztr_cam = 1200.*mm; 
    G4double zmintr_cam = -1200.*mm; 
    G4ThreeVector trcam0(zmintr_cam,0.,0.);
    G4ThreeVector trcam1(zmintr_cam,0.,-354*mm);
    G4ThreeVector trcam2(zmintr_cam,0.,354*mm);
    G4ThreeVector trcam3(zmintr_cam,-354*mm,0.);
    G4ThreeVector trcam4(zmintr_cam,-354*mm,-354*mm);
    G4ThreeVector trcam5(zmintr_cam,-354*mm,354*mm);
    G4ThreeVector trcam6(zmintr_cam,354*mm,0.);
    G4ThreeVector trcam7(zmintr_cam,354*mm,-354*mm);
    G4ThreeVector trcam8(zmintr_cam,354*mm,354*mm);

    G4ThreeVector trlens0(zmintr_cam+0.5*x_camera_body+0.5*z_camera_lens,0.,0.);
    G4ThreeVector trlens1(zmintr_cam+0.5*x_camera_body+0.5*z_camera_lens,0.,-354*mm);
    G4ThreeVector trlens2(zmintr_cam+0.5*x_camera_body+0.5*z_camera_lens,0.,354*mm);
    G4ThreeVector trlens3(zmintr_cam+0.5*x_camera_body+0.5*z_camera_lens,-354*mm,0.);
    G4ThreeVector trlens4(zmintr_cam+0.5*x_camera_body+0.5*z_camera_lens,-354*mm,-354*mm);
    G4ThreeVector trlens5(zmintr_cam+0.5*x_camera_body+0.5*z_camera_lens,-354*mm,354*mm);
    G4ThreeVector trlens6(zmintr_cam+0.5*x_camera_body+0.5*z_camera_lens,354*mm,0.);
    G4ThreeVector trlens7(zmintr_cam+0.5*x_camera_body+0.5*z_camera_lens,354*mm,-354*mm);
    G4ThreeVector trlens8(zmintr_cam+0.5*x_camera_body+0.5*z_camera_lens,354*mm,354*mm);
    
    G4RotationMatrix* rotcam0 = new G4RotationMatrix();
    G4RotationMatrix* rotlens0 = new G4RotationMatrix();
    rotlens0->rotateY(90*deg);
    camera_phys = new G4PVPlacement(rotcam0,trcam0,camera_log,"camera",AirBox_log, false, 0, true);	
    camera_phys = new G4PVPlacement(rotcam0,trcam1,camera_log,"camera",AirBox_log, false, 1, true);	
    camera_phys = new G4PVPlacement(rotcam0,trcam2,camera_log,"camera",AirBox_log, false, 2, true);	
    camera_phys = new G4PVPlacement(rotcam0,trcam3,camera_log,"camera",AirBox_log, false, 3, true);	
    camera_phys = new G4PVPlacement(rotcam0,trcam4,camera_log,"camera",AirBox_log, false, 4, true);	
    camera_phys = new G4PVPlacement(rotcam0,trcam5,camera_log,"camera",AirBox_log, false, 5, true);	
    camera_phys = new G4PVPlacement(rotcam0,trcam6,camera_log,"camera",AirBox_log, false, 6, true);	
    camera_phys = new G4PVPlacement(rotcam0,trcam7,camera_log,"camera",AirBox_log, false, 7, true);	
    camera_phys = new G4PVPlacement(rotcam0,trcam8,camera_log,"camera",AirBox_log, false, 8, true);	
    
    camera_lens_phys = new G4PVPlacement(rotlens0,trlens0,camera_lens_log,"camera_lens",AirBox_log, false, 0, true);	
    camera_lens_phys = new G4PVPlacement(rotlens0,trlens1,camera_lens_log,"camera_lens",AirBox_log, false, 1, true);	
    camera_lens_phys = new G4PVPlacement(rotlens0,trlens2,camera_lens_log,"camera_lens",AirBox_log, false, 2, true);	
    camera_lens_phys = new G4PVPlacement(rotlens0,trlens3,camera_lens_log,"camera_lens",AirBox_log, false, 3, true);	
    camera_lens_phys = new G4PVPlacement(rotlens0,trlens4,camera_lens_log,"camera_lens",AirBox_log, false, 4, true);	
    camera_lens_phys = new G4PVPlacement(rotlens0,trlens5,camera_lens_log,"camera_lens",AirBox_log, false, 5, true);	
    camera_lens_phys = new G4PVPlacement(rotlens0,trlens6,camera_lens_log,"camera_lens",AirBox_log, false, 6, true);	
    camera_lens_phys = new G4PVPlacement(rotlens0,trlens7,camera_lens_log,"camera_lens",AirBox_log, false, 7, true);	
    camera_lens_phys = new G4PVPlacement(rotlens0,trlens8,camera_lens_log,"camera_lens",AirBox_log, false, 8, true);	
    
    G4ThreeVector  trcam9(ztr_cam,0.,0.);
    G4ThreeVector trcam10(ztr_cam,0.,-354*mm);
    G4ThreeVector trcam11(ztr_cam,0.,354*mm);
    G4ThreeVector trcam12(ztr_cam,-354*mm,0.);
    G4ThreeVector trcam13(ztr_cam,-354*mm,-354*mm);
    G4ThreeVector trcam14(ztr_cam,-354*mm,354*mm);
    G4ThreeVector trcam15(ztr_cam,354*mm,0.);
    G4ThreeVector trcam16(ztr_cam,354*mm,-354*mm);
    G4ThreeVector trcam17(ztr_cam,354*mm,354*mm);

    G4ThreeVector trlens9(ztr_cam -0.5*x_camera_body-0.5*z_camera_lens,0.,0.);
    G4ThreeVector trlens10(ztr_cam-0.5*x_camera_body-0.5*z_camera_lens,0.,-354*mm);
    G4ThreeVector trlens11(ztr_cam-0.5*x_camera_body-0.5*z_camera_lens,0.,354*mm);
    G4ThreeVector trlens12(ztr_cam-0.5*x_camera_body-0.5*z_camera_lens,-354*mm,0.);
    G4ThreeVector trlens13(ztr_cam-0.5*x_camera_body-0.5*z_camera_lens,-354*mm,-354*mm);
    G4ThreeVector trlens14(ztr_cam-0.5*x_camera_body-0.5*z_camera_lens,-354*mm,354*mm);
    G4ThreeVector trlens15(ztr_cam-0.5*x_camera_body-0.5*z_camera_lens,354*mm,0.);
    G4ThreeVector trlens16(ztr_cam-0.5*x_camera_body-0.5*z_camera_lens,354*mm,-354*mm);
    G4ThreeVector trlens17(ztr_cam-0.5*x_camera_body-0.5*z_camera_lens,354*mm,354*mm);
   
    G4RotationMatrix* rotcam1 = new G4RotationMatrix();
    //rotcam1->rotateY(180*deg);
    G4RotationMatrix* rotlens1 = new G4RotationMatrix();
    rotlens1->rotateY(90.*deg);
    camera_phys = new G4PVPlacement(rotcam1,trcam9,camera_log,"camera", AirBox_log, false, 9, true);	
    camera_phys = new G4PVPlacement(rotcam1,trcam10,camera_log,"camera",AirBox_log, false, 10, true);	
    camera_phys = new G4PVPlacement(rotcam1,trcam11,camera_log,"camera",AirBox_log, false, 11, true);	
    camera_phys = new G4PVPlacement(rotcam1,trcam12,camera_log,"camera",AirBox_log, false, 12, true);	
    camera_phys = new G4PVPlacement(rotcam1,trcam13,camera_log,"camera",AirBox_log, false, 13, true);	
    camera_phys = new G4PVPlacement(rotcam1,trcam14,camera_log,"camera",AirBox_log, false, 14, true);	
    camera_phys = new G4PVPlacement(rotcam1,trcam15,camera_log,"camera",AirBox_log, false, 15, true);	
    camera_phys = new G4PVPlacement(rotcam1,trcam16,camera_log,"camera",AirBox_log, false, 16, true);	
    camera_phys = new G4PVPlacement(rotcam1,trcam17,camera_log,"camera",AirBox_log, false, 17, true);	
    
    camera_lens_phys = new G4PVPlacement(rotlens1,trlens9,camera_lens_log,"camera_lens",AirBox_log, false, 9, true);	
    camera_lens_phys = new G4PVPlacement(rotlens1,trlens10,camera_lens_log,"camera_lens",AirBox_log, false, 10, true);	
    camera_lens_phys = new G4PVPlacement(rotlens1,trlens11,camera_lens_log,"camera_lens",AirBox_log, false, 11, true);	
    camera_lens_phys = new G4PVPlacement(rotlens1,trlens12,camera_lens_log,"camera_lens",AirBox_log, false, 12, true);	
    camera_lens_phys = new G4PVPlacement(rotlens1,trlens13,camera_lens_log,"camera_lens",AirBox_log, false, 13, true);	
    camera_lens_phys = new G4PVPlacement(rotlens1,trlens14,camera_lens_log,"camera_lens",AirBox_log, false, 14, true);	
    camera_lens_phys = new G4PVPlacement(rotlens1,trlens15,camera_lens_log,"camera_lens",AirBox_log, false, 15, true);	
    camera_lens_phys = new G4PVPlacement(rotlens1,trlens16,camera_lens_log,"camera_lens",AirBox_log, false, 16, true);	
    camera_lens_phys = new G4PVPlacement(rotlens1,trlens17,camera_lens_log,"camera_lens",AirBox_log, false, 17, true);	
    
    camera_shield_phys = new G4PVPlacement(0,G4ThreeVector(zmintr_cam+0.5*x_camera_body+0.5*z_camera_lens,0.,0.),camera_shield_log,"camera_shield",AirBox_log,false, 0, true) ;
    camera_shield_phys = new G4PVPlacement(0,G4ThreeVector(ztr_cam-0.5*x_camera_body-0.5*z_camera_lens,0.,0.),camera_shield_log,"camera_shield",AirBox_log,false, 1, true) ;

    
    //
    //**********************************************************************
    // GLOBAL TRANSLATIONS ***************************************
    G4cout<<"Placement of Laboratory in the World started"<<G4endl;
    
    tr_Rock=-1*(tr_Laboratory+rot_Laboratory*(tr_cad+rot_cad*(tr_InsideVolume+rot_InsideVolume*(tr_Shielding))));//The shift of Rock_log in the world volume to make the origin be the center of the detector

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
    CYGNOProperties->AddVolumeNameMassAndDensity(cad_shell_logical);
    CYGNOProperties->AddVolumeNameMassAndDensity(cad_camera_carter_logical);
    CYGNOProperties->AddVolumeNameMassAndDensity(cad_cameras_all_logical);
    CYGNOProperties->AddVolumeNameMassAndDensity(camera_log);
    CYGNOProperties->AddVolumeNameMassAndDensity(camera_lens_log);
    CYGNOProperties->AddVolumeNameMassAndDensity(camera_shield_log);
    CYGNOProperties->AddVolumeNameMassAndDensity(TPC_log);
    CYGNOProperties->AddVolumeNameMassAndDensity(CYGNO_log);
    CYGNOProperties->AddVolumeNameMassAndDensity(cad_fc_support_logical);
    CYGNOProperties->AddVolumeNameMassAndDensity(cad_turns_support_logical);
    CYGNOProperties->AddVolumeNameMassAndDensity(cad_gem_support_logical);
    CYGNOProperties->AddVolumeNameMassAndDensity(cad_gem_structure_logical);
    CYGNOProperties->AddVolumeNameMassAndDensity(cad_gem_logical);
    CYGNOProperties->AddVolumeNameMassAndDensity(cad_cathode_frame_logical);
    CYGNOProperties->AddVolumeNameMassAndDensity(cad_cathode_logical);
    CYGNOProperties->AddVolumeNameMassAndDensity(cad_field_cage_logical);
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

  cad_shell_logical=0;
  cad_camera_carter_logical=0;
  cad_cameras_all_logical=0;
  TPC_log=0;
  CYGNO_log=0;
  cad_fc_support_logical=0;
  cad_turns_support_logical=0;
  cad_gem_support_logical=0;
  cad_gem_logical=0;
  cad_cathode_frame_logical=0;
  cad_cathode_logical=0;
  cad_field_cage_logical=0;
  camera_log=0;
  camera_lens_log=0;
  camera_shield_log=0;


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
