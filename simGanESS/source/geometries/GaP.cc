#include "GaP.hh"

#include <G4Box.hh>
#include <G4Tubs.hh>
#include <G4SubtractionSolid.hh>
#include <G4UnionSolid.hh>
#include <G4MultiUnion.hh>
#include <G4OpticalSurface.hh>
#include <G4LogicalSkinSurface.hh>
#include "G4LogicalVolume.hh"
#include "G4NistManager.hh"
#include "G4PVPlacement.hh"
#include "G4UnitsTable.hh"
#include <G4UserLimits.hh>
#include <G4SDManager.hh>

#include "nexus/FactoryBase.h"
#include "nexus/UniformElectricDriftField.h"
#include "nexus/OpticalMaterialProperties.h"
#include "nexus/XenonProperties.h"
#include "nexus/IonizationSD.h"

#include <iostream>
#include <cmath>
#include <vector>


using namespace nexus;

REGISTER_CLASS(GaP, GeometryBase)

GaP::GaP():
    GeometryBase(),
    msg_ (nullptr),
    gas_ (nullptr),
    mesh_mat_ (nullptr),

    vessel_out_rad_    (288./2  *mm),
    vessel_out_length_ (46.679  *cm),
    vessel_rad_        (276./2  *mm),
    vessel_length_     (38.599  *cm), // Adjusted length so that the gas volume is centered. Original length (38.639  *cm),

    mesh_rad_          (104./2  *mm),
    mesh_thickn_       (0.075   *mm),
    mesh_transparency_ (0.95),

    meshBracket_rad_      (180./2  *mm),
    meshBracket_thickn_   (6.      *mm),
    anodeBracket_rad_     (160./2  *mm),
    anodeBracket_thickn_  (6.975   *mm),

    pmt_rad_ (25.4/2  *mm),

    enclosure_pmt_rad_        (120./2  *mm),
    enclosure_pmt_thickn_     (8.5     *mm),
    enclosure_pmt_length_     (113.5   *mm),
    enclosurevac_pmt_length_  (110.5   *mm),

    plate_pmt_rad_        (105./2  *mm),
    plate_pmt_thickn_     (105./2  *mm),
    plate_pmt_length_     (10      *mm),
    plateUp_pmt_length_   (15      *mm),
    plateUp_pmt_thickn_   (21.5    *mm),

    pmtHolder_rad_        (115./2  *mm),
    pmtHolder_length_     (9.      *mm),

    quartz_window_rad_    (108./2  *mm),
    quartz_window_thickn_ (3       *mm),
    tpb_coating_thickn_   (3       *micrometer),

    photoe_prob_       (0.),

    pressure_          (10.* bar),
    temperature_       (293. * kelvin),
    sc_yield_          (22222 * 1/MeV), // Wsc = 45 eV, fr
    elifetime_         (1e6* ms),
    drift_vel_         (1. * mm/microsecond),
    drift_transv_diff_ (1. * mm/sqrt(cm)),
    drift_long_diff_   (.3 * mm/sqrt(cm)),
    el_field_          (16.0 * kilovolt/cm),
    el_vel_            (3. * mm/microsecond),
    el_transv_diff_    (1. * mm/sqrt(cm)),
    el_long_diff_      (.3 * mm/sqrt(cm)),
    specific_vertex_{}
{
  // Messenger
  msg_ = new G4GenericMessenger(this, "/Geometry/GaP/",
                                "Control commands of the GaP geometry.");
  // Parametrized dimensions
  DefineConfigurationParameters();
}

GaP::~GaP()
{
  delete msg_;

  delete drift_gen_;
  delete el_gen_;

}

void GaP::Construct()
{
    //Materials
    gas_   = materials::GXe(pressure_, temperature_);
    gas_->SetMaterialPropertiesTable(opticalprops::GXe(pressure_,
                                                      temperature_,
                                                      sc_yield_,
                                                      elifetime_));
    // Mesh materials (cathode, anode and gate)
    mesh_mat_ = materials::FakeDielectric(gas_, "mesh_mat");
    mesh_mat_->SetMaterialPropertiesTable(opticalprops::FakeGrid(pressure_,
              temperature_, mesh_transparency_, mesh_thickn_,
              sc_yield_, elifetime_, photoe_prob_));

    steel_ = materials::Steel();
    steel_->SetMaterialPropertiesTable(new G4MaterialPropertiesTable());

    //PEEK for the holders
    peek_ = materials::PEEK();
    peek_->SetMaterialPropertiesTable(new G4MaterialPropertiesTable());

    // quartz (SiO2, crystalline) with evaporated TPB
    quartz_ = materials::FusedSilica();
    quartz_->SetMaterialPropertiesTable(opticalprops::FusedSilica());

    tpb_ = materials::TPB();
    tpb_->SetMaterialPropertiesTable(opticalprops::TPB());

    vacuum_ = G4NistManager::Instance()->FindOrBuildMaterial("G4_Galactic");


    //Cylinder, acting as the vessel
    G4Tubs          *solid_vessel_steel = new G4Tubs("Vessel", 0, vessel_out_rad_, vessel_out_length_/2 , 0., 360.*deg);
    G4LogicalVolume *logic_vessel_steel = new G4LogicalVolume(solid_vessel_steel, steel_, "Vessel");
    this->SetLogicalVolume(logic_vessel_steel);

    //Build inside detector
    BuildTPC(gas_, mesh_mat_, steel_, peek_, vacuum_, quartz_, tpb_, logic_vessel_steel);
}

G4ThreeVector GaP::GenerateVertex(const G4String& region) const
{
  G4ThreeVector vertex;

  if (region == "AD_HOC") {
    vertex = specific_vertex_;
  }

  //// Gas regions
  else if (
    (region == "GasEL") ||
    (region == "GasDrift")) {
    vertex = GenerateVertexGas(region);
  }

  else {
    G4Exception("[GaP]", "GenerateVertex()", FatalException,
      "Unknown vertex generation region!");
  }

  return vertex;
}

G4ThreeVector GaP::GenerateVertexGas(const G4String& region) const
{
    G4ThreeVector vertex;

    if     (region == "GasEL")     {vertex = el_gen_->GenerateVertex("VOLUME");}
    else if(region == "GasDrift")  {vertex = drift_gen_->GenerateVertex("VOLUME");}
    else{G4Exception("[GaP]", "GenerateVertex()", FatalException,
                "Unknown vertex generation region!");}
    return vertex;
}


void GaP::DefineConfigurationParameters()
{
  // Gas pressure
  G4GenericMessenger::Command& pressure_cmd =
    msg_->DeclareProperty("pressure", pressure_,
                          "Pressure of the gas.");
  pressure_cmd.SetUnitCategory("Pressure");
  pressure_cmd.SetParameterName("pressure", false);
  pressure_cmd.SetRange("pressure>0.");

  // Gas temperature
  G4GenericMessenger::Command& temperature_cmd =
    msg_->DeclareProperty("temperature", temperature_,
                          "Temperature of the gas.");
  temperature_cmd.SetUnitCategory("Temperature");
  temperature_cmd.SetParameterName("temperature", false);
  temperature_cmd.SetRange("temperature>0.");

  // e- lifetime
  G4GenericMessenger::Command& e_lifetime_cmd =
    msg_->DeclareProperty("elifetime", elifetime_,
                          "Electron lifetime in gas.");
  e_lifetime_cmd.SetParameterName("elifetime", false);
  e_lifetime_cmd.SetUnitCategory("Time");
  e_lifetime_cmd.SetRange("elifetime>0.");

  // Drift velocity in drift region
  G4GenericMessenger::Command& drift_vel_cmd =
    msg_->DeclareProperty("drift_vel", drift_vel_,
                          "Electron drift velocity in the drift region.");
  drift_vel_cmd.SetParameterName("drift_vel", false);
  drift_vel_cmd.SetUnitCategory("Velocity");
  drift_vel_cmd.SetRange("drift_vel>0.");

  // Transverse diffusion in drift region
  new G4UnitDefinition("mm/sqrt(cm)", "mm/sqrt(cm)", "Diffusion", mm/sqrt(cm));
  G4GenericMessenger::Command& drift_transv_diff_cmd =
    msg_->DeclareProperty("drift_transv_diff", drift_transv_diff_,
                          "Tranvsersal diffusion in the drift region");
  drift_transv_diff_cmd.SetParameterName("drift_transv_diff", false);
  drift_transv_diff_cmd.SetUnitCategory("Diffusion");
  drift_transv_diff_cmd.SetRange("drift_transv_diff>0.");

  // Longitudinal diffusion in drift region
  G4GenericMessenger::Command& drift_long_diff_cmd =
    msg_->DeclareProperty("drift_long_diff", drift_long_diff_,
                          "Longitudinal diffusion in the drift region");
  drift_long_diff_cmd.SetParameterName("drift_long_diff", false);
  drift_long_diff_cmd.SetUnitCategory("Diffusion");
  drift_long_diff_cmd.SetRange("drift_long_diff>0.");

  // Scintillation yield (for S1)
  new G4UnitDefinition("1/MeV","1/MeV", "1/Energy", 1/MeV);
  G4GenericMessenger::Command& sc_yield_cmd =
    msg_->DeclareProperty("sc_yield", sc_yield_,
        "Set scintillation yield for gas. It is in photons/MeV");
  sc_yield_cmd.SetParameterName("sc_yield", true);
  sc_yield_cmd.SetUnitCategory("1/Energy");

  // Drift velocity in EL region
  G4GenericMessenger::Command& el_vel_cmd =
    msg_->DeclareProperty("el_vel", el_vel_,
                          "Electron drift velocity in the EL region.");
  el_vel_cmd.SetParameterName("el_vel", false);
  el_vel_cmd.SetUnitCategory("Velocity");
  el_vel_cmd.SetRange("el_vel>0.");

  // Transverse diffusion in EL region
  G4GenericMessenger::Command& el_transv_diff_cmd =
    msg_->DeclareProperty("el_transv_diff", el_transv_diff_,
                          "Tranvsersal diffusion in the EL region");
  el_transv_diff_cmd.SetParameterName("el_transv_diff", false);
  el_transv_diff_cmd.SetUnitCategory("Diffusion");
  el_transv_diff_cmd.SetRange("el_transv_diff>0.");

  // Longitudinal diffusion in EL region
  G4GenericMessenger::Command& el_long_diff_cmd =
    msg_->DeclareProperty("el_long_diff", el_long_diff_,
                          "Longitudinal diffusion in the EL region");
  el_long_diff_cmd.SetParameterName("el_long_diff", false);
  el_long_diff_cmd.SetUnitCategory("Diffusion");
  el_long_diff_cmd.SetRange("el_long_diff>0.");

  // EL field
  new G4UnitDefinition("kilovolt/cm", "kV/cm", "Electric field", kilovolt/cm);
  G4GenericMessenger::Command& el_field_cmd =
    msg_->DeclareProperty("el_field", el_field_,
                          "EL electric field intensity");
  el_field_cmd.SetParameterName("el_field", false);
  el_field_cmd.SetUnitCategory("Electric field");

  // Photoelectric probability
  msg_->DeclareProperty("photoe_prob", photoe_prob_,
                        "Probability of optical photon to ie- conversion");


  // Specific vertex in case region to shoot from is AD_HOC
  msg_->DeclarePropertyWithUnit("specific_vertex", "mm",  specific_vertex_, "Set generation vertex.");
}

void GaP::BuildTPC(G4Material* gas, G4Material* mesh_mat, G4Material* steel, G4Material* peek, G4Material* vacuum, G4Material* quartz, G4Material* tpb, G4LogicalVolume* logic_vessel_steel)
{
    //Gas
    G4double drift_length_  = 19.825*mm - mesh_thickn_ ;
    G4double el_length_     = 10.775*mm + mesh_thickn_;

    G4Tubs          *solid_vessel = new G4Tubs("GasVessel", 0, vessel_rad_, vessel_length_/2 , 0., 360.*deg);
    G4LogicalVolume *logic_vessel = new G4LogicalVolume(solid_vessel, gas, "GasVessel");
    new G4PVPlacement(0, G4ThreeVector(), logic_vessel, "GasVessel", logic_vessel_steel, false, 0, true);

    //// Cathode
    G4Tubs          *solid_cathode = new G4Tubs("Cathode", 0., mesh_rad_, (mesh_thickn_)/2, 0., 360.*deg);
    G4LogicalVolume *logic_cathode = new G4LogicalVolume(solid_cathode, mesh_mat, "Cathode");
    G4double cathode_z = 4.505*mm + mesh_thickn_/2;  //cathode center from vessel center
    new G4PVPlacement(0, G4ThreeVector(0., 0., cathode_z), logic_cathode, "Cathode", logic_vessel, false, 0, true);

    // Cathode Bracket
    G4Tubs          *solid_cathBracket = new G4Tubs("CathodeBracket", mesh_rad_, meshBracket_rad_, (meshBracket_thickn_)/2, 0., 360.*deg);
    G4LogicalVolume *logic_cathBracket = new G4LogicalVolume(solid_cathBracket, steel, "CathodeBracket");
    G4double cathBracket_z = 8.005*mm - meshBracket_thickn_/2;
    new G4PVPlacement(0, G4ThreeVector(0., 0., cathBracket_z), logic_cathBracket, "CathodeBracket", logic_vessel, false, 0, true);

    //// Drift
    G4Tubs          *solid_gas_drift = new G4Tubs("GasDrift", 0., mesh_rad_, (drift_length_)/2, 0., 360.*deg);
    G4LogicalVolume *logic_gas_drift = new G4LogicalVolume(solid_gas_drift, gas, "GasDrift");

    G4double drift_z = cathode_z - mesh_thickn_/2 - drift_length_/2;
    G4VPhysicalVolume* drift_phys_ = new G4PVPlacement(0, G4ThreeVector(0., 0., drift_z), logic_gas_drift, "GasDrift", logic_vessel, false, 0, true);
    drift_gen_  = new CylinderPointSampler2020(drift_phys_);

    // Define the drift field
    UniformElectricDriftField* drift_field = new UniformElectricDriftField();
    drift_field->SetCathodePosition(drift_z + drift_length_/2);
    drift_field->SetAnodePosition  (drift_z - drift_length_/2);
    drift_field->SetDriftVelocity(drift_vel_);
    drift_field->SetTransverseDiffusion(drift_transv_diff_);
    drift_field->SetLongitudinalDiffusion(drift_long_diff_);
    G4Region* drift_region = new G4Region("DRIFT");
    drift_region->SetUserInformation(drift_field);
    drift_region->AddRootLogicalVolume(logic_gas_drift);

    logic_gas_drift->SetUserLimits(new G4UserLimits(1.*mm));

    // Set the DRIFT volume as an ionization sensitive detector
    IonizationSD* active_sd = new IonizationSD("/GaP/DRIFT");
    logic_gas_drift->SetSensitiveDetector(active_sd);
    G4SDManager::GetSDMpointer()->AddNewDetector(active_sd);

    //// EL gap
    G4Tubs          *solid_gas_el = new G4Tubs("GasEL", 0., mesh_rad_, (el_length_)/2, 0., 360.*deg);
    G4LogicalVolume *logic_gas_el = new G4LogicalVolume(solid_gas_el, gas, "GasEL");

    G4double el_z = drift_z - drift_length_/2 - el_length_/2;
    G4VPhysicalVolume* el_phys_ = new G4PVPlacement(0, G4ThreeVector(0., 0., el_z), logic_gas_el, "GasEL", logic_vessel, false, 0, true);
    el_gen_  = new CylinderPointSampler2020(el_phys_);

    //// Gate
    G4Tubs          *solid_gate = new G4Tubs("Gate", 0., mesh_rad_, (mesh_thickn_)/2, 0., 360.*deg);
    G4LogicalVolume *logic_gate = new G4LogicalVolume(solid_gate, mesh_mat, "Gate");

    G4double gate_z = el_length_/2 - mesh_thickn_/2;
    new G4PVPlacement(0, G4ThreeVector(0., 0., gate_z), logic_gate, "Gate", logic_gas_el, false, 0, true);

    // Gate Bracket
    G4Tubs          *solid_gateBracket = new G4Tubs("GateBracket", mesh_rad_, meshBracket_rad_, (meshBracket_thickn_)/2, 0., 360.*deg);
    G4LogicalVolume *logic_gateBracket = new G4LogicalVolume(solid_gateBracket, steel, "GateBracket");

    G4double gateBracket_z = 12.745*mm + meshBracket_thickn_/2;
    new G4PVPlacement(0, G4ThreeVector(0., 0., -gateBracket_z), logic_gateBracket, "GateBracket", logic_vessel, false, 0, true);

    /// Peek Mesh Holder (holds cathode and gate)
    G4double meshHolder_length_    = 36.75     *mm;
    G4double meshHolder_width_     = 21.035    *mm;
    G4double meshHolder_thickn_    = 24        *mm;

    G4double meshHolder_hole_length_   = 15.75  *mm;
    G4double meshHolder_holeUp_length_ = 6      *mm;
    G4double meshHolder_hole_thickn_   = 11.667 *mm;

    G4Box *solid_meshHolder_block = new G4Box("MeshHolder", meshHolder_width_/2, meshHolder_thickn_/2 , meshHolder_length_/2);
    G4Box *solid_meshHolder_hole = new G4Box("MeshHolderHoleDown", meshHolder_width_, meshHolder_hole_thickn_ , meshHolder_hole_length_/2);
    G4Box *solid_meshHolder_holeUp = new G4Box("MeshHolderHoleUp", meshHolder_width_, meshHolder_hole_thickn_ , meshHolder_holeUp_length_/2);
    G4VSolid        *solidSub_meshHolder0 = new G4SubtractionSolid("MeshHolderHole_SubDown", solid_meshHolder_block, solid_meshHolder_hole, 0, G4ThreeVector(0., meshHolder_thickn_/2, meshHolder_length_/2-meshHolder_hole_length_/2-5*mm) );
    G4VSolid        *solidSub_meshHolder = new G4SubtractionSolid("MeshHolderHole_SubUp", solidSub_meshHolder0, solid_meshHolder_holeUp, 0, G4ThreeVector(0., meshHolder_thickn_/2, -(meshHolder_length_/2-meshHolder_holeUp_length_/2-5*mm)) );
    G4LogicalVolume *logic_meshHolder = new G4LogicalVolume(solidSub_meshHolder, peek, "MeshHolder");

    G4double meshHolder_x = 5*mm + meshBracket_rad_ * cos(45*deg);
    G4double meshHolder_y = 5*mm + meshBracket_rad_ * sin(45*deg);
    G4double meshHolder_z = -13.005*mm + meshHolder_length_/2 ;

    G4RotationMatrix* Rot45 = new G4RotationMatrix();
    Rot45->rotateZ(45*deg);
    G4RotationMatrix* Rot_45 = new G4RotationMatrix();
    Rot_45->rotateZ(-45*deg);
    G4RotationMatrix* Rot135 = new G4RotationMatrix();
    Rot135->rotateZ(135*deg);
    G4RotationMatrix* Rot_135 = new G4RotationMatrix();
    Rot_135->rotateZ(-135*deg);

    new G4PVPlacement(Rot45,   G4ThreeVector(-meshHolder_x, -meshHolder_y, -meshHolder_z), logic_meshHolder, "MeshHolderA", logic_vessel, false, 0, true);
    new G4PVPlacement(Rot135,  G4ThreeVector(-meshHolder_x, meshHolder_y, -meshHolder_z), logic_meshHolder, "MeshHolderB", logic_vessel, false, 1, true);
    new G4PVPlacement(Rot_45,  G4ThreeVector(meshHolder_x, -meshHolder_y, -meshHolder_z), logic_meshHolder, "MeshHolderC", logic_vessel, false, 2, true);
    new G4PVPlacement(Rot_135, G4ThreeVector(meshHolder_x, meshHolder_y, -meshHolder_z), logic_meshHolder, "MeshHolderD", logic_vessel, false, 3, true);

    //Steel Bar joining Mesh holder and PMT clad
    G4double meshHolderBar_rad_     = 9./2   *mm;
    G4double meshHolderBar_length_  = 35.75  *mm;

    G4Tubs          *solid_meshHolderBar = new G4Tubs("MeshHolderBar", 0., meshHolderBar_rad_, (meshHolderBar_length_)/2, 0., 360.*deg);
    G4LogicalVolume *logic_meshHolderBar = new G4LogicalVolume(solid_meshHolderBar, steel, "MeshHolderBar");

    G4double meshHolderBar_xy = 73.769*mm;
    G4double meshHolderBar_z  = meshHolder_z + meshHolder_length_/2 + meshHolderBar_length_/2 ;
    new G4PVPlacement(0, G4ThreeVector(meshHolderBar_xy, meshHolderBar_xy, -meshHolderBar_z), logic_meshHolderBar, "MeshHolderBarA", logic_vessel, false, 0, true);
    new G4PVPlacement(0, G4ThreeVector(-meshHolderBar_xy, meshHolderBar_xy, -meshHolderBar_z), logic_meshHolderBar, "MeshHolderBarB", logic_vessel, false, 1, true);
    new G4PVPlacement(0, G4ThreeVector(meshHolderBar_xy, -meshHolderBar_xy, -meshHolderBar_z), logic_meshHolderBar, "MeshHolderBarC", logic_vessel, false, 2, true);
    new G4PVPlacement(0, G4ThreeVector(-meshHolderBar_xy, -meshHolderBar_xy, -meshHolderBar_z), logic_meshHolderBar, "MeshHolderBarD", logic_vessel, false, 3, true);

    /// Define EL electric field
    G4double yield = XenonELLightYield(el_field_, gas->GetPressure());
    UniformElectricDriftField* el_field = new UniformElectricDriftField();
    el_field->SetCathodePosition(el_z + el_length_/2.);
    el_field->SetAnodePosition  (el_z - el_length_/2.);
    el_field->SetDriftVelocity        (el_vel_);
    el_field->SetTransverseDiffusion  (el_transv_diff_);
    el_field->SetLongitudinalDiffusion(el_long_diff_);
    el_field->SetLightYield(yield);
    G4Region* el_region = new G4Region("EL");
    el_region->SetUserInformation(el_field);
    el_region->AddRootLogicalVolume(logic_gas_el);

    //// Anode
    G4Tubs          *solid_anode = new G4Tubs("Anode", 0., mesh_rad_, (mesh_thickn_)/2, 0., 360.*deg);
    G4LogicalVolume *logic_anode = new G4LogicalVolume(solid_anode, mesh_mat, "Anode");

    G4double anode_z = - el_length_/2 + mesh_thickn_/2;
    new G4PVPlacement(0, G4ThreeVector(0., 0., anode_z), logic_anode, "Anode", logic_gas_el, false, 0, true);

    // Anode Bracket
    G4Tubs          *solid_anodeBracket = new G4Tubs("AnodeBracket", mesh_rad_, anodeBracket_rad_, (anodeBracket_thickn_)/2, 0., 360.*deg);
    G4LogicalVolume *logic_anodeBracket = new G4LogicalVolume(solid_anodeBracket, steel, "AnodeBracket");

    G4double anodeBracket_z = gateBracket_z + meshBracket_thickn_/2 + 3.775*mm + anodeBracket_thickn_/2;
    new G4PVPlacement(0, G4ThreeVector(0., 0., -anodeBracket_z), logic_anodeBracket, "AnodeBracket", logic_vessel, false, 0, true);

    // Peek Anode Holder
    G4double anodeHolder_length_  = 30        *mm;
    G4double anodeHolder_rad_     = 145.001/2 *mm;
    G4double anodeHolder_thickn_  = 17.       *mm;
    G4double anodeHolder_angle_   = 7.964     *deg;

    G4double anodeHolder_hole_length_ = 20    *mm;
    G4double anodeHolder_hole_rad_    = anodeHolder_rad_ + 8*mm;
    G4double anodeHolder_hole_thickn_ = 10    *mm;
    G4double anodeHolder_hole_angle_  = 7.965 *deg;

    G4Tubs *solid_anodeHolder_block = new G4Tubs("AnodeHolder", anodeHolder_rad_, anodeHolder_rad_+anodeHolder_thickn_, anodeHolder_length_/2, -anodeHolder_angle_/2, anodeHolder_angle_);
    G4Tubs *solid_anodeHolder_hole = new G4Tubs("AnodeHolderHole", anodeHolder_hole_rad_, anodeHolder_hole_rad_+anodeHolder_hole_thickn_, anodeHolder_hole_length_/2, -anodeHolder_hole_angle_/2, anodeHolder_hole_angle_);
    G4VSolid        *solidSub_anodeHolder = new G4SubtractionSolid("AnodeHolderHole_Sub", solid_anodeHolder_block, solid_anodeHolder_hole, 0, G4ThreeVector(0,0,anodeHolder_length_/2-anodeHolder_hole_length_/2) );
    G4LogicalVolume *logic_anodeHolder = new G4LogicalVolume(solidSub_anodeHolder, peek, "AnodeHolder");

    G4double anodeHolder_z = 29.495*mm + anodeHolder_length_/2;
    new G4PVPlacement(Rot45,   G4ThreeVector(0., 0., -anodeHolder_z), logic_anodeHolder, "AnodeHolderA", logic_vessel, false, 0, true);
    new G4PVPlacement(Rot135,  G4ThreeVector(0., .0, -anodeHolder_z), logic_anodeHolder, "AnodeHolderB", logic_vessel, false, 1, true);
    new G4PVPlacement(Rot_45,  G4ThreeVector(0., 0., -anodeHolder_z), logic_anodeHolder, "AnodeHolderC", logic_vessel, false, 2, true);
    new G4PVPlacement(Rot_135, G4ThreeVector(0., 0., -anodeHolder_z), logic_anodeHolder, "AnodeHolderD", logic_vessel, false, 3, true);

    G4cout << "* GATE Z position: " << el_z + gate_z << G4endl;
    G4cout << "* GATE Volt position: " << el_z + el_length_/2. << G4endl;
    G4cout << "* ANODE Z position: " << el_z + anode_z << G4endl;
    G4cout << "* ANODE Volt position: " << el_z - el_length_/2. << G4endl;
    G4cout << "* EL_GAP Z positions: " << el_z - el_length_/2. <<
              " to " << el_z + el_length_/2. << G4endl;

    G4cout << "* EL field intensity (kV/cm): " << el_field_ / (kilovolt/cm)
            << "  ->  EL Light yield (photons/ie-/cm): " << yield / (1/cm) << G4endl;
    G4cout << "* EL Light yield (photons/ie-): " << yield * el_length_ << G4endl;

    //// Quartz Window
    G4Tubs          *solid_quartz_window = new G4Tubs("QuartzWindow", 0, quartz_window_rad_, quartz_window_thickn_/2, 0, 360*deg);
    G4LogicalVolume *logic_quartz_window = new G4LogicalVolume(solid_quartz_window, quartz, "QuartzWindow");

    G4double quartz_window_z = 35.495*mm + quartz_window_thickn_/2 ;
    new G4PVPlacement(0, G4ThreeVector(0., 0., -quartz_window_z ), logic_quartz_window, "QuartzWindow", logic_vessel, false, 0, true);

    /// Evaporated TPB
    G4Tubs          *solid_tpb_coating = new G4Tubs("CoatingTPB", 0, quartz_window_rad_, tpb_coating_thickn_/2, 0, 360*deg);
    G4LogicalVolume *logic_tpb_coating = new G4LogicalVolume(solid_tpb_coating, tpb, "CoatingTPB");

    G4double tpb_coating_z = quartz_window_z + quartz_window_thickn_/2 + tpb_coating_thickn_/2 ;
    new G4PVPlacement(0, G4ThreeVector(0., 0., -tpb_coating_z ), logic_tpb_coating, "CoatingTPB", logic_vessel, false, 0, true);

    // Optical surface between gas and TPB to model the latter's roughness
    G4OpticalSurface* tpb_coating_surf =
      new G4OpticalSurface("TPB_COATING_OPSURF", glisur, ground,
                           dielectric_dielectric, .01);
    new G4LogicalSkinSurface("TPB_COATING_OPSURF", logic_tpb_coating, tpb_coating_surf);

    // Peek Quartz Window Holder
    G4double windowHolder_angle_   = 10      *deg;
    G4double windowHolder_length1_ = 4.       *mm;
    G4double windowHolder_length2_ = 19.      *mm;
    G4double windowHolder_length3_ = 5.       *mm;
    G4double windowHolder_rad0_    = 100.4/2  *mm;
    G4double windowHolder_rad1_    = 118.4/2  *mm;
    G4double windowHolder_rad2_    = 122.4/2  *mm;
    G4double windowHolder_rad3_    = 140.4/2  *mm;

    G4double windowHolder_block2_z = windowHolder_length1_/2 + windowHolder_length2_/2 ;
    G4double windowHolder_block3_z = windowHolder_length3_/2 + windowHolder_length2_/2 + windowHolder_block2_z ;

    G4Tubs *solid_windowHolder_block1 = new G4Tubs("QuartzWindowHolder1", windowHolder_rad0_, windowHolder_rad2_, windowHolder_length1_/2, -windowHolder_angle_/2, windowHolder_angle_);
    G4Tubs *solid_windowHolder_block2 = new G4Tubs("QuartzWindowHolder2", windowHolder_rad1_, windowHolder_rad2_, windowHolder_length2_/2, -windowHolder_angle_/2, windowHolder_angle_);
    G4Tubs *solid_windowHolder_block3 = new G4Tubs("QuartzWindowHolder3", windowHolder_rad1_, windowHolder_rad3_, windowHolder_length3_/2, -windowHolder_angle_/2, windowHolder_angle_);
    G4VSolid        *solidUni_windowHolder0 = new G4UnionSolid("WindowHolder_Union1", solid_windowHolder_block1, solid_windowHolder_block2, 0, G4ThreeVector(0,0,-windowHolder_block2_z) );
    G4VSolid        *solidUni_windowHolder = new G4UnionSolid("WindowHolder_Union2", solidUni_windowHolder0, solid_windowHolder_block3, 0, G4ThreeVector(0,0,-windowHolder_block3_z) );
    G4LogicalVolume *logic_windowHolder = new G4LogicalVolume(solidUni_windowHolder, peek, "QuartzWindowHolder");

    G4double quartz_windowHolder_z = quartz_window_z - quartz_window_thickn_/2 -  windowHolder_length1_/2;
    new G4PVPlacement(Rot45,   G4ThreeVector(0., 0., -quartz_windowHolder_z), logic_windowHolder, "QuartzWindowHolderA", logic_vessel, false, 0, true);
    new G4PVPlacement(Rot135,  G4ThreeVector(0., 0., -quartz_windowHolder_z), logic_windowHolder, "QuartzWindowHolderB", logic_vessel, false, 1, true);
    new G4PVPlacement(Rot_45,  G4ThreeVector(0., 0., -quartz_windowHolder_z), logic_windowHolder, "QuartzWindowHolderC", logic_vessel, false, 2, true);
    new G4PVPlacement(Rot_135, G4ThreeVector(0., 0., -quartz_windowHolder_z), logic_windowHolder, "QuartzWindowHolderD", logic_vessel, false, 3, true);

    G4double windowHolderTop_length1_ = 3.  *mm;
    G4double windowHolderTop_length2_ = 4.  *mm;
    G4double windowHolderTop_rad0_    = windowHolder_rad0_ ;
    G4double windowHolderTop_rad1_    = windowHolderTop_rad0_ + 4.*mm ;
    G4double windowHolderTop_rad2_    = windowHolder_rad1_ ;

    G4double windowHolderTop_block2_z = windowHolderTop_length1_/2 + windowHolderTop_length2_/2 ;

    G4Tubs *solid_windowHolderTop_block1 = new G4Tubs("QuartzWindowHolderTop1", windowHolderTop_rad1_, windowHolderTop_rad2_, windowHolderTop_length1_/2, -windowHolder_angle_/2, windowHolder_angle_);
    G4Tubs *solid_windowHolderTop_block2 = new G4Tubs("QuartzWindowHolderTop2", windowHolderTop_rad0_, windowHolderTop_rad2_, windowHolderTop_length2_/2, -windowHolder_angle_/2, windowHolder_angle_);
    G4VSolid        *solidUni_windowHolderTop = new G4UnionSolid("WindowHolderTop_Union", solid_windowHolderTop_block1, solid_windowHolderTop_block2, 0, G4ThreeVector(0,0,-windowHolderTop_block2_z) );
    G4LogicalVolume *logic_windowHolderTop = new G4LogicalVolume(solidUni_windowHolderTop, peek, "QuartzWindowHolderTop");

    G4double quartz_windowHolderTop_z = quartz_window_z - quartz_window_thickn_/2 +  windowHolderTop_length1_/2;
    new G4PVPlacement(Rot45,   G4ThreeVector(0., 0., -quartz_windowHolderTop_z), logic_windowHolderTop, "QuartzWindowHolderTopA", logic_vessel, false, 0, true);
    new G4PVPlacement(Rot135,  G4ThreeVector(0., 0., -quartz_windowHolderTop_z), logic_windowHolderTop, "QuartzWindowHolderTopB", logic_vessel, false, 1, true);
    new G4PVPlacement(Rot_45,  G4ThreeVector(0., 0., -quartz_windowHolderTop_z), logic_windowHolderTop, "QuartzWindowHolderTopC", logic_vessel, false, 2, true);
    new G4PVPlacement(Rot_135, G4ThreeVector(0., 0., -quartz_windowHolderTop_z), logic_windowHolderTop, "QuartzWindowHolderTopD", logic_vessel, false, 3, true);

    /// Build PMT
    pmt_.Construct();
    G4LogicalVolume* logic_pmt = pmt_.GetLogicalVolume();
    G4double pmt_length_ = pmt_.Length();
    G4double pmt_z  = 42.495*mm + pmt_length_/2;

    G4Tubs *solid_pmt = new G4Tubs("SolidPMT", 0., pmt_rad_, pmt_length_/2, 0., 360.*deg); // Hamamatsu pmt length: 43*mm | STEP pmt gap length: 57.5*mm

    // Position pairs (x,Y) for PMTs
    std::vector <float> pmt_PsX={-15.573, 20.68, -36.253, 0., 36.253, -20.68, 15.573};
    std::vector <float> pmt_PsY={-32.871, -29.922, -2.949, 0., 2.949, 29.922, 32.871};

    // PMT clad
    G4VSolid *solid_enclosure_pmt = new G4Tubs("EnclosurePMT", 0, enclosure_pmt_rad_+ enclosure_pmt_thickn_ , (enclosure_pmt_length_)/2, 0., 360.*deg);
    G4double enclosure_pmt_z = vessel_length_/2 - enclosure_pmt_length_/2;
    G4double relative_pmt_z  = enclosure_pmt_z - pmt_z;
    // Vacuum inside the pmt enclosure
    G4double enclosurevac_pmt_z = vessel_length_/2 - enclosurevac_pmt_length_/2;
    G4double relativevac_pmt_z  = enclosurevac_pmt_z - pmt_z;
    G4VSolid *solid_enclosurevac_pmt = new G4Tubs("EnclosureVacPMT", 0, enclosure_pmt_rad_, (enclosurevac_pmt_length_)/2, 0., 360.*deg);
    // PMT Holder
    G4double pmtHolder_z = enclosurevac_pmt_length_/2 - pmtHolder_length_/2;
    G4VSolid *solid_pmtHolder = new G4Tubs("PMTHolder", 0, pmtHolder_rad_, (pmtHolder_length_)/2, 0., 360.*deg);
    // Steel plate enclosing the pmt tube
    G4double plate1_pmt_z = enclosure_pmt_z - enclosure_pmt_length_/2 - plate_pmt_length_/2;
    G4VSolid *solid_plate1_pmt = new G4Tubs("PMTplateBottom1", 0, plate_pmt_rad_+plate_pmt_thickn_, plate_pmt_length_/2, 0, 360*deg);

    G4int nHole = 7;
    G4ThreeVector pos_pmt = G4ThreeVector(0, 0, 0);
    G4ThreeVector pos_enclosure_pmt = G4ThreeVector(0, 0, 0);
    G4ThreeVector pos_enclosurevac_pmt = G4ThreeVector(0, 0, 0);
    G4ThreeVector pos = G4ThreeVector(0, 0, 0);

    for (G4int i = 0; i < nHole; i++) {
      pos_pmt = G4ThreeVector(pmt_PsX[i]*mm, pmt_PsY[i]*mm, -pmt_z);
      pos_enclosure_pmt = G4ThreeVector(pmt_PsX[i]*mm, pmt_PsY[i]*mm, relative_pmt_z);
      pos_enclosurevac_pmt = G4ThreeVector(pmt_PsX[i]*mm, pmt_PsY[i]*mm, relativevac_pmt_z);
      pos = G4ThreeVector(pmt_PsX[i]*mm, pmt_PsY[i]*mm, 0.);

      solid_enclosure_pmt = new G4SubtractionSolid("EnclosurePMT_Sub", solid_enclosure_pmt, solid_pmt,  0, pos_enclosure_pmt);
      solid_enclosurevac_pmt = new G4SubtractionSolid("EnclosureVacPMT_Sub", solid_enclosurevac_pmt, solid_pmt,  0, pos_enclosurevac_pmt);
      solid_pmtHolder = new G4SubtractionSolid("PMTHolder_Sub", solid_pmtHolder, solid_pmt,  0, pos);
      solid_plate1_pmt = new G4SubtractionSolid("PMTplateBottom1_Sub", solid_plate1_pmt, solid_pmt,  0, pos);

      new G4PVPlacement(0, pos_pmt, logic_pmt, "PMT", logic_vessel, false, i, true);
    }

    G4LogicalVolume *logic_enclosure_pmt = new G4LogicalVolume(solid_enclosure_pmt, steel, "EnclosurePMT");
    new G4PVPlacement(0, G4ThreeVector(0., 0., -enclosure_pmt_z), logic_enclosure_pmt, "EnclosurePMT", logic_vessel, false, 0, true);

    G4LogicalVolume *logic_enclosurevac_pmt = new G4LogicalVolume(solid_enclosurevac_pmt, vacuum, "EnclosureVacPMT");
    new G4PVPlacement(0, G4ThreeVector(0., 0, enclosure_pmt_z - enclosurevac_pmt_z), logic_enclosurevac_pmt, "EnclosureVacPMT", logic_enclosure_pmt, false, 0, true);

    G4LogicalVolume *logic_pmtHolder = new G4LogicalVolume(solid_pmtHolder, steel, "PMTHolder");
    new G4PVPlacement(0, G4ThreeVector(0., 0, pmtHolder_z), logic_pmtHolder, "PMTHolder", logic_enclosurevac_pmt, false, 0, true);

    G4LogicalVolume *logic_plate1_pmt = new G4LogicalVolume(solid_plate1_pmt, steel, "PMTplateBottom1");
    new G4PVPlacement(0, G4ThreeVector(0., 0., -plate1_pmt_z), logic_plate1_pmt, "PMTplateBottom1", logic_vessel, false, 0, true);

    // Steel plate attached where the peek holders are attached
    G4Tubs          *solid_plate0_pmt = new G4Tubs("PMTplateBottom0", plate_pmt_rad_, plate_pmt_rad_+plate_pmt_thickn_, plate_pmt_length_/2, 0, 360*deg);
    G4LogicalVolume *logic_plate0_pmt = new G4LogicalVolume(solid_plate0_pmt, steel, "PMTplateBottom0");

    G4double plate0_pmt_z = plate1_pmt_z - plate_pmt_length_;
    new G4PVPlacement(0, G4ThreeVector(0., 0., -plate0_pmt_z), logic_plate0_pmt, "PMTplateBottom0", logic_vessel, false, 0, true);

    // Upper steel plate at the pmt clad
    G4double plateUp_pmt_rad_ = enclosure_pmt_rad_ + enclosure_pmt_thickn_;

    G4Tubs          *solid_plateUp_pmt = new G4Tubs("PMTplateUp", plateUp_pmt_rad_, plateUp_pmt_rad_+plateUp_pmt_thickn_, plateUp_pmt_length_/2, 0, 360*deg);
    G4LogicalVolume *logic_plateUp_pmt = new G4LogicalVolume(solid_plateUp_pmt, steel, "PMTplateUp");

    G4double plateUp_pmt_z = vessel_length_/2 - plateUp_pmt_length_/2 ;
    new G4PVPlacement(0, G4ThreeVector(0., 0., -plateUp_pmt_z), logic_plateUp_pmt, "PMTplateUp", logic_vessel, false, 0, true);

}
