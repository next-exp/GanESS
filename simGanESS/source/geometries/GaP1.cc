#include "GaP1.hh"

#include <G4Tubs.hh>
#include "G4LogicalVolume.hh"
#include "G4NistManager.hh"
#include "G4PVPlacement.hh"
#include <G4UserLimits.hh>
#include <G4SDManager.hh>

#include "nexus/FactoryBase.h"
#include "nexus/UniformElectricDriftField.h"
#include "nexus/OpticalMaterialProperties.h"
#include "nexus/XenonProperties.h"
#include "nexus/IonizationSD.h"


using namespace nexus;

REGISTER_CLASS(GaP1, GeometryBase)

GaP1::GaP1(): 
    GeometryBase(),
    msg_ (nullptr),
    gas_ (nullptr),
    mesh_mat_ (nullptr),
    vessel_rad_        (104./2  *mm),
    vessel_thickn_     (0.5 *cm),
    vessel_length_     (10.  *cm),
    mesh_thickn_       (0.01  *mm),
    mesh_transparency_ (0.95),
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
  msg_ = new G4GenericMessenger(this, "/Geometry/GaP1/",
                                "Control commands of the GaP1 geometry.");
  // Parametrized dimensions
  DefineConfigurationParameters();
}

GaP1::~GaP1()
{
  delete msg_;

  delete buffer_gen_;
  delete drift_gen_;
  delete gas_pmt_gen_;
  delete el_gen_;

}

void GaP1::Construct()
{
    //Materials
    //auto nistManager = G4NistManager::Instance(); // Nist manager to retrieve materials.
    //auto air   = nistManager->FindOrBuildMaterial("G4_AIR");
    gas_   = materials::GXe(pressure_, temperature_);
    gas_->SetMaterialPropertiesTable(opticalprops::GXe(pressure_,
                                                      temperature_,
                                                      sc_yield_,
                                                      elifetime_));

    auto steel = materials::Steel();
    steel->SetMaterialPropertiesTable(new G4MaterialPropertiesTable());
    // Meshes materials
    G4double photoe_prob_ = 0.;
    mesh_mat_ = materials::FakeDielectric(gas_, "mesh_mat");
    mesh_mat_->SetMaterialPropertiesTable(opticalprops::FakeGrid(pressure_,
              temperature_, mesh_transparency_, mesh_thickn_,
              sc_yield_, elifetime_, photoe_prob_));

    //Cylinder, acting as the vessel
    G4Tubs          *solid_vessel = new G4Tubs("Vessel", 0., vessel_rad_ + vessel_thickn_, vessel_length_/2 + vessel_thickn_, 0., 360.*deg);
    G4LogicalVolume *logic_vessel = new G4LogicalVolume(solid_vessel, steel, "Vessel");
    this->SetLogicalVolume(logic_vessel);

    //Build inside detector
    BuildTPC(gas_, mesh_mat_, logic_vessel);
}

G4ThreeVector GaP1::GenerateVertex(const G4String& region) const
{
  G4ThreeVector vertex;

  if (region == "AD_HOC") {
    vertex = specific_vertex_;
  }

  //// Gas regions
  else if (
    (region == "GasPMT") ||
    (region == "GasEL") ||
    (region == "GasDrift") ||
    (region == "GasBuffer")) {
    vertex = GenerateVertexGas(region);
  }

  else {
    G4Exception("[GaP1]", "GenerateVertex()", FatalException,
      "Unknown vertex generation region!");
  }

  return vertex;
}

G4ThreeVector GaP1::GenerateVertexGas(const G4String& region) const
{
    G4ThreeVector vertex;

    if     (region == "GasEL")     {vertex = el_gen_->GenerateVertex("VOLUME");}
    else if(region == "GasBuffer") {vertex = buffer_gen_->GenerateVertex("VOLUME");}
    else if(region == "GasDrift")  {vertex = drift_gen_->GenerateVertex("VOLUME");}
    else if(region == "GasPMT")    {vertex = gas_pmt_gen_->GenerateVertex("VOLUME");}
    else{G4Exception("[GaP1]", "GenerateVertex()", FatalException,
                "Unknown vertex generation region!");}
    return vertex;
}


void GaP1::DefineConfigurationParameters()
{
  // Gas pressure
  G4GenericMessenger::Command& pressure_cmd =
    msg_->DeclareProperty("pressure", pressure_,
                          "Pressure of the gas.");
  pressure_cmd.SetUnitCategory("Pressure");
  pressure_cmd.SetParameterName("pressure", false);
  pressure_cmd.SetRange("pressure>0.");

  // Specific vertex in case region to shoot from is AD_HOC
  msg_->DeclarePropertyWithUnit("specific_vertex", "mm",  specific_vertex_, "Set generation vertex.");
}

void GaP1::BuildTPC(G4Material* gas, G4Material* mesh_mat, G4LogicalVolume* logic_vessel)
{
    //Gas
    G4double drift_length_  = 19.825*mm;
    G4double el_length_     = 10.775*mm;
    G4double pmt_gap_       = 23.2  *mm;
    G4double buffer_        = 5     *mm; // Cap to cathode buffer

    //// Buffer
    G4Tubs *solid_gas_buffer = new G4Tubs("GasBuffer", 0., vessel_rad_, (buffer_)/2, 0., 360.*deg);
    G4LogicalVolume *logic_gas_buffer = new G4LogicalVolume(solid_gas_buffer, gas, "GasBuffer");

    G4double buffer_z = (vessel_length_/2) - buffer_/2;
    G4VPhysicalVolume* buffer_phys_ = new G4PVPlacement(0, G4ThreeVector(0., 0, buffer_z), logic_gas_buffer, "GasBuffer", logic_vessel, false, 0, true);
    buffer_gen_  = new CylinderPointSampler2020(buffer_phys_);

    //// Cathode 
    G4Tubs *solid_cathode = new G4Tubs("Cathode", 0., vessel_rad_, (mesh_thickn_)/2, 0., 360.*deg);
    G4LogicalVolume *logic_cathode = new G4LogicalVolume(solid_cathode, mesh_mat, "Cathode");

    G4double cathode_z = - buffer_/2 + mesh_thickn_/2;
    new G4PVPlacement(0, G4ThreeVector(0., 0, cathode_z), logic_cathode, "Cathode", logic_gas_buffer, false, 0, true);

    //// Drift
    G4Tubs *solid_gas_drift = new G4Tubs("GasDrift", 0., vessel_rad_, (drift_length_)/2, 0., 360.*deg);
    G4LogicalVolume *logic_gas_drift = new G4LogicalVolume(solid_gas_drift, gas, "GasDrift");

    G4double drift_z = buffer_z - buffer_/2 - drift_length_/2;
    G4VPhysicalVolume* drift_phys_ = new G4PVPlacement(0, G4ThreeVector(0., 0, drift_z), logic_gas_drift, "GasDrift", logic_vessel, false, 0, true);
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
    IonizationSD* active_sd = new IonizationSD("/GaP1/DRIFT");
    logic_gas_drift->SetSensitiveDetector(active_sd);
    G4SDManager::GetSDMpointer()->AddNewDetector(active_sd);

    //// EL gap
    G4Tubs *solid_gas_el = new G4Tubs("GasEL", 0., vessel_rad_, (el_length_)/2, 0., 360.*deg);
    G4LogicalVolume *logic_gas_el = new G4LogicalVolume(solid_gas_el, gas, "GasEL");

    G4double el_z = drift_z - drift_length_/2 - el_length_/2;
    G4VPhysicalVolume* el_phys_ = new G4PVPlacement(0, G4ThreeVector(0., 0, el_z), logic_gas_el, "GasEL", logic_vessel, false, 0, true);
    el_gen_  = new CylinderPointSampler2020(el_phys_);

    //// Gate 
    G4Tubs *solid_gate = new G4Tubs("Gate", 0., vessel_rad_, (mesh_thickn_)/2, 0., 360.*deg);
    G4LogicalVolume *logic_gate = new G4LogicalVolume(solid_gate, mesh_mat, "Gate");

    G4double gate_z = el_length_/2 - mesh_thickn_/2;
    new G4PVPlacement(0, G4ThreeVector(0., 0, gate_z), logic_gate, "Gate", logic_gas_el, false, 0, true);

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
    G4Tubs *solid_anode = new G4Tubs("Anode", 0., vessel_rad_, (mesh_thickn_)/2, 0., 360.*deg);
    G4LogicalVolume *logic_anode = new G4LogicalVolume(solid_anode, mesh_mat, "Anode");

    G4double anode_z = -el_length_/2 + mesh_thickn_/2;
    new G4PVPlacement(0, G4ThreeVector(0., 0, anode_z), logic_anode, "Anode", logic_gas_el, false, 0, true);

    G4cout << "* GATE Z position: " << el_z + gate_z << G4endl;
    G4cout << "* GATE Volt position: " << el_z + el_length_/2. << G4endl;
    G4cout << "* ANODE Z position: " << el_z + anode_z << G4endl;
    G4cout << "* ANODE Volt position: " << el_z - el_length_/2. << G4endl;
    G4cout << "* EL_GAP Z positions: " << el_z - el_length_/2. <<
              " to " << el_z + el_length_/2. << G4endl;

    G4cout << "* EL field intensity (kV/cm): " << el_field_ / (kilovolt/cm)
            << "  ->  EL Light yield (photons/ie-/cm): " << yield / (1/cm) << G4endl;
    G4cout << "* EL Light yield (photons/ie-): " << yield * el_length_ << G4endl;

    //// PMT gap
    G4Tubs *solid_gas_pmt = new G4Tubs("GasPMT", 0., vessel_rad_, (pmt_gap_)/2, 0., 360.*deg);
    G4LogicalVolume *logic_gas_pmt = new G4LogicalVolume(solid_gas_pmt, gas, "GasPMT");

    G4double pmt_gap_z = el_z - el_length_/2 - pmt_gap_/2;
    G4VPhysicalVolume* gas_pmt_phys_ = new G4PVPlacement(0, G4ThreeVector(0., 0, pmt_gap_z), logic_gas_pmt, "GasPMT", logic_vessel, false, 0, true);   
    gas_pmt_gen_  = new CylinderPointSampler2020(gas_pmt_phys_);

    //Build PMT
    pmt_.Construct();
    G4LogicalVolume* logic_pmt = pmt_.GetLogicalVolume();
    G4double pmt_length_ = pmt_.Length();
    G4double pmt_z       = pmt_gap_z - pmt_gap_/2 - pmt_length_/2;
    new G4PVPlacement(0, G4ThreeVector(0.,0., pmt_z),
		              logic_pmt, "PMT",
		              logic_vessel, false, 0, true);
}