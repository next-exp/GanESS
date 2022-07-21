#include "GaP1.hh"

#include <G4Tubs.hh>
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


using namespace nexus;

REGISTER_CLASS(GaP1, GeometryBase)

GaP1::GaP1(): 
    GeometryBase(),
    msg_ (nullptr),
    gas_ (nullptr),
    cath_mat_ (nullptr),
    gate_mat_ (nullptr),
    anode_mat_(nullptr),
    vessel_rad_        (104./2  *mm),
    vessel_thickn_     (0.5 *cm),
    vessel_length_     (10.  *cm),
    cath_thickn_       (0.01  *mm),
    gate_thickn_       (0.01  *mm),
    anode_thickn_      (0.01  *mm),
    cath_transparency_ (0.95),
    gate_transparency_ (0.95),
    anode_transparency_(0.95),
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
    // Cathode material
    cath_mat_ = materials::FakeDielectric(gas_, "cath_mat");
    cath_mat_->SetMaterialPropertiesTable(opticalprops::FakeGrid(pressure_,
              temperature_, cath_transparency_, cath_thickn_,
              sc_yield_, elifetime_, photoe_prob_));
    // Gate material
    gate_mat_ = materials::FakeDielectric(gas_, "gate_mat");
    gate_mat_->SetMaterialPropertiesTable(opticalprops::FakeGrid(pressure_,
              temperature_, gate_transparency_, gate_thickn_,
              sc_yield_, elifetime_, photoe_prob_));
    // Anode material
    anode_mat_ = materials::FakeDielectric(gas_, "anode_mat");
    anode_mat_->SetMaterialPropertiesTable(opticalprops::FakeGrid(pressure_,
              temperature_, anode_transparency_, anode_thickn_,
              sc_yield_, elifetime_, photoe_prob_));


    //Cylinder, acting as the vessel
    G4Tubs          *solid_vessel = new G4Tubs("Vessel", 0., vessel_rad_ + vessel_thickn_, vessel_length_/2 + vessel_thickn_, 0., 360.*deg);
    G4LogicalVolume *logic_vessel = new G4LogicalVolume(solid_vessel, steel, "Vessel");
    this->SetLogicalVolume(logic_vessel);

    //Build inside detector
    BuildTPC(gas_, cath_mat_, gate_mat_, anode_mat_, logic_vessel);
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

void GaP1::BuildTPC(G4Material* gas, G4Material* cath_mat, G4Material* gate_mat, G4Material* anode_mat, G4LogicalVolume* logic_vessel)
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
    G4Tubs *solid_cathode = new G4Tubs("Cathode", 0., vessel_rad_, (cath_thickn_)/2, 0., 360.*deg);
    G4LogicalVolume *logic_cathode = new G4LogicalVolume(solid_cathode, cath_mat, "Cathode");

    G4double cathode_z = - buffer_/2 + cath_thickn_/2;
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
    G4Tubs *solid_gate = new G4Tubs("Gate", 0., vessel_rad_, (gate_thickn_)/2, 0., 360.*deg);
    G4LogicalVolume *logic_gate = new G4LogicalVolume(solid_gate, gate_mat, "Gate");

    G4double gate_z = el_length_/2 - gate_thickn_/2;
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
    G4Tubs *solid_anode = new G4Tubs("Anode", 0., vessel_rad_, (anode_thickn_)/2, 0., 360.*deg);
    G4LogicalVolume *logic_anode = new G4LogicalVolume(solid_anode, anode_mat, "Anode");

    G4double anode_z = -el_length_/2 + anode_thickn_/2;
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