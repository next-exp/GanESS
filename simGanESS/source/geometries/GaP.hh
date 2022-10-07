#ifndef GaP_H
#define GaP_H

#include "nexus/GeometryBase.h"
#include "nexus/CylinderPointSampler2020.h"
#include "nexus/MaterialsList.h"

#include <G4ThreeVector.hh>
#include <G4GenericMessenger.hh>
#include "nexus/PmtR7378A.h"

using namespace nexus;

class GaP : public GeometryBase
{

    public:
    // Constructor
        GaP();
    //Destructor
        ~GaP();

          G4ThreeVector GenerateVertex   (const G4String& region) const;
          G4ThreeVector GenerateVertexGas(const G4String& region) const;

    private:
        void DefineConfigurationParameters();
        void Construct();
        void BuildTPC(G4Material* gas, G4Material* mesh_mat, G4LogicalVolume* logic_vessel);

    private:
        G4GenericMessenger* msg_;

        // Materials
        G4Material* gas_;
        G4Material* mesh_mat_;

        // Vessel parameters
        G4double vessel_out_rad_   ;
        G4double vessel_out_length_;
        G4double vessel_rad_   ;
        G4double vessel_length_;

        // Meshes
        G4double mesh_rad_   ;
        G4double mesh_thickn_;
        G4double mesh_transparency_;

        G4double photoe_prob_;

        // Gas parameters
        G4double pressure_;
        G4double temperature_;
        G4double sc_yield_;
        G4double elifetime_;
        G4double drift_vel_;
        G4double drift_transv_diff_;
        G4double drift_long_diff_;

        // EL field (in kV/cm)
        G4double el_field_;
        G4double el_vel_;
        G4double el_transv_diff_;
        G4double el_long_diff_;

        // PMTs
        PmtR7378A pmt_;

        // Vertex generation
        G4ThreeVector specific_vertex_;

        CylinderPointSampler2020* buffer_gen_;
        CylinderPointSampler2020* drift_gen_;
        CylinderPointSampler2020* el_gen_;
        CylinderPointSampler2020* gas_pmt_gen_;

};

#endif
