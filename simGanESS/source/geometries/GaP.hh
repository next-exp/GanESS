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
        void BuildTPC(G4Material* gas, G4Material* mesh_mat, G4Material* steel, G4Material* peek, G4Material* vacuum, G4Material* quartz, G4Material* tpb, G4LogicalVolume* logic_vessel);

    private:
        G4GenericMessenger* msg_;

        // Materials
        G4Material* gas_;
        G4Material* mesh_mat_;
        G4Material* steel_;
        G4Material* peek_;
        G4Material* vacuum_;
        G4Material* quartz_;
        G4Material* tpb_;

        // Gas being used.
        G4String gas_element_;

        // Vessel parameters
        G4double vessel_out_rad_   ;
        G4double vessel_out_length_;
        G4double vessel_rad_   ;
        G4double vessel_length_;

        // Mesh
        G4double mesh_rad_   ;
        G4double mesh_thickn_;
        G4double mesh_transparency_;

        // Mesh Bracket
        G4double meshBracket_rad_;
        G4double meshBracket_thickn_ ;
        G4double anodeBracket_rad_ ;
        G4double anodeBracket_thickn_ ;

        G4double pmt_rad_;

        // Pmt enclosure
        G4double enclosure_pmt_rad_   ;
        G4double enclosure_pmt_thickn_;
        G4double enclosure_pmt_length_;
        G4double enclosurevac_pmt_length_;

        G4double plate_pmt_rad_;
        G4double plate_pmt_thickn_;
        G4double plate_pmt_length_;
        G4double plateUp_pmt_length_;
        G4double plateUp_pmt_thickn_;

        G4double pmtHolder_rad_;
        G4double pmtHolder_length_;

        //Quartz window with tpb
        G4double quartz_window_rad_;
        G4double quartz_window_thickn_;
        G4double tpb_coating_thickn_;

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
        G4double yield_;
        G4double el_vel_;
        G4double el_transv_diff_;
        G4double el_long_diff_;

        // PMTs
        PmtR7378A pmt_;

        // Vertex generation
        G4ThreeVector specific_vertex_;

        CylinderPointSampler2020* drift_gen_;
        CylinderPointSampler2020* el_gen_;


};

#endif
