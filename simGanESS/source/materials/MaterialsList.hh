#ifndef MATERIALS_LIST_H
#define MATERIALS_LIST_H

#include <globals.hh>

#include <CLHEP/Units/SystemOfUnits.h>
#include <CLHEP/Units/PhysicalConstants.h>

class G4Material;

namespace materials {
    using namespace CLHEP;
    // Argon
    G4Material* GAr(G4double pressure=STP_Pressure,
			        G4double temperature=STP_Temperature);
    // Krypton
    G4Material* GKr(G4double pressure=STP_Pressure,
			        G4double temperature=STP_Temperature);
}

#endif