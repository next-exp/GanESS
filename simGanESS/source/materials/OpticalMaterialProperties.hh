#ifndef OPTICAL_MATERIAL_PROPERTIES_H
#define OPTICAL_MATERIAL_PROPERTIES_H

#include <globals.hh>

#include <CLHEP/Units/SystemOfUnits.h>
#include <CLHEP/Units/PhysicalConstants.h>

class G4MaterialPropertiesTable;


namespace opticalprops {

  using namespace CLHEP;
  G4MaterialPropertiesTable* GAr(G4double sc_yield,
                                G4double e_lifetime=1000.*ms);

  constexpr G4double optPhotMinE_ =  0.2  * eV;
  constexpr G4double optPhotMaxE_ = 11.5  * eV;
  constexpr G4double noAbsLength_ = 1.e8  * m;

  // Constant that allows to convert nm to eV:
  // nm_to_eV_ / wavelength (nm) = energy (eV)
  constexpr G4double nm_to_eV_ = h_Planck * c_light * 1.e6;
}

#endif