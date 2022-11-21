// ----------------------------------------------------------------------------
// GanESS | OpticalMaterialProperties_GanESS.h
//
// Optical properties of relevant materials.
//
// ----------------------------------------------------------------------------

#ifndef OPTICAL_MATERIAL_PROPERTIES_GANESS_H
#define OPTICAL_MATERIAL_PROPERTIES_GANESS_H

#include <globals.hh>

#include <CLHEP/Units/SystemOfUnits.h>
#include <CLHEP/Units/PhysicalConstants.h>

//using namespace nexus;

class G4MaterialPropertiesTable;

namespace opticalprops_GanESS {

  using namespace CLHEP;

  G4MaterialPropertiesTable* Steel();

  constexpr G4double optPhotMinE_GanESS_ =  6.2 * eV;
  constexpr G4double optPhotMaxE_GanESS_ =  11.3 * eV;

  // Constant that allows to convert nm to eV:
  // nm_to_eV_ / wavelength (nm) = energy (eV)
  constexpr G4double nm_to_eV_ = h_Planck * c_light * 1.e6;


} // end namespace opticalprops_GanESS

#endif
