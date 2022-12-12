// ----------------------------------------------------------------------------
// GanESS | OpticalMaterialProperties_GanESS.cc
//
// Optical properties of relevant materials.
//
// ----------------------------------------------------------------------------

#include "nexus/OpticalMaterialProperties.h"

#include "OpticalMaterialProperties_GanESS.h"

#include <G4MaterialPropertiesTable.hh>

#include <assert.h>

using namespace CLHEP;

namespace opticalprops_GanESS{

  /// Steel (== Stainless Steel 340L) ///
  G4MaterialPropertiesTable* Steel()
  {
    G4MaterialPropertiesTable* mpt = new G4MaterialPropertiesTable();

    // REFLECTIVITY between 100-200[nm]
    std::vector<G4double> ENERGIES = {
      optPhotMinE_GanESS_,
      h_Planck * c_light / (191. * nm), h_Planck * c_light / (181. * nm),
      h_Planck * c_light / (171. * nm), h_Planck * c_light / (161. * nm),
      h_Planck * c_light / (151. * nm), h_Planck * c_light / (141. * nm),
      h_Planck * c_light / (131. * nm), h_Planck * c_light / (121. * nm),
      optPhotMaxE_GanESS_
    };
    std::vector<G4double> REFLECTIVITY = {
      .46,
      .43, .40,
      .37, .34,
      .30, .27,
      .24, .20,
      .19
    };
    mpt->AddProperty("REFLECTIVITY", ENERGIES, REFLECTIVITY);

    // REFLEXION BEHAVIOR
    std::vector<G4double> ENERGIES_2    = {optPhotMinE_GanESS_, optPhotMaxE_GanESS_};
    // Specular reflection about the normal to a microfacet.
    // Such a vector is chosen according to a gaussian distribution with
    // sigma = SigmaAlhpa (in rad) and centered in the average normal.
    std::vector<G4double> specularlobe  = {.19, .07};
    // specular reflection about the average normal
    std::vector<G4double> specularspike = {0., 0.};
    // 180 degrees reflection.
    std::vector<G4double> backscatter   = {0., 0.};
    // 1 - the sum of these three last parameters is the percentage of Lambertian reflection

    mpt->AddProperty("SPECULARLOBECONSTANT", ENERGIES_2, specularlobe);
    mpt->AddProperty("SPECULARSPIKECONSTANT",ENERGIES_2, specularspike);
    mpt->AddProperty("BACKSCATTERCONSTANT",  ENERGIES_2, backscatter);

    return mpt;
  }

}
