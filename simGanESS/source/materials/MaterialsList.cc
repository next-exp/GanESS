#include "MaterialsList.hh"

#include <G4Material.hh>
#include <G4Element.hh>
#include <G4NistManager.hh>

namespace materials {

G4Material* GAr(G4double pressure, G4double temperature)
{
  G4String name = "GAr";
  G4Material* mat = G4Material::GetMaterial(name, false);
  
  if (mat == 0) {
    G4double density = 1.60279 * (pressure / bar) * (300 * kelvin) / temperature * kg/m3; // 1.60729 is the pressure at 1 bar and 300K
    G4NistManager* nist = G4NistManager::Instance();

    mat = new G4Material(name, density, 1, kStateGas, temperature, pressure);

    G4Element* Ar = nist->FindOrBuildElement("Ar");

    mat->AddElement(Ar,1);
  }

  return mat;
}

G4Material* GKr(G4double pressure, G4double temperature)
{
  G4String name = "GKr";
  G4Material* mat = G4Material::GetMaterial(name, false);

  if (mat == 0) {
    G4double density = 3.749 * (pressure / bar) * (273.15 * kelvin) / temperature * kg/m3; // 3.749 is the pressure at 1 bar and 273.15K
    G4NistManager* nist = G4NistManager::Instance();

    mat = new G4Material(name, density, 1, kStateGas, temperature, pressure);

    G4Element* Kr = nist->FindOrBuildElement("Kr");

    mat->AddElement(Kr,1);
  }

  return mat;
}
}