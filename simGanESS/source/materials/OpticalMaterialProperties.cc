
#include "OpticalMaterialProperties.hh"
#include <G4MaterialPropertiesTable.hh>

namespace opticalprops {
  /// Gaseous Argon ///
  G4MaterialPropertiesTable* GAr(G4double sc_yield,
                                G4double e_lifetime)
  {
    // An argon gas proportional scintillation counter with UV avalanche photodiode scintillation
    // readout C.M.B. Monteiro, J.A.M. Lopes, P.C.P.S. Simoes, J.M.F. dos Santos, C.A.N. Conde
    G4MaterialPropertiesTable* mpt = new G4MaterialPropertiesTable();

    // REFRACTIVE INDEX
    const G4int ri_entries = 200;
    G4double eWidth = (optPhotMaxE_ - optPhotMinE_) / ri_entries;

    std::vector<G4double> ri_energy;
    for (int i=0; i<ri_entries; i++) {
      ri_energy.push_back(optPhotMinE_ + i * eWidth);
    }

    std::vector<G4double> rIndex;
    for (int i=0; i<ri_entries; i++) {
      G4double wl = h_Planck * c_light / ri_energy[i] * 1000; // in micron
      // From refractiveindex.info
      rIndex.push_back(1 + 0.012055*(0.2075*pow(wl,2)/(91.012*pow(wl,2)-1) +
                                     0.0415*pow(wl,2)/(87.892*pow(wl,2)-1) +
                                     4.3330*pow(wl,2)/(214.02*pow(wl,2)-1)));
      //G4cout << "* GAr rIndex:  " << std::setw(5) << ri_energy[i]/eV
      //       << " eV -> " << rIndex[i] << G4endl;
    }
    mpt->AddProperty("RINDEX", ri_energy, rIndex);

    // ABSORPTION LENGTH
    std::vector<G4double> abs_energy = {optPhotMinE_, optPhotMaxE_};
    std::vector<G4double> absLength  = {noAbsLength_, noAbsLength_};
    mpt->AddProperty("ABSLENGTH", abs_energy, absLength);

    // EMISSION SPECTRUM
//    G4double Wavelength_peak  = 128.000 * nm;
    G4double Wavelength_peak  = 128.000 * nm; // Xe, to be changed back
    G4double Wavelength_sigma =   2.929 * nm;
    G4double Energy_peak  = (h_Planck*c_light / Wavelength_peak);
    G4double Energy_sigma = (h_Planck*c_light * Wavelength_sigma / pow(Wavelength_peak,2));
    //G4cout << "*** GAr Energy_peak: " << Energy_peak/eV << " eV   Energy_sigma: "
    //       << Energy_sigma/eV << " eV" << G4endl;

    // Sampling from ~110 nm to 150 nm <----> from ~11.236 eV to 8.240 eV
    const G4int sc_entries = 380;
    std::vector<G4double> sc_energy;
    std::vector<G4double> intensity;
    for (int i=0; i<sc_entries; i++){
      sc_energy.push_back(8.240*eV + 0.008*i*eV);
      intensity.push_back(exp(-pow(Energy_peak/eV-sc_energy[i]/eV,2) /
                              (2*pow(Energy_sigma/eV, 2)))/(Energy_sigma/eV*sqrt(pi*2.)));
      //G4cout << "* GAr energy: " << std::setw(6) << sc_energy[i]/eV << " eV  ->  "
      //       << std::setw(6) << intensity[i] << G4endl;
    }
    mpt->AddProperty("SCINTILLATIONCOMPONENT1", sc_energy, intensity);
    mpt->AddProperty("SCINTILLATIONCOMPONENT2", sc_energy, intensity);
    mpt->AddProperty("ELSPECTRUM"             , sc_energy, intensity, 1);

    // CONST PROPERTIES
    mpt->AddConstProperty("SCINTILLATIONYIELD", sc_yield);
    mpt->AddConstProperty("SCINTILLATIONTIMECONSTANT1",   6.*ns); // From https://dspace.mit.edu/bitstream/handle/1721.1/129347/1903.06706.pdf?sequence=2&isAllowed=y
    mpt->AddConstProperty("SCINTILLATIONTIMECONSTANT2",   3480.*ns); // From https://dspace.mit.edu/bitstream/handle/1721.1/129347/1903.06706.pdf?sequence=2&isAllowed=y
    mpt->AddConstProperty("SCINTILLATIONYIELD1", .136); // From https://dspace.mit.edu/bitstream/handle/1721.1/129347/1903.06706.pdf?sequence=2&isAllowed=y
    mpt->AddConstProperty("SCINTILLATIONYIELD2", .864); // From https://dspace.mit.edu/bitstream/handle/1721.1/129347/1903.06706.pdf?sequence=2&isAllowed=y
    
    mpt->AddConstProperty("ELTIMECONSTANT1",   6.*ns, 1); // From https://dspace.mit.edu/bitstream/handle/1721.1/129347/1903.06706.pdf?sequence=2&isAllowed=y
    mpt->AddConstProperty("ELTIMECONSTANT2",   3480.*ns, 1); // From https://dspace.mit.edu/bitstream/handle/1721.1/129347/1903.06706.pdf?sequence=2&isAllowed=y
    mpt->AddConstProperty("ELTIMECONSTANT3",  0, 1); // From https://dspace.mit.edu/bitstream/handle/1721.1/129347/1903.06706.pdf?sequence=2&isAllowed=y
    mpt->AddConstProperty("ELYIELD1", .136, 1); // From https://dspace.mit.edu/bitstream/handle/1721.1/129347/1903.06706.pdf?sequence=2&isAllowed=y
    mpt->AddConstProperty("ELYIELD2", .864, 1); // From https://dspace.mit.edu/bitstream/handle/1721.1/129347/1903.06706.pdf?sequence=2&isAllowed=y   
    mpt->AddConstProperty("ELYIELD3", 0, 1); // From https://dspace.mit.edu/bitstream/handle/1721.1/129347/1903.06706.pdf?sequence=2&isAllowed=y   
    mpt->AddConstProperty("RESOLUTIONSCALE",    1.0);
    mpt->AddConstProperty("ATTACHMENT",         e_lifetime, 1);

    return mpt;
  }
}