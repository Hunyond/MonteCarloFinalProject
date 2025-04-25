#pragma once
#include <map>
#include <string>


struct InteractionData {
    float incoherentScatter; // Incoherent scattering cross-section (cm²/g)
    float photoelAbsorb;    // Photoelectric absorption cross-section (cm²/g)
    float pairProd;    // Nuclear pair production cross-section (cm²/g)
};


InteractionData getCrossSectionsAtEnergy (const std::map<float, InteractionData>& dataMap, const float targetEnergy);
std::map<float, InteractionData> loadPhotonDataToMap(const std::string& filename, const float density);

