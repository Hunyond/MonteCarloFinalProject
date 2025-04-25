#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <cctype> // for std::isspace
#include <iterator> // for std::find_if
#include <utility> // for std::pair
#include <vector>
#include <iterator>
#include "crosssections.hpp"






InteractionData getCrossSectionsAtEnergy (const std::map<float, InteractionData>& dataMap, const float targetEnergy)
{
    if (dataMap.empty ()) {
        throw std::out_of_range ("Empty data map");
    }
    
    // Check for exact match first
    auto exact = dataMap.find (targetEnergy);
    if (exact != dataMap.end ()) {
        return exact->second;
    }
    
    // Check if below minimum energy
    if (targetEnergy < dataMap.begin ()->first) {
        std::cout << "Warning: Energy out of range, treating particle as absorbed!!!" << std::endl;
        return  {0.0f, 1.0f, 0.0f}; 
    }
    
    // Check if above maximum energy
    if (targetEnergy > dataMap.rbegin ()->first) {
        std::cout << "Warning: Energy out of range, returning maximum value!!!" << std::endl;
        return dataMap.rbegin ()->second;
    }
    
    // Normal interpolation case
    auto upper = dataMap.lower_bound (targetEnergy);
    auto lower = std::prev (upper);
    const float x0 = lower->first;
    const float x1 = upper->first;
    const InteractionData& y0 = lower->second;
    const InteractionData& y1 = upper->second;
    const float delta = x1 - x0;
    const float t = (targetEnergy - x0) / delta;

    return {
        std::lerp (y0.incoherentScatter, y1.incoherentScatter, t),
        std::lerp (y0.photoelAbsorb, y1.photoelAbsorb, t),
        std::lerp (y0.pairProd, y1.pairProd, t)};
}



std::map<float, InteractionData> loadPhotonDataToMap (const std::string& filename, const float density)
{
    std::map<float, InteractionData> dataMap;
    std::ifstream file(filename);
    
    if (!file.is_open ()) {
        std::cerr << "Error: Could not open file " << filename << std::endl;
        return dataMap;
    }

    std::string line;
    
    // Skip header lines (first 5 lines)
    for (int i = 0; i < 3; ++i) {
        if (!std::getline (file, line)) {
            std::cerr << "Error: File has fewer than 3 header lines" << std::endl;
            return dataMap;
        }
    }
    
    int identicalEnergy = 0;
    while (std::getline (file, line)) {
        // Clean the line
        line.erase (line.begin (), std::find_if (line.begin (), line.end (), [] (int ch) { return !std::isspace (ch); }));
        line.erase (std::find_if (line.rbegin (), line.rend (), [] (int ch) { return !std::isspace (ch); }).base (), line.end ());
        
        if (line.empty()) continue;
        
        // Replace 'E' with 'e' for scientific notation parsing
        std::transform(line.begin (), line.end (), line.begin (), [] (char c) { return c == 'E' ? 'e' : c; });
        
        std::istringstream iss (line);
        float energy;
        float Nuclear;
        float Electron;
        float compton;
        float photoelAbsorb;
        InteractionData entry;
        
        if (iss >> energy >> compton 
                >> photoelAbsorb >> Nuclear 
                >> Electron) {
            // Insert or update the map entry
            entry.pairProd = (Nuclear + Electron) * density; // Nuclear pair production cross-section (cm²/g)
            entry.incoherentScatter = compton * density; // Incoherent scattering cross-section (cm²/g)
            entry.photoelAbsorb = photoelAbsorb * density; // Photoelectric absorption cross-section (cm²/g)
            if (dataMap.find (energy) != dataMap.end ()) {
                identicalEnergy += 1;
                energy += identicalEnergy * 1e-6f; // Adjust energy to avoid duplicates
            } else {
                identicalEnergy = 0; // Reset if no duplicates found
            }

            dataMap[energy] = entry;
        } else {
            std::cerr << "Warning: Failed to parse line: " << line << std::endl;
        }
    }
    
    return dataMap;
}