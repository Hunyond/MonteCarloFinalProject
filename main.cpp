#include <iostream>
#include <vector>
#include <thread>
#include <mutex>
#include "geometry.hpp"
#include <chrono>
#include <random>
#include <map>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <string>
#include "interactions.hpp"
#include <cctype> // for std::isspace
#include <iterator> // for std::find_if
// #include <utility> // for std::pair  


std::map<float, InteractionData> loadPhotonDataToMap(const std::string& filename, const float density) {
    std::map<float, InteractionData> dataMap;
    std::ifstream file(filename);
    
    if (!file.is_open()) {
        std::cerr << "Error: Could not open file " << filename << std::endl;
        return dataMap;
    }

    std::string line;
    
    // Skip header lines (first 5 lines)
    for (int i = 0; i < 3; ++i) {
        if (!std::getline(file, line)) {
            std::cerr << "Error: File has fewer than 3 header lines" << std::endl;
            return dataMap;
        }
    }
    
    int identicalEnergy = 0;
    while (std::getline(file, line)) {
        // Clean the line
        line.erase(line.begin(), std::find_if(line.begin(), line.end(), [](int ch) { return !std::isspace(ch); }));
        line.erase(std::find_if(line.rbegin(), line.rend(), [](int ch) { return !std::isspace(ch); }).base(), line.end());
        
        if (line.empty()) continue;
        
        // Replace 'E' with 'e' for scientific notation parsing
        std::transform(line.begin(), line.end(), line.begin(), [](char c) { return c == 'E' ? 'e' : c; });
        
        std::istringstream iss(line);
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
            if (dataMap.find(energy) != dataMap.end()) {
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



void RunMonteCarloSimulation (int32_t seed, int32_t numberOfNeutrons, const Vector& source, const std::map<float, InteractionData>& crossSections, std::vector<float>& results , std::mutex& resultMutex, const float E, const float R, const float H, const float Ro)
{
    std::mt19937 generator(seed);
    std::uniform_real_distribution<float> randomNumber(0.0f, 1.0f);
    auto getRandomNumber = [&]() {
        return randomNumber(generator);
    };
    const float rg = std::sqrt (R * R + (H * H) / 4.0f);
    const float os = std::sqrt (source.x * source.x + source.y * source.y + source.z * source.z);
    const float alpha = std::atan2 (rg, os);
    Vector direction = GetIsotropicDirectionInAngle (alpha,getRandomNumber);
    const auto res = HitsCylinder (source, direction, R, H/2.0f, -H/2.0f);
    if (!res.first) {
        return;
    }
    Vector startingPosition = res.second;
    TrackPhoton (generator, source, direction, E, crossSections, results, resultMutex, R, H, Ro);
    void TrackPhoton (getRandomNumber, startingPosition, direction, E, crossSections, std::vector<float>& results, resultMutex, R, H, Ro)
}





int main (void)
{



    std::vector<float> results;
    std::mutex results_mutex;
    std::vector<std::thread> threads;
    std::map<float, InteractionData> crossSections = loadPhotonDataToMap ("corsssections.txt");
    const Vector source = {0.0f, 0.0f, 0.0f};
    const float E = 1.0f; // Energy in MeV
    const float R = 1.0f; // Radius of the cylinder in cm
    const float H = 10.0f; // Height of the cylinder in cm
    const float Ro = 0.0f; // Bottom of the cylinder in cm
    const int32_t numberOfNeutrons = 1000; // Number of neutrons to simulate
    unsigned int num_threads = std::thread::hardware_concurrency();
    threads.reserve(num_threads);

    std::random_device rd;
    std::vector<std::uint32_t> seed_values;
    for (unsigned int i = 0; i < num_threads; ++i) {
        seed_values.push_back (rd ());
    }


    const unsigned int numPerThread = numberOfNeutrons / num_threads;
    const unsigned int remainder = numberOfNeutrons % num_threads;
    for (unsigned int i = 0; i < num_threads; ++i) {
        if (i == num_threads - 1) {
            numPerThread += remainder; // Last thread takes the remainder
        }
        threads.emplace_back(RunMonteCarloSimulation, seed_values[i], numPerThread, source, crossSections, results, resultMutex, E, R, H, Ro);
    }


    for (auto& thread : threads) {
        thread.join ();
    }
    for (const auto& result : results) {
        std::cout << result << std::endl;
    }
    std::cout << "Number of results: " << results.size() << std::endl;
}