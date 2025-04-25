#pragma once
#include <mutex>
#include <vector>
#include <map>

#include <cmath>
#include <thread>
#include <iostream>
#include <algorithm>
#include <random>
#include "geometry.hpp"
#include "crosssections.hpp"

template<RandomNumberGenerator GEN>
void TrackPhoton (GEN& getRandomNumber, const Vector& position, const Vector&  direction,const float energy_in, const std::map<float, InteractionData>& corssSections, std::vector<float>&  results, std::mutex& results_mutex,const float R, const float H, const float Ro);

template<RandomNumberGenerator GEN>
std::pair<float, float> PhotonAngleAndEnergy (GEN& getRandomNumber, float energy_in);


template<RandomNumberGenerator GEN>
std::pair<Vector, float> ComptonScatter (GEN& getRandomNumber, const Vector& direction, const float energy_in);


template<RandomNumberGenerator GEN>
void PairProduction (GEN& getRandomNumber, const Vector& position, const Vector&  direction,const float energy_in, const std::map<float, InteractionData>& crossSections, std::vector<float>&  results, std::mutex& results_mutex,const float R, const float H, const float Ro);



template<RandomNumberGenerator GEN>
void TrackPhoton (GEN& getRandomNumber, const Vector& position, const Vector&  direction,const float energy_in, const std::map<float, InteractionData>& corssSections, std::vector<float>&  results, std::mutex& results_mutex,const float R, const float H, const float Ro)
{
    Vector currenctPosition = position;
    Vector currentDirection = direction;
    float energy = energy_in;
    float energyDeposit = 0.0f;

    bool isPhotonAlive = true;
    float distanceTravelled = 0.0f;
    float distanceToCylinder = 0.0f;

    InteractionData currentCrossSection = getCrossSectionsAtEnergy (corssSections, energy_in);
    float sigma = currentCrossSection.incoherentScatter + currentCrossSection.photoelAbsorb + currentCrossSection.pairProd; // Total cross-section
    std::vector<std::pair<uint16_t, float>> interactionList = {
        {1, currentCrossSection.incoherentScatter},
        {2, currentCrossSection.photoelAbsorb},
        {3, currentCrossSection.pairProd}
    };
    std::sort (interactionList.begin (), interactionList.end (), [](const auto& a, const auto& b) { return a.second < b.second; });

    std::cout << "Photon energy: " << energy_in << " MeV" << std::endl;

    while (isPhotonAlive) {
        distanceTravelled = -std::log ((float)getRandomNumber ()) / sigma; 
        distanceToCylinder = GetDistanceToCylinderIn (currenctPosition, currentDirection, R, H, Ro); // Get the distance to the cylinder
        if (distanceToCylinder < distanceTravelled) {
            isPhotonAlive = false; // exiits cylinder
        }
        currenctPosition.x += currentDirection.x * distanceToCylinder;
        currenctPosition.y += currentDirection.y * distanceToCylinder;
        currenctPosition.z += currentDirection.z * distanceToCylinder;

        float rand = getRandomNumber ();
        std::pair<Vector, float> res = {Vector{0.0f, 0.0f, 0.0f}, 0.0f};
        
        currentDirection = GetIsotropicDirectionMarsaglia (getRandomNumber);
        for (const auto& interaction : interactionList) {
            if (rand < interaction.second / sigma) {
                switch (interaction.first) {
                    case 1: // Compton scattering
                        res = ComptonScatter (getRandomNumber, currentDirection, energy);
                        currentDirection = res.first; // Update direction after scattering
                        energy = res.second; // Update energy after scattering
                        energyDeposit += energy_in - res.second; // Energy deposited in the material
                        currentCrossSection = getCrossSectionsAtEnergy (corssSections, energy);
                        sigma = currentCrossSection.incoherentScatter + currentCrossSection.photoelAbsorb + currentCrossSection.pairProd; // Total cross-section
                        interactionList = {
                            {1, currentCrossSection.incoherentScatter},
                            {2, currentCrossSection.photoelAbsorb},
                            {3, currentCrossSection.pairProd}
                        };
                        std::sort (interactionList.begin (), interactionList.end (), [](const auto& a, const auto& b) { return a.second < b.second; });
                        break;
                    case 2: // Photoelectric absorption
                        energyDeposit += energy; // Energy deposited in the material
                        isPhotonAlive = false; // Photon is absorbed
                        break;
                    case 3: // Pair production
                        energyDeposit += 1.022; // Energy deposited in the material
                        isPhotonAlive = false; // Photon is absorbed
                        PairProduction (getRandomNumber, currenctPosition, corssSections, results, results_mutex, R, H, Ro); // Pair production
                        break;
                }
                break;
            }
        }
    }
    std::cout << "fianlly here" << energyDeposit << " MeV" << std::endl;
    std::lock_guard<std::mutex> lock (results_mutex);
    results.push_back (energyDeposit); // Store the energy deposit in the results vector

}

template<RandomNumberGenerator GEN>
std::pair<float, float> PhotonAngleAndEnergy (GEN& getRandomNumber, float energy_in)
{
    constexpr float electron_rest_energy = 0.511f; // MeV
    const float a = energy_in / electron_rest_energy;
    const float b = 1.0f / a;
    const float c = 1.0f + 2.0f * a;
    const float d = c / (9.0f + 2.0f * a);
    
    float f = 0.0f;
    float g = 0.0f;
    while (true) {
        const float r1 = getRandomNumber ();
        const float r2 = getRandomNumber ();
        const float r3 = getRandomNumber ();
        
        const float e = 1.0f + 2.0f * a * r1;
        
        if (r2 < d) {
            f = e;
            const float inv_f = 1.0f / f;
            g = 4.0f * (inv_f - inv_f * inv_f);
        } else {
            f = c / e;
            const float term = 1.0f + b - b * f;
            g = 0.5f * (term * term + 1.0f / f);
        }
        
        if (r3 < g) break;
    }
    
    const float angle = std::acos (1.0f + b - b * f);
    const float energy_out = energy_in / f;
    
    return {angle, energy_out};
}

template<RandomNumberGenerator GEN>
Vector DirectionInComptonScatter (GEN& getRandomNumber, const float angle)
{
    const float nz = std::cos (angle);
    const float rho = std::sin (angle);
    const float phi = 2 * M_PI * getRandomNumber ();

    return {rho * std::cos (phi), rho * std::sin (phi), nz};
}

template<RandomNumberGenerator GEN>
std::pair<Vector, float> ComptonScatter (GEN& getRandomNumber, const Vector& direction, const float energy_in)
{
    const auto [angle, energy_out] = PhotonAngleAndEnergy (getRandomNumber, energy_in);
    const Vector newDirection = DirectionInComptonScatter (getRandomNumber, angle);
    const Vector newDirectionInParticlesCoordinateSystem = TransfromDirection (newDirection, direction); // Transform to the original coordinate system
    return {newDirectionInParticlesCoordinateSystem, energy_out};
}

template<RandomNumberGenerator GEN>
void PairProduction (GEN& getRandomNumber, const Vector& position, const std::map<float, InteractionData>& crossSections, std::vector<float>&  results, std::mutex& results_mutex,const float R, const float H, const float Ro)
{
    Vector direction1 = GetIsotropicDirectionMarsaglia (getRandomNumber);
    Vector direction2 = {-1 * direction1.x, -1 * direction1.y, -1 * direction1.z};
    TrackPhoton (getRandomNumber, position, direction1, 0.511f, crossSections, results, results_mutex, R, H, Ro);
    TrackPhoton (getRandomNumber, position, direction2, 0.511f, crossSections, results, results_mutex, R, H, Ro);
}