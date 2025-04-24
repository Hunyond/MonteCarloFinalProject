
#include <cmath>
#include <random>
#include <utility>
#include "geometry.hpp"

struct InteractionData {
    float incoherentScatter; // Incoherent scattering cross-section (cm²/g)
    float photoelAbsorb;    // Photoelectric absorption cross-section (cm²/g)
    float pairProd;    // Nuclear pair production cross-section (cm²/g)
};


template<RandomNumberGenerator GEN>
float TrackPhotonAfterPairProction (GEN& getRandomNumber, const Vector& position, const Vector&  direction,const float energy_in, const float R, const float H, const float Ro, const float sigmaCompton, const float sigmaAbsorb, const float sigmaPairProd,)
{
    const float sigma = sigmaCompton + sigmaPairProd + sigmaAbsorb;
    Vector currenctPosition = position;
    Vector currentDirection = direction;
    float energy = energy_in;
    float energyDeposit = 0.0f;

    bool isPhotonAlive = true;
    float distanceTravelled = 0.0f;
    float distanceToCylinder = 0.0f;
    while (isPhotonAlive) {
        distanceTravelled = -std::log(getRandomNumber()) / sigma; 
        distanceToCylinder = GetDistanceToCylinderIn (currenctPosition, currentDirection, R, H, Ro); // Get the distance to the cylinder
        if (distanceToCylinder < distanceTravelled) {
            isPhotonAlive = false; // exiits cylinder
        }
        currenctPosition.x += currentDirection.x * distanceToCylinder;
        currenctPosition.y += currentDirection.y * distanceToCylinder;
        currenctPosition.z += currentDirection.z * distanceToCylinder;
        float rand = getRandomNumber();
        
        currentDirection = GetIsotropicDirectionMarsaglia (GEN& getRandomNumber);
    }

}
template<RandomNumberGenerator GEN>
std::pair<float, float> PhotonAngleAndEnergy(GEN& getRandomNumber, float energy_in) {
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
    
    const float angle = std::acos(1.0f + b - b * f);
    const float energy_out = energy_in / f;
    
    return {angle, energy_out};
}

Vector DirectionInComptonScatter (const float angle)
{
    const float nz = std::cos(angle);
    const float rho = std::sin(angle);

    return {rho * std::cos (phi), rho * std::sin (phi), nz};
}

template<RandomNumberGenerator GEN>
std::pair<Vector, float> ComptonScatter (GEN& getRandomNumber, const Vector& direction, const float energy_in)
{
    const auto [angle, energy_out] = PhotonAngleAndEnergy (getRandomNumber, energy_in);
    const Vector newDirection = DirectionInComptonScatter (angle);
    const Vector newDirectionInParticlesCoordinateSystem = TransfromDirection (newDirection, direction); // Transform to the original coordinate system
    return {newDirectionInParticlesCoordinateSystem, energy_out};
}

template<RandomNumberGenerator GEN>
float PairProduction (GEN& getRandomNumber, const Vector& direction, const float energy_in)
{
    float ToTalEnergy = energy_in - 1.022; // MeV

}