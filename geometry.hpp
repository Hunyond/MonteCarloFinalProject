#pragma once
#include <cmath>
#include <concepts>
#include <utility>

template<typename T>
concept RandomNumberGenerator = requires(T t) {
    { t() } -> std::convertible_to<double>;
};

struct Vector
{
    float x;
    float y;
    float z;
};

constexpr float myMPI = 3.1415927f;

template<RandomNumberGenerator GEN>
const Vector GetIsotropicDirectionMarsaglia (GEN& generateRandomNumber)
{
    float length_squared = 2.0f;
    float u = 0;
    float v = 0;
    while (length_squared > 1) {
        u = 2 * generateRandomNumber () - 1;
        v = 2 * generateRandomNumber () - 1;
        length_squared = u * u + v * v;
    }
    const float factor = std::sqrt(1 - length_squared);
    return {1 - 2 * length_squared, 2 * u * factor, 2 * v * factor};
}

template<RandomNumberGenerator GEN>
const Vector GetIsotropicDirectionInAngle (const float alpha, GEN& generateRandomNumber)
{
    const float nz = std::cos (alpha) + (1 - std::cos (alpha)) * generateRandomNumber (); // cos(theta)
    const float theta = std::acos (nz); // polar angle
    const float beta = 2 * myMPI * generateRandomNumber (); // azimuthal angle
    return  {std::sin (theta) * std::cos (beta), std::sin (theta) * std::sin (beta), nz};
}


Vector TransfromDirection (const Vector& direction, const Vector& axis);
//Can only be used for planes parallel to the x-y plane
float GetDistanceToPlane (const Vector& point, const Vector& dir, const float ZCoord);

std::pair<float, float> inline GetDistanceToCylinderMantle (const Vector& point, const Vector& dir, const float R);

float GetDistanceToCylinderIn (const Vector& point, const Vector& dir, const float R, const float topOfCyl, const float botOfCyl);

float GetDistanceToCylinderOut (const Vector& point, const Vector& dir, const float R, const float topOfCyl, const float botOfCyl);

std::pair<bool, Vector> HitsCylinder (const Vector& point, const Vector& dir, const float R, const float topOfCyl, const float botOfCyl);