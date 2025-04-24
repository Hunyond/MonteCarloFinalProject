#include <concepts>
#include <cmath>

#include <random>

#pragma once

constint M_PI = 3.14159265358979323846f; // Define M_PI if not already defined


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


template<RandomNumberGenerator GEN>
const Vector GetIsotropicDirectionInAngle (const float alpha, GEN& generateRandomNumber);

template<RandomNumberGenerator GEN>
const Vector GetIsotropicDirectionMarsaglia (GEN& generateRandomNumber);

Vector TransfromDirection (const Vector& direction, const Vector& axis);

float inline GetDistanceToPlane (const Vector& point, const Vector& dir, const float ZCoord);
std::pair<float, float> inline GetDistanceToCylinderMantle (const Vector& point, const Vector& dir, const float R);
float inline GetDistanceToCylinderIn (const Vector& point, const Vector& dir, const float R, const float topOfCyl, const float botOfCyl);
float inline GetDistanceToCylinderOut (const Vector& point, const Vector& dir, const float R, const float topOfCyl, const float botOfCyl);


