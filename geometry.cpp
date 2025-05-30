#include <random>
#include "geometry.hpp"

Vector TransfromDirection (const Vector& direction, const Vector& axis)
{
    float ax = axis.x;
    float ay = axis.y;
    float az = axis.z;
    float len = std::sqrt (ax * ax + ay * ay + az * az);

    if (len == 0.0f) {
        len = 0.000000001f;
    }
    
    ax = ax / len;
    ay = ay / len;
    az = az / len;

    const float s_squared = (ax * ax + ay * ay);
    if (s_squared < 1e-6f) {
        return az > 0.0f ? direction : Vector{ -direction.x, -direction.y, -direction.z };
    }

    const float s = std::sqrt(s_squared);
    float inv_s = 1.0f / s;
    if (std::isnan(inv_s)) { inv_s = 0.0f; }

    const float t11 = ay * inv_s;
    const float t12 = ax * az * inv_s;
    const float t13 = ax;
    const float t21 = -ax * inv_s;
    const float t22 = ay * az * inv_s;
    const float t23 = ay;
    const float t31 = 0.0f;
    const float t32 = -s;
    const float t33 = az;


    return {t11 * direction.x + t12 * direction.y + t13 * direction.z,
            t21 * direction.x + t22 * direction.y + t23 * direction.z,
            t31 * direction.x + t32 * direction.y + t33 * direction.z};
}


//Can only be used for planes parallel to the x-y plane
float GetDistanceToPlane (const Vector& point, const Vector& dir, const float ZCoord)
{
    const float pz = point.z;
    const float nz = dir.z;
    if (std::abs (nz) < 1e-4f) {
        return INFINITY;
    }
    return (ZCoord - pz) / (nz);
}



std::pair<float, float> GetDistanceToCylinderMantle (const Vector& point, const Vector& dir, const float R)
{
    const float px = point.x;
    const float py = point.y;
    const float dx = dir.x;
    const float dy = dir.y;
    const float a = dx * dx + dy * dy;
    const float b = 2 * (px * dx + py * dy);
    const float c = px * px + py * py - R * R;
    const float discriminant = b * b - 4 * a * c;
    if (discriminant < 0) {
        return {INFINITY, INFINITY};
    }
    const float sqrt_discriminant = std::sqrt(discriminant);
    return { (-b + sqrt_discriminant) / (2 * a), (-b - sqrt_discriminant) / (2 * a)};
}

float GetDistanceToCylinderIn (const Vector& point, const Vector& dir, const float R, const float topOfCyl, const float botOfCyl)
{
    auto [d1, d2] = GetDistanceToCylinderMantle (point, dir, R);
    const float dtop = GetDistanceToPlane (point, dir, topOfCyl);
    const float dbot = GetDistanceToPlane (point, dir, botOfCyl);
    return std::min (std::max (d1, d2), std::max (dtop, dbot));
}

float GetDistanceToCylinderOut (const Vector& point, const Vector& dir, const float R, const float topOfCyl, const float botOfCyl)
{
    const float px = point.x;
    const float py = point.y;
    const float pz = point.z;

    const float dx = dir.x;
    const float dy = dir.y;
    const float dz = dir.z;

    auto [d1, d2] = GetDistanceToCylinderMantle (point, dir, R);
    const float dMinCyl = std::min (d1, d2);
    const float dMaxCyl = std::max (d1, d2);

    if (dMinCyl == INFINITY || dMaxCyl < 0 || (pz > topOfCyl && dz > 0) || (pz < botOfCyl && dz < 0)) {
        return INFINITY;
    }

    float dCyl = 0;

    if (d1 * d2 > 0) {
        dCyl = dMinCyl;
    } else {
        dCyl = dMaxCyl;
    }

    const float zCyl = pz + dCyl * dz;
    if (zCyl < botOfCyl || zCyl > topOfCyl) {
        dCyl = INFINITY;
    }

    float dTop = GetDistanceToPlane (point, dir, topOfCyl);
    const float xTop = px + dTop * dx;
    const float yTop = py + dTop * dy;
    if (xTop * xTop + yTop * yTop > R * R) {
        dTop = INFINITY;
    }

    float dBot = GetDistanceToPlane (point, dir, botOfCyl);
    const float xBot = px + dBot * dx;
    const float yBot = py + dBot * dy;
    if (xBot * xBot + yBot * yBot > R * R) {
        dBot = INFINITY;
    }
    return std::min (dCyl, std::min (dTop, dBot));
}

std::pair<bool, Vector> HitsCylinder (const Vector& point, const Vector& dir, const float R, const float topOfCyl, const float botOfCyl)
{
    const auto res = GetDistanceToCylinderOut (point, dir, R, topOfCyl, botOfCyl);
    if (res == INFINITY) {
        return {false, {0.0f, 0.0f, 0.0f}};
    }
    return {true, {point.x + dir.x * res, point.y + dir.y * res, point.z + dir.z * res}};
}