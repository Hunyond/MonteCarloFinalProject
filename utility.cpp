
#include <cmath>
#include <fstream>
#include <iostream>
#include <numeric>
#include <random>
#include "utility.hpp"



std::vector<float> linspace (double start, double end, size_t num_points)
{
    std::vector<float> result;
    result.reserve (num_points);
    
    if (num_points == 0) return result;
    if (num_points == 1) {
        result.push_back (start);
        return result;
    }
    
    for (size_t i = 0; i < num_points; ++i) {
        double t = static_cast<double>(i) / (num_points - 1);
        result.push_back (std::lerp (start, end, t)); 
    }
    
    return result;
}

std::vector<coordinate> linspace3D (const coordinate start, const coordinate end, const size_t num_points)
{
    std::vector<coordinate> result;
    result.reserve(num_points);
    
    if (num_points == 0) return result;
    if (num_points == 1) {
        result.push_back(start);
        return result;
    }
    
    const float step_x = (end.x - start.x) / (num_points - 1);
    const float step_y = (end.y - start.y) / (num_points - 1);
    const float step_z = (end.z - start.z) / (num_points - 1);
    
    for (size_t i = 0; i < num_points; ++i) {
        coordinate point;
        point.x = start.x + i * step_x;
        point.y = start.y + i * step_y;
        point.z = start.z + i * step_z;
        result.push_back (point);
    }
    
    return result;
}

std::map<float, int> CreateHistogram( std::vector<float>& data, double min, double max, int num_bins)
{
    std::map<float, int> histogram;
    double bin_width = (max - min) / num_bins;

    for (const auto& value : data) {
        if (value >= min && value <= max) {
            int bin_index = static_cast<int>((value - min) / bin_width);
            if (bin_index == num_bins) {
                bin_index--;
            }
            float bin_start = min + bin_index * bin_width;
            histogram[bin_start]++;
        }
    }

    return histogram;
}

void WriteHistogramToFile( std::map<float, int>& histogram, const std::string& filename)
{
    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Could not open file " << filename << std::endl;
        return;
    }

    for (const auto& [bin_start, count] : histogram) {
        file << bin_start << ";" << count << "\n";
    }

    file.close();
}

float GetStatisticalUncertainty (const std::vector<float>& data)
{
    const auto mean = std::accumulate(data.begin (), data.end (), 0.0f) / data.size ();
    const auto meanSquare = std::accumulate(data.begin (), data.end (), 0.0f, [] (float sum, float value) {
        return sum + value * value;
    }) / data.size ();
    const auto variance = meanSquare - mean * mean;
    return variance;
}

void ApplyFWHM (std::vector<float>& data, float fwhm)
{
    const auto sigma = fwhm / 2.3548;
    const auto mean = std::accumulate(data.begin (), data.end (), 0.0f) / data.size ();
    std::random_device rd;
    std::mt19937 gen (rd ());
    std::normal_distribution<> dist (mean, sigma);

    for (auto& value : data) {
        value += dist (gen);
    }
}