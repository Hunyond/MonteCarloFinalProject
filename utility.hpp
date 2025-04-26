#include <vector>
#include <map>
#include <string>

struct coordinate {
    float x;
    float y;
    float z;
};

void ApplyFWHM (std::vector<float>& data, float fwhm);
float GetStatisticalUncertainty (const std::vector<float>& data);
std::vector<float> linspace (double start, double end, size_t num_points);
std::vector<coordinate> linspace3D (const coordinate start, const coordinate end, const size_t num_points);
std::map<float, int> CreateHistogram( std::vector<float>& data, double min, double max, int num_bins);
void WriteHistogramToFile( std::map<float, int>& histogram, const std::string& filename);