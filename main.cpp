#include <iterator>
#include "interactions.hpp"

void RunMonteCarloSimulation (int32_t seed, int32_t numberOfNeutrons, const Vector& source, const std::map<float, InteractionData>& crossSections, std::vector<float>& results , std::mutex& resultMutex, const float E, const float R, const float H, const float Ro)
{
    std::cout << "Thread ID: " << std::this_thread::get_id() << " Seed: " << seed << std::endl;
    std::mt19937 generator(seed);
    std::uniform_real_distribution<float> randomNumber(0.0f, 1.0f);
    auto getRandomNumber = [&generator, &randomNumber]() { 
        return randomNumber(generator); 
    };
    const float rg = std::sqrt (R * R + (H * H) / 4.0f);
    const float os = std::sqrt (source.x * source.x + source.y * source.y + source.z * source.z);
    const float alpha = std::atan2 (rg, os);
    for (int32_t i = 0; i < numberOfNeutrons; ++i) {
        
        Vector direction = GetIsotropicDirectionInAngle (alpha, getRandomNumber);
        const auto res = HitsCylinder (source, direction, R, H/2.0f, -H/2.0f);
        if (!res.first) {
            return;
        }
        Vector startingPosition = res.second;
        TrackPhoton (generator, startingPosition, direction, E, crossSections, results, resultMutex, R, H, Ro);
    }
}



int main (void)
{
    std::vector<float> results;
    std::mutex results_mutex;

    std::vector<std::thread> threads;

    const float density = 1.0f; // Density in g/cm³
    std::map<float, InteractionData> crossSections = loadPhotonDataToMap ("corsssections.txt", density);

    const Vector source = {0.0f, 0.0f, 0.0f};
    const float E = 1.0f; // Energy in MeV
    const float R = 1.0f; // Radius of the cylinder in cm
    const float H = 10.0f; // Height of the cylinder in cm
    const float Ro = 0.0f; // Bottom of the cylinder in cm
    const int32_t numberOfNeutrons = 10; // Number of neutrons to simulate
    unsigned int num_threads = std::thread::hardware_concurrency ();
    threads.reserve (num_threads);


    std::random_device rd;
    std::vector<std::uint32_t> seed_values;
    for (unsigned int i = 0; i < num_threads; ++i) {
        seed_values.push_back (rd ());
    }


    unsigned int baseNumPerThread = numberOfNeutrons / num_threads;
    unsigned int remainder = numberOfNeutrons % num_threads;
    for (unsigned int i = 0; i < num_threads; ++i) {
        unsigned int numNeutronsForThread = baseNumPerThread + (i < remainder ? 1 : 0);
        threads.emplace_back(RunMonteCarloSimulation, seed_values[i], numNeutronsForThread, std::ref(source), std::ref(crossSections), std::ref(results), std::ref(results_mutex), E, R, H, Ro);
    }


    for (auto& thread : threads) {
        thread.join ();
    }
    std::cout << "Number of results: " << results.size () << std::endl;
}