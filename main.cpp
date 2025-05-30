#include <iterator>
#include "interactions.hpp"
#include "utility.hpp"
#include <chrono>


void RunMonteCarloSimulation (int32_t seed, long long numberOfNeutrons, const Vector& source, const std::map<float, InteractionData>& crossSections, std::vector<float>& results , std::mutex& resultMutex, std::atomic<int>& misses, const float E, const float R, const float H, const float alpha)
{
    std::mt19937 generator(seed);
    std::uniform_real_distribution<float> randomNumber(0.0f, 1.0f);
    auto getRandomNumber = [&generator, &randomNumber]() { 
        return randomNumber(generator); 
    };

    for (int32_t i = 0; i < numberOfNeutrons; ++i) {
        Vector direction = GetIsotropicDirectionInAngle (alpha, getRandomNumber);
        direction = TransfromDirection (direction, {-source.x, -source.y, -source.z}); // Transform to the original coordinate system
        const auto res = HitsCylinder (source, direction, R, H/2.0f, -H/2.0f);
        if (!res.first) {
            misses++;
            continue; // Missed the cylinder;
        }
        Vector startingPosition = res.second;
        TrackPhoton (getRandomNumber, startingPosition, direction, E, crossSections, results, resultMutex, R, H);
    }
}


std::pair<float, float> PrepareSimulation (const int simId, const int numPhotons, const Vector& source, const std::map<float, InteractionData>& crossSections, const float E, const float R, const float H, const float FWHM)
{
    std::mutex results_mutex;
    std::vector<std::thread> threads;
    std::atomic<int> misses = 0;


    unsigned int num_threads = std::thread::hardware_concurrency ();
    threads.reserve (num_threads);
    std::vector<float> results;
    results.reserve (numPhotons / 2);


    std::random_device rd;
    std::vector<std::uint32_t> seed_values;
    for (unsigned int i = 0; i < num_threads; ++i) {
        seed_values.push_back (rd ());
    }

    const float rg = std::sqrt (R * R + (H * H) / 4.0f);
    const float os = std::sqrt (source.x * source.x + source.y * source.y + source.z * source.z);
    const float alpha = std::atan2 (rg, os);


    std::cout << "----------------------------------------------------------------------" <<
                 std::endl <<"starting simulation with " << num_threads << " threads" << std::endl;
    std::chrono::high_resolution_clock::time_point start = std::chrono::high_resolution_clock::now ();
    unsigned int baseNumPerThread = numPhotons / num_threads;
    unsigned int remainder = numPhotons % num_threads;
    for (unsigned int i = 0; i < num_threads; ++i) {
        unsigned long long numNPhotonsForThread = baseNumPerThread + (i < remainder ? 1 : 0);
        threads.emplace_back (RunMonteCarloSimulation, seed_values[i], numNPhotonsForThread, std::cref(source), std::cref(crossSections), std::ref(results), std::ref(results_mutex), std::ref(misses), E, R, H, alpha);
    }
    for (auto& thread : threads) {
        thread.join ();
    }
    std::chrono::high_resolution_clock::time_point end = std::chrono::high_resolution_clock::now ();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds> (end - start);
    std::cout << "Finished simulation " << std::endl;
    std::cout << "Time taken (ms): " << duration.count () << std::endl << "----------------------------------------------------------------------" << std::endl;

    const float totalEnergyDeposited = std::accumulate (results.begin (), results.end (), 0.0f);
    const float totalEnergyEmitted = numPhotons * E * 2 / (1 - std::cos (alpha));
    const float totalEnergyReached = numPhotons - misses * E;

    const float totalEfficiency = totalEnergyDeposited / totalEnergyEmitted * 100.0f;
    const float interactionEfficiency = totalEnergyDeposited / totalEnergyReached * 100.0f;


    ApplyFWHM (results, FWHM);
    const auto var = GetStatisticalUncertainty (results);

    std::cout << "Statistical uncertainty: " << var << std::endl;


    auto histogram = CreateHistogram (results, 0.0f, E + 0.1 * E, 1024);
    std::string filename = "histogram_" + std::to_string (simId) + ".csv";
    WriteHistogramToFile (histogram, filename);
    return {totalEfficiency, interactionEfficiency};
}



/*int main (void)
{

    const Vector source = {3.0f, -3.0f, 2.0f};
    const float E = 0.6617f; // Energy in MeV
    const float R = 2.5f; // Radius of the cylinder in cm
    const float H = 3.0f; // Height of the cylinder in cm
    const float Ro = 3.67f;
    const float FWHM = 6.0 / 1000.0f; // FWHM in MeV
    const long long numberOfNeutrons = 100000000; // Number of neutrons to simulate
    std::map<float, InteractionData> crossSections = loadPhotonDataToMap ("corsssections.txt", Ro); // Load the cross-section data from a file
    std::cout << "Loaded cross-section data." << std::endl;
    std::cout << "Number of neutrons: " << numberOfNeutrons << std::endl;
    std::cout << "Source position: (" << source.x << ", " << source.y << ", " << source.z << ")" << std::endl;
    std::cout << "Energy: " << E << " MeV" << std::endl;

    auto efficecnies = PrepareSimulation (0, numberOfNeutrons, source, crossSections, E, R, H, FWHM);
    std::cout << "Total efficiency: " << efficecnies.first << "%" << std::endl;
    std::cout << "Interaction efficiency: " << efficecnies.second << "%" << std::endl;

}*/


/*int main (void)
{

    const Vector source = {4.0f, 4.0f, 0.0f};
    const float E = 1.3325f; // Energy in MeV
    const float R = 3.0f; // Radius of the cylinder in cm
    const float H = 5.0f; // Height of the cylinder in cm
    const float Ro = 3.67f;
    const float FWHM = 8.0 / 1000.0f; // FWHM in MeV
    const long long numberOfNeutrons = 100000000; // Number of neutrons to simulate
    std::map<float, InteractionData> crossSections = loadPhotonDataToMap ("corsssections.txt", Ro); // Load the cross-section data from a file
    std::cout << "Loaded cross-section data." << std::endl;
    std::cout << "Number of neutrons: " << numberOfNeutrons << std::endl;
    std::cout << "Source position: (" << source.x << ", " << source.y << ", " << source.z << ")" << std::endl;
    std::cout << "Energy: " << E << " MeV" << std::endl;

    auto efficecnies = PrepareSimulation (0, numberOfNeutrons, source, crossSections, E, R, H, FWHM);
    std::cout << "Total efficiency: " << efficecnies.first << "%" << std::endl;
    std::cout << "Interaction efficiency: " << efficecnies.second << "%" << std::endl;

}*/

/*int main (void)
{

    const auto sources = linspace3D (coordinate{1.0f, 3.5f, 2.0f}, coordinate{-4.0, -1.5, 2.0}, 11);
    const float E = 0.6617f; // Energy in MeV
    const float R = 2.5f; // Radius of the cylinder in cm
    const float H = 3.0f; // Height of the cylinder in cm
    const float Ro = 3.67f;
    const float FWHM = 6.0 / 1000.0f; // FWHM in MeV
    const long long numberOfNeutrons = 100000000; // Number of neutrons to simulate
    std::map<float, InteractionData> crossSections = loadPhotonDataToMap ("corsssections.txt", Ro); // Load the cross-section data from a file


    std::vector<float> totelEfficiencies;
    std::vector<float> interactionEfficiencies;

    int simId = 0;
    for (const auto& source : sources) {
        std::cout << "----------------------------------------------------------------------" << std::endl;
        std::cout << "Simulation ID: " << simId << std::endl;
        std::cout << "----------------------------------------------------------------------" << std::endl;
        std::cout << "Source position: (" << source.x << ", " << source.y << ", " << source.z << ")" << std::endl;
        std::cout << "Energy: " << E << " MeV" << std::endl;
        auto efficecnies = PrepareSimulation (simId, numberOfNeutrons, {source.x, source.y, source.z}, crossSections, E, R, H, FWHM);
        std::cout << "Total efficiency: " << efficecnies.first << "%" << std::endl;
        totelEfficiencies.push_back (efficecnies.first);
        interactionEfficiencies.push_back (efficecnies.second);
        std::cout << "Interaction efficiency: " << efficecnies.second << "%" << std::endl;
        simId++;
    }
    std::cout << "----------------------------------------------------------------------" << std::endl;
    std::cout << "Total efficiencies: ";
    for (const auto& eff : totelEfficiencies) {
        std::cout << eff;
    }
    std::cout << std::endl;
    std::cout << "Interaction efficiencies: ";
    for (const auto& eff : interactionEfficiencies) {
        std::cout << eff;
    }

}*/

int main (void)
{

    const Vector source = {4.0f, 4.0f, 0.0f};
    std::vector Energies = linspace(0.4, 4.0, 10);
    const float R = 3.0f; // Radius of the cylinder in cm
    const float H = 5.0f; // Height of the cylinder in cm
    const float Ro = 3.67f;
    const float FWHM = 8.0 / 1000.0f; // FWHM in MeV
    const long long numberOfNeutrons = 100000000; // Number of neutrons to simulate
    std::map<float, InteractionData> crossSections = loadPhotonDataToMap ("corsssections.txt", Ro); // Load the cross-section data from a file

    int cnt = 0;
    std::vector<float> totelEfficiencies;
    std::vector<float> interactionEfficiencies;
    for (const auto& E : Energies) {
        std::cout << "----------------------------------------------------------------------" << std::endl;
        std::cout << "Simulation for energy: " << E << " MeV" << std::endl;
        std::cout << "----------------------------------------------------------------------" << std::endl;
        std::cout << "Source position: (" << source.x << ", " << source.y << ", " << source.z << ")" << std::endl;
        auto efficecnies = PrepareSimulation (cnt, numberOfNeutrons, source, crossSections, E, R, H, FWHM);
        totelEfficiencies.push_back (efficecnies.first);
        interactionEfficiencies.push_back (efficecnies.second);
        cnt++;
    }

    std::cout << "----------------------------------------------------------------------" << std::endl;
    std::cout << "Total efficiencies: ";
    for (const auto& eff : totelEfficiencies) {
        std::cout << eff << ", ";
    }
    std::cout << std::endl;
    std::cout << "Interaction efficiencies: ";
    for (const auto& eff : interactionEfficiencies) {
        std::cout << eff << ", ";
    }

}



