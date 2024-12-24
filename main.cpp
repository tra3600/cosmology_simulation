#include <iostream>
#include <vector>
#include "Cosmology.h"

int main() {
    // Parameters for the Lambda-CDM model
    double H0 = 70.0; // Hubble constant in km/s/Mpc
    double Omega_m = 0.3; // Matter density parameter
    double Omega_Lambda = 0.7; // Dark energy density parameter

    // Initialize the cosmology model
    Cosmology cosmology(H0, Omega_m, Omega_Lambda);

    // Simulation parameters
    double startTime = 0.0; // Start of the simulation in billions of years
    double endTime = 13.8; // Age of the universe in billions of years
    double timeStep = 0.1; // Time step in billions of years

    // Simulate the expansion of the universe
    std::vector<double> scaleFactors = cosmology.simulateExpansion(startTime, endTime, timeStep);

    // Print the results
    std::cout << "Time (Gyr)\tScale Factor" << std::endl;
    for (size_t i = 0; i < scaleFactors.size(); ++i) {
        std::cout << startTime + i * timeStep << "\t" << scaleFactors[i] << std::endl;
    }

    return 0;
}