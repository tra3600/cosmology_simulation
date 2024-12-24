#include "Cosmology.h"
#include <cmath>

Cosmology::Cosmology(double H0, double Omega_m, double Omega_Lambda)
    : H0(H0), Omega_m(Omega_m), Omega_Lambda(Omega_Lambda) {}

double Cosmology::scaleFactor(double time) const {
    // Simple scale factor evolution assuming a matter-dominated universe
    return std::pow((3.0 / 2.0) * H0 * std::sqrt(Omega_m) * time, 2.0 / 3.0);
}

double Cosmology::hubbleParameter(double time) const {
    double a = scaleFactor(time);
    return H0 * std::sqrt(Omega_m / std::pow(a, 3) + Omega_Lambda);
}

double Cosmology::friedmannEquation(double a) const {
    return H0 * std::sqrt(Omega_m / std::pow(a, 3) + Omega_Lambda);
}

std::vector<double> Cosmology::simulateExpansion(double startTime, double endTime, double timeStep) const {
    std::vector<double> scaleFactors;
    for (double t = startTime; t <= endTime; t += timeStep) {
        scaleFactors.push_back(scaleFactor(t));
    }
    return scaleFactors;
}