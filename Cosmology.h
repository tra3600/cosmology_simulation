#ifndef COSMOLOGY_H
#define COSMOLOGY_H

#include <vector>

class Cosmology {
public:
    Cosmology(double H0, double Omega_m, double Omega_Lambda);
    double scaleFactor(double time) const; // Calculates the scale factor at a given time
    double hubbleParameter(double time) const; // Calculates the Hubble parameter at a given time
    std::vector<double> simulateExpansion(double startTime, double endTime, double timeStep) const;

private:
    double H0; // Hubble constant at present time (km/s/Mpc)
    double Omega_m; // Density parameter for matter
    double Omega_Lambda; // Density parameter for dark energy
    double friedmannEquation(double a) const; // Friedmann equation for the scale factor
};

#endif // COSMOLOGY_H