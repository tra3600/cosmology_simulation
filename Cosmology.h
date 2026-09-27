#ifndef COSMOLOGY_H
#define COSMOLOGY_H

#include <string>
#include <vector>

// ---------------------------------------------------------------------------
//  Constantes physiques et conversions d'unités
// ---------------------------------------------------------------------------
namespace cosmo {

constexpr double kSpeedOfLight = 299792.458;   // c en km/s
constexpr double kHubbleToGyr  = 1.0227121650e-3; // 1 km/s/Mpc = 1.0227e-3 Gyr^-1
constexpr double kPi           = 3.14159265358979323846;
constexpr double kGravitation  = 6.67430e-11;  // G en m^3 kg^-1 s^-2
constexpr double kMpcInMeters  = 3.0856775814913673e22;

// Conversion d'une constante de Hubble en km/s/Mpc vers des Gyr^-1
inline double hubbleToInverseGyr(double H) { return H * kHubbleToGyr; }

// ---------------------------------------------------------------------------
//  Paramètres d'un modèle cosmologique (Lambda-CDM étendu w0-wa, courbure)
// ---------------------------------------------------------------------------
struct CosmoParams {
    std::string name        = "planck2018";
    std::string description = "Planck 2018 (TT,TE,EE+lowE+lensing), LambdaCDM plat";
    double H0          = 67.36;   // constante de Hubble (km/s/Mpc)
    double Omega_m     = 0.3153;  // matière (baryons + matière noire froide)
    double Omega_Lambda= 0.6847;  // énergie sombre aujourd'hui
    double Omega_r     = -1.0;    // rayonnement ; < 0 => calculé depuis T_cmb et N_eff
    double T_cmb       = 2.7255;  // température du fond diffus cosmologique (K)
    double N_eff       = 3.046;   // nombre effectif d'espèces de neutrinos
    double w0          = -1.0;    // équation d'état de l'énergie sombre, paramétrisation
    double wa          = 0.0;     //   CPL : w(a) = w0 + wa (1 - a)
    double sigma8      = 0.8111;  // amplitude des fluctuations de matière aujourd'hui
    bool   flat        = true;    // si vrai, Omega_Lambda = 1 - Omega_m - Omega_r
};

// Modèles prédéfinis (Planck, tension de Hubble, DESI, Big Rip, Big Crunch...)
const std::vector<CosmoParams>& presets();
bool findPreset(const std::string& name, CosmoParams& out);

// Un point de la trajectoire a(t) obtenue par intégration numérique
struct ExpansionState {
    double t;      // temps cosmique depuis le Big Bang (Gyr)
    double a;      // facteur d'échelle (a = 1 aujourd'hui)
    double adot;   // da/dt (Gyr^-1)
    double H;      // paramètre de Hubble (km/s/Mpc), négatif en contraction
    double q;      // paramètre de décélération
};

// Destin de l'univers déterminé par l'intégration dans le futur
enum class Fate { EternalAcceleration, EternalDeceleration, EternalCoasting, BigCrunch, BigRip };
std::string fateToString(Fate fate);

struct ExpansionHistory {
    std::vector<ExpansionState> states;
    Fate   fate;
    double fateTime = -1.0; // temps (Gyr) du Big Crunch / Big Rip, sinon -1
};

// ---------------------------------------------------------------------------
//  Modèle cosmologique : équations de Friedmann et observables
// ---------------------------------------------------------------------------
class Cosmology {
public:
    explicit Cosmology(const CosmoParams& params);
    // Compatibilité avec l'ancienne interface (univers sans rayonnement)
    Cosmology(double H0, double Omega_m, double Omega_Lambda);

    const CosmoParams& params() const { return p; }
    double H0() const { return p.H0; }
    double omegaM() const { return p.Omega_m; }
    double omegaR() const { return Omega_r; }
    double omegaDE() const { return p.Omega_Lambda; }
    double omegaK() const { return Omega_k; }

    // --- Dynamique ---------------------------------------------------------
    double w(double a) const;                  // équation d'état de l'énergie sombre
    double darkEnergyScaling(double a) const;  // rho_DE(a) / rho_DE(1)
    double E2(double a) const;                 // (H/H0)^2
    double E(double a) const;                  // H/H0
    double hubble(double a) const;             // H(a) en km/s/Mpc
    double deceleration(double a) const;       // q(a) = -a a'' / a'^2
    double omegaMatter(double a) const;        // fractions de densité à l'époque a
    double omegaRadiation(double a) const;
    double omegaDarkEnergy(double a) const;
    double omegaCurvature(double a) const;

    // --- Temps (Gyr) -------------------------------------------------------
    double hubbleTime() const;                 // 1/H0
    double age(double a) const;                // âge de l'univers quand le facteur d'échelle vaut a
    double ageToday() const { return age(1.0); }
    double ageAtRedshift(double z) const { return age(1.0 / (1.0 + z)); }
    double lookbackTime(double z) const;

    // Ancienne interface : facteur d'échelle et H au temps t (Gyr), par intégration
    double scaleFactor(double time) const;
    double hubbleParameter(double time) const;
    std::vector<double> simulateExpansion(double startTime, double endTime, double timeStep) const;

    // --- Distances (Mpc) ---------------------------------------------------
    double comovingDistance(double z) const;
    double transverseComovingDistance(double z) const;
    double luminosityDistance(double z) const;
    double angularDiameterDistance(double z) const;
    double distanceModulus(double z) const;    // mu = 5 log10(D_L / 10 pc)
    double particleHorizon(double a) const;    // horizon des particules (Mpc comobiles)

    // --- Croissance des structures (régime linéaire) -----------------------
    double growthFactor(double a) const;       // D(a), normalisé à D(1) = 1
    double growthRate(double a) const;         // f = dlnD/dlna
    double fSigma8(double z) const;            // f(z) sigma8(z), mesuré par les relevés de galaxies

    // --- Époques clés (redshift ; -1 si inexistante) -----------------------
    double redshiftMatterRadiationEquality() const;
    double redshiftMatterDarkEnergyEquality() const;
    double redshiftAccelerationOnset() const;

    // --- Divers -------------------------------------------------------------
    double criticalDensity() const;            // rho_c aujourd'hui (kg/m^3)
    double cmbTemperature(double z) const { return p.T_cmb * (1.0 + z); }
    bool   hasBigBang(std::string* why = nullptr) const;

    // Intègre a(t) (équation de Friedmann du second ordre, RK4 adaptatif)
    // depuis les premiers instants jusqu'à t = tMax (Gyr) ou la singularité finale.
    ExpansionHistory evolve(double tMax) const;

private:
    CosmoParams p;
    double Omega_r;
    double Omega_k;
    double H0inv;  // H0 en Gyr^-1

    // Table de croissance D(a) calculée une fois à la construction
    std::vector<double> growthLnA, growthD, growthF;
    void computeGrowth();
    double interpolateGrowth(const std::vector<double>& table, double a) const;
    double acceleration(double a) const;       // a'' en Gyr^-2
};

// Facteur d'échelle a(t) interpolé linéairement sur une trajectoire intégrée
double interpolateScaleFactor(const ExpansionHistory& h, double t);

} // namespace cosmo

#endif // COSMOLOGY_H
