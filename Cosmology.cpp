#include "Cosmology.h"

#include <algorithm>
#include <cmath>
#include <limits>

namespace cosmo {

// ---------------------------------------------------------------------------
//  Outils numériques
// ---------------------------------------------------------------------------
namespace {

template <class F>
double simpsonStep(F& f, double a, double b, double fa, double fm, double fb,
                   double whole, double eps, int depth) {
    const double m = 0.5 * (a + b);
    const double lm = 0.5 * (a + m), rm = 0.5 * (m + b);
    const double flm = f(lm), frm = f(rm);
    const double left = (m - a) / 6.0 * (fa + 4.0 * flm + fm);
    const double right = (b - m) / 6.0 * (fm + 4.0 * frm + fb);
    const double delta = left + right - whole;
    if (depth <= 0 || std::fabs(delta) <= 15.0 * eps)
        return left + right + delta / 15.0;
    return simpsonStep(f, a, m, fa, flm, fm, left, 0.5 * eps, depth - 1) +
           simpsonStep(f, m, b, fm, frm, fb, right, 0.5 * eps, depth - 1);
}

// Intégration de Simpson adaptative sur [a, b]
template <class F>
double integrate(F f, double a, double b, double eps = 1e-9, int depth = 30) {
    if (b <= a) return 0.0;
    const double fa = f(a), fb = f(b), fm = f(0.5 * (a + b));
    const double whole = (b - a) / 6.0 * (fa + 4.0 * fm + fb);
    return simpsonStep(f, a, b, fa, fm, fb, whole, eps, depth);
}

// Racine de g sur [lo, hi] par dichotomie (g(lo) et g(hi) de signes opposés)
template <class G>
double bisect(G g, double lo, double hi) {
    double glo = g(lo);
    for (int i = 0; i < 100; ++i) {
        const double mid = 0.5 * (lo + hi);
        const double gm = g(mid);
        if ((gm > 0) == (glo > 0)) { lo = mid; glo = gm; }
        else hi = mid;
    }
    return 0.5 * (lo + hi);
}

constexpr double kTiny = 1e-12;

} // namespace

// ---------------------------------------------------------------------------
//  Modèles prédéfinis
// ---------------------------------------------------------------------------
const std::vector<CosmoParams>& presets() {
    static const std::vector<CosmoParams> list = [] {
        std::vector<CosmoParams> v;
        CosmoParams planck;  // valeurs par défaut de CosmoParams
        v.push_back(planck);

        CosmoParams shoes = planck;
        shoes.name = "sh0es";
        shoes.description = "Tension de Hubble : H0 = 73.04 (SH0ES, céphéides + SN Ia), reste Planck";
        shoes.H0 = 73.04;
        v.push_back(shoes);

        CosmoParams desi = planck;
        desi.name = "desi";
        desi.description = "Énergie sombre dynamique w0waCDM (~DESI DR2 + CMB + SN, valeurs arrondies)";
        desi.H0 = 66.7; desi.Omega_m = 0.319; desi.w0 = -0.75; desi.wa = -0.86;
        v.push_back(desi);

        CosmoParams phantom = planck;
        phantom.name = "phantom";
        phantom.description = "Énergie fantôme w = -1.5 : l'univers finit déchiré (Big Rip)";
        phantom.w0 = -1.5;
        v.push_back(phantom);

        CosmoParams eds = planck;
        eds.name = "eds";
        eds.description = "Einstein-de Sitter : matière seule, plat (modèle historique)";
        eds.Omega_m = 1.0; eds.Omega_Lambda = 0.0; eds.Omega_r = 0.0; eds.flat = false;
        v.push_back(eds);

        CosmoParams open = planck;
        open.name = "open";
        open.description = "Univers ouvert sans énergie sombre (Omega_m = 0.3, courbure négative)";
        open.Omega_m = 0.3; open.Omega_Lambda = 0.0; open.flat = false;
        v.push_back(open);

        CosmoParams closed = planck;
        closed.name = "closed";
        closed.description = "Univers fermé très dense (Omega_m = 3) : il s'effondre (Big Crunch)";
        closed.Omega_m = 3.0; closed.Omega_Lambda = 0.0; closed.flat = false;
        v.push_back(closed);

        CosmoParams milne = planck;
        milne.name = "milne";
        milne.description = "Univers de Milne : vide, expansion à vitesse constante a = t / t0";
        milne.Omega_m = 0.0; milne.Omega_Lambda = 0.0; milne.Omega_r = 0.0; milne.flat = false;
        v.push_back(milne);
        return v;
    }();
    return list;
}

bool findPreset(const std::string& name, CosmoParams& out) {
    for (const auto& p : presets()) {
        if (p.name == name) { out = p; return true; }
    }
    return false;
}

std::string fateToString(Fate fate) {
    switch (fate) {
        case Fate::EternalAcceleration: return "expansion accélérée éternelle (Big Freeze)";
        case Fate::EternalDeceleration: return "expansion éternelle qui ralentit";
        case Fate::EternalCoasting:     return "expansion éternelle à vitesse constante";
        case Fate::BigCrunch:           return "effondrement final (Big Crunch)";
        case Fate::BigRip:              return "déchirement final (Big Rip)";
    }
    return "?";
}

// ---------------------------------------------------------------------------
//  Construction
// ---------------------------------------------------------------------------
Cosmology::Cosmology(const CosmoParams& params) : p(params) {
    if (p.Omega_r < 0.0) {
        // Photons : Omega_gamma h^2 = 2.4728e-5 (T/2.7255)^4 ; neutrinos : x 0.2271 N_eff
        const double h = p.H0 / 100.0;
        const double omegaGamma = 2.4728e-5 * std::pow(p.T_cmb / 2.7255, 4) / (h * h);
        Omega_r = omegaGamma * (1.0 + 0.2271 * p.N_eff);
    } else {
        Omega_r = p.Omega_r;
    }
    if (p.flat) p.Omega_Lambda = 1.0 - p.Omega_m - Omega_r;
    Omega_k = 1.0 - p.Omega_m - Omega_r - p.Omega_Lambda;
    H0inv = hubbleToInverseGyr(p.H0);
    if (hasBigBang()) computeGrowth();
}

Cosmology::Cosmology(double H0, double Omega_m, double Omega_Lambda)
    : Cosmology([&] {
          CosmoParams c;
          c.name = "custom";
          c.description = "Modèle personnalisé";
          c.H0 = H0; c.Omega_m = Omega_m; c.Omega_Lambda = Omega_Lambda;
          c.Omega_r = 0.0; c.flat = false;
          return c;
      }()) {}

// ---------------------------------------------------------------------------
//  Équation de Friedmann
// ---------------------------------------------------------------------------
double Cosmology::w(double a) const { return p.w0 + p.wa * (1.0 - a); }

double Cosmology::darkEnergyScaling(double a) const {
    // rho_DE ∝ a^{-3(1+w0+wa)} exp(-3 wa (1-a)) pour la paramétrisation CPL
    if (p.w0 == -1.0 && p.wa == 0.0) return 1.0;
    return std::pow(a, -3.0 * (1.0 + p.w0 + p.wa)) * std::exp(-3.0 * p.wa * (1.0 - a));
}

double Cosmology::E2(double a) const {
    return Omega_r / (a * a * a * a) + p.Omega_m / (a * a * a) + Omega_k / (a * a) +
           p.Omega_Lambda * darkEnergyScaling(a);
}

double Cosmology::E(double a) const { return std::sqrt(std::max(E2(a), 0.0)); }

double Cosmology::hubble(double a) const { return p.H0 * E(a); }

double Cosmology::deceleration(double a) const {
    const double num = Omega_r / std::pow(a, 4) + 0.5 * p.Omega_m / std::pow(a, 3) +
                       0.5 * p.Omega_Lambda * darkEnergyScaling(a) * (1.0 + 3.0 * w(a));
    return num / E2(a);
}

double Cosmology::omegaMatter(double a) const { return p.Omega_m / std::pow(a, 3) / E2(a); }
double Cosmology::omegaRadiation(double a) const { return Omega_r / std::pow(a, 4) / E2(a); }
double Cosmology::omegaDarkEnergy(double a) const {
    return p.Omega_Lambda * darkEnergyScaling(a) / E2(a);
}
double Cosmology::omegaCurvature(double a) const { return Omega_k / (a * a) / E2(a); }

double Cosmology::acceleration(double a) const {
    // a''/a = -(H0^2/2) Σ Omega_i(a) (1 + 3 w_i)  (la courbure n'y contribue pas)
    const double src = 2.0 * Omega_r / std::pow(a, 4) + p.Omega_m / std::pow(a, 3) +
                       p.Omega_Lambda * darkEnergyScaling(a) * (1.0 + 3.0 * w(a));
    return -0.5 * H0inv * H0inv * a * src;
}

bool Cosmology::hasBigBang(std::string* why) const {
    if (p.H0 <= 0.0) {
        if (why) *why = "H0 doit être strictement positif";
        return false;
    }
    if (p.Omega_m <= 0.0 && Omega_r <= 0.0 && Omega_k <= 0.0) {
        if (why) *why = "sans matière, rayonnement ni courbure négative, il n'y a pas de Big Bang";
        return false;
    }
    for (int i = 0; i <= 800; ++i) {
        const double a = std::pow(10.0, -8.0 + 8.0 * i / 800.0);
        if (E2(a) <= 0.0) {
            if (why) *why = "H^2 < 0 dans le passé : l'univers aurait « rebondi » (pas de Big Bang)";
            return false;
        }
    }
    return true;
}

// ---------------------------------------------------------------------------
//  Temps
// ---------------------------------------------------------------------------
double Cosmology::hubbleTime() const { return 1.0 / H0inv; }

double Cosmology::age(double a) const {
    // t(a) = ∫0^a da'/(a' H(a')) ; changement de variable a' = s^2 pour lisser l'intégrande
    auto f = [this](double s) {
        s = std::max(s, kTiny);
        return 2.0 / (s * E(s * s));
    };
    return integrate(f, 0.0, std::sqrt(a)) / H0inv;
}

double Cosmology::lookbackTime(double z) const { return ageToday() - ageAtRedshift(z); }

double Cosmology::scaleFactor(double time) const {
    if (time <= 0.0) return 0.0;
    return interpolateScaleFactor(evolve(time), time);
}

double Cosmology::hubbleParameter(double time) const {
    const ExpansionHistory h = evolve(time);
    return h.states.empty() ? 0.0 : h.states.back().H;
}

std::vector<double> Cosmology::simulateExpansion(double startTime, double endTime,
                                                 double timeStep) const {
    std::vector<double> scaleFactors;
    const ExpansionHistory h = evolve(endTime);
    const int n = static_cast<int>(std::floor((endTime - startTime) / timeStep + 1e-9));
    for (int i = 0; i <= n; ++i)
        scaleFactors.push_back(interpolateScaleFactor(h, startTime + i * timeStep));
    return scaleFactors;
}

ExpansionHistory Cosmology::evolve(double tMax) const {
    ExpansionHistory hist;
    hist.fate = Fate::EternalAcceleration;
    if (!hasBigBang()) return hist;

    const double aInit = 1e-5;
    double t = age(aInit);
    double a = aInit;
    double v = a * H0inv * E(a);  // da/dt

    auto record = [&] {
        ExpansionState s;
        s.t = t; s.a = a; s.adot = v;
        s.H = (v / a) / kHubbleToGyr;
        s.q = (std::fabs(v) > 1e-300) ? -a * acceleration(a) / (v * v)
                                      : std::numeric_limits<double>::infinity();
        hist.states.push_back(s);
    };
    // Point initial a = 0 à t = 0 pour les tracés
    hist.states.push_back({0.0, 0.0, std::numeric_limits<double>::infinity(),
                           std::numeric_limits<double>::infinity(), deceleration(aInit)});
    record();

    const double dtMax = 0.02;   // Gyr
    const double eta = 2e-3;     // fraction du temps d'expansion a/|a'| par pas
    while (t < tMax) {
        double dt = std::min(dtMax, eta * a / std::max(std::fabs(v), 1e-300));
        dt = std::min(dt, tMax - t);
        // Runge-Kutta d'ordre 4 sur (a, a')
        const double k1a = v,                  k1v = acceleration(a);
        const double k2a = v + 0.5 * dt * k1v, k2v = acceleration(std::max(a + 0.5 * dt * k1a, kTiny));
        const double k3a = v + 0.5 * dt * k2v, k3v = acceleration(std::max(a + 0.5 * dt * k2a, kTiny));
        const double k4a = v + dt * k3v,       k4v = acceleration(std::max(a + dt * k3a, kTiny));
        a += dt / 6.0 * (k1a + 2.0 * k2a + 2.0 * k3a + k4a);
        v += dt / 6.0 * (k1v + 2.0 * k2v + 2.0 * k3v + k4v);
        t += dt;

        if (a <= aInit) {  // retour à une densité infinie
            a = std::max(a, 0.0);
            record();
            hist.fate = Fate::BigCrunch;
            hist.fateTime = t;
            return hist;
        }
        record();
        if (a > 1e8 && w(a) < -1.0) {  // expansion super-exponentielle : singularité imminente
            hist.fate = Fate::BigRip;
            hist.fateTime = t;
            return hist;
        }
        if (a > 1e12) break;
    }
    const ExpansionState& last = hist.states.back();
    if (last.adot < 0.0) {
        hist.fate = Fate::BigCrunch;  // en contraction : l'effondrement viendra plus tard
    } else if (w(last.a) < -1.0 && p.Omega_Lambda > 0.0) {
        hist.fate = Fate::BigRip;
    } else {
        hist.fate = std::fabs(last.q) < 1e-6 ? Fate::EternalCoasting
                  : last.q < 0.0            ? Fate::EternalAcceleration
                                            : Fate::EternalDeceleration;
    }
    return hist;
}

double interpolateScaleFactor(const ExpansionHistory& h, double t) {
    const auto& s = h.states;
    if (s.empty()) return 0.0;
    if (t <= s.front().t) return s.front().a;
    if (t >= s.back().t) return s.back().a;
    auto it = std::lower_bound(s.begin(), s.end(), t,
                               [](const ExpansionState& st, double tt) { return st.t < tt; });
    const ExpansionState& hi = *it;
    const ExpansionState& lo = *(it - 1);
    const double u = (t - lo.t) / (hi.t - lo.t);
    return lo.a + u * (hi.a - lo.a);
}

// ---------------------------------------------------------------------------
//  Distances
// ---------------------------------------------------------------------------
double Cosmology::comovingDistance(double z) const {
    auto f = [this](double zz) { return 1.0 / E(1.0 / (1.0 + zz)); };
    return kSpeedOfLight / p.H0 * integrate(f, 0.0, z, 1e-10);
}

double Cosmology::transverseComovingDistance(double z) const {
    const double DH = kSpeedOfLight / p.H0;
    const double DC = comovingDistance(z);
    if (std::fabs(Omega_k) < 1e-8) return DC;
    const double sk = std::sqrt(std::fabs(Omega_k));
    if (Omega_k > 0.0) return DH / sk * std::sinh(sk * DC / DH);
    return DH / sk * std::sin(sk * DC / DH);
}

double Cosmology::luminosityDistance(double z) const {
    return (1.0 + z) * transverseComovingDistance(z);
}

double Cosmology::angularDiameterDistance(double z) const {
    return transverseComovingDistance(z) / (1.0 + z);
}

double Cosmology::distanceModulus(double z) const {
    return 5.0 * std::log10(luminosityDistance(z) * 1e6 / 10.0);
}

double Cosmology::particleHorizon(double a) const {
    // chi = c ∫0^t dt/a = c/H0 ∫0^a da/(a^2 E) ; a' = s^2
    auto f = [this](double s) {
        s = std::max(s, 1e-8);
        return 2.0 / (s * s * s * E(s * s));
    };
    return kSpeedOfLight / p.H0 * integrate(f, 0.0, std::sqrt(a), 1e-10);
}

// ---------------------------------------------------------------------------
//  Croissance linéaire des structures
//  D'' + (1 - q) D' = (3/2) Omega_m(a) D, dérivées par rapport à ln a
// ---------------------------------------------------------------------------
void Cosmology::computeGrowth() {
    const double aStart = 1e-5;
    const int n = 4000;
    const double x0 = std::log(aStart), x1 = 0.0, h = (x1 - x0) / n;
    // Solution exacte matière + rayonnement (Meszaros) : D ∝ a + 2 a_eq / 3
    const double aEq = p.Omega_m > 0.0 ? Omega_r / p.Omega_m : 0.0;
    double D = aStart + 2.0 * aEq / 3.0, dD = aStart;

    auto acc = [this](double x, double y, double dy) {
        const double a = std::exp(x);
        return -(1.0 - deceleration(a)) * dy + 1.5 * omegaMatter(a) * y;
    };
    growthLnA.assign(1, x0); growthD.assign(1, D); growthF.assign(1, dD / D);
    for (int i = 0; i < n; ++i) {
        const double x = x0 + i * h;
        const double k1D = dD,                 k1v = acc(x, D, dD);
        const double k2D = dD + 0.5 * h * k1v, k2v = acc(x + 0.5 * h, D + 0.5 * h * k1D, k2D);
        const double k3D = dD + 0.5 * h * k2v, k3v = acc(x + 0.5 * h, D + 0.5 * h * k2D, k3D);
        const double k4D = dD + h * k3v,       k4v = acc(x + h, D + h * k3D, k4D);
        D  += h / 6.0 * (k1D + 2.0 * k2D + 2.0 * k3D + k4D);
        dD += h / 6.0 * (k1v + 2.0 * k2v + 2.0 * k3v + k4v);
        growthLnA.push_back(x + h);
        growthD.push_back(D);
        growthF.push_back(dD / D);
    }
    const double D1 = growthD.back();
    for (double& d : growthD) d /= D1;
}

double Cosmology::interpolateGrowth(const std::vector<double>& table, double a) const {
    if (growthLnA.empty()) return 0.0;
    const double x = std::log(std::max(a, kTiny));
    if (x <= growthLnA.front()) return table.front();
    if (x >= growthLnA.back()) return table.back();
    const double step = growthLnA[1] - growthLnA[0];
    const size_t i = static_cast<size_t>((x - growthLnA.front()) / step);
    const size_t j = std::min(i + 1, growthLnA.size() - 1);
    const double u = (x - growthLnA[i]) / step;
    return table[i] + u * (table[j] - table[i]);
}

double Cosmology::growthFactor(double a) const { return interpolateGrowth(growthD, a); }
double Cosmology::growthRate(double a) const { return interpolateGrowth(growthF, a); }

double Cosmology::fSigma8(double z) const {
    const double a = 1.0 / (1.0 + z);
    return growthRate(a) * p.sigma8 * growthFactor(a);
}

// ---------------------------------------------------------------------------
//  Époques clés
// ---------------------------------------------------------------------------
double Cosmology::redshiftMatterRadiationEquality() const {
    if (Omega_r <= 0.0 || p.Omega_m <= 0.0) return -1.0;
    return p.Omega_m / Omega_r - 1.0;
}

double Cosmology::redshiftMatterDarkEnergyEquality() const {
    if (p.Omega_Lambda <= 0.0 || p.Omega_m <= 0.0) return -1.0;
    auto g = [this](double lna) {
        const double a = std::exp(lna);
        return p.Omega_m / std::pow(a, 3) - p.Omega_Lambda * darkEnergyScaling(a);
    };
    // Balayage de a = 100 (futur) vers le passé : première égalité rencontrée
    const int n = 600;
    double hi = std::log(100.0), ghi = g(hi);
    for (int i = 1; i <= n; ++i) {
        const double lo = std::log(100.0) - i * (std::log(100.0) - std::log(1e-4)) / n;
        const double glo = g(lo);
        if ((glo > 0) != (ghi > 0)) return 1.0 / std::exp(bisect(g, lo, hi)) - 1.0;
        hi = lo; ghi = glo;
    }
    return -1.0;
}

double Cosmology::redshiftAccelerationOnset() const {
    if (deceleration(1.0) >= 0.0) return -1.0;
    auto g = [this](double lna) { return deceleration(std::exp(lna)); };
    const int n = 600;
    double hi = 0.0, ghi = g(hi);
    for (int i = 1; i <= n; ++i) {
        const double lo = -i * (-std::log(1e-4)) / n;
        const double glo = g(lo);
        if ((glo > 0) != (ghi > 0)) return 1.0 / std::exp(bisect(g, lo, hi)) - 1.0;
        hi = lo; ghi = glo;
    }
    return -1.0;
}

double Cosmology::criticalDensity() const {
    const double H = p.H0 * 1000.0 / kMpcInMeters;  // s^-1
    return 3.0 * H * H / (8.0 * kPi * kGravitation);
}

} // namespace cosmo
