// Tests de validation : on compare la simulation à des solutions analytiques connues.
#include "../Cosmology.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <string>

using namespace cosmo;

static int failures = 0;

static void check(const std::string& what, double got, double expected, double relTol) {
    const double err = std::fabs(got - expected) / std::max(std::fabs(expected), 1e-300);
    const bool ok = err <= relTol;
    std::printf("[%s] %-58s obtenu %-12.6g attendu %-12.6g (écart %.2e)\n",
                ok ? " OK " : "FAIL", what.c_str(), got, expected, err);
    if (!ok) ++failures;
}

static CosmoParams preset(const char* name) {
    CosmoParams p;
    findPreset(name, p);
    return p;
}

int main() {
    // Einstein-de Sitter : a = (t/t0)^(2/3), t0 = 2/(3 H0), q = 1/2, D = a
    {
        Cosmology c(preset("eds"));
        const double t0 = 2.0 / (3.0 * hubbleToInverseGyr(c.H0()));
        check("EdS : âge = 2/(3 H0)", c.ageToday(), t0, 1e-6);
        check("EdS : a(t0/8) = 1/4", c.scaleFactor(t0 / 8.0), 0.25, 1e-4);
        check("EdS : q = 1/2", c.deceleration(0.3), 0.5, 1e-12);
        check("EdS : facteur de croissance D(a) = a", c.growthFactor(0.25), 0.25, 1e-5);
        check("EdS : taux de croissance f = 1", c.growthRate(0.5), 1.0, 1e-6);
        const double DH = kSpeedOfLight / c.H0();
        check("EdS : D_C(z=1) = 2 c/H0 (1 - 1/sqrt 2)", c.comovingDistance(1.0),
              2.0 * DH * (1.0 - 1.0 / std::sqrt(2.0)), 1e-8);
        check("EdS : horizon des particules = 2 c/H0", c.particleHorizon(1.0), 2.0 * DH, 1e-5);
    }
    // Milne : univers vide, a = H0 t, D_L = c/H0 z (1 + z/2)
    {
        Cosmology c(preset("milne"));
        check("Milne : âge = 1/H0", c.ageToday(), c.hubbleTime(), 1e-6);
        const double DH = kSpeedOfLight / c.H0();
        check("Milne : D_L(z=2) = c/H0 z(1+z/2)", c.luminosityDistance(2.0), DH * 2.0 * 2.0, 1e-8);
    }
    // Planck 2018 : âge publié 13.797 ± 0.023 Gyr, z_eq ≈ 3400, transition vers l'accélération z ≈ 0.63
    {
        Cosmology c(preset("planck2018"));
        check("Planck 2018 : âge de l'univers (Gyr)", c.ageToday(), 13.797, 2e-3);
        check("Planck 2018 : égalité matière-rayonnement z_eq", c.redshiftMatterRadiationEquality(), 3402, 1e-2);
        check("Planck 2018 : début de l'accélération", c.redshiftAccelerationOnset(),
              std::cbrt(2.0 * c.omegaDE() / c.omegaM()) - 1.0, 1e-2);
        check("Planck 2018 : somme des Omega = 1 (univers plat)",
              c.omegaM() + c.omegaR() + c.omegaDE(), 1.0, 1e-12);
        // La trajectoire intégrée doit donner a = 1 à l'âge calculé par intégrale
        check("Planck 2018 : a(t0) = 1 par intégration RK4", c.scaleFactor(c.ageToday()), 1.0, 1e-4);
        check("Planck 2018 : H(t0) = H0 par intégration RK4", c.hubbleParameter(c.ageToday()), c.H0(), 1e-4);
        // Solution analytique ΛCDM plat sans rayonnement : f(Omega_m) ≈ Omega_m^0.55
        check("Planck 2018 : taux de croissance f ≈ Omega_m^0.55", c.growthRate(1.0),
              std::pow(c.omegaMatter(1.0), 0.55), 1e-2);
    }
    // Univers fermé sans Λ : cycloïde, t_crunch = π Omega_m / (H0 (Omega_m - 1)^(3/2))
    {
        Cosmology c(preset("closed"));
        const double Om = c.omegaM();
        const double tc = kPi * Om / (hubbleToInverseGyr(c.H0()) * std::pow(Om - 1.0, 1.5));
        const ExpansionHistory h = c.evolve(200.0);
        check("Fermé : Big Crunch détecté", h.fate == Fate::BigCrunch ? 1.0 : 0.0, 1.0, 0.0);
        check("Fermé : instant du Big Crunch (cycloïde)", h.fateTime, tc, 2e-3);
    }
    // Énergie fantôme : Big Rip ≈ t0 + 2 / (3 |1+w| H0 sqrt(Omega_DE)) (approximation tardive)
    {
        Cosmology c(preset("phantom"));
        const ExpansionHistory h = c.evolve(200.0);
        const double approx = 2.0 / (3.0 * 0.5 * hubbleToInverseGyr(c.H0()) * std::sqrt(c.omegaDE()));
        check("Fantôme : Big Rip détecté", h.fate == Fate::BigRip ? 1.0 : 0.0, 1.0, 0.0);
        check("Fantôme : délai avant le Big Rip (Gyr)", h.fateTime - c.ageToday(), approx, 5e-2);
    }

    std::printf("\n%s (%d échec(s))\n", failures ? "ÉCHEC" : "TOUS LES TESTS PASSENT", failures);
    return failures ? 1 : 0;
}
