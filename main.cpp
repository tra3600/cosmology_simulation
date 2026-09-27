// Simulation de l'expansion de l'univers (modèle ΛCDM et extensions)
//
// Exemples :
//   ./cosmology_simulation                         modèle Planck 2018, sortie console
//   ./cosmology_simulation --preset desi           énergie sombre dynamique (DESI)
//   ./cosmology_simulation --H0 73 --Om 0.3        modèle personnalisé (plat)
//   ./cosmology_simulation --compare               compare tous les modèles prédéfinis
//   ./cosmology_simulation --svg figures           génère les illustrations + rapport HTML
#include "Cosmology.h"
#include "Report.h"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>

using namespace cosmo;

namespace {

struct Options {
    CosmoParams params;
    double futureGyr = 50.0;
    double dt = 1.0;
    std::string csvPath, svgDir;
    bool compare = false, plot = true, table = true;
};

void printUsage() {
    std::cout <<
        "Usage : cosmology_simulation [options]\n\n"
        "Modèle :\n"
        "  --preset NOM     modèle prédéfini (voir --list), défaut : planck2018\n"
        "  --H0 X           constante de Hubble (km/s/Mpc)\n"
        "  --Om X           densité de matière Omega_m\n"
        "  --OL X           densité d'énergie sombre (sinon : univers plat, OL = 1 - Om - Or)\n"
        "  --Or X           densité de rayonnement (défaut : calculée depuis T_CMB et N_eff)\n"
        "  --w0 X --wa X    équation d'état de l'énergie sombre w(a) = w0 + wa (1 - a)\n"
        "  --sigma8 X       amplitude des fluctuations de matière\n\n"
        "Simulation et sorties :\n"
        "  --future GYR     durée simulée dans le futur (défaut 50)\n"
        "  --dt GYR         pas de la table a(t) (défaut 1)\n"
        "  --csv FICHIER    exporte a(t), H, q, Omega_i... en CSV\n"
        "  --svg DOSSIER    génère les illustrations SVG et un rapport HTML (index.html)\n"
        "  --compare        tableau comparatif de tous les modèles prédéfinis\n"
        "  --no-table       n'affiche pas la table a(t)\n"
        "  --no-plot        n'affiche pas le graphique ASCII\n"
        "  --list           liste les modèles prédéfinis\n"
        "  -h, --help       cette aide\n";
}

void listPresets() {
    std::cout << "Modèles prédéfinis :\n";
    for (const auto& p : presets())
        std::cout << "  " << std::left << std::setw(12) << p.name << p.description << '\n';
}

bool parseArgs(int argc, char** argv, Options& o) {
    auto number = [&](int& i, const std::string& flag, double& out) {
        if (i + 1 >= argc) { std::cerr << "Valeur manquante pour " << flag << '\n'; return false; }
        char* end = nullptr;
        out = std::strtod(argv[++i], &end);
        if (end == argv[i] || *end) { std::cerr << "Nombre invalide pour " << flag << " : " << argv[i] << '\n'; return false; }
        return true;
    };
    // Le preset est lu en premier pour que les autres options puissent le modifier
    for (int i = 1; i + 1 < argc; ++i) {
        if (std::string(argv[i]) == "--preset" && !findPreset(argv[i + 1], o.params)) {
            std::cerr << "Modèle inconnu : " << argv[i + 1] << "\n";
            listPresets();
            return false;
        }
    }
    bool custom = false;
    for (int i = 1; i < argc; ++i) {
        const std::string a = argv[i];
        double v = 0;
        if (a == "-h" || a == "--help") { printUsage(); std::exit(0); }
        else if (a == "--list") { listPresets(); std::exit(0); }
        else if (a == "--preset") ++i;
        else if (a == "--H0") { if (!number(i, a, v)) return false; o.params.H0 = v; custom = true; }
        else if (a == "--Om") { if (!number(i, a, v)) return false; o.params.Omega_m = v; custom = true; }
        else if (a == "--OL") { if (!number(i, a, v)) return false; o.params.Omega_Lambda = v; o.params.flat = false; custom = true; }
        else if (a == "--Or") { if (!number(i, a, v)) return false; o.params.Omega_r = v; custom = true; }
        else if (a == "--w0") { if (!number(i, a, v)) return false; o.params.w0 = v; custom = true; }
        else if (a == "--wa") { if (!number(i, a, v)) return false; o.params.wa = v; custom = true; }
        else if (a == "--sigma8") { if (!number(i, a, v)) return false; o.params.sigma8 = v; }
        else if (a == "--future") { if (!number(i, a, o.futureGyr)) return false; }
        else if (a == "--dt") { if (!number(i, a, o.dt)) return false; }
        else if (a == "--csv" && i + 1 < argc) o.csvPath = argv[++i];
        else if (a == "--svg" && i + 1 < argc) o.svgDir = argv[++i];
        else if (a == "--compare") o.compare = true;
        else if (a == "--no-plot") o.plot = false;
        else if (a == "--no-table") o.table = false;
        else { std::cerr << "Option inconnue ou incomplète : " << a << "\n\n"; printUsage(); return false; }
    }
    if (custom) {
        o.params.name = "custom";
        o.params.description = "Modèle personnalisé (dérivé de " + o.params.description + ")";
    }
    if (o.dt <= 0 || o.futureGyr < 0) { std::cerr << "--dt doit être > 0 et --future >= 0\n"; return false; }
    return true;
}

std::string dominantComponent(const Cosmology& c, double a) {
    const double r = c.omegaRadiation(a), m = c.omegaMatter(a), d = c.omegaDarkEnergy(a),
                 k = c.omegaCurvature(a);
    if (r >= m && r >= d && r >= k) return "rayonnement";
    if (m >= d && m >= k) return "matière";
    if (d >= k) return "énergie sombre";
    return "courbure";
}

void printSummary(const Cosmology& c, const ExpansionHistory& h) {
    const CosmoParams& p = c.params();
    std::cout << std::fixed;
    std::cout << "══════════════════════════════════════════════════════════════════════\n"
              << "  Simulation cosmologique — modèle « " << p.name << " »\n  " << p.description << '\n'
              << "══════════════════════════════════════════════════════════════════════\n";
    std::cout << std::setprecision(4)
              << "  H0 = " << p.H0 << " km/s/Mpc   Ωm = " << c.omegaM() << "   Ωr = " << std::scientific
              << std::setprecision(3) << c.omegaR() << std::fixed << std::setprecision(4) << "   ΩDE = " << c.omegaDE()
              << "   Ωk = " << c.omegaK() << "\n  w(a) = " << p.w0 << " + " << p.wa << " (1 - a)   σ8 = " << p.sigma8
              << "   T_CMB = " << p.T_cmb << " K\n\n";

    std::cout << std::setprecision(3);
    std::cout << "  Âge de l'univers              : " << c.ageToday() << " milliards d'années (Gyr)\n"
              << "  Temps de Hubble 1/H0          : " << c.hubbleTime() << " Gyr\n"
              << "  Paramètre de décélération q0  : " << c.deceleration(1.0)
              << (c.deceleration(1.0) < 0 ? "  (l'expansion accélère)" : "  (l'expansion ralentit)") << '\n';
    const double zEq = c.redshiftMatterRadiationEquality(), zAcc = c.redshiftAccelerationOnset(),
                 zDE = c.redshiftMatterDarkEnergyEquality();
    if (zEq > 0) std::cout << "  Égalité matière-rayonnement   : z = " << std::setprecision(0) << zEq << std::setprecision(3)
                           << "  (t = " << c.ageAtRedshift(zEq) * 1e6 << " kyr)\n";
    if (zAcc > 0) std::cout << "  Début de l'accélération       : z = " << zAcc << "  (t = " << c.ageAtRedshift(zAcc)
                            << " Gyr, il y a " << c.lookbackTime(zAcc) << " Gyr)\n";
    if (zDE > 0) std::cout << "  Égalité matière-énergie sombre: z = " << zDE << "  (t = " << c.ageAtRedshift(zDE) << " Gyr)\n";
    else if (zDE > -1) std::cout << "  Égalité matière-énergie sombre: dans le futur, a = " << 1.0 / (1.0 + zDE) << '\n';
    std::cout << "  Rayon de l'univers observable : " << std::setprecision(1) << c.particleHorizon(1.0) * 3.26156e-3
              << " milliards d'années-lumière (" << c.particleHorizon(1.0) / 1000.0 << " Gpc)\n"
              << "  Densité critique              : " << std::scientific << std::setprecision(3) << c.criticalDensity()
              << " kg/m³ ≈ " << std::fixed << std::setprecision(2) << c.criticalDensity() / 1.6735575e-27
              << " atomes d'hydrogène par m³\n";
    std::cout << "  Destin de l'univers           : " << fateToString(h.fate);
    if (h.fateTime > 0)
        std::cout << " dans " << std::setprecision(1) << h.fateTime - c.ageToday() << " Gyr";
    std::cout << "\n\n";
}

void printTimeline(const Cosmology& c) {
    std::cout << "  Chronologie cosmique\n  "
              << std::left << std::setw(44) << "événement" << std::right << std::setw(12) << "z"
              << std::setw(14) << "âge" << std::setw(14) << "T_CMB (K)" << '\n';
    for (const auto& e : cosmicEvents(c)) {
        const double t = c.ageAtRedshift(e.z);
        std::ostringstream age;
        age << std::fixed;
        if (t < 1e-9) age << std::setprecision(0) << t * 3.15576e16 << " s";
        else if (t < 1e-3) age << std::setprecision(0) << t * 1e9 << " ans";
        else if (t < 1.0) age << std::setprecision(0) << t * 1e3 << " Myr";
        else age << std::setprecision(2) << t << " Gyr";
        // Les accents (UTF-8) occupent 2 octets : on compense pour aligner les colonnes
        int extra = 0;
        for (unsigned char ch : e.name) if ((ch & 0xC0) == 0x80) ++extra;
        std::cout << "  " << std::left << std::setw(44 + extra) << e.name << std::right << std::setw(12)
                  << std::setprecision(e.z >= 100 ? 0 : 2) << (e.z >= 1e5 ? std::scientific : std::fixed) << e.z
                  << std::fixed << std::setw(14) << age.str() << std::setw(14) << std::setprecision(1);
        if (c.cmbTemperature(e.z) >= 1e6) std::cout << std::scientific << c.cmbTemperature(e.z) << std::fixed;
        else std::cout << c.cmbTemperature(e.z);
        std::cout << '\n';
    }
    std::cout << '\n';
}

void printEvolutionTable(const Cosmology& c, const ExpansionHistory& h, double dt) {
    std::cout << "  Évolution temporelle (intégration RK4 de l'équation de Friedmann)\n"
              << "    t (Gyr)           a           z   H (km/s/Mpc)         q   composante dominante\n";
    const double t0 = c.ageToday(), tEnd = h.states.back().t;
    // Grille régulière + une ligne exacte pour aujourd'hui
    std::vector<double> times;
    for (double t = 0.0; t <= tEnd + 1e-9; t += dt) times.push_back(t);
    times.push_back(t0);
    std::sort(times.begin(), times.end());
    size_t k = 1;
    for (double t : times) {
        while (k + 1 < h.states.size() && h.states[k].t < t) ++k;
        const double a = interpolateScaleFactor(h, t);
        const bool now = t == t0;
        std::cout << (now ? "  ► " : "    ") << std::fixed << std::setprecision(2) << std::setw(7) << t;
        if (a <= 0 || t == 0.0) {
            std::cout << "  Big Bang\n";
            continue;
        }
        if (now) std::cout << std::setprecision(4) << std::setw(12) << 1.0 << std::setprecision(3) << std::setw(12) << 0.0;  // a = 1 et z = 0 par définition
        else std::cout << std::setprecision(4) << std::setw(12) << a << std::setw(12) << std::setprecision(3) << 1.0 / a - 1.0;
        const double aa = now ? 1.0 : a;
        const double H = (h.states[k].adot < 0 ? -1.0 : 1.0) * c.hubble(aa);
        std::cout << std::setw(15) << std::setprecision(2) << H << std::setw(10) << std::setprecision(3)
                  << c.deceleration(aa) << "   " << dominantComponent(c, aa) << (now ? "   ← aujourd'hui" : "") << '\n';
    }
    if (h.fateTime > 0)
        std::cout << "    " << std::setprecision(2) << std::setw(7) << h.fateTime << "  "
                  << (h.fate == Fate::BigRip ? "Big Rip : a → ∞" : "Big Crunch : a → 0") << '\n';
    std::cout << '\n';
}

void printDistances(const Cosmology& c) {
    std::cout << "  Distances (Mpc) et observables\n"
              << "        z    D_comobile    D_lumineuse   D_angulaire    μ (mag)   âge (Gyr)   regard en arrière (Gyr)\n";
    for (double z : {0.01, 0.1, 0.5, 1.0, 2.0, 3.0, 6.0, 10.0, 1089.9}) {
        std::cout << std::fixed << std::setw(9) << std::setprecision(z < 0.1 ? 2 : 1) << z << std::setprecision(1)
                  << std::setw(14) << c.transverseComovingDistance(z) << std::setw(15) << c.luminosityDistance(z)
                  << std::setw(14) << c.angularDiameterDistance(z) << std::setw(11) << std::setprecision(2)
                  << c.distanceModulus(z) << std::setw(12) << std::setprecision(3) << c.ageAtRedshift(z) << std::setw(14)
                  << c.lookbackTime(z) << '\n';
    }
    std::cout << '\n';
}

void printComparison(double futureGyr) {
    std::cout << "  Comparaison des modèles prédéfinis\n    " << std::left << std::setw(12) << "modèle" << std::right
              << std::setw(8) << "H0" << std::setw(8) << "Ωm" << std::setw(8) << "ΩDE" << std::setw(8) << "w0"
              << std::setw(11) << "âge (Gyr)" << std::setw(8) << "q0" << std::setw(8) << "z_acc" << "   destin\n";
    for (const auto& p : presets()) {
        Cosmology c(p);
        if (!c.hasBigBang()) continue;
        const ExpansionHistory h = c.evolve(c.ageToday() + futureGyr);
        const double zAcc = c.redshiftAccelerationOnset();
        std::cout << "    " << std::left << std::setw(12) << p.name << std::right << std::fixed << std::setprecision(2)
                  << std::setw(8) << p.H0 << std::setw(8) << c.omegaM() << std::setw(8) << c.omegaDE() << std::setw(8)
                  << p.w0 << std::setw(11) << c.ageToday() << std::setw(8) << c.deceleration(1.0) << std::setw(8);
        if (zAcc > 0) std::cout << zAcc; else std::cout << "—";
        std::cout << "   " << fateToString(h.fate);
        if (h.fateTime > 0) std::cout << " dans " << std::setprecision(1) << h.fateTime - c.ageToday() << " Gyr";
        std::cout << '\n';
    }
    std::cout << '\n';
}

} // namespace

int main(int argc, char** argv) {
    Options opt;
    if (!parseArgs(argc, argv, opt)) return 1;

    Cosmology cosmology(opt.params);
    std::string why;
    if (!cosmology.hasBigBang(&why)) {
        std::cerr << "Modèle non physique pour cette simulation : " << why << '\n';
        return 1;
    }
    const double t0 = cosmology.ageToday();
    const ExpansionHistory history = cosmology.evolve(t0 + opt.futureGyr);

    printSummary(cosmology, history);
    printTimeline(cosmology);
    if (opt.table) printEvolutionTable(cosmology, history, opt.dt);
    printDistances(cosmology);
    if (opt.plot) {
        std::cout << "  Facteur d'échelle a(t)\n" << asciiPlot(history, t0) << '\n';
    }
    if (opt.compare) printComparison(opt.futureGyr);

    if (!opt.csvPath.empty()) {
        if (writeCsv(cosmology, history, std::min(opt.dt, 0.1), opt.csvPath))
            std::cout << "  CSV écrit : " << opt.csvPath << '\n';
        else
            std::cerr << "  Impossible d'écrire " << opt.csvPath << '\n';
    }
    if (!opt.svgDir.empty()) {
        std::error_code ec;
        std::filesystem::create_directories(opt.svgDir, ec);
        std::vector<CosmoParams> models = presets();
        bool isPreset = false;
        for (const auto& p : models) isPreset = isPreset || p.name == opt.params.name;
        if (!isPreset) models.insert(models.begin(), opt.params);
        const auto figures = writeFigures(cosmology, models, opt.futureGyr, opt.svgDir);
        writeHtmlReport(cosmology, history, figures, opt.svgDir);
        std::cout << "  " << figures.size() << " illustrations SVG + rapport écrits dans " << opt.svgDir
                  << "/ (ouvrez " << opt.svgDir << "/index.html)\n";
    }
    return 0;
}
