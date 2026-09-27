#include "Report.h"
#include "SvgPlot.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <sstream>

namespace cosmo {

namespace {

constexpr double kSecondsPerGyr = 3.15576e16;
constexpr double kSolarSystemAgeGyr = 4.567;

// Âge lisible : secondes, minutes, années, millions ou milliards d'années
std::string formatAge(double gyr) {
    std::ostringstream o;
    o << std::setprecision(3);
    const double s = gyr * kSecondsPerGyr;
    if (s < 120) o << s << " s";
    else if (s < 7200) o << s / 60 << " min";
    else if (gyr < 1e-9) o << s / 86400 << " jours";
    else if (gyr < 1e-3) o << gyr * 1e9 << " ans";
    else if (gyr < 1.0) o << gyr * 1e3 << " Myr";
    else o << gyr << " Gyr";
    return o.str();
}

std::string formatRedshift(double z) {
    std::ostringstream o;
    if (z >= 1e5) o << std::scientific << std::setprecision(1) << z;
    else if (z >= 100) o << std::fixed << std::setprecision(0) << z;
    else o << std::fixed << std::setprecision(2) << z;
    return o.str();
}

// Downsampling d'une trajectoire pour des SVG légers
std::vector<size_t> strideIndices(size_t n, size_t maxPoints) {
    std::vector<size_t> idx;
    const size_t step = std::max<size_t>(1, n / maxPoints);
    for (size_t i = 0; i < n; i += step) idx.push_back(i);
    if (n && idx.back() != n - 1) idx.push_back(n - 1);
    return idx;
}

std::vector<double> linspace(double a, double b, int n) {
    std::vector<double> v(n);
    for (int i = 0; i < n; ++i) v[i] = a + (b - a) * i / (n - 1);
    return v;
}

std::vector<double> logspace(double a, double b, int n) {
    std::vector<double> v(n);
    for (int i = 0; i < n; ++i) v[i] = std::pow(10.0, std::log10(a) + (std::log10(b) - std::log10(a)) * i / (n - 1));
    return v;
}

struct Model {
    CosmoParams params;
    Cosmology cosmo;
    std::string color;
    bool highlighted;
};

// Frise chronologique verticale (temps en échelle logarithmique)
bool writeTimeline(const Cosmology& c, const std::string& path) {
    const std::vector<CosmicEvent> events = cosmicEvents(c);
    const double W = 920, H = 760, axisX = 170, top = 90, bottom = H - 40;
    const double lt0 = 0.0;                                           // 1 seconde
    const double lt1 = std::log10(c.ageToday() * kSecondsPerGyr) + 0.1;
    auto yOf = [&](double gyr) {
        const double lt = std::log10(std::max(gyr * kSecondsPerGyr, 1.0));
        return top + (lt - lt0) / (lt1 - lt0) * (bottom - top);
    };

    std::ostringstream o;
    o << "<svg xmlns=\"http://www.w3.org/2000/svg\" viewBox=\"0 0 " << W << ' ' << H << "\" width=\"" << W
      << "\" height=\"" << H << "\" font-family=\"system-ui, -apple-system, 'Segoe UI', Roboto, sans-serif\">\n"
      << "<defs><linearGradient id=\"heat\" x1=\"0\" y1=\"0\" x2=\"0\" y2=\"1\">"
      << "<stop offset=\"0\" stop-color=\"#fff7c2\"/><stop offset=\"0.25\" stop-color=\"#ffb74d\"/>"
      << "<stop offset=\"0.55\" stop-color=\"#e5533d\"/><stop offset=\"0.75\" stop-color=\"#7b3fa0\"/>"
      << "<stop offset=\"1\" stop-color=\"#1b2a6b\"/></linearGradient></defs>\n"
      << "<rect width=\"100%\" height=\"100%\" fill=\"#0b1026\"/>\n"
      << "<text x=\"40\" y=\"38\" fill=\"#e6e9f5\" font-size=\"21\" font-weight=\"600\">"
      << "Chronologie de l'univers — modèle " << xmlEscape(c.params().name) << "</text>\n"
      << "<text x=\"40\" y=\"62\" fill=\"#9aa4c7\" font-size=\"13.5\">Temps depuis le Big Bang en échelle "
         "logarithmique ; la couleur suit la température du rayonnement fossile, T = T0 (1+z)</text>\n";

    o << "<rect x=\"" << axisX - 7 << "\" y=\"" << top << "\" width=\"14\" height=\"" << bottom - top
      << "\" rx=\"7\" fill=\"url(#heat)\"/>\n";
    // Graduations : 1 s, 1 min, 1 an, 1000 ans, 1 Myr, 100 Myr, 1 Gyr, 10 Gyr
    const std::vector<std::pair<double, std::string>> ticks = {
        {1.0, "1 seconde"}, {60.0, "1 minute"}, {3600.0, "1 heure"}, {86400.0 * 365.25, "1 an"},
        {3.15576e10, "1000 ans"}, {3.15576e13, "1 Myr"}, {3.15576e15, "100 Myr"},
        {3.15576e16, "1 Gyr"}, {3.15576e17, "10 Gyr"}};
    for (const auto& [sec, label] : ticks) {
        const double y = yOf(sec / kSecondsPerGyr);
        o << "<line x1=\"" << axisX - 16 << "\" x2=\"" << axisX - 8 << "\" y1=\"" << y << "\" y2=\"" << y
          << "\" stroke=\"#9aa4c7\"/><text x=\"" << axisX - 22 << "\" y=\"" << y + 4
          << "\" fill=\"#9aa4c7\" font-size=\"12\" text-anchor=\"end\">" << label << "</text>\n";
    }

    // Placement des étiquettes sans chevauchement (espacement minimal vertical)
    const double gap = 46;
    std::vector<double> ly;
    for (const auto& e : events) ly.push_back(yOf(c.ageAtRedshift(e.z)));
    std::vector<double> target = ly;
    for (size_t i = 1; i < ly.size(); ++i) ly[i] = std::max(ly[i], ly[i - 1] + gap);
    if (!ly.empty() && ly.back() > bottom) {
        ly.back() = bottom;
        for (size_t i = ly.size() - 1; i-- > 0;) ly[i] = std::min(ly[i], ly[i + 1] - gap);
    }
    for (size_t i = 0; i < events.size(); ++i) {
        const auto& e = events[i];
        const double t = c.ageAtRedshift(e.z);
        const double labelX = 300;
        o << "<circle cx=\"" << axisX << "\" cy=\"" << target[i] << "\" r=\"4.5\" fill=\"#ffffff\"/>\n"
          << "<path d=\"M" << axisX + 6 << ',' << target[i] << " L" << axisX + 60 << ',' << target[i] << " L"
          << labelX - 12 << ',' << ly[i] << " L" << labelX - 4 << ',' << ly[i]
          << "\" fill=\"none\" stroke=\"#5b6488\" stroke-width=\"1.2\"/>\n"
          << "<text x=\"" << labelX << "\" y=\"" << ly[i] - 3 << "\" fill=\"#e6e9f5\" font-size=\"14.5\" "
          << "font-weight=\"600\">" << xmlEscape(e.name) << "</text>\n"
          << "<text x=\"" << labelX << "\" y=\"" << ly[i] + 15 << "\" fill=\"#9aa4c7\" font-size=\"12.5\">t = "
          << formatAge(t) << "  ·  z = " << formatRedshift(e.z) << "  ·  T = "
          << formatNumber(std::round(c.cmbTemperature(e.z) * 10) / 10) << " K  ·  " << xmlEscape(e.note)
          << "</text>\n";
    }
    o << "</svg>\n";
    std::ofstream f(path);
    f << o.str();
    return static_cast<bool>(f);
}

} // namespace

// ---------------------------------------------------------------------------
std::vector<CosmicEvent> cosmicEvents(const Cosmology& c) {
    std::vector<CosmicEvent> ev;
    const double T0 = c.params().T_cmb;
    ev.push_back({"Nucléosynthèse primordiale", "formation de l'hélium et du deutérium", 8e8 / T0 - 1.0});
    const double zEq = c.redshiftMatterRadiationEquality();
    if (zEq > 0) ev.push_back({"Égalité matière-rayonnement", "la matière commence à dominer", zEq});
    ev.push_back({"Recombinaison : fond diffus cosmologique", "l'univers devient transparent", 1089.9});
    ev.push_back({"Premières étoiles (aube cosmique)", "fin des âges sombres, estimation", 20.0});
    ev.push_back({"Galaxie JADES-GS-z14-0", "parmi les plus lointaines vues par JWST", 14.32});
    ev.push_back({"Milieu de la réionisation", "l'hydrogène intergalactique est ionisé", 7.7});
    ev.push_back({"« Midi cosmique »", "pic de formation des étoiles", 2.0});
    const double zAcc = c.redshiftAccelerationOnset();
    if (zAcc > 0) ev.push_back({"Début de l'accélération", "q devient négatif", zAcc});
    // Système solaire : il y a 4,567 Gyr
    const double tSolar = c.ageToday() - kSolarSystemAgeGyr;
    if (tSolar > 0) {
        double lo = 1e-4, hi = 1.0;
        for (int i = 0; i < 80; ++i) {
            const double mid = 0.5 * (lo + hi);
            (c.age(mid) < tSolar ? lo : hi) = mid;
        }
        ev.push_back({"Formation du Système solaire", "il y a 4,57 Gyr", 1.0 / lo - 1.0});
    }
    const double zDE = c.redshiftMatterDarkEnergyEquality();
    if (zDE > 0) ev.push_back({"Égalité matière-énergie sombre", "l'énergie sombre prend le dessus", zDE});
    ev.push_back({"Aujourd'hui", "", 0.0});
    std::stable_sort(ev.begin(), ev.end(), [](const CosmicEvent& a, const CosmicEvent& b) { return a.z > b.z; });
    return ev;
}

bool writeCsv(const Cosmology& c, const ExpansionHistory& h, double dt, const std::string& path) {
    std::ofstream f(path);
    if (!f) return false;
    f << "t_Gyr,a,z,H_km_s_Mpc,q,Omega_m,Omega_r,Omega_de,Omega_k,T_cmb_K,growth_D\n";
    f << std::setprecision(8);
    const double tEnd = h.states.empty() ? 0.0 : h.states.back().t;
    size_t k = 1;
    for (double t = dt; t <= tEnd + 1e-9; t += dt) {
        while (k + 1 < h.states.size() && h.states[k].t < t) ++k;
        const double a = interpolateScaleFactor(h, t);
        if (a <= 0) continue;
        const ExpansionState& s = h.states[k];
        const double sign = s.adot < 0 ? -1.0 : 1.0;
        f << t << ',' << a << ',' << 1.0 / a - 1.0 << ',' << sign * c.hubble(a) << ','
          << c.deceleration(a) << ',' << c.omegaMatter(a) << ',' << c.omegaRadiation(a) << ','
          << c.omegaDarkEnergy(a) << ',' << c.omegaCurvature(a) << ',' << c.cmbTemperature(1.0 / a - 1.0)
          << ',' << (a <= 1.0 ? c.growthFactor(a) : std::nan("")) << '\n';
    }
    return static_cast<bool>(f);
}

std::vector<std::string> writeFigures(const Cosmology& main, const std::vector<CosmoParams>& models,
                                      double futureGyr, const std::string& dir) {
    std::vector<std::string> files;
    std::vector<Model> ms;
    const auto& pal = SvgPlot::palette();
    for (size_t i = 0; i < models.size(); ++i) {
        Cosmology c(models[i]);
        if (!c.hasBigBang()) continue;
        const bool hl = models[i].name == main.params().name;
        ms.push_back({models[i], c, pal[i % pal.size()], hl});
    }
    auto save = [&](const SvgPlot& p, const std::string& name) {
        if (p.save(dir + "/" + name)) files.push_back(name);
    };

    // 1. Facteur d'échelle en fonction du temps (passé et futur) --------------
    {
        SvgPlot p("Expansion de l'univers : facteur d'échelle a(t)", "temps depuis aujourd'hui (milliards d'années)",
                  "facteur d'échelle a  (a = 1 aujourd'hui)");
        p.setSubtitle("Même point de départ aujourd'hui, destins différents : Big Freeze, Big Crunch ou Big Rip");
        p.setXRange(-16.0, futureGyr);
        p.setYRange(0.0, 5.0);
        for (const auto& m : ms) {
            const double t0 = m.cosmo.ageToday();
            const ExpansionHistory h = m.cosmo.evolve(t0 + futureGyr);
            Series s{m.params.name, {}, {}, m.color, false, m.highlighted ? 3.4 : 2.0};
            for (size_t i : strideIndices(h.states.size(), 1500)) {
                s.x.push_back(h.states[i].t - t0);
                s.y.push_back(h.states[i].a);
            }
            p.addSeries(s);
        }
        p.addVerticalMarker({0.0, "aujourd'hui"});
        p.addHorizontalLine(1.0);
        save(p, "expansion.svg");
    }
    // 2. Vitesse d'expansion a' = H/(1+z) : le minimum marque le début de l'accélération
    {
        SvgPlot p("Vitesse d'expansion : ȧ = H(z) / (1+z)", "redshift z", "ȧ  (km/s/Mpc)");
        p.setSubtitle("Une courbe qui remonte vers z = 0 signe une expansion qui accélère");
        p.setXRange(0.0, 4.0);
        p.setLegendRight(true);
        const auto zs = linspace(0.0, 4.0, 300);
        for (const auto& m : ms) {
            Series s{m.params.name, {}, {}, m.color, false, m.highlighted ? 3.4 : 2.0};
            for (double z : zs) { s.x.push_back(z); s.y.push_back(m.cosmo.hubble(1.0 / (1.0 + z)) / (1.0 + z)); }
            p.addSeries(s);
        }
        save(p, "expansion_rate.svg");
    }
    // 3. Paramètre de décélération q(z)
    {
        SvgPlot p("Paramètre de décélération q(z)", "redshift z", "q = -a ä / ȧ²");
        p.setSubtitle("q > 0 : l'expansion freine ; q < 0 : elle accélère (énergie sombre)");
        p.setXRange(0.0, 5.0);
        p.setLegendRight(true);
        const auto zs = linspace(0.0, 5.0, 300);
        for (const auto& m : ms) {
            Series s{m.params.name, {}, {}, m.color, false, m.highlighted ? 3.4 : 2.0};
            for (double z : zs) { s.x.push_back(z); s.y.push_back(m.cosmo.deceleration(1.0 / (1.0 + z))); }
            p.addSeries(s);
        }
        p.addHorizontalLine(0.0, "#e6e9f5");
        save(p, "deceleration.svg");
    }
    // 4. Composition de l'univers au cours du temps (modèle principal)
    {
        const Cosmology& c = main;
        SvgPlot p("Qui domine l'univers ? Fractions de densité Ω_i(a)", "facteur d'échelle a (log)", "Ω_i(a)");
        p.setSubtitle("Modèle " + c.params().name + " : ère du rayonnement → ère de la matière → ère de l'énergie sombre");
        p.setLogX(true);
        p.setXRange(1e-6, 100.0);
        p.setYRange(0.0, 1.05);
        const auto as = logspace(1e-6, 100.0, 500);
        Series r{"rayonnement", {}, {}, "#ffb74d"}, mm{"matière (noire + ordinaire)", {}, {}, "#4fc3f7"},
            de{"énergie sombre", {}, {}, "#ba68c8"}, k{"courbure", {}, {}, "#81c784", true};
        for (double a : as) {
            r.x.push_back(a); r.y.push_back(c.omegaRadiation(a));
            mm.x.push_back(a); mm.y.push_back(c.omegaMatter(a));
            de.x.push_back(a); de.y.push_back(c.omegaDarkEnergy(a));
            k.x.push_back(a); k.y.push_back(c.omegaCurvature(a));
        }
        p.addSeries(r); p.addSeries(mm); p.addSeries(de);
        if (std::fabs(c.omegaK()) > 1e-6) p.addSeries(k);
        const double zEq = c.redshiftMatterRadiationEquality();
        if (zEq > 0) p.addVerticalMarker({1.0 / (1.0 + zEq), "égalité matière-rayonnement"});
        p.addVerticalMarker({1.0 / 1090.9, "fond diffus cosmologique", "#ffb74d"});
        p.addVerticalMarker({1.0, "aujourd'hui", "#e6e9f5"});
        p.setLegendAt(0.43, 0.40);
        save(p, "densities.svg");
    }
    // 5. Distances cosmologiques (modèle principal)
    {
        const Cosmology& c = main;
        SvgPlot p("Les distances en cosmologie", "redshift z (log)", "distance (Mpc)");
        p.setSubtitle("D_A passe par un maximum : au-delà de z ≈ 1.6, les objets lointains paraissent plus grands !");
        p.setLogX(true);
        p.setLogY(true);
        p.setXRange(0.01, 1100.0);
        p.setYRange(20.0, 3e7);
        const auto zs = logspace(0.01, 1100.0, 220);
        Series dl{"distance lumineuse D_L (supernovae)", {}, {}, "#ffb74d"},
            dc{"distance comobile D_C", {}, {}, "#4fc3f7"},
            da{"distance angulaire D_A (taille apparente)", {}, {}, "#f06292"},
            lb{"c × temps de regard en arrière", {}, {}, "#81c784", true};
        for (double z : zs) {
            dl.x.push_back(z); dl.y.push_back(c.luminosityDistance(z));
            dc.x.push_back(z); dc.y.push_back(c.transverseComovingDistance(z));
            da.x.push_back(z); da.y.push_back(c.angularDiameterDistance(z));
            lb.x.push_back(z); lb.y.push_back(c.lookbackTime(z) * 306.601);  // Gyr → Mpc (c = 306.6 Mpc/Gyr)
        }
        p.addSeries(dl); p.addSeries(dc); p.addSeries(da); p.addSeries(lb);
        p.addHorizontalLine(c.particleHorizon(1.0), "#e6e9f5");
        save(p, "distances.svg");
    }
    // 6. Diagramme de Hubble différentiel (ce que mesurent les supernovae Ia)
    {
        const Cosmology ref(presets().front());
        SvgPlot p("Diagramme de Hubble : écart de module de distance", "redshift z",
                  "Δμ = μ(modèle) − μ(Planck 2018)  (magnitudes)");
        p.setSubtitle("Les supernovae Ia mesurent μ à ~0.1 mag près : c'est ainsi qu'on a découvert l'accélération (1998)");
        p.setXRange(0.0, 2.5);
        p.setYRange(-0.8, 0.8);
        p.setLegendRight(true);
        const auto zs = linspace(0.01, 2.5, 200);
        for (const auto& m : ms) {
            Series s{m.params.name, {}, {}, m.color, false, m.highlighted ? 3.4 : 2.0};
            for (double z : zs) { s.x.push_back(z); s.y.push_back(m.cosmo.distanceModulus(z) - ref.distanceModulus(z)); }
            p.addSeries(s);
        }
        p.addHorizontalLine(0.0, "#e6e9f5");
        save(p, "hubble_diagram.svg");
    }
    // 7. Croissance des structures f σ8(z)
    {
        SvgPlot p("Croissance des grandes structures : f σ8(z)", "redshift z", "f σ8");
        p.setSubtitle("Mesuré par les distorsions dans l'espace des redshifts (relevés BOSS, eBOSS, DESI)");
        p.setXRange(0.0, 2.0);
        p.setLegendRight(true);
        const auto zs = linspace(0.0, 2.0, 200);
        for (const auto& m : ms) {
            if (m.params.Omega_m <= 0) continue;
            Series s{m.params.name, {}, {}, m.color, false, m.highlighted ? 3.4 : 2.0};
            for (double z : zs) { s.x.push_back(z); s.y.push_back(m.cosmo.fSigma8(z)); }
            p.addSeries(s);
        }
        save(p, "growth.svg");
    }
    // 8. Âge de l'univers vs H0 : la tension de Hubble
    {
        SvgPlot p("Tension de Hubble : âge de l'univers selon H0", "H0 (km/s/Mpc)", "âge (milliards d'années)");
        p.setSubtitle("Planck (CMB) : 67.4 ± 0.5 · SH0ES (céphéides + SN Ia) : 73.0 ± 1.0 · les plus vieilles étoiles : ≳ 12-13 Gyr");
        p.setXRange(60.0, 80.0);
        const auto hs = linspace(60.0, 80.0, 81);
        Series lcdm{"ΛCDM plat, Ωm = 0.315", {}, {}, "#4fc3f7"}, eds{"matière seule (Einstein-de Sitter)", {}, {}, "#ba68c8"};
        for (double H : hs) {
            CosmoParams a = presets().front(); a.H0 = H;
            CosmoParams b; findPreset("eds", b); b.H0 = H;
            lcdm.x.push_back(H); lcdm.y.push_back(Cosmology(a).ageToday());
            eds.x.push_back(H); eds.y.push_back(Cosmology(b).ageToday());
        }
        p.addSeries(lcdm); p.addSeries(eds);
        p.addVerticalMarker({67.36, "Planck 2018", "#4fc3f7"});
        p.addVerticalMarker({73.04, "SH0ES 2022", "#ffb74d"});
        p.setLegendAt(0.0, -1.0);
        save(p, "hubble_tension.svg");
    }
    // 9. Frise chronologique
    if (writeTimeline(main, dir + "/timeline.svg")) files.push_back("timeline.svg");
    return files;
}

bool writeHtmlReport(const Cosmology& c, const ExpansionHistory& h, const std::vector<std::string>& figures,
                     const std::string& dir) {
    std::ofstream f(dir + "/index.html");
    if (!f) return false;
    const CosmoParams& p = c.params();
    f << std::setprecision(4);
    f << "<!doctype html>\n<html lang=\"fr\"><head><meta charset=\"utf-8\">"
         "<meta name=\"viewport\" content=\"width=device-width, initial-scale=1\">"
         "<title>Simulation cosmologique</title><style>"
         ":root{--bg:#0b1026;--panel:#11183a;--line:#1f2a55;--text:#e6e9f5;--muted:#9aa4c7;--accent:#4fc3f7}"
         "body{margin:0;background:var(--bg);color:var(--text);font:16px/1.55 system-ui,-apple-system,'Segoe UI',Roboto,sans-serif}"
         "main{max-width:960px;margin:0 auto;padding:32px 16px}h1{font-size:30px;margin:0 0 4px}"
         "h2{margin-top:40px;border-bottom:1px solid var(--line);padding-bottom:6px}.muted{color:var(--muted)}"
         ".cards{display:grid;grid-template-columns:repeat(auto-fit,minmax(170px,1fr));gap:12px;margin:20px 0}"
         ".card{background:var(--panel);border:1px solid var(--line);border-radius:12px;padding:14px}"
         ".card b{display:block;font-size:22px;color:var(--accent)}table{border-collapse:collapse;width:100%;font-size:14.5px}"
         "td,th{border-bottom:1px solid var(--line);padding:6px 8px;text-align:left}th{color:var(--muted);font-weight:500}"
         "img{width:100%;height:auto;border-radius:12px;border:1px solid var(--line);margin:12px 0}"
         ".scroll{overflow-x:auto}</style></head><body><main>\n";
    f << "<h1>Simulation de l'expansion de l'univers</h1><p class=\"muted\">Modèle <b>" << xmlEscape(p.name) << "</b> — "
      << xmlEscape(p.description) << "</p>\n";
    f << "<div class=\"cards\">"
      << "<div class=\"card\"><b>" << c.ageToday() << " Gyr</b>âge de l'univers</div>"
      << "<div class=\"card\"><b>" << c.hubbleTime() << " Gyr</b>temps de Hubble 1/H0</div>"
      << "<div class=\"card\"><b>" << c.particleHorizon(1.0) * 3.26156e-3 << " Gal</b>rayon de l'univers observable (milliards d'années-lumière)</div>"
      << "<div class=\"card\"><b>" << c.criticalDensity() / 1.6735575e-27 << "</b>atomes d'hydrogène par m³ (densité critique)</div>"
      << "<div class=\"card\"><b>" << xmlEscape(fateToString(h.fate)) << "</b>destin"
      << (h.fateTime > 0 ? " dans " + formatAge(h.fateTime - c.ageToday()) : std::string()) << "</div></div>\n";
    f << "<h2>Paramètres</h2><div class=\"scroll\"><table><tr><th>H0</th><th>Ωm</th><th>Ωr</th><th>ΩDE</th><th>Ωk</th>"
         "<th>w0</th><th>wa</th><th>σ8</th><th>T_CMB</th></tr><tr>"
      << "<td>" << p.H0 << " km/s/Mpc</td><td>" << c.omegaM() << "</td><td>" << c.omegaR() << "</td><td>" << c.omegaDE()
      << "</td><td>" << c.omegaK() << "</td><td>" << p.w0 << "</td><td>" << p.wa << "</td><td>" << p.sigma8 << "</td><td>"
      << p.T_cmb << " K</td></tr></table></div>\n";
    f << "<h2>Chronologie</h2><div class=\"scroll\"><table><tr><th>Événement</th><th>redshift z</th><th>âge</th>"
         "<th>il y a</th><th>T du CMB</th></tr>";
    for (const auto& e : cosmicEvents(c)) {
        f << "<tr><td>" << xmlEscape(e.name) << "</td><td>" << formatRedshift(e.z) << "</td><td>"
          << formatAge(c.ageAtRedshift(e.z)) << "</td><td>" << formatAge(std::max(c.lookbackTime(e.z), 0.0))
          << "</td><td>" << formatNumber(std::round(c.cmbTemperature(e.z) * 10) / 10) << " K</td></tr>";
    }
    f << "</table></div>\n<h2>Illustrations</h2>\n";
    for (const auto& fig : figures) f << "<img src=\"" << fig << "\" alt=\"" << fig << "\" loading=\"lazy\">\n";
    f << "<p class=\"muted\">Généré par cosmology_simulation.</p></main></body></html>\n";
    return static_cast<bool>(f);
}

std::string asciiPlot(const ExpansionHistory& h, double t0, int width, int height) {
    if (h.states.size() < 2) return "";
    const double tEnd = h.states.back().t;
    double aMax = 0;
    for (const auto& s : h.states) aMax = std::max(aMax, s.a);
    aMax = std::min(aMax, 5.0);
    std::vector<std::string> grid(height, std::string(width, ' '));
    for (int col = 0; col < width; ++col) {
        const double t = tEnd * col / (width - 1);
        const double a = interpolateScaleFactor(h, t);
        const int row = static_cast<int>(std::round((1.0 - std::min(a, aMax) / aMax) * (height - 1)));
        grid[row][col] = a > aMax ? '^' : '*';
    }
    const int todayCol = static_cast<int>(std::round(t0 / tEnd * (width - 1)));
    const int oneRow = static_cast<int>(std::round((1.0 - 1.0 / aMax) * (height - 1)));
    for (int r = 0; r < height; ++r)
        if (todayCol >= 0 && todayCol < width && grid[r][todayCol] == ' ') grid[r][todayCol] = '|';
    for (int col = 0; col < width; ++col)
        if (oneRow >= 0 && oneRow < height && grid[oneRow][col] == ' ') grid[oneRow][col] = '.';

    std::ostringstream o;
    o << std::fixed << std::setprecision(2);
    for (int r = 0; r < height; ++r) {
        const double a = aMax * (1.0 - static_cast<double>(r) / (height - 1));
        if (r % 4 == 0 || r == height - 1) o << std::setw(6) << a << " ┤";
        else o << "       │";
        o << grid[r] << '\n';
    }
    o << "       └" << std::string(width, '-') << '\n';
    o << "        0" << std::string(std::max(0, width / 2 - 10), ' ') << "temps (Gyr)"
      << std::string(std::max(0, width / 2 - 8), ' ') << std::setprecision(1) << tEnd << '\n';
    o << "        | = aujourd'hui   . = a = 1   ^ = hors échelle\n";
    return o.str();
}

} // namespace cosmo
