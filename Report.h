#ifndef REPORT_H
#define REPORT_H

#include "Cosmology.h"

#include <string>
#include <vector>

namespace cosmo {

// Un événement marquant de l'histoire cosmique, repéré par son redshift
struct CosmicEvent {
    std::string name;
    std::string note;
    double z;
};

// Chronologie de l'univers pour un modèle donné (triée du plus ancien au plus récent)
std::vector<CosmicEvent> cosmicEvents(const Cosmology& c);

// Table a(t), H, q, fractions de densité... au format CSV
bool writeCsv(const Cosmology& c, const ExpansionHistory& h, double dt, const std::string& path);

// Génère toutes les illustrations SVG + un rapport HTML dans `dir`.
// `models` sert aux figures comparatives ; `main` est le modèle mis en avant.
std::vector<std::string> writeFigures(const Cosmology& main, const std::vector<CosmoParams>& models,
                                      double futureGyr, const std::string& dir);
bool writeHtmlReport(const Cosmology& main, const ExpansionHistory& h,
                     const std::vector<std::string>& figures, const std::string& dir);

// Graphique ASCII de a(t) pour le terminal
std::string asciiPlot(const ExpansionHistory& h, double t0, int width = 72, int height = 20);

} // namespace cosmo

#endif // REPORT_H
