#ifndef SVG_PLOT_H
#define SVG_PLOT_H

#include <string>
#include <vector>

// Petit moteur de graphiques SVG sans dépendance : courbes, axes linéaires ou
// logarithmiques, légende, repères verticaux. Thème sombre « ciel profond ».
namespace cosmo {

struct Series {
    std::string label;
    std::vector<double> x, y;
    std::string color;
    bool dashed = false;
    double width = 2.2;
};

struct Marker {
    double x;
    std::string label;
    std::string color = "#9aa4c7";
};

class SvgPlot {
public:
    SvgPlot(std::string title, std::string xLabel, std::string yLabel);

    void setSubtitle(const std::string& s) { subtitle = s; }
    void setLogX(bool v) { logX = v; }
    void setLogY(bool v) { logY = v; }
    void setXRange(double lo, double hi) { xMin = lo; xMax = hi; fixedX = true; }
    void setYRange(double lo, double hi) { yMin = lo; yMax = hi; fixedY = true; }
    void setLegendRight(bool v) { legendFx = v ? -1.0 : 0.0; legendFy = 0.0; }
    // Position du coin supérieur gauche de la légende, en fraction de la zone de tracé
    // (fx < 0 : collée à droite ; fy < 0 : collée en bas)
    void setLegendAt(double fx, double fy) { legendFx = fx; legendFy = fy; }
    void addSeries(const Series& s) { series.push_back(s); }
    void addVerticalMarker(const Marker& m) { vMarkers.push_back(m); }
    void addHorizontalLine(double y, const std::string& color = "#5b6488") { hLines.push_back({y, "", color}); }

    std::string render() const;
    bool save(const std::string& path) const;

    // Palette utilisée pour les modèles cosmologiques
    static const std::vector<std::string>& palette();

private:
    std::string title, subtitle, xLabel, yLabel;
    std::vector<Series> series;
    std::vector<Marker> vMarkers, hLines;
    bool logX = false, logY = false, fixedX = false, fixedY = false;
    double legendFx = 0.0, legendFy = 0.0;
    double xMin = 0, xMax = 1, yMin = 0, yMax = 1;
};

std::string xmlEscape(const std::string& s);
std::string formatNumber(double v);

} // namespace cosmo

#endif // SVG_PLOT_H
