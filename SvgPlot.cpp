#include "SvgPlot.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <limits>
#include <sstream>

namespace cosmo {

namespace {

constexpr double kWidth = 920, kHeight = 560;
constexpr double kLeft = 86, kRight = 30, kTop = 80, kBottom = 64;

const char* kBackground = "#0b1026";
const char* kPanel = "#11183a";
const char* kGrid = "#1f2a55";
const char* kText = "#e6e9f5";
const char* kMuted = "#9aa4c7";

// Pas « rond » (1, 2, 5 x 10^n) pour les graduations linéaires
double niceStep(double range, int targetTicks) {
    const double raw = range / std::max(targetTicks, 1);
    const double mag = std::pow(10.0, std::floor(std::log10(raw)));
    const double r = raw / mag;
    const double nice = r < 1.5 ? 1 : r < 3 ? 2 : r < 7 ? 5 : 10;
    return nice * mag;
}

std::vector<double> linearTicks(double lo, double hi) {
    std::vector<double> t;
    const double step = niceStep(hi - lo, 6);
    for (double v = std::ceil(lo / step - 1e-9) * step; v <= hi + 1e-9 * step; v += step)
        t.push_back(std::fabs(v) < 1e-12 * step ? 0.0 : v);
    return t;
}

std::vector<double> logTicks(double lo, double hi) {
    std::vector<double> t;
    const int d0 = static_cast<int>(std::floor(std::log10(lo) + 1e-9));
    const int d1 = static_cast<int>(std::ceil(std::log10(hi) - 1e-9));
    const int decades = d1 - d0;
    const int stride = decades > 12 ? 2 : 1;
    for (int d = d0; d <= d1; d += stride) {
        const double base = std::pow(10.0, d);
        if (decades <= 2) {
            for (double m : {1.0, 2.0, 5.0})
                if (base * m >= lo * (1 - 1e-9) && base * m <= hi * (1 + 1e-9)) t.push_back(base * m);
        } else if (base >= lo * (1 - 1e-9) && base <= hi * (1 + 1e-9)) {
            t.push_back(base);
        }
    }
    return t;
}

// Exposant en caractères Unicode : 10^-6 -> 10⁻⁶
std::string superscript(int e) {
    static const char* digits[] = {"⁰", "¹", "²", "³", "⁴", "⁵", "⁶", "⁷", "⁸", "⁹"};
    std::string out = e < 0 ? "⁻" : "";
    for (char c : std::to_string(std::abs(e))) out += digits[c - '0'];
    return out;
}

// Nombre de caractères affichés d'une chaîne UTF-8
size_t displayLength(const std::string& s) {
    size_t n = 0;
    for (unsigned char c : s) n += (c & 0xC0) != 0x80;
    return n;
}

std::string fmt(double v) {
    std::ostringstream os;
    os.setf(std::ios::fixed);
    os.precision(1);
    os << v;
    return os.str();
}

} // namespace

std::string xmlEscape(const std::string& s) {
    std::string out;
    for (char c : s) {
        switch (c) {
            case '&': out += "&amp;"; break;
            case '<': out += "&lt;"; break;
            case '>': out += "&gt;"; break;
            case '"': out += "&quot;"; break;
            default: out += c;
        }
    }
    return out;
}

std::string formatNumber(double v) {
    if (v == 0.0) return "0";
    const double av = std::fabs(v);
    std::ostringstream os;
    if (av >= 1e5 || av < 1e-3) {
        const int e = static_cast<int>(std::floor(std::log10(av) + 1e-9));
        const double m = v / std::pow(10.0, e);
        if (std::fabs(m - std::round(m)) < 1e-6 && std::fabs(std::round(m)) == 1.0)
            os << (m < 0 ? "-" : "") << "10" << superscript(e);
        else
            os << std::round(m * 10) / 10 << "×10" << superscript(e);
        return os.str();
    }
    os.precision(6);
    os << v;
    return os.str();
}

const std::vector<std::string>& SvgPlot::palette() {
    static const std::vector<std::string> colors = {
        "#4fc3f7", "#ffb74d", "#81c784", "#f06292", "#ba68c8", "#fff176", "#4db6ac", "#e0e0e0"};
    return colors;
}

SvgPlot::SvgPlot(std::string t, std::string xl, std::string yl)
    : title(std::move(t)), xLabel(std::move(xl)), yLabel(std::move(yl)) {}

std::string SvgPlot::render() const {
    // --- Détermination des bornes --------------------------------------------
    double x0 = xMin, x1 = xMax, y0 = yMin, y1 = yMax;
    if (!fixedX || !fixedY) {
        double ax0 = std::numeric_limits<double>::infinity(), ax1 = -ax0, ay0 = ax0, ay1 = -ax0;
        for (const auto& s : series) {
            for (size_t i = 0; i < s.x.size() && i < s.y.size(); ++i) {
                const double x = s.x[i], y = s.y[i];
                if (!std::isfinite(x) || !std::isfinite(y)) continue;
                if ((logX && x <= 0) || (logY && y <= 0)) continue;
                if (fixedX && (x < xMin || x > xMax)) continue;
                ax0 = std::min(ax0, x); ax1 = std::max(ax1, x);
                ay0 = std::min(ay0, y); ay1 = std::max(ay1, y);
            }
        }
        if (!std::isfinite(ax0)) { ax0 = 0; ax1 = 1; ay0 = 0; ay1 = 1; }
        if (!fixedX) { x0 = ax0; x1 = ax1; }
        if (!fixedY) {
            if (logY) { y0 = ay0 / 1.5; y1 = ay1 * 1.5; }
            else {
                const double pad = 0.06 * std::max(ay1 - ay0, 1e-12);
                y0 = ay0 - pad; y1 = ay1 + pad;
            }
        }
    }
    if (x1 <= x0) x1 = x0 + 1;
    if (y1 <= y0) y1 = y0 + 1;

    const double pw = kWidth - kLeft - kRight, ph = kHeight - kTop - kBottom;
    auto tx = [&](double x) {
        const double u = logX ? (std::log10(x) - std::log10(x0)) / (std::log10(x1) - std::log10(x0))
                              : (x - x0) / (x1 - x0);
        return kLeft + std::clamp(u, -20.0, 20.0) * pw;
    };
    auto ty = [&](double y) {
        const double u = logY ? (std::log10(y) - std::log10(y0)) / (std::log10(y1) - std::log10(y0))
                              : (y - y0) / (y1 - y0);
        return kTop + ph - std::clamp(u, -20.0, 20.0) * ph;
    };

    std::ostringstream o;
    o << "<svg xmlns=\"http://www.w3.org/2000/svg\" viewBox=\"0 0 " << kWidth << ' ' << kHeight
      << "\" width=\"" << kWidth << "\" height=\"" << kHeight
      << "\" font-family=\"system-ui, -apple-system, 'Segoe UI', Roboto, sans-serif\">\n";
    o << "<rect width=\"100%\" height=\"100%\" fill=\"" << kBackground << "\"/>\n";
    // Quelques étoiles décoratives (positions pseudo-aléatoires déterministes)
    unsigned seed = 12345;
    for (int i = 0; i < 70; ++i) {
        seed = seed * 1103515245u + 12345u; const double sx = (seed >> 8) % 920;
        seed = seed * 1103515245u + 12345u; const double sy = (seed >> 8) % 560;
        seed = seed * 1103515245u + 12345u; const double r = 0.4 + ((seed >> 8) % 10) / 12.0;
        o << "<circle cx=\"" << sx << "\" cy=\"" << sy << "\" r=\"" << fmt(r)
          << "\" fill=\"#ffffff\" opacity=\"0.25\"/>";
    }
    o << "\n<text x=\"" << kLeft << "\" y=\"34\" fill=\"" << kText
      << "\" font-size=\"21\" font-weight=\"600\">" << xmlEscape(title) << "</text>\n";
    if (!subtitle.empty())
        o << "<text x=\"" << kLeft << "\" y=\"58\" fill=\"" << kMuted << "\" font-size=\"13.5\">"
          << xmlEscape(subtitle) << "</text>\n";

    o << "<rect x=\"" << kLeft << "\" y=\"" << kTop << "\" width=\"" << pw << "\" height=\"" << ph
      << "\" fill=\"" << kPanel << "\" opacity=\"0.85\" stroke=\"" << kGrid << "\"/>\n";
    o << "<clipPath id=\"plot\"><rect x=\"" << kLeft << "\" y=\"" << kTop << "\" width=\"" << pw
      << "\" height=\"" << ph << "\"/></clipPath>\n";

    // --- Grille et graduations -------------------------------------------------
    for (double v : logX ? logTicks(x0, x1) : linearTicks(x0, x1)) {
        const double px = tx(v);
        o << "<line x1=\"" << fmt(px) << "\" y1=\"" << kTop << "\" x2=\"" << fmt(px) << "\" y2=\""
          << kTop + ph << "\" stroke=\"" << kGrid << "\"/>";
        o << "<text x=\"" << fmt(px) << "\" y=\"" << kTop + ph + 20 << "\" fill=\"" << kMuted
          << "\" font-size=\"12.5\" text-anchor=\"middle\">" << formatNumber(v) << "</text>\n";
    }
    for (double v : logY ? logTicks(y0, y1) : linearTicks(y0, y1)) {
        const double py = ty(v);
        o << "<line x1=\"" << kLeft << "\" y1=\"" << fmt(py) << "\" x2=\"" << kLeft + pw << "\" y2=\""
          << fmt(py) << "\" stroke=\"" << kGrid << "\"/>";
        o << "<text x=\"" << kLeft - 8 << "\" y=\"" << fmt(py + 4) << "\" fill=\"" << kMuted
          << "\" font-size=\"12.5\" text-anchor=\"end\">" << formatNumber(v) << "</text>\n";
    }
    o << "<text x=\"" << kLeft + pw / 2 << "\" y=\"" << kHeight - 18 << "\" fill=\"" << kText
      << "\" font-size=\"14\" text-anchor=\"middle\">" << xmlEscape(xLabel) << "</text>\n";
    o << "<text transform=\"translate(22," << kTop + ph / 2 << ") rotate(-90)\" fill=\"" << kText
      << "\" font-size=\"14\" text-anchor=\"middle\">" << xmlEscape(yLabel) << "</text>\n";

    o << "<g clip-path=\"url(#plot)\">\n";
    for (const auto& h : hLines) {
        if (logY && h.x <= 0) continue;
        o << "<line x1=\"" << kLeft << "\" x2=\"" << kLeft + pw << "\" y1=\"" << fmt(ty(h.x))
          << "\" y2=\"" << fmt(ty(h.x)) << "\" stroke=\"" << h.color << "\" stroke-width=\"1.2\"/>\n";
    }
    for (const auto& m : vMarkers) {
        if ((logX && m.x <= 0) || m.x < std::min(x0, x1) || m.x > std::max(x0, x1)) continue;
        o << "<line x1=\"" << fmt(tx(m.x)) << "\" x2=\"" << fmt(tx(m.x)) << "\" y1=\"" << kTop
          << "\" y2=\"" << kTop + ph << "\" stroke=\"" << m.color
          << "\" stroke-width=\"1.3\" stroke-dasharray=\"5 5\"/>\n";
    }
    // --- Courbes ---------------------------------------------------------------
    for (const auto& s : series) {
        std::ostringstream path;
        bool penDown = false;
        for (size_t i = 0; i < s.x.size() && i < s.y.size(); ++i) {
            const double x = s.x[i], y = s.y[i];
            const bool valid = std::isfinite(x) && std::isfinite(y) && !(logX && x <= 0) && !(logY && y <= 0);
            if (!valid) { penDown = false; continue; }
            path << (penDown ? 'L' : 'M') << fmt(tx(x)) << ',' << fmt(ty(y)) << ' ';
            penDown = true;
        }
        o << "<path d=\"" << path.str() << "\" fill=\"none\" stroke=\"" << s.color
          << "\" stroke-width=\"" << s.width << "\" stroke-linejoin=\"round\" stroke-linecap=\"round\""
          << (s.dashed ? " stroke-dasharray=\"8 6\"" : "") << "/>\n";
    }
    o << "</g>\n";
    for (const auto& m : vMarkers) {
        if ((logX && m.x <= 0) || m.x < std::min(x0, x1) || m.x > std::max(x0, x1)) continue;
        o << "<text x=\"" << fmt(tx(m.x) + 5) << "\" y=\"" << kTop + 16 << "\" fill=\"" << m.color
          << "\" font-size=\"12\" stroke=\"" << kBackground << "\" stroke-width=\"4\" paint-order=\"stroke\" transform=\"rotate(90 " << fmt(tx(m.x) + 5) << ' ' << kTop + 16
          << ")\">" << xmlEscape(m.label) << "</text>\n";
    }

    // --- Légende ---------------------------------------------------------------
    std::vector<const Series*> labelled;
    for (const auto& s : series) if (!s.label.empty()) labelled.push_back(&s);
    if (!labelled.empty()) {
        size_t maxLen = 0;
        for (const auto* s : labelled) maxLen = std::max(maxLen, displayLength(s->label));
        const double lw = 56 + 7.3 * static_cast<double>(maxLen), lh = 12 + 21.0 * labelled.size();
        const double lx = legendFx < 0 ? kLeft + pw - lw - 12 : kLeft + std::max(12.0, legendFx * pw);
        const double ly = legendFy < 0 ? kTop + ph - lh - 12 : kTop + std::max(12.0, legendFy * ph);
        o << "<rect x=\"" << fmt(lx) << "\" y=\"" << ly << "\" width=\"" << fmt(lw) << "\" height=\""
          << fmt(lh) << "\" rx=\"8\" fill=\"" << kBackground << "\" opacity=\"0.82\" stroke=\"" << kGrid << "\"/>\n";
        for (size_t i = 0; i < labelled.size(); ++i) {
            const double yy = ly + 18 + 21.0 * i;
            o << "<line x1=\"" << fmt(lx + 10) << "\" x2=\"" << fmt(lx + 36) << "\" y1=\"" << fmt(yy - 4)
              << "\" y2=\"" << fmt(yy - 4) << "\" stroke=\"" << labelled[i]->color << "\" stroke-width=\"3\""
              << (labelled[i]->dashed ? " stroke-dasharray=\"6 4\"" : "") << "/>";
            o << "<text x=\"" << fmt(lx + 44) << "\" y=\"" << fmt(yy) << "\" fill=\"" << kText
              << "\" font-size=\"12.5\">" << xmlEscape(labelled[i]->label) << "</text>\n";
        }
    }
    o << "</svg>\n";
    return o.str();
}

bool SvgPlot::save(const std::string& path) const {
    std::ofstream f(path);
    if (!f) return false;
    f << render();
    return static_cast<bool>(f);
}

} // namespace cosmo
