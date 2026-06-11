#include <volt/plugin/plugin_entry.h>
#include <volt/steinhardt_service.h>

#include <sstream>
#include <string>
#include <vector>

using namespace Volt;
using namespace Volt::Plugin;
using S = SteinhardtService;

namespace {

// Parse the --qlist CSV ("4,6,8,10,12") into the orders l. Mirrors the legacy
// parseQList: split on ',', trim, drop negatives/malformed; fall back to the
// default set when nothing valid remains so compute() never sees an empty list.
std::vector<int> parseQList(const std::string& raw) {
    std::vector<int> result;
    std::stringstream ss(raw);
    std::string token;
    while (std::getline(ss, token, ',')) {
        const std::size_t start = token.find_first_not_of(" \t\r\n");
        const std::size_t end = token.find_last_not_of(" \t\r\n");
        if (start == std::string::npos) continue;
        token = token.substr(start, end - start + 1);
        try {
            const int l = std::stoi(token);
            if (l >= 0) result.push_back(l);
        } catch (...) {
            // ignore malformed entry
        }
    }
    if (result.empty()) result = {4, 6, 8, 10, 12};
    return result;
}

} // namespace

static const std::vector<OptionBinding<S>> bindings = {
    // Custom-lambda overload: setQList takes a vector<int>, not a primitive.
    opt<S>(CliOption{"--qlist", "string", "Comma-separated orders l to compute", "4,6,8,10,12"},
        [](S& s, const OptsMap& o) {
            s.setQList(parseQList(CLI::getString(o, "--qlist", "4,6,8,10,12")));
        }),
    opt("--nnn",          "Nearest neighbors (0 => use --cutoff)",     12,    &S::setNnn),
    opt("--cutoff",       "Cutoff radius (used when --nnn 0)",         0.0,   &S::setCutoff),
    opt("--wl",           "Also compute Wigner-3j invariant w_l",      false, &S::setWlFlag),
    opt("--wlHat",        "Also compute normalized invariant w_l_hat", false, &S::setWlHatFlag),
    opt("--components",   "Emit Q_l^m components for this l (-1 off)",  -1,    &S::setComponents),
    opt("--onlySelected", "Only classify selected atoms",              false, &S::setOnlySelected),
};

VOLT_SERVICE_PLUGIN("volt-steinhardt", "Steinhardt Order Parameters", S, bindings)
