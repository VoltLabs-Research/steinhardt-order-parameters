#include <volt/steinhardt_engine.h>

#include <volt/analysis/cutoff_neighbor_finder.h>
#include <volt/analysis/nearest_neighbor_finder.h>

#include <oneapi/tbb/blocked_range.h>
#include <oneapi/tbb/enumerable_thread_specific.h>
#include <oneapi/tbb/parallel_for.h>
#include <spdlog/spdlog.h>

#include <algorithm>
#include <array>
#include <atomic>
#include <cfloat>
#include <cmath>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <vector>

namespace Volt{

namespace {

// LAMMPS `MY_EPSILON` (10 * DBL_EPSILON) and `QEPSILON` (1e-6) are required
// for bit-for-bit agreement with the reference implementation.
constexpr double MY_EPSILON = 10.0 * DBL_EPSILON;
constexpr double QEPSILON   = 1.0e-6;
constexpr double MY_4PI     = 4.0 * M_PI;

// q_6-based structure thresholds (see class docstring).
constexpr double Q6_LIQUID_MAX   = 0.30;
constexpr double Q6_INTERFACE_MAX = 0.45;

// Heuristic: among the requested orders pick the entry whose l is closest
// to 6; q_6 is the canonical crystallinity discriminator.
int pickReferenceIndex(const std::vector<int>& qlist){
    int best = -1;
    int bestDist = std::numeric_limits<int>::max();
    for(std::size_t i = 0; i < qlist.size(); ++i){
        const int d = std::abs(qlist[i] - 6);
        if(d < bestDist){
            bestDist = d;
            best = static_cast<int>(i);
        }
    }
    return best;
}

int classifyFromQ6(double q6){
    if(q6 < Q6_LIQUID_MAX)   return SteinhardtEngine::STRUCTURE_LIQUID;
    if(q6 < Q6_INTERFACE_MAX) return SteinhardtEngine::STRUCTURE_INTERFACE;
    return SteinhardtEngine::STRUCTURE_CRYSTAL;
}

} // namespace

// ---------------------------------------------------------------------------
// Math helpers ported verbatim from LAMMPS compute_orientorder_atom.cpp.
// ---------------------------------------------------------------------------

double SteinhardtEngine::factorial(int n){
    // LAMMPS ships a precomputed factorial table in MathSpecial; here we
    // compute it directly in `double`. Any n beyond ~170 overflows, but the
    // Wigner-3j table only ever calls factorial() with n <= 3*lmax + 1, and
    // `lmax` is bounded by the public `components <= 20` cap and
    // qlist entries in practice rarely exceed 12.
    double result = 1.0;
    for(int i = 2; i <= n; ++i){
        result *= static_cast<double>(i);
    }
    return result;
}

// sign convention: sign(Y_ll(0,0)) = (-1)^l
double SteinhardtEngine::associatedLegendre(int l, int m, double x){
    if(l < m) return 0.0;

    double p(1.0), pm1(0.0), pm2(0.0);

    if(m != 0){
        const double msqx = -std::sqrt(1.0 - x * x);
        for(int i = 1; i < m + 1; ++i){
            p *= static_cast<double>(2 * i - 1) * msqx;
        }
    }

    for(int i = m + 1; i < l + 1; ++i){
        pm2 = pm1;
        pm1 = p;
        p = (static_cast<double>(2 * i - 1) * x * pm1
            - static_cast<double>(i + m - 1) * pm2)
            / static_cast<double>(i - m);
    }

    return p;
}

// Y_l^m (theta, phi) = polar_prefactor(l, m, cos theta) * exp(i m phi)
double SteinhardtEngine::polarPrefactor(int l, int m, double costheta){
    const int mabs = std::abs(m);

    double prefactor = 1.0;
    for(int i = l - mabs + 1; i < l + mabs + 1; ++i){
        prefactor *= static_cast<double>(i);
    }

    prefactor = std::sqrt(static_cast<double>(2 * l + 1) / (MY_4PI * prefactor))
        * associatedLegendre(l, mabs, costheta);

    if((m < 0) && (m % 2)) prefactor = -prefactor;
    return prefactor;
}

double SteinhardtEngine::triangleCoeff(int a, int b, int c){
    return factorial(a + b - c) * factorial(a - b + c) * factorial(-a + b + c)
        / factorial(a + b + c + 1);
}

double SteinhardtEngine::wigner3j(int lmax, int j1, int j2, int j3){
    const int a = lmax, b = lmax, c = lmax;
    const int alpha = j1, beta = j2, gamma = j3;

    struct Denom{
        double operator()(int a, int b, int c, int alpha, int beta, int t) const{
            return factorial(t) * factorial(c - b + t + alpha)
                * factorial(c - a + t - beta) * factorial(a + b - c - t)
                * factorial(a - t - alpha) * factorial(b - t + beta);
        }
    } x;

    const double sgn = 1.0 - 2.0 * static_cast<double>((a - b - gamma) & 1);
    const double g = std::sqrt(triangleCoeff(lmax, lmax, lmax))
        * std::sqrt(
            factorial(a + alpha) * factorial(a - alpha)
            * factorial(b + beta) * factorial(b - beta)
            * factorial(c + gamma) * factorial(c - gamma)
        );

    double s = 0.0;
    int t = 0;
    while(c - b + t + alpha < 0 || c - a + t - beta < 0) ++t;
    while(true){
        if(a + b - c - t < 0) break;
        if(a - t - alpha < 0) break;
        if(b - t + beta < 0) break;
        const double m1t = 1.0 - 2.0 * static_cast<double>(t & 1);
        s += m1t / x(lmax, lmax, lmax, alpha, beta, t);
        ++t;
    }
    return sgn * g * s;
}

// ---------------------------------------------------------------------------
// Engine implementation.
// ---------------------------------------------------------------------------

SteinhardtEngine::SteinhardtEngine(
    ParticleProperty* positions,
    const SimulationCell& simCell,
    const Parameters& parameters,
    bool onlySelected,
    std::shared_ptr<ParticleProperty> selection
)
    : _positions(positions),
      _simCell(simCell),
      _parameters(parameters),
      _onlySelected(onlySelected),
      _selection(std::move(selection))
{
    if(!_positions){
        throw std::invalid_argument("SteinhardtEngine: positions is null");
    }
    if(_parameters.qlist.empty()){
        throw std::invalid_argument("SteinhardtEngine: qlist is empty");
    }

    // Normalization factors for q_l and w_l_hat (matches LAMMPS line 160-162).
    _qnormfac.resize(_parameters.qlist.size());
    _qnormfac2.resize(_parameters.qlist.size());
    _qmax = 0;
    for(std::size_t i = 0; i < _parameters.qlist.size(); ++i){
        const int l = _parameters.qlist[i];
        _qnormfac[i]  = std::sqrt(MY_4PI / (2.0 * l + 1.0));
        _qnormfac2[i] = std::sqrt(2.0 * l + 1.0);
        if(l > _qmax) _qmax = l;
    }

    // Build column layout: [q | w? | w_hat? | components?].
    const int nq = static_cast<int>(_parameters.qlist.size());
    int offset = nq;
    if(_parameters.wlFlag){
        _wColumnOffset = offset;
        offset += nq;
    }
    if(_parameters.wlHatFlag){
        _wHatColumnOffset = offset;
        offset += nq;
    }
    if(_parameters.components >= 0){
        _componentsL = _parameters.components;
        _componentsColumnCount = 2 * (2 * _componentsL + 1);
        _componentsColumnOffset = offset;
        offset += _componentsColumnCount;
    }
    _columnsPerAtom = offset;

    // Precompute the Wigner-3j table if w_l / w_l_hat are requested.
    if(_parameters.wlFlag || _parameters.wlHatFlag){
        initWigner3jTable();
    }

    _q6Index = pickReferenceIndex(_parameters.qlist);

    // Allocate packed output storage.
    _orderParameters = std::make_shared<ParticleProperty>(
        _positions->size(),
        DataType::Double,
        static_cast<std::size_t>(_columnsPerAtom),
        0,
        true
    );

    _structureLabels.assign(_positions->size(), STRUCTURE_LIQUID);
    _structureCounts.fill(0);
}

void SteinhardtEngine::initWigner3jTable(){
    // First pass: count entries per l to build offsets.
    _w3jOffsets.assign(_parameters.qlist.size(), 0);
    std::size_t total = 0;
    for(std::size_t il = 0; il < _parameters.qlist.size(); ++il){
        _w3jOffsets[il] = total;
        const int l = _parameters.qlist[il];
        for(int m1 = -l; m1 <= 0; ++m1){
            for(int m2 = 0; m2 <= ((-m1) >> 1); ++m2){
                ++total;
            }
        }
    }
    _w3jlist.assign(total, 0.0);

    // Second pass: evaluate Wigner-3j and scale by permutation multiplicity
    // (matches the `pfac` logic in LAMMPS init_wigner3j).
    for(std::size_t il = 0; il < _parameters.qlist.size(); ++il){
        const int l = _parameters.qlist[il];
        std::size_t widx = _w3jOffsets[il];
        for(int m1 = -l; m1 <= 0; ++m1){
            for(int m2 = 0; m2 <= ((-m1) >> 1); ++m2){
                const int m3 = -(m1 + m2);
                int pfac;
                if(m1 == 0){
                    pfac = 1;
                }else if(m2 == 0 || m2 == m3){
                    pfac = 6;
                }else{
                    pfac = 12;
                }
                _w3jlist[widx++] = wigner3j(l, m1, m2, m3) * static_cast<double>(pfac);
            }
        }
    }
}

// Core per-atom kernel. `rlist` holds the (x,y,z) bond vectors for the
// selected neighbors, `qn` is the row slice in the packed output where the
// results must land, and `qnm` is a scratch buffer of shape
// (nqlist, qmax+1) storing the per-l, per-m Q_l^m accumulator (real, imag).
void SteinhardtEngine::computeBondOrderParameters(
    const std::vector<std::array<double, 3>>& rlist,
    int ncount,
    double* qn,
    std::vector<std::array<double, 2>>& qnm
) const {
    const int nqlist = static_cast<int>(_parameters.qlist.size());
    const int stride = _qmax + 1;

    // Reset Q_l^m accumulator for this atom.
    for(int il = 0; il < nqlist; ++il){
        const int l = _parameters.qlist[il];
        for(int m = 0; m <= l; ++m){
            qnm[il * stride + m][0] = 0.0;
            qnm[il * stride + m][1] = 0.0;
        }
    }

    // Early-exit: LAMMPS zeroes qn[] when any bond is degenerate.
    for(int ineigh = 0; ineigh < ncount; ++ineigh){
        const auto& r = rlist[ineigh];
        const double rmag = std::sqrt(r[0] * r[0] + r[1] * r[1] + r[2] * r[2]);
        if(rmag <= MY_EPSILON){
            for(int jj = 0; jj < _columnsPerAtom; ++jj) qn[jj] = 0.0;
            return;
        }

        const double costheta = r[2] / rmag;
        double expphi_r = r[0];
        double expphi_i = r[1];
        const double rxymag = std::sqrt(expphi_r * expphi_r + expphi_i * expphi_i);
        if(rxymag <= MY_EPSILON){
            expphi_r = 1.0;
            expphi_i = 0.0;
        }else{
            const double inv = 1.0 / rxymag;
            expphi_r *= inv;
            expphi_i *= inv;
        }

        for(int il = 0; il < nqlist; ++il){
            const int l = _parameters.qlist[il];
            qnm[il * stride + 0][0] += polarPrefactor(l, 0, costheta);
            double expphim_r = expphi_r;
            double expphim_i = expphi_i;
            for(int m = 1; m <= l; ++m){
                const double prefactor = polarPrefactor(l, m, costheta);
                qnm[il * stride + m][0] += prefactor * expphim_r;
                qnm[il * stride + m][1] += prefactor * expphim_i;
                // Advance e^{i(m+1)phi} via one complex multiply.
                const double tmp_r = expphim_r * expphi_r - expphim_i * expphi_i;
                const double tmp_i = expphim_r * expphi_i + expphim_i * expphi_r;
                expphim_r = tmp_r;
                expphim_i = tmp_i;
            }
        }
    }

    // Convert sums to averages.
    const double facn = 1.0 / static_cast<double>(ncount);
    for(int il = 0; il < nqlist; ++il){
        const int l = _parameters.qlist[il];
        for(int m = 0; m <= l; ++m){
            qnm[il * stride + m][0] *= facn;
            qnm[il * stride + m][1] *= facn;
        }
    }

    // Compute q_l = qnormfac[il] * sqrt(q_l^0^2 + 2 * sum_{m>0} |q_l^m|^2).
    int jj = 0;
    for(int il = 0; il < nqlist; ++il){
        const int l = _parameters.qlist[il];
        const double q0 = qnm[il * stride + 0][0];
        double qm_sum = q0 * q0;
        for(int m = 1; m <= l; ++m){
            const double r0 = qnm[il * stride + m][0];
            const double r1 = qnm[il * stride + m][1];
            qm_sum += 2.0 * (r0 * r0 + r1 * r1);
        }
        qn[jj++] = _qnormfac[il] * std::sqrt(qm_sum);
    }

    // Compute w_l using the symmetry-reduced (-l <= m1 <= 0 <= m2 <= m3) loop
    // and the precomputed Wigner-3j table.
    int nterms = 0;
    if(_parameters.wlFlag || _parameters.wlHatFlag){
        for(int il = 0; il < nqlist; ++il){
            const int l = _parameters.qlist[il];
            std::size_t widx = _w3jOffsets[il];
            double wlsum = 0.0;
            for(int m1 = -l; m1 <= 0; ++m1){
                const int sgn = 1 - 2 * (m1 & 1); // (-1)^m1
                for(int m2 = 0; m2 <= ((-m1) >> 1); ++m2){
                    const int m3 = -(m1 + m2);
                    // Q_l^{-m1} = (-1)^{m1} * conj(Q_l^{m1})
                    const double q_m1_r = qnm[il * stride + (-m1)][0];
                    const double q_m1_i = qnm[il * stride + (-m1)][1];
                    const double q_m2_r = qnm[il * stride + m2][0];
                    const double q_m2_i = qnm[il * stride + m2][1];
                    const double q_m3_r = qnm[il * stride + m3][0];
                    const double q_m3_i = qnm[il * stride + m3][1];

                    const double Q1Q2_r = (q_m1_r * q_m2_r + q_m1_i * q_m2_i) * sgn;
                    const double Q1Q2_i = (q_m1_r * q_m2_i - q_m1_i * q_m2_r) * sgn;
                    const double Q1Q2Q3 = Q1Q2_r * q_m3_r - Q1Q2_i * q_m3_i;
                    const double c = _w3jlist[widx++];
                    wlsum += Q1Q2Q3 * c;
                }
            }
            qn[jj++] = wlsum / _qnormfac2[il];
            ++nterms;
        }
    }

    // w_l_hat = w_l * (qnormfac[il] / q_l)^3 * qnormfac2[il] — matches
    // LAMMPS line 540-543 so downstream comparisons agree.
    if(_parameters.wlHatFlag){
        const int jptr = jj - nterms; // start of w_l block
        if(!_parameters.wlFlag) jj = jptr;
        for(int il = 0; il < nqlist; ++il){
            if(qn[il] < QEPSILON){
                qn[jj++] = 0.0;
            }else{
                const double qnfac = _qnormfac[il] / qn[il];
                qn[jj++] = qn[jptr + il] * (qnfac * qnfac * qnfac) * _qnormfac2[il];
            }
        }
    }

    // Optionally emit 2*(2l+1) real components of Q_l^m / |Q_l|.
    if(_parameters.components >= 0 && _componentsColumnOffset >= 0){
        const int target_l = _componentsL;
        // Find il matching target_l; if not present, zero-fill.
        int ilMatch = -1;
        for(int il = 0; il < nqlist; ++il){
            if(_parameters.qlist[il] == target_l){ ilMatch = il; break; }
        }
        if(ilMatch < 0 || qn[ilMatch] < QEPSILON){
            for(int k = 0; k < _componentsColumnCount; ++k){
                qn[_componentsColumnOffset + k] = 0.0;
            }
        }else{
            const double qnfac = _qnormfac[ilMatch] / qn[ilMatch];
            int out = _componentsColumnOffset;
            for(int m = -target_l; m < 0; ++m){
                const int sgn = 1 - 2 * (m & 1); // (-1)^m
                qn[out++] =  qnm[ilMatch * stride + (-m)][0] * qnfac * sgn;
                qn[out++] = -qnm[ilMatch * stride + (-m)][1] * qnfac * sgn;
            }
            for(int m = 0; m <= target_l; ++m){
                qn[out++] = qnm[ilMatch * stride + m][0] * qnfac;
                qn[out++] = qnm[ilMatch * stride + m][1] * qnfac;
            }
        }
    }
}

// ---------------------------------------------------------------------------
// Top-level per-atom driver. Picks nearest-neighbor or cutoff mode and
// launches a TBB parallel_for. Per-thread scratch buffers avoid allocation
// in the hot path.
// ---------------------------------------------------------------------------
void SteinhardtEngine::perform(){
    const std::size_t particleCount = _positions->size();
    if(particleCount == 0){
        spdlog::warn("SteinhardtEngine: no atoms in frame.");
        return;
    }

    const int nqlist = static_cast<int>(_parameters.qlist.size());
    const int stride = _qmax + 1;
    const int* selectionData = _selection ? _selection->constDataInt() : nullptr;

    double* outputBase = _orderParameters->dataDouble();

    const bool useNearest = (_parameters.nnn > 0);
    const int nnn = _parameters.nnn;
    const double cutoff = _parameters.cutoff;
    const double cutsq = cutoff * cutoff;

    if(useNearest && nnn > MAX_NEAREST_NEIGHBORS){
        throw std::invalid_argument(
            "SteinhardtEngine: nnn exceeds MAX_NEAREST_NEIGHBORS (64)"
        );
    }

    if(!useNearest && cutoff <= 0.0){
        throw std::invalid_argument(
            "SteinhardtEngine: cutoff must be > 0 when nnn == 0"
        );
    }

    // Prepare whichever neighbor finder we need. The two code paths share
    // the same per-atom kernel.
    std::unique_ptr<NearestNeighborFinder> nearestFinder;
    std::unique_ptr<CutoffNeighborFinder>  cutoffFinder;

    if(useNearest){
        nearestFinder = std::make_unique<NearestNeighborFinder>(nnn);
        if(!nearestFinder->prepare(_positions, _simCell)){
            throw std::runtime_error("SteinhardtEngine: nearest neighbor finder prepare() failed");
        }
        spdlog::info(
            "Steinhardt: using {} nearest neighbors per atom (nqlist={})",
            nnn, nqlist
        );
    }else{
        cutoffFinder = std::make_unique<CutoffNeighborFinder>();
        if(!cutoffFinder->prepare(cutoff, _positions, _simCell)){
            throw std::runtime_error("SteinhardtEngine: cutoff neighbor finder prepare() failed");
        }
        spdlog::info(
            "Steinhardt: using cutoff = {} (nqlist={})",
            cutoff, nqlist
        );
    }

    struct ThreadScratch{
        std::vector<std::array<double, 2>> qnm;      // (nqlist * stride) complex accumulators
        std::vector<std::array<double, 3>> rlist;    // neighbor bond vectors
        std::array<std::int64_t, 3> counts;
    };

    oneapi::tbb::enumerable_thread_specific<ThreadScratch> scratch(
        [&]{
            ThreadScratch s;
            s.qnm.resize(static_cast<std::size_t>(nqlist) * static_cast<std::size_t>(stride));
            s.rlist.reserve(64);
            s.counts.fill(0);
            return s;
        }
    );

    std::atomic<std::size_t> processedCount{0};
    const std::size_t reportInterval = std::max<std::size_t>(1, particleCount / 20);

    oneapi::tbb::parallel_for(
        oneapi::tbb::blocked_range<std::size_t>(0, particleCount, 1024),
        [&](const oneapi::tbb::blocked_range<std::size_t>& range){
            auto& s = scratch.local();

            for(std::size_t i = range.begin(); i != range.end(); ++i){
                double* row = outputBase + i * static_cast<std::size_t>(_columnsPerAtom);
                for(int k = 0; k < _columnsPerAtom; ++k) row[k] = 0.0;
                _structureLabels[i] = STRUCTURE_LIQUID;

                if(_onlySelected && selectionData && selectionData[i] == 0){
                    s.counts[STRUCTURE_LIQUID]++;
                    continue;
                }

                s.rlist.clear();

                if(useNearest){
                    NearestNeighborFinder::Query<MAX_NEAREST_NEIGHBORS> query(*nearestFinder);
                    query.findNeighbors(i);
                    const auto& results = query.results();
                    if(static_cast<int>(results.size()) < nnn){
                        // Under-coordinated atom: leave qn = 0, classify as liquid.
                        s.counts[STRUCTURE_LIQUID]++;
                    }else{
                        // Collect, sort by distance, keep the `nnn` nearest.
                        struct Neigh{ double d; std::array<double, 3> delta; };
                        std::array<Neigh, MAX_NEAREST_NEIGHBORS> tmp{};
                        int nres = 0;
                        for(const auto& r : results){
                            if(nres >= MAX_NEAREST_NEIGHBORS) break;
                            tmp[nres].d = static_cast<double>(r.distanceSq);
                            tmp[nres].delta = {
                                static_cast<double>(r.delta.x()),
                                static_cast<double>(r.delta.y()),
                                static_cast<double>(r.delta.z())
                            };
                            ++nres;
                        }
                        std::sort(tmp.begin(), tmp.begin() + nres,
                            [](const Neigh& a, const Neigh& b){ return a.d < b.d; });
                        const int take = std::min(nnn, nres);
                        for(int k = 0; k < take; ++k){
                            s.rlist.push_back(tmp[k].delta);
                        }
                        computeBondOrderParameters(s.rlist, take, row, s.qnm);

                        const int label = (_q6Index >= 0)
                            ? classifyFromQ6(row[_q6Index])
                            : STRUCTURE_LIQUID;
                        _structureLabels[i] = label;
                        s.counts[label]++;
                    }
                }else{
                    for(CutoffNeighborFinder::Query nq(*cutoffFinder, i); !nq.atEnd(); nq.next()){
                        if(nq.distanceSquared() >= cutsq) continue;
                        s.rlist.push_back({
                            static_cast<double>(nq.delta().x()),
                            static_cast<double>(nq.delta().y()),
                            static_cast<double>(nq.delta().z())
                        });
                    }
                    const int ncount = static_cast<int>(s.rlist.size());
                    if(ncount == 0){
                        s.counts[STRUCTURE_LIQUID]++;
                    }else{
                        computeBondOrderParameters(s.rlist, ncount, row, s.qnm);
                        const int label = (_q6Index >= 0)
                            ? classifyFromQ6(row[_q6Index])
                            : STRUCTURE_LIQUID;
                        _structureLabels[i] = label;
                        s.counts[label]++;
                    }
                }

                const std::size_t p = processedCount.fetch_add(1, std::memory_order_relaxed) + 1;
                if((p % reportInterval) == 0){
                    spdlog::info("Steinhardt: processed {}/{} atoms", p, particleCount);
                }
            }
        }
    );

    _structureCounts.fill(0);
    for(const auto& s : scratch){
        for(int k = 0; k < 3; ++k) _structureCounts[k] += s.counts[k];
    }

    spdlog::info(
        "Steinhardt: liquid={}, interface={}, crystal={}",
        _structureCounts[STRUCTURE_LIQUID],
        _structureCounts[STRUCTURE_INTERFACE],
        _structureCounts[STRUCTURE_CRYSTAL]
    );
}

const char* SteinhardtEngine::structureLabelName(int label){
    switch(label){
        case STRUCTURE_LIQUID:    return "Liquid-like";
        case STRUCTURE_INTERFACE: return "Interface";
        case STRUCTURE_CRYSTAL:   return "Crystal-like";
        default:                  return "Unknown";
    }
}

}
