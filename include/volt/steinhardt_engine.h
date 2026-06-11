#pragma once

#include <volt/core/volt.h>
#include <volt/core/simulation_cell.h>
#include <volt/core/particle_property.h>

#include <array>
#include <cstddef>
#include <memory>
#include <vector>

namespace Volt{

using namespace Volt::Particles;

/**
 * Per-atom Steinhardt bond-orientational order parameters.
 *
 * Ports LAMMPS `compute orientorder/atom` (Thompson & Kohlmeyer) which
 * implements the q_l, w_l and w_l_hat family of order parameters
 * (Steinhardt-Nelson-Ronchetti 1983; Lechner-Dellago 2008).
 *
 * For each atom the engine selects either `nnn` nearest neighbors or all
 * neighbors inside `cutoff`, then averages spherical harmonics Y_lm over
 * those bond directions. The result is one row per atom packed into a
 * single `ParticleProperty` of doubles whose column layout is:
 *   [ q(l_0..l_k)  w(l_0..l_k)?  w_hat(l_0..l_k)?  re/im Q_l^m components? ]
 */
class SteinhardtEngine{
public:
    /** Maximum supported nearest-neighbor count (template bound on Query). */
    static constexpr int MAX_NEAREST_NEIGHBORS = 64;

    /** Heuristic structure labels derived from q_6 (used by the exporter). */
    enum StructureLabel{
        STRUCTURE_LIQUID    = 0,
        STRUCTURE_INTERFACE = 1,
        STRUCTURE_CRYSTAL   = 2
    };

    struct Parameters{
        std::vector<int> qlist = {4, 6, 8, 10, 12};
        int nnn = 12;            ///< 0 = use cutoff mode
        double cutoff = 0.0;     ///< only used when nnn == 0
        bool wlFlag = false;
        bool wlHatFlag = false;
        int components = -1;     ///< emit Q_l^m / |Q_l| for this l (-1 disables)
    };

    SteinhardtEngine(
        ParticleProperty* positions,
        const SimulationCell& simCell,
        const Parameters& parameters,
        bool onlySelected = false,
        std::shared_ptr<ParticleProperty> selection = nullptr
    );

    /// Runs the computation on all atoms in parallel.
    void perform();

    ParticleProperty* positions() const { return _positions; }
    const SimulationCell& cell() const { return _simCell; }
    const Parameters& parameters() const { return _parameters; }

    /// Per-atom packed result; see class docstring for column layout.
    std::shared_ptr<ParticleProperty> orderParameters() const { return _orderParameters; }

    /// Number of columns packed per atom in orderParameters().
    int columnsPerAtom() const { return _columnsPerAtom; }

    /// Column offset at which the q_l block starts (always 0).
    int qColumnOffset() const { return 0; }

    /// Column offset at which the w_l block starts (-1 if disabled).
    int wColumnOffset() const { return _wColumnOffset; }

    /// Column offset at which the w_l_hat block starts (-1 if disabled).
    int wHatColumnOffset() const { return _wHatColumnOffset; }

    /// Column offset at which the Q_l^m / |Q_l| block starts (-1 if disabled).
    int componentsColumnOffset() const { return _componentsColumnOffset; }

    /// Pre-classified structure label per atom (derived from q_6).
    const std::vector<int>& structureLabels() const { return _structureLabels; }

    /// Histogram of structure labels (size 3).
    const std::array<std::int64_t, 3>& structureCounts() const { return _structureCounts; }

    /// Index into qlist that produced q_6 (or closest to l=6) or -1.
    int q6ColumnIndex() const { return _q6Index; }

    static const char* structureLabelName(int label);

private:
    // Precomputed table of Wigner-3j coefficients (scaled by permutation
    // multiplicity), flat-packed across all (l, m1, m2) triples.
    void initWigner3jTable();

    // Core per-atom kernel: computes `ncount` bond directions' spherical
    // harmonic averages and packs the result into the `qn` row.
    void computeBondOrderParameters(
        const std::vector<std::array<double, 3>>& rlist,
        int ncount,
        double* qn,
        std::vector<std::array<double, 2>>& qnm
    ) const;

private:
    // --- LAMMPS-compatible math helpers (ported verbatim) ---
    static double polarPrefactor(int l, int m, double costheta);
    static double associatedLegendre(int l, int m, double x);
    static double triangleCoeff(int a, int b, int c);
    static double wigner3j(int lmax, int j1, int j2, int j3);
    static double factorial(int n);

    ParticleProperty* _positions;
    SimulationCell _simCell;
    Parameters _parameters;
    bool _onlySelected;
    std::shared_ptr<ParticleProperty> _selection;

    // Per-l normalization constants (size = qlist.size()).
    std::vector<double> _qnormfac;
    std::vector<double> _qnormfac2;

    // Flat Wigner-3j table, scaled by permutation multiplicity.
    std::vector<double> _w3jlist;

    // Per-l offset into _w3jlist where that l's block begins.
    std::vector<std::size_t> _w3jOffsets;

    // Layout bookkeeping for the output ParticleProperty.
    int _columnsPerAtom = 0;
    int _wColumnOffset = -1;
    int _wHatColumnOffset = -1;
    int _componentsColumnOffset = -1;
    int _componentsL = -1;
    int _componentsColumnCount = 0;
    int _q6Index = -1;

    // Max l across qlist (used to size the per-thread Q_lm scratch buffer).
    int _qmax = 0;

    // Output storage and structure classification side-product.
    std::shared_ptr<ParticleProperty> _orderParameters;
    std::vector<int> _structureLabels;
    std::array<std::int64_t, 3> _structureCounts{};
};

}
