#pragma once

#include <volt/core/volt.h>
#include <volt/core/lammps_parser.h>
#include <volt/core/particle_property.h>
#include <volt/steinhardt_engine.h>

#include <nlohmann/json.hpp>

#include <string>
#include <vector>

namespace Volt{

using json = nlohmann::json;

/**
 * Orchestrates the per-frame Steinhardt order parameter analysis:
 * validates the input frame, drives the engine, writes the canonical
 * `_steinhardt.parquet` summary plus the AtomisticExporter-compatible
 * `_atoms.parquet` with q_6-based structure grouping.
 */
class SteinhardtService{
public:
    SteinhardtService();

    void setQList(const std::vector<int>& qlist);
    void setNnn(int nnn);
    void setCutoff(double cutoff);
    void setWlFlag(bool enabled);
    void setWlHatFlag(bool enabled);
    void setComponents(int l);
    void setOnlySelected(bool enabled);

    const SteinhardtEngine::Parameters& parameters() const { return _parameters; }

    json compute(
        const LammpsParser::Frame& frame,
        const std::string& outputBase = ""
    );

private:
    SteinhardtEngine::Parameters _parameters;
    bool _onlySelected;
};

}
