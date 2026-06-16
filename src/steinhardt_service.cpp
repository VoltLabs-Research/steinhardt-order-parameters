#include <volt/steinhardt_service.h>

#include <volt/core/analysis_result.h>
#include <volt/core/frame_adapter.h>
#include <volt/utilities/json_utils.h>
#include <volt/utilities/parquet_atom_writer.h>

#include <spdlog/spdlog.h>

#include <chrono>
#include <cstddef>
#include <memory>
#include <string>
#include <vector>

namespace Volt{

using namespace Volt::Particles;

SteinhardtService::SteinhardtService()
    : _parameters(),
      _onlySelected(false){}

void SteinhardtService::setQList(const std::vector<int>& qlist){
    _parameters.qlist = qlist;
}

void SteinhardtService::setNnn(int nnn){
    _parameters.nnn = nnn;
}

void SteinhardtService::setCutoff(double cutoff){
    _parameters.cutoff = cutoff;
}

void SteinhardtService::setWlFlag(bool enabled){
    _parameters.wlFlag = enabled;
}

void SteinhardtService::setWlHatFlag(bool enabled){
    _parameters.wlHatFlag = enabled;
}

void SteinhardtService::setComponents(int l){
    _parameters.components = l;
}

void SteinhardtService::setOnlySelected(bool enabled){
    _onlySelected = enabled;
}

json SteinhardtService::compute(const LammpsParser::Frame& frame, const std::string& outputBase){
    const auto startTime = std::chrono::high_resolution_clock::now();

    if(frame.natoms <= 0){
        return AnalysisResult::failure("Invalid number of atoms");
    }

    if(!FrameAdapter::validateSimulationCell(frame.simulationCell)){
        return AnalysisResult::failure("Invalid simulation cell");
    }

    auto positions = FrameAdapter::createPositionPropertyShared(frame);
    if(!positions){
        return AnalysisResult::failure("Failed to create position property");
    }

    if(_parameters.qlist.empty()){
        return AnalysisResult::failure("qlist is empty");
    }

    if(_parameters.nnn == 0 && _parameters.cutoff <= 0.0){
        return AnalysisResult::failure("--cutoff must be > 0 when --nnn 0 is used");
    }

    // Mirror the convention used by the rest of VOLT: look for a
    // Selection column only when the user asked for onlySelected.
    std::shared_ptr<ParticleProperty> selection;
    if(_onlySelected){
        for(const char* col : {"Selection", "selection", "v_selection"}){
            selection = FrameAdapter::createIntPropertyShared(frame, col);
            if(selection) break;
        }
        if(!selection){
            spdlog::warn("--onlySelected requested but no Selection column found; classifying all atoms.");
        }
    }

    spdlog::info(
        "Starting Steinhardt analysis (nqlist={}, nnn={}, cutoff={}, wl={}, wlHat={}, components={}, onlySelected={})",
        _parameters.qlist.size(),
        _parameters.nnn,
        _parameters.cutoff,
        _parameters.wlFlag,
        _parameters.wlHatFlag,
        _parameters.components,
        _onlySelected
    );

    SteinhardtEngine engine(
        positions.get(),
        frame.simulationCell,
        _parameters,
        _onlySelected,
        selection
    );

    try{
        engine.perform();
    }catch(const std::exception& ex){
        return AnalysisResult::failure(std::string("Engine failed: ") + ex.what());
    }

    auto orderProp = engine.orderParameters();
    const double* base = orderProp ? orderProp->constDataDouble() : nullptr;
    const int cols = engine.columnsPerAtom();
    const auto& labels = engine.structureLabels();
    const int q6col = engine.q6ColumnIndex();        // index into qlist (qColumnOffset is 0)
    const int qOff = engine.qColumnOffset();
    const int wOff = engine.wColumnOffset();
    const int whOff = engine.wHatColumnOffset();
    const int cOff = engine.componentsColumnOffset();
    const auto& ql = _parameters.qlist;
    const int nq = static_cast<int>(ql.size());

    // --- Build the canonical _steinhardt.parquet summary ------------------
    json result;
    result["main_listing"] = {
        {"total_atoms", frame.natoms},
        {"qlist", ql},
        {"nnn", _parameters.nnn},
        {"cutoff", _parameters.cutoff},
        {"wl", _parameters.wlFlag},
        {"wl_hat", _parameters.wlHatFlag},
        {"components", _parameters.components},
        {"liquid_like", engine.structureCounts()[SteinhardtEngine::STRUCTURE_LIQUID]},
        {"interface", engine.structureCounts()[SteinhardtEngine::STRUCTURE_INTERFACE]},
        {"crystal_like", engine.structureCounts()[SteinhardtEngine::STRUCTURE_CRYSTAL]}
    };

    json structRows = json::array();
    for(int s = 0; s < 3; ++s){
        const std::int64_t count = engine.structureCounts()[s];
        if(count <= 0) continue;
        structRows.push_back({
            {"structure_id", s},
            {"structure_name", SteinhardtEngine::structureLabelName(s)},
            {"atom_count", count}
        });
    }
    result["sub_listings"] = {{"structures", structRows}};

    if(!outputBase.empty()){
        const std::string summaryPath = outputBase + "_steinhardt.parquet";
        if(JsonUtils::writeJsonToParquet(result, summaryPath)){
            spdlog::info("Steinhardt summary parquet written to {}", summaryPath);
        }else{
            spdlog::warn("Could not write Steinhardt summary parquet: {}", summaryPath);
        }

        // Group by the q_6-derived structure label so the AtomisticExporter GLB
        // splits Liquid-like / Interface / Crystal-like.
        auto bucketResolver = [&labels](std::size_t i) -> std::string {
            const int label = i < labels.size() ? labels[i] : SteinhardtEngine::STRUCTURE_LIQUID;
            return std::string(SteinhardtEngine::structureLabelName(label));
        };
        auto perAtom = [&](ColumnarAtomWriter& w, std::size_t i){
            const double* row = base ? base + i * static_cast<std::size_t>(cols) : nullptr;

            // Per-l scalar columns (q_4, q_6, ...) for direct coloring/filtering,
            // plus the whole q_l set as a single list<double> for downstream analysis.
            std::vector<double> qVec(static_cast<std::size_t>(nq), 0.0);
            for(int k = 0; k < nq; ++k){
                const double v = row ? row[qOff + k] : 0.0;
                w.field("q_" + std::to_string(ql[k]), v);
                qVec[static_cast<std::size_t>(k)] = v;
            }
            // The classifying scalar stays a scalar column.
            w.field("q6", (row && q6col >= 0) ? row[qOff + q6col] : 0.0);
            w.field("q", qVec);

            if(_parameters.wlFlag){
                for(int k = 0; k < nq; ++k){
                    w.field("w_" + std::to_string(ql[k]), (row && wOff >= 0) ? row[wOff + k] : 0.0);
                }
            }
            if(_parameters.wlHatFlag){
                for(int k = 0; k < nq; ++k){
                    w.field("w_hat_" + std::to_string(ql[k]), (row && whOff >= 0) ? row[whOff + k] : 0.0);
                }
            }
            if(_parameters.components >= 0){
                const int n = 2 * (2 * _parameters.components + 1);
                std::vector<double> comp(static_cast<std::size_t>(n), 0.0);
                if(row && cOff >= 0){
                    for(int j = 0; j < n; ++j) comp[static_cast<std::size_t>(j)] = row[cOff + j];
                }
                w.field("components", comp);
            }
        };

        // Steinhardt is not structural identification — the q6 label drives only the
        // bucket/q6 property column, not structure_id (which would clobber an upstream
        // classifier in a pipeline). structure_id/structure_name stay off (default).
        Volt::streamAtomsToParquet(outputBase + "_atoms.parquet", frame, bucketResolver, perAtom);
    }

    const auto endTime = std::chrono::high_resolution_clock::now();
    const auto elapsedMs = std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime).count();
    spdlog::info("Steinhardt analysis finished in {} ms", elapsedMs);

    return result;
}

}
