#include <volt/steinhardt_service.h>

#include <volt/cli/common.h>
#include <volt/core/analysis_result.h>
#include <volt/core/frame_adapter.h>
#include <volt/core/progress.h>
#include <volt/plugin/output_files.h>
#include <volt/steinhardt_log.h>
#include <volt/utilities/json_utils.h>
#include <volt/utilities/parquet_atom_writer.h>

#include <spdlog/spdlog.h>

#include <cstddef>
#include <exception>
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

std::vector<int> SteinhardtService::parseQList(const std::string& csv){
    std::vector<int> qlist;
    for(const std::string& token : CLI::splitList(csv, ',')){
        try{
            const int l = std::stoi(token);
            if(l >= 0) qlist.push_back(l);
        }catch(...){
        }
    }
    if(qlist.empty()) qlist = {4, 6, 8, 10, 12};
    return qlist;
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
    SteinhardtLog log;

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
        selection = FrameAdapter::findColumnShared(frame, {"Selection", "selection", "v_selection"}, DataType::Int);
        if(!selection){
            spdlog::warn("--onlySelected requested but no Selection column found; classifying all atoms.");
        }
    }

    SteinhardtEngine engine(
        positions.get(),
        frame.simulationCell,
        _parameters,
        _onlySelected,
        selection
    );

    // Weighted by how long each step took on a 100k-atom BCC frame.
    ProgressSteps steps({
        {"Computing Steinhardt order parameters", 1.0},
        {"Writing outputs", 2.0},
    });
    try{
        engine.perform();
    }catch(const std::exception& ex){
        return AnalysisResult::failure(std::string("Engine failed: ") + ex.what());
    }
    log.dump(engine);
    steps.next();

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
        WrittenTable summary{"Steinhardt summary", outputPath(outputBase, "steinhardt")};
        try{
            JsonUtils::writeJsonToParquet(result, summary.path);
        }catch(const std::exception& error){
            summary.error = error.what();
        }
        log.dump(summary);

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
        Volt::streamAtomsToParquet(outputPath(outputBase, "atoms"), frame, bucketResolver, perAtom);
    }

    log.dump(*this);
    return result;
}

}
