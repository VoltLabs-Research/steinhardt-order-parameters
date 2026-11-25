#include <volt/steinhardt_log.h>

#include <spdlog/spdlog.h>

namespace Volt{

void SteinhardtLog::dump(const SteinhardtEngine& engine){
    const auto& counts = engine.structureCounts();
    spdlog::info(
        "Steinhardt: liquid={}, interface={}, crystal={}",
        counts[SteinhardtEngine::STRUCTURE_LIQUID],
        counts[SteinhardtEngine::STRUCTURE_INTERFACE],
        counts[SteinhardtEngine::STRUCTURE_CRYSTAL]
    );
}

void SteinhardtLog::dump(const SteinhardtService&){
    spdlog::info("Steinhardt analysis finished in {} ms", runMilliseconds());
}

}
