#pragma once

#include <volt/core/run_log.h>
#include <volt/steinhardt_engine.h>

namespace Volt{

class SteinhardtService;

// The execution log of one analysis: each phase hands it what it built.
class SteinhardtLog : public RunLog{
public:
    using RunLog::dump;

    void dump(const SteinhardtEngine& engine);
    // The service is done with the run: closes it with its total time.
    void dump(const SteinhardtService& service);
};

}
