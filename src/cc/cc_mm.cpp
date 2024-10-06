#include "cc/cc_mm.hpp"
#include "cc_modules.hpp"

namespace cc {

void load_modules(pluginplay::ModuleManager& mm) { 
    // cc::load_modules<double>(mm); 
    mm.add_module<CCSDEnergy>("CCSD Energy");
}

} // namespace cc