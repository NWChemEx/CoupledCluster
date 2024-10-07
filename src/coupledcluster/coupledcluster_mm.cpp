#include "coupledcluster/coupledcluster_mm.hpp"
#include "coupledcluster_modules.hpp"

namespace coupledcluster {

void load_modules(pluginplay::ModuleManager& mm) { 
    // coupledcluster::load_modules<double>(mm); 
    mm.add_module<CCSDEnergy>("CCSD Energy");
}

} // namespace coupledcluster