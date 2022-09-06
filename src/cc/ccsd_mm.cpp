#include "energy.hpp"
#include "cc/ccsd_mm.hpp"

template<typename T>
void load_modules_(pluginplay::ModuleManager& mm) {
    ccsd::load_energy_modules<T>(mm);
}

namespace ccsd{
    
void load_modules(pluginplay::ModuleManager& mm) {
    load_modules_<double>(mm);
}

} // namespace ccsd