#include "cc/cc_mm.hpp"
#include "cc_modules.hpp"

namespace ccsd {

void load_modules(pluginplay::ModuleManager& mm) { ccsd::load_modules<double>(mm); }

} // namespace ccsd