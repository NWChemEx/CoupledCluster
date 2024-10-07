#pragma once
#include <pluginplay/pluginplay.hpp>

namespace coupledcluster {

// template<typename T>
DECLARE_MODULE(CCSDEnergy);

// inline void set_defaults(pluginplay::ModuleManager& mm) {
//   // mm.change_submod("CCSD", "Energy", "Energy");
// }

// template<typename T>
// inline void load_modules(pluginplay::ModuleManager& mm) {
//   mm.add_module<CCSD<T>>("CCSD");
//   set_defaults(mm);
// }

// extern template class CCSD<double>;

} // namespace coupledcluster
