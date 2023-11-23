#pragma once
#include <simde/simde.hpp>

namespace ccsd {

template<typename T>
DECLARE_MODULE(CCSD);

inline void set_defaults(pluginplay::ModuleManager& mm) {
  // mm.change_submod("CCSD", "Energy", "Energy");
}

template<typename T>
inline void load_modules(pluginplay::ModuleManager& mm) {
  mm.add_module<CCSD<T>>("CCSD");
  set_defaults(mm);
}

extern template class CCSD<double>;

} // namespace ccsd
