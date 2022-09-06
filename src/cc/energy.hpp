#pragma once
#include <simde/simde.hpp>

namespace ccsd {

template<typename T>
DECLARE_MODULE(CCSD);

extern template class CCSD<double>;

template<typename T>
inline void load_energy_modules(pluginplay::ModuleManager& mm) {
    mm.add_module<CCSD<T>>("CCSD");
}

} // namespace ccsd
