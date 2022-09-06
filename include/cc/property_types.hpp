#pragma once
// #include <simde/correlation_energy.hpp>
#include <simde/simde.hpp>
#include "tamm/tamm.hpp"

namespace coupledcluster {

using BaseType = simde::CorrelationEnergy<simde::type::canonical_reference,simde::type::canonical_reference>;

template<typename ElementType>
DECLARE_DERIVED_TEMPLATED_PROPERTY_TYPE(CCSDResults, BaseType, ElementType);

template<typename ElementType>
TEMPLATED_PROPERTY_TYPE_INPUTS(CCSDResults,ElementType) {
    auto rv = pluginplay::declare_input();
    return rv;
}

template<typename ElementType>
TEMPLATED_PROPERTY_TYPE_RESULTS(CCSDResults,ElementType) {
    using tensor_type = tamm::Tensor<ElementType>;

    auto rv = pluginplay::declare_result()
                .add_field<tensor_type>("Fock MO")
                .template add_field<tensor_type>("t1 amplitude")
                .template add_field<tensor_type>("t2 amplitude");
    return rv;
}


}