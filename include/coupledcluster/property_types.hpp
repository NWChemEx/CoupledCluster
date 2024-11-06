/*
 * Copyright 2024 NWChemEx-Project
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#pragma once
#include "simde/simde.hpp"
#include <tamm/tamm.hpp>

namespace coupledcluster {

template<typename BraType, typename OpType, typename KetType>
DECLARE_TEMPLATED_PROPERTY_TYPE(BraKet, BraType, OpType, KetType);

template<typename BraType, typename OpType, typename KetType>
TEMPLATED_PROPERTY_TYPE_INPUTS(BraKet, BraType, OpType, KetType) {
  using op_t = const OpType&;

  auto rv = pluginplay::declare_input()
              .add_field<const BraType&>("Bra")
              .template add_field<op_t>("Operator")
              .template add_field<const KetType&>("Ket");
  return rv;
}

template<typename BraType, typename OpType, typename KetType>
TEMPLATED_PROPERTY_TYPE_RESULTS(BraKet, BraType, OpType, KetType) {
  auto rv = pluginplay::declare_result()
              .add_field<tamm::Tensor<double>>("Fock MO")
              .template add_field<tamm::Tensor<double>>("T1 amplitude")
              .template add_field<tamm::Tensor<double>>("T2 amplitude");
  return rv;
}

template<typename BraType, typename KetType>
using CorrelationEnergy = BraKet<BraType, simde::type::els_hamiltonian, KetType>;

template<typename BraType>
using ElectronicEnergy = BraKet<BraType, simde::type::els_hamiltonian, BraType>;

template<typename BraType>
using TotalEnergy = BraKet<BraType, simde::type::hamiltonian, BraType>;

// CCSD_T

template<typename T>
DECLARE_TEMPLATED_PROPERTY_TYPE(CCSD_T_PT, T);

template<typename T>
TEMPLATED_PROPERTY_TYPE_INPUTS(CCSD_T_PT, T) {
  auto rv = pluginplay::declare_input()
              .add_field<tamm::Tensor<T>>("Fock MO")
              .template add_field<tamm::Tensor<T>>("T1 amplitude")
              .template add_field<tamm::Tensor<T>>("T2 amplitude");
  return rv;
}

template<typename T>
TEMPLATED_PROPERTY_TYPE_RESULTS(CCSD_T_PT, T) {
  auto rv = pluginplay::declare_result().add_field<T>("Energy");
  return rv;
}

} // namespace coupledcluster
