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

#include <catch2/catch_approx.hpp>
#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>
#include <chemcache/chemcache.hpp>
#include <iostream>
#include "cc/cc.hpp"

using namespace mokup;

using pt     = simde::CanonicalCorrelationEnergy;
using eri_pt = simde::TransformedERI4;

TEST_CASE("CCSD") { 
    std::cout << "TBD: Canonical CCSD Test" << std::endl; 
}
