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

#include <iostream>

#include <coupledcluster/coupledcluster.hpp>
#include <chemcache/chemcache.hpp>
#include <tamm/tamm.hpp>

int main(int argc, char** argv) {

    tamm::initialize(argc, argv);

    // Populate modules
    pluginplay::ModuleManager mm;
    chemcache::load_modules(mm);
    coupledcluster::load_modules(mm);

    // Create ChemicalSystem
    std::string mol_name = "water";
    auto mol = mm.at("NWX Molecules").run_as<simde::MoleculeFromString>(mol_name);
    simde::type::chemical_system cs(mol);

    mm.change_input("CCSD Energy", "molecule_name", mol_name);

    // Create BasisSet
    std::string basis_name = "sto-3g"; // This is the only supported basis in ChemCache
    auto aos = mm.at(basis_name).run_as<simde::MolecularBasisSet>(mol);

    // Run module
    auto E = mm.at("CCSD Energy").run_as<simde::AOEnergy>(aos, cs);
    std::cout << std::endl << "CCSD Energy = " << std::fixed << std::setprecision(12) << E << " Hartree" << std::endl;
    
    tamm::finalize();
    
    return 0;
}
