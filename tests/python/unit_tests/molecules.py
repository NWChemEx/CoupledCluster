# Copyright 2023 NWChemEx-Project
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
# http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

from chemist import Atom, Molecule, ChemicalSystem


def make_h2():
    mol = Molecule()
    mol.push_back(Atom('H', 1, 1837.15264648179, 0.0, 0.0, 0.0))
    mol.push_back(Atom('H', 1, 1837.15264648179, 0.0, 0.0, 1.68185))

    return ChemicalSystem(mol)
