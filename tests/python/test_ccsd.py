#!/Users/jwaldrop/venvs/nwx/bin/python
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

from pluginplay import ModuleManager
from cc import load_modules, tamm_finalize, tamm_initialize
from simde import AOEnergy, MoleculeFromString, MolecularBasisSet
import chemcache as ccache
from molecules import make_h2
import parallelzone as pz
from chemist import ChemicalSystem
import unittest
import os
import sys


class TestSCF(unittest.TestCase):

    def test_ccsd(self):
        mol_name = "water"
        mol = self.mm.run_as(MoleculeFromString(), "NWX Molecules", mol_name)
        cs = ChemicalSystem(mol)

        basis_name = "sto-3g"
        aos = self.mm.run_as(MolecularBasisSet(), basis_name, mol)

        # key = 'SCF Energy'
        # self.mm.change_input(key, 'molecule_name', mol_name)
        # egy = self.mm.run_as(AOEnergy(), key, aos, cs)
        # self.assertAlmostEqual(egy, -74.3670617803483, places=6)

    def setUp(self):
        self.mm = ModuleManager()
        ccache.load_modules(self.mm)
        load_modules(self.mm)


if __name__ == '__main__':
    tamm_initialize(sys.argv)
    rv = pz.runtime.RuntimeView()

    my_dir = os.path.dirname(os.path.realpath(__file__))

    loader = unittest.TestLoader()
    tests = loader.discover(my_dir)
    testrunner = unittest.runner.TextTestRunner()
    ret = not testrunner.run(tests).wasSuccessful()

    tamm_finalize()

    sys.exit(ret)
