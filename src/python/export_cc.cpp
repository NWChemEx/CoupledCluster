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

#include <cc/cc.hpp>
#include <pluginplay/plugin/plugin.hpp>
#include <pybind11/pybind11.h>
#include <tamm/tamm.hpp>

namespace cc {

EXPORT_PLUGIN(scf, m) {
    m.def("tamm_initialize", [](pybind11::list py_args) {
        std::vector<std::string> args;
        for (const auto& arg : py_args) 
            args.push_back(arg.cast<std::string>());

        std::vector<char*> argv;
        for (const auto& arg : args) 
            argv.push_back(const_cast<char*>(arg.c_str()));

        int argc = static_cast<int>(argv.size());        
        tamm::initialize(argc, argv.data());
    });
    m.def("tamm_finalize", []() {
        tamm::finalize();
    });
}

} // namespace cc
