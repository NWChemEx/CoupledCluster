#pragma once

#include "tamm/tamm.hpp"

using namespace tamm;

inline auto sum_tensor_sizes = [](auto&&... t) {
  return ((compute_tensor_size(t) + ...) * 8) / (1024 * 1024 * 1024.0);
};

void check_memory_requirements(ExecutionContext& ec, double calc_mem) {
  auto minfo = ec.mem_info();
  if(calc_mem > static_cast<double>(minfo.total_cpu_mem)) {
    ec.print_mem_info();
    std::string err_msg = "ERROR: Insufficient CPU memory, required = " + std::to_string(calc_mem) +
                          "GiB, available = " + std::to_string(minfo.total_cpu_mem) + " GiB";
    tamm_terminate(err_msg);
  }
}
