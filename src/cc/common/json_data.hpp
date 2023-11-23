#pragma once

#include <cctype>
#include <iostream>
#include <regex>
#include <string>
#include <vector>

#include "tamm/eigen_utils.hpp"
#include "tamm/tamm.hpp"

using namespace tamm;

// Libint Gaussian integrals library
#include <libint2.hpp>
#include <libint2/basis.h>
#include <libint2/chemistry/sto3g_atomic_density.h>

#include "ga/ga-mpi.h"
#include "ga/ga.h"

#include <nlohmann/json.hpp>
using json = nlohmann::ordered_json;

using libint2::Atom;
using std::cerr;
using std::cout;
using std::endl;
using std::string;

#include "utils.hpp"

#include <filesystem>
namespace fs = std::filesystem;

struct SystemData {
  OptionsMap options_map;
  int        n_occ_alpha{};
  int        n_vir_alpha{};
  int        n_occ_beta{};
  int        n_vir_beta{};
  int        n_lindep;
  int        ndf{};
  int        nbf{};
  int        nbf_orig{};
  int        nelectrons{};
  int        nelectrons_alpha{};
  int        nelectrons_beta{};
  int        n_frozen_core{};
  int        n_frozen_virtual{};
  int        nmo{};
  int        nocc{};
  int        nvir{};
  int        focc{};
  bool       ediis{};
  bool       is_restricted{};
  bool       is_unrestricted{};
  bool       is_restricted_os{};
  bool       is_ks{};

  std::string basis;
  std::string scf_type_string;
  std::string input_molecule;
  std::string output_file_prefix;

  // output data
  double scf_energy{};
  int    num_chol_vectors{};
  double ccsd_corr_energy{};

  // SCF
  int AO_tilesize{30};
  int dfAO_tilesize{30};

  // json data
  json results;

  void print() {
    std::cout << std::endl << "----------------------------" << std::endl;
    std::cout << "scf_type = " << scf_type_string << std::endl;
    if(is_restricted) std::cout << "Closed-Shell SCF" << std::endl;
    if(is_unrestricted) std::cout << "Open-Shell SCF" << std::endl;
    if(is_restricted_os) std::cout << "Restricted Open-Shell SCF" << std::endl;
    if(is_ks) std::cout << "KS-DFT Enabled" << std::endl;

    std::cout << "nbf = " << nbf << std::endl;
    std::cout << "nbf_orig = " << nbf_orig << std::endl;
    std::cout << "n_lindep = " << n_lindep << std::endl;

    std::cout << "focc = " << focc << std::endl;
    std::cout << "nmo = " << nmo << std::endl;
    std::cout << "nocc = " << nocc << std::endl;
    std::cout << "nvir = " << nvir << std::endl;

    std::cout << "n_occ_alpha = " << n_occ_alpha << std::endl;
    std::cout << "n_vir_alpha = " << n_vir_alpha << std::endl;
    std::cout << "n_occ_beta = " << n_occ_beta << std::endl;
    std::cout << "n_vir_beta = " << n_vir_beta << std::endl;

    std::cout << "nelectrons = " << nelectrons << std::endl;
    std::cout << "nelectrons_alpha = " << nelectrons_alpha << std::endl;
    std::cout << "nelectrons_beta = " << nelectrons_beta << std::endl;
    std::cout << "n_frozen_core = " << n_frozen_core << std::endl;
    std::cout << "n_frozen_virtual = " << n_frozen_virtual << std::endl;
    std::cout << "----------------------------" << std::endl;
  }

  void update() {
    EXPECTS(nbf == n_occ_alpha + n_vir_alpha); // lin-deps
    // EXPECTS(nbf_orig == n_occ_alpha + n_vir_alpha + n_lindep + n_frozen_core + n_frozen_virtual);
    nocc = n_occ_alpha + n_occ_beta;
    nvir = n_vir_alpha + n_vir_beta;
    // EXPECTS(nelectrons == n_occ_alpha + n_occ_beta);
    // EXPECTS(nelectrons == nelectrons_alpha+nelectrons_beta);
    nmo = n_occ_alpha + n_vir_alpha + n_occ_beta + n_vir_beta; // lin-deps
  }

  SystemData(OptionsMap options_map_, const std::string scf_type_string):
    options_map(options_map_), scf_type_string(scf_type_string) {
    results          = json::object();
    is_restricted    = false;
    is_unrestricted  = false;
    is_restricted_os = false;
    is_ks            = false;
    if(scf_type_string == "restricted") {
      focc          = 1;
      is_restricted = true;
    }
    else if(scf_type_string == "unrestricted") {
      focc            = 2;
      is_unrestricted = true;
    }
    else if(scf_type_string == "restricted_os") {
      focc             = -1;
      is_restricted_os = true;
    }
    else
      tamm_terminate("ERROR: unrecognized scf_type [" + scf_type_string + "] provided");
    // if(!options_map_.scf_options.xc_type.empty()) { is_ks = true; }
  }
};

void check_json(std::string filename) {
  std::string get_ext = fs::path(filename).extension();
  const bool  is_json = (get_ext == ".json");
  if(!is_json) tamm_terminate("ERROR: Input file provided [" + filename + "] must be a json file");
}

std::string getfilename(std::string filename) {
  size_t lastindex = filename.find_last_of(".");
  auto   fname     = filename.substr(0, lastindex);
  return fname.substr(fname.find_last_of("/") + 1, fname.length());
}

void write_json_data(SystemData& sys_data, const std::string module) {
  auto options = sys_data.options_map;
  auto cd      = options.cd_options;
  auto ccsd    = options.ccsd_options;

  json& results = sys_data.results;

  auto str_bool = [=](const bool val) {
    if(val) return "true";
    return "false";
  };

  results["input"]["molecule"]["name"]  = sys_data.input_molecule;
  results["input"]["molecule"]["basis"] = sys_data.basis;
  // results["input"]["molecule"]["basis_sphcart"] = scf.sphcart;
  // results["input"]["molecule"]["geometry_units"] = scf.geom_units;
  // SCF options
  // results["input"]["SCF"]["tol_int"] = scf.tol_int;
  // results["input"]["SCF"]["tol_lindep"] = scf.tol_lindep;
  // results["input"]["SCF"]["conve"] = scf.conve;
  // results["input"]["SCF"]["convd"] = scf.convd;
  // results["input"]["SCF"]["diis_hist"] = scf.diis_hist;
  results["input"]["SCF"]["AO_tilesize"] = sys_data.AO_tilesize;
  // results["input"]["SCF"]["force_tilesize"] = str_bool(scf.force_tilesize);
  results["input"]["SCF"]["scf_type"] = sys_data.scf_type_string;
  // results["input"]["SCF"]["multiplicity"] = scf.multiplicity;

  if(module == "CD" || module == "CCSD") {
    // CD options
    results["input"]["CD"]["diagtol"]          = cd.diagtol;
    results["input"]["CD"]["itilesize"]        = cd.itilesize;
    results["input"]["CD"]["write_cv"]         = cd.write_cv;
    results["input"]["CD"]["write_vcount"]     = cd.write_vcount;
    results["input"]["CD"]["max_cvecs_factor"] = cd.max_cvecs_factor;
  }

  if(module == "CCSD") {
    // CCSD options
    results["input"]["CCSD"]["threshold"]     = ccsd.threshold;
    results["input"]["CCSD"]["tilesize"]      = ccsd.tilesize;
    results["input"]["CCSD"]["ndiis"]         = ccsd.ndiis;
    results["input"]["CCSD"]["readt"]         = str_bool(ccsd.readt);
    results["input"]["CCSD"]["writet"]        = str_bool(ccsd.writet);
    results["input"]["CCSD"]["ccsd_maxiter"]  = ccsd.ccsd_maxiter;
    results["input"]["CCSD"]["balance_tiles"] = str_bool(ccsd.balance_tiles);
    // results["input"]["CCSD"]["cache_size"]         = ccsd.cache_size;
  }

  std::string l_module = module;
  std::transform(l_module.begin(), l_module.end(), l_module.begin(), ::tolower);

  std::string out_fp = sys_data.output_file_prefix + "." + sys_data.options_map.ccsd_options.basis;
  std::string files_dir    = out_fp + "_files/" + sys_data.scf_type_string;
  std::string files_prefix = files_dir + "/" + out_fp;
  std::string json_file    = files_prefix + "." + l_module + ".json";
  bool        json_exists  = std::filesystem::exists(json_file);
  if(json_exists) {
    // std::ifstream jread(json_file);
    // jread >> results;
    std::filesystem::remove(json_file);
  }

  // std::cout << std::endl << std::endl << results.dump() << std::endl;
  std::ofstream res_file(json_file);
  res_file << std::setw(2) << results << std::endl;
}

std::tuple<std::vector<size_t>, std::vector<Tile>, std::vector<Tile>>
compute_AO_tiles(const ExecutionContext& ec, const SystemData& sys_data, libint2::BasisSet& shells,
                 const bool is_df = false) {
  auto rank      = ec.pg().rank();
  int  tile_size = sys_data.AO_tilesize;
  if(is_df) tile_size = sys_data.dfAO_tilesize;

  std::vector<Tile> AO_tiles;
  for(auto s: shells) AO_tiles.push_back(s.size());
  // if(rank == 0) cout << "Number of AO tiles = " << AO_tiles.size() << endl;

  tamm::Tile          est_ts = 0;
  std::vector<Tile>   AO_opttiles;
  std::vector<size_t> shell_tile_map;
  for(size_t s = 0; s < shells.size(); s++) {
    est_ts += shells[s].size();
    if(est_ts >= (size_t) tile_size) {
      AO_opttiles.push_back(est_ts);
      shell_tile_map.push_back(s); // shell id specifying tile boundary
      est_ts = 0;
    }
  }
  if(est_ts > 0) {
    AO_opttiles.push_back(est_ts);
    shell_tile_map.push_back(shells.size() - 1);
  }

  return std::make_tuple(shell_tile_map, AO_tiles, AO_opttiles);
}
