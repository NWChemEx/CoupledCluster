// clang-format off
#include "cd_ccsd_cs_ann.hpp"
#include "cc_modules.hpp"
#include <libint2.hpp>
#include <simde/simde.hpp>
#include "ccsd_t/ccsd_t_fused_driver.hpp"
#include "cc/property_types.hpp"
// clang-format on

namespace ccsd {

using f_pt   = simde::Fock;
using eri_pt = simde::TransformedERI4;
using ee_pt  = simde::CanonicalElectronicEnergy;

using ce_pt =
  simde::CorrelationEnergy<simde::type::canonical_reference, simde::type::canonical_reference>;

using ccsd_pt = coupledcluster::CorrelationEnergy<simde::type::canonical_reference,
                                                  simde::type::canonical_reference>;

inline libint2::BasisSet ccsd_make_basis(const simde::type::ao_basis_set& bs) {
  /// Typedefs for everything
  using atom_t          = libint2::Atom;
  using shell_t         = libint2::Shell;
  using basis_t         = libint2::BasisSet;
  using cont_t          = libint2::Shell::Contraction;
  using svec_d_t        = libint2::svector<double>;
  using conts_t         = libint2::svector<cont_t>;
  using centers_t       = std::vector<atom_t>;
  using atom_bases_t    = std::vector<shell_t>;
  using element_bases_t = std::vector<atom_bases_t>;

  /// Inputs for BasisSet constructor
  centers_t       centers{};
  element_bases_t element_bases{};

  /// Atom doesn't have a value ctor, so here's a stand in
  auto atom_ctor = [](int Z, double x, double y, double z) {
    atom_t atom{};
    atom.atomic_number = Z;
    atom.x             = x;
    atom.y             = y;
    atom.z             = z;
    return atom;
  };

  /// Origin for shell construction
  std::array<double, 3> origin = {0.0, 0.0, 0.0};

  /// Convert centers and their shells to libint equivalents.
  for(auto abs_i = 0; abs_i < bs.size(); ++abs_i) {
    /// Add current center to atoms list
    const auto& abs = bs[abs_i];
    centers.push_back(atom_ctor(abs_i, abs.center().x(), abs.center().y(), abs.center().z()));

    /// Gather shells for this center and add them to element_bases
    atom_bases_t atom_bases{};
    for(const auto&& shelli: abs) {
      const auto nprims = shelli.n_primitives();
      const auto prim0  = shelli.primitive(0);
      const auto primN  = shelli.primitive(nprims - 1);
      const bool pure   = shelli.pure() == chemist::ShellType::pure;
      const int  l      = shelli.l();

      svec_d_t alphas(&prim0.exponent(), &primN.exponent() + 1);
      svec_d_t coefs(&prim0.coefficient(), &primN.coefficient() + 1);
      conts_t  conts{cont_t{l, pure, coefs}};
      /// Use origin for position, because BasisSet moves shells to center
      atom_bases.push_back(shell_t(alphas, conts, origin));
    }
    element_bases.push_back(atom_bases);
  }

  /// Return the new basis set
  return basis_t(centers, element_bases);
}

// DECLARE_MODULE(CCSD);

template<typename T>
TEMPLATED_MODULE_CTOR(CCSD, T) {
  description("Computes CCSD correlation energy");
  satisfies_property_type<ce_pt>();
  satisfies_property_type<ccsd_pt>();
  satisfies_property_type<simde::TotalCanonicalEnergy>();

  add_submodule<f_pt>("Fock Builder");
  add_submodule<eri_pt>("Transformed ERIs");
  add_submodule<ee_pt>("Electronic Energy");

  add_input<bool>("debug").set_default(false).set_description("Debugging flag");

  add_input<double>("diagtol").set_default(1.0e-5).set_description(
    "Cholesky Decomposition Threshold");

  // write to disk after every count number of vectors are computed.
  // enabled only if write_cv=true and nbf>1000
  add_input<bool>("write_cv").set_default(false).set_description("write chol vecs to disk");
  add_input<int>("write_vcount")
    .set_default(5000)
    .set_description("write to disk after every count number of vectors are computed");

  add_input<int>("max_cvecs_factor")
    .set_default(12)
    .set_description("Limit Max. number of cholesky vectors to 12*N");

  add_input<double>("printtol")
    .set_default(0.05)
    .set_description("Write T1,T2 amplitudes above a certain threshold to a file");

  add_input<double>("threshold").set_default(1.0e-6).set_description("CCSD Threshold");

  add_input<int>("tilesize")
    .set_default(50)
    .set_description("Tilesize for the MO space. Will be reset automatically "
                     "based on MO size");

  add_input<bool>("force_tilesize").set_default(false).set_description("Force tilesize specified");

  add_input<int>("itilesize")
    .set_default(1000)
    .set_description("Tilesize for the Cholesky Dimension");

  add_input<int>("ndiis").set_default(5).set_description("number of diis entries");

  add_input<double>("lshift").set_default(0.0).set_description("Level Shift");

  add_input<int>("ccsd_maxiter").set_default(50).set_description("Maximum number of iterations");

  add_input<int>("freeze_core")
    .set_default(0)
    .set_description("Specify number of core orbitals to freeze");

  add_input<int>("freeze_virtual")
    .set_default(0)
    .set_description("Specify number of virtuals to freeze");

  add_input<bool>("balance_tiles")
    .set_default(true)
    .set_description("Balanced tiling scheme, determined automatically based "
                     "on CC module being run");

  add_input<bool>("profile_ccsd")
    .set_default(false)
    .set_description("Write profiling information to csv file");

  add_input<bool>("writet").set_default(false).set_description(
    "Write Fock, 2e integral, T1, T2 amplitude tensors to disk");

  add_input<int>("writet_iter")
    .set_default(5)
    .set_description("Write Fock, 2e integral, T1, T2 amplitude tensors to "
                     "disk after every 5 iterations");

  add_input<bool>("readt").set_default(false).set_description(
    "Read Fock, 2e integral, T1, T2 amplitude tensors to disk. Not required "
    "when writet=true");

  add_input<bool>("writev").set_default(false).set_description(
    "Write the 4D 2e integral tensor to disk");

  add_input<bool>("computeTData")
    .set_default(true)
    .set_description("Compute and write data needed for (T) calculation to disk");

  add_input<int>("cache_size").set_default(8).set_description("cache size for (T)");
  add_input<int>("ccsdt_tilesize").set_default(32).set_description("MO tilesize for (T)");

  // add_result<SystemData>("CCSD System Data").set_description("CCSD System Data");
}

template<typename T>
TEMPLATED_MODULE_RUN(CCSD, T) {
  using Matrix = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;

  const auto& [bra, H_e, ket] = ce_pt::unwrap_inputs(inputs);
  const auto& i               = bra.basis_set().occupied_orbitals();
  const auto& a               = bra.basis_set().virtual_orbitals();

  const auto& bra_aos = bra.basis_set().occupied_orbitals().from_space();

  auto& f_mod = submods.at("Fock Builder");

  const auto& f_hat     = bra.basis_set().fock_operator();
  const auto& f_wrapper = f_mod.run_as<f_pt>(bra_aos, f_hat, bra_aos);

  auto& ee_mod    = submods.at("Electronic Energy");
  auto  hf_energy = ee_mod.run_as<ee_pt>(bra, H_e, ket);

  const auto& C_occ  = i.C();
  const auto& C_virt = a.C();

  Matrix f_ao_eig   = tensor_wrapper_to_eigen(f_wrapper);
  Matrix C_occ_eig  = tensor_wrapper_to_eigen(C_occ);
  Matrix C_virt_eig = tensor_wrapper_to_eigen(C_virt);

  auto nwx_shells     = bra.basis_set().occupied_orbitals().from_space().basis_set();
  auto nbf            = nwx_shells.n_aos();
  auto basis_set_name = nwx_shells[0].basis_set_name();

  libint2::BasisSet li_shells = ccsd_make_basis(nwx_shells);

  // Setup system data
  OptionsMap options_map;
  SystemData sys_data(options_map, "restricted");
  sys_data.nbf_orig = sys_data.nbf = C_occ_eig.rows();
  sys_data.n_occ_alpha = sys_data.n_occ_beta = C_occ_eig.cols();
  sys_data.n_vir_alpha = sys_data.n_vir_beta = C_virt_eig.cols();
  // Check: Compute Northo
  sys_data.n_lindep = sys_data.nbf_orig - sys_data.n_occ_alpha - sys_data.n_vir_alpha;
  sys_data.nbf -= sys_data.n_lindep;
  
  sys_data.options_map.cd_options.diagtol          = inputs.at("diagtol").value<double>();
  sys_data.options_map.cd_options.write_cv         = inputs.at("write_cv").value<bool>();
  sys_data.options_map.cd_options.write_vcount     = inputs.at("write_vcount").value<int>();
  sys_data.options_map.cd_options.max_cvecs_factor = inputs.at("max_cvecs_factor").value<int>();
  sys_data.options_map.ccsd_options.debug          = inputs.at("debug").value<bool>();
  sys_data.options_map.ccsd_options.printtol       = inputs.at("printtol").value<double>();
  sys_data.options_map.ccsd_options.threshold      = inputs.at("threshold").value<double>();
  sys_data.options_map.ccsd_options.force_tilesize = inputs.at("force_tilesize").value<bool>();
  sys_data.options_map.ccsd_options.tilesize       = inputs.at("tilesize").value<int>();
  sys_data.options_map.ccsd_options.itilesize      = inputs.at("itilesize").value<int>();
  sys_data.options_map.ccsd_options.ndiis          = inputs.at("ndiis").value<int>();
  sys_data.options_map.ccsd_options.lshift         = inputs.at("lshift").value<double>();
  sys_data.options_map.ccsd_options.ccsd_maxiter   = inputs.at("ccsd_maxiter").value<int>();
  sys_data.options_map.ccsd_options.freeze_core    = inputs.at("freeze_core").value<int>();
  sys_data.options_map.ccsd_options.freeze_virtual = inputs.at("freeze_virtual").value<int>();
  sys_data.options_map.ccsd_options.balance_tiles  = inputs.at("balance_tiles").value<bool>();
  sys_data.options_map.ccsd_options.profile_ccsd   = inputs.at("profile_ccsd").value<bool>();
  sys_data.options_map.ccsd_options.writet         = inputs.at("writet").value<bool>();
  sys_data.options_map.ccsd_options.writev         = inputs.at("writev").value<bool>();
  sys_data.options_map.ccsd_options.writet_iter    = inputs.at("writet_iter").value<int>();
  sys_data.options_map.ccsd_options.readt          = inputs.at("readt").value<bool>();
  sys_data.options_map.ccsd_options.computeTData   = inputs.at("computeTData").value<bool>();

  sys_data.options_map.ccsd_options.cache_size     = inputs.at("cache_size").value<int>();
  sys_data.options_map.ccsd_options.ccsdt_tilesize = inputs.at("ccsdt_tilesize").value<int>();

  // TODO: get
  sys_data.input_molecule                 = "h2o";
  sys_data.output_file_prefix             = sys_data.input_molecule;
  sys_data.basis                          = basis_set_name.value();
  sys_data.options_map.ccsd_options.basis = sys_data.basis;
  sys_data.scf_energy                     = hf_energy;

  sys_data.update();

  // MPI is already initialized by TA, can we terminate madworld explicitly at
  // this point ?
  if(!GA_Initialized()) {
    GA_Initialize();
    (void) ProcGroup::self_ga_pgroup(true);
  }

  ProcGroup        pg = ProcGroup::create_coll(GA_MPI_Comm());
  ExecutionContext ec{pg, DistributionKind::nw, MemoryManagerKind::ga};
  auto             rank = ec.pg().rank();
  std::cout << std::defaultfloat;

  // Setup SCF data in TAMM
  const auto scf_conv = true;

  std::vector<size_t>     shell_tile_map; //not used
  std::vector<tamm::Tile> AO_tiles, AO_opttiles;

  auto       ao_ts    = static_cast<tamm::Tile>(std::ceil(sys_data.nbf_orig * 0.05));
  if(sys_data.AO_tilesize > ao_ts) ao_ts = sys_data.AO_tilesize;
  else sys_data.AO_tilesize = ao_ts;

  std::tie(shell_tile_map, AO_tiles, AO_opttiles) =
    compute_AO_tiles(ec, sys_data, li_shells);
  
  TiledIndexSpace AO_opt{IndexSpace{range(0, (size_t) (sys_data.nbf))}, AO_opttiles};
  TiledIndexSpace AO_ortho{IndexSpace{range(0, (size_t) (sys_data.nbf_orig))}, ao_ts};

  tamm::Tensor<T>     C_AO{AO_opt, AO_ortho}, F_AO{AO_opt, AO_opt};
  tamm::Tensor<T>     C_beta_AO{AO_opt, AO_ortho}, F_beta_AO{AO_opt, AO_opt}; // not used
  tamm::Tensor<T>::allocate(&ec, C_AO, F_AO);
  if(rank == 0) {
    Matrix C_AO_eig(sys_data.nbf_orig, sys_data.nbf); // NxNortho
    C_AO_eig << C_occ_eig, C_virt_eig;
    tamm::eigen_to_tamm_tensor(C_AO, C_AO_eig);
    tamm::eigen_to_tamm_tensor(F_AO, f_ao_eig);
  }
  C_occ_eig.resize(0, 0);
  C_virt_eig.resize(0, 0);
  f_ao_eig.resize(0, 0);
  ec.pg().barrier();

  // Start CCSD
  CCSDOptions& ccsd_options = sys_data.options_map.ccsd_options;
  const auto   debug        = ccsd_options.debug;
  if(rank == 0) {
    ccsd_options.print();
    std::cout << std::endl
              << "#occupied, #virtual = " << sys_data.nocc << ", " << sys_data.nvir << std::endl;
  }

  auto [MO, total_orbitals] = setupMOIS(sys_data);

  std::string out_fp       = sys_data.output_file_prefix + "." + ccsd_options.basis;
  std::string files_dir    = out_fp + "_files/" + sys_data.scf_type_string;
  std::string files_prefix = /*out_fp;*/ files_dir + "/" + out_fp;
  std::string f1file       = files_prefix + ".f1_mo";
  std::string t1file       = files_prefix + ".t1amp";
  std::string t2file       = files_prefix + ".t2amp";
  std::string v2file       = files_prefix + ".cholv2";
  std::string cholfile     = files_prefix + ".cholcount";
  std::string ccsdstatus   = files_prefix + ".ccsdstatus";

  const bool is_rhf = sys_data.is_restricted;
  if(!fs::exists(files_dir)) fs::create_directories(files_dir);

  bool ccsd_restart = ccsd_options.readt || ((fs::exists(t1file) && fs::exists(t2file) &&
                                              fs::exists(f1file) && fs::exists(v2file)));

  auto [cholVpr, d_f1_, lcao, chol_count, max_cvecs, CI] =
    cd_svd_driver<T>(sys_data, ec, MO, AO_opt, C_AO, F_AO, C_beta_AO, F_beta_AO, li_shells,
                        shell_tile_map, ccsd_restart, cholfile);

  Tensor<T> d_f1 = d_f1_;

  free_tensors(lcao);

  if(ccsd_options.writev) ccsd_options.writet = true;

  TiledIndexSpace N = MO("all");

  std::vector<T>         p_evl_sorted;
  Tensor<T>              d_r1, d_r2, d_t1, d_t2;
  std::vector<Tensor<T>> d_r1s, d_r2s, d_t1s, d_t2s;

  if(is_rhf)
    std::tie(p_evl_sorted, d_t1, d_t2, d_r1, d_r2, d_r1s, d_r2s, d_t1s, d_t2s) = setupTensors_cs(
      ec, MO, d_f1, ccsd_options.ndiis, ccsd_restart && fs::exists(ccsdstatus) && scf_conv);
  else
    std::tie(p_evl_sorted, d_t1, d_t2, d_r1, d_r2, d_r1s, d_r2s, d_t1s, d_t2s) = setupTensors(
      ec, MO, d_f1, ccsd_options.ndiis, ccsd_restart && fs::exists(ccsdstatus) && scf_conv);

  if(ccsd_restart) {
    read_from_disk(d_f1, f1file);
    if(fs::exists(t1file) && fs::exists(t2file)) {
      read_from_disk(d_t1, t1file);
      read_from_disk(d_t2, t2file);
    }
    read_from_disk(cholVpr, v2file);
    ec.pg().barrier();
    p_evl_sorted = tamm::diagonal(d_f1);
  }

  else if(ccsd_options.writet) {
    // fs::remove_all(files_dir);
    if(!fs::exists(files_dir)) fs::create_directories(files_dir);

    write_to_disk(d_f1, f1file);
    write_to_disk(cholVpr, v2file);

    if(rank == 0) {
      std::ofstream out(cholfile, std::ios::out);
      if(!out) cerr << "Error opening file " << cholfile << endl;
      out << chol_count << std::endl;
      out.close();
    }
  }

  if(rank == 0 && debug) {
    print_vector(p_evl_sorted, files_prefix + ".eigen_values.txt");
    cout << "Eigen values written to file: " << files_prefix + ".eigen_values.txt" << endl << endl;
  }

  ec.pg().barrier();

  auto cc_t1 = std::chrono::high_resolution_clock::now();

  ccsd_restart = ccsd_restart && fs::exists(ccsdstatus) && scf_conv;

  std::string fullV2file = files_prefix + ".fullV2";
  t1file                 = files_prefix + ".fullT1amp";
  t2file                 = files_prefix + ".fullT2amp";

  bool computeTData = ccsd_options.computeTData;
  if(ccsd_options.writev)
    computeTData = computeTData && !fs::exists(fullV2file) && !fs::exists(t1file) &&
                   !fs::exists(t2file);

  if(computeTData && is_rhf) setup_full_t1t2(ec, MO, dt1_full, dt2_full);

  double residual = 0, corr_energy = 0;

  if(is_rhf)
    std::tie(residual, corr_energy) =
      cd_ccsd_cs_driver<T>(sys_data, ec, MO, CI, d_t1, d_t2, d_f1, d_r1, d_r2, d_r1s, d_r2s, d_t1s,
                           d_t2s, p_evl_sorted, cholVpr, ccsd_restart, files_prefix, computeTData);

  if(computeTData && is_rhf) {
    if(ccsd_options.writev) {
      write_to_disk(dt1_full, t1file);
      write_to_disk(dt2_full, t2file);
      free_tensors(dt1_full, dt2_full);
    }
  }

  ccsd_stats(ec, sys_data.scf_energy, residual, corr_energy, ccsd_options.threshold);

  if(ccsd_options.writet && !fs::exists(ccsdstatus)) {
    // write_to_disk(d_t1,t1file);
    // write_to_disk(d_t2,t2file);
    if(rank == 0) {
      std::ofstream out(ccsdstatus, std::ios::out);
      if(!out) cerr << "Error opening file " << ccsdstatus << endl;
      out << 1 << std::endl;
      out.close();
    }
  }

  auto   cc_t2 = std::chrono::high_resolution_clock::now();
  double ccsd_time =
    std::chrono::duration_cast<std::chrono::duration<double>>((cc_t2 - cc_t1)).count();
  if(rank == 0) {
    if(is_rhf)
      std::cout << std::endl
                << "Time taken for Closed Shell Cholesky CCSD: " << ccsd_time << " secs"
                << std::endl;
    else
      std::cout << std::endl
                << "Time taken for Open Shell Cholesky CCSD: " << ccsd_time << " secs" << std::endl;
  }

  double printtol = ccsd_options.printtol;
  if(rank == 0 && debug) {
    std::cout << std::endl << "Threshold for printing amplitudes set to: " << printtol << std::endl;
    std::cout << "T1, T2 amplitudes written to files: " << files_prefix + ".print_t1amp.txt"
              << ", " << files_prefix + ".print_t2amp.txt" << std::endl
              << std::endl;
    print_max_above_threshold(d_t1, printtol, files_prefix + ".print_t1amp.txt");
    print_max_above_threshold(d_t2, printtol, files_prefix + ".print_t2amp.txt");
  }

  if(!ccsd_restart) {
    free_tensors(d_r1, d_r2);
    free_vec_tensors(d_r1s, d_r2s, d_t1s, d_t2s);
  }

  if(is_rhf) free_tensors(d_t1, d_t2);
  ec.flush_and_sync();

#if 1

  // double ccsdt_s1_t1_GetTime;
  // double ccsdt_s1_v2_GetTime;
  // double ccsd_t_data_per_rank;
  // double ccsdt_d1_t2_GetTime;
  // double ccsdt_d1_v2_GetTime;
  // double ccsdt_d2_t2_GetTime;
  // double ccsdt_d2_v2_GetTime;

  // else { //skip ccsd
  //     d_f1 = {{N,N},{1,1}};
  //     Tensor<T>::allocate(&ec,d_f1);
  // }

  auto [MO1, total_orbitals1] = setupMOIS(sys_data, true);
  TiledIndexSpace N1          = MO1("all");
  TiledIndexSpace O1          = MO1("occ");
  TiledIndexSpace V1          = MO1("virt");

  // Tensor<T> d_v2{{N,N,N,N},{2,2}};
  // Tensor<T> t_d_f1{{N1,N1},{1,1}};
  Tensor<T>    t_d_t1{{V1, O1}, {1, 1}};
  Tensor<T>    t_d_t2{{V1, V1, O1, O1}, {2, 2}};
  Tensor<T>    t_d_cv2{{N1, N1, CI}, {1, 1}};
  V2Tensors<T> v2tensors({"ijab", "ijka", "iabc"});

  bool skip_ccsd{false};

  T            ccsd_t_mem{};
  const double gib   = (1024 * 1024 * 1024.0);
  const double Osize = MO("occ").max_num_indices();
  const double Vsize = MO("virt").max_num_indices();
  // const double Nsize = N.max_num_indices();
  // const double cind_size = CI.max_num_indices();

  ccsd_t_mem = sum_tensor_sizes(d_f1, t_d_t1, t_d_t2);
  if(!skip_ccsd) {
    // auto v2_setup_mem = sum_tensor_sizes(d_f1,t_d_v2,t_d_cv2);
    // auto cv2_retile = (Nsize*Nsize*cind_size*8)/gib + sum_tensor_sizes(d_f1,cholVpr,t_d_cv2);
    if(is_rhf)
      ccsd_t_mem += sum_tensor_sizes(dt1_full, dt2_full);
    else
      ccsd_t_mem += sum_tensor_sizes(d_t1, d_t2);

    // retiling allocates full GA versions of the t1,t2 tensors.
    ccsd_t_mem += (Osize * Vsize + Vsize * Vsize * Osize * Osize) * 8 / gib;
  }

  // const auto ccsd_t_mem_old = ccsd_t_mem + sum_tensor_sizes(t_d_v2);
  ccsd_t_mem += v2tensors.tensor_sizes(MO1);

  Index noab           = MO1("occ").num_tiles();
  Index nvab           = MO1("virt").num_tiles();
  Index cache_size     = ccsd_options.cache_size;
  auto  ccsdt_tilesize = ccsd_options.ccsdt_tilesize;
  auto  ex_hw          = ec.exhw();

  {
    Index noa    = MO1("occ_alpha").num_tiles();
    Index nva    = MO1("virt_alpha").num_tiles();
    auto  nranks = ec.pg().size().value();

    auto                mo_tiles = MO1.input_tile_sizes();
    std::vector<size_t> k_range;
    for(auto x: mo_tiles) k_range.push_back(x);
    size_t max_pdim = 0;
    size_t max_hdim = 0;
    for(size_t t_p4b = noab; t_p4b < noab + nvab; t_p4b++)
      max_pdim = std::max(max_pdim, k_range[t_p4b]);
    for(size_t t_h1b = 0; t_h1b < noab; t_h1b++) max_hdim = std::max(max_hdim, k_range[t_h1b]);

    size_t max_d1_kernels_pertask = 9 * noa;
    size_t max_d2_kernels_pertask = 9 * nva;
    size_t size_T_s1_t1           = 9 * (max_pdim) * (max_hdim);
    size_t size_T_s1_v2           = 9 * (max_pdim * max_pdim) * (max_hdim * max_hdim);
    size_t size_T_d1_t2 = max_d1_kernels_pertask * (max_pdim * max_pdim) * (max_hdim * max_hdim);
    size_t size_T_d1_v2 = max_d1_kernels_pertask * (max_pdim) * (max_hdim * max_hdim * max_hdim);
    size_t size_T_d2_t2 = max_d2_kernels_pertask * (max_pdim * max_pdim) * (max_hdim * max_hdim);
    size_t size_T_d2_v2 = max_d2_kernels_pertask * (max_pdim * max_pdim * max_pdim) * (max_hdim);

    double extra_buf_mem_per_rank =
      size_T_s1_t1 + size_T_s1_v2 + size_T_d1_t2 + size_T_d1_v2 + size_T_d2_t2 + size_T_d2_v2;
    extra_buf_mem_per_rank     = extra_buf_mem_per_rank * 8 / gib;
    double total_extra_buf_mem = extra_buf_mem_per_rank * nranks;

    size_t cache_buf_size =
      ccsdt_tilesize * ccsdt_tilesize * ccsdt_tilesize * ccsdt_tilesize * 8; // bytes
    double cache_mem_per_rank =
      (ccsdt_tilesize * ccsdt_tilesize * 8 + cache_buf_size) * cache_size; // s1 t1+v2
    cache_mem_per_rank += (noab + nvab) * 2 * cache_size * cache_buf_size; // d1,d2 t2+v2
    cache_mem_per_rank     = cache_mem_per_rank / gib;
    double total_cache_mem = cache_mem_per_rank * nranks; // GiB

    double total_ccsd_t_mem = ccsd_t_mem + total_extra_buf_mem + total_cache_mem;
    if(rank == 0) {
      std::cout << std::string(70, '-') << std::fixed << std::setprecision(2) << std::endl;
      std::cout << "Total CPU memory required for (T) calculation = " << total_ccsd_t_mem << " GiB"
                << std::endl;
      std::cout << " -- memory required for the input tensors: " << ccsd_t_mem << " GiB"
                << std::endl;
      std::cout << " -- memory required for intermediate buffers: " << total_extra_buf_mem << " GiB"
                << std::endl;
      std::string cache_msg = " -- memory required for caching t1,t2,v2 blocks";
      if(total_cache_mem > (ccsd_t_mem + total_extra_buf_mem) / 2.0)
        cache_msg += " (set cache_size option in the input file to a lower value to reduce this "
                     "memory requirement further)";
      std::cout << cache_msg << ": " << total_cache_mem << " GiB" << std::endl;
      // std::cout << "***** old memory requirement was "
      //           << ccsd_t_mem_old + total_extra_buf_mem + total_cache_mem
      //           << " GiB (old v2 = " << sum_tensor_sizes(t_d_v2)
      //           << " GiB, new v2 = " << v2tensors.tensor_sizes(MO1) << " GiB)" << std::endl;
      std::cout << std::string(70, '-') << std::endl;
    }
    check_memory_requirements(ec, total_ccsd_t_mem);
  }

  if(computeTData && !skip_ccsd) {
    Tensor<T>::allocate(&ec, t_d_cv2);
    retile_tamm_tensor(cholVpr, t_d_cv2, "CholV2");
    free_tensors(cholVpr);

    v2tensors = setupV2Tensors<T>(ec, t_d_cv2, ex_hw, v2tensors.get_blocks());
    if(ccsd_options.writev) {
      v2tensors.write_to_disk(files_prefix);
      v2tensors.deallocate();
    }
    free_tensors(t_d_cv2);
  }

  double energy1 = 0, energy2 = 0;

  if(rank == 0) {
    auto mo_tiles = MO.input_tile_sizes();
    cout << endl << "CCSD MO Tiles = " << mo_tiles << endl;
  }

  Tensor<T>::allocate(&ec, t_d_t1, t_d_t2);
  if(skip_ccsd) v2tensors.allocate(ec, MO1);

  bool ccsd_t_restart = fs::exists(t1file) && fs::exists(t2file) && fs::exists(f1file) &&
                        v2tensors.exist_on_disk(files_prefix);

  if(!ccsd_t_restart && !skip_ccsd) {
    if(!is_rhf) {
      dt1_full = d_t1;
      dt2_full = d_t2;
    }
    if(rank == 0) { cout << endl << "Retile T1,T2 tensors ... " << endl; }

    Scheduler{ec}(t_d_t1() = 0)(t_d_t2() = 0).execute();

    TiledIndexSpace O = MO("occ");
    TiledIndexSpace V = MO("virt");

    if(ccsd_options.writev) {
      Tensor<T> wd_t1{{V, O}, {1, 1}};
      Tensor<T> wd_t2{{V, V, O, O}, {2, 2}};

      read_from_disk(t_d_t1, t1file, false, wd_t1);
      read_from_disk(t_d_t2, t2file, false, wd_t2);

      ec.pg().barrier();
      write_to_disk(t_d_t1, t1file);
      write_to_disk(t_d_t2, t2file);
    }

    else {
      retile_tamm_tensor(dt1_full, t_d_t1);
      retile_tamm_tensor(dt2_full, t_d_t2);
      if(is_rhf) free_tensors(dt1_full, dt2_full);
    }
  }
  else if(ccsd_options.writev && !skip_ccsd) {
    // read_from_disk(t_d_f1,f1file);
    read_from_disk(t_d_t1, t1file);
    read_from_disk(t_d_t2, t2file);
    v2tensors.read_from_disk(files_prefix);
  }

  if(!is_rhf && !skip_ccsd) free_tensors(d_t1, d_t2);

  p_evl_sorted = tamm::diagonal(d_f1);

  // cc_t1 = std::chrono::high_resolution_clock::now();

  bool is_restricted = is_rhf;

  if(rank == 0) {
    if(is_restricted)
      cout << endl << "Running Closed Shell CCSD(T) calculation" << endl;
    else
      cout << endl << "Running Open Shell CCSD(T) calculation" << endl;
  }

  bool                            seq_h3b = true;
  LRUCache<Index, std::vector<T>> cache_s1t{cache_size};
  LRUCache<Index, std::vector<T>> cache_s1v{cache_size};
  LRUCache<Index, std::vector<T>> cache_d1t{cache_size * noab};
  LRUCache<Index, std::vector<T>> cache_d1v{cache_size * noab};
  LRUCache<Index, std::vector<T>> cache_d2t{cache_size * nvab};
  LRUCache<Index, std::vector<T>> cache_d2v{cache_size * nvab};

  if(rank == 0 && seq_h3b) cout << "running seq h3b loop variant..." << endl;

  std::vector<int> k_spin;
  for(tamm::Index x = 0; x < noab / 2; x++) k_spin.push_back(1);
  for(tamm::Index x = noab / 2; x < noab; x++) k_spin.push_back(2);
  for(tamm::Index x = 0; x < nvab / 2; x++) k_spin.push_back(1);
  for(tamm::Index x = nvab / 2; x < nvab; x++) k_spin.push_back(2);

  double ccsd_t_time = 0, total_t_time = 0;
  // cc_t1 = std::chrono::high_resolution_clock::now();
  std::tie(energy1, energy2, ccsd_t_time, total_t_time) = ccsd_t_fused_driver_new<T>(
    sys_data, ec, k_spin, MO1, t_d_t1, t_d_t2, v2tensors, p_evl_sorted, hf_energy + corr_energy,
    is_restricted, cache_s1t, cache_s1v, cache_d1t, cache_d1v, cache_d2t, cache_d2v, seq_h3b);

  energy1 = ec.pg().reduce(&energy1, ReduceOp::sum, 0);
  energy2 = ec.pg().reduce(&energy2, ReduceOp::sum, 0);

  if(rank == 0 && !skip_ccsd) {
    std::cout.precision(15);
    cout << "CCSD[T] correction energy / hartree  = " << energy1 << endl;
    cout << "CCSD[T] correlation energy / hartree = " << corr_energy + energy1 << endl;
    cout << "CCSD[T] total energy / hartree       = " << hf_energy + corr_energy + energy1 << endl;

    cout << "CCSD(T) correction energy / hartree  = " << energy2 << endl;
    cout << "CCSD(T) correlation energy / hartree = " << corr_energy + energy2 << endl;
    cout << "CCSD(T) total energy / hartree       = " << hf_energy + corr_energy + energy2 << endl;

    sys_data.results["output"]["CCSD(T)"]["[T]Energies"]["correction"]  = energy1;
    sys_data.results["output"]["CCSD(T)"]["[T]Energies"]["correlation"] = corr_energy + energy1;
    sys_data.results["output"]["CCSD(T)"]["[T]Energies"]["total"] =
      hf_energy + corr_energy + energy1;
    sys_data.results["output"]["CCSD(T)"]["(T)Energies"]["correction"]  = energy2;
    sys_data.results["output"]["CCSD(T)"]["(T)Energies"]["correlation"] = corr_energy + energy2;
    sys_data.results["output"]["CCSD(T)"]["(T)Energies"]["total"] =
      hf_energy + corr_energy + energy2;
  }

  long double total_num_ops = 0;
  //
  if(rank == 0) {
    // std::cout << "--------------------------------------------------------------------" <<
    // std::endl;
    ccsd_t_fused_driver_calculator_ops<T>(sys_data, ec, k_spin, MO1, p_evl_sorted,
                                          hf_energy + corr_energy, is_restricted, total_num_ops,
                                          seq_h3b);
    // std::cout << "--------------------------------------------------------------------" <<
    // std::endl;
  }

  ec.pg().barrier();

  auto nranks = ec.pg().size().value();

  auto print_profile_stats = [&](const std::string& timer_type, const double g_tval,
                                 const double tval_min, const double tval_max) {
    const double tval = g_tval / nranks;
    std::cout.precision(3);
    std::cout << "   -> " << timer_type << ": " << tval << "s (" << tval * 100.0 / total_t_time
              << "%), (min,max) = (" << tval_min << "," << tval_max << ")" << std::endl;
  };

  auto comm_stats = [&](const std::string& timer_type, const double ctime) {
    double g_getTime     = ec.pg().reduce(&ctime, ReduceOp::sum, 0);
    double g_min_getTime = ec.pg().reduce(&ctime, ReduceOp::min, 0);
    double g_max_getTime = ec.pg().reduce(&ctime, ReduceOp::max, 0);

    if(rank == 0) print_profile_stats(timer_type, g_getTime, g_min_getTime, g_max_getTime);
    return g_getTime / nranks;
  };

  if(rank == 0) {
    std::cout << std::endl << "------CCSD(T) Performance------" << std::endl;
    std::cout << "Total CCSD(T) Time: " << total_t_time << std::endl;
  }
  ccsd_t_time = comm_stats("CCSD(T) Avg. Work Time", ccsd_t_time);
  if(rank == 0) {
    const double n_gflops = total_num_ops / (total_t_time * 1e9);
    const double load_imb = (1.0 - ccsd_t_time / total_t_time);
    std::cout << std::scientific << "   -> Total Number of Operations: " << total_num_ops
              << std::endl;
    std::cout << std::fixed << "   -> GFLOPS: " << n_gflops << std::endl;
    std::cout << std::fixed << "   -> Load imbalance: " << load_imb << std::endl;

    sys_data.results["output"]["CCSD(T)"]["performance"]["total_time"]     = total_t_time;
    sys_data.results["output"]["CCSD(T)"]["performance"]["gflops"]         = n_gflops;
    sys_data.results["output"]["CCSD(T)"]["performance"]["total_num_ops"]  = total_num_ops;
    sys_data.results["output"]["CCSD(T)"]["performance"]["load_imbalance"] = load_imb;
    write_json_data(sys_data, "CCSD_T");
  }

  // comm_stats("S1-T1 GetTime", ccsdt_s1_t1_GetTime);
  // comm_stats("S1-V2 GetTime", ccsdt_s1_v2_GetTime);
  // comm_stats("D1-T2 GetTime", ccsdt_d1_t2_GetTime);
  // comm_stats("D1-V2 GetTime", ccsdt_d1_v2_GetTime);
  // comm_stats("D2-T2 GetTime", ccsdt_d2_t2_GetTime);
  // comm_stats("D2-V2 GetTime", ccsdt_d2_v2_GetTime);

  // ccsd_t_data_per_rank          = (ccsd_t_data_per_rank * 8.0) / (1024 * 1024.0 * 1024); // GB
  // double g_ccsd_t_data_per_rank = ec.pg().reduce(&ccsd_t_data_per_rank, ReduceOp::sum, 0);
  // if(rank == 0)
  //   std::cout << "   -> Data Transfer (GB): " << g_ccsd_t_data_per_rank / nranks << std::endl;

  ec.pg().barrier();

  free_tensors(t_d_t1, t_d_t2, d_f1);
  v2tensors.deallocate();

  ec.flush_and_sync();
#endif

  // GA Terminate
  // GA_Terminate();

  auto rv = results();
  rv      = ce_pt::wrap_results(rv, corr_energy);
  rv      = ccsd_pt::wrap_results(rv, d_f1, d_t1, d_t2);
  // rv.at("CCSD System Data").change(sys_data);
  return rv;
}

// Instantiations
template class CCSD<double>;

} // namespace ccsd