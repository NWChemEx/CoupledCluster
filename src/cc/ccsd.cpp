#include "cd_ccsd_cs_ann.hpp"
#include "energy.hpp"
#include <libint2.hpp>
#include <simde/simde.hpp>

namespace ccsd {

using f_pt   = simde::Fock;
using eri_pt = simde::TransformedERI4;
using ee_pt  = simde::CanonicalElectronicEnergy;

using ce_pt =
  simde::CorrelationEnergy<simde::type::canonical_reference, simde::type::canonical_reference>;

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
    centers_t centers{};
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
    for(auto center_i = 0; center_i < bs.size(); ++center_i) {
        /// Add current center to atoms list
        auto& center = bs[center_i];
        centers.push_back(
          atom_ctor(center_i, center.x(), center.y(), center.z()));

        /// Gather shells for this center and add them to element_bases
        atom_bases_t atom_bases{};
        for(const auto&& shelli : center) {
            const auto nprims = shelli.n_unique_primitives();
            const auto prim0  = shelli.unique_primitive(0);
            const auto primN  = shelli.unique_primitive(nprims - 1);
            const bool pure   = shelli.pure() == chemist::ShellType::pure;
            const int l       = shelli.l();

            svec_d_t alphas(&prim0.exponent(), &primN.exponent() + 1);
            svec_d_t coefs(&prim0.coefficient(), &primN.coefficient() + 1);
            conts_t conts{cont_t{l, pure, coefs}};
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
  satisfies_property_type<simde::TotalCanonicalEnergy>();

  add_submodule<f_pt>("Fock Builder");
  add_submodule<eri_pt>("Transformed ERIs");
  add_submodule<ee_pt>("Electronic Energy");

  add_input<bool>("debug").set_default(false).set_description("Debugging flag");

  add_input<double>("diagtol").set_default(1.0e-6).set_description(
    "Cholesky Decomposition Threshold");

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
    .set_default(false)
    .set_description("Compute and write data needed for (T) calculation to disk");

  add_result<tamm::Tensor<T>>("MO Fock Matrix").set_description("The Fock matrix in MO");
}

template<typename T>
TEMPLATED_MODULE_RUN(CCSD, T) {
  using Matrix = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;

  const auto& [bra, H_e, ket] = ce_pt::unwrap_inputs(inputs);
  const auto& i                 = bra.basis_set().occupied_orbitals();
  const auto& a                 = bra.basis_set().virtual_orbitals();

  const auto& bra_aos = bra.basis_set().occupied_orbitals().from_space();

  auto& f_mod = submods.at("Fock Builder");

  const auto& f_hat       = bra.basis_set().fock_operator();
  const auto& [f_wrapper] = f_mod.run_as<f_pt>(bra_aos, f_hat, bra_aos);

  auto&                        ee_mod = submods.at("Electronic Energy");
  auto [scf_energy]                   = ee_mod.run_as<ee_pt>(bra, H_e, ket);

  const auto& C_occ  = i.C();
  const auto& C_virt = a.C();

  Matrix f_ao_eig   = tensor_wrapper_to_eigen(f_wrapper);
  Matrix C_occ_eig  = tensor_wrapper_to_eigen(C_occ);
  Matrix C_virt_eig = tensor_wrapper_to_eigen(C_virt);

  auto              nwx_shells     = bra.basis_set().occupied_orbitals().from_space().basis_set();
  auto              nbf            = nwx_shells.n_aos();
  const std::string basis_set_name = nwx_shells[0].basis_set_name();

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

  // TODO: get
  sys_data.input_molecule     = "h2o";
  sys_data.output_file_prefix = sys_data.input_molecule;
  sys_data.basis = sys_data.options_map.ccsd_options.basis = basis_set_name;
  sys_data.scf_energy                                      = scf_energy;

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
  auto       ao_ts    = static_cast<tamm::Tile>(std::ceil(sys_data.nbf_orig * 0.05));
  if(sys_data.AO_tilesize > ao_ts) ao_ts = sys_data.AO_tilesize;
  TiledIndexSpace AO_opt{IndexSpace{range(0, (size_t) (sys_data.nbf))}, ao_ts};
  TiledIndexSpace AO_ortho{IndexSpace{range(0, (size_t) (sys_data.nbf_orig))}, ao_ts};

  std::vector<size_t> shell_tile_map; // not used
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
    cd_svd_ga_driver<T>(sys_data, ec, MO, AO_opt, C_AO, F_AO, C_beta_AO, F_beta_AO, li_shells,
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

  ExecutionHW ex_hw = ExecutionHW::CPU;
#ifdef USE_DPCPP
  ex_hw = ExecutionHW::GPU;
#endif
#ifdef USE_TALSH
  ex_hw              = ExecutionHW::GPU;
  const bool has_gpu = ec.has_gpu();
  TALSH      talsh_instance;
  if(has_gpu) talsh_instance.initialize(ec.gpu_devid(), rank.value());
#endif

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

  // if(is_rhf)
  //     free_tensors(d_t1, d_t2, d_f1);
  // else
  //     free_tensors(d_f1);

  Tensor<T> d_v2;
  if(computeTData) {
    d_v2 = setupV2<T>(ec, MO, CI, cholVpr, chol_count, ex_hw);
    if(ccsd_options.writev) {
      write_to_disk(d_v2, fullV2file, true);
      Tensor<T>::deallocate(d_v2);
    }
  }

  free_tensors(cholVpr);

#ifdef USE_TALSH
  // talshStats();
  if(has_gpu) talsh_instance.shutdown();
#endif

  if(computeTData) {
    if(!is_rhf) {
      dt1_full = d_t1;
      dt2_full = d_t2;
    }
    if(rank == 0) { cout << endl << "Retile T1,T2,V2 ... " << endl; }

    auto [MO1, total_orbitals1] = setupMOIS(sys_data, true);
    TiledIndexSpace N1          = MO1("all");
    TiledIndexSpace O1          = MO1("occ");
    TiledIndexSpace V1          = MO1("virt");

    // Tensor<T> t_d_f1{{N1,N1},{1,1}};
    Tensor<T> t_d_t1{{V1, O1}, {1, 1}};
    Tensor<T> t_d_t2{{V1, V1, O1, O1}, {2, 2}};
    Tensor<T> t_d_v2{{N1, N1, N1, N1}, {2, 2}};
    Tensor<T>::allocate(&ec, t_d_t1, t_d_t2, t_d_v2);

    Scheduler{ec} // (t_d_f1() = 0)
    (t_d_t1() = 0)(t_d_t2() = 0)(t_d_v2() = 0)
      .execute();

    TiledIndexSpace O = MO("occ");
    TiledIndexSpace V = MO("virt");

    if(ccsd_options.writev) {
      // Tensor<T> wd_f1{{N,N},{1,1}};
      Tensor<T> wd_t1{{V, O}, {1, 1}};
      Tensor<T> wd_t2{{V, V, O, O}, {2, 2}};
      Tensor<T> wd_v2{{N, N, N, N}, {2, 2}};

      // read_from_disk(t_d_f1,f1file,false,wd_f1);
      read_from_disk(t_d_t1, t1file, false, wd_t1);
      read_from_disk(t_d_t2, t2file, false, wd_t2);
      read_from_disk(t_d_v2, fullV2file, false, wd_v2);

      ec.pg().barrier();
      // write_to_disk(t_d_f1,f1file);
      write_to_disk(t_d_t1, t1file);
      write_to_disk(t_d_t2, t2file);
      write_to_disk(t_d_v2, fullV2file);
    }

    else {
      retile_tamm_tensor(dt1_full, t_d_t1);
      retile_tamm_tensor(dt2_full, t_d_t2);
      if(is_rhf) free_tensors(dt1_full, dt2_full);
      retile_tamm_tensor(d_v2, t_d_v2, "V2");
      free_tensors(d_v2);
    }

    free_tensors(t_d_t1, t_d_t2, t_d_v2);
  }

  if(!is_rhf) free_tensors(d_t1, d_t2);

  ec.flush_and_sync();

  // GA Terminate
  // GA_Terminate();

  auto rv = results();
  rv      = ce_pt::wrap_results(rv, corr_energy);
  rv.at("MO Fock Matrix").change(d_f1);
  return rv;
}

// Instantiations
template class CCSD<double>;

} // namespace ccsd