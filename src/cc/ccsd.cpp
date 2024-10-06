// clang-format off

#include "cc_modules.hpp"
#include <libint2.hpp>
#include <simde/simde.hpp>
#include "exachem/common/chemenv.hpp"
#include "exachem/scf/scf_main.hpp"
#include "exachem/common/initialize_system_data.hpp"
#include "exachem/cc/ccsd/cd_ccsd_os_ann.hpp"
// #include "exachem/cc/ccsd_t/ccsd_t_fused_driver.hpp"

// clang-format on

namespace cc {

using energy_pt = simde::AOEnergy;

inline libint2::BasisSet make_libint_basis(const simde::type::ao_basis_set& bs) {
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

MODULE_CTOR(CCSDEnergy) {
  satisfies_property_type<energy_pt>();

  // add_submodule<f_pt>("Fock Builder");
  // add_submodule<eri_pt>("Transformed ERIs");
  // add_submodule<ee_pt>("Electronic Energy");

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

MODULE_RUN(CCSDEnergy) {

    const auto       rank = ProcGroup::world_rank();
    ProcGroup        pg   = ProcGroup::create_world_coll();
    ExecutionContext ec{pg, DistributionKind::nw, MemoryManagerKind::ga};

    const auto& [aos, cs] = energy_pt::unwrap_inputs(inputs);

    const double angstrom_to_bohr = 1.8897259878858;

    libint2::BasisSet li_shells = make_libint_basis(aos);

    ChemEnv chem_env;
    chem_env.input_file = inputs.at("molecule_name").value<std::string>();
    auto mol = cs.molecule();
    chem_env.ec_atoms.resize(mol.size());
    chem_env.atoms.resize(mol.size());

    CommonOptions& coptions = chem_env.ioptions.common_options;
    //parse common options
    coptions.geom_units = inputs.at("units").value<std::string>();
    const double convert_units = (coptions.geom_units == "angstrom") ? angstrom_to_bohr : 1.0;    
    coptions.basis = aos[0].basis_set_name().value_or("sto-3g");

    std::cout << std::endl;
    for (int i = 0; i < mol.size(); i++) {
      auto atom_i = mol[i];
      chem_env.atoms[i] = {(int)atom_i.Z(), atom_i.x() * convert_units, atom_i.y() * convert_units, atom_i.z() * convert_units};
      chem_env.ec_atoms[i].atom    = chem_env.atoms[i];
      chem_env.ec_atoms[i].esymbol = atom_i.name();    
      chem_env.ec_atoms[i].basis = coptions.basis;    

      std::cout << std::setw(3) << std::left << chem_env.ec_atoms[i].esymbol << " " << std::right << std::setw(14)
                << std::fixed << std::setprecision(10) << chem_env.atoms[i].x << " " << std::right
                << std::setw(14) << std::fixed << std::setprecision(10) << chem_env.atoms[i].y << " "
                << std::right << std::setw(14) << std::fixed << std::setprecision(10) << chem_env.atoms[i].z << "\n";      
    }

    chem_env.sys_data.input_molecule = chem_env.input_file;

    if(chem_env.ioptions.common_options.file_prefix.empty()) {
      chem_env.ioptions.common_options.file_prefix = chem_env.sys_data.input_molecule;
    }

    chem_env.sys_data.output_file_prefix =
      chem_env.ioptions.common_options.file_prefix + "." + chem_env.ioptions.common_options.basis;
    chem_env.workspace_dir = chem_env.sys_data.output_file_prefix + "_files/";

    // Set SCF options
    SCFOptions& scf = chem_env.ioptions.scf_options;
    scf.charge = inputs.at("charge").value<int>();
    scf.multiplicity = inputs.at("multiplicity").value<int>();
    scf.lshift = inputs.at("lshift").value<double>();
    scf.tol_int = inputs.at("tol_int").value<double>();
    scf.tol_sch = inputs.at("tol_sch").value<double>();
    scf.tol_lindep = inputs.at("tol_lindep").value<double>();
    scf.conve = inputs.at("conve").value<double>();
    scf.conve = inputs.at("convd").value<double>();
    scf.diis_hist = inputs.at("diis_hist").value<int>();
    scf.damp = inputs.at("damp").value<int>();
    scf.writem = inputs.at("writem").value<int>();
    scf.debug = inputs.at("debug").value<bool>();
    scf.restart = inputs.at("restart").value<bool>();
    scf.noscf = inputs.at("noscf").value<bool>();
    scf.scf_type = inputs.at("scf_type").value<std::string>();
    scf.direct_df = inputs.at("direct_df").value<bool>();

    // DFT
    scf.snK = inputs.at("snK").value<bool>();
    scf.xc_type = inputs.at("xc_type").value<std::vector<std::string>>();

    scf.xc_grid_type = inputs.at("xc_grid_type").value<std::string>();
    scf.xc_pruning_scheme = inputs.at("xc_pruning_scheme").value<std::string>();
    scf.xc_rad_quad = inputs.at("xc_rad_quad").value<std::string>();
    scf.xc_weight_scheme = inputs.at("xc_weight_scheme").value<std::string>();
    scf.xc_exec_space = inputs.at("xc_exec_space").value<std::string>();

    scf.xc_basis_tol = inputs.at("xc_basis_tol").value<double>();
    scf.xc_batch_size = inputs.at("xc_batch_size").value<int>();
    scf.xc_snK_etol = inputs.at("xc_snK_etol").value<double>();
    scf.xc_snK_ktol = inputs.at("xc_snK_ktol").value<double>();

    scf.xc_lb_kernel = inputs.at("xc_lb_kernel" ).value<std::string>();
    scf.xc_mw_kernel = inputs.at("xc_mw_kernel" ).value<std::string>();
    scf.xc_int_kernel = inputs.at("xc_int_kernel").value<std::string>();
    scf.xc_red_kernel = inputs.at("xc_red_kernel").value<std::string>();
    scf.xc_lwd_kernel = inputs.at("xc_lwd_kernel").value<std::string>();    

    IniSystemData    ini_sys_data(chem_env);
    SCFOptions& scf_options   = chem_env.ioptions.scf_options;
    chem_env.ec_basis         = ECBasis(ec, scf_options.basis, scf_options.basisfile,
                                        scf_options.gaussian_type, chem_env.atoms, chem_env.ec_atoms);
    chem_env.shells           = chem_env.ec_basis.shells;
    chem_env.sys_data.has_ecp = chem_env.ec_basis.has_ecp;

  // sys_data.options_map.cd_options.diagtol          = inputs.at("diagtol").value<double>();
  // sys_data.options_map.cd_options.itilesize        = inputs.at("itilesize").value<int>();
  // sys_data.options_map.cd_options.write_cv         = inputs.at("write_cv").value<bool>();
  // sys_data.options_map.cd_options.write_vcount     = inputs.at("write_vcount").value<int>();
  // sys_data.options_map.cd_options.max_cvecs_factor = inputs.at("max_cvecs_factor").value<int>();
  // sys_data.options_map.ccsd_options.debug          = inputs.at("debug").value<bool>();
  // sys_data.options_map.ccsd_options.printtol       = inputs.at("printtol").value<double>();
  // sys_data.options_map.ccsd_options.threshold      = inputs.at("threshold").value<double>();
  // sys_data.options_map.ccsd_options.force_tilesize = inputs.at("force_tilesize").value<bool>();
  // sys_data.options_map.ccsd_options.tilesize       = inputs.at("tilesize").value<int>();
  // sys_data.options_map.ccsd_options.ndiis          = inputs.at("ndiis").value<int>();
  // sys_data.options_map.ccsd_options.lshift         = inputs.at("lshift").value<double>();
  // sys_data.options_map.ccsd_options.ccsd_maxiter   = inputs.at("ccsd_maxiter").value<int>();
  // sys_data.options_map.ccsd_options.freeze_core    = inputs.at("freeze_core").value<int>();
  // sys_data.options_map.ccsd_options.freeze_virtual = inputs.at("freeze_virtual").value<int>();
  // sys_data.options_map.ccsd_options.balance_tiles  = inputs.at("balance_tiles").value<bool>();
  // sys_data.options_map.ccsd_options.profile_ccsd   = inputs.at("profile_ccsd").value<bool>();
  // sys_data.options_map.ccsd_options.writet         = inputs.at("writet").value<bool>();
  // sys_data.options_map.ccsd_options.writev         = inputs.at("writev").value<bool>();
  // sys_data.options_map.ccsd_options.writet_iter    = inputs.at("writet_iter").value<int>();
  // sys_data.options_map.ccsd_options.readt          = inputs.at("readt").value<bool>();
  // sys_data.options_map.ccsd_options.computeTData   = inputs.at("computeTData").value<bool>();

  // sys_data.options_map.ccsd_options.cache_size     = inputs.at("cache_size").value<int>();
  // sys_data.options_map.ccsd_options.ccsdt_tilesize = inputs.at("ccsdt_tilesize").value<int>();

  // exachem::scf::scf(ec, chem_env);
  exachem::cc::ccsd::cd_ccsd(ec, chem_env);

  double E0 = chem_env.hf_energy; // This is a total energy in Hartree
  auto rv = results();
  return energy_pt::wrap_results(rv, E0);

  return rv;
}

// Instantiations
// template class CCSD<double>;

} // namespace cc