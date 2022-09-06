#pragma once

#include "json_data.hpp"

using TensorType = double;

size_t nbasis(const std::vector<libint2::Shell>& shells) {
    size_t n = 0;
    for(const auto& shell : shells) n += shell.size();
    return n;
}

std::vector<size_t> map_shell_to_basis_function(
  const std::vector<libint2::Shell>& shells) {
    std::vector<size_t> result;
    result.reserve(shells.size());

    size_t n = 0;
    for(auto shell : shells) {
        result.push_back(n);
        n += shell.size();
    }

    return result;
}

std::vector<size_t> map_basis_function_to_shell(
  const std::vector<libint2::Shell>& shells) {
    std::vector<size_t> result(nbasis(shells));

    auto shell2bf = map_shell_to_basis_function(shells);
    for(size_t s1 = 0; s1 != shells.size(); ++s1) {
        auto bf1_first = shell2bf[s1]; // first basis function in this shell
        auto n1        = shells[s1].size();
        for(size_t f1 = 0; f1 != n1; ++f1) {
            const auto bf1 = f1 + bf1_first;
            result[bf1]    = s1;
        }
    }
    return result;
}

void print_bool(std::string str, bool val) {
    if(val)
        cout << str << " = true" << endl;
    else
        cout << str << " = false" << endl;
}

class Options {
public:
    Options() {
        maxiter            = 50;
        debug              = false;
        basis              = "sto-3g";
        dfbasis            = "";
        geom_units         = "bohr";
        sphcart            = "spherical";
        output_file_prefix = "";
        ext_data_path      = "";
    }

    bool debug;
    int maxiter;
    std::string basis;
    std::string dfbasis;
    std::string sphcart;
    std::string geom_units;
    std::string output_file_prefix;
    std::string ext_data_path;
    std::vector<libint2::Atom> atoms;

    void print() {
        std::cout << std::defaultfloat;
        cout << endl << "Common Options" << endl;
        cout << "{" << endl;
        cout << " maxiter    = " << maxiter << endl;
        cout << " basis      = " << basis << " ";
        cout << sphcart;
        cout << endl;
        if(!dfbasis.empty()) cout << " dfbasis    = " << dfbasis << endl;
        cout << " geom_units = " << geom_units << endl;
        print_bool(" debug     ", debug);
        if(!output_file_prefix.empty())
            cout << " output_file_prefix    = " << output_file_prefix << endl;
        cout << "}" << endl;
    }
};

class CDOptions : public Options {
public:
    CDOptions() = default;
    CDOptions(Options o) : Options(o) {
        diagtol = 1e-6;
        // At most 8*ao CholVec's. For vast majority cases, this is way
        // more than enough. For very large basis, it can be increased.
        max_cvecs_factor = 12;
    }

    double diagtol;
    int max_cvecs_factor;

    void print() {
        std::cout << std::defaultfloat;
        cout << endl << "CD Options" << endl;
        cout << "{" << endl;
        cout << " diagtol          = " << std::scientific
             << std::setprecision(2) << diagtol << endl;
        cout << " max_cvecs_factor = " << max_cvecs_factor << endl;
        print_bool(" debug           ", debug);
        cout << "}" << endl;
    }
};

class CCSDOptions : public Options {
public:
    CCSDOptions() = default;
    CCSDOptions(Options o) : Options(o) {
        printtol       = 0.05;
        threshold      = 1e-6;
        force_tilesize = false;
        tilesize       = 50;
        itilesize      = 1000;
        ndiis          = 5;
        lshift         = 0;
        ccsd_maxiter   = 50;
        freeze_core    = 0;
        freeze_virtual = 0;
        balance_tiles  = true;
        profile_ccsd   = false;

        writet       = false;
        writev       = false;
        writet_iter  = ndiis;
        readt        = false;
        computeTData = false;

        localize      = false;
        skip_dlpno    = false;
        keep_npairs   = 1;
        max_pnos      = 1;
        dlpno_dfbasis = "";
        TCutEN        = 0.97;
        TCutPNO       = 0.00;
        TCutPre       = -1.0;
        TCutPairs     = 0.00;
        TCutDO        = 1e-2;
        TCutDOij      = 1e-5;
        TCutDOPre     = 3e-2;

        ngpu           = 0;
        ccsdt_tilesize = 28;
    }

    int tilesize;
    int itilesize;
    bool force_tilesize;
    int ndiis;
    int writet_iter;
    bool readt, writet, writev, balance_tiles, computeTData;
    bool profile_ccsd;
    double lshift;
    double printtol;
    double threshold;

    int ccsd_maxiter;
    int freeze_core;
    int freeze_virtual;

    // CCSD(T)
    int ngpu;
    int ccsdt_tilesize;

    // DLPNO
    bool localize;
    bool skip_dlpno;
    int max_pnos;
    size_t keep_npairs;
    std::string dlpno_dfbasis;
    double TCutEN;
    double TCutPNO;
    double TCutPre;
    double TCutPairs;
    double TCutDO;
    double TCutDOij;
    double TCutDOPre;
    std::vector<int> doubles_opt_eqns;

    void print() {
        std::cout << std::defaultfloat;
        cout << endl << "CCSD Options" << endl;
        cout << "{" << endl;
        // if(ngpu > 0) {
        //   cout << " ngpu                 = " << ngpu          << endl;
        //   cout << " ccsdt_tilesize       = " << ccsdt_tilesize << endl;
        // }
        cout << " ndiis                = " << ndiis << endl;
        cout << " printtol             = " << std::scientific
             << std::setprecision(2) << printtol << endl;
        cout << " threshold            = " << std::scientific
             << std::setprecision(2) << threshold << endl;
        cout << " tilesize             = " << tilesize << endl;
        cout << " ccsd_maxiter         = " << ccsd_maxiter << endl;
        cout << " freeze_core          = " << freeze_core << endl;
        cout << " freeze_virtual       = " << freeze_virtual << endl;
        cout << " itilesize            = " << itilesize << endl;
        if(lshift != 0) cout << " lshift               = " << lshift << endl;
        print_bool(" readt               ", readt);
        print_bool(" writet              ", writet);
        print_bool(" writev              ", writev);
        // print_bool(" computeTData        ", computeTData);
        cout << " writet_iter          = " << writet_iter << endl;
        print_bool(" profile_ccsd        ", profile_ccsd);
        print_bool(" balance_tiles       ", balance_tiles);

        if(!dlpno_dfbasis.empty())
            cout << " dlpno_dfbasis        = " << dlpno_dfbasis << endl;
        if(!doubles_opt_eqns.empty()) {
            cout << " doubles_opt_eqns        = [";
            for(auto x : doubles_opt_eqns) cout << x << ",";
            cout << "]" << endl;
        }

        if(!ext_data_path.empty()) {
            cout << " ext_data_path   = " << ext_data_path << endl;
        }

        print_bool(" debug               ", debug);
        cout << "}" << endl;
    }
};

class OptionsMap {
public:
    OptionsMap() = default;
    Options options;
    CDOptions cd_options;
    CCSDOptions ccsd_options;
};
