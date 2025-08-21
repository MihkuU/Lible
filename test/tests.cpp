#include <tests.hpp>
#include <available_basis_sets.hpp>

#include <lible/ints/ints.hpp>

#include <vector>

namespace ltests = lible::tests;
namespace lints = lible::ints;

namespace lible::tests
{
    const static double tol = 1e-12;

    // H2O
    std::vector<int> atomic_nrs_h2o{8, 1, 1};
    std::vector<double> coords_h2o{0.0, 0.0, 0.0,
                                   0.757, 0.586, 0.0,
                                   -0.757, 0.586, 0.0};

    // Ethylene
    std::vector<int> atomic_nrs_c2h4{6, 6, 1, 1, 1, 1};
    std::vector<double> coords_c2h4{3.402, 0.773, -9.252,
                                    4.697, 0.791, -8.909,
                                    2.933, -0.150, -9.521,
                                    2.837, 1.682, -9.258,
                                    5.262, -0.118, -8.904,
                                    5.167, 1.714, -8.641};

    // CO2-
    std::vector<int> atomic_nrs_co2{6, 8, 8};
    std::vector<double> coords_co2{0.00000, 0.00000, 0.32942,
                                   0.00000, 1.15362, -0.12353,
                                   0.00000, -1.15362, -0.12353};
}

bool ltests::numCartesians()
{
    const int correct_answer = 21;

    int l = 5;
    int num_cartesians = lints::numCartesians(l);

    if (num_cartesians == correct_answer)
        return true;
    else
        return false;
}

bool ltests::numSphericals()
{
    const int correct_answer = 13;

    int l = 6;
    int num_sphericals = lints::numSphericals(l);

    if (num_sphericals == correct_answer)
        return true;
    else
        return false;
}

bool ltests::eri2Diagonal()
{
    const double correct_answer = 830.788976485982;

    lints::Structure structure("def2-SVP", "def2-universal-jkfit", atomic_nrs_h2o, coords_h2o);

    std::vector<double> eri2_diag = lints::eri2Diagonal(structure);

    double sum_eri2_diag = 0;
    for (size_t idx = 0; idx < eri2_diag.size(); idx++)
        sum_eri2_diag += std::fabs(eri2_diag[idx]);

    if (std::fabs(sum_eri2_diag - correct_answer) < tol)
        return true;
    else
        return false;
}

bool ltests::overlap()
{
    const double correct_answer = 226.099184393591;

    lints::Structure structure("cc-pvdz", atomic_nrs_c2h4, coords_c2h4);

    lible::vec2d sints = lints::overlap(structure);

    double sum_sints = 0;
    for (size_t i = 0; i < sints.size(); i++)
        sum_sints += std::fabs(sints[i]);

    if (std::fabs(sum_sints - correct_answer) < tol)
        return true;
    else
        return false;
}

bool ltests::overlapKernel()
{
    const double correct_answer = 137.049592196795;

    lints::Structure structure("cc-pvdz", atomic_nrs_c2h4, coords_c2h4);
    std::vector<lints::ShellPairData> shell_pair_datas = lints::shellPairDatas(true, structure);

    double sum_sints = 0;
    for (const auto &sp_data : shell_pair_datas)
        for (int ipair = 0; ipair < sp_data.n_pairs; ipair++)
        {
            lible::vec2d sints_batch = lints::overlapKernel(ipair, sp_data);

            for (size_t a = 0; a < sints_batch.dim<0>(); a++)
                for (size_t b = 0; b < sints_batch.dim<1>(); b++)
                    sum_sints += std::fabs(sints_batch(a, b));
        }

    if (std::fabs(sum_sints - correct_answer) < tol)
        return true;
    else
        return false;
}

bool ltests::overlapD1Kernel()
{
    const double correct_answer = 712.697494155785;

    lints::Structure structure("cc-pvdz", atomic_nrs_c2h4, coords_c2h4);
    std::vector<lints::ShellPairData> shell_pair_datas = lints::shellPairDatas(false, structure);

    double sum_sints = 0;
    for (const auto &sp_data : shell_pair_datas)
        for (int ipair = 0; ipair < sp_data.n_pairs; ipair++)
        {
            std::array<lible::vec2d, 6> sints_batch = lints::overlapD1Kernel(ipair, sp_data);

            for (int icart = 0; icart < 6; icart++)
                for (size_t a = 0; a < sints_batch[icart].dim<0>(); a++)
                    for (size_t b = 0; b < sints_batch[icart].dim<1>(); b++)
                        sum_sints += std::fabs(sints_batch[icart](a, b));
        }

    if (std::fabs(sum_sints - correct_answer) < tol)
        return true;
    else
        return false;
}

bool ltests::kineticEnergy()
{
    const double correct_answer = 252.755644635021;

    lints::Structure structure("cc-pvdz", atomic_nrs_c2h4, coords_c2h4);

    lible::vec2d kints = lints::kineticEnergy(structure);

    double sum_kints = 0;
    for (size_t i = 0; i < kints.size(); i++)
        sum_kints += std::fabs(kints[i]);

    if (std::fabs(sum_kints - correct_answer) < tol)
        return true;
    else
        return false;
}

bool ltests::kineticEnergyKernel()
{
    const double correct_answer = 252.755644635020;

    lints::Structure structure("cc-pvdz", atomic_nrs_c2h4, coords_c2h4);
    std::vector<lints::ShellPairData> shell_pair_datas = lints::shellPairDatas(false, structure);

    double sum_tints = 0;
    for (const auto &sp_data : shell_pair_datas)
        for (int ipair = 0; ipair < sp_data.n_pairs; ipair++)
        {
            lible::vec2d tints_batch = lints::kineticEnergyKernel(ipair, sp_data);

            for (size_t a = 0; a < tints_batch.dim<0>(); a++)
                for (size_t b = 0; b < tints_batch.dim<1>(); b++)
                    sum_tints += std::fabs(tints_batch(a, b));
        }

    if (std::fabs(sum_tints - correct_answer) < tol)
        return true;
    else
        return false;
}

bool ltests::kineticEnergyD1Kernel()
{
    const double correct_answer = 1153.592041307788;

    lints::Structure structure("cc-pvdz", atomic_nrs_c2h4, coords_c2h4);
    std::vector<lints::ShellPairData> shell_pair_datas = lints::shellPairDatas(false, structure);

    double sum_tints = 0;
    for (const auto &sp_data : shell_pair_datas)
        for (int ipair = 0; ipair < sp_data.n_pairs; ipair++)
        {
            std::array<lible::vec2d, 6> tints_batch = lints::kineticEnergyD1Kernel(ipair, sp_data);

            for (int icart = 0; icart < 6; icart++)
                for (size_t a = 0; a < tints_batch[icart].dim<0>(); a++)
                    for (size_t b = 0; b < tints_batch[icart].dim<1>(); b++)
                        sum_tints += std::fabs(tints_batch[icart](a, b));
        }

    if (std::fabs(sum_tints - correct_answer) < tol)
        return true;
    else
        return false;
}

bool ltests::nuclearAttraction()
{
    const double correct_answer = 2141.927982294058;

    lints::Structure structure("6-31+g", atomic_nrs_co2, coords_co2);

    lible::vec2d nints = lints::nuclearAttraction(structure);

    double sum_nints = 0;
    for (size_t i = 0; i < nints.size(); i++)
        sum_nints += std::fabs(nints[i]);

    if (std::fabs(sum_nints - correct_answer) < tol)
        return true;
    else
        return false;
}

bool ltests::nuclearAttractionErf()
{
    const double correct_answer = 257.523987741588;

    lints::Structure structure("6-31+g", atomic_nrs_co2, coords_co2);

    std::vector<double> omegas = {0.0035, -0.5, 0.6};
    lible::vec2d nints = lints::nuclearAttractionErf(structure, omegas);

    double sum_nints = 0;
    for (size_t i = 0; i < nints.size(); i++)
        sum_nints += std::fabs(nints[i]);

    if (std::fabs(sum_nints - correct_answer) < tol)
        return true;
    else
        return false;
}

bool ltests::externalCharges()
{
    const double correct_answer = 2141.927982294058;

    lints::Structure structure("6-31+g", atomic_nrs_co2, coords_co2);

    std::vector<std::array<double, 4>> charges = structure.getZs();

    lible::vec2d ecints = lints::externalCharges(charges, structure);

    double sum_ecints = 0;
    for (size_t i = 0; i < ecints.size(); i++)
        sum_ecints += std::fabs(ecints[i]);

    if (std::fabs(sum_ecints - correct_answer) < tol)
        return true;
    else
        return false;
}

bool ltests::externalChargesErf()
{
    const double correct_answer = 419.197077258843;

    lints::Structure structure("6-31+g", atomic_nrs_co2, coords_co2);

    std::vector<std::array<double, 4>> charges = structure.getZs();
    std::vector<double> omegas = {1.2, -10.5, 0.031};

    lible::vec2d ecints = lints::externalChargesErf(charges, omegas, structure);

    double sum_ecints = 0;
    for (size_t i = 0; i < ecints.size(); i++)
        sum_ecints += std::fabs(ecints[i]);

    if (std::fabs(sum_ecints - correct_answer) < tol)
        return true;
    else
        return false;
}

bool ltests::externalChargesKernel()
{
    const double correct_answer = 1338.606063918550;

    lints::Structure structure("6-31+g", atomic_nrs_co2, coords_co2);

    std::vector<lints::ShellPairData> shell_pair_datas = lints::shellPairDatas(true, structure);

    std::vector<std::array<double, 4>> charges = structure.getZs();

    double sum_ecints = 0;
    for (const lints::ShellPairData &sp_data : shell_pair_datas)
    {
        auto [la, lb] = sp_data.getLPair();
        int lab = la + lb;
        lints::BoysGrid boys_grid(lab);

        for (int ipair = 0; ipair < sp_data.n_pairs; ipair++)
        {
            lible::vec2d ecints = lints::externalChargesKernel(ipair, charges, boys_grid, sp_data);

            for (size_t a = 0; a < ecints.dim<0>(); a++)
                for (size_t b = 0; b < ecints.dim<1>(); b++)
                    sum_ecints += std::fabs(ecints(a, b));
        }
    }

    if (std::fabs(sum_ecints - correct_answer) < tol)
        return true;
    else
        return false;
}

bool ltests::externalChargesErfKernel()
{
    const double correct_answer = 419.197077258844;

    lints::Structure structure("6-31+g", atomic_nrs_co2, coords_co2);

    std::vector<lints::ShellPairData> shell_pair_datas = lints::shellPairDatas(false, structure);

    std::vector<std::array<double, 4>> charges = structure.getZs();
    std::vector<double> omegas = {1.2, -10.5, 0.031};

    double sum_ecints = 0;
    for (const lints::ShellPairData &sp_data : shell_pair_datas)
    {
        auto [la, lb] = sp_data.getLPair();
        int lab = la + lb;
        lints::BoysGrid boys_grid(lab);

        for (int ipair = 0; ipair < sp_data.n_pairs; ipair++)
        {
            lible::vec2d ecints = lints::externalChargesErfKernel(ipair, charges, omegas,
                                                                  boys_grid, sp_data);

            for (size_t a = 0; a < ecints.dim<0>(); a++)
                for (size_t b = 0; b < ecints.dim<1>(); b++)
                    sum_ecints += std::fabs(ecints(a, b));
        }
    }

    if (std::fabs(sum_ecints - correct_answer) < tol)
        return true;
    else
        return false;
}

bool ltests::externalChargesD1Kernel()
{
    const double correct_answer = 3274.419345348967;

    lints::Structure structure("6-31+g", atomic_nrs_co2, coords_co2);

    std::vector<lints::ShellPairData> shell_pair_datas = lints::shellPairDatas(true, structure);

    std::vector<std::array<double, 4>> charges = structure.getZs();

    double sum_ecints = 0;
    for (const lints::ShellPairData &sp_data : shell_pair_datas)
    {
        auto [la, lb] = sp_data.getLPair();
        int lab = la + lb;
        lints::BoysGrid boys_grid(lab);

        for (int ipair = 0; ipair < sp_data.n_pairs; ipair++)
        {
            std::array<lible::vec2d, 6> ecints = lints::externalChargesD1Kernel(ipair, charges,
                                                                                boys_grid, sp_data);

            for (int icart = 0; icart < 6; icart++)
                for (size_t a = 0; a < ecints[icart].dim<0>(); a++)
                    for (size_t b = 0; b < ecints[icart].dim<1>(); b++)
                        sum_ecints += std::fabs(ecints[icart](a, b));
        }
    }

    if (std::fabs(sum_ecints - correct_answer) < tol)
        return true;
    else
        return false;
}

bool ltests::externalChargesOperatorD1Kernel()
{
    const double correct_answer = 1067.752096834278;

    lints::Structure structure("6-31+g", atomic_nrs_co2, coords_co2);

    std::vector<lints::ShellPairData> shell_pair_datas = lints::shellPairDatas(true, structure);

    std::vector<std::array<double, 4>> charges = structure.getZs();

    double sum_ecints = 0;
    for (const lints::ShellPairData &sp_data : shell_pair_datas)
    {
        auto [la, lb] = sp_data.getLPair();
        int lab = la + lb;
        lints::BoysGrid boys_grid(lab);

        for (int ipair = 0; ipair < sp_data.n_pairs; ipair++)
        {
            std::vector<std::array<lible::vec2d, 3>> ecints =
                lints::externalChargesOperatorD1Kernel(ipair, charges, boys_grid, sp_data);

            for (const auto &ecints_at_charge : ecints)
                for (int icart = 0; icart < 3; icart++)
                    for (size_t a = 0; a < ecints_at_charge[icart].dim<0>(); a++)
                        for (size_t b = 0; b < ecints_at_charge[icart].dim<1>(); b++)
                            sum_ecints += std::fabs(ecints_at_charge[icart](a, b));
        }
    }

    if (std::fabs(sum_ecints - correct_answer) < tol)
        return true;
    else
        return false;
}

bool ltests::potentialAtExternalChargesKernel()
{
    const double correct_answer = 1357.142718060126;

    lints::Structure structure("6-31+g", atomic_nrs_co2, coords_co2);

    std::vector<lints::ShellPairData> shell_pair_datas = lints::shellPairDatas(true, structure);

    std::vector<std::array<double, 4>> charges = structure.getZs();

    double sum_ecints = 0;
    for (const lints::ShellPairData &sp_data : shell_pair_datas)
    {
        auto [la, lb] = sp_data.getLPair();
        int lab = la + lb;
        lints::BoysGrid boys_grid(lab);

        for (int ipair = 0; ipair < sp_data.n_pairs; ipair++)
        {
            std::vector<lible::vec2d> ecints =
                lints::potentialAtExternalChargesKernel(ipair, charges, boys_grid, sp_data);

            for (const auto &ecints_at_charge : ecints)            
                for (size_t a = 0; a < ecints_at_charge.dim<0>(); a++)
                    for (size_t b = 0; b < ecints_at_charge.dim<1>(); b++)
                        sum_ecints += std::fabs(ecints_at_charge(a, b));            
        }
    }

    if (std::fabs(sum_ecints - correct_answer) < tol)
        return true;
    else
        return false;
}

bool ltests::potentialAtExternalChargesErfKernel()
{
    const double correct_answer = 66.950210612294;

    lints::Structure structure("6-31+g", atomic_nrs_co2, coords_co2);

    std::vector<lints::ShellPairData> shell_pair_datas = lints::shellPairDatas(true, structure);

    std::vector<std::array<double, 4>> charges{{0.0, -0.025, 0.05, 1.0},
                                               {1.0, -0.025, 0.05, -1.0},
                                               {-1.0, 0.5, 0.6, 0.25},
                                               {0.0, -0.5, 0.24, -0.3}};
    std::vector<double> omegas{0.1, -0.5, 0.61, 0.000015};

    double sum_ecints = 0;
    for (const lints::ShellPairData &sp_data : shell_pair_datas)
    {
        auto [la, lb] = sp_data.getLPair();
        int lab = la + lb;
        lints::BoysGrid boys_grid(lab);

        for (int ipair = 0; ipair < sp_data.n_pairs; ipair++)
        {
            std::vector<lible::vec2d> ecints =
                lints::potentialAtExternalChargesErfKernel(ipair, charges, omegas, boys_grid,
                                                           sp_data);

            for (const auto &ecints_at_charge : ecints)
                for (size_t a = 0; a < ecints_at_charge.dim<0>(); a++)
                    for (size_t b = 0; b < ecints_at_charge.dim<1>(); b++)
                        sum_ecints += std::fabs(ecints_at_charge(a, b));
        }
    }

    if (std::fabs(sum_ecints - correct_answer) < tol)
        return true;
    else
        return false;
}

bool ltests::availableBasisSets()
{
    std::set<std::string> available_basis_sets = lints::availableBasisSets();
    
    for (const std::string &basis_set : available_basis_sets)
        if (basis_sets_table.contains(basis_set) == false)
            return false;
        
    return true;
}

bool ltests::availableBasisSetsAux()
{
    std::set<std::string> available_basis_sets_aux = lints::availableBasisSetsAux();

    for (const std::string &basis_set : available_basis_sets_aux)
        if (auxbasis_sets_table.contains(basis_set) == false)
            return false;
        
    return true;
}

bool ltests::purePrimitiveNorm()
{
    const double correct_answer = 8.063341937228;    

    double exp = 1.70325241;
    int l = 4;
    double norm = lints::purePrimitiveNorm(exp, l);
    
    if (std::fabs(norm - correct_answer) < tol)
        return true;
    else
        return false;    
}

bool ltests::numHermites()
{
    const int correct_answer = 35;

    int l = 4;
    int num_hermites = lints::numHermites(l);

    return (num_hermites == correct_answer);
}

bool ltests::cartExps()
{
    const std::vector<std::array<int, 3>> correct_answer{
        {5, 0, 0},
        {4, 1, 0},
        {4, 0, 1},
        {3, 2, 0},
        {3, 1, 1},
        {3, 0, 2},
        {2, 3, 0},
        {2, 2, 1},
        {2, 1, 2},
        {2, 0, 3},
        {1, 4, 0},
        {1, 3, 1},
        {1, 2, 2},
        {1, 1, 3},
        {1, 0, 4},
        {0, 5, 0},
        {0, 4, 1},
        {0, 3, 2},
        {0, 2, 3},
        {0, 1, 4},
        {0, 0, 5}};

    int l = 5;
    std::vector<std::array<int, 3>> cart_exps = lints::cartExps(l);
    
    return (cart_exps == correct_answer);
}

bool ltests::getLPairsSymm()
{
    const std::vector<std::pair<int, int>> correct_answer{
        {0, 0},
        {1, 0},
        {1, 1},
        {2, 0},
        {2, 1},
        {2, 2},
        {3, 0},
        {3, 1},
        {3, 2},
        {3, 3},
        {4, 0},
        {4, 1},
        {4, 2},
        {4, 3},
        {4, 4},
        {5, 0},
        {5, 1},
        {5, 2},
        {5, 3},
        {5, 4},
        {5, 5}};

    int l_max = 5;
    std::vector<std::pair<int, int>> l_pairs = lints::getLPairsSymm(l_max);

    return (l_pairs == correct_answer);
}

bool ltests::getLPairsNoSymm()
{
    const std::vector<std::pair<int, int>> correct_answer{
        {0, 0},
        {0, 1},
        {0, 2},
        {0, 3},
        {0, 4},
        {1, 0},
        {1, 1},
        {1, 2},
        {1, 3},
        {1, 4},
        {2, 0},
        {2, 1},
        {2, 2},
        {2, 3},
        {2, 4},
        {3, 0},
        {3, 1},
        {3, 2},
        {3, 3},
        {3, 4},
        {4, 0},
        {4, 1},
        {4, 2},
        {4, 3},
        {4, 4}};

    int l_max = 4;
    std::vector<std::pair<int, int>> l_pairs = lints::getLPairsNoSymm(l_max);

    return (l_pairs == correct_answer);
}