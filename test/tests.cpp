#include <tests.hpp>
#include <available_basis_sets.hpp>

#include <lible/ints/ints.hpp>

#include <vector>

namespace ltests = lible::tests;
namespace lints = lible::ints;

namespace lible::tests
{
    const static double tol = 1e-12;

    // CO2-
    std::vector<int> atomic_nrs_co2{6, 8, 8};
    std::vector<std::array<double, 3>> coords_co2{
        {0.00000, 0.00000, 0.32942},
        {0.00000, 1.15362, -0.12353},
        {0.00000, -1.15362, -0.12353}
    };

    // Ethane
    std::vector<int> atomic_nrs_c2h6{1, 6, 1, 1, 6, 1, 1, 1};
    std::vector<std::array<double, 3>> coords_c2h6{
        {1.1851, -0.0039, 0.9875},
        {0.7516, -0.0225, -0.0209},
        {1.1669, 0.8330, -0.5693},
        {1.1155, -0.9329, -0.5145},
        {-0.7516, 0.0225, 0.0209},
        {-1.1669, -0.8334, 0.5687},
        {-1.1157, 0.9326, 0.5151},
        {-1.1850, 0.0044, -0.9875}
    };

    // FeH6
    std::vector<int> atomic_nrs_feh6{26, 1, 1, 1, 1, 1, 1};
    std::vector<std::array<double, 3>> coords_feh6{
        {-0.00000000700783, -0.00000000961045, -0.00000001699045},
        {-1.90000000000000, -0.06531984963340, 0.06531987277077},
        {1.90000000000000, 0.06531978816537, -0.06531984034895},
        {-0.06531980735085, -1.69016064020970, 0.06531920821976},
        {0.06531974281206, 1.69016062544695, -0.06531916850398},
        {0.06531985909338, 0.06531923539693, -1.69016066429686},
        {-0.06531977461945, -0.06531914955570, 1.69016060914971}
    };

    // Ethylene
    std::vector<int> atomic_nrs_c2h4{6, 6, 1, 1, 1, 1};
    std::vector<std::array<double, 3>> coords_c2h4{
        {3.402, 0.773, -9.252},
        {4.697, 0.791, -8.909},
        {2.933, -0.150, -9.521},
        {2.837, 1.682, -9.258},
        {5.262, -0.118, -8.904},
        {5.167, 1.714, -8.641}
    };

    // H2O
    std::vector<int> atomic_nrs_h2o{8, 1, 1};
    std::vector<std::array<double, 3>> coords_h2o{
        {0.0, 0.0, 0.0},
        {0.757, 0.586, 0.0},
        {-0.757, 0.586, 0.0}
    };

    // Ozone
    std::vector<int> atomic_nrs_o3{8, 8, 8};
    std::vector<std::array<double, 3>> coords_o3{
        {0, 0, 0},
        {0, 1, 1},
        {0, 0, 2}
    };
}

bool ltests::numCartesians()
{
    const int correct_answer = 21;

    int l = 5;
    int num_cartesians = lints::numCartesians(l);

    if (num_cartesians == correct_answer)
        return true;

    return false;
}

bool ltests::numSphericals()
{
    const int correct_answer = 13;

    int l = 6;
    int num_sphericals = lints::numSphericals(l);

    if (num_sphericals == correct_answer)
        return true;

    return false;
}

bool ltests::eri2Diagonal()
{
    const double correct_answer = 830.788976485982;

    lints::Structure structure("def2-SVP", "def2-universal-jkfit", atomic_nrs_h2o, coords_h2o);

    std::vector<double> eri2_diag = lints::eri2Diagonal(structure);

    double sum_eri2_diag = 0;
    for (double idx : eri2_diag)
        sum_eri2_diag += std::fabs(idx);

    if (std::fabs(sum_eri2_diag - correct_answer) < tol)
        return true;

    return false;
}

bool ltests::overlap()
{
    const double correct_answer = 226.099184393591;

    lints::Structure structure("cc-pvdz", atomic_nrs_c2h4, coords_c2h4);

    vec2d sints = lints::overlap(structure);

    double sum_sints = 0;
    for (double sint : sints)
        sum_sints += std::fabs(sint);

    if (std::fabs(sum_sints - correct_answer) < tol)
        return true;

    return false;
}

bool ltests::overlapKernel()
{
    constexpr double correct_answer = 137.049592196795;

    lints::Structure structure("cc-pvdz", atomic_nrs_c2h4, coords_c2h4);
    std::vector<lints::ShellPairData> shell_pair_datas = lints::shellPairData(true, structure);

    double sum_sints = 0;
    for (const auto &sp_data : shell_pair_datas)
    {
        for (size_t ipair = 0; ipair < sp_data.n_pairs_; ipair++)
        {
            vec2d sints_batch = lints::overlapKernel(ipair, sp_data);

            for (size_t a = 0; a < sints_batch.dim<0>(); a++)
                for (size_t b = 0; b < sints_batch.dim<1>(); b++)
                    sum_sints += std::fabs(sints_batch(a, b));
        }
    }

    if (std::fabs(sum_sints - correct_answer) < tol)
        return true;

    return false;
}

bool ltests::overlapD1Kernel()
{
    const double correct_answer = 712.697494155785;

    lints::Structure structure("cc-pvdz", atomic_nrs_c2h4, coords_c2h4);
    std::vector<lints::ShellPairData> shell_pair_datas = lints::shellPairData(false, structure);

    double sum_sints = 0;
    for (const auto &sp_data : shell_pair_datas)
        for (size_t ipair = 0; ipair < sp_data.n_pairs_; ipair++)
        {
            std::array<vec2d, 6> sints_batch = lints::overlapD1Kernel(ipair, sp_data);

            for (int icart = 0; icart < 6; icart++)
                for (size_t a = 0; a < sints_batch[icart].dim<0>(); a++)
                    for (size_t b = 0; b < sints_batch[icart].dim<1>(); b++)
                        sum_sints += std::fabs(sints_batch[icart](a, b));
        }

    if (std::fabs(sum_sints - correct_answer) < tol)
        return true;

    return false;
}

bool ltests::kineticEnergy()
{
    const double correct_answer = 252.755644635021;

    lints::Structure structure("cc-pvdz", atomic_nrs_c2h4, coords_c2h4);

    vec2d kints = lints::kineticEnergy(structure);

    double sum_kints = 0;
    for (double kint : kints)
        sum_kints += std::fabs(kint);

    if (std::fabs(sum_kints - correct_answer) < tol)
        return true;

    return false;
}

bool ltests::kineticEnergyKernel()
{
    const double correct_answer = 252.755644635020;

    lints::Structure structure("cc-pvdz", atomic_nrs_c2h4, coords_c2h4);
    std::vector<lints::ShellPairData> shell_pair_datas = lints::shellPairData(false, structure);

    double sum_tints = 0;
    for (const auto &sp_data : shell_pair_datas)
        for (size_t ipair = 0; ipair < sp_data.n_pairs_; ipair++)
        {
            vec2d tints_batch = lints::kineticEnergyKernel(ipair, sp_data);

            for (size_t a = 0; a < tints_batch.dim<0>(); a++)
                for (size_t b = 0; b < tints_batch.dim<1>(); b++)
                    sum_tints += std::fabs(tints_batch(a, b));
        }

    if (std::fabs(sum_tints - correct_answer) < tol)
        return true;

    return false;
}

bool ltests::kineticEnergyD1Kernel()
{
    const double correct_answer = 1153.592041307799;

    lints::Structure structure("cc-pvdz", atomic_nrs_c2h4, coords_c2h4);
    std::vector<lints::ShellPairData> shell_pair_datas = lints::shellPairData(false, structure);

    double sum_tints = 0;
    for (const auto &sp_data : shell_pair_datas)
        for (size_t ipair = 0; ipair < sp_data.n_pairs_; ipair++)
        {
            std::array<vec2d, 6> tints_batch = lints::kineticEnergyD1Kernel(ipair, sp_data);

            for (int icart = 0; icart < 6; icart++)
                for (size_t a = 0; a < tints_batch[icart].dim<0>(); a++)
                    for (size_t b = 0; b < tints_batch[icart].dim<1>(); b++)
                        sum_tints += std::fabs(tints_batch[icart](a, b));
        }

    if (std::fabs(sum_tints - correct_answer) < tol)
        return true;

    return false;
}

bool ltests::nuclearAttraction()
{
    const double correct_answer = 2141.927982294058;

    lints::Structure structure("6-31+g", atomic_nrs_co2, coords_co2);

    vec2d nints = lints::nuclearAttraction(structure);

    double sum_nints = 0;
    for (double nint : nints)
        sum_nints += std::fabs(nint);

    if (std::fabs(sum_nints - correct_answer) < tol)
        return true;

    return false;
}

bool ltests::nuclearAttractionErf()
{
    const double correct_answer = 257.523987741588;

    lints::Structure structure("6-31+g", atomic_nrs_co2, coords_co2);

    std::vector<double> omegas = {0.0035, -0.5, 0.6};
    vec2d nints = lints::nuclearAttractionErf(structure, omegas);

    double sum_nints = 0;
    for (double nint : nints)
        sum_nints += std::fabs(nint);

    if (std::fabs(sum_nints - correct_answer) < tol)
        return true;

    return false;
}

bool ltests::externalCharges()
{
    const double correct_answer = 2141.927982294058;

    lints::Structure structure("6-31+g", atomic_nrs_co2, coords_co2);

    std::vector<std::array<double, 4>> charges = structure.getZs();

    vec2d ecints = lints::externalCharges(charges, structure);

    double sum_ecints = 0;
    for (double ecint : ecints)
        sum_ecints += std::fabs(ecint);

    if (std::fabs(sum_ecints - correct_answer) < tol)
        return true;

    return false;
}

bool ltests::externalChargesErf()
{
    const double correct_answer = 419.197077258843;

    lints::Structure structure("6-31+g", atomic_nrs_co2, coords_co2);

    std::vector<std::array<double, 4>> charges = structure.getZs();
    std::vector<double> omegas = {1.2, -10.5, 0.031};

    vec2d ecints = lints::externalChargesErf(charges, omegas, structure);

    double sum_ecints = 0;
    for (double ecint : ecints)
        sum_ecints += std::fabs(ecint);

    if (std::fabs(sum_ecints - correct_answer) < tol)
        return true;

    return false;
}

bool ltests::externalChargesKernel()
{
    const double correct_answer = 1338.606063918550;

    lints::Structure structure("6-31+g", atomic_nrs_co2, coords_co2);

    std::vector<lints::ShellPairData> shell_pair_datas = lints::shellPairData(true, structure);

    std::vector<std::array<double, 4>> charges = structure.getZs();

    double sum_ecints = 0;
    for (const lints::ShellPairData &sp_data : shell_pair_datas)
    {
        auto [la, lb] = sp_data.getLPair();
        int lab = la + lb;
        lints::BoysGrid boys_grid(lab);

        for (size_t ipair = 0; ipair < sp_data.n_pairs_; ipair++)
        {
            vec2d ecints = lints::externalChargesKernel(ipair, charges, boys_grid, sp_data);

            for (size_t a = 0; a < ecints.dim<0>(); a++)
                for (size_t b = 0; b < ecints.dim<1>(); b++)
                    sum_ecints += std::fabs(ecints(a, b));
        }
    }

    if (std::fabs(sum_ecints - correct_answer) < tol)
        return true;

    return false;
}

bool ltests::externalChargesErfKernel()
{
    const double correct_answer = 419.197077258844;

    lints::Structure structure("6-31+g", atomic_nrs_co2, coords_co2);

    std::vector<lints::ShellPairData> shell_pair_datas = lints::shellPairData(false, structure);

    std::vector<std::array<double, 4>> charges = structure.getZs();
    std::vector<double> omegas = {1.2, -10.5, 0.031};

    double sum_ecints = 0;
    for (const lints::ShellPairData &sp_data : shell_pair_datas)
    {
        auto [la, lb] = sp_data.getLPair();
        int lab = la + lb;
        lints::BoysGrid boys_grid(lab);

        for (size_t ipair = 0; ipair < sp_data.n_pairs_; ipair++)
        {
            vec2d ecints = lints::externalChargesErfKernel(ipair, charges, omegas,
                                                           boys_grid, sp_data);

            for (size_t a = 0; a < ecints.dim<0>(); a++)
                for (size_t b = 0; b < ecints.dim<1>(); b++)
                    sum_ecints += std::fabs(ecints(a, b));
        }
    }

    if (std::fabs(sum_ecints - correct_answer) < tol)
        return true;

    return false;
}

bool ltests::externalChargesD1Kernel()
{
    const double correct_answer = 3274.419345348968;

    lints::Structure structure("6-31+g", atomic_nrs_co2, coords_co2);

    std::vector<lints::ShellPairData> shell_pair_datas = lints::shellPairData(true, structure);

    std::vector<std::array<double, 4>> charges = structure.getZs();

    double sum_ecints = 0;
    for (const lints::ShellPairData &sp_data : shell_pair_datas)
    {
        auto [la, lb] = sp_data.getLPair();
        int lab = la + lb;
        lints::BoysGrid boys_grid(lab);

        for (size_t ipair = 0; ipair < sp_data.n_pairs_; ipair++)
        {
            std::array<vec2d, 6> ecints = lints::externalChargesD1Kernel(ipair, charges,
                                                                         boys_grid, sp_data);

            for (int icart = 0; icart < 6; icart++)
                for (size_t a = 0; a < ecints[icart].dim<0>(); a++)
                    for (size_t b = 0; b < ecints[icart].dim<1>(); b++)
                        sum_ecints += std::fabs(ecints[icart](a, b));
        }
    }

    if (std::fabs(sum_ecints - correct_answer) < tol)
        return true;

    return false;
}

bool ltests::externalChargesOperatorD1Kernel()
{
    const double correct_answer = 1067.752096834278;

    lints::Structure structure("6-31+g", atomic_nrs_co2, coords_co2);

    std::vector<lints::ShellPairData> shell_pair_datas = lints::shellPairData(true, structure);

    std::vector<std::array<double, 4>> charges = structure.getZs();

    double sum_ecints = 0;
    for (const lints::ShellPairData &sp_data : shell_pair_datas)
    {
        auto [la, lb] = sp_data.getLPair();
        int lab = la + lb;
        lints::BoysGrid boys_grid(lab);

        for (size_t ipair = 0; ipair < sp_data.n_pairs_; ipair++)
        {
            std::vector<std::array<vec2d, 3>> ecints =
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

    return false;
}

bool ltests::potentialAtExternalChargesKernel()
{
    const double correct_answer = 1357.142718060126;

    lints::Structure structure("6-31+g", atomic_nrs_co2, coords_co2);

    std::vector<lints::ShellPairData> shell_pair_datas = lints::shellPairData(true, structure);

    std::vector<std::array<double, 4>> charges = structure.getZs();

    double sum_ecints = 0;
    for (const lints::ShellPairData &sp_data : shell_pair_datas)
    {
        auto [la, lb] = sp_data.getLPair();
        int lab = la + lb;
        lints::BoysGrid boys_grid(lab);

        for (size_t ipair = 0; ipair < sp_data.n_pairs_; ipair++)
        {
            std::vector<vec2d> ecints =
                    lints::potentialAtExternalChargesKernel(ipair, charges, boys_grid, sp_data);

            for (const auto &ecints_at_charge : ecints)
                for (size_t a = 0; a < ecints_at_charge.dim<0>(); a++)
                    for (size_t b = 0; b < ecints_at_charge.dim<1>(); b++)
                        sum_ecints += std::fabs(ecints_at_charge(a, b));
        }
    }

    if (std::fabs(sum_ecints - correct_answer) < tol)
        return true;

    return false;
}

bool ltests::potentialAtExternalChargesErfKernel()
{
    const double correct_answer = 66.950210612294;

    lints::Structure structure("6-31+g", atomic_nrs_co2, coords_co2);

    std::vector<lints::ShellPairData> shell_pair_datas = lints::shellPairData(true, structure);

    std::vector<std::array<double, 4>> charges{
        {0.0, -0.025, 0.05, 1.0},
        {1.0, -0.025, 0.05, -1.0},
        {-1.0, 0.5, 0.6, 0.25},
        {0.0, -0.5, 0.24, -0.3}
    };
    std::vector<double> omegas{0.1, -0.5, 0.61, 0.000015};

    double sum_ecints = 0;
    for (const lints::ShellPairData &sp_data : shell_pair_datas)
    {
        auto [la, lb] = sp_data.getLPair();
        int lab = la + lb;
        lints::BoysGrid boys_grid(lab);

        for (int ipair = 0; ipair < sp_data.n_pairs_; ipair++)
        {
            std::vector<vec2d> ecints =
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

bool ltests::dipoleMoment()
{
    const double correct_answer = 1205.074870756679;

    lints::Structure structure("cc-pvdz", atomic_nrs_c2h6, coords_c2h6);

    std::array<lible::vec2d, 3> dipole_moment_ints = lints::dipoleMoment({0, 0, 0}, structure);

    double sum_dipole_moment = 0;
    for (int icart = 0; icart < 3; icart++)
        for (size_t i = 0; i < dipole_moment_ints[icart].size(); i++)
            sum_dipole_moment += std::fabs(dipole_moment_ints[icart][i]);

    if (std::fabs(sum_dipole_moment - correct_answer) < tol)
        return true;
    else
        return false;
}

bool ltests::dipoleMomentKernel()
{
    const double correct_answer = 1205.074870756674;

    lints::Structure structure("cc-pvdz", atomic_nrs_c2h6, coords_c2h6);

    std::vector<lints::ShellPairData> shell_pair_datas = lints::shellPairData(false, structure);

    double sum_dipole_moment = 0;
    for (const lints::ShellPairData &sp_data : shell_pair_datas)
        for (size_t ipair = 0; ipair < sp_data.n_pairs_; ipair++)
        {
            std::array<vec2d, 3> dipole_moment_ints =
                    lints::dipoleMomentKernel(ipair, {0, 0, 0}, sp_data);

            for (int icart = 0; icart < 3; icart++)
                for (double val : dipole_moment_ints[icart])
                    sum_dipole_moment += std::fabs(val);
        }

    if (std::fabs(sum_dipole_moment - correct_answer) < tol)
        return true;

    return false;
}

bool ltests::spinOrbitCoupling1El()
{
    const double correct_answer = 86318.844650089319;

    lints::Structure structure("ahlrichs-vdz", atomic_nrs_feh6, coords_feh6);

    std::array<vec2d, 3> soc_ints = lints::spinOrbitCoupling1El(structure);

    double sum_soc_ints = 0;
    for (int icart = 0; icart < 3; icart++)
        for (double val : soc_ints[icart])
            sum_soc_ints += std::fabs(val);

    if (std::fabs(sum_soc_ints - correct_answer) < tol)
        return true;

    return false;
}

bool ltests::spinOrbitCoupling1ElKernel()
{
    const double correct_answer = 86318.844650089552;

    lints::Structure structure("ahlrichs-vdz", atomic_nrs_feh6, coords_feh6);

    std::vector<lints::ShellPairData> shell_pair_datas = lints::shellPairData(false, structure);

    std::vector<std::array<double, 4>> charges = structure.getZs();

    double sum_soc_ints = 0;
    for (const lints::ShellPairData &sp_data : shell_pair_datas)
    {
        auto [la, lb] = sp_data.getLPair();
        lints::BoysGrid boys_grid(la + lb + 1);

        for (size_t ipair = 0; ipair < sp_data.n_pairs_; ipair++)
        {
            std::array<vec2d, 3> soc_ints =
                    lints::spinOrbitCoupling1ElKernel(ipair, charges, boys_grid, sp_data);

            for (int icart = 0; icart < 3; icart++)
                for (double val : soc_ints[icart])
                    sum_soc_ints += std::fabs(val);
        }
    }

    if (std::fabs(sum_soc_ints - correct_answer) < tol)
        return true;

    return false;
}

bool ltests::eri2()
{
    const double correct_answer = 8801.334703460838;

    lints::Structure structure("def2-SVP", "def2-qzvp-rifit", atomic_nrs_h2o, coords_h2o);

    vec2d eri2 = lints::eri2(structure);

    double eri2_sum = 0;
    for (double val : eri2)
        eri2_sum += std::fabs(val);

    if (std::fabs(eri2_sum - correct_answer) < tol)
        return true;

    return false;
}

bool ltests::eri4Diagonal()
{
    const double correct_answer = 306.806419250122;

    lints::Structure structure("def2-tzvp", atomic_nrs_o3, coords_o3);

    vec2d eri4_diag = lints::eri4Diagonal(structure);

    double eri4_diag_sum = 0;
    for (double val : eri4_diag)
        eri4_diag_sum += std::fabs(val);

    if (std::fabs(eri4_diag_sum - correct_answer) < tol)
        return true;

    return false;
}

bool ltests::eri3()
{
    const double correct_answer = 25854.105851585416531;

    lints::Structure structure("def2-tzvp", "def2-qzvp-rifit", atomic_nrs_o3, coords_o3);

    vec3d eri3 = lints::eri3(structure);

    double eri3_sum = 0;
    for (double val : eri3)
        eri3_sum += std::fabs(val);

    if (std::fabs(eri3_sum - correct_answer) < tol)
        return true;

    return false;
}

bool ltests::eri4()
{
    const double correct_answer = 8246.711763196334;

    lints::Structure structure("def2-svp", atomic_nrs_o3, coords_o3);

    vec4d eri4 = lints::eri4(structure);

    double eri4_sum = 0;
    for (double val : eri4)
        eri4_sum += std::fabs(val);

    if (std::fabs(eri4_sum - correct_answer) < tol)
        return true;

    return false;
}

bool ltests::basisForAtom()
{
    const double correct_answer = 107960921.709905013442;

    lints::BasisAtom basis_atom = lints::basisForAtom(34, "ano-rcc-vtz");

    double sum_data = 0;
    for (const auto &[l, exps, coeffs] : basis_atom.basis_shells_)
    {
        for (double exp : exps)
            sum_data += std::fabs(exp);

        for (double coeff : coeffs)
            sum_data += std::fabs(coeff);
    }

    if (std::fabs(sum_data - correct_answer) < tol)
        return true;

    return false;
}

bool ltests::basisForAtomAux()
{
    const double correct_answer = 89.733725300000;

    lints::BasisAtom basis_atom = lints::basisForAtomAux(43, "def2-svpd-rifit");

    double sum_data = 0;
    for (const lints::BasisShell &basis_shell : basis_atom.basis_shells_)
    {
        for (double exp : basis_shell.exps_)
            sum_data += std::fabs(exp);

        for (double coeff : basis_shell.coeffs_)
            sum_data += std::fabs(coeff);
    }

    if (std::fabs(sum_data - correct_answer) < tol)
        return true;

    return false;
}

bool ltests::basisForAtoms()
{
    const double correct_answer = 55347.164895499984;

    lints::basis_atoms_t basis_atoms = lints::basisForAtoms({14, 6, 17}, "6-31g(3df,3pd)");

    double sum_data = 0;
    for (const lints::BasisAtom &basis_atom : basis_atoms)
        for (const auto& [l, exps, coeffs] : basis_atom.basis_shells_)
        {
            for (double exp : exps)
                sum_data += std::fabs(exp);

            for (double coeff : coeffs)
                sum_data += std::fabs(coeff);
        }

    if (std::fabs(sum_data - correct_answer) < tol)
        return true;

    return false;
}

bool ltests::basisForAtomsAux()
{
    const double correct_answer = 338.915714100950;

    lints::basis_atoms_t basis_atoms = lints::basisForAtomsAux({76, 19, 1}, "def2-svpd-rifit");

    double sum_data = 0;
    for (const lints::BasisAtom &basis_atom : basis_atoms)
        for (const auto &[l, exps, coeffs] : basis_atom.basis_shells_)
            {
                for (double exp : exps)
                    sum_data += std::fabs(exp);

                for (double coeff : coeffs)
                    sum_data += std::fabs(coeff);
            }

    if (std::fabs(sum_data - correct_answer) < tol)
        return true;

    return false;
}

bool ltests::sphericalTrafo()
{
    const double correct_answer = 7005.923248093791;

    int l = 8;
    std::vector<std::tuple<int, int, double>> trafo = lints::sphericalTrafo(l);

    double sum_trafo = 0;
    for (const auto &[i, j, val] : trafo)
        sum_trafo += std::fabs(val) + i + j;

    if (std::fabs(sum_trafo - correct_answer) < tol)
        return true;

    return false;
}

bool ltests::deployERI4Kernel()
{
    const double correct_answer = 8246.711763197583;

    lints::Structure structure("def2-svp", atomic_nrs_o3, coords_o3);

    std::vector<lints::ShellPairData> shell_pair_datas = lints::shellPairData(false, structure);

    double eri4_sum = 0;
    for (const lints::ShellPairData &sp_data_ab : shell_pair_datas)
        for (const lints::ShellPairData &sp_data_cd : shell_pair_datas)
        {
            lints::ERI4Kernel eri4_kernel(sp_data_ab, sp_data_cd);
            for (size_t ipair1 = 0; ipair1 < sp_data_ab.n_pairs_; ipair1++)
                for (size_t ipair2 = 0; ipair2 < sp_data_cd.n_pairs_; ipair2++)
                {
                    vec4d eri4_batch = eri4_kernel(ipair1, ipair2, sp_data_ab, sp_data_cd);

                    for (double val : eri4_batch)
                        eri4_sum += std::fabs(val);
                }
        }

    if (std::fabs(eri4_sum - correct_answer) < tol)
        return true;

    return false;
}

bool ltests::deployERI3Kernel()
{
    const double correct_answer = 12350.130354189232;

    lints::Structure structure("def2-svp", "def2-universal-jkfit", atomic_nrs_o3, coords_o3);

    std::vector<lints::ShellPairData> shell_pair_datas = lints::shellPairData(false, structure);
    std::vector<lints::ShellData> shell_datas = lints::shellDataAux(structure);

    double eri3_sum = 0;
    for (const lints::ShellPairData &sp_data_ab : shell_pair_datas)
        for (const lints::ShellData &sh_data_c : shell_datas)
        {
            lints::ERI3Kernel eri3_kernel(sp_data_ab, sh_data_c);
            for (size_t ipair = 0; ipair < sp_data_ab.n_pairs_; ipair++)
                for (size_t ishell = 0; ishell < sh_data_c.n_shells_; ishell++)
                {
                    vec3d eri3_batch = eri3_kernel(ipair, ishell, sp_data_ab, sh_data_c);

                    for (double val : eri3_batch)
                        eri3_sum += std::fabs(val);
                }
        }

    if (std::fabs(eri3_sum - correct_answer) < tol)
        return true;

    return false;
}

bool ltests::deployERI2Kernel()
{
    const double correct_answer = 21014.512841176023;

    lints::Structure structure("def2-svp", "def2-universal-jkfit", atomic_nrs_o3, coords_o3);

    std::vector<lints::ShellData> shell_datas = lints::shellDataAux(structure);

    double eri2_sum = 0;
    for (const lints::ShellData &sh_data_a : shell_datas)
        for (const lints::ShellData &sh_data_b : shell_datas)
        {
            lints::ERI2Kernel eri2_kernel(sh_data_a, sh_data_b);
            for (int ishell_a = 0; ishell_a < sh_data_a.n_shells_; ishell_a++)
                for (int ishell_b = 0; ishell_b < sh_data_b.n_shells_; ishell_b++)
                {
                    vec2d eri2_batch = eri2_kernel(ishell_a, ishell_b, sh_data_a, sh_data_b);

                    for (double val : eri2_batch)
                        eri2_sum += std::fabs(val);
                }
        }

    if (std::fabs(eri2_sum - correct_answer) < tol)
        return true;

    return false;
}

bool ltests::deployERI4D1Kernel()
{
    const double correct_answer = 30546.780073140635;

    lints::Structure structure("def2-svp", atomic_nrs_o3, coords_o3);

    std::vector<lints::ShellPairData> shell_pair_datas = lints::shellPairData(true, structure);

    double eri4_sum = 0;
    for (const lints::ShellPairData &sp_data_ab : shell_pair_datas)
        for (const lints::ShellPairData &sp_data_cd : shell_pair_datas)
        {
            lints::ERI4D1Kernel eri4d1_kernel(sp_data_ab, sp_data_cd);
            for (size_t ipair1 = 0; ipair1 < sp_data_ab.n_pairs_; ipair1++)
                for (size_t ipair2 = 0; ipair2 < sp_data_cd.n_pairs_; ipair2++)
                {
                    std::array<vec4d, 12> eri4_batch =
                            eri4d1_kernel(ipair1, ipair2, sp_data_ab, sp_data_cd);

                    for (int icart = 0; icart < 12; icart++)
                        for (double val : eri4_batch[icart])
                            eri4_sum += std::fabs(val);
                }
        }

    if (std::fabs(eri4_sum - correct_answer) < tol)
        return true;

    return false;
}

bool ltests::deployERI3D1Kernel()
{
    const double correct_answer = 78466.360848033844;

    lints::Structure structure("def2-svp", "def2-universal-jkfit", atomic_nrs_o3, coords_o3);

    std::vector<lints::ShellPairData> shell_pair_datas = lints::shellPairData(false, structure);
    std::vector<lints::ShellData> shell_datas = lints::shellDataAux(structure);

    double eri3_sum = 0;
    for (const lints::ShellPairData &sp_data_ab : shell_pair_datas)
        for (const lints::ShellData &sh_data_c : shell_datas)
        {
            lints::ERI3D1Kernel eri3d1_kernel(sp_data_ab, sh_data_c);
            for (size_t ipair = 0; ipair < sp_data_ab.n_pairs_; ipair++)
                for (size_t ishell = 0; ishell < sh_data_c.n_shells_; ishell++)
                {
                    std::array<vec3d, 9> eri3_batch =
                            eri3d1_kernel(ipair, ishell, sp_data_ab, sh_data_c);

                    for (int icart = 0; icart < 9; icart++)
                        for (double val : eri3_batch[icart])
                            eri3_sum += std::fabs(val);
                }
        }

    if (std::fabs(eri3_sum - correct_answer) < tol)
        return true;

    return false;
}

bool ltests::deployERI2D1Kernel()
{
    const double correct_answer = 39949.218421684251;

    lints::Structure structure("def2-svp", "def2-universal-jkfit", atomic_nrs_o3, coords_o3);

    std::vector<lints::ShellData> shell_datas = lints::shellDataAux(structure);

    double eri2_sum = 0;
    for (const lints::ShellData &sh_data_a : shell_datas)
        for (const lints::ShellData &sh_data_b : shell_datas)
        {
            lints::ERI2D1Kernel eri2d1_kernel(sh_data_a, sh_data_b);
            for (size_t ishell_a = 0; ishell_a < sh_data_a.n_shells_; ishell_a++)
                for (size_t ishell_b = 0; ishell_b < sh_data_b.n_shells_; ishell_b++)
                {
                    std::array<vec2d, 6> eri2_batch =
                            eri2d1_kernel(ishell_a, ishell_b, sh_data_a, sh_data_b);

                    for (int icart = 0; icart < 6; icart++)
                        for (double val : eri2_batch[icart])
                            eri2_sum += std::fabs(val);
                }
        }

    if (std::fabs(eri2_sum - correct_answer) < tol)
        return true;

    return false;
}

bool ltests::deployERI2D2Kernel()
{
    const double correct_answer = 175567.635699242965;

    lints::Structure structure("def2-svp", "def2-universal-jkfit", atomic_nrs_o3, coords_o3);

    std::vector<lints::ShellData> shell_datas = lints::shellDataAux(structure);

    double eri2_sum = 0;
    for (const lints::ShellData &sh_data_a : shell_datas)
        for (const lints::ShellData &sh_data_b : shell_datas)
        {
            lints::ERI2D2Kernel eri2d2_kernel(sh_data_a, sh_data_b);
            for (size_t ishell_a = 0; ishell_a < sh_data_a.n_shells_; ishell_a++)
                for (size_t ishell_b = 0; ishell_b < sh_data_b.n_shells_; ishell_b++)
                {
                    arr2d<vec2d, 6, 6> eri2_batch =
                            eri2d2_kernel(ishell_a, ishell_b, sh_data_a, sh_data_b);

                    for (int icart = 0; icart < 6; icart++)
                        for (int jcart = 0; jcart < 6; jcart++)
                            for (double i : eri2_batch[icart][jcart])
                                eri2_sum += std::fabs(i);
                }
        }

    if (std::fabs(eri2_sum - correct_answer) < tol)
        return true;

    return false;
}

bool ltests::deployERI4SOCKernel()
{
    const double correct_answer = 12940.953099658431;

    lints::Structure structure("def2-svp", atomic_nrs_o3, coords_o3);

    std::vector<lints::ShellPairData> shell_pair_datas = lints::shellPairData(false, structure);

    double eri4_sum = 0;
    for (const lints::ShellPairData &sp_data_ab : shell_pair_datas)
        for (const lints::ShellPairData &sp_data_cd : shell_pair_datas)
        {
            lints::ERI4SOCKernel eri4soc_kernel(sp_data_ab, sp_data_cd);
            for (size_t ipair1 = 0; ipair1 < sp_data_ab.n_pairs_; ipair1++)
                for (size_t ipair2 = 0; ipair2 < sp_data_cd.n_pairs_; ipair2++)
                {
                    std::array<vec4d, 3> eri4_batch =
                            eri4soc_kernel(ipair1, ipair2, sp_data_ab, sp_data_cd);

                    for (int icart = 0; icart < 3; icart++)
                        for (double val : eri4_batch[icart])
                            eri4_sum += std::fabs(val);
                }
        }

    if (std::fabs(eri4_sum - correct_answer) < tol)
        return true;

    return false;
}

bool ltests::deployERI3SOCKernel()
{
    const double correct_answer = 18124.713152990007;

    lints::Structure structure("def2-svp", "def2-universal-jkfit", atomic_nrs_o3, coords_o3);

    std::vector<lints::ShellPairData> shell_pair_datas = lints::shellPairData(false, structure);
    std::vector<lints::ShellData> shell_datas = lints::shellDataAux(structure);

    double eri3_sum = 0;
    for (const lints::ShellPairData &sp_data_ab : shell_pair_datas)
        for (const lints::ShellData &sh_data_c : shell_datas)
        {
            lints::ERI3SOCKernel eri3soc_kernel(sp_data_ab, sh_data_c);
            for (size_t ipair = 0; ipair < sp_data_ab.n_pairs_; ipair++)
                for (size_t ishell = 0; ishell < sh_data_c.n_shells_; ishell++)
                {
                    std::array<vec3d, 3> eri3_batch =
                            eri3soc_kernel(ipair, ishell, sp_data_ab, sh_data_c);

                    for (int icart = 0; icart < 3; icart++)
                        for (double val : eri3_batch[icart])
                            eri3_sum += std::fabs(val);
                }
        }

    if (std::fabs(eri3_sum - correct_answer) < tol)
        return true;

    return false;
}

bool ltests::ecoeffsRecurrence1()
{
    const double correct_answer = 181021.807437296084;

    vec2d ecoeffs = lints::ecoeffsRecurrence1(19.0034512342, 4);

    double sum_ecoeffs = 0;
    for (double ecoeff : ecoeffs)
        sum_ecoeffs += std::fabs(ecoeff);

    if (std::fabs(sum_ecoeffs - correct_answer) < tol)
        return true;

    return false;
}

bool ltests::ecoeffsRecurrence2()
{
    const double correct_answer = 0.218834692700;

    vec3d ecoeffs = lints::ecoeffsRecurrence2(0.0023445, 1.123621, 2, 3, 0.5113241, 1.2345123,
                                              0.003451234);

    double sum_ecoeffs = 0;
    for (double ecoeff : ecoeffs)
        sum_ecoeffs += std::fabs(ecoeff);

    if (std::fabs(sum_ecoeffs - correct_answer) < tol)
        return true;

    return false;
}

bool ltests::ecoeffsRecurrence2_n1()
{
    const double correct_answer = 0.053908505303;

    double a = 0.0023445;
    double b = 1.123621;
    int la = 2;
    int lb = 3;
    double PA = 0.5113241;
    double PB = 1.2345123;
    double Kab = 0.003451234;
    double A = -0.571235891;
    double B = 1.612357121;

    vec3d ecoeffs = lints::ecoeffsRecurrence2(a, b, la, lb, PA, PB, Kab);
    vec3d ecoeffs1 = lints::ecoeffsRecurrence2_n1(a, b, la, lb, A, B, ecoeffs);

    double sum_ecoeffs = 0;
    for (double ecoeff : ecoeffs1)
        sum_ecoeffs += std::fabs(ecoeff);

    if (std::fabs(sum_ecoeffs - correct_answer) < tol)
        return true;

    return false;
}

bool ltests::ecoeffsPrimitive()
{
    const double correct_answer = 3.002442449382;

    std::array<vec2d, 3> ecoeffs = lints::ecoeffsPrimitive(1230.02342162, 4);

    double sum_ecoeffs = 0;
    for (int icart = 0; icart < 3; icart++)
        for (size_t i = 0; i < ecoeffs[icart].size(); i++)
            sum_ecoeffs += std::fabs(ecoeffs[icart][i]);

    if (std::fabs(sum_ecoeffs - correct_answer) < tol)
        return true;

    return false;
}

bool ltests::ecoeffsPrimitivePair()
{
    const double correct_answer = 16.100457827896;

    std::array<double, 3> xyz_a{0, 1.5, 0};
    std::array<double, 3> xyz_b{1.0, 0, -1.0};

    std::array<vec3d, 3> ecoeffs =
            lints::ecoeffsPrimitivePair(1230.02342162, 0.0023445, 4, 3, xyz_a.data(), xyz_b.data());

    double sum_ecoeffs = 0;
    for (int icart = 0; icart < 3; icart++)
        for (double ecoeff : ecoeffs[icart])
            sum_ecoeffs += std::fabs(ecoeff);

    if (std::fabs(sum_ecoeffs - correct_answer) < tol)
        return true;

    return false;
}

bool ltests::ecoeffsPrimitivePair_n1()
{
    const double correct_answer = 21.712096988812;

    std::array<double, 3> xyz_a{0, 1.5, 0};
    std::array<double, 3> xyz_b{1.0, 0, -1.0};

    double a = 1230.02342162;
    double b = 0.0023445;
    int la = 4;
    int lb = 3;
    double A[3] = {0.571235891, 1.612357121, -0.234512341};
    double B[3] = {1.2345123, -0.5113241, 0.003451234};

    std::array<vec3d, 3> ecoeffs0 =
            lints::ecoeffsPrimitivePair(a, b, la, lb, xyz_a.data(), xyz_b.data());
    std::array<vec3d, 3> ecoeffs1 =
            lints::ecoeffsPrimitivePair_n1(a, b, la, lb, A, B, ecoeffs0);

    double sum_ecoeffs = 0;
    for (int icart = 0; icart < 3; icart++)
        for (double ecoeff : ecoeffs1[icart])
            sum_ecoeffs += std::fabs(ecoeff);

    if (std::fabs(sum_ecoeffs - correct_answer) < tol)
        return true;

    return false;
}

bool ltests::ecoeffsShell()
{
    const double correct_answer = 112231.357645556767;

    std::vector<std::vector<double>> ecoeffs =
            lints::ecoeffsShell(5, {0.95556865, 0.36628794, 0.14651517});

    double sum_ecoeffs = 0;
    for (const auto &ecoeffs_primpair : ecoeffs)
        for (double ecoeff : ecoeffs_primpair)
            sum_ecoeffs += std::fabs(ecoeff);

    if (std::fabs(sum_ecoeffs - correct_answer) < tol)
        return true;

    return false;
}

bool ltests::calcBoysF()
{
    const double correct_answer = 0.728561124180;

    int lab = 5;
    lints::BoysGrid boys_grid(lab);

    std::vector<double> fnx = lints::calcBoysF(lab, 2.5, boys_grid);

    double sum_fnx = 0;
    for (double val : fnx)
        sum_fnx += val;

    if (std::fabs(sum_fnx - correct_answer) < tol)
        return true;

    return false;
}

bool ltests::calcRInts3D()
{
    const double correct_answer = 2.814391945942;

    int lab = 5;
    lints::BoysGrid boys_grid(lab);

    std::vector<double> fnx = lints::calcBoysF(lab, 2.5, boys_grid);

    std::array<double, 3> xyz_ab{0.5, -1.2, 3.4};
    vec3d rints = lints::calcRInts3D(lab, 0.342412347, xyz_ab.data(), fnx.data());

    double sum_rints = 0;
    for (double rint : rints)
        sum_rints += std::fabs(rint);

    if (std::fabs(sum_rints - correct_answer) < tol)
        return true;

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
    const double correct_answer = 4.813408368164;

    double exp = 1.70325241;
    int l = 4;
    double norm = lints::purePrimitiveNorm(l, exp);

    if (std::fabs(norm - correct_answer) < tol)
        return true;

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
        {0, 0, 5}
    };

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
        {5, 5}
    };

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
        {4, 4}
    };

    int l_max = 4;
    std::vector<std::pair<int, int>> l_pairs = lints::getLPairsNoSymm(l_max);

    return (l_pairs == correct_answer);
}
