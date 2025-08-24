#include <lible/utils.hpp>
#include <lible/ints/ecoeffs.hpp>
#include <lible/ints/rints.hpp>
#include <lible/ints/spherical_trafo.hpp>
#include <lible/ints/utils.hpp>
#include <lible/ints/twoel/eri_kernels.hpp>

#include <chrono>
#include <cstring>
#include <format>
#include <map>

namespace LI = lible::ints;

using std::array, std::map, std::pair, std::vector;

namespace lible::ints
{
    vec2d eri4Diagonal(const Structure &structure);

    vec4d eri4(const Structure &structure);

    void eri4Benchmark(const Structure &structure);

    void transferIntsERI4Diag(const int ipair_ab, const ShellPairData &sp_data_ab,
                              const vec4d &eri4_batch, vec2d &eri4_diagonal)
    {
        int dim_a = numSphericals(sp_data_ab.la);
        int dim_b = numSphericals(sp_data_ab.lb);
        int ofs_a = sp_data_ab.offsets_sph[2 * ipair_ab];
        int ofs_b = sp_data_ab.offsets_sph[2 * ipair_ab + 1];
        for (int a = 0; a < dim_a; a++)
            for (int b = 0; b < dim_b; b++)
            {
                double integral = eri4_batch(a, b, a, b);

                int mu = ofs_a + a;
                int nu = ofs_b + b;

                eri4_diagonal(mu, nu) = integral;
                eri4_diagonal(nu, mu) = integral;
            }
    }
}

lible::vec4d LI::eri4(const Structure &structure)
{
    vector<ShellPairData> sp_data = shellPairData(true, structure);

    int dim_ao = structure.getDimAO();
    vec4d eri4(Fill(0), dim_ao);
    for (size_t ispdata_ab = 0; ispdata_ab < sp_data.size(); ispdata_ab++)
        for (size_t isp_data_cd = 0; isp_data_cd <= ispdata_ab; isp_data_cd++)
        {
            const ShellPairData &sp_data_ab = sp_data[ispdata_ab];
            const ShellPairData &sp_data_cd = sp_data[isp_data_cd];

            int n_sph_a = numSphericals(sp_data_ab.la);
            int n_sph_b = numSphericals(sp_data_ab.lb);
            int n_sph_c = numSphericals(sp_data_cd.la);
            int n_sph_d = numSphericals(sp_data_cd.lb);

            ERI4Kernel eri4_kernel = deployERI4Kernel(sp_data_ab, sp_data_cd);

            for (int ipair_ab = 0; ipair_ab < sp_data_ab.n_pairs; ipair_ab++)
            {
                int bound_cd = (ispdata_ab == isp_data_cd) ? ipair_ab + 1 : sp_data_cd.n_pairs;
                for (int ipair_cd = 0; ipair_cd < bound_cd; ipair_cd++)
                {
                    vec4d eri4_batch = eri4_kernel(ipair_ab, ipair_cd, sp_data_ab, sp_data_cd);

                    int ofs_a = sp_data_ab.offsets_sph[2 * ipair_ab];
                    int ofs_b = sp_data_ab.offsets_sph[2 * ipair_ab + 1];
                    int ofs_c = sp_data_cd.offsets_sph[2 * ipair_cd];
                    int ofs_d = sp_data_cd.offsets_sph[2 * ipair_cd + 1];

                    for (int ia = 0; ia < n_sph_a; ia++)
                        for (int ib = 0; ib < n_sph_b; ib++)
                            for (int ic = 0; ic < n_sph_c; ic++)
                                for (int id = 0; id < n_sph_d; id++)
                                {
                                    int mu = ofs_a + ia;
                                    int nu = ofs_b + ib;
                                    int ka = ofs_c + ic;
                                    int ta = ofs_d + id;

                                    double integral = eri4_batch(ia, ib, ic, id);
                                    eri4(mu, nu, ka, ta) = integral;
                                    eri4(mu, nu, ta, ka) = integral;
                                    eri4(nu, mu, ka, ta) = integral;
                                    eri4(nu, mu, ta, ka) = integral;
                                    eri4(ka, ta, mu, nu) = integral;
                                    eri4(ka, ta, nu, mu) = integral;
                                    eri4(ta, ka, mu, nu) = integral;
                                    eri4(ta, ka, nu, mu) = integral;
                                }
                }
            }
        }

    return eri4;
}

void LI::eri4Benchmark(const Structure &structure)
{
    palPrint(std::format("Lible::{:<40}\n", "ERI4 benchmark..."));

    auto start{std::chrono::steady_clock::now()};

    vector<ShellPairData> sp_data = shellPairData(true, structure);

    double sum_eri4 = 0;
    for (size_t ispdata_ab = 0; ispdata_ab < sp_data.size(); ispdata_ab++)
        for (size_t ispdata_cd = 0; ispdata_cd <= ispdata_ab; ispdata_cd++)
        {
            auto start{std::chrono::steady_clock::now()};

            const auto &sp_data_ab = sp_data[ispdata_ab];
            const auto &sp_data_cd = sp_data[ispdata_cd];

            int n_pairs_ab = sp_data_ab.n_pairs;
            int n_pairs_cd = sp_data_cd.n_pairs;

            ERI4Kernel eri4_kernel = deployERI4Kernel(sp_data_ab, sp_data_cd);

            size_t n_shells_abcd = 0;
            for (int ipair_ab = 0; ipair_ab < n_pairs_ab; ipair_ab++)
            {
                int bound_cd = (ispdata_ab == ispdata_cd) ? ipair_ab + 1 : n_pairs_cd;
                for (int ipair_cd = 0; ipair_cd < bound_cd; ipair_cd++)
                {
                    vec4d eri4_batch = eri4_kernel(ipair_ab, ipair_cd, sp_data_ab, sp_data_cd);

                    for (double x : eri4_batch)
                        sum_eri4 += std::fabs(x);

                    n_shells_abcd++;
                }
            }

            auto end{std::chrono::steady_clock::now()};
            std::chrono::duration<double> duration{end - start};

            int la = sp_data_ab.la;
            int lb = sp_data_ab.lb;
            int lc = sp_data_cd.la;
            int ld = sp_data_cd.lb;
            palPrint(std::format("   {} {} {} {} ; {:10} ; {:.2e} s\n", la, lb, lc, ld,
                                 n_shells_abcd, duration.count()));
        }

    palPrint(std::format("   sum_eri4 = {:16.12f}\n", sum_eri4));

    auto end{std::chrono::steady_clock::now()};
    std::chrono::duration<double> duration{end - start};
    palPrint(std::format("done {:.2e} s\n", duration.count()));
}

lible::vec2d LI::eri4Diagonal(const Structure &structure)
{
    vector<ShellPairData> sp_data = shellPairData(true, structure);

    int dim_ao = structure.getDimAO();
    vec2d eri4_diagonal(Fill(0), dim_ao, dim_ao);
    for (size_t ispdata = 0; ispdata < sp_data.size(); ispdata++)
    {
        const auto &sp_data_ab = sp_data[ispdata];

        ERI4Kernel eri4_kernel = deployERI4Kernel(sp_data_ab, sp_data_ab);

        for (int ipair_ab = 0; ipair_ab < sp_data_ab.n_pairs; ipair_ab++)
        {
            vec4d eri4_batch = eri4_kernel(ipair_ab, ipair_ab, sp_data_ab, sp_data_ab);

            transferIntsERI4Diag(ipair_ab, sp_data_ab, eri4_batch, eri4_diagonal);
        }
    }

    return eri4_diagonal;
}