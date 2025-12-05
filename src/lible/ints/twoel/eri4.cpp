#include <lible/utils.hpp>
#include <lible/ints/ints.hpp>
#include <lible/ints/twoel/eri_kernels.hpp>

#include <chrono>
#include <cstring>
#include <format>

namespace lints = lible::ints;

namespace lible::ints
{
    /// Copies the integrals from the shell batch to the target container.
    void transferIntsERI4Diag(int ipair_ab, const ShellPairData &sp_data_ab,
                              const vec4d &eri4_batch, vec2d &eri4_diagonal);
}

void lints::transferIntsERI4Diag(const int ipair_ab, const ShellPairData &sp_data_ab,
                                 const vec4d &eri4_batch, vec2d &eri4_diagonal)
{
    int dim_a = numSphericals(sp_data_ab.la_);
    int dim_b = numSphericals(sp_data_ab.lb_);
    size_t ofs_a = sp_data_ab.offsets_sph_[2 * ipair_ab];
    size_t ofs_b = sp_data_ab.offsets_sph_[2 * ipair_ab + 1];
    for (int a = 0; a < dim_a; a++)
        for (int b = 0; b < dim_b; b++)
        {
            double integral = eri4_batch(a, b, a, b);

            size_t mu = ofs_a + a;
            size_t nu = ofs_b + b;

            eri4_diagonal(mu, nu) = integral;
            eri4_diagonal(nu, mu) = integral;
        }
}

lible::vec4d lints::eri4(const Structure &structure)
{
    std::vector<ShellPairData> sp_data = shellPairData(true, structure);

    size_t dim_ao = structure.getDimAO();
    vec4d eri4(Fill(0), dim_ao);
    for (size_t ispdata_ab = 0; ispdata_ab < sp_data.size(); ispdata_ab++)
        for (size_t isp_data_cd = 0; isp_data_cd <= ispdata_ab; isp_data_cd++)
        {
            const ShellPairData &sp_data_ab = sp_data[ispdata_ab];
            const ShellPairData &sp_data_cd = sp_data[isp_data_cd];

            int n_sph_a = numSphericals(sp_data_ab.la_);
            int n_sph_b = numSphericals(sp_data_ab.lb_);
            int n_sph_c = numSphericals(sp_data_cd.la_);
            int n_sph_d = numSphericals(sp_data_cd.lb_);

            ERI4Kernel eri4_kernel(sp_data_ab, sp_data_cd);

#pragma omp parallel for
            for (size_t ipair_ab = 0; ipair_ab < sp_data_ab.n_pairs_; ipair_ab++)
            {
                size_t bound_cd = (ispdata_ab == isp_data_cd) ? ipair_ab + 1 : sp_data_cd.n_pairs_;
                for (size_t ipair_cd = 0; ipair_cd < bound_cd; ipair_cd++)
                {
                    vec4d eri4_batch = eri4_kernel(ipair_ab, ipair_cd, sp_data_ab, sp_data_cd);

                    size_t ofs_a = sp_data_ab.offsets_sph_[2 * ipair_ab];
                    size_t ofs_b = sp_data_ab.offsets_sph_[2 * ipair_ab + 1];
                    size_t ofs_c = sp_data_cd.offsets_sph_[2 * ipair_cd];
                    size_t ofs_d = sp_data_cd.offsets_sph_[2 * ipair_cd + 1];

                    for (int ia = 0; ia < n_sph_a; ia++)
                        for (int ib = 0; ib < n_sph_b; ib++)
                            for (int ic = 0; ic < n_sph_c; ic++)
                                for (int id = 0; id < n_sph_d; id++)
                                {
                                    size_t mu = ofs_a + ia;
                                    size_t nu = ofs_b + ib;
                                    size_t ka = ofs_c + ic;
                                    size_t ta = ofs_d + id;

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

void lints::eri4Benchmark(const Structure &structure)
{
    palPrint(std::format("Lible::{:<40}\n", "ERI4 benchmark..."));

    auto start_total{std::chrono::steady_clock::now()};

    std::vector<ShellPairData> sp_data = shellPairData(true, structure);

    double sum_eri4 = 0;
    for (size_t ispdata_ab = 0; ispdata_ab < sp_data.size(); ispdata_ab++)
        for (size_t ispdata_cd = 0; ispdata_cd <= ispdata_ab; ispdata_cd++)
        {
            auto start{std::chrono::steady_clock::now()};

            const auto &sp_data_ab = sp_data[ispdata_ab];
            const auto &sp_data_cd = sp_data[ispdata_cd];

            size_t n_pairs_ab = sp_data_ab.n_pairs_;
            size_t n_pairs_cd = sp_data_cd.n_pairs_;

            ERI4Kernel eri4_kernel(sp_data_ab, sp_data_cd);

            size_t n_shells_abcd = 0;
            for (size_t ipair_ab = 0; ipair_ab < n_pairs_ab; ipair_ab++)
            {
                size_t bound_cd = (ispdata_ab == ispdata_cd) ? ipair_ab + 1 : n_pairs_cd;
                for (size_t ipair_cd = 0; ipair_cd < bound_cd; ipair_cd++)
                {
                    vec4d eri4_batch = eri4_kernel(ipair_ab, ipair_cd, sp_data_ab, sp_data_cd);

                    for (double x : eri4_batch)
                        sum_eri4 += std::fabs(x);

                    n_shells_abcd++;
                }
            }

            auto end{std::chrono::steady_clock::now()};
            std::chrono::duration<double> duration{end - start};

            int la = sp_data_ab.la_;
            int lb = sp_data_ab.lb_;
            int lc = sp_data_cd.la_;
            int ld = sp_data_cd.lb_;
            palPrint(std::format("   {} {} {} {} ; {:10} ; {:.2e} s\n", la, lb, lc, ld,
                                 n_shells_abcd, duration.count()));
        }

    palPrint(std::format("   sum_eri4 = {:16.12f}\n", sum_eri4));

    const auto end_total{std::chrono::steady_clock::now()};
    std::chrono::duration<double> duration{end_total - start_total};
    palPrint(std::format("done {:.2e} s\n", duration.count()));
}

lible::vec2d lints::eri4Diagonal(const Structure &structure)
{
    std::vector<ShellPairData> sp_data = shellPairData(true, structure);

    size_t dim_ao = structure.getDimAO();
    vec2d eri4_diagonal(Fill(0), dim_ao, dim_ao);
    for (size_t ispdata = 0; ispdata < sp_data.size(); ispdata++)
    {
        const auto &sp_data_ab = sp_data[ispdata];

        ERI4Kernel eri4_kernel(sp_data_ab, sp_data_ab);

#pragma omp parallel for
        for (size_t ipair_ab = 0; ipair_ab < sp_data_ab.n_pairs_; ipair_ab++)
        {
            vec4d eri4_batch = eri4_kernel(ipair_ab, ipair_ab, sp_data_ab, sp_data_ab);

            transferIntsERI4Diag(ipair_ab, sp_data_ab, eri4_batch, eri4_diagonal);
        }
    }

    return eri4_diagonal;
}
