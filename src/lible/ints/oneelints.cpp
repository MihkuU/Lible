#include <lible/ints/cart_exps.hpp>
#include <lible/ints/ints.hpp>
#include <lible/ints/rints.hpp>
#include <lible/ints/spherical_trafo.hpp>

namespace lints = lible::ints;

namespace lible::ints
{
    /// A helper function for 'kineticEnergyD1Kernel'.
    double kineticEKernelKernel(double b, double b2, double fac, const vec3d &ecoeffs_x,
                                const vec3d &ecoeffs_y, const vec3d &ecoeffs_z,
                                const std::array<int, 3> &ijk, const std::array<int, 3> &i_j_k_);
}

double lints::kineticEKernelKernel(const double b, const double b2, const double fac,
                                   const vec3d &ecoeffs_x, const vec3d &ecoeffs_y,
                                   const vec3d &ecoeffs_z, const std::array<int, 3> &ijk,
                                   const std::array<int, 3> &i_j_k_)
{
    auto [i, j, k] = ijk;
    auto [i_, j_, k_] = i_j_k_;

    double Tx, Ty, Tz;
    if (i_ < 2)
        Tx = -2 * b2 * ecoeffs_x(i, i_ + 2, 0) +
             b * (2 * i_ + 1) * ecoeffs_x(i, i_, 0);
    else
        Tx = -2 * b2 * ecoeffs_x(i, i_ + 2, 0) +
             b * (2 * i_ + 1) * ecoeffs_x(i, i_, 0) -
             0.5 * i_ * (i_ - 1) * ecoeffs_x(i, i_ - 2, 0);

    if (j_ < 2)
        Ty = -2 * b2 * ecoeffs_y(j, j_ + 2, 0) +
             b * (2 * j_ + 1) * ecoeffs_y(j, j_, 0);
    else
        Ty = -2 * b2 * ecoeffs_y(j, j_ + 2, 0) +
             b * (2 * j_ + 1) * ecoeffs_y(j, j_, 0) -
             0.5 * j_ * (j_ - 1) * ecoeffs_y(j, j_ - 2, 0);

    if (k_ < 2)
        Tz = -2 * b2 * ecoeffs_z(k, k_ + 2, 0) +
             b * (2 * k_ + 1) * ecoeffs_z(k, k_, 0);
    else
        Tz = -2 * b2 * ecoeffs_z(k, k_ + 2, 0) +
             b * (2 * k_ + 1) * ecoeffs_z(k, k_, 0) -
             0.5 * k_ * (k_ - 1) * ecoeffs_z(k, k_ - 2, 0);

    double integral = 0;
    integral += fac * Tx * ecoeffs_y(j, j_, 0) * ecoeffs_z(k, k_, 0);
    integral += fac * ecoeffs_x(i, i_, 0) * Ty * ecoeffs_z(k, k_, 0);
    integral += fac * ecoeffs_x(i, i_, 0) * ecoeffs_y(j, j_, 0) * Tz;

    return integral;
}

lible::vec2d lints::overlapKernel(const size_t ipair, const ShellPairData &sp_data)
{
    size_t ofs_prim = sp_data.offsets_primitives_[ipair];
    const double *exps = &sp_data.exps_[ofs_prim];
    const double *coeffs = &sp_data.coeffs_[ofs_prim];
    const double *xyz_a = &sp_data.coords_[6 * ipair + 0];
    const double *xyz_b = &sp_data.coords_[6 * ipair + 3];

    auto [la, lb] = sp_data.getLPair();
    const auto &cart_exps_a = cart_exps[la];
    const auto &cart_exps_b = cart_exps[lb];

    vec2d ints_cart(Fill(0), numCartesians(la), numCartesians(lb));
    for (size_t iab = 0; iab < sp_data.nrs_ppairs_[ipair]; iab++)
    {
        double a = exps[iab * 2];
        double b = exps[iab * 2 + 1];
        double da = coeffs[iab * 2];
        double db = coeffs[iab * 2 + 1];

        double p = a + b;
        double dadb = da * db;
        double fac = dadb * std::pow(M_PI / p, 1.5);

        auto [Ex, Ey, Ez] = ecoeffsPrimitivePair(a, b, la, lb, xyz_a, xyz_b);

        for (const auto &[i, j, k, mu] : cart_exps_a)
            for (const auto &[i_, j_, k_, nu] : cart_exps_b)
                ints_cart(mu, nu) += fac * Ex(i, i_, 0) * Ey(j, j_, 0) * Ez(k, k_, 0);
    }

    vec2d ints_sph = trafo2Spherical(la, lb, ints_cart);

    size_t ofs_norm_a = sp_data.offsets_norms_[2 * ipair + 0];
    size_t ofs_norm_b = sp_data.offsets_norms_[2 * ipair + 1];
    for (size_t mu = 0; mu < ints_sph.dim<0>(); mu++)
        for (size_t nu = 0; nu < ints_sph.dim<1>(); nu++)
        {
            double norm_a = sp_data.norms_[ofs_norm_a + mu];
            double norm_b = sp_data.norms_[ofs_norm_b + nu];
            ints_sph(mu, nu) *= norm_a * norm_b;
        }

    return ints_sph;
}

std::array<lible::vec2d, 6> lints::overlapD1Kernel(const size_t ipair, const ShellPairData &sp_data)
{
    size_t ofs_prim = sp_data.offsets_primitives_[ipair];
    const double *exps = &sp_data.exps_[ofs_prim];
    const double *coeffs = &sp_data.coeffs_[ofs_prim];
    const double *xyz_a = &sp_data.coords_[6 * ipair + 0];
    const double *xyz_b = &sp_data.coords_[6 * ipair + 3];

    int la = sp_data.la_;
    int lb = sp_data.lb_;
    const auto &cart_exps_a = cart_exps[la];
    const auto &cart_exps_b = cart_exps[lb];

    std::array<vec2d, 6> ints_cart;
    for (int ideriv = 0; ideriv < 6; ideriv++)
        ints_cart[ideriv] = vec2d(Fill(0), numCartesians(la), numCartesians(lb));

    for (size_t iab = 0; iab < sp_data.nrs_ppairs_[ipair]; iab++)
    {
        double a = exps[iab * 2];
        double b = exps[iab * 2 + 1];
        double da = coeffs[iab * 2];
        double db = coeffs[iab * 2 + 1];

        double p = a + b;
        double dadb = da * db;
        double fac = dadb * std::pow(M_PI / p, 1.5);

        auto [Ex, Ey, Ez] = ecoeffsPrimitivePair(a, b, la, lb, xyz_a, xyz_b);

        auto [E1x, E1y, E1z] = ecoeffsPrimitivePair_n1(a, b, la, lb, xyz_a, xyz_b, {Ex, Ey, Ez});

        for (const auto &[i, j, k, mu] : cart_exps_a)
            for (const auto &[i_, j_, k_, nu] : cart_exps_b)
            {
                // d/dA
                ints_cart[0](mu, nu) += fac * E1x(i, i_, 0) * Ey(j, j_, 0) * Ez(k, k_, 0);
                ints_cart[1](mu, nu) += fac * Ex(i, i_, 0) * E1y(j, j_, 0) * Ez(k, k_, 0);
                ints_cart[2](mu, nu) += fac * Ex(i, i_, 0) * Ey(j, j_, 0) * E1z(k, k_, 0);

                // d/dB
                ints_cart[3](mu, nu) -= fac * E1x(i, i_, 0) * Ey(j, j_, 0) * Ez(k, k_, 0);
                ints_cart[4](mu, nu) -= fac * Ex(i, i_, 0) * E1y(j, j_, 0) * Ez(k, k_, 0);
                ints_cart[5](mu, nu) -= fac * Ex(i, i_, 0) * Ey(j, j_, 0) * E1z(k, k_, 0);
            }
    }

    std::array<vec2d, 6> ints_sph;
    for (int ideriv = 0; ideriv < 6; ideriv++)
        ints_sph[ideriv] = trafo2Spherical(la, lb, ints_cart[ideriv]);

    size_t ofs_norm_a = sp_data.offsets_norms_[2 * ipair + 0];
    size_t ofs_norm_b = sp_data.offsets_norms_[2 * ipair + 1];
    for (int ideriv = 0; ideriv < 6; ideriv++)
        for (size_t mu = 0; mu < ints_sph[ideriv].dim<0>(); mu++)
            for (size_t nu = 0; nu < ints_sph[ideriv].dim<1>(); nu++)
            {
                double norm_a = sp_data.norms_[ofs_norm_a + mu];
                double norm_b = sp_data.norms_[ofs_norm_b + nu];
                ints_sph[ideriv](mu, nu) *= norm_a * norm_b;
            }

    return ints_sph;
}


lible::vec2d lints::overlap(const Structure &structure)
{
    int l_max = structure.getMaxL();
    size_t dim_ao = structure.getDimAO();

    vec2d ints(Fill(0), dim_ao, dim_ao);
    for (int la = l_max; la >= 0; la--)
        for (int lb = la; lb >= 0; lb--)
        {
            ShellPairData sp_data(true, la, lb, structure.getShellsL(la), structure.getShellsL(lb));

#pragma omp parallel for
            for (size_t ipair = 0; ipair < sp_data.n_pairs_; ipair++)
            {
                vec2d ints_ipair = overlapKernel(ipair, sp_data);

                size_t ofs_a = sp_data.offsets_sph_[2 * ipair + 0];
                size_t ofs_b = sp_data.offsets_sph_[2 * ipair + 1];
                for (size_t mu = 0; mu < ints_ipair.dim<0>(); mu++)
                    for (size_t nu = 0; nu < ints_ipair.dim<1>(); nu++)
                    {
                        ints(ofs_a + mu, ofs_b + nu) = ints_ipair(mu, nu);
                        ints(ofs_b + nu, ofs_a + mu) = ints_ipair(mu, nu);
                    }
            }
        }

    return ints;
}

lible::vec2d lints::kineticEnergyKernel(const size_t ipair, const ShellPairData &sp_data)
{
    // Formula taken from https://gqcg-res.github.io/knowdes/the-mcmurchie-davidson-integral-scheme.html.

    size_t ofs_prim = sp_data.offsets_primitives_[ipair];
    const double *exps = &sp_data.exps_[ofs_prim];
    const double *coeffs = &sp_data.coeffs_[ofs_prim];
    const double *xyz_a = &sp_data.coords_[6 * ipair + 0];
    const double *xyz_b = &sp_data.coords_[6 * ipair + 3];

    int la = sp_data.la_;
    int lb = sp_data.lb_;
    const auto &cart_exps_a = cart_exps[la];
    const auto &cart_exps_b = cart_exps[lb];

    vec2d ints_cart(Fill(0), numCartesians(la), numCartesians(lb));
    for (size_t iab = 0; iab < sp_data.nrs_ppairs_[ipair]; iab++)
    {
        double a = exps[iab * 2];
        double b = exps[iab * 2 + 1];
        double b2 = b * b;
        double da = coeffs[iab * 2];
        double db = coeffs[iab * 2 + 1];

        double p = a + b;
        double dadb = da * db;
        double fac = dadb * std::pow(M_PI / p, 1.5);

        auto [Ex, Ey, Ez] = ecoeffsPrimitivePair(a, b, la, lb + 2, xyz_a, xyz_b);

        for (const auto &[i, j, k, mu] : cart_exps_a)
            for (const auto &[i_, j_, k_, nu] : cart_exps_b)
            {
                double Tx, Ty, Tz;
                if (i_ < 2)
                    Tx = -2 * b2 * Ex(i, i_ + 2, 0) +
                         b * (2 * i_ + 1) * Ex(i, i_, 0);
                else
                    Tx = -2 * b2 * Ex(i, i_ + 2, 0) +
                         b * (2 * i_ + 1) * Ex(i, i_, 0) -
                         0.5 * i_ * (i_ - 1) * Ex(i, i_ - 2, 0);

                if (j_ < 2)
                    Ty = -2 * b2 * Ey(j, j_ + 2, 0) +
                         b * (2 * j_ + 1) * Ey(j, j_, 0);
                else
                    Ty = -2 * b2 * Ey(j, j_ + 2, 0) +
                         b * (2 * j_ + 1) * Ey(j, j_, 0) -
                         0.5 * j_ * (j_ - 1) * Ey(j, j_ - 2, 0);

                if (k_ < 2)
                    Tz = -2 * b2 * Ez(k, k_ + 2, 0) +
                         b * (2 * k_ + 1) * Ez(k, k_, 0);
                else
                    Tz = -2 * b2 * Ez(k, k_ + 2, 0) +
                         b * (2 * k_ + 1) * Ez(k, k_, 0) -
                         0.5 * k_ * (k_ - 1) * Ez(k, k_ - 2, 0);

                ints_cart(mu, nu) += fac * Tx * Ey(j, j_, 0) * Ez(k, k_, 0);
                ints_cart(mu, nu) += fac * Ex(i, i_, 0) * Ty * Ez(k, k_, 0);
                ints_cart(mu, nu) += fac * Ex(i, i_, 0) * Ey(j, j_, 0) * Tz;
            }
    }

    vec2d ints_sph = trafo2Spherical(la, lb, ints_cart);

    size_t ofs_norm_a = sp_data.offsets_norms_[2 * ipair + 0];
    size_t ofs_norm_b = sp_data.offsets_norms_[2 * ipair + 1];
    for (size_t mu = 0; mu < ints_sph.dim<0>(); mu++)
        for (size_t nu = 0; nu < ints_sph.dim<1>(); nu++)
        {
            double norm_a = sp_data.norms_[ofs_norm_a + mu];
            double norm_b = sp_data.norms_[ofs_norm_b + nu];
            ints_sph(mu, nu) *= norm_a * norm_b;
        }

    return ints_sph;
}

std::array<lible::vec2d, 6> lints::kineticEnergyD1Kernel(const size_t ipair, const ShellPairData &sp_data)
{
    size_t ofs_prim = sp_data.offsets_primitives_[ipair];
    const double *exps = &sp_data.exps_[ofs_prim];
    const double *coeffs = &sp_data.coeffs_[ofs_prim];
    const double *xyz_a = &sp_data.coords_[6 * ipair + 0];
    const double *xyz_b = &sp_data.coords_[6 * ipair + 3];

    auto [la, lb] = sp_data.getLPair();
    const auto &cart_exps_a = cart_exps[la];
    const auto &cart_exps_b = cart_exps[lb];

    std::array<vec2d, 6> ints_cart;
    for (int ideriv = 0; ideriv < 6; ideriv++)
        ints_cart[ideriv] = vec2d(Fill(0), numCartesians(la), numCartesians(lb));

    for (size_t iab = 0; iab < sp_data.nrs_ppairs_[ipair]; iab++)
    {
        double a = exps[iab * 2];
        double b = exps[iab * 2 + 1];
        double b2 = b * b;
        double da = coeffs[iab * 2];
        double db = coeffs[iab * 2 + 1];

        double p = a + b;
        double dadb = da * db;
        double fac = dadb * std::pow(M_PI / p, 1.5);

        auto [Ex, Ey, Ez] = ecoeffsPrimitivePair(a, b, la, lb + 2, xyz_a, xyz_b);

        auto [E1x, E1y, E1z] = ecoeffsPrimitivePair_n1(a, b, la, lb + 2, xyz_a, xyz_b,
                                                       {Ex, Ey, Ez});

        for (const auto &[i, j, k, mu] : cart_exps_a)
            for (const auto &[i_, j_, k_, nu] : cart_exps_b)
            {
                double kin_x = kineticEKernelKernel(b, b2, fac, E1x, Ey, Ez, {i, j, k}, {i_, j_, k_});
                double kin_y = kineticEKernelKernel(b, b2, fac, Ex, E1y, Ez, {i, j, k}, {i_, j_, k_});
                double kin_z = kineticEKernelKernel(b, b2, fac, Ex, Ey, E1z, {i, j, k}, {i_, j_, k_});

                // d/dA
                ints_cart[0](mu, nu) += kin_x;
                ints_cart[1](mu, nu) += kin_y;
                ints_cart[2](mu, nu) += kin_z;

                // d/dB
                ints_cart[3](mu, nu) -= kin_x;
                ints_cart[4](mu, nu) -= kin_y;
                ints_cart[5](mu, nu) -= kin_z;
            }
    }

    std::array<vec2d, 6> ints_sph;
    for (int ideriv = 0; ideriv < 6; ideriv++)
        ints_sph[ideriv] = trafo2Spherical(la, lb, ints_cart[ideriv]);

    size_t ofs_norm_a = sp_data.offsets_norms_[2 * ipair + 0];
    size_t ofs_norm_b = sp_data.offsets_norms_[2 * ipair + 1];
    for (int ideriv = 0; ideriv < 6; ideriv++)
        for (size_t mu = 0; mu < ints_sph[ideriv].dim<0>(); mu++)
            for (size_t nu = 0; nu < ints_sph[ideriv].dim<1>(); nu++)
            {
                double norm_a = sp_data.norms_[ofs_norm_a + mu];
                double norm_b = sp_data.norms_[ofs_norm_b + nu];
                ints_sph[ideriv](mu, nu) *= norm_a * norm_b;
            }

    return ints_sph;
}

lible::vec2d lints::kineticEnergy(const Structure &structure)
{
    int l_max = structure.getMaxL();
    size_t dim_ao = structure.getDimAO();

    vec2d ints(Fill(0), dim_ao, dim_ao);
    for (int la = l_max; la >= 0; la--)
        for (int lb = la; lb >= 0; lb--)
        {
            ShellPairData sp_data(true, la, lb, structure.getShellsL(la), structure.getShellsL(lb));

#pragma omp parallel for
            for (size_t ipair = 0; ipair < sp_data.n_pairs_; ipair++)
            {
                vec2d ints_ipair = kineticEnergyKernel(ipair, sp_data);

                size_t ofs_a = sp_data.offsets_sph_[2 * ipair + 0];
                size_t ofs_b = sp_data.offsets_sph_[2 * ipair + 1];
                for (size_t mu = 0; mu < ints_ipair.dim<0>(); mu++)
                    for (size_t nu = 0; nu < ints_ipair.dim<1>(); nu++)
                    {
                        ints(ofs_a + mu, ofs_b + nu) = ints_ipair(mu, nu);
                        ints(ofs_b + nu, ofs_a + mu) = ints_ipair(mu, nu);
                    }
            }
        }

    return ints;
}

std::array<lible::vec2d, 3>
lints::dipoleMomentKernel(const size_t ipair, const std::array<double, 3> &origin,
                          const ShellPairData &sp_data)
{
    size_t ofs_prim = sp_data.offsets_primitives_[ipair];
    const double *exps = &sp_data.exps_[ofs_prim];
    const double *coeffs = &sp_data.coeffs_[ofs_prim];
    const double *xyz_a = &sp_data.coords_[6 * ipair + 0];
    const double *xyz_b = &sp_data.coords_[6 * ipair + 3];

    auto [la, lb] = sp_data.getLPair();
    int lab = la + lb;
    const auto &cart_exps_a = cart_exps[la];
    const auto &cart_exps_b = cart_exps[lb];

    std::array<vec2d, 3> ints_cart;
    for (int icart = 0; icart < 3; icart++)
        ints_cart[icart] = vec2d(Fill(0), numCartesians(la), numCartesians(lb));

    for (size_t iab = 0; iab < sp_data.nrs_ppairs_[ipair]; iab++)
    {
        double a = exps[iab * 2];
        double b = exps[iab * 2 + 1];
        double da = coeffs[iab * 2];
        double db = coeffs[iab * 2 + 1];

        double p = a + b;
        double dadb = da * db;
        double fac = dadb * std::pow(M_PI / p, 1.5);

        auto [Ex, Ey, Ez] = ecoeffsPrimitivePair(a, b, la, lb, xyz_a, xyz_b);

        std::array<double, 3> xyz_p{
            (a * xyz_a[0] + b * xyz_b[0]) / p,
            (a * xyz_a[1] + b * xyz_b[1]) / p,
            (a * xyz_a[2] + b * xyz_b[2]) / p
        };

        std::array<double, 3> xyz_po{
            xyz_p[0] - origin[0],
            xyz_p[1] - origin[1],
            xyz_p[2] - origin[2]
        };

        for (const auto &[i, j, k, mu] : cart_exps_a)
            for (const auto &[i_, j_, k_, nu] : cart_exps_b)
            {
                double dx, dy, dz;
                if (lab > 0)
                {
                    dx = fac * (Ex(i, i_, 1) + xyz_po[0] * Ex(i, i_, 0)) *
                         Ey(j, j_, 0) * Ez(k, k_, 0);

                    dy = fac * (Ex(i, i_, 0)) *
                         (Ey(j, j_, 1) + xyz_po[1] * Ey(j, j_, 0)) *
                         Ez(k, k_, 0);

                    dz = fac * (Ex(i, i_, 0)) * Ey(j, j_, 0) *
                         (Ez(k, k_, 1) + xyz_po[2] * Ez(k, k_, 0));
                }
                else
                {
                    dx = fac * (xyz_po[0] * Ex(i, i_, 0)) *
                         Ey(j, j_, 0) * Ez(k, k_, 0);

                    dy = fac * (Ex(i, i_, 0)) *
                         (xyz_po[1] * Ey(j, j_, 0)) *
                         Ez(k, k_, 0);

                    dz = fac * (Ex(i, i_, 0)) * Ey(j, j_, 0) *
                         (xyz_po[2] * Ez(k, k_, 0));
                }

                ints_cart[0](mu, nu) += dx;
                ints_cart[1](mu, nu) += dy;
                ints_cart[2](mu, nu) += dz;
            }
    }

    std::array<vec2d, 3> ints_sph;
    for (int icart = 0; icart < 3; icart++)
        ints_sph[icart] = trafo2Spherical(la, lb, ints_cart[icart]);

    size_t ofs_norm_a = sp_data.offsets_norms_[2 * ipair + 0];
    size_t ofs_norm_b = sp_data.offsets_norms_[2 * ipair + 1];
    for (int icart = 0; icart < 3; icart++)
        for (size_t mu = 0; mu < ints_sph[icart].dim<0>(); mu++)
            for (size_t nu = 0; nu < ints_sph[icart].dim<1>(); nu++)
            {
                double norm_a = sp_data.norms_[ofs_norm_a + mu];
                double norm_b = sp_data.norms_[ofs_norm_b + nu];
                ints_sph[icart](mu, nu) *= norm_a * norm_b;
            }

    return ints_sph;
}

std::array<lible::vec2d, 3> lints::dipoleMoment(const std::array<double, 3> &origin,
                                                const Structure &structure)
{
    int l_max = structure.getMaxL();
    size_t dim_ao = structure.getDimAO();

    vec2d filler(Fill(0), dim_ao, dim_ao);
    std::array<vec2d, 3> ints{filler, filler, filler};
    for (int la = l_max; la >= 0; la--)
        for (int lb = la; lb >= 0; lb--)
        {
            ShellPairData sp_data(true, la, lb, structure.getShellsL(la), structure.getShellsL(lb));

#pragma omp parallel for
            for (size_t ipair = 0; ipair < sp_data.n_pairs_; ipair++)
            {
                std::array<vec2d, 3> ints_ipair = dipoleMomentKernel(ipair, origin, sp_data);

                size_t ofs_a = sp_data.offsets_sph_[2 * ipair + 0];
                size_t ofs_b = sp_data.offsets_sph_[2 * ipair + 1];
                for (size_t mu = 0; mu < ints_ipair[0].dim<0>(); mu++)
                    for (size_t nu = 0; nu < ints_ipair[0].dim<1>(); nu++)
                    {
                        for (int icart = 0; icart < 3; icart++)
                        {
                            ints[icart](ofs_a + mu, ofs_b + nu) = ints_ipair[icart](mu, nu);
                            ints[icart](ofs_b + nu, ofs_a + mu) = ints_ipair[icart](mu, nu);
                        }
                    }
            }
        }

    return ints;
}

lible::vec2d
lints::externalChargesKernel(const size_t ipair, const std::vector<std::array<double, 4>> &charges,
                             const BoysGrid &boys_grid, const ShellPairData &sp_data)
{
    size_t ofs_prim = sp_data.offsets_primitives_[ipair];
    const double *exps = &sp_data.exps_[ofs_prim];
    const double *coeffs = &sp_data.coeffs_[ofs_prim];
    const double *xyz_a = &sp_data.coords_[6 * ipair + 0];
    const double *xyz_b = &sp_data.coords_[6 * ipair + 3];

    auto [la, lb] = sp_data.getLPair();
    int lab = la + lb;
    const auto &cart_exps_a = cart_exps[la];
    const auto &cart_exps_b = cart_exps[lb];

    vec2d ints_cart(Fill(0), numCartesians(la), numCartesians(lb));
    for (size_t iab = 0; iab < sp_data.nrs_ppairs_[ipair]; iab++)
    {
        double a = exps[iab * 2];
        double b = exps[iab * 2 + 1];
        double da = coeffs[iab * 2];
        double db = coeffs[iab * 2 + 1];

        double p = a + b;
        double dadb = da * db;
        double fac = 2 * (M_PI / p) * dadb;

        auto [Ex, Ey, Ez] = ecoeffsPrimitivePair(a, b, la, lb, xyz_a, xyz_b);

        std::array<double, 3> xyz_p{
            (a * xyz_a[0] + b * xyz_b[0]) / p,
            (a * xyz_a[1] + b * xyz_b[1]) / p,
            (a * xyz_a[2] + b * xyz_b[2]) / p
        };

        vec3d rints_sum(Fill(0), lab + 1);
        for (auto [xc, yc, zc, charge] : charges)
        {
            std::array<double, 3> xyz_pc{xyz_p[0] - xc, xyz_p[1] - yc, xyz_p[2] - zc};

            double xx{xyz_pc[0]}, xy{xyz_pc[1]}, xz{xyz_pc[2]};
            double xyz_pc_dot = xx * xx + xy * xy + xz * xz;
            double x = p * xyz_pc_dot;

            std::vector<double> fnx = calcBoysF(lab, x, boys_grid);

            vec3d rints = calcRInts3D(lab, p, &xyz_pc[0], &fnx[0]);

            for (int t = 0; t <= lab; t++)
                for (int u = 0; u <= lab; u++)
                    for (int v = 0; v <= lab; v++)
                        rints_sum(t, u, v) += charge * rints(t, u, v);
        }

        for (const auto &[i, j, k, mu] : cart_exps_a)
            for (const auto &[i_, j_, k_, nu] : cart_exps_b)
                for (int t = 0; t <= i + i_; t++)
                    for (int u = 0; u <= j + j_; u++)
                        for (int v = 0; v <= k + k_; v++)
                            ints_cart(mu, nu) += (-1) * fac * // -1 = charge of electron
                                    Ex(i, i_, t) * Ey(j, j_, u) * Ez(k, k_, v) *
                                    rints_sum(t, u, v);
    }

    vec2d ints_sph = trafo2Spherical(la, lb, ints_cart);

    size_t ofs_norm_a = sp_data.offsets_norms_[2 * ipair + 0];
    size_t ofs_norm_b = sp_data.offsets_norms_[2 * ipair + 1];
    for (size_t mu = 0; mu < ints_sph.dim<0>(); mu++)
        for (size_t nu = 0; nu < ints_sph.dim<1>(); nu++)
        {
            double norm_a = sp_data.norms_[ofs_norm_a + mu];
            double norm_b = sp_data.norms_[ofs_norm_b + nu];
            ints_sph(mu, nu) *= norm_a * norm_b;
        }

    return ints_sph;
}

lible::vec2d lints::externalChargesErfKernel(const size_t ipair,
                                             const std::vector<std::array<double, 4>> &charges,
                                             const std::vector<double> &omegas,
                                             const BoysGrid &boys_grid,
                                             const ShellPairData &sp_data)
{
    if (omegas.size() != charges.size())
        throw std::runtime_error("externalChargesErfKernel(): number of omega values must match "
            "number of charges");

    auto [la, lb] = sp_data.getLPair();
    int lab = la + lb;

    size_t ofs_prim = sp_data.offsets_primitives_[ipair];
    const double *exps = &sp_data.exps_[ofs_prim];
    const double *coeffs = &sp_data.coeffs_[ofs_prim];
    const double *xyz_a = &sp_data.coords_[6 * ipair + 0];
    const double *xyz_b = &sp_data.coords_[6 * ipair + 3];

    const auto &cart_exps_a = cart_exps[la];
    const auto &cart_exps_b = cart_exps[lb];

    vec2d ints_cart(Fill(0), numCartesians(la), numCartesians(lb));
    for (size_t iab = 0; iab < sp_data.nrs_ppairs_[ipair]; iab++)
    {
        double a = exps[iab * 2];
        double b = exps[iab * 2 + 1];
        double da = coeffs[iab * 2];
        double db = coeffs[iab * 2 + 1];

        double p = a + b;
        double dadb = da * db;
        double fac = 2 * (M_PI / p) * dadb;

        auto [Ex, Ey, Ez] = ecoeffsPrimitivePair(a, b, la, lb, xyz_a, xyz_b);

        std::array<double, 3> xyz_p{
            (a * xyz_a[0] + b * xyz_b[0]) / p,
            (a * xyz_a[1] + b * xyz_b[1]) / p,
            (a * xyz_a[2] + b * xyz_b[2]) / p
        };

        vec3d rints_sum(Fill(0), lab + 1);
        for (size_t icharge = 0; icharge < charges.size(); icharge++)
        {
            auto [xc, yc, zc, charge] = charges[icharge];

            std::array<double, 3> xyz_pc{xyz_p[0] - xc, xyz_p[1] - yc, xyz_p[2] - zc};

            double xx{xyz_pc[0]}, xy{xyz_pc[1]}, xz{xyz_pc[2]};
            double xyz_pc_dot = xx * xx + xy * xy + xz * xz;

            double omega_factor = omegas[icharge] / std::pow(omegas[icharge] * omegas[icharge] + p, 0.5);

            double omega_squared = omegas[icharge] * omegas[icharge];

            // second part of (52) in https://doi.org/10.1039/B605188J
            double x = p * xyz_pc_dot * omega_squared / (omega_squared + p);

            std::vector<double> fnx = calcBoysF(lab, x, boys_grid);

            vec3d rints = calcRInts3DErf(lab, p, omegas[icharge], &xyz_pc[0], &fnx[0]);

            for (int t = 0; t <= lab; t++)
                for (int u = 0; u <= lab; u++)
                    for (int v = 0; v <= lab; v++)
                        rints_sum(t, u, v) += charge * omega_factor // erf-related factors
                                * rints(t, u, v);
        }

        for (const auto &[i, j, k, mu] : cart_exps_a)
            for (const auto &[i_, j_, k_, nu] : cart_exps_b)
                for (int t = 0; t <= i + i_; t++)
                    for (int u = 0; u <= j + j_; u++)
                        for (int v = 0; v <= k + k_; v++)
                            ints_cart(mu, nu) += (-1) * fac * // -1 = charge of electron
                                    Ex(i, i_, t) * Ey(j, j_, u) * Ez(k, k_, v) *
                                    rints_sum(t, u, v);
    }

    vec2d ints_sph = trafo2Spherical(la, lb, ints_cart);

    size_t ofs_norm_a = sp_data.offsets_norms_[2 * ipair + 0];
    size_t ofs_norm_b = sp_data.offsets_norms_[2 * ipair + 1];
    for (size_t mu = 0; mu < ints_sph.dim<0>(); mu++)
        for (size_t nu = 0; nu < ints_sph.dim<1>(); nu++)
        {
            double norm_a = sp_data.norms_[ofs_norm_a + mu];
            double norm_b = sp_data.norms_[ofs_norm_b + nu];
            ints_sph(mu, nu) *= norm_a * norm_b;
        }

    return ints_sph;
}

std::array<lible::vec2d, 6>
lints::externalChargesD1Kernel(const size_t ipair, const std::vector<std::array<double, 4>> &charges,
                               const BoysGrid &boys_grid, const ShellPairData &sp_data)
{
    auto [la, lb] = sp_data.getLPair();
    const int lab = la + lb;

    size_t ofs_prim = sp_data.offsets_primitives_[ipair];
    const double *exps = &sp_data.exps_[ofs_prim];
    const double *coeffs = &sp_data.coeffs_[ofs_prim];
    const double *xyz_a = &sp_data.coords_[6 * ipair + 0];
    const double *xyz_b = &sp_data.coords_[6 * ipair + 3];

    const auto &cart_exps_a = cart_exps[la];
    const auto &cart_exps_b = cart_exps[lb];

    std::array<vec2d, 6> ints_cart;
    for (int ideriv = 0; ideriv < 6; ideriv++)
        ints_cart[ideriv] = vec2d(Fill(0), numCartesians(la), numCartesians(lb));

    for (size_t iab = 0; iab < sp_data.nrs_ppairs_[ipair]; iab++)
    {
        double a = exps[iab * 2];
        double b = exps[iab * 2 + 1];
        double da = coeffs[iab * 2];
        double db = coeffs[iab * 2 + 1];

        double p = a + b;
        double dadb = da * db;
        double fac = 2 * (M_PI / p) * dadb;

        auto [Ex, Ey, Ez] = ecoeffsPrimitivePair(a, b, la, lb, xyz_a, xyz_b);

        auto [E1x, E1y, E1z] = ecoeffsPrimitivePair_n1(a, b, la, lb, xyz_a, xyz_b,
                                                       {Ex, Ey, Ez});

        std::array<double, 3> xyz_p{
            (a * xyz_a[0] + b * xyz_b[0]) / p,
            (a * xyz_a[1] + b * xyz_b[1]) / p,
            (a * xyz_a[2] + b * xyz_b[2]) / p
        };

        vec3d rints_sum(Fill(0), lab + 2);
        for (auto [xc, yc, zc, charge] : charges)
        {
            std::array<double, 3> xyz_pc{xyz_p[0] - xc, xyz_p[1] - yc, xyz_p[2] - zc};

            double xx{xyz_pc[0]}, xy{xyz_pc[1]}, xz{xyz_pc[2]};
            double xyz_pc_dot = xx * xx + xy * xy + xz * xz;
            double x = p * xyz_pc_dot;

            std::vector<double> fnx = calcBoysF(lab + 1, x, boys_grid);

            vec3d rints = calcRInts3D(lab + 1, p, &xyz_pc[0], &fnx[0]);

            for (int t = 0; t <= lab + 1; t++)
                for (int u = 0; u <= lab + 1; u++)
                    for (int v = 0; v <= lab + 1; v++)
                        rints_sum(t, u, v) += charge * rints(t, u, v);
        }

        for (const auto &[i, j, k, mu] : cart_exps_a)
            for (const auto &[i_, j_, k_, nu] : cart_exps_b)
                for (int t = 0; t <= i + i_; t++)
                    for (int u = 0; u <= j + j_; u++)
                        for (int v = 0; v <= k + k_; v++)
                        {
                            double dpx = fac * Ex(i, i_, t) * Ey(j, j_, u) * Ez(k, k_, v) * rints_sum(t + 1, u, v);
                            double dpy = fac * Ex(i, i_, t) * Ey(j, j_, u) * Ez(k, k_, v) * rints_sum(t, u + 1, v);
                            double dpz = fac * Ex(i, i_, t) * Ey(j, j_, u) * Ez(k, k_, v) * rints_sum(t, u, v + 1);

                            double drx = fac * E1x(i, i_, t) * Ey(j, j_, u) * Ez(k, k_, v) * rints_sum(t, u, v);
                            double dry = fac * Ex(i, i_, t) * E1y(j, j_, u) * Ez(k, k_, v) * rints_sum(t, u, v);
                            double drz = fac * Ex(i, i_, t) * Ey(j, j_, u) * E1z(k, k_, v) * rints_sum(t, u, v);

                            // d/dA
                            ints_cart[0](mu, nu) += -1 * ((a / p) * dpx + drx); // -1 = charge of electron
                            ints_cart[1](mu, nu) += -1 * ((a / p) * dpy + dry);
                            ints_cart[2](mu, nu) += -1 * ((a / p) * dpz + drz);

                            // d/dB
                            ints_cart[3](mu, nu) += -1 * ((b / p) * dpx - drx);
                            ints_cart[4](mu, nu) += -1 * ((b / p) * dpy - dry);
                            ints_cart[5](mu, nu) += -1 * ((b / p) * dpz - drz);
                        }
    }

    std::array<vec2d, 6> ints_sph;
    for (int ideriv = 0; ideriv < 6; ideriv++)
        ints_sph[ideriv] = trafo2Spherical(la, lb, ints_cart[ideriv]);

    size_t ofs_norm_a = sp_data.offsets_norms_[2 * ipair + 0];
    size_t ofs_norm_b = sp_data.offsets_norms_[2 * ipair + 1];
    for (int ideriv = 0; ideriv < 6; ideriv++)
        for (size_t mu = 0; mu < ints_sph[ideriv].dim<0>(); mu++)
            for (size_t nu = 0; nu < ints_sph[ideriv].dim<1>(); nu++)
            {
                double norm_a = sp_data.norms_[ofs_norm_a + mu];
                double norm_b = sp_data.norms_[ofs_norm_b + nu];
                ints_sph[ideriv](mu, nu) *= norm_a * norm_b;
            }

    return ints_sph;
}

std::vector<std::array<lible::vec2d, 3>>
lints::externalChargesOperatorD1Kernel(const size_t ipair, const std::vector<std::array<double, 4>> &charges,
                                       const BoysGrid &boys_grid, const ShellPairData &sp_data)
{
    auto [la, lb] = sp_data.getLPair();
    int lab = la + lb;

    size_t ofs_prim = sp_data.offsets_primitives_[ipair];
    const double *exps = &sp_data.exps_[ofs_prim];
    const double *coeffs = &sp_data.coeffs_[ofs_prim];
    const double *xyz_a = &sp_data.coords_[6 * ipair + 0];
    const double *xyz_b = &sp_data.coords_[6 * ipair + 3];

    const auto &cart_exps_a = cart_exps[la];
    const auto &cart_exps_b = cart_exps[lb];

    size_t n_charges = charges.size();
    std::vector<std::array<vec2d, 3>> ints_cart(n_charges);
    for (size_t icharge = 0; icharge < n_charges; icharge++)
        for (int icoord = 0; icoord < 3; icoord++)
            ints_cart[icharge][icoord] = vec2d(Fill(0), numCartesians(la), numCartesians(lb));

    for (size_t iab = 0; iab < sp_data.nrs_ppairs_[ipair]; iab++)
    {
        double a = exps[iab * 2];
        double b = exps[iab * 2 + 1];
        double da = coeffs[iab * 2];
        double db = coeffs[iab * 2 + 1];

        double p = a + b;
        double dadb = da * db;
        double fac = 2 * (M_PI / p) * dadb;

        auto [Ex, Ey, Ez] = ecoeffsPrimitivePair(a, b, la, lb, xyz_a, xyz_b);

        auto [E1x, E1y, E1z] = ecoeffsPrimitivePair_n1(a, b, la, lb, xyz_a, xyz_b,
                                                       {Ex, Ey, Ez});

        std::array<double, 3> xyz_p{
            (a * xyz_a[0] + b * xyz_b[0]) / p,
            (a * xyz_a[1] + b * xyz_b[1]) / p,
            (a * xyz_a[2] + b * xyz_b[2]) / p
        };

        for (size_t icharge = 0; icharge < n_charges; icharge++)
        {
            auto [xc, yc, zc, charge] = charges[icharge];

            std::array<double, 3> xyz_pc{xyz_p[0] - xc, xyz_p[1] - yc, xyz_p[2] - zc};

            double xx{xyz_pc[0]}, xy{xyz_pc[1]}, xz{xyz_pc[2]};
            double xyz_pc_dot = xx * xx + xy * xy + xz * xz;
            double x = p * xyz_pc_dot;

            std::vector<double> fnx = calcBoysF(lab + 1, x, boys_grid);

            vec3d rints = calcRInts3D(lab + 1, p, &xyz_pc[0], &fnx[0]);

            for (const auto &[i, j, k, mu] : cart_exps_a)
                for (const auto &[i_, j_, k_, nu] : cart_exps_b)
                    for (int t = 0; t <= i + i_; t++)
                        for (int u = 0; u <= j + j_; u++)
                            for (int v = 0; v <= k + k_; v++)
                            {
                                double Exyz = Ex(i, i_, t) * Ey(j, j_, u) * Ez(k, k_, v);
                                ints_cart[icharge][0](mu, nu) += charge * fac * Exyz * rints(t + 1, u, v);
                                ints_cart[icharge][1](mu, nu) += charge * fac * Exyz * rints(t, u + 1, v);
                                ints_cart[icharge][2](mu, nu) += charge * fac * Exyz * rints(t, u, v + 1);
                            }
        }
    }

    std::vector<std::array<vec2d, 3>> ints_sph(n_charges);
    for (size_t icharge = 0; icharge < n_charges; icharge++)
    {
        ints_sph[icharge][0] = trafo2Spherical(la, lb, ints_cart[icharge][0]);
        ints_sph[icharge][1] = trafo2Spherical(la, lb, ints_cart[icharge][1]);
        ints_sph[icharge][2] = trafo2Spherical(la, lb, ints_cart[icharge][2]);
    }

    size_t ofs_norm_a = sp_data.offsets_norms_[2 * ipair + 0];
    size_t ofs_norm_b = sp_data.offsets_norms_[2 * ipair + 1];
    for (size_t icharge = 0; icharge < n_charges; icharge++)
        for (size_t mu = 0; mu < ints_sph[icharge][0].dim<0>(); mu++)
            for (size_t nu = 0; nu < ints_sph[icharge][0].dim<1>(); nu++)
            {
                double norm_a = sp_data.norms_[ofs_norm_a + mu];
                double norm_b = sp_data.norms_[ofs_norm_b + nu];

                ints_sph[icharge][0](mu, nu) *= norm_a * norm_b;
                ints_sph[icharge][1](mu, nu) *= norm_a * norm_b;
                ints_sph[icharge][2](mu, nu) *= norm_a * norm_b;
            }

    return ints_sph;
}

std::vector<lible::vec2d>
lints::potentialAtExternalChargesKernel(const size_t ipair,
                                        const std::vector<std::array<double, 4>> &charges,
                                        const BoysGrid &boys_grid, const ShellPairData &sp_data)
{
    auto [la, lb] = sp_data.getLPair();
    int lab = la + lb;

    size_t ofs_prim = sp_data.offsets_primitives_[ipair];
    const double *exps = &sp_data.exps_[ofs_prim];
    const double *coeffs = &sp_data.coeffs_[ofs_prim];
    const double *xyz_a = &sp_data.coords_[6 * ipair + 0];
    const double *xyz_b = &sp_data.coords_[6 * ipair + 3];

    const auto &cart_exps_a = cart_exps[la];
    const auto &cart_exps_b = cart_exps[lb];

    const size_t n_ext_points = charges.size();

    std::vector<vec2d> ints_cart(n_ext_points, vec2d(Fill(0), numCartesians(la), numCartesians(lb)));
    for (size_t iab = 0; iab < sp_data.nrs_ppairs_[ipair]; iab++)
    {
        double a = exps[iab * 2];
        double b = exps[iab * 2 + 1];
        double da = coeffs[iab * 2];
        double db = coeffs[iab * 2 + 1];

        double p = a + b;
        double dadb = da * db;
        double fac = 2 * (M_PI / p) * dadb;

        auto [Ex, Ey, Ez] = ecoeffsPrimitivePair(a, b, la, lb, xyz_a, xyz_b);

        std::array<double, 3> xyz_p{
            (a * xyz_a[0] + b * xyz_b[0]) / p,
            (a * xyz_a[1] + b * xyz_b[1]) / p,
            (a * xyz_a[2] + b * xyz_b[2]) / p
        };

        std::vector<vec3d> rints_sum(n_ext_points, vec3d(Fill(0), lab + 1));

        for (size_t icharge = 0; icharge < charges.size(); icharge++)
        {
            auto [xc, yc, zc, charge] = charges[icharge];

            std::array<double, 3> xyz_pc{xyz_p[0] - xc, xyz_p[1] - yc, xyz_p[2] - zc};

            double xx{xyz_pc[0]}, xy{xyz_pc[1]}, xz{xyz_pc[2]};
            double xyz_pc_dot = xx * xx + xy * xy + xz * xz;
            double x = p * xyz_pc_dot;

            std::vector<double> fnx = calcBoysF(lab, x, boys_grid);

            vec3d rints = calcRInts3D(lab, p, &xyz_pc[0], &fnx[0]);

            for (int t = 0; t <= lab; t++)
                for (int u = 0; u <= lab; u++)
                    for (int v = 0; v <= lab; v++)
                        rints_sum[icharge](t, u, v) += charge * rints(t, u, v);

            for (const auto &[i, j, k, mu] : cart_exps_a)
                for (const auto &[i_, j_, k_, nu] : cart_exps_b)
                    for (int t = 0; t <= i + i_; t++)
                        for (int u = 0; u <= j + j_; u++)
                            for (int v = 0; v <= k + k_; v++)
                                ints_cart[icharge](mu, nu) += (-1) * fac * // -1 = charge of electron
                                        Ex(i, i_, t) * Ey(j, j_, u) * Ez(k, k_, v) *
                                        rints_sum[icharge](t, u, v);
        }
    }

    std::vector<vec2d> ints_sph(n_ext_points, vec2d(Fill(0), numSphericals(la), numSphericals(lb)));

    size_t ofs_norm_a = sp_data.offsets_norms_[2 * ipair + 0];
    size_t ofs_norm_b = sp_data.offsets_norms_[2 * ipair + 1];

    for (size_t icharge = 0; icharge < charges.size(); icharge++)
    {
        ints_sph[icharge] = trafo2Spherical(la, lb, ints_cart[icharge]);

        for (size_t mu = 0; mu < ints_sph[icharge].dim<0>(); mu++)
            for (size_t nu = 0; nu < ints_sph[icharge].dim<1>(); nu++)
            {
                double norm_a = sp_data.norms_[ofs_norm_a + mu];
                double norm_b = sp_data.norms_[ofs_norm_b + nu];
                ints_sph[icharge](mu, nu) *= norm_a * norm_b;
            }
    }

    return ints_sph;
}

std::vector<lible::vec2d>
lints::potentialAtExternalChargesErfKernel(const size_t ipair,
                                           const std::vector<std::array<double, 4>> &charges,
                                           const std::vector<double> &omegas,
                                           const BoysGrid &boys_grid, const ShellPairData &sp_data)
{
    if (charges.size() != omegas.size())
        throw std::runtime_error("potentialAtExternalChargesErfKernel(): number of omega values "
            "must match number of charges.");

    auto [la, lb] = sp_data.getLPair();
    int lab = la + lb;

    size_t ofs_prim = sp_data.offsets_primitives_[ipair];
    const double *exps = &sp_data.exps_[ofs_prim];
    const double *coeffs = &sp_data.coeffs_[ofs_prim];
    const double *xyz_a = &sp_data.coords_[6 * ipair + 0];
    const double *xyz_b = &sp_data.coords_[6 * ipair + 3];

    const auto &cart_exps_a = cart_exps[la];
    const auto &cart_exps_b = cart_exps[lb];

    size_t n_ext_points = charges.size();

    std::vector<vec2d> ints_cart(n_ext_points, vec2d(Fill(0), numCartesians(la), numCartesians(lb)));

    for (size_t iab = 0; iab < sp_data.nrs_ppairs_[ipair]; iab++)
    {
        double a = exps[iab * 2];
        double b = exps[iab * 2 + 1];
        double da = coeffs[iab * 2];
        double db = coeffs[iab * 2 + 1];

        double p = a + b;
        double dadb = da * db;
        double fac = 2 * (M_PI / p) * dadb;

        auto [Ex, Ey, Ez] = ecoeffsPrimitivePair(a, b, la, lb, xyz_a, xyz_b);

        std::array<double, 3> xyz_p{
            (a * xyz_a[0] + b * xyz_b[0]) / p,
            (a * xyz_a[1] + b * xyz_b[1]) / p,
            (a * xyz_a[2] + b * xyz_b[2]) / p
        };

        std::vector<vec3d> rints_sum(n_ext_points, vec3d(Fill(0), lab + 1));

        for (size_t icharge = 0; icharge < charges.size(); icharge++)
        {
            auto [xc, yc, zc, charge] = charges[icharge];

            std::array<double, 3> xyz_pc{xyz_p[0] - xc, xyz_p[1] - yc, xyz_p[2] - zc};

            double xx{xyz_pc[0]}, xy{xyz_pc[1]}, xz{xyz_pc[2]};
            double xyz_pc_dot = xx * xx + xy * xy + xz * xz;

            double omega_factor = omegas[icharge] /
                                  std::pow(omegas[icharge] * omegas[icharge] + p, 0.5);

            double omega_squared = omegas[icharge] * omegas[icharge];

            // second part of (52) in https://doi.org/10.1039/B605188J
            double x = p * xyz_pc_dot * omega_squared / (omega_squared + p);

            std::vector<double> fnx = calcBoysF(lab, x, boys_grid);

            vec3d rints = calcRInts3DErf(lab, p, omegas[icharge], &xyz_pc[0], &fnx[0]);

            for (int t = 0; t <= lab; t++)
                for (int u = 0; u <= lab; u++)
                    for (int v = 0; v <= lab; v++)
                        rints_sum[icharge](t, u, v) += charge * omega_factor * rints(t, u, v); // erf-related factor

            for (const auto &[i, j, k, mu] : cart_exps_a)
                for (const auto &[i_, j_, k_, nu] : cart_exps_b)
                    for (int t = 0; t <= i + i_; t++)
                        for (int u = 0; u <= j + j_; u++)
                            for (int v = 0; v <= k + k_; v++)
                                ints_cart[icharge](mu, nu) += (-1) * fac * // -1 = charge of electron
                                        Ex(i, i_, t) * Ey(j, j_, u) * Ez(k, k_, v) *
                                        rints_sum[icharge](t, u, v);
        }
    }

    std::vector<vec2d> ints_sph(n_ext_points, vec2d(Fill(0), numSphericals(la), numSphericals(lb)));

    size_t ofs_norm_a = sp_data.offsets_norms_[2 * ipair + 0];
    size_t ofs_norm_b = sp_data.offsets_norms_[2 * ipair + 1];

    for (size_t icharge = 0; icharge < charges.size(); icharge++)
    {
        ints_sph[icharge] = trafo2Spherical(la, lb, ints_cart[icharge]);

        for (size_t mu = 0; mu < ints_sph[icharge].dim<0>(); mu++)
            for (size_t nu = 0; nu < ints_sph[icharge].dim<1>(); nu++)
            {
                double norm_a = sp_data.norms_[ofs_norm_a + mu];
                double norm_b = sp_data.norms_[ofs_norm_b + nu];
                ints_sph[icharge](mu, nu) *= norm_a * norm_b;
            }
    }

    return ints_sph;
}

lible::vec2d lints::externalCharges(const std::vector<std::array<double, 4>> &charges,
                                    const Structure &structure)
{
    int l_max = structure.getMaxL();
    size_t dim_ao = structure.getDimAO();

    vec2d ints(Fill(0), dim_ao, dim_ao);
    for (int la = l_max; la >= 0; la--)
        for (int lb = la; lb >= 0; lb--)
        {
            ShellPairData sp_data(true, la, lb, structure.getShellsL(la), structure.getShellsL(lb));

            BoysGrid boys_grid(la + lb);

#pragma omp parallel for
            for (size_t ipair = 0; ipair < sp_data.n_pairs_; ipair++)
            {
                vec2d ints_ipair = externalChargesKernel(ipair, charges, boys_grid, sp_data);

                size_t ofs_a = sp_data.offsets_sph_[2 * ipair + 0];
                size_t ofs_b = sp_data.offsets_sph_[2 * ipair + 1];
                for (int a = 0; a < 2 * la + 1; a++)
                    for (int b = 0; b < 2 * lb + 1; b++)
                    {
                        size_t mu = ofs_a + a;
                        size_t nu = ofs_b + b;
                        ints(mu, nu) = ints_ipair(a, b);
                        ints(nu, mu) = ints_ipair(a, b);
                    }
            }
        }

    return ints;
}

lible::vec2d lints::externalChargesErf(const std::vector<std::array<double, 4>> &charges,
                                       const std::vector<double> &omegas,
                                       const Structure &structure)
{
    if (charges.size() != omegas.size())
        throw std::runtime_error("Number of charges and omegas must match.");

    int l_max = structure.getMaxL();
    size_t dim_ao = structure.getDimAO();

    vec2d ints(Fill(0), dim_ao, dim_ao);
    for (int la = l_max; la >= 0; la--)
        for (int lb = la; lb >= 0; lb--)
        {
            ShellPairData sp_data(true, la, lb, structure.getShellsL(la), structure.getShellsL(lb));

            int lab = la + lb;
            BoysGrid boys_grid(lab);

#pragma omp parallel for
            for (size_t ipair = 0; ipair < sp_data.n_pairs_; ipair++)
            {
                vec2d ints_ipair = externalChargesErfKernel(ipair, charges, omegas, boys_grid, sp_data);

                int n_sph_a = 2 * la + 1;
                int n_sph_b = 2 * lb + 1;
                size_t ofs_a = sp_data.offsets_sph_[2 * ipair + 0];
                size_t ofs_b = sp_data.offsets_sph_[2 * ipair + 1];
                for (int a = 0; a < n_sph_a; a++)
                    for (int b = 0; b < n_sph_b; b++)
                    {
                        size_t mu = ofs_a + a;
                        size_t nu = ofs_b + b;
                        ints(mu, nu) = ints_ipair(a, b);
                        ints(nu, mu) = ints_ipair(a, b);
                    }
            }
        }

    return ints;
}

std::array<lible::vec2d, 3> lints::spinOrbitCoupling1El(const Structure &structure)
{
    int l_max = structure.getMaxL();
    size_t dim_ao = structure.getDimAO();

    std::vector<std::array<double, 4>> charges = structure.getZs();

    std::array<vec2d, 3> ints;
    for (int i = 0; i < 3; i++)
        ints[i] = vec2d(Fill(0), dim_ao, dim_ao);

    for (int la = 0; la <= l_max; la++)
        for (int lb = 0; lb <= l_max; lb++)
        {
            ShellPairData sp_data(false, la, lb, structure.getShellsL(la), structure.getShellsL(lb));

            int lab = la + lb;
            BoysGrid boys_grid(lab + 1);

#pragma omp parallel for
            for (size_t ipair = 0; ipair < sp_data.n_pairs_; ipair++)
            {
                std::array<vec2d, 3> ints_batch =
                        spinOrbitCoupling1ElKernel(ipair, charges, boys_grid, sp_data);

                int n_sph_a = 2 * la + 1;
                int n_sph_b = 2 * lb + 1;
                size_t ofs_a = sp_data.offsets_sph_[2 * ipair];
                size_t ofs_b = sp_data.offsets_sph_[2 * ipair + 1];
                for (int a = 0; a < n_sph_a; a++)
                    for (int b = 0; b < n_sph_b; b++)
                    {
                        size_t mu = ofs_a + a;
                        size_t nu = ofs_b + b;
                        ints[0](mu, nu) = ints_batch[0](a, b);
                        ints[1](mu, nu) = ints_batch[1](a, b);
                        ints[2](mu, nu) = ints_batch[2](a, b);
                    }
            }
        }

    return ints;
}

std::array<lible::vec2d, 3>
lints::spinOrbitCoupling1ElKernel(const size_t ipair,
                                  const std::vector<std::array<double, 4>> &charges,
                                  const BoysGrid &boys_grid,
                                  const ShellPairData &sp_data)
{
    auto [la, lb] = sp_data.getLPair();
    int lab = la + lb;

    size_t ofs_prim = sp_data.offsets_primitives_[ipair];
    const double *exps = &sp_data.exps_[ofs_prim];
    const double *coeffs = &sp_data.coeffs_[ofs_prim];
    const double *xyz_a = &sp_data.coords_[6 * ipair + 0];
    const double *xyz_b = &sp_data.coords_[6 * ipair + 3];

    const auto &cart_exps_a = cart_exps[la];
    const auto &cart_exps_b = cart_exps[lb];

    std::array<vec2d, 3> ints_cart;
    for (int ideriv = 0; ideriv < 3; ideriv++)
        ints_cart[ideriv] = vec2d(Fill(0), numCartesians(la), numCartesians(lb));

    for (size_t iab = 0; iab < sp_data.nrs_ppairs_[ipair]; iab++)
    {
        double a = exps[iab * 2];
        double b = exps[iab * 2 + 1];
        double da = coeffs[iab * 2];
        double db = coeffs[iab * 2 + 1];

        auto [Ex, Ey, Ez] = ecoeffsPrimitivePair(a, b, la, lb, xyz_a, xyz_b);

        auto [E1x, E1y, E1z] = ecoeffsPrimitivePair_n1(a, b, la, lb, xyz_a, xyz_b,
                                                       {Ex, Ey, Ez});

        // R integrals
        double p = a + b;
        std::array<double, 3> xyz_p{
            (a * xyz_a[0] + b * xyz_b[0]) / p,
            (a * xyz_a[1] + b * xyz_b[1]) / p,
            (a * xyz_a[2] + b * xyz_b[2]) / p
        };

        vec3d rints_sum(Fill(0), lab + 2);
        for (auto [xc, yc, zc, charge] : charges)
        {
            std::array<double, 3> xyz_pc{xyz_p[0] - xc, xyz_p[1] - yc, xyz_p[2] - zc};

            double xx{xyz_pc[0]}, xy{xyz_pc[1]}, xz{xyz_pc[2]};
            double xyz_pc_dot = xx * xx + xy * xy + xz * xz;
            double x = p * xyz_pc_dot;

            std::vector<double> fnx = calcBoysF(lab + 1, x, boys_grid);

            vec3d rints = calcRInts3D(lab + 1, p, &xyz_pc[0], &fnx[0]);

            for (int t = 0; t <= lab + 1; t++)
                for (int u = 0; u <= lab + 1; u++)
                    for (int v = 0; v <= lab + 1; v++)
                        rints_sum(t, u, v) += charge * rints(t, u, v);
        }

        // (\nabla_A x \nabla_B) integrals
        double dadb = da * db;
        double fac = (2 * M_PI / p) * dadb;
        for (const auto &[i, j, k, mu] : cart_exps_a)
            for (const auto &[i_, j_, k_, nu] : cart_exps_b)
                for (int t = 0; t <= i + i_; t++)
                    for (int u = 0; u <= j + j_; u++)
                        for (int v = 0; v <= k + k_; v++)
                        {
                            double r100 = rints_sum(t + 1, u, v);
                            double r010 = rints_sum(t, u + 1, v);
                            double r001 = rints_sum(t, u, v + 1);
                            double e100 = E1x(i, i_, t) * Ey(j, j_, u) * Ez(k, k_, v);
                            double e010 = Ex(i, i_, t) * E1y(j, j_, u) * Ez(k, k_, v);
                            double e001 = Ex(i, i_, t) * Ey(j, j_, u) * E1z(k, k_, v);

                            double pr_yz = r010 * e001;
                            double pr_zy = r001 * e010;
                            double pr_zx = r001 * e100;
                            double pr_xz = r100 * e001;
                            double pr_xy = r100 * e010;
                            double pr_yx = r010 * e100;

                            // -1 = charge of electron
                            ints_cart[0](mu, nu) += -1.0 * fac * (pr_zy - pr_yz);
                            ints_cart[1](mu, nu) += -1.0 * fac * (pr_xz - pr_zx);
                            ints_cart[2](mu, nu) += -1.0 * fac * (pr_yx - pr_xy);
                        }
    }

    std::array<vec2d, 3> ints_sph;
    for (int icoord = 0; icoord < 3; icoord++)
        ints_sph[icoord] = trafo2Spherical(la, lb, ints_cart[icoord]);

    size_t ofs_norm_a = sp_data.offsets_norms_[2 * ipair + 0];
    size_t ofs_norm_b = sp_data.offsets_norms_[2 * ipair + 1];
    for (int icoord = 0; icoord < 3; icoord++)
        for (size_t mu = 0; mu < ints_sph[icoord].dim<0>(); mu++)
            for (size_t nu = 0; nu < ints_sph[icoord].dim<1>(); nu++)
            {
                double norm_a = sp_data.norms_[ofs_norm_a + mu];
                double norm_b = sp_data.norms_[ofs_norm_b + nu];
                ints_sph[icoord](mu, nu) *= norm_a * norm_b;
            }

    return ints_sph;
}

lible::arr2d<lible::vec2d, 3, 3>
lints::pVpKernel(const size_t ipair, const std::vector<std::array<double, 4>> &charges,
                 const BoysGrid &boys_grid, const ShellPairData &sp_data)
{
    auto [la, lb] = sp_data.getLPair();
    int lab = la + lb;

    size_t ofs_prim = sp_data.offsets_primitives_[ipair];
    const double *exps = &sp_data.exps_[ofs_prim];
    const double *coeffs = &sp_data.coeffs_[ofs_prim];
    const double *xyz_a = &sp_data.coords_[6 * ipair + 0];
    const double *xyz_b = &sp_data.coords_[6 * ipair + 3];

    const auto &cart_exps_a = cart_exps[la];
    const auto &cart_exps_b = cart_exps[lb];

    arr2d<vec2d, 3, 3> ints_cart;
    for (int id = 0; id < 3; id++)
        for (int jd = 0; jd < 3; jd++)
            ints_cart[id][jd] = vec2d(Fill(0), numCartesians(la), numCartesians(lb));

    for (size_t iab = 0; iab < sp_data.nrs_ppairs_[ipair]; iab++)
    {
        double a = exps[iab * 2];
        double b = exps[iab * 2 + 1];
        double da = coeffs[iab * 2];
        double db = coeffs[iab * 2 + 1];

        const auto [Ex, Ey, Ez] = ecoeffsPrimitivePair(a, b, la, lb, xyz_a, xyz_b);

        const auto [E1x, E1y, E1z] = ecoeffsPrimitivePair_n1(a, b, la, lb, xyz_a, xyz_b,
                                                             {Ex, Ey, Ez});

        const auto [E2x, E2y, E2z] = ecoeffsPrimitivePair_n2(a, b, la, lb, xyz_a, xyz_b,
                                                             {Ex, Ey, Ez}, {E1x, E1y, E1z});

        // R integrals
        const double p = a + b;
        std::array<double, 3> xyz_p{
            (a * xyz_a[0] + b * xyz_b[0]) / p,
            (a * xyz_a[1] + b * xyz_b[1]) / p,
            (a * xyz_a[2] + b * xyz_b[2]) / p
        };

        vec3d rints_sum(Fill(0), lab + 3);
        for (const auto &[xc, yc, zc, charge] : charges)
        {
            std::array<double, 3> xyz_pc{xyz_p[0] - xc, xyz_p[1] - yc, xyz_p[2] - zc};

            const double xx{xyz_pc[0]}, xy{xyz_pc[1]}, xz{xyz_pc[2]};
            const double xyz_pc_dot = xx * xx + xy * xy + xz * xz;
            const double x = p * xyz_pc_dot;

            std::vector<double> fnx = calcBoysF(lab + 2, x, boys_grid);

            vec3d rints = calcRInts3D(lab + 2, p, &xyz_pc[0], &fnx[0]);

            for (int t = 0; t <= lab + 2; t++)
                for (int u = 0; u <= lab + 2; u++)
                    for (int v = 0; v <= lab + 2; v++)
                        rints_sum(t, u, v) += charge * rints(t, u, v);
        }

        //
        double ab = a * b;
        double p2 = p * p;
        double dadb = da * db;
        double fac = (2 * M_PI / p) * dadb;
        for (const auto &[i, j, k, mu] : cart_exps_a)
            for (const auto &[i_, j_, k_, nu] : cart_exps_b)
                for (int t = 0; t <= i + i_; t++)
                    for (int u = 0; u <= j + j_; u++)
                        for (int v = 0; v <= k + k_; v++)
                        {
                            arr2d<double, 3, 3> pp{};
                            arr2d<double, 3, 3> pr{};
                            arr2d<double, 3, 3> rr{};

                            // PP
                            double r200 = rints_sum(t + 2, u, v);
                            double r110 = rints_sum(t + 1, u + 1, v);
                            double r101 = rints_sum(t + 1, u, v + 1);
                            double r020 = rints_sum(t, u + 2, v);
                            double r011 = rints_sum(t, u + 1, v + 1);
                            double r002 = rints_sum(t, u, v + 2);
                            pp[0][0] = r200;
                            pp[0][1] = r110;
                            pp[0][2] = r101;
                            pp[1][0] = r110;
                            pp[1][1] = r020;
                            pp[1][2] = r011;
                            pp[2][0] = r101;
                            pp[2][1] = r011;
                            pp[2][2] = r002;

                            // PR
                            double r100 = rints_sum(t + 1, u, v);
                            double r010 = rints_sum(t, u + 1, v);
                            double r001 = rints_sum(t, u, v + 1);
                            double e100 = E1x(i, i_, t) * Ey(j, j_, u) * Ez(k, k_, v);
                            double e010 = Ex(i, i_, t) * E1y(j, j_, u) * Ez(k, k_, v);
                            double e001 = Ex(i, i_, t) * Ey(j, j_, u) * E1z(k, k_, v);
                            pr[0][0] = r100 * e100;
                            pr[0][1] = r100 * e010;
                            pr[0][2] = r100 * e001;
                            pr[1][0] = r010 * e100;
                            pr[1][1] = r010 * e010;
                            pr[1][2] = r010 * e001;
                            pr[2][0] = r001 * e100;
                            pr[2][1] = r001 * e010;
                            pr[2][2] = r001 * e001;

                            // RR
                            double e200 = E2x(i, i_, t) * Ey(j, j_, u) * Ez(k, k_, v);
                            double e110 = E1x(i, i_, t) * E1y(j, j_, u) * Ez(k, k_, v);
                            double e101 = E1x(i, i_, t) * Ey(j, j_, u) * E1z(k, k_, v);
                            double e020 = Ex(i, i_, t) * E2y(j, j_, u) * Ez(k, k_, v);
                            double e011 = Ex(i, i_, t) * E1y(j, j_, u) * E1z(k, k_, v);
                            double e002 = Ex(i, i_, t) * Ey(j, j_, u) * E2z(k, k_, v);
                            rr[0][0] = e200;
                            rr[0][1] = e110;
                            rr[0][2] = e101;
                            rr[1][0] = e110;
                            rr[1][1] = e020;
                            rr[1][2] = e011;
                            rr[2][0] = e101;
                            rr[2][1] = e011;
                            rr[2][2] = e002;

                            for (int id = 0; id < 3; id++)
                                for (int jd = 0; jd < 3; jd++)
                                {
                                    ints_cart[id][jd](mu, nu) +=
                                            fac * ((ab / p2) * pp[id][jd] -
                                                   (a / p) * pr[id][jd] +
                                                   (b / p) * pr[jd][id] - rr[id][jd]);
                                }
                        }
    }

    arr2d<vec2d, 3, 3> ints_sph;
    for (int id = 0; id < 3; id++)
        for (int jd = 0; jd < 3; jd++)
            ints_sph[id][jd] = trafo2Spherical(la, lb, ints_cart[id][jd]);

    size_t ofs_norm_a = sp_data.offsets_norms_[2 * ipair + 0];
    size_t ofs_norm_b = sp_data.offsets_norms_[2 * ipair + 1];
    for (int id = 0; id < 3; id++)
        for (int jd = 0; jd < 3; jd++)
            for (size_t mu = 0; mu < ints_sph[id][jd].dim<0>(); mu++)
                for (size_t nu = 0; nu < ints_sph[id][jd].dim<1>(); nu++)
                {
                    double norm_a = sp_data.norms_[ofs_norm_a + mu];
                    double norm_b = sp_data.norms_[ofs_norm_b + nu];
                    ints_sph[id][jd](mu, nu) *= norm_a * norm_b;
                }

    return ints_sph;
}

std::array<lible::vec2d, 3> lints::momentumKernel(const size_t ipair, const ShellPairData &sp_data)
{
    auto [la, lb] = sp_data.getLPair();

    size_t ofs_prim = sp_data.offsets_primitives_[ipair];
    const double *exps = &sp_data.exps_[ofs_prim];
    const double *coeffs = &sp_data.coeffs_[ofs_prim];
    const double *xyz_a = &sp_data.coords_[6 * ipair + 0];
    const double *xyz_b = &sp_data.coords_[6 * ipair + 3];

    const auto &cart_exps_a = cart_exps[la];
    const auto &cart_exps_b = cart_exps[lb];

    std::array<vec2d, 3> ints_cart;
    for (int ideriv = 0; ideriv < 3; ideriv++)
        ints_cart[ideriv] = vec2d(Fill(0), numCartesians(la), numCartesians(lb));

    for (size_t iab = 0; iab < sp_data.nrs_ppairs_[ipair]; iab++)
    {
        double a = exps[iab * 2];
        double b = exps[iab * 2 + 1];
        double da = coeffs[iab * 2];
        double db = coeffs[iab * 2 + 1];

        const double p = a + b;
        const double dadb = da * db;
        const double fac = dadb * std::pow(M_PI / p, 1.5);

        const auto [Ex, Ey, Ez] = ecoeffsPrimitivePair(a, b, la, lb + 1, xyz_a, xyz_b);

        for (const auto &[i, j, k, mu] : cart_exps_a)
            for (const auto &[i_, j_, k_, nu] : cart_exps_b)
            {
                double D1x = -2 * b * Ex(i, i_ + 1, 0);
                if (i_ > 0)
                    D1x += i_ * Ex(i, i_ - 1, 0);

                double D1y = -2 * b * Ey(j, j_ + 1, 0);
                if (j_ > 0)
                    D1y += j_ * Ey(j, j_ - 1, 0);

                double D1z = -2 * b * Ez(k, k_ + 1, 0);
                if (k_ > 0)
                    D1z += k_ * Ez(k, k_ - 1, 0);

                ints_cart[0](mu, nu) -= fac * D1x * Ey(j, j_, 0) * Ez(k, k_, 0);
                ints_cart[1](mu, nu) -= fac * Ex(i, i_, 0) * D1y * Ez(k, k_, 0);
                ints_cart[2](mu, nu) -= fac * Ex(i, i_, 0) * Ey(j, j_, 0) * D1z;
            }
    }

    std::array<vec2d, 3> ints_sph;
    for (int ideriv = 0; ideriv < 3; ideriv++)
        ints_sph[ideriv] = trafo2Spherical(la, lb, ints_cart[ideriv]);

    size_t ofs_norm_a = sp_data.offsets_norms_[2 * ipair + 0];
    size_t ofs_norm_b = sp_data.offsets_norms_[2 * ipair + 1];
    for (int ideriv = 0; ideriv < 3; ideriv++)
        for (size_t mu = 0; mu < ints_sph[ideriv].dim<0>(); mu++)
            for (size_t nu = 0; nu < ints_sph[ideriv].dim<1>(); nu++)
                ints_sph[ideriv](mu, nu) *=
                        sp_data.norms_[ofs_norm_a + mu] * sp_data.norms_[ofs_norm_b + nu];

    return ints_sph;
}

std::array<lible::vec2d, 3>
lints::angularMomentumKernel(const size_t ipair, const ShellPairData &sp_data)
{
    auto [la, lb] = sp_data.getLPair();

    size_t ofs_prim = sp_data.offsets_primitives_[ipair];
    const double *exps = &sp_data.exps_[ofs_prim];
    const double *coeffs = &sp_data.coeffs_[ofs_prim];

    const double *xyz_a = &sp_data.coords_[6 * ipair + 0];
    const double *xyz_b = &sp_data.coords_[6 * ipair + 3];

    const auto &cart_exps_a = cart_exps[la];
    const auto &cart_exps_b = cart_exps[lb];

    std::array<vec2d, 3> ints_cart;
    for (int ideriv = 0; ideriv < 3; ideriv++)
        ints_cart[ideriv] = vec2d(Fill(0), numCartesians(la), numCartesians(lb));

    for (size_t iab = 0; iab < sp_data.nrs_ppairs_[ipair]; iab++)
    {
        double a = exps[iab * 2];
        double b = exps[iab * 2 + 1];
        double da = coeffs[iab * 2];
        double db = coeffs[iab * 2 + 1];

        const double p = a + b;
        const double dadb = da * db;
        const double fac = dadb * std::pow(M_PI / p, 1.5);

        const auto [Ex, Ey, Ez] = ecoeffsPrimitivePair(a, b, la, lb + 1, xyz_a, xyz_b);

        std::array<double, 3> xyz_p{
            (a * xyz_a[0] + b * xyz_b[0]) / p,
            (a * xyz_a[1] + b * xyz_b[1]) / p,
            (a * xyz_a[2] + b * xyz_b[2]) / p
        };

        for (const auto &[i, j, k, mu] : cart_exps_a)
            for (const auto &[i_, j_, k_, nu] : cart_exps_b)
            {
                double D1x = -2 * b * Ex(i, i_ + 1, 0);
                if (i_ > 0)
                    D1x += i_ * Ex(i, i_ - 1, 0);

                double D1y = -2 * b * Ey(j, j_ + 1, 0);
                if (j_ > 0)
                    D1y += j_ * Ey(j, j_ - 1, 0);

                double D1z = -2 * b * Ez(k, k_ + 1, 0);
                if (k_ > 0)
                    D1z += k_ * Ez(k, k_ - 1, 0);

                double S1x = xyz_p[0] * Ex(i, i_, 0);
                if (i + i_ > 0)
                    S1x += Ex(i, i_, 1);

                double S1y = xyz_p[1] * Ey(j, j_, 0);
                if (j + j_ > 0)
                    S1y += Ey(j, j_, 1);

                double S1z = xyz_p[2] * Ez(k, k_, 0);
                if (k + k_ > 0)
                    S1z += Ez(k, k_, 1);

                ints_cart[0](mu, nu) -= fac * Ex(i, i_, 0) * (S1y * D1z - D1y * S1z);
                ints_cart[1](mu, nu) -= fac * Ex(j, j_, 0) * (D1x * S1z - S1x * D1z);
                ints_cart[2](mu, nu) -= fac * Ez(k, k_, 0) * (S1x * D1y - D1x * S1y);
            }
    }

    std::array<vec2d, 3> ints_sph;
    for (int ideriv = 0; ideriv < 3; ideriv++)
        ints_sph[ideriv] = trafo2Spherical(la, lb, ints_cart[ideriv]);

    size_t ofs_norm_a = sp_data.offsets_norms_[2 * ipair + 0];
    size_t ofs_norm_b = sp_data.offsets_norms_[2 * ipair + 1];
    for (int ideriv = 0; ideriv < 3; ideriv++)
        for (size_t mu = 0; mu < ints_sph[ideriv].dim<0>(); mu++)
            for (size_t nu = 0; nu < ints_sph[ideriv].dim<1>(); nu++)
                ints_sph[ideriv](mu, nu) *=
                        sp_data.norms_[ofs_norm_a + mu] * sp_data.norms_[ofs_norm_b + nu];

    return ints_sph;
}

lible::vec2d lints::nuclearAttraction(const Structure &structure)
{
    int l_max = structure.getMaxL();
    size_t dim_ao = structure.getDimAO();

    std::vector<std::array<double, 4>> charges(structure.getNAtoms());
    for (size_t iatom = 0; iatom < structure.getNAtoms(); iatom++)
    {
        std::array<double, 3> coords = structure.getCoordsAtom(iatom);
        std::array<double, 3> xyz_c{coords[0], coords[1], coords[2]};

        const double Z = structure.getZ(iatom);
        charges[iatom] = {xyz_c[0], xyz_c[1], xyz_c[2], Z};
    }

    vec2d ints(Fill(0), dim_ao, dim_ao);
    for (int la = l_max; la >= 0; la--)
        for (int lb = la; lb >= 0; lb--)
        {
            ShellPairData sp_data(true, la, lb, structure.getShellsL(la), structure.getShellsL(lb));

            int lab = la + lb;
            BoysGrid boys_grid(lab);

#pragma omp parallel for
            for (size_t ipair = 0; ipair < sp_data.n_pairs_; ipair++)
            {
                vec2d ints_ipair = externalChargesKernel(ipair, charges, boys_grid, sp_data);

                size_t ofs_a = sp_data.offsets_sph_[2 * ipair + 0];
                size_t ofs_b = sp_data.offsets_sph_[2 * ipair + 1];
                for (size_t mu = 0; mu < ints_ipair.dim<0>(); mu++)
                    for (size_t nu = 0; nu < ints_ipair.dim<1>(); nu++)
                    {
                        ints(ofs_a + mu, ofs_b + nu) = ints_ipair(mu, nu);
                        ints(ofs_b + nu, ofs_a + mu) = ints_ipair(mu, nu);
                    }
            }
        }

    return ints;
}

lible::vec2d lints::nuclearAttractionErf(const Structure &structure,
                                         const std::vector<double> &omegas)
{
    if (omegas.size() != structure.getNAtoms())
        throw std::runtime_error("Number of omega values must match number of atoms.");

    int l_max = structure.getMaxL();
    size_t dim_ao = structure.getDimAO();

    std::vector<std::array<double, 4>> charges(structure.getNAtoms());
    for (size_t iatom = 0; iatom < structure.getNAtoms(); iatom++)
    {
        const std::array<double, 3> coords = structure.getCoordsAtom(iatom);
        const std::array<double, 3> xyz_c{coords[0], coords[1], coords[2]};

        const double Z = structure.getZ(iatom);
        charges[iatom] = {xyz_c[0], xyz_c[1], xyz_c[2], Z};
    }

    vec2d ints(Fill(0), dim_ao, dim_ao);
    for (int la = l_max; la >= 0; la--)
        for (int lb = la; lb >= 0; lb--)
        {
            ShellPairData sp_data(true, la, lb, structure.getShellsL(la), structure.getShellsL(lb));

            int lab = la + lb;
            BoysGrid boys_grid(lab);

#pragma omp parallel for
            for (size_t ipair = 0; ipair < sp_data.n_pairs_; ipair++)
            {
                vec2d ints_ipair = externalChargesErfKernel(ipair, charges, omegas, boys_grid, sp_data);

                size_t ofs_a = sp_data.offsets_sph_[2 * ipair + 0];
                size_t ofs_b = sp_data.offsets_sph_[2 * ipair + 1];
                for (size_t mu = 0; mu < ints_ipair.dim<0>(); mu++)
                    for (size_t nu = 0; nu < ints_ipair.dim<1>(); nu++)
                    {
                        ints(ofs_a + mu, ofs_b + nu) = ints_ipair(mu, nu);
                        ints(ofs_b + nu, ofs_a + mu) = ints_ipair(mu, nu);
                    }
            }
        }

    return ints;
}


lible::arr2d<lible::vec2d, 3, 3> lints::pVpIntegrals(const Structure &structure)
{
    int l_max = structure.getMaxL();
    size_t dim_ao = structure.getDimAO();

    std::vector<std::array<double, 4>> charges(structure.getNAtoms());
    for (size_t iatom = 0; iatom < structure.getNAtoms(); iatom++)
    {
        std::array<double, 3> coords = structure.getCoordsAtom(iatom);
        std::array<double, 3> xyz_c{coords[0], coords[1], coords[2]};

        const double Z = structure.getZ(iatom);
        charges[iatom] = {xyz_c[0], xyz_c[1], xyz_c[2], Z};
    }

    arr2d<vec2d, 3, 3> ints;
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            ints[i][j] = vec2d(Fill(0), dim_ao, dim_ao);

    for (int la = l_max; la >= 0; la--)
        for (int lb = la; lb >= 0; lb--)
        {
            ShellPairData sp_data(true, la, lb, structure.getShellsL(la), structure.getShellsL(lb));

            int lab = la + lb;
            BoysGrid boys_grid(lab + 2);

#pragma omp parallel for
            for (size_t ipair = 0; ipair < sp_data.n_pairs_; ipair++)
            {
                arr2d<vec2d, 3, 3> ints_ipair = pVpKernel(ipair, charges, boys_grid, sp_data);

                size_t ofs_a = sp_data.offsets_sph_[2 * ipair + 0];
                size_t ofs_b = sp_data.offsets_sph_[2 * ipair + 1];

                for (int i = 0; i < 3; i++)
                    for (int j = 0; j < 3; j++)
                        for (size_t mu = 0; mu < ints_ipair[i][j].dim<0>(); mu++)
                            for (size_t nu = 0; nu < ints_ipair[i][j].dim<1>(); nu++)
                                ints[i][j](ofs_a + mu, ofs_b + nu) = ints_ipair[i][j](mu, nu);
            }
        }

    return ints;
}
