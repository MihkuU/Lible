#include <lible/ints/ints.hpp>
#include <lible/ints/shell_pair_data.hpp>
#include <lible/ints/utils.hpp>

namespace lints = lible::ints;

lints::ShellData::ShellData(const int l, const std::vector<Shell> &shells)
    : l_(l)
{
    if (omp_in_parallel() == true)
        throw std::runtime_error("ShellData(): cannot be initialized inside a parallel region");

    n_shells_ = shells.size();
    n_primitives_ = 0;
    for (const Shell &shell : shells)
        n_primitives_ += shell.coeffs_.size();

    int n_sph = numSphericals(l);
    int n_herm = numHermites(l);
    int n_sph_ecoeffs = n_sph * n_herm;
    size_t n_norms = n_sph * n_shells_;
    size_t n_coords = 3 * n_shells_;

    atomic_idxs_.resize(n_shells_);
    cdepths_.resize(n_shells_);
    coffsets_.resize(n_shells_);
    offsets_ecoeffs_.resize(n_shells_);
    offsets_norms_.resize(n_shells_);
    offsets_sph_.resize(n_shells_);
    shell_idxs_.resize(n_shells_);

    coeffs_.resize(n_primitives_);
    exps_.resize(n_primitives_);
    coords_.resize(n_coords);
    norms_.resize(n_norms);

    size_t ofs_coeffs{0}, ofs_ecoeffs{0}, ofs_norms{0};
    for (size_t ishell = 0; ishell < shells.size(); ishell++)
    {
        const Shell &shell = shells[ishell];
        if (shell.l_ != l)
            throw std::runtime_error("ShellData(): given angular momentum and shell angular "
                "momentum don't match");

        size_t cdepth = shell.coeffs_.size();
        for (size_t i = 0; i < cdepth; i++)
        {
            coeffs_[ofs_coeffs + i] = shell.coeffs_[i] * shell.norms_prim_[i]; // include primitive norms
            exps_[ofs_coeffs + i] = shell.exps_[i];
        }

        coords_[3 * ishell + 0] = shell.xyz_coords_[0];
        coords_[3 * ishell + 1] = shell.xyz_coords_[1];
        coords_[3 * ishell + 2] = shell.xyz_coords_[2];

        for (int i = 0; i < n_sph; i++)
            norms_[ofs_norms + i] = shell.norms_[i];

        atomic_idxs_[ishell] = shell.idx_atom_;
        cdepths_[ishell] = cdepth;
        coffsets_[ishell] = ofs_coeffs;
        offsets_ecoeffs_[ishell] = ofs_ecoeffs;
        offsets_norms_[ishell] = ofs_norms;
        offsets_sph_[ishell] = shell.ofs_sph_;
        shell_idxs_[ishell] = shell.idx_;

        ofs_coeffs += cdepth;
        ofs_ecoeffs += cdepth * n_sph_ecoeffs;
        ofs_norms += n_sph;
    }
}

lints::ShellPairData::ShellPairData(const bool use_symm, const int la, const int lb,
                                    const std::vector<Shell> &shells_a,
                                    const std::vector<Shell> &shells_b,
                                    const double primitives_thrs)
    : uses_symm_(use_symm), primitives_thrs_(primitives_thrs), la_(la), lb_(lb)
{
    if (omp_in_parallel() == true)
        throw std::runtime_error("ShellPairData(): cannot be initialized inside a parallel region");

    countPairs(shells_a, shells_b, n_pairs_, n_pairs_total_, n_ppairs_, n_ppairs_total_);

    // Write the shell pair data based on how many significant shell and primitive pairs are there.
    int n_sph_a = numSphericals(la_);
    int n_sph_b = numSphericals(lb_);
    int n_herm = numHermites(la + lb);
    int n_sph_ecoeffs = n_sph_a * n_sph_b * n_herm;
    size_t n_norms = (n_sph_a + n_sph_b) * n_pairs_;
    size_t n_shells_ab = 2 * n_pairs_;
    size_t n_primitives = 2 * n_ppairs_;

    coeffs_.resize(n_primitives);
    exps_.resize(n_primitives);
    coords_.resize(6 * n_pairs_);
    norms_.resize(n_norms);

    nrs_ppairs_.resize(n_pairs_);
    offsets_primitives_.resize(n_pairs_);
    offsets_sph_.resize(n_shells_ab);
    offsets_cart_.resize(n_shells_ab);
    offsets_norms_.resize(n_shells_ab);
    atomic_idxs_.resize(n_shells_ab);
    shell_idxs_.resize(n_shells_ab);
    offsets_ecoeffs_.resize(n_pairs_);
    offsets_ecoeffs_deriv1_.resize(n_pairs_);
    offsets_ecoeffs_deriv2_.resize(n_pairs_);

    size_t ofs_primitives{0}, ofs_norms{0}, ofs_ecoeffs{0}, ofs_ecoeffs_d1{0}, ofs_ecoeffs_d2{0};
    for (size_t ishell = 0, ipair = 0; ishell < shells_a.size(); ishell++)
    {
        size_t bound_shell_b;
        if (uses_symm_ == true && la_ == lb_)
            bound_shell_b = ishell + 1;
        else
            bound_shell_b = shells_b.size();

        for (size_t jshell = 0; jshell < bound_shell_b; jshell++)
        {
            const auto &shell_a = shells_a[ishell];
            const auto &shell_b = shells_b[jshell];

            const auto &xyz_a = shell_a.xyz_coords_;
            const auto &xyz_b = shell_b.xyz_coords_;
            double RAB2 = std::pow(xyz_a[0] - xyz_b[0], 2) +
                          std::pow(xyz_a[1] - xyz_b[1], 2) +
                          std::pow(xyz_a[2] - xyz_b[2], 2);

            size_t cdepth_a = shells_a[ishell].coeffs_.size();
            size_t cdepth_b = shells_b[jshell].coeffs_.size();

            // Find how many survivors are there.
            size_t n_ppairs_survived = 0;
            for (size_t ia = 0, iab = 0; ia < cdepth_a; ia++)
                for (size_t ib = 0; ib < cdepth_b; ib++)
                {
                    double a = shell_a.exps_[ia];
                    double b = shell_b.exps_[ib];
                    double da = shell_a.norms_prim_[ia] * shell_a.coeffs_[ia];
                    double db = shell_b.norms_prim_[ib] * shell_b.coeffs_[ib];
                    double dadb = da * db;

                    double mu = a * b / (a + b);
                    double Kab = std::exp(-mu * RAB2);
                    double screen_val = std::fabs(dadb * Kab);
                    if (screen_val >= primitives_thrs)
                    {
                        exps_[ofs_primitives + iab * 2 + 0] = a;
                        exps_[ofs_primitives + iab * 2 + 1] = b;
                        coeffs_[ofs_primitives + iab * 2 + 0] = da;
                        coeffs_[ofs_primitives + iab * 2 + 1] = db;

                        n_ppairs_survived++;
                        iab++;
                    }
                }

            // Write data.
            if (n_ppairs_survived > 0)
            {
                coords_[6 * ipair + 0] = xyz_a[0];
                coords_[6 * ipair + 1] = xyz_a[1];
                coords_[6 * ipair + 2] = xyz_a[2];
                coords_[6 * ipair + 3] = xyz_b[0];
                coords_[6 * ipair + 4] = xyz_b[1];
                coords_[6 * ipair + 5] = xyz_b[2];

                nrs_ppairs_[ipair] = n_ppairs_survived;
                offsets_primitives_[ipair] = ofs_primitives;
                offsets_cart_[2 * ipair + 0] = shell_a.ofs_cart_;
                offsets_cart_[2 * ipair + 1] = shell_b.ofs_cart_;
                offsets_sph_[2 * ipair + 0] = shell_a.ofs_sph_;
                offsets_sph_[2 * ipair + 1] = shell_b.ofs_sph_;

                offsets_norms_[2 * ipair + 0] = ofs_norms;
                for (int mu = 0; mu < n_sph_a; mu++)
                    norms_[ofs_norms + mu] = shell_a.norms_[mu];
                ofs_norms += n_sph_a;

                offsets_norms_[2 * ipair + 1] = ofs_norms;
                for (int nu = 0; nu < n_sph_b; nu++)
                    norms_[ofs_norms + nu] = shell_b.norms_[nu];
                ofs_norms += n_sph_b;

                offsets_ecoeffs_[ipair] = ofs_ecoeffs;
                offsets_ecoeffs_deriv1_[ipair] = ofs_ecoeffs_d1;
                offsets_ecoeffs_deriv2_[ipair] = ofs_ecoeffs_d2;

                atomic_idxs_[2 * ipair + 0] = shell_a.idx_atom_;
                atomic_idxs_[2 * ipair + 1] = shell_b.idx_atom_;
                shell_idxs_[2 * ipair + 0] = shell_a.idx_;
                shell_idxs_[2 * ipair + 1] = shell_b.idx_;

                ipair++;
                ofs_primitives += 2 * n_ppairs_survived;
                ofs_ecoeffs += n_ppairs_survived * n_sph_ecoeffs;
                ofs_ecoeffs_d1 += 3 * n_ppairs_survived * n_sph_ecoeffs;
                ofs_ecoeffs_d2 += 6 * n_ppairs_survived * n_sph_ecoeffs;
            }
        }
    }
}

void lints::ShellPairData::countPairs(const std::vector<Shell> &shells_a,
                                      const std::vector<Shell> &shells_b, size_t &n_pairs,
                                      size_t &n_pairs_total, size_t &n_ppairs,
                                      size_t &n_ppairs_total) const
{
    size_t n_shells_a = shells_a.size();
    size_t n_shells_b = shells_b.size();

    for (size_t ishell = 0; ishell < n_shells_a; ishell++)
    {
        size_t bound_shell_b;
        if (uses_symm_ == true && la_ == lb_)
            bound_shell_b = ishell + 1;
        else
            bound_shell_b = n_shells_b;

        for (size_t jshell = 0; jshell < bound_shell_b; jshell++)
        {
            const auto &shell_a = shells_a[ishell];
            const auto &shell_b = shells_b[jshell];

            if (shell_a.l_ != la_ || shell_b.l_ != lb_)
                throw std::runtime_error("ShellPairData(): given angular momentum and shell "
                    "angular momentum don't match");

            const auto &xyz_a = shell_a.xyz_coords_;
            const auto &xyz_b = shell_b.xyz_coords_;
            double RAB2 = std::pow(xyz_a[0] - xyz_b[0], 2) +
                          std::pow(xyz_a[1] - xyz_b[1], 2) +
                          std::pow(xyz_a[2] - xyz_b[2], 2);

            size_t cdepth_a = shells_a[ishell].coeffs_.size();
            size_t cdepth_b = shells_b[jshell].coeffs_.size();

            size_t n_ppairs_survived = 0;
            for (size_t ia = 0; ia < cdepth_a; ia++)
                for (size_t ib = 0; ib < cdepth_b; ib++)
                {
                    double a = shell_a.exps_[ia];
                    double b = shell_b.exps_[ib];
                    double da = shell_a.norms_prim_[ia] * shell_a.coeffs_[ia];
                    double db = shell_b.norms_prim_[ib] * shell_b.coeffs_[ib];
                    double dadb = da * db;

                    double mu = a * b / (a + b);
                    double Kab = std::exp(-mu * RAB2);
                    double screen_val = std::fabs(dadb * Kab);
                    if (screen_val >= primitives_thrs_)
                        n_ppairs_survived++;
                }

            n_pairs_total++;
            n_ppairs_total += cdepth_a * cdepth_b;

            if (n_ppairs_survived > 0)
            {
                n_pairs++;
                n_ppairs += n_ppairs_survived;
            }
        }
    }
}

std::vector<lible::ints::ShellData> lible::ints::shellData(const std::vector<Shell> &shells)
{
    std::map<int, std::vector<Shell>> shells_map;
    for (const Shell &shell : shells)
        shells_map[shell.l_].push_back(shell);

    std::vector<ShellData> sh_data;
    for (const auto &[l, shells_l] : shells_map)
        sh_data.emplace_back(l, shells_l);

    return sh_data;
}

std::vector<lints::ShellData> lints::shellDataAux(const Structure &structure)
{
    if (structure.getUseRI() == false)
        throw std::runtime_error("shellDataAux(): RI approximation is not enabled!");

    int l_max_aux = structure.getMaxLAux();

    std::vector<ShellData> sh_data;
    for (int l = 0; l <= l_max_aux; l++)
        sh_data.emplace_back(l, structure.getShellsLAux(l));

    return sh_data;
}

std::vector<lints::ShellPairData>
lints::shellPairData(const bool use_symm, const Structure &structure)
{
    std::vector<std::pair<int, int>> l_pairs;
    if (use_symm == true)
        l_pairs = getLPairsSymm(structure.getMaxL());
    else
        l_pairs = getLPairsNoSymm(structure.getMaxL());

    std::vector<ShellPairData> sp_data;
    for (const auto &[la, lb] : l_pairs)
        sp_data.emplace_back(use_symm, la, lb, structure.getShellsL(la),
                             structure.getShellsL(lb));

    return sp_data;
}

std::vector<lints::ShellPairData> lints::shellPairData(const std::vector<Shell> &shells_a,
                                                       const std::vector<Shell> &shells_b)
{
    std::map<int, std::vector<Shell>> shells_a_map;
    for (const Shell &shell : shells_a)
        shells_a_map[shell.l_].push_back(shell);

    std::map<int, std::vector<Shell>> shells_b_map;
    for (const Shell &shell : shells_b)
        shells_b_map[shell.l_].push_back(shell);

    std::vector<ShellPairData> sp_data;
    for (const auto &[la, shells_la] : shells_a_map)
        for (const auto &[lb, shells_lb] : shells_b_map)
            sp_data.emplace_back(false, la, lb, shells_la, shells_lb);

    return sp_data;
}