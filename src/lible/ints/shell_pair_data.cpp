#include <lible/ints/shell_pair_data.hpp>
#include <lible/ints/utils.hpp>

namespace lints = lible::ints;

using std::pair, std::vector;

lints::ShellData::ShellData(const int l, const std::vector<Shell> &shells)
    : l(l)
{
    n_shells = shells.size();
    n_primitives = 0;
    for (const Shell &shell : shells)
        n_primitives += shell.coeffs.size();

    const int n_coords = 3 * n_shells;
    const int n_sph = numSphericals(l);
    const int n_herm = numHermites(l);
    const int n_sph_ecoeffs = n_sph * n_herm;
    const int n_norms = n_sph * n_shells;

    atomic_idxs.resize(n_shells);
    cdepths.resize(n_shells);
    coffsets.resize(n_shells);
    offsets_ecoeffs.resize(n_shells);
    offsets_norms.resize(n_shells);
    offsets_sph.resize(n_shells);
    shell_idxs.resize(n_shells);

    coeffs.resize(n_primitives);
    exps.resize(n_primitives);
    coords.resize(n_coords);
    norms.resize(n_norms);

    int ofs_coeffs{0}, ofs_ecoeffs{0}, ofs_norms{0};
    for (size_t ishell = 0; ishell < shells.size(); ishell++)
    {
        const Shell &shell = shells[ishell];
        if (shell.l != l)
            throw std::runtime_error("ShellData::ShellData(): given angular momentum and shell "
                                     "angular momentum don't match");

        int cdepth = shell.coeffs.size();
        for (int i = 0; i < cdepth; i++)
        {
            coeffs[ofs_coeffs + i] = shell.coeffs[i] * shell.norms_prim[i]; // include primitive norms
            exps[ofs_coeffs + i] = shell.exps[i];
        }

        coords[3 * ishell + 0] = shell.xyz_coords[0];
        coords[3 * ishell + 1] = shell.xyz_coords[1];
        coords[3 * ishell + 2] = shell.xyz_coords[2];

        for (int i = 0; i < n_sph; i++)
            norms[ofs_norms + i] = shell.norms[i];

        atomic_idxs[ishell] = shell.atom_idx;
        cdepths[ishell] = cdepth;
        coffsets[ishell] = ofs_coeffs;
        offsets_ecoeffs[ishell] = ofs_ecoeffs;
        offsets_norms[ishell] = ofs_norms;
        offsets_sph[ishell] = shell.ofs_sph;
        shell_idxs[ishell] = shell.idx;

        ofs_coeffs += cdepth;
        ofs_ecoeffs += cdepth * n_sph_ecoeffs;
        ofs_norms += n_sph;
    }
}

lints::ShellPairData::ShellPairData(const bool use_symm, const int la, const int lb,
                                    const vector<Shell> &shells_a, const vector<Shell> &shells_b)
    : uses_symm(use_symm), la(la), lb(lb)
{
    const int n_shells_a = shells_a.size();
    const int n_shells_b = shells_b.size();

    if (uses_symm == true)
    {
        if (la == lb)
            n_pairs = n_shells_a * (n_shells_a + 1) / 2;
        else
            n_pairs = n_shells_a * n_shells_b;
    }
    else
        n_pairs = n_shells_a * n_shells_b;

    vector<pair<int, int>> shell_pair_idxs(n_pairs);
    for (int ishell = 0, ipair = 0; ishell < n_shells_a; ishell++)
    {
        int bound_shell_b;
        if ((uses_symm == true) && (la == lb))
            bound_shell_b = ishell + 1;
        else
            bound_shell_b = n_shells_b;

        for (int jshell = 0; jshell < bound_shell_b; jshell++, ipair++)
            shell_pair_idxs[ipair] = {ishell, jshell};
    }

    n_prim_pairs = 0;
    int n_coeffs = 0;
    for (const auto &[ishell, jshell] : shell_pair_idxs)
    {
        n_prim_pairs += shells_a[ishell].coeffs.size() * shells_b[jshell].coeffs.size();
        n_coeffs += shells_a[ishell].coeffs.size() + shells_b[jshell].coeffs.size();
    }

    const int n_shells_ab = n_pairs * 2;
    const int n_sph_a = numSphericals(la);
    const int n_sph_b = numSphericals(lb);
    const int n_herm = numHermites(la + lb);
    const int n_sph_ecoeffs = n_sph_a * n_sph_b * n_herm;
    const int n_norms = (n_sph_a + n_sph_b) * n_pairs;

    coeffs.resize(n_coeffs);
    coords.resize(3 * n_shells_ab);
    exps.resize(n_coeffs);
    norms.resize(n_norms);

    atomic_idxs.resize(n_shells_ab);
    cdepths.resize(n_shells_ab);
    coffsets.resize(n_shells_ab);
    offsets_cart.resize(n_shells_ab);
    offsets_ecoeffs.resize(n_pairs);
    offsets_ecoeffs_deriv1.resize(n_pairs);
    offsets_ecoeffs_deriv2.resize(n_pairs);
    offsets_norms.resize(n_shells_ab);
    offsets_sph.resize(n_shells_ab);
    shell_idxs.resize(n_shells_ab);

    int ofs_coeffs{0}, ofs_norms{0}, ofs_ecoeffs{0}, ofs_ecoeffs_d1{0}, ofs_ecoeffs_d2{0};
    for (int ipair = 0; ipair < n_pairs; ipair++)
    {
        const auto &[ishell, jshell] = shell_pair_idxs[ipair];

        const Shell &shell_a = shells_a[ishell];
        const Shell &shell_b = shells_b[jshell];

        if (shell_a.l != la)
            throw std::runtime_error("ShellPairData::ShellPairData(): given la and shell_a angular "
                                     "momentum don't match");
        if (shell_b.l != lb)
            throw std::runtime_error("ShellPairData::ShellPairData(): given lb and shell_b angular "
                                     "momentum don't match");

        const int cdepth_a = shell_a.coeffs.size();
        const int cdepth_b = shell_b.coeffs.size();
        const int cdepth_ab = cdepth_a * cdepth_b;

        // printf("ipair = %d, cdepth_a = %d, cdepth_b = %d\n", ipair, cdepth_a, cdepth_b);

        atomic_idxs[2 * ipair + 0] = shell_a.atom_idx;
        atomic_idxs[2 * ipair + 1] = shell_b.atom_idx;
        cdepths[2 * ipair + 0] = cdepth_a;
        cdepths[2 * ipair + 1] = cdepth_b;

        coords[6 * ipair + 0] = shell_a.xyz_coords[0];
        coords[6 * ipair + 1] = shell_a.xyz_coords[1];
        coords[6 * ipair + 2] = shell_a.xyz_coords[2];
        coords[6 * ipair + 3] = shell_b.xyz_coords[0];
        coords[6 * ipair + 4] = shell_b.xyz_coords[1];
        coords[6 * ipair + 5] = shell_b.xyz_coords[2];

        coffsets[2 * ipair + 0] = ofs_coeffs;
        for (int a = 0; a < cdepth_a; a++, ofs_coeffs++)
        {
            coeffs[ofs_coeffs] = shell_a.coeffs[a] * shell_a.norms_prim[a]; // include primitive norms
            exps[ofs_coeffs] = shell_a.exps[a];
        }

        coffsets[2 * ipair + 1] = ofs_coeffs;
        for (int b = 0; b < cdepth_b; b++, ofs_coeffs++)
        {
            coeffs[ofs_coeffs] = shell_b.coeffs[b] * shell_b.norms_prim[b]; // include primitive norms
            exps[ofs_coeffs] = shell_b.exps[b];
        }

        offsets_cart[2 * ipair + 0] = shell_a.ofs_cart;
        offsets_cart[2 * ipair + 1] = shell_b.ofs_cart;
        offsets_sph[2 * ipair + 0] = shell_a.ofs_sph;
        offsets_sph[2 * ipair + 1] = shell_b.ofs_sph;
        offsets_ecoeffs[ipair] = ofs_ecoeffs;
        offsets_ecoeffs_deriv1[ipair] = ofs_ecoeffs_d1;
        offsets_ecoeffs_deriv2[ipair] = ofs_ecoeffs_d2;

        offsets_norms[2 * ipair + 0] = ofs_norms;
        for (int a = 0; a < n_sph_a; a++, ofs_norms++)
            norms[ofs_norms] = shell_a.norms[a];

        offsets_norms[2 * ipair + 1] = ofs_norms;
        for (int b = 0; b < n_sph_b; b++, ofs_norms++)
            norms[ofs_norms] = shell_b.norms[b];

        offsets_sph[2 * ipair + 0] = shell_a.ofs_sph;
        offsets_sph[2 * ipair + 1] = shell_b.ofs_sph;
        shell_idxs[2 * ipair + 0] = shell_a.idx;
        shell_idxs[2 * ipair + 1] = shell_b.idx;

        ofs_ecoeffs += cdepth_ab * n_sph_ecoeffs;
        ofs_ecoeffs_d1 += 3 * cdepth_ab * n_sph_ecoeffs;
        ofs_ecoeffs_d2 += 6 * cdepth_ab * n_sph_ecoeffs;
    }
}

vector<lints::ShellData> lints::shellDataAux(const Structure &structure)
{
    if (structure.getUseRI() == false)
        throw std::runtime_error("shellDataAux(): RI approximation is not enabled!");

    int l_max_aux = structure.getMaxLAux();

    vector<ShellData> sh_data;
    for (int l = 0; l <= l_max_aux; l++)
        sh_data.emplace_back(ShellData(l, structure.getShellsLAux(l)));

    return sh_data;
}

vector<lints::ShellPairData> lints::shellPairData(const bool use_symm, const Structure &structure)
{
    vector<pair<int, int>> l_pairs;
    if (use_symm == true)
        l_pairs = getLPairsSymm(structure.getMaxL());
    else
        l_pairs = getLPairsNoSymm(structure.getMaxL());

    vector<ShellPairData> sp_data;
    for (const auto &[la, lb] : l_pairs)
        sp_data.emplace_back(ShellPairData(use_symm, la, lb, structure.getShellsL(la),
                                           structure.getShellsL(lb)));

    return sp_data;
}