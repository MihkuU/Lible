#include <lible/ints/shell_pair_data.hpp>
#include <lible/ints/utils.hpp>

namespace LI = lible::ints;

using std::pair, std::vector;

LI::ShellPairData_new LI::constructShellPairData(const int la, const int lb,
                                                 const Structure &structure)
{
    const auto &shells_a = structure.getShellsL(la);
    const auto &shells_b = structure.getShellsL(lb);
    int n_shells_a = shells_a.size();
    int n_shells_b = shells_b.size();

    int n_pairs;
    if (la == lb)
        n_pairs = n_shells_a * (n_shells_a + 1) / 2;
    else
        n_pairs = n_shells_a * n_shells_b;

    vector<pair<int, int>> shell_pair_idxs(n_pairs);
    if (la == lb)
        for (int ishell = 0, idx = 0; ishell < n_shells_a; ishell++)
            for (int jshell = 0; jshell <= ishell; jshell++, idx++)
                shell_pair_idxs[idx] = std::make_pair(ishell, jshell);
    else
        for (int ishell = 0, idx = 0; ishell < n_shells_a; ishell++)
            for (int jshell = 0; jshell < n_shells_b; jshell++, idx++)
                shell_pair_idxs[idx] = std::make_pair(ishell, jshell);

    int dim_sph_a = dimSphericals(la);
    int dim_sph_b = dimSphericals(lb);
    int dim_herm_gauss = dimHermiteGaussians(la + lb);
    int n_sph_ecoeffs = dim_sph_a * dim_sph_b * dim_herm_gauss;

    int n_coords = 6 * n_pairs;
    int n_norms = (dim_sph_a + dim_sph_b) * n_pairs;

    int n_coeffs{0};
    for (auto &[ishell, jshell] : shell_pair_idxs)
    {
        int cdepth_a = shells_a[ishell].coeffs.size();
        int cdepth_b = shells_b[jshell].coeffs.size();
        n_coeffs += cdepth_a + cdepth_b;
    }

    int n_prim_pairs{0};

    std::vector<double> coeffs(n_coeffs);
    std::vector<double> coords(n_coords);
    std::vector<double> exps(n_coeffs);
    std::vector<double> norms(n_norms);

    std::vector<int> cdepths(2 * n_pairs);
    std::vector<int> coffsets(2 * n_pairs);

    std::vector<int> offsets_cart(2 * n_pairs);
    std::vector<int> offsets_ecoeffs(n_pairs); 
    std::vector<int> offsets_norms(2 * n_pairs);
    std::vector<int> offsets_sph(2 * n_pairs);

    int pos_coords{0}, pos_coeffs{0}, pos_norms{0}, pos_ecoeffs{0};
    for (int ipair = 0; ipair < n_pairs; ipair++)
    {
        auto [ishell, jshell] = shell_pair_idxs[ipair];

        const auto &shell_a = shells_a[ishell];
        const auto &shell_b = shells_b[jshell];

        coords[pos_coords + 0] = shell_a.xyz_coords[0];
        coords[pos_coords + 1] = shell_a.xyz_coords[1];
        coords[pos_coords + 2] = shell_a.xyz_coords[2];
        coords[pos_coords + 3] = shell_b.xyz_coords[0];
        coords[pos_coords + 4] = shell_b.xyz_coords[1];
        coords[pos_coords + 5] = shell_b.xyz_coords[2];
        pos_coords += 6;

        int cdepth_a = shell_a.coeffs.size();
        int cdepth_b = shell_b.coeffs.size();

        offsets_ecoeffs[ipair] += pos_ecoeffs;
        pos_ecoeffs += cdepth_a * cdepth_b * n_sph_ecoeffs;
        n_prim_pairs += cdepth_a * cdepth_b;

        coffsets[2 * ipair + 0] = pos_coeffs;

        for (int ia = 0; ia < cdepth_a; ia++)
        {
            coeffs[pos_coeffs] = shell_a.coeffs[ia];
            exps[pos_coeffs] = shell_a.exps[ia];
            pos_coeffs++;
        }

        coffsets[2 * ipair + 1] = pos_coeffs;

        for (int ib = 0; ib < cdepth_b; ib++)
        {
            coeffs[pos_coeffs] = shell_b.coeffs[ib];
            exps[pos_coeffs] = shell_b.exps[ib];
            pos_coeffs++;
        }

        cdepths[2 * ipair + 0] = cdepth_a;
        cdepths[2 * ipair + 1] = cdepth_b;

        offsets_cart[2 * ipair + 0] = shell_a.pos_cart;
        offsets_cart[2 * ipair + 1] = shell_b.pos_cart;
        offsets_sph[2 * ipair + 0] = shell_a.pos;
        offsets_sph[2 * ipair + 1] = shell_b.pos;

        const auto &norms_a = shell_a.norms;
        const auto &norms_b = shell_b.norms;

        offsets_norms[2 * ipair + 0] = pos_norms;
        
        for (int mu = 0; mu < dim_sph_a; mu++)
        {
            norms[pos_norms] = norms_a[mu];
            pos_norms++;
        }            

        offsets_norms[2 * ipair + 1] = pos_norms;

        for (int mu = 0; mu < dim_sph_b; mu++)
        {
            norms[pos_norms] = norms_b[mu];
            pos_norms++;
        }
    }

    int n_atoms = structure.getNAtoms();
    vector<double> atomic_coords(3 * n_atoms);
    vector<int> atomic_nrs(n_atoms);
    for (int iatom = 0; iatom < n_atoms; iatom++)
    {
        auto [x, y, z] = structure.getCoordsAtom(iatom);
        atomic_coords[3 * iatom] = x;
        atomic_coords[3 * iatom + 1] = y;
        atomic_coords[3 * iatom + 2] = z;
        atomic_nrs[iatom] = structure.getZ(iatom);
    }

    ShellPairData_new shell_pair_data(la, lb, n_atoms, n_pairs, n_prim_pairs, atomic_coords,
                                      coeffs, coords, exps, norms, atomic_nrs, cdepths, coffsets,
                                      offsets_cart, offsets_ecoeffs, offsets_norms, offsets_sph);

    return shell_pair_data;
}