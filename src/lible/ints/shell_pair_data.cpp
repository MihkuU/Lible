#include <lible/ints/shell_pair_data.hpp>
#include <lible/ints/utils.hpp>

namespace LI = lible::ints;

using std::pair, std::vector;

LI::ShellPairData::ShellPairData()
{
}

LI::ShellPairData::ShellPairData(const int la, const int lb, const Structure &structure)
    : la(la), lb(lb), structure(&structure)
{
    int lab = la + lb;
    int n_ecoeffs_sph = (2 * la + 1) * (2 * lb + 1) * (lab + 1) * (lab + 2) * (lab + 3) / 6;

    if (la == lb)
    {
        auto shells_a = structure.getShellsL(la);

        size_t size_a = shells_a.size();

        n_pairs = size_a * (size_a + 1) / 2;

        offsets.reserve(n_pairs);
        offsets_cart.reserve(n_pairs);
        coords.reserve(n_pairs);
        ccoeffs.reserve(n_pairs);
        exps.reserve(n_pairs);

        n_prim_pairs = 0;
        for (size_t ishell = 0; ishell < size_a; ishell++)
            for (size_t jshell = 0; jshell <= ishell; jshell++)
            {
                offsets.push_back(std::make_pair(shells_a[ishell].pos,
                                                 shells_a[jshell].pos));

                offsets_cart.push_back(std::make_pair(shells_a[ishell].pos_cart,
                                                      shells_a[jshell].pos_cart));

                coords.push_back(std::make_pair(shells_a[ishell].xyz_coords,
                                                shells_a[jshell].xyz_coords));

                ccoeffs.push_back(std::make_pair(shells_a[ishell].coeffs,
                                                 shells_a[jshell].coeffs));

                exps.push_back(std::make_pair(shells_a[ishell].exps,
                                              shells_a[jshell].exps));

                norms.push_back(std::make_pair(shells_a[ishell].norms,
                                               shells_a[jshell].norms));

                for (size_t i = 0; i < shells_a[ishell].exps.size(); i++)
                    for (size_t j = 0; j < shells_a[jshell].exps.size(); j++)
                        n_prim_pairs++;
            }

        offsets_prims.resize(n_pairs);
        offsets_ecoeffs.resize(n_prim_pairs);

        size_t ipair = 0, iprim_pair = 0, offset = 0;
        for (size_t ishell = 0; ishell < size_a; ishell++)
            for (size_t jshell = 0; jshell <= ishell; jshell++)
            {
                offsets_prims[ipair] = iprim_pair;
                for (size_t i = 0; i < shells_a[ishell].exps.size(); i++)
                    for (size_t j = 0; j < shells_a[jshell].exps.size(); j++)
                    {
                        offsets_ecoeffs[iprim_pair] = offset;
                        offset += n_ecoeffs_sph;
                        iprim_pair++;
                    }
                ipair++;
            }
    }
    else
    {
        auto shells_a = structure.getShellsL(la);
        auto shells_b = structure.getShellsL(lb);

        size_t size_a = shells_a.size();
        size_t size_b = shells_b.size();

        n_pairs = size_a * size_b;

        offsets.reserve(n_pairs);
        offsets_cart.reserve(n_pairs);
        coords.reserve(n_pairs);
        ccoeffs.reserve(n_pairs);
        exps.reserve(n_pairs);

        n_prim_pairs = 0;
        for (size_t ishell = 0; ishell < size_a; ishell++)
            for (size_t jshell = 0; jshell < size_b; jshell++)
            {
                offsets.push_back(std::make_pair(shells_a[ishell].pos,
                                                 shells_b[jshell].pos));

                offsets_cart.push_back(std::make_pair(shells_a[ishell].pos_cart,
                                                      shells_b[jshell].pos_cart));

                coords.push_back(std::make_pair(shells_a[ishell].xyz_coords,
                                                shells_b[jshell].xyz_coords));

                ccoeffs.push_back(std::make_pair(shells_a[ishell].coeffs,
                                                 shells_b[jshell].coeffs));

                exps.push_back(std::make_pair(shells_a[ishell].exps,
                                              shells_b[jshell].exps));

                norms.push_back(std::make_pair(shells_a[ishell].norms,
                                               shells_b[jshell].norms));

                for (size_t i = 0; i < shells_a[ishell].exps.size(); i++)
                    for (size_t j = 0; j < shells_b[jshell].exps.size(); j++)
                        n_prim_pairs++;
            }

        offsets_prims.resize(n_pairs);
        offsets_ecoeffs.resize(n_prim_pairs);

        size_t ipair = 0, iprim_pair = 0, offset = 0;
        for (size_t ishell = 0; ishell < size_a; ishell++)
            for (size_t jshell = 0; jshell < size_b; jshell++)
            {
                offsets_prims[ipair] = iprim_pair;
                for (size_t i = 0; i < shells_a[ishell].exps.size(); i++)
                    for (size_t j = 0; j < shells_b[jshell].exps.size(); j++)
                    {
                        offsets_ecoeffs[iprim_pair] = offset;
                        offset += n_ecoeffs_sph;
                        iprim_pair++;
                    }
                ipair++;
            }
    }
}

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

    int n_coords = 6 * n_pairs;
    int n_norms = (dim_sph_a + dim_sph_b) * n_pairs;

    int n_coeffs{0};
    for (auto &[ishell, jshell] : shell_pair_idxs)
    {
        int cdepth_a = shells_a[ishell].coeffs.size();
        int cdepth_b = shells_b[jshell].coeffs.size();

        n_coeffs += cdepth_a + cdepth_b;
    }

    int n_prim_pairs{}; // TODO

    std::vector<double> coeffs(n_coeffs);
    std::vector<double> coords(n_coords);
    std::vector<double> exps(n_coeffs);
    std::vector<double> norms(n_norms);

    std::vector<int> cdepths(2 * n_pairs);
    std::vector<int> coffsets(2 * n_pairs);

    std::vector<int> offsets_cart(2 * n_pairs);
    std::vector<int> offsets_ecoeffs; // TODO
    std::vector<int> offsets_norms(2 * n_pairs);
    std::vector<int> offsets_primitives; // TODO
    std::vector<int> offsets_sph(2 * n_pairs);

    int pos_coords{0}, pos_coeffs{0}, pos_norms{0};
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
                                      offsets_cart, offsets_ecoeffs, offsets_norms,
                                      offsets_primitives, offsets_sph);

    return shell_pair_data;
}