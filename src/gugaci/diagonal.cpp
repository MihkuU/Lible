#include <lible/guga_sci.h>
#include <lible/util.h>

#ifdef _USE_MPI_
// #include <boost/mpi.hpp>
// namespace mpi = boost::mpi;
#endif

using namespace lible;
using namespace lible::guga;
using namespace lible::guga::util;

using std::string;
using std::vector;

vector<double> GCI::Impl::calcCCXDiag(const int &p, const int &q, const CFG *cfg)
{
    vector<double> ccx(cfg->getNumCSFs(), 0);

    char dp, dq, dk;
    int b;
    arma::ivec b_vals(n_orbs, arma::fill::zeros);
    for (size_t icsf = 0; icsf < cfg->getNumCSFs(); icsf++)
    {
        string csf = cfg->getCSF(icsf);
        dp = csf[p];
        dq = csf[q];
        int b_val = 0;
        for (size_t i = 0; i < n_orbs; i++)
        {
            b_val += determineStepb(csf[i]);
            b_vals[i] = b_val;
        }

        double contrib;
        if (dp == '1' and dq == '1')
        {
            b = b_vals[p];
            contrib = A(b, 2, 0);
            for (int k = p + 1; k < q; k++)
            {
                b = b_vals[k];
                dk = csf[k];
                contrib *= f(dk, b);
            }
            b = b_vals[q];
            contrib *= A(b, -1, 1);
            ccx[icsf] = contrib;
        }
        else if (dp == '1' and dq == '2')
        {
            b = b_vals[p];
            contrib = -A(b, 2, 0);
            for (int k = p + 1; k < q; k++)
            {
                b = b_vals[k];
                dk = csf[k];
                contrib *= f(dk, b);
            }
            b = b_vals[q];
            contrib *= A(b, 3, 1);
            ccx[icsf] = contrib;
        }
        else if (dp == '2' and dq == '1')
        {
            b = b_vals[p];
            contrib = -A(b, 0, 2);
            for (int k = p + 1; k < q; k++)
            {
                b = b_vals[k];
                dk = csf[k];
                contrib *= f(dk, b);
            }
            b = b_vals[q];
            contrib *= A(b, -1, 1);
            ccx[icsf] = contrib;
        }
        else if (dp == '2' and dq == '2')
        {
            b = b_vals[p];
            contrib = A(b, 0, 2);
            for (int k = p + 1; k < q; k++)
            {
                b = b_vals[k];
                dk = csf[k];
                contrib *= f(dk, b);
            }
            b = b_vals[q];
            contrib *= A(b, 3, 1);
            ccx[icsf] = contrib;
        }
    }

    return ccx;
}

vector<double> GCI::Impl::calcDiag()
{
    arma::dvec diag(wave_function->getNumCSFs(), arma::fill::zeros);

#pragma omp parallel
    {

#ifdef _USE_MPI_
        // int rank_total = returnTotalRank(world);
        // int size_total = returnTotalSize(world);
        int rank_total = omp_get_thread_num(); // TMP
        int size_total = omp_get_num_threads(); // TMP
#else
        int rank_total = omp_get_thread_num();
        int size_total = omp_get_num_threads();
#endif

        arma::dvec diag_omp(wave_function->getNumCSFs(), arma::fill::zeros);
        for (size_t icfg = 0; icfg < wave_function->getNumCFGs(); icfg++)
        {
            if (icfg % size_total != rank_total)
                continue;

            CFG *conf_p = wave_function->getCFGPtr(icfg);

            size_t dim = wave_function->getDim(icfg);
            size_t pos = wave_function->getPos(icfg);
            string onv = conf_p->getONV();

            for (size_t p = 0; p < n_orbs; p++)
            {
                if (onv[p] == '0')
                    continue;

                int p_occ = onv[p] - '0';
                for (size_t mu = 0; mu < dim; mu++)
                    diag_omp(pos + mu) += p_occ * one_el_ints(p, p);
            }

            for (size_t p = 0; p < n_orbs; p++)
            {
                if (onv[p] != '2')
                    continue;

                for (size_t mu = 0; mu < dim; mu++)
                    diag_omp(pos + mu) += two_el_ints(p, p, p, p);
            }

            for (size_t p = 0; p < n_orbs; p++)
            {
                if (onv[p] == '0')
                    continue;

                int p_occ = onv[p] - '0';
                for (size_t q = p + 1; q < n_orbs; q++)
                {
                    if (onv[q] == '0')
                        continue;
                    int q_occ = onv[q] - '0';

                    for (size_t mu = 0; mu < dim; mu++)
                        diag_omp(pos + mu) += p_occ * q_occ * two_el_ints(p, p, q, q);
                }
            }

            for (size_t p = 0; p < n_orbs; p++)
            {
                if (onv[p] == '0')
                    continue;

                int p_occ = onv[p] - '0';
                for (size_t q = p + 1; q < n_orbs; q++)
                {
                    int q_occ = onv[q] - '0';

                    vector<double> ccx_diag = calcCCXDiag(p, q, conf_p);
                    for (size_t mu = 0; mu < dim; mu++)
                        diag_omp(pos + mu) -= 0.5 * two_el_ints(p, q, q, p) *
                                              (p_occ * q_occ + ccx_diag[mu]);
                }
            }
        }

#pragma omp critical
        {
            diag += diag_omp;
        }
    }

#ifdef _USE_MPI_
    // mpi::all_reduce(mpi::communicator(world, mpi::comm_duplicate),
    //                 diag, diag, std::plus<arma::dvec>());
#endif

    return arma::conv_to<vector<double>>::from(diag);
}

vector<double> GCI::Impl::calcDiag(const WaveFunction *wave_function)
{
    vector<double> diag(wave_function->getNumCSFs(), 0);

    for (size_t icfg = 0; icfg < wave_function->getNumCFGs(); icfg++)
    {
        const CFG *cfg = wave_function->getCFGPtr(icfg);

        size_t dim = wave_function->getDim(icfg);
        size_t pos = wave_function->getPos(icfg);
        string onv = cfg->getONV();

        for (size_t p = 0; p < n_orbs; p++)
        {
            if (onv[p] == '0')
                continue;

            int p_occ = onv[p] - '0';
            for (size_t mu = 0; mu < dim; mu++)
                diag[pos + mu] += p_occ * one_el_ints(p, p);
        }

        for (size_t p = 0; p < n_orbs; p++)
        {
            if (onv[p] != '2')
                continue;

            for (size_t mu = 0; mu < dim; mu++)
                diag[pos + mu] += two_el_ints(p, p, p, p);
        }

        for (size_t p = 0; p < n_orbs; p++)
        {
            if (onv[p] == '0')
                continue;

            int p_occ = onv[p] - '0';
            for (size_t q = p + 1; q < n_orbs; q++)
            {
                if (onv[q] == '0')
                    continue;
                int q_occ = onv[q] - '0';

                for (size_t mu = 0; mu < dim; mu++)
                    diag[pos + mu] += p_occ * q_occ * two_el_ints(p, p, q, q);
            }
        }

        for (size_t p = 0; p < n_orbs; p++)
        {
            if (onv[p] == '0')
                continue;

            int p_occ = onv[p] - '0';
            for (size_t q = p + 1; q < n_orbs; q++)
            {
                int q_occ = onv[q] - '0';

                vector<double> ccx_diag = calcCCXDiag(p, q, cfg);
                for (size_t mu = 0; mu < dim; mu++)
                    diag[pos + mu] -= 0.5 * two_el_ints(p, q, q, p) * (p_occ * q_occ + ccx_diag[mu]);
            }
        }
    }

    return diag;
}