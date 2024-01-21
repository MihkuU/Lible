#include <lible/mcmurchie_davidson.h>

#include <any> // TMP
#include <chrono>
#include <functional> //TMP
#include <map>     // TMP
#include <utility> // TMP

#include <armadillo>

namespace LIMD = lible::ints::MD;

using namespace lible;

void LIMD::coeffs(const double &a, const double &b, const double &PA, const double &PB,
                  const double one_o_2p, const int &la, const int &lb, vec3d &E)
{
    for (int i = 1; i <= la; i++)
    {
        E(i, 0, 0) = PA * E(i - 1, 0, 0) + E(i - 1, 0, 1);

        for (int t = 1; t < i; t++)
            E(i, 0, t) = one_o_2p * E(i - 1, 0, t - 1) +
                         PA * E(i - 1, 0, t) +
                         (t + 1) * E(i - 1, 0, t + 1);

        E(i, 0, i) = one_o_2p * E(i - 1, 0, i - 1) + PA * E(i - 1, 0, i);
    }

    for (int j = 1; j <= lb; j++)
        for (int i = 0; i <= la; i++)
        {
            E(i, j, 0) = PB * E(i, j - 1, 0) + E(i, j - 1, 1);

            for (int t = 1; t < i + j; t++)
                E(i, j, t) = one_o_2p * E(i, j - 1, t - 1) +
                             PB * E(i, j - 1, t) +
                             (t + 1) * E(i, j - 1, t + 1);

            E(i, j, i + j) = one_o_2p * E(i, j - 1, i + j - 1) + PB * E(i, j - 1, i + j);
        }
}

void LIMD::coeffs(const double &a, const double &b, const int &la, const int &lb,
                  const std::array<double, 3> &A, const std::array<double, 3> &B,
                  vec3d &E_x, vec3d &E_y, vec3d &E_z)
{
    E_x.set(0);
    E_y.set(0);
    E_z.set(0);

    E_x(0, 0, 0) = 1;
    E_y(0, 0, 0) = 1;
    E_z(0, 0, 0) = 1;

    double p = a + b;
    double one_o_2p = 1.0 / (2 * p);

    arma::vec::fixed<3> RA{A[0], A[1], A[2]};
    arma::vec::fixed<3> RB{B[0], B[1], B[2]};

    arma::vec::fixed<3> P = (a * RA + b * RB) / p;
    arma::vec::fixed<3> PA = P - RA;
    arma::vec::fixed<3> PB = P - RB;

    coeffs(a, b, PA[0], PB[0], one_o_2p, la, lb, E_x);
    coeffs(a, b, PA[1], PB[1], one_o_2p, la, lb, E_y);
    coeffs(a, b, PA[2], PB[2], one_o_2p, la, lb, E_z);
}

/*
 * TODO: Test the following:
 *   1) Fastor
 *   2) XTensor
 *   3) Templates
 *   4) Recursive templates
 *   5) One array for all three
 *   6) Nested-pointer 3D array
 *   ...
 */

void coeffs4D(const double one_o_2p, const int &la, const int &lb,
              const arma::vec::fixed<3> &PA, const arma::vec::fixed<3> &PB,
              vec4d &E)
{
    for (int x = 0; x < 3; x++)
    {
        for (int i = 1; i <= la; i++)
        {
            E(x, i, 0, 0) = PA[x] * E(x, i - 1, 0, 0) + E(x, i - 1, 0, 1);

            for (int t = 1; t < i; t++)
                E(x, i, 0, t) = one_o_2p * E(x, i - 1, 0, t - 1) +
                                PA[x] * E(x, i - 1, 0, t) +
                                (t + 1) * E(x, i - 1, 0, t + 1);

            E(x, i, 0, i) = one_o_2p * E(x, i - 1, 0, i - 1) + PA[x] * E(x, i - 1, 0, i);
        }

        for (int j = 1; j <= lb; j++)
            for (int i = 0; i <= la; i++)
            {
                E(x, i, j, 0) = PB[x] * E(x, i, j - 1, 0) + E(x, i, j - 1, 1);

                for (int t = 1; t < i + j; t++)
                    E(x, i, j, t) = one_o_2p * E(x, i, j - 1, t - 1) +
                                    PB[x] * E(x, i, j - 1, t) +
                                    (t + 1) * E(x, i, j - 1, t + 1);

                E(x, i, j, i + j) = one_o_2p * E(x, i, j - 1, i + j - 1) + PB[x] * E(x, i, j - 1, i + j);
            }
    }
}

void testCoeffs4D(const double &a, const double &b, const int &la, const int &lb,
                  const std::array<double, 3> &A, const std::array<double, 3> &B,
                  vec4d &E_xyz)
{
    E_xyz.set(0);
    E_xyz(0, 0, 0, 0) = 1;
    E_xyz(1, 0, 0, 0) = 1;
    E_xyz(2, 0, 0, 0) = 1;

    double p = a + b;
    double one_o_2p = 1.0 / (2 * p);

    arma::vec::fixed<3> RA{A[0], A[1], A[2]};
    arma::vec::fixed<3> RB{B[0], B[1], B[2]};

    arma::vec::fixed<3> P = (a * RA + b * RB) / p;
    arma::vec::fixed<3> PA = P - RA;
    arma::vec::fixed<3> PB = P - RB;

    coeffs4D(one_o_2p, la, lb, PA, PB, E_xyz);
}

void coeffs3DRolled(const double one_o_2p, const int &la, const int &lb,
                    const arma::vec::fixed<3> &PA, const arma::vec::fixed<3> &PB,
                    vec3d &E_x, vec3d &E_y, vec3d &E_z)
{
    for (int i = 1; i <= la; i++)
    {
        E_x(i, 0, 0) = PA[0] * E_x(i - 1, 0, 0) + E_x(i - 1, 0, 1);
        E_y(i, 0, 0) = PA[1] * E_y(i - 1, 0, 0) + E_y(i - 1, 0, 1);
        E_z(i, 0, 0) = PA[2] * E_z(i - 1, 0, 0) + E_z(i - 1, 0, 1);

        for (int t = 1; t < i; t++)
        {
            E_x(i, 0, t) = one_o_2p * E_x(i - 1, 0, t - 1) +
                           PA[0] * E_x(i - 1, 0, t) +
                           (t + 1) * E_x(i - 1, 0, t + 1);

            E_y(i, 0, t) = one_o_2p * E_y(i - 1, 0, t - 1) +
                           PA[1] * E_y(i - 1, 0, t) +
                           (t + 1) * E_y(i - 1, 0, t + 1);

            E_z(i, 0, t) = one_o_2p * E_z(i - 1, 0, t - 1) +
                           PA[2] * E_z(i - 1, 0, t) +
                           (t + 1) * E_z(i - 1, 0, t + 1);
        }

        E_x(i, 0, i) = one_o_2p * E_x(i - 1, 0, i - 1) + PA[0] * E_x(i - 1, 0, i);
        E_y(i, 0, i) = one_o_2p * E_y(i - 1, 0, i - 1) + PA[1] * E_y(i - 1, 0, i);
        E_z(i, 0, i) = one_o_2p * E_z(i - 1, 0, i - 1) + PA[2] * E_z(i - 1, 0, i);
    }

    for (int j = 1; j <= lb; j++)
        for (int i = 0; i <= la; i++)
        {
            E_x(i, j, 0) = PB[0] * E_x(i, j - 1, 0) + E_x(i, j - 1, 1);
            E_y(i, j, 0) = PB[1] * E_y(i, j - 1, 0) + E_y(i, j - 1, 1);
            E_z(i, j, 0) = PB[2] * E_z(i, j - 1, 0) + E_z(i, j - 1, 1);

            for (int t = 1; t < i + j; t++)
            {
                E_x(i, j, t) = one_o_2p * E_x(i, j - 1, t - 1) +
                               PB[0] * E_x(i, j - 1, t) +
                               (t + 1) * E_x(i, j - 1, t + 1);

                E_y(i, j, t) = one_o_2p * E_y(i, j - 1, t - 1) +
                               PB[1] * E_y(i, j - 1, t) +
                               (t + 1) * E_y(i, j - 1, t + 1);

                E_z(i, j, t) = one_o_2p * E_z(i, j - 1, t - 1) +
                               PB[2] * E_z(i, j - 1, t) +
                               (t + 1) * E_z(i, j - 1, t + 1);
            }

            E_x(i, j, i + j) = one_o_2p * E_x(i, j - 1, i + j - 1) + PB[0] * E_x(i, j - 1, i + j);
            E_y(i, j, i + j) = one_o_2p * E_y(i, j - 1, i + j - 1) + PB[1] * E_y(i, j - 1, i + j);
            E_z(i, j, i + j) = one_o_2p * E_z(i, j - 1, i + j - 1) + PB[2] * E_z(i, j - 1, i + j);
        }
}

void testCoeffs3DRolled(const double &a, const double &b, const int &la, const int &lb,
                        const std::array<double, 3> &A, const std::array<double, 3> &B,
                        vec3d &E_x, vec3d &E_y, vec3d &E_z)
{
    E_x.set(0);
    E_y.set(0);
    E_z.set(0);
    E_x(0, 0, 0) = 1;
    E_y(0, 0, 0) = 1;
    E_z(0, 0, 0) = 1;

    double p = a + b;
    double one_o_2p = 1.0 / (2 * p);

    arma::vec::fixed<3> RA{A[0], A[1], A[2]};
    arma::vec::fixed<3> RB{B[0], B[1], B[2]};

    arma::vec::fixed<3> P = (a * RA + b * RB) / p;
    arma::vec::fixed<3> PA = P - RA;
    arma::vec::fixed<3> PB = P - RB;

    coeffs3DRolled(one_o_2p, la, lb, PA, PB, E_x, E_y, E_z);
}

template <int la, int lb>
void coeffs3DRolledTemplated(const double one_o_2p, const arma::vec::fixed<3> &PA,
                             const arma::vec::fixed<3> &PB,
                             vec3d &E_x, vec3d &E_y, vec3d &E_z)
{
    for (int i = 1; i <= la; i++)
    {
        E_x(i, 0, 0) = PA[0] * E_x(i - 1, 0, 0) + E_x(i - 1, 0, 1);
        E_y(i, 0, 0) = PA[1] * E_y(i - 1, 0, 0) + E_y(i - 1, 0, 1);
        E_z(i, 0, 0) = PA[2] * E_z(i - 1, 0, 0) + E_z(i - 1, 0, 1);

        for (int t = 1; t < i; t++)
        {
            E_x(i, 0, t) = one_o_2p * E_x(i - 1, 0, t - 1) +
                           PA[0] * E_x(i - 1, 0, t) +
                           double(t + 1) * E_x(i - 1, 0, t + 1);

            E_y(i, 0, t) = one_o_2p * E_y(i - 1, 0, t - 1) +
                           PA[1] * E_y(i - 1, 0, t) +
                           double(t + 1) * E_y(i - 1, 0, t + 1);

            E_z(i, 0, t) = one_o_2p * E_z(i - 1, 0, t - 1) +
                           PA[2] * E_z(i - 1, 0, t) +
                           double(t + 1) * E_z(i - 1, 0, t + 1);
        }

        E_x(i, 0, i) = one_o_2p * E_x(i - 1, 0, i - 1) + PA[0] * E_x(i - 1, 0, i);
        E_y(i, 0, i) = one_o_2p * E_y(i - 1, 0, i - 1) + PA[1] * E_y(i - 1, 0, i);
        E_z(i, 0, i) = one_o_2p * E_z(i - 1, 0, i - 1) + PA[2] * E_z(i - 1, 0, i);
    }

    for (int j = 1; j <= lb; j++)
        for (int i = 0; i <= la; i++)
        {
            E_x(i, j, 0) = PB[0] * E_x(i, j - 1, 0) + E_x(i, j - 1, 1);
            E_y(i, j, 0) = PB[1] * E_y(i, j - 1, 0) + E_y(i, j - 1, 1);
            E_z(i, j, 0) = PB[2] * E_z(i, j - 1, 0) + E_z(i, j - 1, 1);

            for (int t = 1; t < i + j; t++)
            {
                E_x(i, j, t) = one_o_2p * E_x(i, j - 1, t - 1) +
                               PB[0] * E_x(i, j - 1, t) +
                               double(t + 1) * E_x(i, j - 1, t + 1);

                E_y(i, j, t) = one_o_2p * E_y(i, j - 1, t - 1) +
                               PB[1] * E_y(i, j - 1, t) +
                               double(t + 1) * E_y(i, j - 1, t + 1);

                E_z(i, j, t) = one_o_2p * E_z(i, j - 1, t - 1) +
                               PB[2] * E_z(i, j - 1, t) +
                               double(t + 1) * E_z(i, j - 1, t + 1);
            }

            E_x(i, j, i + j) = one_o_2p * E_x(i, j - 1, i + j - 1) + PB[0] * E_x(i, j - 1, i + j);
            E_y(i, j, i + j) = one_o_2p * E_y(i, j - 1, i + j - 1) + PB[1] * E_y(i, j - 1, i + j);
            E_z(i, j, i + j) = one_o_2p * E_z(i, j - 1, i + j - 1) + PB[2] * E_z(i, j - 1, i + j);
        }
}

template <int la, int lb>
void testCoeffs3DRolledTemplated(const double &a, const double &b,
                                 const std::array<double, 3> &A,
                                 const std::array<double, 3> &B,
                                 vec3d &E_x, vec3d &E_y, vec3d &E_z)
{
    E_x.set(0);
    E_y.set(0);
    E_z.set(0);
    E_x(0, 0, 0) = 1;
    E_y(0, 0, 0) = 1;
    E_z(0, 0, 0) = 1;

    double p = a + b;
    double one_o_2p = 1.0 / (2 * p);

    arma::vec::fixed<3> RA{A[0], A[1], A[2]};
    arma::vec::fixed<3> RB{B[0], B[1], B[2]};

    arma::vec::fixed<3> P = (a * RA + b * RB) / p;
    arma::vec::fixed<3> PA = P - RA;
    arma::vec::fixed<3> PB = P - RB;

    coeffs3DRolledTemplated<la, lb>(one_o_2p, PA, PB, E_x, E_y, E_z);
}

using kernel_t = std::function<void(const double &, const double &,
                                    const std::array<double, 3> &,
                                    const std::array<double, 3> &,
                                    vec3d &, vec3d &, vec3d &)>;

// std::map<std::pair<int, int>, kernel_t> assss{{{1, 1}, testCoeffs3DRolledTemplated<1, 1>}};

std::map<std::pair<int, int>, kernel_t> assss{{{7, 7}, testCoeffs3DRolledTemplated<7, 7>},
                                              {{7, 6}, testCoeffs3DRolledTemplated<7, 6>},
                                              {{7, 5}, testCoeffs3DRolledTemplated<7, 5>},
                                              {{7, 4}, testCoeffs3DRolledTemplated<7, 4>},
                                              {{7, 3}, testCoeffs3DRolledTemplated<7, 3>},
                                              {{7, 2}, testCoeffs3DRolledTemplated<7, 2>},
                                              {{7, 1}, testCoeffs3DRolledTemplated<7, 1>},
                                              {{7, 0}, testCoeffs3DRolledTemplated<7, 0>},
                                              {{6, 6}, testCoeffs3DRolledTemplated<6, 6>},
                                              {{6, 5}, testCoeffs3DRolledTemplated<6, 5>},
                                              {{6, 4}, testCoeffs3DRolledTemplated<6, 4>},
                                              {{6, 3}, testCoeffs3DRolledTemplated<6, 3>},
                                              {{6, 2}, testCoeffs3DRolledTemplated<6, 2>},
                                              {{6, 1}, testCoeffs3DRolledTemplated<6, 1>},
                                              {{6, 0}, testCoeffs3DRolledTemplated<6, 0>},
                                              {{5, 5}, testCoeffs3DRolledTemplated<5, 5>},
                                              {{5, 4}, testCoeffs3DRolledTemplated<5, 4>},
                                              {{5, 3}, testCoeffs3DRolledTemplated<5, 3>},
                                              {{5, 2}, testCoeffs3DRolledTemplated<5, 2>},
                                              {{5, 1}, testCoeffs3DRolledTemplated<5, 1>},
                                              {{5, 0}, testCoeffs3DRolledTemplated<5, 0>},
                                              {{4, 4}, testCoeffs3DRolledTemplated<4, 4>},
                                              {{4, 3}, testCoeffs3DRolledTemplated<4, 3>},
                                              {{4, 2}, testCoeffs3DRolledTemplated<4, 2>},
                                              {{4, 1}, testCoeffs3DRolledTemplated<4, 1>},
                                              {{4, 0}, testCoeffs3DRolledTemplated<4, 0>},
                                              {{3, 3}, testCoeffs3DRolledTemplated<3, 3>},
                                              {{3, 2}, testCoeffs3DRolledTemplated<3, 2>},
                                              {{3, 1}, testCoeffs3DRolledTemplated<3, 1>},
                                              {{3, 0}, testCoeffs3DRolledTemplated<3, 0>},
                                              {{2, 2}, testCoeffs3DRolledTemplated<2, 2>},
                                              {{2, 1}, testCoeffs3DRolledTemplated<2, 1>},
                                              {{2, 0}, testCoeffs3DRolledTemplated<2, 0>},
                                              {{1, 1}, testCoeffs3DRolledTemplated<1, 1>},
                                              {{1, 0}, testCoeffs3DRolledTemplated<1, 0>},
                                              {{0, 0}, testCoeffs3DRolledTemplated<0, 0>}};

// using kernel_t = std::function<void(const double &, const double &,
//                                     const std::array<double, 3> &,
//                                     const std::array<double, 3> &,
//                                     vec3d &, vec3d &, vec3d &)>;

std::map<std::pair<int, int>, kernel_t> assss2;
// assss2[std::make_pair(7, 7)] = [](const double &a, const double &b,
//                                   const std::array<double, 3> &A,
//                                   const std::array<double, 3> &B,
//                                   vec3d &E_x, vec3d &E_y, vec3d &E_z)
// {
//     testCoeffs3DRolledTemplated<7, 7>(a, b, A, B, E_x, E_y, E_z);
// };

inline constexpr int calcIdx3D(const int &i, const int &j, const int &k,
              const int &idim, const int &jdim, const int &kdim)
{
    return i * jdim * kdim + j * kdim + k;
}

void coeffsFlat(const int &la, const int &lb,
                const double one_o_2p, const arma::vec::fixed<3> &PA,
                const arma::vec::fixed<3> &PB,
                std::vector<double> &E_x, std::vector<double> &E_y, std::vector<double> &E_z)
{
    int idim = la + 1;
    int jdim = lb + 1;
    int kdim = la + lb + 1;
    for (int i = 1; i <= la; i++)
    {
        // E_x(i, 0, 0) = PA[0] * E_x(i - 1, 0, 0) + E_x(i - 1, 0, 1);
        // E_y(i, 0, 0) = PA[1] * E_y(i - 1, 0, 0) + E_y(i - 1, 0, 1);
        // E_z(i, 0, 0) = PA[2] * E_z(i - 1, 0, 0) + E_z(i - 1, 0, 1);
        E_x[calcIdx3D(i, 0, 0, idim, jdim, kdim)] = PA[0] * E_x[calcIdx3D(i - 1, 0, 0, idim, jdim, kdim)] +
                                                    E_x[calcIdx3D(i - 1, 0, 1, idim, jdim, kdim)];

        E_y[calcIdx3D(i, 0, 0, idim, jdim, kdim)] = PA[1] * E_y[calcIdx3D(i - 1, 0, 0, idim, jdim, kdim)] +
                                                    E_y[calcIdx3D(i - 1, 0, 1, idim, jdim, kdim)];

        E_z[calcIdx3D(i, 0, 0, idim, jdim, kdim)] = PA[2] * E_z[calcIdx3D(i - 1, 0, 0, idim, jdim, kdim)] +
                                                    E_z[calcIdx3D(i - 1, 0, 1, idim, jdim, kdim)];

        for (int t = 1; t < i; t++)
        {
            // E_x(i, 0, t) = one_o_2p * E_x(i - 1, 0, t - 1) +
            //                PA[0] * E_x(i - 1, 0, t) +
            //                double(t + 1) * E_x(i - 1, 0, t + 1);

            // E_y(i, 0, t) = one_o_2p * E_y(i - 1, 0, t - 1) +
            //                PA[1] * E_y(i - 1, 0, t) +
            //                double(t + 1) * E_y(i - 1, 0, t + 1);

            // E_z(i, 0, t) = one_o_2p * E_z(i - 1, 0, t - 1) +
            //                PA[2] * E_z(i - 1, 0, t) +
            //                double(t + 1) * E_z(i - 1, 0, t + 1);

            // E_x[calcIdx3D(i, 0, t, idim, jdim, kdim)] = one_o_2p * E_x[calcIdx3D(i - 1, 0, t - 1, idim, jdim, kdim)] +
            //                                             PA[0] * E_x[calcIdx3D(i - 1, 0, t, idim, jdim, kdim)] +
            //                                             double(t + 1) * E_x[calcIdx3D(i - 1, 0, t + 1, idim, jdim, kdim)];

            // E_y[calcIdx3D(i, 0, t, idim, jdim, kdim)] = one_o_2p * E_y[calcIdx3D(i - 1, 0, t - 1, idim, jdim, kdim)] +
            //                                             PA[1] * E_y[calcIdx3D(i - 1, 0, t, idim, jdim, kdim)] +
            //                                             double(t + 1) * E_y[calcIdx3D(i - 1, 0, t + 1, idim, jdim, kdim)];

            // E_z[calcIdx3D(i, 0, t, idim, jdim, kdim)] = one_o_2p * E_z[calcIdx3D(i - 1, 0, t - 1, idim, jdim, kdim)] +
            //                                             PA[2] * E_z[calcIdx3D(i - 1, 0, t, idim, jdim, kdim)] +
            //                                             double(t + 1) * E_z[calcIdx3D(i - 1, 0, t + 1, idim, jdim, kdim)];
        }

        // E_x(i, 0, i) = one_o_2p * E_x(i - 1, 0, i - 1) + PA[0] * E_x(i - 1, 0, i);
        // E_y(i, 0, i) = one_o_2p * E_y(i - 1, 0, i - 1) + PA[1] * E_y(i - 1, 0, i);
        // E_z(i, 0, i) = one_o_2p * E_z(i - 1, 0, i - 1) + PA[2] * E_z(i - 1, 0, i);

        // E_x[calcIdx3D(i, 0, i, idim, jdim, kdim)] = one_o_2p * E_x[calcIdx3D(i - 1, 0, i - 1, idim, jdim, kdim)] +
        //                                             PA[0] * E_x[calcIdx3D(i - 1, 0, i, idim, jdim, kdim)];

        // E_y[calcIdx3D(i, 0, i, idim, jdim, kdim)] = one_o_2p * E_y[calcIdx3D(i - 1, 0, i - 1, idim, jdim, kdim)] +
        //                                             PA[1] * E_y[calcIdx3D(i - 1, 0, i, idim, jdim, kdim)];

        // E_z[calcIdx3D(i, 0, i, idim, jdim, kdim)] = one_o_2p * E_z[calcIdx3D(i - 1, 0, i - 1, idim, jdim, kdim)] +
        //                                             PA[2] * E_z[calcIdx3D(i - 1, 0, i, idim, jdim, kdim)];                                                                                                        
    }

    // for (int j = 1; j <= lb; j++)
    //     for (int i = 0; i <= la; i++)
    //     {
    //         // E_x(i, j, 0) = PB[0] * E_x(i, j - 1, 0) + E_x(i, j - 1, 1);
    //         // E_y(i, j, 0) = PB[1] * E_y(i, j - 1, 0) + E_y(i, j - 1, 1);
    //         // E_z(i, j, 0) = PB[2] * E_z(i, j - 1, 0) + E_z(i, j - 1, 1);

    //         E_x[calcIdx3D(i, j, 0, idim, jdim, kdim)] = PB[0] * E_x[calcIdx3D(i, j - 1, 0, idim, jdim, kdim)] +
    //                                                     E_x[calcIdx3D(i, j - 1, 1, idim, jdim, kdim)];

    //         E_y[calcIdx3D(i, j, 0, idim, jdim, kdim)] = PB[1] * E_y[calcIdx3D(i, j - 1, 0, idim, jdim, kdim)] +
    //                                                     E_y[calcIdx3D(i, j - 1, 1, idim, jdim, kdim)];

    //         E_z[calcIdx3D(i, j, 0, idim, jdim, kdim)] = PB[2] * E_z[calcIdx3D(i, j - 1, 0, idim, jdim, kdim)] +
    //                                                     E_z[calcIdx3D(i, j - 1, 1, idim, jdim, kdim)];                                                                                                                

    //         for (int t = 1; t < i + j; t++)
    //         {
    //             // E_x(i, j, t) = one_o_2p * E_x(i, j - 1, t - 1) +
    //             //                PB[0] * E_x(i, j - 1, t) +
    //             //                double(t + 1) * E_x(i, j - 1, t + 1);

    //             // E_y(i, j, t) = one_o_2p * E_y(i, j - 1, t - 1) +
    //             //                PB[1] * E_y(i, j - 1, t) +
    //             //                double(t + 1) * E_y(i, j - 1, t + 1);

    //             // E_z(i, j, t) = one_o_2p * E_z(i, j - 1, t - 1) +
    //             //                PB[2] * E_z(i, j - 1, t) +
    //             //                double(t + 1) * E_z(i, j - 1, t + 1);

    //             E_x[calcIdx3D(i, j, t, idim, jdim, kdim)] = one_o_2p * E_x[calcIdx3D(i, j - 1, t - 1, idim, jdim, kdim)] +
    //                                                         PB[0] * E_x[calcIdx3D(i, j - 1, t, idim, jdim, kdim)] +
    //                                                         double(t + 1) * E_x[calcIdx3D(i, j - 1, t + 1, idim, jdim, kdim)];

    //             E_y[calcIdx3D(i, j, t, idim, jdim, kdim)] = one_o_2p * E_y[calcIdx3D(i, j - 1, t - 1, idim, jdim, kdim)] +
    //                                                         PB[1] * E_y[calcIdx3D(i, j - 1, t, idim, jdim, kdim)] +
    //                                                         double(t + 1) * E_y[calcIdx3D(i, j - 1, t + 1, idim, jdim, kdim)];

    //             E_z[calcIdx3D(i, j, t, idim, jdim, kdim)] = one_o_2p * E_z[calcIdx3D(i, j - 1, t - 1, idim, jdim, kdim)] +
    //                                                         PB[2] * E_z[calcIdx3D(i, j - 1, t, idim, jdim, kdim)] +
    //                                                         double(t + 1) * E_z[calcIdx3D(i, j - 1, t + 1, idim, jdim, kdim)];                                                                                                                        
    //         }

    //         // E_x(i, j, i + j) = one_o_2p * E_x(i, j - 1, i + j - 1) + PB[0] * E_x(i, j - 1, i + j);
    //         // E_y(i, j, i + j) = one_o_2p * E_y(i, j - 1, i + j - 1) + PB[1] * E_y(i, j - 1, i + j);
    //         // E_z(i, j, i + j) = one_o_2p * E_z(i, j - 1, i + j - 1) + PB[2] * E_z(i, j - 1, i + j);

    //         E_x[calcIdx3D(i, j, i + j, idim, jdim, kdim)] = one_o_2p * E_x[calcIdx3D(i, j - 1, i + j - 1, idim, jdim, kdim)] +
    //                                                         PB[0] * E_x[calcIdx3D(i, j - 1, i + j, idim, jdim, kdim)];

    //         E_y[calcIdx3D(i, j, i + j, idim, jdim, kdim)] = one_o_2p * E_y[calcIdx3D(i, j - 1, i + j - 1, idim, jdim, kdim)] +
    //                                                         PB[1] * E_y[calcIdx3D(i, j - 1, i + j, idim, jdim, kdim)];

    //         E_z[calcIdx3D(i, j, i + j, idim, jdim, kdim)] = one_o_2p * E_z[calcIdx3D(i, j - 1, i + j - 1, idim, jdim, kdim)] +
    //                                                         PB[2] * E_z[calcIdx3D(i, j - 1, i + j, idim, jdim, kdim)];                                                                                                                        
    //     }
}

void testCoeffsFlat(const int &la, const int &lb, 
                    const double &a, const double &b,
                    const std::array<double, 3> &A,
                    const std::array<double, 3> &B,
                    std::vector<double> &E_x, std::vector<double> &E_y, std::vector<double> &E_z)
{
    // std::fill(E_x.begin(), E_x.end(), 0);
    // std::fill(E_y.begin(), E_y.end(), 0);
    // std::fill(E_z.begin(), E_z.end(), 0);
    // memset(&E_x[0], 0, sizeof(E_x[0]) * E_x.size());
    // memset(&E_y[0], 0, sizeof(E_y[0]) * E_y.size());
    // memset(&E_z[0], 0, sizeof(E_z[0]) * E_z.size());
    E_x[0] = 1;
    E_y[0] = 1;
    E_z[0] = 1;

    double p = a + b;
    double one_o_2p = 1.0 / (2 * p);

    arma::vec::fixed<3> RA{A[0], A[1], A[2]};
    arma::vec::fixed<3> RB{B[0], B[1], B[2]};

    arma::vec::fixed<3> P = (a * RA + b * RB) / p;
    arma::vec::fixed<3> PA = P - RA;
    arma::vec::fixed<3> PB = P - RB;

    coeffsFlat(la, lb, one_o_2p, PA, PB, E_x, E_y, E_z);
}

template <int la>
void coeffsMeta1(const int &idx, const int &idim, const int &jdim, const int &kdim,
                 const arma::vec::fixed<3> &PA,
                 std::vector<double> &E)
{
    coeffsMeta1<la - 1>(idx, idim, jdim, kdim, PA, E);
    E[calcIdx3D(la, 0, 0, idim, jdim, kdim)] = PA[idx] * E[calcIdx3D(la - 1, 0, 0, idim, jdim, kdim)] +
                                               E[calcIdx3D(la - 1, 0, 1, idim, jdim, kdim)];
}

template <>
void coeffsMeta1<1>(const int &idx, const int &idim, const int &jdim, const int &kdim,
                   const arma::vec::fixed<3> &PA,
                   std::vector<double> &E)
{
    E[calcIdx3D(1, 0, 0, idim, jdim, kdim)] = PA[idx] * E[0];
}

template <>
void coeffsMeta1<0>(const int &idx, const int &idim, const int &jdim, const int &kdim,
                   const arma::vec::fixed<3> &PA,
                   std::vector<double> &E)
{
    E[calcIdx3D(1, 0, 0, idim, jdim, kdim)] = PA[idx] * E[0];
}

template <int la, int lb>
void coeffsMeta(const int &idim, const int &jdim, const int &kdim,                
                const double one_o_2p, 
                const arma::vec::fixed<3> &PA,
                const arma::vec::fixed<3> &PB,
                std::vector<double> &E_x, std::vector<double> &E_y, std::vector<double> &E_z)
{
    if (la > 1)
    {
        coeffsMeta1<la>(0, idim, jdim, kdim, PA, E_x);
        coeffsMeta1<la>(1, idim, jdim, kdim, PA, E_y);
        coeffsMeta1<la>(2, idim, jdim, kdim, PA, E_z);
    }
    else if (la == 1)
    {
        coeffsMeta1<1>(0, idim, jdim, kdim, PA, E_x);
        coeffsMeta1<1>(1, idim, jdim, kdim, PA, E_y);
        coeffsMeta1<1>(2, idim, jdim, kdim, PA, E_z);
    }
}

template <int la, int lb>
void coeffsMetaTest(const double &a, const double &b,
                    const std::array<double, 3> &A,
                    const std::array<double, 3> &B,
                    std::vector<double> &E_x, std::vector<double> &E_y, std::vector<double> &E_z)
{
    constexpr int idim = la + 1;
    constexpr int jdim = lb + 1;
    constexpr int kdim = la + lb + 1;

    E_x[0] = 1;
    E_y[0] = 1;
    E_z[0] = 1;

    double p = a + b;
    double one_o_2p = 1.0 / (2 * p);

    arma::vec::fixed<3> RA{A[0], A[1], A[2]};
    arma::vec::fixed<3> RB{B[0], B[1], B[2]};

    arma::vec::fixed<3> P = (a * RA + b * RB) / p;
    arma::vec::fixed<3> PA = P - RA;
    arma::vec::fixed<3> PB = P - RB;

    coeffsMeta<la, lb>(idim, jdim, kdim, one_o_2p, PA, PB, E_x, E_y, E_z);
}

using meta_kernel_t = std::function<void(const double &, const double &,
                                         const std::array<double, 3> &,
                                         const std::array<double, 3> &,
                                         std::vector<double> &, std::vector<double> &, std::vector<double> &)>;

// std::map<std::pair<int, int>, meta_kernel_t> metass{{{0, 0}, coeffsMetaTest<0, 0>}};

std::map<std::pair<int, int>, meta_kernel_t> metass{{{7, 7}, coeffsMetaTest<7, 7>},
                                                    {{7, 6}, coeffsMetaTest<7, 6>},
                                                    {{7, 5}, coeffsMetaTest<7, 5>},
                                                    {{7, 4}, coeffsMetaTest<7, 4>},
                                                    {{7, 3}, coeffsMetaTest<7, 3>},
                                                    {{7, 2}, coeffsMetaTest<7, 2>},
                                                    {{7, 1}, coeffsMetaTest<7, 1>},
                                                    {{7, 0}, coeffsMetaTest<7, 0>},
                                                    {{6, 6}, coeffsMetaTest<6, 6>},
                                                    {{6, 5}, coeffsMetaTest<6, 5>},
                                                    {{6, 4}, coeffsMetaTest<6, 4>},
                                                    {{6, 3}, coeffsMetaTest<6, 3>},
                                                    {{6, 2}, coeffsMetaTest<6, 2>},
                                                    {{6, 1}, coeffsMetaTest<6, 1>},
                                                    {{6, 0}, coeffsMetaTest<6, 0>},
                                                    {{5, 5}, coeffsMetaTest<5, 5>},
                                                    {{5, 4}, coeffsMetaTest<5, 4>},
                                                    {{5, 3}, coeffsMetaTest<5, 3>},
                                                    {{5, 2}, coeffsMetaTest<5, 2>},
                                                    {{5, 1}, coeffsMetaTest<5, 1>},
                                                    {{5, 0}, coeffsMetaTest<5, 0>},
                                                    {{4, 4}, coeffsMetaTest<4, 4>},
                                                    {{4, 3}, coeffsMetaTest<4, 3>},
                                                    {{4, 2}, coeffsMetaTest<4, 2>},
                                                    {{4, 1}, coeffsMetaTest<4, 1>},
                                                    {{4, 0}, coeffsMetaTest<4, 0>},
                                                    {{3, 3}, coeffsMetaTest<3, 3>},
                                                    {{3, 2}, coeffsMetaTest<3, 2>},
                                                    {{3, 1}, coeffsMetaTest<3, 1>},
                                                    {{3, 0}, coeffsMetaTest<3, 0>},
                                                    {{2, 2}, coeffsMetaTest<2, 2>},
                                                    {{2, 1}, coeffsMetaTest<2, 1>},
                                                    {{2, 0}, coeffsMetaTest<2, 0>},
                                                    {{1, 1}, coeffsMetaTest<1, 1>},
                                                    {{1, 0}, coeffsMetaTest<1, 0>},
                                                    {{0, 0}, coeffsMetaTest<0, 0>}};

void LIMD::test()
{
    typedef std::chrono::nanoseconds ns;
    typedef std::chrono::microseconds mus;
    printf("Testing LIMD::coeffs\n");

    double a = 2.0;
    double b = 3.0;
    std::array<double, 3> A{0.5, 1.1, -0.4};
    std::array<double, 3> B{1.5, -1.2, 0.7};

    int l_max = 7;
    for (int la = l_max; la >= 0; la--)
        for (int lb = la; lb >= 0; lb--)
        {
            // vec3d E_x(la + 1, lb + 1, la + lb + 1, 0);
            // vec3d E_y(la + 1, lb + 1, la + lb + 1, 0);
            // vec3d E_z(la + 1, lb + 1, la + lb + 1, 0);
            size_t dim_E = (la + 1) * (lb + 1) * (la + lb + 1);
            std::vector<double> E_x(dim_E, 0);
            std::vector<double> E_y(dim_E, 0);
            std::vector<double> E_z(dim_E, 0);

            // vec4d E_xyz(3, la + 1, lb + 1, la + lb + 1, 0);

            auto start{std::chrono::steady_clock::now()};

            // kernel_t kernel = testCoeffs3DRolledTemplated<, 1>;
            // kernel_t kernel = assss.at(std::make_pair(la, lb));
            meta_kernel_t kernel = metass.at(std::make_pair(la, lb));

            // printf("\nla = %d, lb = %d", la, lb);
            // auto kernel = std::any_cast<kernel_t>(assss.at(std::make_pair(la, lb)));

            for (int i = 0; i < 1000; i++)
                kernel(a, b, A, B, E_x, E_y, E_z);
                // testCoeffsFlat(la, lb, a, b, A, B, E_x, E_y, E_z);                

                // testCoeffsFlat(la, lb, a, b, A, B, E_x, E_y, E_z);

                // testCoeffsFlat(la, lb, a, b, A, B, E_x, E_y, E_z);
                // kernel(a, b, A, B, E_x, E_y, E_z);
                // testCoeffs3DRolled(a, b, la, lb, A, B, E_x, E_y, E_z);
                // coeffs(a, b, la, lb, A, B, E_x, E_y, E_z);                
                // kernel(a, b, A, B, E_x, E_y, E_z);

            // coeffs(a, b, la, lb, A, B, E_x, E_y, E_z);
            // testCoeffs4D(a, b, la, lb, A, B, E_xyz);
            // testCoeffs3DRolled(a, b, la, lb, A, B, E_x, E_y, E_z);

            auto end{std::chrono::steady_clock::now()};
            std::chrono::duration<double> duration{end - start};

            double sum_coeffs = 0;
            for (size_t i = 0; i < dim_E; i++)
            {
                sum_coeffs += E_x[i];
                sum_coeffs += E_y[i];
                sum_coeffs += E_z[i];
            }

            // for (int i = 0; i <= la; i++)
            //     for (int j = 0; j <= lb; j++)
            //         for (int t = 0; t <= (la + lb); t++)
            //             sum_coeffs += E_x(i, j, t);

            // for (int i = 0; i <= la; i++)
            //     for (int j = 0; j <= lb; j++)
            //         for (int t = 0; t <= (la + lb); t++)
            //             sum_coeffs += E_y(i, j, t);

            // for (int i = 0; i <= la; i++)
            //     for (int j = 0; j <= lb; j++)
            //         for (int t = 0; t <= (la + lb); t++)
            //             sum_coeffs += E_z(i, j, t);

            // for (int x = 0; x < 3; x++)
            //     for (int i = 0; i <= la; i++)
            //         for (int j = 0; j <= lb; j++)
            //             for (int t = 0; t <= (la + lb); t++)
            //                 sum_coeffs += E_xyz(x, i, j, t);

            auto duration_ns = std::chrono::duration_cast<ns>(duration);
            auto duration_mus = std::chrono::duration_cast<mus>(duration);
            std::cout << la << ", " << lb << ", " << duration_ns.count() << ", " << sum_coeffs << "\n";
        }
}