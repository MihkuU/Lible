#pragma once

#include <lible/types.hpp>
#include <lible/ints/shell_pair_data.hpp>

#include <vector>

namespace lible
{
    namespace ints
    {
        /** */
        vec3d ecoeffsRecurrence2(const double a, const double b, const int la, const int lb,
                                 const double PA, const double PB, const double Kab);

        /** */
        vec3d ecoeffsRecurrence2_n1(const double a, const double b, const int la, const int lb,
                                    const double A, const double B, const vec3d &ecoeffs);

        /** */
        vec2d ecoeffsRecurrence1(const double one_o_2a, const int l);

        /** */
        std::array<vec3d, 3> ecoeffsPrimitivePair(const double a, const double b, const int la,
                                                  const int lb, const double *xyz_a,
                                                  const double *xyz_b);

        /** */
        std::array<vec3d, 3> ecoeffsPrimitivePair_n1(const double a, const double b, const int la,
                                                     const int lb, const double *xyz_a,
                                                     const double *xyz_b,
                                                     const std::array<lible::vec3d, 3> &ecoeffs);

        /** */
        std::array<lible::vec2d, 3> ecoeffsPrimitive(const double a, const int l);

        /** */
        std::vector<lible::vec3d>
        ecoeffsShellPair_Eij0(const int la, const int lb, const int cdepth_a, const int cdepth_b,
                              const double *exps_a, const double *exps_b, const double *xyz_a,
                              const double *xyz_b);

        /** */
        std::vector<lible::vec4d>
        ecoeffsShellPair_Eijt(const int la, const int lb, const int cdepth_a, const int cdepth_b,
                              const double *exps_a, const double *exps_b, const double *xyz_a,
                              const double *xyz_b);

        /** TODO: remove */
        std::vector<lible::vec4d>
        ecoeffsShellPair_Eijt_Debug(const int la, const int lb, const int cdepth_a, const int cdepth_b,
                                    const double *exps_a, const double *exps_b, const double *xyz_a,
                                    const double *xyz_b);

        /** */
        std::vector<std::vector<double>> ecoeffsShell(const int l, const std::vector<double> &exps);

        /**
         * Calculates the Hermite expansion coefficients { E(r, i, j, 0) | r = x,y,z } for each
         * pair of Gaussian primitives in each shell pair.
         */
        std::vector<std::vector<vec3d>>
        ecoeffsSPData_Eij0(const ShellPairData &sp_data);

        /**
         * Calculates the Hermite expansion coefficients { E(r, i, j, 0) | r = x,y,z } for each
         * pair of Gaussian primitives in each shell pair.
         */
        std::vector<std::vector<vec4d>>
        ecoeffsSPData_Eijt(const ShellPairData &sp_data);

        /** TODO: mention that c-coeffs are inside the e-coeffs!!!! */
        std::vector<double>
        ecoeffsSphericalSPData_Bra(const ShellPairData &sp_data);

        /** TODO: mention that c-coeffs are inside the e-coeffs!!!! */
        std::vector<double>
        ecoeffsSphericalSPData_Bra_Deriv1(const ShellPairData &sp_data);

        /** TODO: mention that c-coeffs are inside the e-coeffs!!!! */
        std::pair<std::vector<double>, std::vector<double>>
        ecoeffsSphericalSPData_BraKet(const ShellPairData &sp_data);

        /** TODO: mention that c-coeffs are inside the e-coeffs!!!! */
        std::pair<std::vector<double>, std::vector<double>>
        ecoeffsSphericalSPData_BraKet_Deriv1(const ShellPairData &sp_data);

        /** TODO: mention that c-coeffs are inside the e-coeffs!!!! */
        std::vector<double>
        ecoeffsSphericalShellData_Bra(const ShellData &sh_data);

        /** TODO: mention that c-coeffs are inside the e-coeffs!!!! */
        std::pair<std::vector<double>, std::vector<double>>
        ecoeffsSphericalShellData_BraKet(const ShellData &sh_data);

        /** TODO: mention that c-coeffs are inside the e-coeffs!!!! */
        std::vector<std::vector<double>>
        ecoeffsSphericalSPDatas_Bra(const std::vector<ShellPairData> &sp_datas);

        /** TODO: mention that c-coeffs are inside the e-coeffs!!!! */
        std::pair<std::vector<std::vector<double>>, std::vector<std::vector<double>>>
        ecoeffsSphericalSPDatas_BraKet(const std::vector<ShellPairData> &sp_datas);

        ////////////////////////////////
        std::pair<std::vector<double>, std::vector<double>>
        ecoeffsSphericalSPData_BraKet(const ShellPairData &sp_data);

        ////////////////////////////////
        std::vector<double>
        ecoeffsSHARKBig(const ShellPairData &sp_data);

        std::vector<double>
        ecoeffsSHARKBig(const ShellData &sh_data);

        ///// only this >>>>>>>>>>>>>>>>>>>>>

        std::vector<double>
        ecoeffsSHARKSmall(const ShellPairData &sp_data, const bool transpose = false);

        std::vector<double>
        ecoeffsSHARKSmall(const ShellData &sh_data, const bool transpose = false);

        std::vector<double>
        ecoeffsD1SHARKSmall(const ShellPairData &sp_data, const bool transpose = false);

        ///// only this <<<<<<<<<<<<<<<<<<<<<

        std::vector<double> 
        ecoeffsD1SHARK(const ShellPairData &sp_data);

        std::vector<double>
        ecoeffsSHARKShellPair(const int ipair, const ShellPairData &sp_data);

        std::vector<double>
        ecoeffsD1SHARKShellPair(const int ipair, const ShellPairData &sp_data);        

        std::vector<double>
        ecoeffsSHARKBigShell(const int ishell, const ShellData &sh_data);

    }
}