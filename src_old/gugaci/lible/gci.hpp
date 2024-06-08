#pragma once

#include <memory>
#include <tuple>
#include <utility>
#include <vector>

#include <lible/davidson_settings.hpp>
#include <lible/gci_settings.hpp>
#include <lible/types.hpp>

namespace lible
{
    namespace guga
    {
        class GCI;

        GCI run(const int &n_orbs, const int &n_els,
                const int &n_roots, const int &multiplicity,
                const vec2d &one_el_ints, const vec4d &two_el_ints,
                std::vector<double> &ci_energies_out,
                std::vector<std::vector<double>> &ci_vectors_out,
                const double &core_energy = 0);

        GCI runFromCSFs(const int &n_orbs, const int &n_els,
                        const int &n_roots, const int &multiplicity,
                        const vec2d &one_el_ints, const vec4d &two_el_ints,
                        const std::vector<std::string> &csfs,
                        std::vector<double> &ci_energies_out,
                        std::vector<std::vector<double>> &ci_vectors_out,
                        const double &core_energy = 0);

        GCI runFromCSFsFile(const int &n_orbs, const int &n_els,
                            const int &n_roots, const int &multiplicity,
                            const vec2d &one_el_ints, const vec4d &two_el_ints,
                            const std::string &csfs_fname,
                            std::vector<double> &ci_energies_out,
                            std::vector<std::vector<double>> &ci_vectors_out,
                            const double &core_energy = 0);

        vec2d calc1RDM(const GCI &gci, const size_t &iroot, const size_t &jroot);

        vec2d calc1SRDM(const GCI &gci, const size_t &iroot);

        vec4d calc2RDM(const GCI &gci, const size_t &iroot, const size_t &jroot);

        std::pair<vec2d, vec4d> calc12RDMs(const GCI &gci, const size_t &iroot,
                                           const size_t &jroot);

        std::vector<double> calcSigma(const GCI &gci);

        std::vector<std::vector<std::tuple<std::string, std::string, double>>>
        returnSignificantCSFs(const GCI &gci);

        class GCI
        {
        public:
            class Impl;

            GCI();
            GCI(std::unique_ptr<Impl> &&impl);
            ~GCI();

            GCI(const GCI &) = delete;
            GCI &operator=(const GCI &) = delete;

            GCI(GCI &&) noexcept;
            GCI &operator=(GCI &&) noexcept;

        private:
            std::unique_ptr<Impl> impl;
        };
    }
}