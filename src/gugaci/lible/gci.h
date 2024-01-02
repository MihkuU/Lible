#pragma once

#include <memory>
#include <utility>
#include <vector>

#include <lible/types.h>

namespace lible
{
    namespace guga
    {
        class GCICalc;

        GCICalc runGCI(const int &n_orbs, const int &n_els, const int &n_roots,
                       const int &multiplicity, const Vec2D<double> &one_el_ints,
                       const Vec4D<double> &two_el_ints);

        vec2d calc1RDM(const GCICalc &gci_calc);

        vec2d calc1SRDM(const GCICalc &gci_calc);

        vec4d calc2RDM(const GCICalc &gci_calc);

        std::pair<vec2d, vec4d> calc12RDMs(const GCICalc &gci_calc);

        std::vector<double> calcSigma();

        class SCI;
        class GCICalc
        {
        public:
            GCICalc() = delete;
            GCICalc(std::unique_ptr<SCI> &&sci);
            ~GCICalc();

            GCICalc(const GCICalc &) = delete;
            GCICalc &operator=(const GCICalc &) = delete;

            GCICalc(GCICalc &&) noexcept;
            GCICalc &operator=(GCICalc &&) noexcept;

        private:
            std::unique_ptr<SCI> sci;
        };
    }
}