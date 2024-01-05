#include <lible/gci.h>

#include <lible/util.h>
#include <lible/guga_sci.h>

namespace LG = lible::guga;

using namespace lible;

LG::GCI::GCI()
{
}

LG::GCI::GCI(std::unique_ptr<GCI::Impl> &&impl) : impl(std::move(impl))
{
}

LG::GCI::~GCI() = default;

LG::GCI::GCI(GCI &&) noexcept = default;
LG::GCI &LG::GCI::operator=(GCI &&) noexcept = default;

LG::GCI LG::run(const int &n_orbs, const int &n_els,
                const int &n_roots, const int &multiplicity,
                const vec2d &one_el_ints, const vec4d &two_el_ints,
                std::vector<double> &ci_energies_out,
                std::vector<std::vector<double>> &ci_vectors_out)
{
    std::unique_ptr<GCI::Impl> impl = std::make_unique<GCI::Impl>(n_orbs, n_els,
                                                                  n_roots, multiplicity,
                                                                  one_el_ints, two_el_ints);

    impl->run(ci_energies_out, ci_vectors_out);

    return GCI(std::move(impl));
}

vec2d LG::calc1RDM(const GCI &gci, const size_t &iroot, const size_t &jroot)
{
}

vec2d LG::calc1SRDM(const GCI &gci, const size_t &iroot)
{
}

vec4d LG::calc2RDM(const GCI &gci, const size_t &iroot, const size_t &jroot)
{
}

std::pair<vec2d, vec4d>
LG::calc12RDMs(const GCI &gci, const size_t &iroot, const size_t &jroot)
{
}

std::vector<double> LG::calcSigma()
{
}

std::vector<std::vector<std::tuple<std::string, std::string, double>>>
LG::returnSignificantCSFs(const GCI &gci)
{
}