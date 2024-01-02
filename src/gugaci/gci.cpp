#include <lible/gci.h>

#include <lible/util.h>

#include <lible/cipsi.h>
#include <lible/coupling_coeffs.h>
#include <lible/connections.h>
#include <lible/prefix_algorithm.h>
#include <lible/guga_sci.h>

using namespace lible;
using namespace lible::guga;

GCICalc::GCICalc(std::unique_ptr<SCI> &&sci) : sci(std::move(sci))
{
    // sci = std::move(sci_);
}

GCICalc::~GCICalc() = default;

GCICalc::GCICalc(GCICalc &&) noexcept = default;
GCICalc &GCICalc::operator=(GCICalc &&) noexcept = default;

GCICalc lible::guga::runGCI(const int &n_orbs, const int &n_els, const int &n_roots,
                            const int &multiplicity, const Vec2D<double> &one_el_ints,
                            const Vec4D<double> &two_el_ints)
{
    std::unique_ptr<SCI> sci = std::make_unique<SCI>(n_orbs, n_els, n_roots, multiplicity,
                                                     std::vector<double>({1, 2}),
                                                     std::vector<double>({1, 2, 3, 4}));

    palPrint("Calling runGCI!\n");                                                    

    return GCICalc(std::move(sci));
}