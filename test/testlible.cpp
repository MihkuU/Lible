
#include <lible>

using namespace Lible;
using std::vector;

void singlePointCalc(double &energy, vector<double> geometry, vector<double> gradient)
{

}

void testLibleGeomOpt()
{
    // GeomOpt geom_opt(vector<double>{});

    // auto singlePointCalculation =
    // std::function<void(double &energy, std::vector<double> &geometry, std::vector<double> &gradient)> singlePointCalculation =
        // [](double &energy, std::vector<double> &geometry, std::vector<double> &gradient) {};
    // std::vector<double> sitt = geom_opt.sitt();
    // vector<double> optimized_geometry = geom_opt.optimize<GeomOpt::Option::KRIGING>(
    //     [](double &energy, vector<double> &geometry, vector<double> &gradient)
    //     { singlePointCalc(energy, geometry, gradient); });
}

void testLibleInts()
{
    Ints ints("def2-svp", vector<double>{0, 0, 0, 0, 0, 1}, vector<std::string>{"H", "H"});

    vector<double> one_el_ints = ints.calcOneElInts<Ints::Option1El::OVERLAP>();
}

int main()
{
    // testLibleGeomOpt();

    testLibleInts();
}
