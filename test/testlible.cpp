
#include <chrono> // Change <> to ""
#include <armadillo>

#include "lible"

using namespace Lible;
using std::vector;
using namespace std::chrono;

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
    // Ints ints("def2-svp", vector<double>{0, 0, 0, 0, 0, 1}, vector<std::string>{"H", "H"});
    Ints ints("def2-qzvp", vector<double>{0, 0, 0, 0, 0, 1}, vector<std::string>{"H", "H"});
    // Ints ints("def2-qzvp", vector<double>{0, 0, 0, 0, 0, 1}, vector<std::string>{"C", "C"});
    // Ints ints("aug-cc-pv6z", vector<double>{0, 0, 0, 0, 0, 1}, vector<std::string>{"H", "H"});
    // Ints ints("aug-cc-pv6z", vector<double>{0, 0, 0, 0, 0, 1}, vector<std::string>{"C", "C"});

    for (int i = 0; i < 10; i++)
    {
        auto start = high_resolution_clock::now();
        // vector<double> one_el_ints = ints.calcOneElInts<Ints::Option1El::OVERLAP>();
        vector<double> one_el_ints = ints.calcOneElInts<Ints::OVERLAP>();        
        auto stop = high_resolution_clock::now();
        auto duration = duration_cast<nanoseconds>(stop - start);
        std::cout << "duration.count() = " << duration.count() << " nanoseconds" << std::endl;

        std::size_t dim = std::sqrt(one_el_ints.size());
        arma::dmat one_el_ints_arma(dim, dim, arma::fill::zeros);
        for (std::size_t i = 0, ij = 0; i < dim; i++)
            for (std::size_t j = 0; j < dim; j++, ij++)
                one_el_ints_arma(i, j) = one_el_ints[ij];

        double sum_one_el_ints = arma::accu(one_el_ints_arma);
        printf("sum_one_el_ints = %16.12lf\n", sum_one_el_ints);
        // std::cout << "sum_one_el_ints = " << sum_one_el_ints << std::endl;
        // std::cout << one_el_ints_arma << std::endl;
    }

    // arma::dmat one_el_ints_arma = arma::conv_to<arma::mat>::from(one_el_ints);
    // std:: cout << "one_el_ints_arma.n_elem = " << one_el_ints_arma.n_elem << std::endl;
    // std::cout << one_el_ints_arma << std::endl;
    
}

int main()
{
    // testLibleGeomOpt();

    testLibleInts();
}
