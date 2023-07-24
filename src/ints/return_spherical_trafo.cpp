#include "spherical_trafo.h"

using namespace Lible;

typedef std::tuple<std::size_t, std::size_t, double> trafo_coeff_tuple;

std::vector<trafo_coeff_tuple> SphericalTrafo::returnSphericalTrafo(const int &angmom)
{
    std::vector<trafo_coeff_tuple> spherical_trafo;

    /*
     * Quick word on conventions used here:
     * The transformation from cartesian to spherical gaussian basis sets is stored
     * in a vector of tuples with three values.
     * Values in tuple:
     *   first integer  - spherical gaussian index
     *   second integer - cartesian gaussian index
     *   third double   - transformation coefficient
     */
    switch (angmom)
    {
    case (0):
    {
        spherical_trafo = std::vector<trafo_coeff_tuple>{
            {0, 0, 1.0}};
        break;
    }

    case (1):
    {
        spherical_trafo = std::vector<trafo_coeff_tuple>{
            {0, 2, 1.0},
            {1, 0, 1.0},
            {2, 1, 1.0}};
        break;
    }

    case (2):
    {
        spherical_trafo = std::vector<trafo_coeff_tuple>{
            {0, 5, 1.0},
            {0, 0, -0.5},
            {0, 3, -0.5},
            {1, 2, 1.7320508075688772},
            {2, 4, 1.7320508075688772},
            {3, 0, 0.8660254037844386},
            {3, 3, -0.8660254037844386},
            {4, 1, 1.7320508075688772},
        };
        break;
    }

    case (3):
    {
        spherical_trafo = std::vector<trafo_coeff_tuple>{
            {0, 9, 1.0},
            {0, 2, -1.5},
            {0, 7, -1.5},
            {1, 5, 2.449489742783178},
            {1, 0, -0.6123724356957945},
            {1, 3, -0.6123724356957945},
            {2, 8, 2.449489742783178},
            {2, 1, -0.6123724356957945},
            {2, 6, -0.6123724356957945},
            {3, 2, 1.9364916731037085},
            {3, 7, -1.9364916731037085},
            {4, 4, 3.872983346207417},
            {5, 0, 0.7905694150420948},
            {5, 3, -2.3717082451262845},
            {6, 1, 2.3717082451262845},
            {6, 6, -0.7905694150420948}};
        break;
    }

    case (4):
    {
        spherical_trafo = std::vector<trafo_coeff_tuple>{
            {0, 14, 1.0},
            {0, 5, -3.0},
            {0, 12, -3.0},
            {0, 0, 0.375},
            {0, 3, 0.75},
            {0, 10, 0.375},
            {1, 9, 3.162277660168379},
            {1, 2, -2.3717082451262845},
            {1, 7, -2.3717082451262845},
            {2, 13, 3.162277660168379},
            {2, 4, -2.3717082451262845},
            {2, 11, -2.3717082451262845},
            {3, 5, 3.3541019662496847},
            {3, 12, -3.3541019662496847},
            {3, 0, -0.5590169943749475},
            {3, 3, 0.5590169943749475},
            {3, 3, -0.5590169943749475},
            {3, 10, 0.5590169943749475},
            {4, 8, 6.708203932499369},
            {4, 1, -1.118033988749895},
            {4, 6, -1.118033988749895},
            {5, 2, 2.0916500663351885},
            {5, 7, -6.274950199005565},
            {6, 4, 6.274950199005565},
            {6, 11, -2.0916500663351885},
            {7, 0, 0.739509972887452},
            {7, 3, -4.437059837324712},
            {7, 10, 0.739509972887452},
            {8, 1, 2.958039891549808},
            {8, 6, -2.958039891549808}};
        break;
    }

    case (5):
    {
        spherical_trafo = std::vector<trafo_coeff_tuple>{
            {0, 20, 1.0},
            {0, 9, -5.0},
            {0, 18, -5.0},
            {0, 2, 1.875},
            {0, 7, 3.75},
            {0, 16, 1.875},
            {1, 14, 3.8729833462074166},
            {1, 5, -5.809475019311125},
            {1, 12, -5.809475019311125},
            {1, 0, 0.4841229182759271},
            {1, 3, 0.9682458365518541},
            {1, 10, 0.4841229182759271},
            {2, 19, 3.8729833462074166},
            {2, 8, -5.809475019311125},
            {2, 17, -5.809475019311125},
            {2, 1, 0.4841229182759271},
            {2, 6, 0.9682458365518541},
            {2, 15, 0.4841229182759271},
            {3, 9, 5.1234753829798},
            {3, 18, -5.1234753829798},
            {3, 2, -2.5617376914899},
            {3, 7, 2.5617376914899},
            {3, 7, -2.5617376914899},
            {3, 16, 2.5617376914899},
            {4, 13, 10.2469507659596},
            {4, 4, -5.1234753829798},
            {4, 11, -5.1234753829798},
            {5, 5, 4.183300132670377},
            {5, 12, -12.549900398011133},
            {5, 0, -0.5229125165837971},
            {5, 3, 1.5687375497513916},
            {5, 3, -0.5229125165837971},
            {5, 10, 1.5687375497513916},
            {6, 8, 12.549900398011133},
            {6, 17, -4.183300132670377},
            {6, 1, -1.5687375497513916},
            {6, 6, 0.5229125165837971},
            {6, 6, -1.5687375497513916},
            {6, 15, 0.5229125165837971},
            {7, 2, 2.218529918662356},
            {7, 7, -13.311179511974137},
            {7, 16, 2.218529918662356},
            {8, 4, 8.874119674649425},
            {8, 11, -8.874119674649425},
            {9, 0, 0.7015607600201139},
            {9, 3, -7.01560760020114},
            {9, 10, 3.50780380010057},
            {10, 1, 3.50780380010057},
            {10, 6, -7.01560760020114},
            {10, 15, 0.7015607600201139}};
        break;
    }

    case (6):
    {
        spherical_trafo = std::vector<trafo_coeff_tuple>{
            {0, 27, 1.0},
            {0, 14, -7.5},
            {0, 25, -7.5},
            {0, 5, 5.625},
            {0, 12, 11.25},
            {0, 23, 5.625},
            {0, 0, -0.3125},
            {0, 3, -0.9375},
            {0, 10, -0.9375},
            {0, 21, -0.3125},
            {1, 20, 4.582575694955841},
            {1, 9, -11.456439237389601},
            {1, 18, -11.456439237389601},
            {1, 2, 2.8641098093474002},
            {1, 7, 5.7282196186948005},
            {1, 16, 2.8641098093474002},
            {2, 26, 4.582575694955841},
            {2, 13, -11.456439237389601},
            {2, 24, -11.456439237389601},
            {2, 4, 2.8641098093474002},
            {2, 11, 5.7282196186948005},
            {2, 22, 2.8641098093474002},
            {3, 14, 7.245688373094719},
            {3, 25, -7.245688373094719},
            {3, 5, -7.245688373094719},
            {3, 12, 7.245688373094719},
            {3, 12, -7.245688373094719},
            {3, 23, 7.245688373094719},
            {3, 0, 0.4528555233184199},
            {3, 3, -0.4528555233184199},
            {3, 3, 0.9057110466368399},
            {3, 10, -0.9057110466368399},
            {3, 10, 0.4528555233184199},
            {3, 21, -0.4528555233184199},
            {4, 19, 14.491376746189438},
            {4, 8, -14.491376746189438},
            {4, 17, -14.491376746189438},
            {4, 1, 0.9057110466368399},
            {4, 6, 1.8114220932736798},
            {4, 15, 0.9057110466368399},
            {5, 9, 7.24568837309472},
            {5, 18, -21.73706511928416},
            {5, 2, -2.71713313991052},
            {5, 7, 8.15139941973156},
            {5, 7, -2.71713313991052},
            {5, 16, 8.15139941973156},
            {6, 13, 21.73706511928416},
            {6, 24, -7.24568837309472},
            {6, 4, -8.15139941973156},
            {6, 11, 2.71713313991052},
            {6, 11, -8.15139941973156},
            {6, 22, 2.71713313991052},
            {7, 5, 4.960783708246107},
            {7, 12, -29.764702249476645},
            {7, 23, 4.960783708246107},
            {7, 0, -0.4960783708246108},
            {7, 3, 2.9764702249476644},
            {7, 10, -0.4960783708246108},
            {7, 3, -0.4960783708246108},
            {7, 10, 2.9764702249476644},
            {7, 21, -0.4960783708246108},
            {8, 8, 19.84313483298443},
            {8, 17, -19.84313483298443},
            {8, 1, -1.984313483298443},
            {8, 6, 1.984313483298443},
            {8, 6, -1.984313483298443},
            {8, 15, 1.984313483298443},
            {9, 2, 2.3268138086232857},
            {9, 7, -23.268138086232856},
            {9, 16, 11.634069043116428},
            {10, 4, 11.634069043116428},
            {10, 11, -23.268138086232856},
            {10, 22, 2.3268138086232857},
            {11, 0, 0.6716932893813962},
            {11, 3, -10.075399340720942},
            {11, 10, 10.075399340720942},
            {11, 21, -0.6716932893813962},
            {12, 1, 4.030159736288377},
            {12, 6, -13.433865787627923},
            {12, 15, 4.030159736288377}};
        break;
    }
    }

    return spherical_trafo;
}