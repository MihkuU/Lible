#pragma once

#include <array>
#include <cassert>
#include <vector>

namespace lible
{
    namespace IntsUtil // TODO: rename
    {
        inline double doubleFactorial(const int &n)
        {
            /*
             * Definition from DOI:10.1002/9781119019572 in eq. (6.5.10).
             */
            assert(n >= 0 or (n % 2 == -1));

            double double_factorial = 1;
            if (n == 0)            
                return 1;        
            else if (n > 0 and (n % 2 == 0))            
                for (int i = n; i >= 2; i -= 2)
                    double_factorial *= i;            
            else if (n > 0 and (n % 2 == 1))            
                for (int i = n; i >= 1; i -= 2)
                    double_factorial *= i;            
            else            
                for (int i = n + 2; i < 1; i += 2)
                    double_factorial *= 1 / i;
                                
            return double_factorial;
        }
    }
}