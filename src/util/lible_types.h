#pragma once

/*
 *
 *
 */

#include <array>
#include <vector>

namespace Lible
{
    template <typename T>
    class Vec2D
    {
    public:
        Vec2D(const size_t &dim1, const size_t &dim2, const T &val);        

        T &operator()(const size_t &i, const size_t &j);

        size_t getDim(const size_t &i) const;

        size_t getSize() const;
        
        T *getData() const;

    private:
        std::array<size_t, 2> dims;
        std::vector<T> data;
    };

    // template<typename T>
    // class Vec3d;

    // template<typename T>
    // class Vec4d;

    /*
     * Add reshaping utilities
     */
}