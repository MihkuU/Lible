#pragma once

/*
 *
 *
 */

#include <lible/vec2d.hpp>
#include <lible/vec3d.hpp>
#include <lible/vec4d.hpp>

#include <array>
#include <vector>

namespace lible
{
    /*
     * TODO: 
     *   1. Implement each class in a separate header file and then include these files here. 
     *      The entire implementation of the class should be entirely within that header file 
     *      to make the template code fully visible. The below template class forward declarations
     *      should be removed, aliases via the 'using' clause are fine and good. 
     * 
     *   2. Test and optimize array access, especially for the 4D array, where speedups could be gained.
     * 
     *   3. Add utility functions/operators that seem missing. For example, add another constructor 
     *      that allows for simpler initialization if the dimensions are equal, e.g:
     *        Vec3D(const size_t &dim, const T &val); // three dimensions all equal
     * 
     *   4. Generate comments/documentation/webpage. 
     * 
     *   5. Feel free to add stuff and have fun with it. 
     * 
     *   6. If you have suggestions for additional types and stuff that the library can provide, let me know. 
     * 
     *   7. Conversions? This is for the future, there could be conversion functions that allow to convert to/from 
     *      std::vector, Armadillo, Eigen etc. quantities. 
     *
     */

    using vec2d = Vec2D<double>;
    using vec3d = Vec3D<double>;
    using vec4d = Vec4D<double>;

    // template <typename T>
    // class Vec2D
    // {
    // public:
    //     // TODO: make a separate file
    //     Vec2D();
    //     Vec2D(const size_t &dim, const T &val);

    //     Vec2D(const size_t &dim1, const size_t &dim2, const T &val);

    //     T operator()(const size_t &i, const size_t &j) const;
    //     T &operator()(const size_t &i, const size_t &j);

    //     size_t getDim(const size_t &i) const;

    //     size_t getSize() const;

    //     const T *getData() const;

    // private:
    //     std::array<size_t, 2> dims;
    //     std::vector<T> data;
    // };

    // template <typename T>
    // class Vec4D
    // {
    // public:
    //     // TODO: make a separate file
    //     Vec4D();
    //     Vec4D(const size_t &dim, const T &val); 

    //     Vec4D(const size_t &dim1, const size_t &dim2, const size_t &dim3, const size_t &dim4,
    //           const T &val);

    //     T operator()(const size_t &i, const size_t &j, const size_t &k, const size_t &l) const;
    //     T &operator()(const size_t &i, const size_t &j, const size_t &k, const size_t &l);

    //     size_t getDim(const size_t &i) const;

    //     size_t getSize() const;

    //     const T *getData() const;

    // private:
    //     std::array<size_t, 4> dims;
    //     std::vector<T> data;
    // };


    /*
     * Add reshaping utilities
     */
}