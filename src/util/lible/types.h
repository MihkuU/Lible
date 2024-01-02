#pragma once

/*
 *
 *
 */

#include <array>
#include <vector>

namespace lible
{
    template <typename T>
    class Vec2D;
    template <typename T>
    class Vec3D;
    template <typename T>
    class Vec4D;

    using vec2d = Vec2D<double>;
    using vec3d = Vec3D<double>;
    using vec4d = Vec4D<double>;

    template <typename T>
    class Vec2D
    {
    public:
        Vec2D(const size_t &dim1, const size_t &dim2, const T &val);

        T &operator()(const size_t &i, const size_t &j);

        size_t getDim(const size_t &i) const;

        size_t getSize() const;

        const T *getData() const;

    private:
        std::array<size_t, 2> dims;
        std::vector<T> data;
    };

    template <typename T>
    class Vec4D
    {
    public:
        Vec4D(const size_t &dim1, const size_t &dim2, const size_t &dim3, const size_t &dim4,
              const T &val);

        T &operator()(const size_t &i, const size_t &j, const size_t &k, const size_t &l);

        size_t getDim(const size_t &i) const;

        size_t getSize() const;

        const T *getData() const;

    private:
        std::array<size_t, 4> dims;
        std::vector<T> data;
    };
    /*
     * Add reshaping utilities
     */
}