#pragma once

#include <array>
#include <vector>

namespace lible
{
    template <typename T>
    class Vec3D
    {
    public:
        Vec3D()
        {
            dims = {0, 0, 0};
            data.resize(0);
        }

        Vec3D(const size_t &dim, const T &val)
        {
            dims = {dim, dim, dim};
            data.resize(dim * dim * dim, val);
        }

        Vec3D(const size_t &dim1, const size_t &dim2, const size_t &dim3, const T &val)
        {
            dims = {dim1, dim2, dim3};
            data.resize(dim1 * dim2 * dim3, val);
        }

        T &operator()(const size_t &i, const size_t &j, const size_t &k)
        {
            return data[i * dims[1] * dims[2] + j * dims[2] + k];
        }

        const T &operator()(const size_t &i, const size_t &j, const size_t &k) const
        {
            return data[i * dims[1] * dims[2] + j * dims[2] + k];
        }

        size_t getDim(const size_t &i) const
        {
            return dims.at(i);
        }

        size_t getSize() const
        {
            return data.size();
        }

        T *getData()
        {
            return data.data();
        }

        const T *getData() const
        {
            return data.data();
        }

        void set(const T &val)
        {
            std::fill(data.begin(), data.end(), val);
        }

    private:
        std::array<size_t, 3> dims;
        std::vector<T> data;
    };
}