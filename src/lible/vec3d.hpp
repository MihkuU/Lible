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

        Vec3D(const std::size_t &dim, const T &val)
        {
            dims = {dim, dim, dim};
            data.resize(dim * dim * dim, val);
        }

        Vec3D(const std::size_t &dim1, const std::size_t &dim2, const std::size_t &dim3,
              const T &val)
        {
            dims = {dim1, dim2, dim3};
            data.resize(dim1 * dim2 * dim3, val);
        }

        T &operator()(const std::size_t &i, const std::size_t &j, const std::size_t &k)
        {
            return data[i * dims[1] * dims[2] + j * dims[2] + k];
        }

        const T &operator()(const std::size_t &i, const std::size_t &j, const std::size_t &k) const
        {
            return data[i * dims[1] * dims[2] + j * dims[2] + k];
        }

        std::size_t getDim(const std::size_t &i) const
        {
            return dims.at(i);
        }

        std::size_t getSize() const
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
        std::array<std::size_t, 3> dims;
        std::vector<T> data;
    };
}