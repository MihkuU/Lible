#pragma once

#include <array>
#include <vector>

namespace lible
{
    template <typename T>
    class Vec4D
    {
    public:
        Vec4D()
        {
            dims = {0, 0, 0, 0};
            data.resize(0);
        }

        Vec4D(const std::size_t &dim, const T &val)
        {
            dims = {dim, dim, dim, dim};
            data.resize(dim * dim * dim * dim, val);
        }

        Vec4D(const std::size_t &dim1, const std::size_t &dim2, const std::size_t &dim3,
              const std::size_t &dim4, const T &val)
        {
            dims = {dim1, dim2, dim3, dim4};
            data.resize(dim1 * dim2 * dim3 * dim4, val);
        }

        T &operator()(const std::size_t &i, const std::size_t &j, const std::size_t &k,
                      const std::size_t &l)
        {
            return data[i * dims[1] * dims[2] * dims[3] + j * dims[2] * dims[3] + k * dims[3] + l];
        }

        const T &operator()(const std::size_t &i, const std::size_t &j, const std::size_t &k,
                            const std::size_t &l) const
        {
            return data[i * dims[1] * dims[2] * dims[3] + j * dims[2] * dims[3] + k * dims[3] + l];
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
        std::array<std::size_t, 4> dims;
        std::vector<T> data;
    };
}