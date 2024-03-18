#pragma once

#include <array>
#include <vector>

namespace lible
{
    template <typename T>
    class Vec2D
    {
    public:
        Vec2D()
        {
            dims = {0, 0};
            data.resize(0);
        }

        Vec2D(const size_t &dim, const T &val)
        {
            dims = {dim, dim};
            data.resize(dim * dim, val);
        }

        Vec2D(const size_t &dim1, const size_t &dim2, const T &val)
        {
            dims = {dim1, dim2};
            data.resize(dim1 * dim2, val);
        }

        T &operator()(const size_t &i, const size_t &j)
        {
            return data[i * dims[1] + j];
        }

        const T &operator()(const size_t &i, const size_t &j) const
        {
            return data[i * dims[1] + j];
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
        std::array<size_t, 2> dims;
        std::vector<T> data;
    };
}