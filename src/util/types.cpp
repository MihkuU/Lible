#include <lible/types.h>

using namespace lible;
using namespace std;

template <typename T>
Vec2D<T>::Vec2D(const size_t &dim1, const size_t &dim2, const T &val)
{
    dims = {dim1, dim2};
    data.resize(dim1 * dim2, val);
}

template <typename T>
T &Vec2D<T>::operator()(const size_t &i, const size_t &j)
{
    return data[i * dims[0] + j];
}

template <typename T>
size_t Vec2D<T>::getDim(const size_t &i) const
{
    return dims.at(i);
}

template <typename T>
size_t Vec2D<T>::getSize() const
{
    return data.size();
}

template <typename T>
const T *Vec2D<T>::getData() const
{
    return data.data();
}

template <typename T>
Vec4D<T>::Vec4D(const size_t &dim1, const size_t &dim2, const size_t &dim3, const size_t &dim4,
                const T &val)
{
    dims = {dim1, dim2, dim3, dim4};
    data.resize(dim1 * dim2 * dim3 * dim4, val);
}

template <typename T>
T &Vec4D<T>::operator()(const size_t &i, const size_t &j, const size_t &k, const size_t &l)
{
    return data[i * dims[0] * dims[1] * dims[2] + j * dims[0] * dims[1] + k * dims[0] + l];
}

template <typename T>
size_t Vec4D<T>::getDim(const size_t &i) const
{
    return dims.at(i);
}

template <typename T>
size_t Vec4D<T>::getSize() const
{
    return data.size();
}

template <typename T>
const T *Vec4D<T>::getData() const
{
    return data.data();
}

namespace lible
{
    template class Vec2D<double>;
    template class Vec4D<double>;
}