#include "lible_types.h"

using namespace Lible;
using namespace std;

template <typename T>
Vec2D<T>::Vec2D(const size_t &dim1, const size_t &dim2, const T &val)
{
    dims = {dim1, dim2};
    data.resize(dim1 * dim2, val);
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
T &Vec2D<T>::operator()(const size_t &i, const size_t &j)
{
    return data[i * dims[0] + j];
}

template <typename T>
T *Vec2D<T>::getData() const
{
    return data.data();
}