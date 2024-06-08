#include <lible/types.hpp>

using namespace lible;
using namespace std;

// template <typename T>
// Vec2D<T>::Vec2D()
// {
//     dims = {0, 0};
//     data.resize(0);
// }

// template <typename T>
// Vec2D<T>::Vec2D(const size_t &dim1, const size_t &dim2, const T &val)
// {
//     dims = {dim1, dim2};
//     data.resize(dim1 * dim2, val);
// }

// template <typename T>
// T &Vec2D<T>::operator()(const size_t &i, const size_t &j)
// {
//     return data[i * dims[1] + j];
// }

// template <typename T>
// T Vec2D<T>::operator()(const size_t &i, const size_t &j) const
// {
//     return data[i * dims[1] + j];
// }

// template <typename T>
// size_t Vec2D<T>::getDim(const size_t &i) const
// {
//     return dims.at(i);
// }

// template <typename T>
// size_t Vec2D<T>::getSize() const
// {
//     return data.size();
// }

// template <typename T>
// const T *Vec2D<T>::getData() const
// {
//     return data.data();
// }

// template <typename T>
// Vec4D<T>::Vec4D()
// {
//     dims = {0, 0, 0, 0};
//     data.resize(0);    
// }

// template <typename T>
// Vec4D<T>::Vec4D(const size_t &dim1, const size_t &dim2, const size_t &dim3, const size_t &dim4,
//                 const T &val)
// {
//     dims = {dim1, dim2, dim3, dim4};
//     data.resize(dim1 * dim2 * dim3 * dim4, val);
// }

// template <typename T>
// T &Vec4D<T>::operator()(const size_t &i, const size_t &j, const size_t &k, const size_t &l)
// {
//     return data[i * dims[1] * dims[2] * dims[3] + j * dims[2] * dims[3] + k * dims[3] + l];
// }

// template <typename T>
// T Vec4D<T>::operator()(const size_t &i, const size_t &j, const size_t &k, const size_t &l) const
// {
//     return data[i * dims[1] * dims[2] * dims[3] + j * dims[2] * dims[3] + k * dims[3] + l];
// }

// template <typename T>
// size_t Vec4D<T>::getDim(const size_t &i) const
// {
//     return dims.at(i);
// }

// template <typename T>
// size_t Vec4D<T>::getSize() const
// {
//     return data.size();
// }

// template <typename T>
// const T *Vec4D<T>::getData() const
// {
//     return data.data();
// }

// namespace lible
// {
//     template class Vec2D<double>;
//     template class Vec2D<int>;
//     template class Vec4D<double>;
// }