#pragma once

#include <array>
#include <concepts>
#include <format>
#include <numeric>
#include <stdexcept>
#include <vector>

#ifndef NDEBUG
#define NODEBUG false
#else
#define NODEBUG true
#endif

#ifndef DEBUG
#define DEBUG !NODEBUG
#endif

namespace lible
{
    /** Object for initializing VectorMD to a given value. */
    template <typename T>
    struct Fill
    {
        Fill(T val) : val(val)
        {
        }

        const T val;
    };

    /**
     * Class for a multidimensional vector with dimension N and type T.
     * The data is stored contiguously in a vector.
     */
    template <typename T, std::size_t N>
    class VectorMD
    {
    public:
        // Constructors

        /** */
        VectorMD()
        {
        }

        /** */
        template <std::integral... Args,
                  typename = typename std::enable_if<sizeof...(Args) == N>::type>
        VectorMD(Args... dims)
        {
            static_assert(sizeof...(dims) == N);

            resize(dims...);
        }

        /** */
        template <typename U, std::integral... Args,
                  typename = typename std::enable_if<sizeof...(Args) == N>::type>
        VectorMD(Fill<U> init_val, Args... dims)
            : VectorMD(dims...)
        {
            std::fill(data.begin(), data.end(), init_val.val);
        }

        /** */
        VectorMD(std::size_t dim)
        {
            resize(dim);
        }

        /** */
        template <typename U>
        VectorMD(Fill<U> init_val, std::size_t dim)
            : VectorMD(dim)
        {
            std::fill(data.begin(), data.end(), init_val.val);
        }

        // Member functions

        /** */
        template<std::size_t idim>
        std::size_t dim() const
        {
            static_assert(idim < N);

            return dimensions[idim];
        }

        /** */
        std::size_t size() const
        {
            return data.size();
        }

        /** */
        const std::vector<T>::const_iterator begin() const 
        {
            return data.begin();
        }

        /** */
        const std::vector<T>::const_iterator end() const 
        {
            return data.end();
        }

        /** */
        T* memptr()
        {
            return data.data();
        }

        /** */
        const T* memptr() const 
        {
            return data.data();
        }

        /** */
        template <std::integral... Args,
                  typename = typename std::enable_if<sizeof...(Args) == N>::type>
        void resize(Args... dims)
        {
            static_assert(sizeof...(dims) == N);

            std::size_t i = 0;
            (void(dimensions[i++] = dims), ...);

            for (std::size_t i = 1; i < N; i++)
                block_sizes[i - 1] = std::accumulate(dimensions.begin() + i, dimensions.end(), 1,
                                                     std::multiplies<T>());

            std::size_t size = std::accumulate(dimensions.begin(), dimensions.end(), 1,
                                               std::multiplies<T>());
            data.resize(size);
        }

        /** */
        void resize(std::size_t dim)
        {
            for (std::size_t i = 0; i < N; i++)
                dimensions[i] = dim;

            for (std::size_t i = 1; i < N; i++)
                block_sizes[i - 1] = std::accumulate(dimensions.begin() + i, dimensions.end(), 1,
                                                     std::multiplies<T>());

            std::size_t size = std::accumulate(dimensions.begin(), dimensions.end(), 1,
                                               std::multiplies<T>());
            data.resize(size);
        }

        /** */
        void set(T val)
        {
            std::fill(data.begin(), data.end(), val);
        }

        // Getters

        std::array<std::size_t, N> getDimensions() const
        {
            return dimensions;
        }

        std::array<std::size_t, N - 1> getBlockSizes() const
        {
            return block_sizes;
        }

        std::vector<T> getData() const
        {
            return data;
        }

        // Operators

        /** */
        template <std::integral... Args>
        T& operator()(Args... idxs)
        {
            static_assert(sizeof...(idxs) == N);

            std::size_t idx = calcIdx(idxs...);

            return data[idx];
        }

        /** */
        template<std::integral... Args>
        const T& operator()(Args... idxs) const 
        {
            static_assert(sizeof...(idxs) == N);

            std::size_t idx = calcIdx(idxs...);

            return data[idx];
        }

        /** */
        T& operator[](std::size_t idx)
        {
            return data[idx];
        }

        /** */
        const T& operator[](std::size_t idx) const 
        {
            return data[idx];
        }

        /** */
        VectorMD operator+(const VectorMD &other) const
        {
#ifdef DEBUG
            if (this->dimensions != other.dimensions)
                throw std::out_of_range("VectorMD::operator+: dimensions don't match");
#endif

            VectorMD out = *this;
            for (std::size_t i = 0; i < data.size(); i++)
                out.data[i] += other.data[i];

            return out;
        }

        /** */
        VectorMD operator-(const VectorMD &other) const
        {
#ifdef DEBUG
            if (this->dimensions != other.dimensions)
                throw std::out_of_range("VectorMD::operator-: dimensions don't match");
#endif

            VectorMD out = *this;
            for (std::size_t i = 0; i < data.size(); i++)
                out.data[i] -= other.data[i];

            return out;
        }

        /** */
        VectorMD &operator+=(const VectorMD &other)
        {
#ifdef DEBUG
            if (this->dimensions != other.dimensions)
                throw std::out_of_range("VectorMD::operator+=: dimensions don't match");
#endif

            for (std::size_t i = 0; i < data.size(); i++)
                this->data[i] += other.data[i];

            return *this;
        }

        /** */
        VectorMD& operator-=(const VectorMD &other)
        {
#ifdef DEBUG
            if (this->dimensions != other.dimensions)
                throw std::out_of_range("VectorMD::operator-=: dimensions don't match");
#endif
            for (std::size_t i = 0; i < data.size(); i++)
                this->data[i] -= other.data[i];

            return *this;
        }

        /** */
        VectorMD operator*(T val) const
        {
            VectorMD out = *this;
            for (std::size_t i = 0; i < data.size(); i++)
                out.data[i] *= val;

            return out;
        }

        /** */
        VectorMD& operator*=(T val)
        {
            for (std::size_t i = 0; i < data.size(); i++)
                this->data[i] *= val;

            return *this;
        }

        /** */
        friend VectorMD operator*(T val, const VectorMD &other)
        {
            return other * val;
        }

        /** */
        friend bool operator==(const VectorMD &lhs, const VectorMD &rhs)
        {
            return ((lhs.getDimensions() == rhs.getDimensions()) &&
                    (lhs.getBlockSizes() == rhs.getBlockSizes()) &&
                    (lhs.getData() == rhs.getData()));
        }

    private:
        std::array<std::size_t, N> dimensions{};      /** */
        std::array<std::size_t, N - 1> block_sizes{}; /** */
        std::vector<T> data{};                        /** */

        /** */
        template <typename U>
        std::size_t calcIdx(U idx) const
        {
#ifdef DEBUG
            // It is assumed that calcIdx() is used by operator() only.
            constexpr std::size_t idim = N - 1;
            if (idx > (dimensions[idim] - 1))
                throw std::out_of_range(std::format("VectorMD::operator(): index value {} along "
                                                    "dimension {} is out of bounds",
                                                    idx, idim));
#endif

            return idx;
        }

        /** */
        template <typename U, typename... Idxs>
        std::size_t calcIdx(U idx, Idxs... idxs) const
        {
            constexpr std::size_t idim = N - sizeof...(idxs) - 1;
#ifdef DEBUG
            // It is assumed that calcIdx() is used by operator() only.
            if (idx > (dimensions[idim] - 1))
                throw std::out_of_range(std::format("VectorMD::operator(): index value {} along "
                                                    "dimension {} is out of bounds",
                                                    idx, idim));
#endif

            return idx * block_sizes[idim] + calcIdx(idxs...);
        }
    };
}