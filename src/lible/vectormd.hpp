#pragma once

#include <array>
#include <concepts>
#include <format>
#include <numeric>
#include <stdexcept>
#include <vector>

namespace lible
{
    /// Object for initializing a VectorMD object to a given value.
    template <typename T>
    struct Fill
    {
        explicit Fill(T val) : val_(val)
        {
        }

        /// Value used for filling `VectorMD` values.
        const T val_;
    };

    /// Class for a multidimensional vector with dimension N and type T.
    /// The data is stored contiguously in a vector in a row-major format.
    template <typename T, std::size_t N>
    class VectorMD
    {
    public:
        // Constructors

        /// Default constructor.
        VectorMD() = default;

        /// Constructor for initializing a multidimensional array with given dimensions
        /// (dim1, dim2, ...). The length of the list must equal the dimension of the array.
        template <std::integral... Args,
            typename = std::enable_if_t<sizeof...(Args) == N>>
        explicit VectorMD(Args... dims)
        {
            static_assert(sizeof...(dims) == N);

            resize(dims...);
        }

        /// Constructor for initializing a multidimensional array with a given initial value,
        /// Fill(init_val), and dimensions (dim1, dim2, ...). The length of the list must
        /// equal the dimension of the array.
        template <typename U, std::integral... Args,
            typename = std::enable_if_t<sizeof...(Args) == N>>
        explicit VectorMD(const Fill<U> init_val, Args... dims)
            : VectorMD(dims...)
        {
            std::fill(data_.begin(), data_.end(), init_val.val_);
        }

        /// Constructor for initializing a multidimensional array with all equal dimensions.
        explicit VectorMD(const std::size_t dim)
        {
            resize(dim);
        }

        /// Constructor for initializing a multidimensional array with a given initial value,
        /// Fill(init_val), and all equal dimensions.
        template <typename U>
        VectorMD(const Fill<U> init_val, const std::size_t dim)
            : VectorMD(dim)
        {
            std::fill(data_.begin(), data_.end(), init_val.val_);
        }

        // Member functions

        /// Returns the value of the `idim` dimension.
        template <std::size_t idim>
        std::size_t dim() const
        {
            static_assert(idim < N);

            return dimensions_[idim];
        }

        /// Returns the total number of elements in the array.
        std::size_t size() const
        {
            return data_.size();
        }

        /// Constant `begin` iterator.
        std::vector<T>::iterator begin()
        {
            return data_.begin();
        }

        /// Constant `end` iterator.
        std::vector<T>::iterator end()
        {
            return data_.end();
        }

        /// Returns a raw pointer to the underlying data (beginning).
        T *memptr()
        {
            return data_.data();
        }

        /// Returns a constant raw pointer to the underlying data (beginning).
        const T *memptr() const
        {
            return data_.data();
        }

        /// Resizes the dimensions and underlying data.
        template <std::integral... Args,
            typename = std::enable_if_t<sizeof...(Args) == N>>
        void resize(Args... dims)
        {
            static_assert(sizeof...(dims) == N);

            std::size_t i = 0;
            (void(dimensions_[i++] = dims), ...);

            for (std::size_t j = 1; j < N; j++)
                block_sizes_[j - 1] = std::accumulate(dimensions_.begin() + j, dimensions_.end(), 1,
                                                      std::multiplies<T>());

            std::size_t size = std::accumulate(dimensions_.begin(), dimensions_.end(), 1,
                                               std::multiplies<T>());
            data_.resize(size);
        }

        /// Resizes the dimensions (all equal) and underlying data.
        void resize(const std::size_t dim)
        {
            for (std::size_t i = 0; i < N; i++)
                dimensions_[i] = dim;

            for (std::size_t i = 1; i < N; i++)
                block_sizes_[i - 1] = std::accumulate(dimensions_.begin() + i, dimensions_.end(), 1,
                                                      std::multiplies<T>());

            std::size_t size = std::accumulate(dimensions_.begin(), dimensions_.end(), 1,
                                               std::multiplies<T>());
            data_.resize(size);
        }

        /// Sets the values of the underlying data.
        void set(T val)
        {
            std::fill(data_.begin(), data_.end(), val);
        }

        // Getters

        /// Returns the dimensions as a static array.
        std::array<std::size_t, N> getDimensions() const
        {
            return dimensions_;
        }

        /// Returns the block sizes.
        std::array<std::size_t, N - 1> getBlockSizes() const
        {
            return block_sizes_;
        }

        /// Returns the underlying data as a vector.
        std::vector<T> getData() const
        {
            return data_;
        }

        // Operators

        /// Returns a reference to the data element at given indices (i, j, ...). Does a bounds
        /// checking in the `DEBUG` mode.
        template <std::integral... Args>
        T &operator()(Args... idxs)
        {
            static_assert(sizeof...(idxs) == N);

            std::size_t idx = calcIdx(idxs...);

            return data_[idx];
        }

        /// Returns a const reference to the data element at given indices (i, j, ...). Does a
        /// bounds checking in `DEBUG` mode.
        template <std::integral... Args>
        const T &operator()(Args... idxs) const
        {
            static_assert(sizeof...(idxs) == N);

            std::size_t idx = calcIdx(idxs...);

            return data_[idx];
        }

        /// Returns a reference to the continuous data at index `idx`.
        T &operator[](const std::size_t idx)
        {
            return data_[idx];
        }

        /// Returns a const reference to the continuous data at index `idx`.
        const T &operator[](const std::size_t idx) const
        {
            return data_[idx];
        }

        /// Adds `other` elements to `this` and returns a copy.
        VectorMD operator+(const VectorMD &other) const
        {
#ifndef NDEBUG
            if (this->dimensions_ != other.dimensions_)
                throw std::out_of_range("VectorMD::operator+: dimensions don't match");
#endif

            VectorMD out = *this;
            for (std::size_t i = 0; i < data_.size(); i++)
                out.data_[i] += other.data_[i];

            return out;
        }

        /// Subtracts `other` elements from `this` and returns a copy.
        VectorMD operator-(const VectorMD &other) const
        {
#ifndef NDEBUG
            if (this->dimensions_ != other.dimensions_)
                throw std::out_of_range("VectorMD::operator-: dimensions don't match");
#endif

            VectorMD out = *this;
            for (std::size_t i = 0; i < data_.size(); i++)
                out.data_[i] -= other.data_[i];

            return out;
        }

        /// Adds the other elements to `this` in place.
        VectorMD &operator+=(const VectorMD &other)
        {
#ifndef NDEBUG
            if (this->dimensions_ != other.dimensions_)
                throw std::out_of_range("VectorMD::operator+=: dimensions don't match");
#endif

            for (std::size_t i = 0; i < data_.size(); i++)
                this->data_[i] += other.data_[i];

            return *this;
        }

        /// Subtracts the other elements from `this` in place.
        VectorMD &operator-=(const VectorMD &other)
        {
#ifndef NDEBUG
            if (this->dimensions_ != other.dimensions_)
                throw std::out_of_range("VectorMD::operator-=: dimensions don't match");
#endif
            for (std::size_t i = 0; i < data_.size(); i++)
                this->data_[i] -= other.data_[i];

            return *this;
        }

        /// Scales the elements and returns a copy.
        VectorMD operator*(T val) const
        {
            VectorMD out = *this;
            for (std::size_t i = 0; i < data_.size(); i++)
                out.data_[i] *= val;

            return out;
        }

        /// Scales the elements in place.
        VectorMD &operator*=(T val)
        {
            for (std::size_t i = 0; i < data_.size(); i++)
                this->data_[i] *= val;

            return *this;
        }

        /// Scales the elements with `val`.
        friend VectorMD operator*(T val, const VectorMD &other)
        {
            return other * val;
        }

        /// equals-operator.
        friend bool operator==(const VectorMD &lhs, const VectorMD &rhs)
        {
            return (lhs.getDimensions() == rhs.getDimensions() &&
                    lhs.getBlockSizes() == rhs.getBlockSizes() &&
                    lhs.getData() == rhs.getData());
        }

    private:
        /// The dimensions, (d1, d2, ...) .
        std::array<std::size_t, N> dimensions_{};
        /// Block sizes for index calculation. For example, for (d1, d2, d3) the block sizes are
        /// (d2 * d3, d3).
        std::array<std::size_t, N - 1> block_sizes_{};
        /// The underlying data in a continuous form.
        std::vector<T> data_{};

        /// Backend function to calculate the index in the contiguous storage.
        template <typename U>
        std::size_t calcIdx(U idx) const
        {
#ifndef NDEBUG
            // It is assumed that calcIdx() is used by operator() only.
            constexpr std::size_t idim = N - 1;
            if (static_cast<std::size_t>(idx) > (dimensions_[idim] - 1))
                throw std::out_of_range(std::format("VectorMD::operator(): index value {} along "
                                                    "dimension {} is out of bounds",
                                                    idx, idim));
#endif

            return idx;
        }

        /// Backend function to calculate the index in the contiguous storage.
        template <typename U, typename... Idxs>
        std::size_t calcIdx(U idx, Idxs... idxs) const
        {
            constexpr std::size_t idim = N - sizeof...(idxs) - 1;
#ifndef NDEBUG
            // It is assumed that calcIdx() is used by operator() only.
            if (static_cast<std::size_t>(idx) > (dimensions_[idim] - 1))
                throw std::out_of_range(std::format("VectorMD::operator(): index value {} along "
                                                    "dimension {} is out of bounds",
                                                    idx, idim));
#endif

            return idx * block_sizes_[idim] + calcIdx(idxs...);
        }
    };
}
