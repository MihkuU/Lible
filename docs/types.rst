.. _Types:

Types
=====

Lible provides some generic types for a more convenient handling of data. Some of the types shown 
here merely type aliases of the common standard template library ones. In addition, a generic 
multidimensional array with a contiguous storage is provided. 

Type aliases
------------

The following types are based on the standard template library:

.. cpp:type:: template<typename T, size_t d1, size_t d2> \
               arr2d = std::array<std::array<T, d2>, d1>

    Two-dimensional static array of type ``T`` and dimensions ``d1`` and ``d2``.
    
    .. code-block:: c++

        #include <lible/types.hpp>
        ...        
        // Initialization of a 2D array of integers, d1 = 2 and d2 = 3.
        std::array<int, 2, 3> data_2d{{{1, 2, 3}, {4, 5, 6}}}; // Note the extra pair of braces!

        // Data access
        int three = data_2d[0][2];
        int five = data_2d[1][1];

.. cpp:type:: template<typename T, size_t d1, size_t d2, size_t d3> \
              arr3d = std::array<std::array<std::array<T, d3>, d2>, d1>

    Three-dimensional static array of type ``T`` and dimensions ``d1``, ``d2`` and ``d3``.

    .. code-block:: c++

        #include <lible/types.hpp>
        ...
        // Initialization of a 3D array of doubles, d1 = d2 = d3 = 2.
        lible::arr3d<double, 2, 2, 2> data_3d{{
            {{
                {1.0, 2.0},
                {3.0, 4.0},
            }},
            {{
                {5.0, 6.0},
                {7.0, 8.0}
            }}
            }};

        // Data access 
        double two = data_3d[0][0][1]
        double seven = data_3d[1][1][0];

    The curly braces usage here is cumbersome, but for common practical usage, it is sufficient to 
    know just the data access. 


.. cpp:type:: template<typename T, size_t d1, size_t d2, size_t d3, size_t d4> \
              arr4d = std::array<std::array<std::array<std::array<T, d4>, d3>, d2>, d1>

    Four-dimensional static array of type ``T`` and dimensions ``d1``, ``d2``, ``d3`` and ``d4``.

.. cpp:type:: template<typename T> \ 
              vecvec = std::vector<std::vector<T>>

    Nested two-dimensional vector. 

The following multidimensional array type aliases are defined:    

.. cpp:type:: vec2d = VectorMD<double, 2>;

    Two-dimensional dynamical array of doubles. 

.. cpp:type:: vec3d = VectorMD<double, 3>;

    Three-dimensional dynamical array of doubles.

.. cpp:type:: vec4d = VectorMD<double, 4>;

    Four-dimensional dynamical array of doubles.

.. cpp:type:: vec3i = VectorMD<int, 3>;

    Three-dimensional dynamical array of integers.

Multidimensional array
----------------------

.. cpp:struct:: template<typename T> \ 
    Fill 

    Helper struct to initialize the multidimensional array with desired values. 

    .. cpp:function:: Fill(T val)

        Constructor that initializes the filling value.

    .. cpp:member:: T val_ 

        Value used for filling the multidimensional array.
    
.. cpp:class:: template<typename T, std::size_t N> \
    VectorMD 

    Class for a multidimensional array of type ``T`` and ``N`` dimensions. 

    .. cpp:function:: template<std::integral... Args, std::enable_if_t<sizeof...(Args) == N>> \ 
        VectorMD(Args... dims)

        Constructor for initializing a multidimensional array with given dimensions ``Args... dims``.
        The values are unitialized.

        .. code-block:: c++

            // Initializes a 2x2 dynamical array of doubles. `vec2d` is a type alias of VectorMD.
            lible::vec2d vec(2, 2); 

    .. cpp:function:: template<typename U, std::integral... Args, std::enable_if_t<sizeof...(Args) == N>> \ 
        VectorMD(const Fill<U> init_val, Args... dims)

        Constructor for initializing a multidimensional array with given dimensions ``Args... dims``.
        The values are initialized to value ``Fill<U> init_val``.

        .. code-block:: c++

            // Initializes a 2x2 dynamical array of doubles with value 1.0.
            lible::vec2d vec(lible::Fill(1), 2, 2);

    .. cpp:function:: VectorMD(const std::size_t dim)

        Constructor for initializing a multidimensional array with all equal dimensions, ``dim``. 

        .. code-block:: c++ 

            // Initializes a 3x3 dynamical array of doubles.
            lible::vec2d vec(3);
          
    .. cpp:function:: template<typename U> \
        VectorMD(const Fill<U> init_val, const std::size_t dim)

        Constructor for initializing a multidimensional array with all equal dimensions, ``dim``.
        The values are initialized to value ``Fill<U> init_val``.
        
        .. code-block:: c++ 

            // Initializes a 3x3 dynamical array of doubles. All values are set to 1.5
            lible::vec2d vec(lible::Fill(1.5) 3); 

    .. cpp:function:: template<std::size_t idim> \ 
        std::size_t dim() const 

        Returns the value of the dimension specified by ``idim``.

        .. code-block:: c++
            
            lible::vec3d vec(2, 3, 4); // Initializes a 2x3x4 3D array.
            size_t dim2 = vec.dim<2>(); // Returns 4.
            size_t dim3 = vec.dim<3>(); // Won't compile!

    .. cpp:function:: std::size_t size() const 

        Returns the total number of elements in the array. 

    .. cpp:function:: std::vector<T>::iterator begin() 

        Begin iterator.

    .. cpp:function:: std::vector<T>::iterator end() 

        End iterator. 

        .. code-block:: c++

            lible::vec2d vec(lible::Fill(0.5), 3, 2);
            for (double &val : vec)
                val *= 0.5; // Scale values with 0.5

            for (double val : vec) // begin and end iterator enables range-based for-loops
                printf("val = %f\n", val);

    .. cpp:function:: T *memptr()

        Returns a raw pointer to the beginning of the underlying data.

    .. cpp:function:: const T *memptr() const 

        Returns a constant raw pointer to the beginning of the underlying data.

    .. cpp:function:: template<std::integral... Args, typename = std::enable_if_t<sizeof...(Args) == N>> \
        void resize(Args... dims)

        Resizes the array with given dimensions, ``Args... dims``.

        .. code-block:: c++

            lible::vec4d vec(2, 3, 2, 3); // A 4D array with dimensions 2x3x2x3.
            printf("(%zu, %zu)", vec.dim<0>(), vec.dim<1>()); // Prints (2, 3)
            vec.resize(3, 2, 3, 2); // Resizes the dimensions to 3x2x3x2.
            printf("(%zu, %zu)", vec.dim<0>(), vec.dim<1>()); // Prints (3, 2)

    .. cpp:function:: void resize(const std::size_t dim)

        Resizes the array with all equal dimensions, ``dim``.

    .. cpp:function:: void set(T val)

        Sets the values to ``val``.

    .. cpp:function:: std::array<std::size_t, N> getDimensions() const

        Returns the dimensions. 

    .. cpp:function:: std::array<std::size_t, N - 1> getBlockSizes() const 

        Returns the block sizes. The blocks are used for calculating the row-major in the contiguous
        data from a given index tuple. For example, in case of a 3D array with dimensions 
        :math:`d_1`, :math:`d_2` and :math:`d_3`, the block sizes are given by :math:`b_1 = d_2d_3`
        and :math:`b_2 = d_3`. The 1D index is then :math:`ijk = ib_1 + jb_2 + k`.

    .. cpp:function:: std::vector<T> getData() const 

        Returns the underlying data as a vector.

    .. cpp:function:: template<std::integral... Args> \ 
        T& operator()(Args... args)

        Returns a reference to the data element at given indices. Does bounds checking in debug 
        mode.

        .. code-block:: c++

            lible::VectorMD<int, 2> vec(4, 4);
            int sum{};
            for (int i = 0, ij = 0; i < 4; i++)
                for (int j = 0; j < 4; j++, ij++)
                {
                    vec(i, j) = ij; // Assign value
                    sum += vec(i, j); // Access assigned value
                }

    .. cpp:function:: template<std::integral... Args> \ 
        const T& operator()(Args... args) const 

        Returns a constant reference to the data element at given indices. Does a bounds checking 
        in debug mode. 

    .. cpp:function:: T &operator[](const std::size_t idx)

        Returns a reference to the continuous data at index ``idx``.

    .. cpp:function:: const T &operator[](const std::size_t idx) const

        Returns a constant reference tothe continuous data at index ``idx``.

    .. cpp:function:: VectorMD operator+(const VectorMD &other) const 

        Adds the elements of two multidimensional arrays. 

        .. code-block:: c++

            lible::vec3d vec1(lible::Fill(0.5), 2, 2, 2);
            lible::vec3d vec2(lible::Fill(0.5), 2, 2, 2);
            lible::vec3d vec_sum = vec1 + vec2; // Sums the elements. Result: 8.0

    .. cpp:function:: VectorMD operator-(const VectorMD &other) const 

        Subtracts the elements of a multidimensional array from another.

    .. cpp:function:: VectorMD operator*(T val) const 

        Scales the elements and returns a copy.

    .. cpp:function:: VectorMD &operator*=(T val)

        Scales the elements in place. 

    .. cpp:function:: friend VectorMD operator*(T val, const VectorMD &other)

        Scales the elements with ``val``.

    .. cpp:function:: friend bool operator==(const VectorMD &lhs, const VectorMD &rhs)

        Comparison operator for two multidimensional arrays.