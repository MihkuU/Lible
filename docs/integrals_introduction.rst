

Introduction
============

Lible provides utilities to calculate molecular integrals over contracted Gaussian type atomic 
orbitals (GTAOs). A contracted GTAO can be written as 

.. math::
   G_{\mu} = N_{\mu} \sum^K_{i=1} d_{\mu,i} g_{\mu,i}

where :math:`g_{\mu,i}` denotes a primitive Gaussian basis function. The contraction depth is 
denoted by :math:`K` and the contraction coefficients by :math:`d_{\mu,i}`. The normalization 
coefficient :math:`N_{\mu}` is usually obtained from the overlap integral

.. math::
   (G_{\mu}|G_{\nu}) = (\mu|\nu) = \int  G_{\mu} G_{\nu} d\mathbf{r}
   \; \Rightarrow \; N_{\mu} = 1.0 / \sqrt{(\mu|\mu)}

The primitive Gaussian basis functions are normalized, and depend on the orbital angular momentum
:math:`l` and its projection, :math:`m = -l, -l + 1, \ldots, l`. Omitting the :math:`\mu`-label, 
the primitive Gaussian function can be written as 

.. math::
   g^{lm}_{i}(\mathbf{r}, \mathbf{a}, \mathbf{A}) = N_{l} (a_i) S_{lm} (\mathbf{r}, \mathbf{A})
   e^{-a_i r^2_A}

The contents of the above expression warrant an explanation: 

1. :math:`\mathbf{r} = (x, y, z)` -- an arbitrary point in space.
2. :math:`\mathbf{a} = (a_1,\ldots,a_K)` -- list of Gaussian primitive exponents.
3. :math:`\mathbf{A} = (A_x, A_y, A_z)` -- point in space where the basis function is centered 
   (i.e., a nucleus).
4. :math:`N_{l} (a_i)` -- pure (harmonic) Gaussian primitive normalization coefficient that can be 
   calculated analytically as

.. math::
   N_{l}(a_i) = \sqrt{\frac{(2a_i/\pi)^{3/2}(4a_i)^l}{(2l - 1)!!}}

5. :math:`S_{lm} (\mathbf{A})` -- a real-valued solic harmonic. The explicit expressions for these 
   shall be not given here. However, what is important is that the solid harmonics can be expressed 
   in terms of Cartesian directions through the transformation :math:`t_{lm;ijk}`,

.. math::
   S_{lm} = \sum_{ijk} t_{lm,ijk} x_A^i y_A^j z_A^k 

6. :math:`\mathbf{r}_A = (x - A_x, y - A_y, z - A_z)`

