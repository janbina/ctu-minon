# ctu-minon

Semetral work for MI-NON course (Nonlinear Continuous Optimization and Numerical Methods) held at [FIT CTU](https://fit.cvut.cz/en) in 2019/2020.
It implements two algorithms for the numerical solution of particular systems of linear equations, [gradient descent](https://en.wikipedia.org/wiki/Gradient_descent)
and [conjugate gradient method](https://en.wikipedia.org/wiki/Conjugate_gradient_method). There are two ways of storing the sparse matrix
provided as input: as a full matrix, storing zeros, or [compressed sparse row](https://en.wikipedia.org/wiki/Sparse_matrix#Compressed_sparse_row_(CSR,_CRS_or_Yale_format)).
