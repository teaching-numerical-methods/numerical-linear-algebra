# Matlab demonstration codes

These demonstration codes focus on a number of core concepts in an advanced undergraduate / beginning graduate Numerial Linear Algebra course.  Each demonstration may be found inside a separate folder, summarized briefly below:

* matvec: performs various formulations for matrix-vector multiplication (row-oriented, column-oriented, native [LAPACK]).  The matrices increase in size, and both the runtimes and approximate error are output to the screen.  The main script is named "driver.m".

* matmat: performs various formulations for matrix-matrix multiplication (various loop orderings, native [LAPACK]).  The matrices increase in size, and both the runtimes and approximate error are output to the screen.  The main script is named "driver.m".

* cholesky: performs two different formulations for the Cholesky factorization (outer-product vs inner-product).  The factorizations are used within column-oriented forward/backward substitution routines to solve linear systems of increasing size; both the runtimes and solution error are output to the screen.  The main script is named "driver.m".

* sp_chol: demonstrates the use of various reordering algorithms (symamd, symrcm, and none) to reduce fill-in when computing the Cholesky factorization of sparse symmetric matrices.  The main script is named "driver.m".

* sp_lu: demonstrates the use of various reordering algorithms (symamd, symrcm, and none) to reduce fill-in when computing the LU factorization of sparse matrices.  The main script is named "driver.m".  This uses the same sparse symmetric matrices as in the "sp_chol" demonstration codes, but with LU instead of Cholesky factorization.  The main script is named "driver.m".

* GE_roundoff: demonstrates the challenges with using floating-point arithmetic for matrix computations.  All use the "roundsd" function for limiting computations to fixed precision.  Three scripts are included:

  - "example1.m" performs step-by-step LU factorization of a 4x4 Hilbert matrix, rounding each entry to only 3 significant digits.  

  - "example2.m" performs step-by-step Gaussian elimination of a 3x3 linear system using 4-digit precision, both with and without row-swapping.

  - "example3.m" performs step-by-step Gaussian elimination of a 3x3 linear system using 4-digit precision without row swapping, but with row scaling to "fix" small pivot values.

* eigenvalues: demonstrations related to Francis's algorithm for eigenvalue computation.  Two "driver" scripts are included:

  - "driver1.m" converts both symmetric and non-symmetric matrices to upper-Hessenberg form via the unitary similarity transformation ``H = Q'*A*Q``, where A is a general square matrix, Q is a unitary transformation matrix, and H is the upper-Hessenberg result.

  - "driver2.m" shows what happens with successive Francis iterations on both symmetric and non-symmetric matrices, using both the Rayleigh quotient shift and the Wilkinson shift.

  
