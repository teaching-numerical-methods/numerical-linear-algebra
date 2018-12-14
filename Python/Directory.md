# Python demonstration codes

These demonstration codes focus on a number of core concepts in an advanced undergraduate / beginning graduate Numerial Linear Algebra course.  Each demonstration may be found inside a separate folder, summarized briefly below:

* matvec (Python v3.5 or higher): performs various formulations for matrix-vector multiplication (row-oriented, column-oriented, native [LAPACK]).  The matrices increase in size, and both the runtimes and approximate error are output to the screen.  The main script is named "driver.py".

* matmat (Python v3.5 or higher): performs various formulations for matrix-matrix multiplication (various loop orderings, native [LAPACK]).  The matrices increase in size, and both the runtimes and approximate error are output to the screen.  The main script is named "driver.py".

* cholesky (Python v3.5 or higher): performs two different formulations for the Cholesky factorization (outer-product vs inner-product).  The factorizations are used within column-oriented forward/backward substitution routines to solve linear systems of increasing size; both the runtimes and solution error are output to the screen.  The main script is named "driver.py".

* sp_lu (Python v3.5 or higher): demonstrates the use of various reordering algorithms (symamd, symrcm, and none) to reduce fill-in when computing the LU factorization of sparse matrices.  The main script is named "driver.m".  This uses the same sparse symmetric matrices as in the "sp_chol" demonstration codes, but with LU instead of Cholesky factorization.  The main script is named "driver.py".

* eigenvalues (Python v3.5 or higher): demonstrations related to Francis's algorithm for eigenvalue computation.  Two "driver" scripts are included:

  - "driver1.py" converts both symmetric and non-symmetric matrices to upper-Hessenberg form via the unitary similarity transformation ``H = Q'*A*Q``, where A is a general square matrix, Q is a unitary transformation matrix, and H is the upper-Hessenberg result.

  - "driver2.py" shows what happens with successive Francis iterations on both symmetric and non-symmetric matrices, using both the Rayleigh quotient shift and the Wilkinson shift.
