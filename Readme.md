# Linear Algebra UDF

A library for linear algebra, developed in AutoIt. This library offers a variety of functions for performing calculations and operations in linear algebra.

## Contents

- [Linear Algebra UDF](#linear-algebra-udf)
	- [Contents](#contents)
	- [Idea](#idea)
	- [Structure](#structure)
	- [Install](#install)
	- [Features](#features)
	- [Documentation](#documentation)
	- [To-Do](#to-do)

## Idea

A UDF for linear algebra in AutoIt.
The aim of this is to be as thematically comprehensive, high-performance and easily accessible as possible.

The widely used software library BLAS/LAPACK serves as the basis for the UDF.
The user should be able to work as intuitively as possible and get by without any major administrative effort.

A particular focus is on extensive functionalities for non-linear adjustment calculations.

## Structure
The UDF is divided into 3 sub-UDFs: `BLAS.au3`, `LAPACK.au3` and `LinearAlgebra.au3`.
The low-level interfaces to the respective BLAS/LAPACK functionalities are implemented in the first two.
`LinearAlgebra.au3`, which is intended as the primary interface for the end user, is built on this basis.
The functions here offer simpler interfaces and access to more complex algorithms such as regressions and adjustment calculations.

In addition to these 3 files, a DLL is also required which implements the BLAS/LAPACK interface.

## Install

* Download the 3 files `BLAS.au3`, `LAPACK.au3` and `LinearAlgebra.au3` (or clone the repository)
* Download a current BLAS/LAPACK DLL:
  * Recommendation: [OpenBLAS](https://github.com/OpenMathLib/OpenBLAS/releases) (theoretically, other BLAS/LAPACK implementations should also work - however, additional adaptations may then be necessary)
  * Download the file `OpenBLAS-x.x.xx-x64.zip` from there and extract the file `libopenblas.dll` into the same folder as the `LinearAlgebra.au3`.
* The sample files in the subfolder `/examples` should now be executable .

## Features
```
---- vector/matrix creation ----
_la_fromArray          - converts a AutoIt array or array define string into a matrix map
_la_fromStruct         - creates a matrix/vector map from a DllStruct as used here in the UDF
_la_createVector       - creates new empty vector
_la_createMatrix       - creates new empty matrix
_la_createIdentity     - create identity matrix/vector
_la_duplicate          - creates an independent copy of a matrix/vector map
_la_fromFile           - reads a matrix or a vector from a file created by _la_toFile()

---- extraction/transforming ----
_la_join               - combines 2 matrices
_la_transpose          - transposes a matrix in-place or out-place and [optional] scaling
_la_ReDim              - changes the shape of a matrix by by changing the number of columns (also matrix <-> vector conversion)
_la_getRow             - extracts a row of a matrix as a vector
_la_getColumn          - extracts a column of a matrix as a vector
_la_getDiag            - extracts the diagonal of a matrix as a vector
_la_getTriangle        - extract upper or lower triangle part of a matrix
_la_VectorToDiag       - creates a diagonal matrix from a vector

---- data output ----
_la_display            - displays a matrix/vector map, similar to _ArrayDisplay
_la_toArray            - converts a matrix/vector map into an AutoIt array
_la_toFile             - write a matrix/vector into a file

---- scalar operations ----
_la_rotate             - applies a plane rotation to coordinate-pairs

---- matrix attributes ----
_la_isPositiveDefinite - checks whether a matrix is positive definite
_la_isSymmetric        - checks whether a matrix is symmetrical
_la_rank               - determines the rank of a matrix
_la_determinant        - calculate the determinant of a matrix
_la_conditionNumber    - determine the condition number of a matrix

---- unary operations ----
_la_inverse            - calculates the inverse of a matrix
_la_pseudoInverse      - calculate the Moore-Penrose pseudo inverse of a matrix
_la_sum                - calculates the sum of the elements of a matrix, vector or parts thereof
_la_asum               - calculate the sum of the absolute(!) values of a matrix/vector
_la_amin               - finds the first element having the minimum absolute(!) value
_la_amax               - finds the first element having the maximum absolute(!) value
_la_norm               - calculate the euclidian norm of a vector
_la_mean               - calculate the mean of a vector or parts of a matrix

---- element wise operations ----
_la_sqrtElements       - calculates the square root of each element of a matrix/vector
_la_squareElements     - calculates the square of each element of a matrix/vector
_la_invElements        - forms the reciprocal (1/x) for each element of the matrix/vector

---- addition subtraction ----
_la_sub                - subtracts a matrix/vector B from matrix/vector A
_la_add                - calculate the sum of a matrix/vector/scalar mA and a matrix/vector/scalar mB

---- multiplication ----
_la_mul                - calculates a multiplication between a matrix/vector/scalar A and a matrix/vector/scalar B
_la_outerproduct       - calculates the outer product ("tensor product") of two vectors
_la_dot                - calculate the "dot product"/"scalar product"/"inner product" of two vectors
_la_scale              - multiplies the elements of a matrix/vector by a scalar value
_la_mulElementWise     - calculates the element-wise ("Hadarmard") product between two matrices/vectors
_la_cross              - calculates the cross product between two 3-element vectors

---- factorization / decomposition ----
_la_LU                 - calculates the LU decomposition of a matrix
_la_QR                 - calculates the QR decomposition of a matrix
_la_SVD                - calculates the singular value decomposition (SVD) of a matrix
_la_cholesky           - calculate the cholesky decomposition of a symmetric, positive definite matrix ( A --> L * Lᵀ or A --> U * Uᵀ )

---- eigenvalues / eigenvectors ----
_la_eigen              - computes for an N-by-N real matrix A, the eigenvalues and the left and/or right eigenvectors.

---- solve linear equation systems ----
_la_solve              - computes the solution to a system of linear equations A * X = B

---- least squares solving ----
_la_lstsq              - solves overdetermined or underdetermined [weighted] linear system

---- regression ----
_la_regression         - calculates an n-dimensional linear or non-linear regression

---- adjustment ----
_la_adjustment         - performs a least-squares adjustment calculation for a system of different [weighted] non-linear equations
_la_adjustment_l1      - performs a adjustment calculation to L1 norm for a system of different [weighted] non-linear equations
_la_adj_addObservation - adds an observation to the adjustment system

---- additional helper functions ----
_la_adj_showResult     - formats the results of _la_adj more clearly and display them in a window
```

## Documentation
The documentation for the individual functions is contained directly in the source code. Each function is provided with a description that explains its parameters and return values. In most cases, a short example is also included here. You will also find detailed explanations in the example files.

## To-Do
* Certain functions (e.g. `_la_add()`, `_la_mul()`, ...) are currently implemented for the general case only. However, these would benefit accordingly if the specific functions for special matrix geometries (symmetric matrices, triangular matrices, band matrices, ...) were used in each case. The basic functions required for this are already implemented in `BLAS.au3` and `LAPACK.au3`.