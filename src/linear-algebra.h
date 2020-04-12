/**************************************************************************
 * FILE: linear-algebra.h
 * AUTHOR: William Stafford Noble
 * CREATE DATE: 2-4-97
 * PROJECT: SVM
 * COPYRIGHT: 1999-2001, Columbia University (see ../doc/copyright.html)
 * VERSION: $Revision: 1.6 $
 * DESCRIPTION: Generic linear algebra routines.
 **************************************************************************/
#ifndef LINEAR_ALGEBRA_H
#define LINEAR_ALGEBRA_H

#include "matrix.h"

/* How many digits of precision do we require? */
#define MATRIX_PRECISION 0.0000001

/***********************************************************************
 * Turn a given square matrix into an identity matrix.
 ***********************************************************************/
void make_identity_matrix
  (MATRIX_T* matrix);

/***********************************************************************
 * Transpose a matrix.
 ***********************************************************************/
MATRIX_T* transpose_matrix
  (MATRIX_T* matrix);

/***********************************************************************
 * Invert a matrix.
 *
 * This function implements the Gauss-Jordan method for computing the
 * inverse.  See pp.43-44 of (Strang 1986).
 *
 * If the given matrix is not invertible, return NULL.
 ***********************************************************************/
MATRIX_T* invert_matrix
  (MATRIX_T* matrix);

/***********************************************************************
 * Multiply two matrices to get a third.
 *
 * Allocates memory for the third matrix, which must be freed by the
 * caller.
 ***********************************************************************/
MATRIX_T* matrix_multiply
  (MATRIX_T* matrix1,
   MATRIX_T* matrix2,
   MATRIX_T* matrix3);

/***********************************************************************
 * Compute the covariance matrix of a given matrix.
 ***********************************************************************/
MATRIX_T* compute_covariance
  (MATRIX_T* matrix);

/***********************************************************************
 * Compute all correlations between pairs of rows or columns in two 
 * matrices.  Allocates and returns a new matrix.
 ***********************************************************************/
MATRIX_T* compute_pairwise_correlations
(BOOLEAN_T use_rows, // else use columns
 MATRIX_T* matrix1,
 MATRIX_T* matrix2);

/***********************************************************************
 * Compute Jacard similarity between pairs of rows in a matrix.
 ***********************************************************************/
MATRIX_T* compute_pairwise_jacard
  (MATRIX_T* matrix);

/**************************************************************************
 * Compute the eigenvectors and the eigenvalues of a symmetric matrix.
 **************************************************************************/
void find_symmetric_eigenvectors
  (MATRIX_T*  matrix,
   ARRAY_T**  eigenvalues,
   MATRIX_T** eigenvectors);

/****************************************************************************
 * Find the largest eigenvector of a given matrix.
 ****************************************************************************/
ARRAY_T* find_first_eigenvector
  (MATRIX_T*  matrix);

/**************************************************************************
 * Find the smallest eigenvector of a given matrix.
 **************************************************************************/
double find_last_eigenvalue
  (MATRIX_T*  matrix);

/****************************************************************************
 * Force a square matrix to be positive definite by adding to the diagonal.
 ****************************************************************************/
void make_positive_definite
  (MATRIX_T* matrix);

/*****************************************************************************
 * Convert a kernel matrix to a Euclidean distance matrix, or vice versa
 *****************************************************************************/
MATRIX_T* kernel_to_distance
(MATRIX_T* matrix);

/*****************************************************************************
 * Compute the intrinsic dimensionality of a Euclidean distance
 * matrix, which is defined as the mean inter-object distance, divided
 * by the corresponding variance.
 *****************************************************************************/
double get_intrinsic_dimensionality
(MATRIX_T* matrix);

/*****************************************************************************
 * Perform a diffusion on a given symmetrix distance matrix.
 *****************************************************************************/
void perform_diffusion
(float     diffusion_constant,
 MATRIX_T* matrix);

/****************************************************************************
 * Sort a matrix by column, according to a given set of sort keys.
 * This is in here, rather than matrix.c, because it requires the
 * transpose function.
 ****************************************************************************/
void sort_matrix_cols
  (BOOLEAN_T reverse_sort,
   ARRAY_T*  keys,
   MATRIX_T* matrix);

/****************************************************************************
 * Compute the matrix exponential.
 * 
 * If X has a full set of eigenvectors V with corresponding
 * eigenvalues D, then [V,D] = EIG(X) and EXPM(X) =
 * V*diag(exp(diag(D)))/V.
 * 
 * This function requires that the given matrix is symmetric along the
 * diagonal.
 ****************************************************************************/
void matrix_exponential
(MATRIX_T* matrix);

/****************************************************************************
 * Compute a diffusion on the given matrix.
 *
 * This function requires a square matrix.
 ****************************************************************************/
MATRIX_T* diffuse_matrix
(MTYPE     diffusion_constant,
 int       iterations,
 MATRIX_T* matrix);

#endif

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
