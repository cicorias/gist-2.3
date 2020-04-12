/*****************************************************************************
 * FILE: kernel-pca.h
 * AUTHOR: William Stafford Noble
 * CREATE DATE: 12/14/99
 * PROJECT: SVM
 * COPYRIGHT: 1999-2001, Columbia University (see ../doc/copyright.html)
 * VERSION: $Revision: 1.1 $
 * DESCRIPTION: Given a set of training examples in RDB format,
 * compute kernel-based eigenvectors.  This program follows the
 * outline given in "Neural Networks" by Haykin (1999), p. 435-436.
 *****************************************************************************/
#ifndef KERNEL_PCA_H
#define KERNEL_PCA_H

#include "rdb-matrix.h"
#include "matrix.h"

/*****************************************************************************
 * Find eigenvectors of a given matrix, normalized so that the dot
 * product of the eigenvector with itself equals the reciprocal of the
 * corresponding eigenvalue.  In the output, the eigenvectors are
 * sorted by increasing magnitude.
 *****************************************************************************/
MATRIX_T* find_normalized_eigenvectors
  (int       num_eigens,          /* Maximum number of eigenvectors to find. */
   double    eigen_threshold,     /* Only print eigenvectors whose 
				     corresponding eigenvalues account for at
				     least this much of the variance. */
   char*     eigenvalue_filename, /* Print the eigenvalues to this file. */
   MATRIX_T* matrix);             /* Matrix to analyze. */

/*****************************************************************************
 * Zero-mean the kernel matrix.  See .c file for full description.
 *****************************************************************************/
void center_kernel_matrix
  (MATRIX_T* kernel_matrix);

/*****************************************************************************
 * A goofy helper function that assigns names to the columns in the
 * eigenvector matrix.
 *****************************************************************************/
void add_eigen_col_names
  (RDB_MATRIX_T* rdb_eigenvectors);

#endif
/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */

