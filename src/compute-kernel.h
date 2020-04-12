/*****************************************************************************
 * FILE: compute-kernel.h
 * AUTHOR: William Stafford Noble
 * CREATE DATE: 3/1/99
 * PROJECT: SVM
 * COPYRIGHT: 1999-2001, Columbia University (see ../doc/copyright.html)
 * VERSION: $Revision: 1.2 $
 * DESCRIPTION: Given an m x n matrix of training data (m examples
 * with n features for each), compute an m x m matrix of kernel values.
 *****************************************************************************/
#ifndef COMPUTE_KERNEL_H
#define COMPUTE_KERNEL_H

#include "kernel.h"
#include "rdb-matrix.h"
#include "matrix.h"
#include "class-array.h"
#include "linear-algebra.h"
#include "utils.h"
#include <stdio.h>
#include <assert.h>

/*****************************************************************************
 * Given a two data matrices, compute a matrix of pairwise dot
 * products.  K[i,j] = M1[i] . M2[j]
 *****************************************************************************/
void compute_base_kernel_matrix
  (MATRIX_T* matrix1,
   MATRIX_T* matrix2,
   MATRIX_T* kernel_matrix);

/*****************************************************************************
 * Perform the following operations on the given kernel matrix:
 *
 *  - add constant to diagonal
 *  - normalize
 *  - polynomialize
 *  - radialize
 *  - add asymmetrically to diagonal
 *  - set features to have zero mean
 *
 *****************************************************************************/
void transform_kernel_matrix
  (BOOLEAN_T      square, // Is this a square matrix?
   KERNEL_T*      kernel,
   CLASS_ARRAY_T* classes,
   ARRAY_T*       self_train,
   ARRAY_T*       self_test,
   MATRIX_T*      kernel_matrix);

/***************************************************************************
 * Do/Undo the kernel diagonal additon. We do this prior to classifying
 * the training examples.
 ***************************************************************************/
void revert_kernel_diagonal
  (CLASS_ARRAY_T* classes, 
   KERNEL_T*      kernel);

/*****************************************************************************
 * Compute all kernel transformations requested in a given file.
 *****************************************************************************/
void compute_kernel_transformations
(char*     kernel_filename,
 MATRIX_T* kernel_matrix);

#endif



/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
