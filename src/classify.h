/****************************************************************************
 * FILE: classify.h
 * AUTHOR: William Stafford Noble
 * CREATE DATE: 1/25/2000
 * PROJECT: SVM
 * COPYRIGHT: 1999-2001, Columbia University (see ../doc/copyright.html)
 * VERSION: $Revision: 1.3 $
 * DESCRIPTION: Classify a set of examples using a kernel function.
 ****************************************************************************/
#ifndef CLASSIFY_H
#define CLASSIFY_H

#include "kernel.h"
#include "matrix.h"
#include "class-array.h"
#include "array.h"

/****************************************************************************
 * Compute the hyperplane coordinates.
 ****************************************************************************/
void compute_hyperplane
(KERNEL_T* kernel,
 ARRAY_T*  weights,
 ARRAY_T** hyperplane);

/****************************************************************************
 * Classify a single example, using the hyperplane coordinates.
 ****************************************************************************/
double fast_classify
(ARRAY_T*  example,
 ARRAY_T*  hyperplane,
 double    sum_of_weights,
 KERNEL_T* kernel);

/****************************************************************************
 * Classify a list of examples.
 ****************************************************************************/
void classify_list
  (int        num_test,
   double     bias,
   ARRAY_T*   weights,
   MATRIX_T*  kernel_matrix,
   ARRAY_T**  classifications,
   ARRAY_T**  discriminants);

/****************************************************************************
 * Classify a list of examples, performing hold-one-out
 * cross-validation.  Note that this only works if the kernel matrix
 * is square.
 ****************************************************************************/
void hold_out_classify_list
  (double         hold_out,
   double         convergence_threshold,
   long int       maxiter,
   double         maxtime,
   double         positive_constraint,
   double         negative_constraint,
   BOOLEAN_T      format_line,
   char*          score_filename,
   KERNEL_T*      kernel,
   CLASS_ARRAY_T* classes,
   ARRAY_T*       weights,
   ARRAY_T*       score_array,
   ARRAY_T**      classifications,
   ARRAY_T**      discriminants);

#endif
