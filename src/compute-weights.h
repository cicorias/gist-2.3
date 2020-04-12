/*****************************************************************************
 * FILE: compute-weights.h
 * AUTHOR: William Stafford Noble
 * CREATE DATE: 1/26/2000
 * PROJECT: SVM
 * COPYRIGHT: 1999-2001, Columbia University (see ../doc/copyright.html)
 * VERSION: $Revision: 1.1 $
 * DESCRIPTION: Given a matrix of pairwise distances for positive and
 * negative examples, compute weights that maximize a quadratic
 * objective function.
 *
 * This program is based upon
 *
 *   "A discriminative framework for detecting remote protein
 *   homologies."  T. Jaakkola, M. Diekhans, D. Haussler.  1998.
 *
 *****************************************************************************/
#ifndef COMPUTE_WEIGHTS_H
#define COMPUTE_WEIGHTS_H
#include "matrix.h"
#include "array.h"
#include "class-array.h"
#include "utils.h"

/*****************************************************************************
 * Randomly return an index between 0 and a given value.
 *
 * On the first call, current_item should be -1.
 *
 * The function allocates the data array and randomly fill it with
 * ascending integers.  On subsequent calls, iterate through the
 * array, returning the next index.  Once we iterate through the whole
 * thing, re-shuffle and start over.
 *****************************************************************************/
int select_random_item
  (int   num_items,
   int*  current_item,
   int** data);

/*****************************************************************************
 * Keep a local copy of the weights array and signal if they've
 * stopped changing.
 *
 * Convergence is reached when the delta is below the convergence
 * threshold.
 *****************************************************************************/
#define RESET_CONVERGENCE() (converged(TRUE, 0, NULL, NULL, NULL))
BOOLEAN_T converged
  (BOOLEAN_T      reset,  /* Reset everything for a new optimization. */
   double         convergence_threshold,
   ARRAY_T*       weights,
   CLASS_ARRAY_T* classes,
   MATRIX_T*      kernel_matrix);

/*****************************************************************************
 * Compute the objective function, equation (7) (but divide by the
 * number of items).
 *****************************************************************************/
double compute_objective
  (ARRAY_T*       weights,
   CLASS_ARRAY_T* classes,
   MATRIX_T*      kernel_matrix);

/*****************************************************************************
 * Optimize the weights on the entire training set, then extract the
 * support vectors and optimize them some more.
 *****************************************************************************/
long int double_optimize_weights
  (double         convergence_threshold,
   long int       maxiter,
   double         maxtime,
   double         positive_constraint,
   double         negative_constraint,
   CLASS_ARRAY_T* classes,
   MATRIX_T*      kernel_matrix,
   ARRAY_T*       weights);

/****************************************************************************
 * Encode the classifications in the weights by multiplying the negative
 * examples by -1.
 ****************************************************************************/
void sign_weights
  (CLASS_ARRAY_T* classes,
   ARRAY_T*       weights);

#endif
