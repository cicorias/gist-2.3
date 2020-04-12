/*****************************************************************************
 * FILE: fselect.h
 * AUTHOR: William Stafford Noble
 * CREATE DATE: 10/24/01
 * PROJECT: SVM
 * COPYRIGHT: 2001, Columbia University
 * VERSION: $Revision: 1.1 $
 * DESCRIPTION: Feature selection by Fisher criterion score.
 *****************************************************************************/
#ifndef FSELECT_H
#define FSELECT_H

#include "kernel.h"
#include "matrix.h"
#include "rdb-matrix.h"
#include "utils.h"
#include "class-array.h"

/*****************************************************************************
 * Reduce a given data set using feature selection.
 *****************************************************************************/
void select_features
  (KERNEL_T*      kernel,        // Kernel parameters.
   CLASS_ARRAY_T* classes,       // Classes.
   STRING_LIST_T* feature_names, // Names of all features.
   MATRIX_T*      train_matrix,  // Data from which to score features.
   MATRIX_T*      test_matrix,   // Second data to remove features from.
   ARRAY_T**      score_array,   // Score array.
   STRING_LIST_T* selected_feature_names, // Names of selected features.
   MATRIX_T**     reduced_train, // Training set with fewer features.
   MATRIX_T**     reduced_test); // Test set with fewer features.

/*****************************************************************************
 * Compute the standard t-test score, or Welch's approximation
 * thereof.
 *****************************************************************************/
double ttest
  (BOOLEAN_T welch,
   double mean1,
   double mean2,
   double var1,
   double var2,
   double num1,
   double num2);

/*****************************************************************************
 * Compute the Fisher criterion score.
 *****************************************************************************/
double fisher_score
  (double mean1,
   double mean2,
   double variance1,
   double variance2);

/*****************************************************************************
 * Mann-whitney u-test.
 *****************************************************************************/
double mannwhitney
  (ARRAY_T* positive_data,
   ARRAY_T* negative_data);

/****************************************************************************
 * TNoM feature selection, implemented according to the algorithm
 * described in "Tissue Classification with Gene Expression Profiles",
 * Ben-Dor et al., JCompBio, vol7,pp559-583.
 ***************************************************************************/
double tnom
  (ARRAY_T* positive_data,
   ARRAY_T* negative_data);

#endif
/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
