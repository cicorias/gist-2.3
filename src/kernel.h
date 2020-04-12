/****************************************************************************
 * FILE: kernel.h
 * AUTHOR: William Stafford Noble
 * CREATE DATE: 10/24/2001
 * PROJECT: SVM
 * COPYRIGHT: 1999-2001, Columbia University (see ../doc/copyright.html)
 * VERSION: $Revision: 1.1 $
 * DESCRIPTION: A data structure to store kernel function parameters.
 ****************************************************************************/
#ifndef KERNEL_H
#define KERNEL_H

#include "rdb-matrix.h"
#include "utils.h"
#include "class-array.h"
#include <stdio.h>

#define DEFAULT_PRECISION 6  // How many digits after decimal.
#define MINIMUM_PRECISION 4  // At least this many digits in output.

// Enumerated type for different ways of setting selection threshold.
typedef enum {INVALID_THRESH, PERCENT_THRESH, NUMBER_THRESH, 
	      VALUE_THRESH} THRESH_T;
extern char* THRESH_STRS[];
extern int NUM_THRESH_T;

// Enumerated type for different feature selection metrics.
typedef enum {INVALID_FSELECT, NONE_FSELECT, FISHER_FSELECT, TTEST_FSELECT, 
	      WELCH_FSELECT, MANNWHITNEY_FSELECT, TNOM_FSELECT, SAM_FSELECT,
	      SAM_FSELECT_FAST} FSELECT_T;
extern char* FSELECT_STRS[];
extern int NUM_FSELECT_T;


typedef struct kernel {

  char*          train_filename;    // Training data.
  FILE*          train_file;	    
  MATRIX_T*      train_matrix;
  RDB_MATRIX_T*  train_matrix_rdb;

  char*          test_filename;     // Test data.
  FILE*          test_file;
  MATRIX_T*      test_matrix;
  RDB_MATRIX_T*  test_matrix_rdb;   // (May be pointer to training matrix.)

  char*          self_train_filename; // Self-kernel values for training set.
  FILE*          self_train_file;
  ARRAY_T*       self_train;

  char*          self_test_filename; // Self-kernel values for test set.
  FILE*          self_test_file;
  ARRAY_T*       self_test;

  BOOLEAN_T      square_matrix;     // Are the rows and columns of the 
                                    // kernel matrix identical?
  BOOLEAN_T      matrix_from_file;  // Read base kernel matrix from a file.
  BOOLEAN_T      zero_mean_row;     // Make each data row have zero mean.
  BOOLEAN_T      zero_mean_col;     // Subtract mean from each feature.
  BOOLEAN_T      variance_one;      // Make each data row have variance one.
  FSELECT_T      feature_select;    // Use t-test (or Fisher score) to select.
  THRESH_T       thresh_type;       // Type of selection threshold.
  double         fthreshold;        // Feature selection threshold.
  BOOLEAN_T      normalize;         // Normalize rows of training matrix?
  double         constant;          // Constant to add to kernel.
  double         coefficient;       // Linear multiplier on kernel.
  double         power;             // Power to raise the kernel to.
  BOOLEAN_T      radial;            // Create radial basis kernel?
  double         width_factor;      // Multiplicative factor on width.
  double         width;             // User-specified width.
  double         two_squared_width; // Width of radial basis.
  double         add_diag;          // Add to the kernel diagonal.
  double         positive_diagonal; // Add to diagonal of kernel matrix.
  double         negative_diagonal; // Add to diagonal of kernel matrix.
  double         diagonal_factor;   // This one is complicated...
  double         regularizer;       // Regularizer to add to matrix diagonal.

  MATRIX_T*      kernel_matrix;     // Kernel values for train vs test.

  BOOLEAN_T      diagonal_isModified; // flag that tells us if the
				      // diagonal factors have been
				      // applied to this matrix

  BOOLEAN_T      kernel_out; // true = kernel matrix is destined for output.

  BOOLEAN_T      fthresholdset; // flag that tells us if the fthreshold has been set.

  RDB_MATRIX_T*  kernel_matrix_rdb;
} KERNEL_T;


/*****************************************************************************
 * Dynamically allocate one kernel object.
 *****************************************************************************/
KERNEL_T* allocate_kernel
  (BOOLEAN_T square_matrix);

/*****************************************************************************
 * Insert an RDB matrix and its corresponding raw matrix.
 *****************************************************************************/
void set_train_matrix
  (RDB_MATRIX_T* train_matrix_rdb,
   KERNEL_T*     kernel);

void set_test_matrix
  (RDB_MATRIX_T* test_matrix_rdb,
   KERNEL_T*     kernel);

void set_kernel_matrix
  (RDB_MATRIX_T* kernel_matrix_rdb,
   KERNEL_T*     kernel);

/*****************************************************************************
 * Set default values for feature selection parameters.
 *****************************************************************************/
void set_default_feature_selection
  (KERNEL_T* kernel);


/*****************************************************************************
 * Set the feature selection threshold
 *****************************************************************************/
void set_fthreshold
  (KERNEL_T* kernel, double fthreshold);

/*****************************************************************************
 * Read train and/or test set files.
 *****************************************************************************/
void read_data_files
  (BOOLEAN_T format_line,
   KERNEL_T* kernel);

/*****************************************************************************
 * Compute an array of kernel values of the form K(x,x).
 *****************************************************************************/
ARRAY_T* compute_self_kernel_values
  (MATRIX_T* data_matrix);

/*****************************************************************************
 * Compute the kernel matrix.
 *****************************************************************************/
void compute_base_kernel
  (KERNEL_T*      kernel,
   CLASS_ARRAY_T* classes,
   ARRAY_T**      score_array);

void free_kernel
  (KERNEL_T* kernel);


#endif
/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
