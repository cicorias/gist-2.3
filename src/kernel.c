/****************************************************************************
 * FILE: kernel.c
 * AUTHOR: William Stafford Noble
 * CREATE DATE: 10/24/2001
 * PROJECT: SVM
 * COPYRIGHT: 1999-2001, Columbia University (see ../doc/copyright.html)
 * VERSION: $Revision: 1.1 $
 * DESCRIPTION: A data structure to store kernel function parameters.
 ****************************************************************************/
#include "compute-kernel.h"
#include "kernel.h"
#include "fselect.h"
#include "utils.h"
#include <string.h>

// Default kernel settings.
#define DEFAULT_DIAG_FACTOR 0.1
#define DEFAULT_KERNEL_CONSTANT 10

// Enumerated type for spacer emission probabilities.
char* THRESH_STRS[] = {"invalid", "percent", "number", "value"};
int NUM_THRESH_T = 4;

// Enumerated type for different feature selection metrics.
char* FSELECT_STRS[] = {"invalid", "none", "fisher", "ttest", "welch", "mannwhitney", "tnom", "sam", "samf"};
int NUM_FSELECT_T = 9;

/***************************************************************************
 * Dynamically allocate space for a kernel object.
 ***************************************************************************/
KERNEL_T* allocate_kernel
  (BOOLEAN_T square_matrix)
{
  KERNEL_T* kernel = mycalloc(1, sizeof(KERNEL_T));

  kernel->square_matrix = square_matrix;
  kernel->matrix_from_file = FALSE;
  kernel->zero_mean_row = FALSE;
  kernel->zero_mean_col = FALSE;
  kernel->variance_one = FALSE;
  kernel->feature_select = NONE_FSELECT;
  kernel->thresh_type = PERCENT_THRESH;
  kernel->fthreshold = 0.0; // This gets set later.
  kernel->normalize = TRUE;
  kernel->constant = DEFAULT_KERNEL_CONSTANT;
  kernel->coefficient = 1.0;
  kernel->power = 1.0;
  kernel->radial = FALSE;
  kernel->width_factor = 1.0;
  kernel->width = 0;
  kernel->add_diag = 0.0;
  kernel->positive_diagonal = 0.0;
  kernel->negative_diagonal = 0.0;
  kernel->diagonal_factor = DEFAULT_DIAG_FACTOR;
  kernel->diagonal_isModified = FALSE;
  kernel->kernel_out = FALSE;
  kernel->regularizer = 1.0;

  return(kernel);
}

/*****************************************************************************
 * Insert an RDB matrix and its corresponding raw matrix.
 *****************************************************************************/
void set_train_matrix
  (RDB_MATRIX_T* train_matrix_rdb,
   KERNEL_T*     kernel)
{
  kernel->train_matrix_rdb = train_matrix_rdb;
  kernel->train_matrix = get_raw_matrix(train_matrix_rdb);
}

void set_test_matrix
  (RDB_MATRIX_T* test_matrix_rdb,
   KERNEL_T*     kernel)
{
  kernel->test_matrix_rdb = test_matrix_rdb;
  if (test_matrix_rdb == NULL) {
    kernel->test_matrix = NULL;
  } else {
    kernel->test_matrix = get_raw_matrix(test_matrix_rdb);
  }
}

void set_kernel_matrix
  (RDB_MATRIX_T* kernel_matrix_rdb,
   KERNEL_T*     kernel)
{
  kernel->kernel_matrix_rdb = kernel_matrix_rdb;
  kernel->kernel_matrix = get_raw_matrix(kernel_matrix_rdb);

}


/*****************************************************************************
 * Set default values for feature selection parameters.
 *****************************************************************************/
void set_default_feature_selection
  (KERNEL_T* kernel)
{
  // Make sure feature selection is on.
  if (kernel->feature_select != NONE_FSELECT) {

    // Make sure the threshold didn't already get set.
    if (kernel->fthreshold == 0.0) {

      // 10% or top 10 features.
      if ((kernel->thresh_type == PERCENT_THRESH)
	  || (kernel->thresh_type == NUMBER_THRESH)) {
	kernel->fthreshold = 10.0;
      }
      else if (kernel->thresh_type == VALUE_THRESH) {
	kernel->fthreshold = 1.0;
      }
    }
  }
}

/*****************************************************************************
 * Set the feature selection threshold
 *****************************************************************************/
void set_fthreshold
  (KERNEL_T* kernel, double fthreshold)
{
  kernel->fthreshold = fthreshold;
  kernel->fthresholdset = TRUE;
}



/*****************************************************************************
 * Read train and/or test set files.
 *****************************************************************************/
void read_data_files
  (BOOLEAN_T format_line,
   KERNEL_T* kernel)
{
  RDB_MATRIX_T* temp_matrix;
  
  if (verbosity > NORMAL_VERBOSE)
    fprintf(stderr, "Reading data\n");

  if (!kernel->matrix_from_file) {

    // Read the training set examples.
    if (open_file(kernel->train_filename, "r", TRUE, "train", "training set",
		  &(kernel->train_file)) == 0)
      exit(1);
    set_train_matrix(read_rdb_matrix(format_line, NULL, kernel->train_file),
		     kernel);
    myassert(1, get_num_cols(kernel->train_matrix) > 0, "No columns in training matrix");
    myassert(1, get_num_rows(kernel->train_matrix) > 0, "No rows in training matrix");

    if (verbosity > HIGH_VERBOSE)
      fprintf(stderr, "Read data from %s\n", kernel->train_filename);

  } else {

    // Read the kernel matrix from a file.
    if (kernel->train_filename != NULL ) {
      if (open_file(kernel->train_filename, "r", TRUE, "kernel",
		    "kernel values", &(kernel->train_file)) == 0)
	exit(1);
      set_kernel_matrix(read_rdb_matrix(format_line, NULL, kernel->train_file),
			kernel);
    } else { /* testing */
      if (open_file(kernel->test_filename, "r", TRUE, "kernel",
		    "kernel values", &(kernel->test_file)) == 0)
	exit(1);
      set_kernel_matrix(read_rdb_matrix(format_line, NULL, kernel->test_file),
			kernel);
    }
    myassert(1, get_num_cols(kernel->kernel_matrix) > 0, "No columns in kernel matrix");
    myassert(1, get_num_rows(kernel->kernel_matrix) > 0, "No rows in kernel matrix");

    if (verbosity > HIGH_VERBOSE)
      fprintf(stderr, "Read kernel matrix from %s\n", kernel->train_filename);
  }

  // Read the test set.
  if (kernel->test_filename != NULL) {
    if (open_file(kernel->test_filename, "r", TRUE, "test", "test values",
		  &(kernel->test_file)) == 0)
      exit(1);
    
    set_test_matrix(read_rdb_matrix(format_line, NULL, kernel->test_file),
		    kernel);
    myassert(1, get_num_cols(kernel->test_matrix) > 0, "No columns in test matrix");
    myassert(1, get_num_rows(kernel->test_matrix) > 0, "No rows in test matrix");
    if (kernel->train_matrix != NULL) {
      myassert(1, get_num_cols(kernel->train_matrix) == get_num_cols(kernel->test_matrix),
	       "Number of features in training data (%d) does not equal that in the test data (%d)", 
	       get_num_cols(kernel->train_matrix), 
	       get_num_cols(kernel->test_matrix));
    }
  } else {
    set_test_matrix(kernel->train_matrix_rdb, kernel);
  }
  
  // Read the self kernel values.
  if (kernel->self_train_filename != NULL) {
    if (open_file(kernel->self_train_filename, "r", TRUE, "self-train",
		  "self-train values", &(kernel->self_train_file)) == 0)
      exit(1);
    temp_matrix = read_rdb_matrix(format_line, NULL, 
				  kernel->self_train_file);
    kernel->self_train = get_matrix_column(0, get_raw_matrix(temp_matrix));
    free_rdb_matrix(temp_matrix);

 /* we are classifying, so the kernel matrix should have the number of
    columns equal to the number of training examples */
    myassert(1, get_array_length(kernel->self_train) == get_num_cols(kernel->kernel_matrix), 
	     "Self train file does not have the same number of items as the kernel matrix has columns (%d != %d).", 
	     get_array_length(kernel->self_train), 
	     get_num_cols(kernel->kernel_matrix));

  }

  if (kernel->self_test_filename != NULL) {
    if (open_file(kernel->self_test_filename, "r", TRUE, "self-test",
		  "self-test values", &(kernel->self_test_file)) == 0)
      exit(1);
    temp_matrix = read_rdb_matrix(format_line, NULL, 
				  kernel->self_test_file);
    kernel->self_test = get_matrix_column(0, get_raw_matrix(temp_matrix));
    free_rdb_matrix(temp_matrix);

    myassert(1, get_array_length(kernel->self_test) == get_num_rows(kernel->kernel_matrix), 
	     "Self test file does not have the same number of items as the kernel matrix has rows (%d != %d).", 
	     get_array_length(kernel->self_test), 
	     get_num_rows(kernel->kernel_matrix));

  }
}

/*****************************************************************************
 * Compute an array of kernel values of the form K(x,x).
 *****************************************************************************/
ARRAY_T* compute_self_kernel_values
  (MATRIX_T* data_matrix)
{
  ARRAY_T* self_values;
  int      num_rows;
  int      i_row;
  ARRAY_T* example;
  double   value;

  if (data_matrix == NULL) {
    return(NULL);
  }

  /* Allocate the kernel values array. */
  num_rows = get_num_rows(data_matrix);
  self_values = allocate_array(num_rows);

  /* Extract the kernel for every example against itself. */
  for (i_row = 0; i_row < num_rows; i_row++) {
    example = get_matrix_row(i_row, data_matrix);

    /* Compute the dot product against itself. */
    value = dot_product(example, example);

    /* Store the resulting distance in the array. */
    set_array_item(i_row, value, self_values);
  }
  return(self_values);
}

/*****************************************************************************
 * Compute the kernel matrix.
 *****************************************************************************/
void compute_base_kernel
  (KERNEL_T*      kernel,
   CLASS_ARRAY_T* classes,
   ARRAY_T**      score_array)
{
  // We may have already read the base kernel matrix
  if (kernel->kernel_matrix == NULL) {

    // Normalize the data.
    if (kernel->zero_mean_row) {
      if (verbosity > NORMAL_VERBOSE) {
	fprintf(stderr, "Subtracting the mean from each row.\n");
      }
      zero_mean_matrix_rows(kernel->train_matrix);
      if (kernel->test_matrix != kernel->train_matrix) {
	zero_mean_matrix_rows(kernel->test_matrix);
      }
    }
    if (kernel->variance_one) {
      if (verbosity > NORMAL_VERBOSE) {
	fprintf(stderr, "Dividing by the standard deviation.\n");
      }
      variance_one_matrix_rows(kernel->train_matrix);
      if (kernel->test_matrix != kernel->train_matrix) {
	variance_one_matrix_rows(kernel->test_matrix);
      }
    }

    // Allocate the kernel matrix.
    set_kernel_matrix(allocate_rdb_matrix(get_num_rows(kernel->test_matrix),
					  get_num_rows(kernel->train_matrix),
					  NULL),
		      kernel);
    set_corner_string(get_corner_string(kernel->train_matrix_rdb),
		      kernel->kernel_matrix_rdb);

    if (kernel->feature_select != NONE_FSELECT) {
      MATRIX_T* reduced_train;
      MATRIX_T* reduced_test;
      assert(score_array != NULL);

      // Select features.
      if (verbosity >= NORMAL_VERBOSE) {
	fprintf(stderr, "Selecting features.\n");
      }
      select_features(kernel,
		      classes,
		      NULL, // Don't care about feature names.
		      kernel->train_matrix,
		      kernel->test_matrix,
		      score_array,
		      NULL,
		      &reduced_train,
		      &reduced_test);

      // Compute the kernel matrix.
      compute_base_kernel_matrix(reduced_test,
				 reduced_train,
				 kernel->kernel_matrix);

      free_matrix(reduced_train);
      if (reduced_train != reduced_test) {
	free_matrix(reduced_test);
      }

    } else {

      // Compute the kernel matrix.
      compute_base_kernel_matrix(kernel->test_matrix,
				 kernel->train_matrix,
				 kernel->kernel_matrix);
    }


    // Add labels to the kernel matrix.
    set_row_names(get_row_names(kernel->test_matrix_rdb), 
		  kernel->kernel_matrix_rdb);
    set_col_names(get_row_names(kernel->train_matrix_rdb), 
		  kernel->kernel_matrix_rdb);

  }

  // We don't always need the diagonals.
  // N.B. Added third clause in 'if' 7-24-03.  -- WSN
  if ((kernel->radial || kernel->normalize || (kernel->add_diag != 0.0))) {

    // During training, we have already computed the diagonals.
    if (kernel->square_matrix) {
      kernel->self_train = extract_diagonal(kernel->kernel_matrix);
      kernel->self_test = extract_diagonal(kernel->kernel_matrix);
    }
    // During testing, we may have already read diagonals from files.
    else if ((kernel->self_train == NULL) || (kernel->self_test == NULL)) {
      kernel->self_train = compute_self_kernel_values(kernel->train_matrix);
      kernel->self_test = compute_self_kernel_values(kernel->test_matrix);
    }
  }

}

/***************************************************************************
 * Free dynamic memory for a kernel object.
 ***************************************************************************/
void free_kernel
  (KERNEL_T* kernel)
{

  if (kernel == NULL) {
    return;
  }

  free_rdb_matrix(kernel->train_matrix_rdb);
  if (kernel->train_matrix_rdb != kernel->test_matrix_rdb) {
    free_rdb_matrix(kernel->test_matrix_rdb);
  }
  free_array(kernel->self_train);
  free_array(kernel->self_test);
  free_rdb_matrix(kernel->kernel_matrix_rdb);
  myfree(kernel);
}


/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
