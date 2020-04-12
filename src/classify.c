/****************************************************************************
 * FILE: classify.c
 * AUTHOR: William Stafford Noble
 * CREATE DATE: 3/3/99
 * PROJECT: SVM
 * COPYRIGHT: 1999-2001, Columbia University (see ../doc/copyright.html)
 * VERSION: $Revision: 1.4 $
 * DESCRIPTION: Classify a set of examples using a kernel function.
 ****************************************************************************/
#include "kernel.h"
#include "fselect.h"
#include "classify.h"
#include "compute-weights.h"
#include "compute-kernel.h"
#include "rdb-matrix.h"
#include "string-list.h"
#include "matrix.h"
#include "array.h"
#include "utils.h"
#include <stdio.h>
#include <string.h>

#ifndef CLASSIFY_DEBUG
#define CLASSIFY_DEBUG 0
#endif

extern VERBOSE_T verbosity;

/****************************************************************************
 * Classify a single example.
 ****************************************************************************/
static double classify
  (int       i_test,
   ARRAY_T*  weights,
   MATRIX_T* kernel_matrix)
{
  double    return_value;
  int       num_train;
  int       i_train;
  double    this_weight;
  double    this_value;

  return_value = 0.0;
  num_train = get_array_length(weights);
  for (i_train = 0; i_train < num_train; i_train++) {

    // Get the current weight.
    this_weight = get_array_item(i_train, weights);

    // If the weight is zero, skip.
    if (this_weight == 0.0) {
      continue;
    }

    // Compute the kernel between the two examples.
    this_value = get_matrix_cell(i_test, i_train, kernel_matrix);

    if (CLASSIFY_DEBUG) {
      fprintf(stderr, "Example %d: %g * %g = %g\n",
	      i_train, this_value, this_weight, this_value * this_weight);
    }

    /* Weight the distance appropriately. This assumes that the
       classification of the training set example is encoded in the
       sign of the weight. */
    this_value *= this_weight;
    return_value += this_value;
  }

  if (CLASSIFY_DEBUG) {
    fprintf(stderr, "Total: %g\n", return_value);
  }
  return(return_value);
}

/****************************************************************************
 * Classify a single example, using the hyperplane coordinates.
 *
 * When a constant C is added to the kernel, this ends up adding a 
 * a term C \sum_i \alpha_i to each discriminant.
 ****************************************************************************/
double fast_classify
(ARRAY_T*  example,
 ARRAY_T*  hyperplane,
 double    sum_of_weights,
 KERNEL_T* kernel)
{
  // Make sure the kernel is linear.
  if (kernel->radial || (kernel->power != 1)) {
    die("Can't compute hyperplane of non-linear kernel.\n");
  }

  // For now, don't deal with coefficients.
  if (kernel->coefficient != 1) {
    die("Can't yet handle power=%g\n", kernel->coefficient);
  }

  // If necessary, normalize the example.
  if (kernel->normalize) {
    two_normalize(example);
  }

  // Compute the scalar product.
  return(dot_product(example, hyperplane) + 
	 (kernel->constant * sum_of_weights));
}



/****************************************************************************
 * Compute the hyperplane coordinates.
 ****************************************************************************/
void compute_hyperplane
(KERNEL_T* kernel,
 ARRAY_T*  weights,
 ARRAY_T** hyperplane)
{
  int num_train;
  int i_train;
  int num_features;
  int i_feature;
  ARRAY_T* example;

  // Make sure the kernel is linear.
  if (kernel->radial || (kernel->power != 1)) {
    die("Can't compute hyperplane of non-linear kernel.\n");
  }

  // Make sure the kernel doesn't have a coefficient.
  if (kernel->coefficient != 1) {
    die("Can't compute hyperplane of kernel with coefficient of %g.\n",
	kernel->coefficient);
  }

  // Get the size of the data set.
  num_train = get_array_length(weights);
  num_features = get_num_cols(kernel->train_matrix);

  // Allocate the hyperplane array.
  *hyperplane = allocate_array(num_features);
  init_array(0.0, *hyperplane);

  // Allocate local storage for the copy of the training example.
  example = allocate_array(num_features);

  // Compute the hyperplane coordinates.
  for (i_train = 0; i_train < num_train; i_train++) {
    ATYPE weight = get_array_item(i_train, weights);

    // If it's not a support vector, skip.
    if (weight == 0.0) {
      continue;
    }

    // Get a copy of this example.
    copy_array(get_matrix_row(i_train, kernel->train_matrix), example);

    // If necessary, normalize it.
    if (kernel->normalize) {
      two_normalize(example);
    }

    // Add to the hyperplane coordinates.
    for (i_feature = 0; i_feature < num_features; i_feature++) {
      incr_array_item(i_feature, 
		      get_array_item(i_feature, example) * weight,
		      *hyperplane);
    }
  }

  // Free dynamically allocated space.
  myfree(example);
}


/****************************************************************************
 * Classify a list of examples.
 ****************************************************************************/
void classify_list
  (int        num_test,
   double     bias,
   ARRAY_T*   weights,
   MATRIX_T*  kernel_matrix,
   ARRAY_T**  classifications,
   ARRAY_T**  discriminants)
{
  int    i_test;
  double this_discriminant;

  // Allocate the output arrays.
  *classifications = allocate_array(num_test);
  *discriminants = allocate_array(num_test);

  for (i_test = 0; i_test < num_test; i_test++) {

    if (verbosity >= HIGHER_VERBOSE) {
      fprintf(stderr, "Classifying example %d.\n", i_test);
    }

    // Compute the discriminant.
    this_discriminant = classify(i_test, weights, kernel_matrix) - bias;

    // Store the classification.
    if (this_discriminant >= 0.0) {
      set_array_item(i_test, 1.0, *classifications);
    } else {
      set_array_item(i_test, -1.0, *classifications);
    }

    // Store the discriminant.
    set_array_item(i_test, this_discriminant, *discriminants);
  }
}

/****************************************************************************
 * Copy an array into a second array, skipping a single value.
 ****************************************************************************/
static void fill_hold_out_array
  (ARRAY_T* array,
   int      skip,
   ARRAY_T* hold_out_array)
{
  int    num_array;
  int    i_array;
  double array_item;
  
  num_array = get_array_length(array);
  assert((num_array - 1) == get_array_length(hold_out_array));

  for (i_array = 0; i_array < num_array; i_array++) {

    // Is this the one to skip?
    if (i_array != skip) {
      array_item = get_array_item(i_array, array);

      // Copy from the array to the hold-out array.
      if (i_array < skip) {
	set_array_item(i_array, array_item, hold_out_array);
      } else {
	set_array_item(i_array - 1, array_item, hold_out_array);
      }
    }
  }
}

/****************************************************************************
 * Copy a class array into a second array, skipping a single value.
 ****************************************************************************/
static void fill_hold_out_class_array
  (CLASS_ARRAY_T* array,
   int            skip,
   CLASS_ARRAY_T* hold_out_array)
{
  int       num_array;
  int       i_array;
  BOOLEAN_T array_item;
  
  num_array = get_num_items(array);
  assert(num_array - 1 == get_num_items(hold_out_array));

  for (i_array = 0; i_array < num_array; i_array++) {

    // Is this the one to skip?
    if (i_array != skip) {
      array_item = get_class(i_array, array);

      // Copy from the array to the hold-out array.
      if (i_array < skip) {
	set_class(i_array, array_item, hold_out_array);
      } else {
	set_class(i_array - 1, array_item, hold_out_array);
      }
    }
  }
}

/****************************************************************************
 * Convert a matrix into two matrices: one containing a single
 * row, and the other containing the rest of the rows.
 ****************************************************************************/
static void fill_hold_out_matrix
  (int       hold_out,       // Index of row to hold out.
   MATRIX_T* full_matrix,    // N by M
   MATRIX_T* held_out_row,   // 1 by M
   MATRIX_T* reduced_matrix) // N-1 by M
{
  ARRAY_T* current_row;
  ARRAY_T* hold_out_row;
  int num_features;
  int num_examples;
  int i_example;

  num_examples = get_num_rows(full_matrix);
  num_features = get_num_cols(full_matrix);
  assert(num_examples == get_num_rows(reduced_matrix) + 1);
  assert(num_features == get_num_cols(held_out_row));
  assert(num_features == get_num_cols(reduced_matrix));

  for (i_example = 0; i_example < num_examples; i_example++) {
    current_row = get_matrix_row(i_example, full_matrix);

    // Copy the non-held out item into the reduced matrix.
    if (i_example < hold_out) {
      hold_out_row = get_matrix_row(i_example, reduced_matrix);
    }
    // Store the held out item in the hold-out matrix.
    else if (i_example == hold_out) {
      hold_out_row = get_matrix_row(0, held_out_row);
    }
    // Copy the non-held out item into the reduced matrix.
    else {
      hold_out_row = get_matrix_row(i_example - 1, reduced_matrix);
    }

    copy_array(current_row, hold_out_row);
  }
}


/****************************************************************************
 * Convert a square matrix into two matrices: one containing a single
 * row, and the other containing the rest of the rows.  The single row
 * is eliminated from both the rows and the columns of the original
 * matrix.
 ****************************************************************************/
static void fill_hold_out_kernel_matrix
  (int       hold_out,       // Index of row to hold out.
   MATRIX_T* full_matrix,    // N by N
   MATRIX_T* held_out_row,   // 1 by N-1
   MATRIX_T* reduced_matrix) // N-1 by N-1
{
  ARRAY_T* current_row;
  ARRAY_T* hold_out_row;
  int num_examples;
  int i_example;

  num_examples = get_num_cols(full_matrix);
  for (i_example = 0; i_example < num_examples; i_example++) {
    current_row = get_matrix_row(i_example, full_matrix);

    // Store the held out item in the hold-out matrix.
    if (i_example == hold_out) {
      fill_hold_out_array(current_row, hold_out, 
			  get_matrix_row(0, held_out_row));
    }
    // Copy the non-held out item into the reduced matrix.
    else {
      if (i_example < hold_out) {
	hold_out_row = get_matrix_row(i_example, reduced_matrix);
      } else {
	hold_out_row = get_matrix_row(i_example - 1, reduced_matrix);
      }
      fill_hold_out_array(current_row, hold_out, hold_out_row);
    }
  }
}


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
   ARRAY_T**      discriminants)
{
  int            num_examples;
  int            num_hold_out;
  int            i_hold_out;
  double         this_discriminant;
  ARRAY_T*       hold_out_weights;
  MATRIX_T*      hold_out_kernel_matrix;
  CLASS_ARRAY_T* hold_out_classes;
  MATRIX_T*      test_kernel_matrix;
  STRING_LIST_T* hold_out_names;
  static int     current_item; // Used by select_random_item.
  static int*    data;         // Used by select_random_item.
  MATRIX_T*  score_matrix = NULL; // Feature quality scores.
  RDB_MATRIX_T* score_matrix_rdb = NULL;

  // Make sure we're testing on the training set.
  num_examples = get_array_length(weights);
  myassert(TRUE, num_examples == get_num_rows(kernel->kernel_matrix),
	   "Performing hold-one-out cross-validation on non-square matrix.");
  myassert(TRUE, num_examples == get_num_cols(kernel->kernel_matrix),
	   "Performing hold-one-out cross-validation on non-square matrix.");

  // Allocate the output arrays and fill them with markers.
  *classifications = allocate_array(num_examples);
  *discriminants = allocate_array(num_examples);
  init_array(NaN(), *classifications);
  init_array(NaN(), *discriminants);

  // Allocate the smaller weights matrix and kernel matrix.
  hold_out_weights = allocate_array(num_examples - 1);
  hold_out_kernel_matrix = allocate_matrix(num_examples - 1, num_examples - 1);
  hold_out_classes = allocate_class_array(num_examples - 1);
  test_kernel_matrix = allocate_matrix(1, num_examples - 1);

  // Decide how many examples to hold out.
  num_hold_out = num_examples * hold_out;
  if (verbosity >= LOW_VERBOSE) {
    fprintf(stderr, "Will hold out %d out of %d examples.\n", num_hold_out,
	    num_examples);
  }

  // Allocate the score matrix, if necessary.
  hold_out_names = new_string_list();
  if (score_filename != NULL) {
    score_matrix = allocate_matrix(get_num_cols(kernel->train_matrix),
				   num_hold_out + 1);
				   
    // Put the first column of scores in.
    set_matrix_column(score_array, 0, score_matrix);
    add_string(convert_enum_type(kernel->feature_select, FSELECT_STRS, 
				NUM_FSELECT_T), hold_out_names);
    free_array(score_array);
  }

  current_item = -1;
  for (i_hold_out = 0; i_hold_out < num_hold_out; i_hold_out++) {
    int rand_item;

    // Randomly select an item to hold out.
    rand_item = select_random_item(num_examples, &current_item, &data);

    if (verbosity >= LOW_VERBOSE) {
      fprintf(stderr, "%d: Holding out example %d.\n", i_hold_out, rand_item);
    }

    // Fill the hold-out class array.
    fill_hold_out_class_array(classes, rand_item, hold_out_classes);

    // Re-compute if this example is a support vector or we're doing
    // feature selection.
    if ((get_array_item(rand_item, weights) != 0.0) || 
	(kernel->feature_select != NONE_FSELECT)) {
      char* hold_out_name 
	= get_nth_string(rand_item,
			 get_row_names(kernel->kernel_matrix_rdb));
      add_string(hold_out_name, hold_out_names);

      // Re-compute the kernel matrix, if necessary.
      if (kernel->feature_select != NONE_FSELECT) {
	int       num_features;
	MATRIX_T* hold_out_train;
	MATRIX_T* hold_out_test;
	MATRIX_T* reduced_train;
	MATRIX_T* reduced_test;
	ARRAY_T*  self_train = NULL;
	ARRAY_T*  self_test = NULL;

	// Create a smaller version of the training and test matrices.
	num_features = get_num_cols(kernel->train_matrix);
	hold_out_train = allocate_matrix(num_examples - 1, num_features);
	hold_out_test = allocate_matrix(1, num_features);
	fill_hold_out_matrix(rand_item,
			     kernel->train_matrix,
			     hold_out_test,
			     hold_out_train);

	// Select features.
	if (verbosity > NORMAL_VERBOSE) {
	  fprintf(stderr, "Selecting features.\n");
	}
	select_features(kernel,
			hold_out_classes,
			get_col_names(kernel->train_matrix_rdb),
			hold_out_train,
			hold_out_test,
			&score_array,
			NULL,
			&reduced_train,
			&reduced_test);

	// Add the scores to the matrix.
	if (score_matrix != NULL) {
	  set_matrix_column(score_array, i_hold_out + 1, score_matrix);
	}
	free_array(score_array);

	// Compute the base kernel matrices.
	if (verbosity > NORMAL_VERBOSE) {
	  fprintf(stderr, "Re-computing kernel matrix.\n");
	}
	compute_base_kernel_matrix(reduced_train,
				   reduced_train,
				   hold_out_kernel_matrix);
	compute_base_kernel_matrix(reduced_test,
				   reduced_train,
				   test_kernel_matrix);
	
	// Extract the diagonal from the hold-out kernel matrix.
	if (((kernel->radial) || (kernel->normalize))) {
	  self_train = extract_diagonal(hold_out_kernel_matrix);
	}

	// Turbo-ize the kernel matrices.
	transform_kernel_matrix(TRUE,
				kernel,
				hold_out_classes,
				self_train,
				self_train,
				hold_out_kernel_matrix);

	// Compute the diagonals.
	if (((kernel->radial) || (kernel->normalize))) {
	  self_train = compute_self_kernel_values(reduced_train);
	  self_test = compute_self_kernel_values(reduced_test);
	}

	transform_kernel_matrix(FALSE, 
				kernel,
				NULL,
				self_train,
				self_test,
				test_kernel_matrix);

	// Free dynamic memory.
	free_matrix(hold_out_train);
	free_matrix(hold_out_test);
	free_matrix(reduced_train);
	free_matrix(reduced_test);
	free_array(self_train);
	free_array(self_test);

      } else { // not feature selecting.
      
	// Fill in the smaller kernel matrix.
	fill_hold_out_kernel_matrix(rand_item,
				    kernel->kernel_matrix,
				    test_kernel_matrix,
				    hold_out_kernel_matrix);
      }

      // Initialize the weights to reasonable values.
      fill_hold_out_array(weights, rand_item, hold_out_weights);


      // Re-optimize the weights.
      RESET_CONVERGENCE();
      double_optimize_weights(convergence_threshold,
			      maxiter,
			      maxtime,
			      positive_constraint,
			      negative_constraint,
			      hold_out_classes,
			      hold_out_kernel_matrix,
			      hold_out_weights);

      // Encode the classifications as the signs of the weights.
      sign_weights(hold_out_classes, hold_out_weights);

      // Compute the discriminant with the re-computed weights.
      this_discriminant = classify(0, hold_out_weights, test_kernel_matrix);

    } else { // not feature selecting or not a support vector.

      // Compute the discriminant with the original weights.
      revert_kernel_diagonal(classes, kernel);
      myassert(1, !kernel->diagonal_isModified, 
	       "Attempting to classify with a modified kernel diagonal.\n");

      this_discriminant = classify(rand_item, weights, kernel->kernel_matrix);

      // Put the kernel diagonal back the way it was.
      revert_kernel_diagonal(classes, kernel);
    }


    // Store the classification.
    if (this_discriminant >= 0.0) {
      set_array_item(rand_item, 1.0, *classifications);
    } else {
      set_array_item(rand_item, -1.0, *classifications);
    }

    // Store the discriminant.
    set_array_item(rand_item, this_discriminant, *discriminants);
  }

  // FIXME: Move this block to train-main.c? -- WSN 3/29/05
  // Print the scores matrix, if requested.
  if (score_filename != NULL) {
    FILE* score_file;

    // Convert the matrix to RDB format.
    score_matrix_rdb = rdbize_matrix("fselect",
				     get_col_names(kernel->train_matrix_rdb),
				     hold_out_names,
				     score_matrix);

    // Open the file.
    if (open_file(score_filename, "w", FALSE, "score", "scores",
		  &score_file) == 0)
      exit(1);
    
    // Print it.
    print_rdb_matrix(score_matrix_rdb, format_line, 8, 6, score_file);
    fclose(score_file);
  }

  // Free dynamic memory.
  free_array(hold_out_weights);
  free_matrix(hold_out_kernel_matrix);
  free_class_array(hold_out_classes);
  free_matrix(test_kernel_matrix);
  free_string_list(hold_out_names);
  free_rdb_matrix(score_matrix_rdb);
}

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
 
