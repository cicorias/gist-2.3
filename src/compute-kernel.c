/*****************************************************************************
 * FILE: compute-kernel.c
 * AUTHOR: William Stafford Noble
 * CREATE DATE: 3/1/99
 * PROJECT: SVM
 * COPYRIGHT: 1999-2001, Columbia University (see ../doc/copyright.html)
 * VERSION: $Revision: 1.10 $
 * DESCRIPTION: Given an m x n matrix of training data (m examples
 *              with n features for each), compute an m x m matrix
 *              of kernel values.
 *****************************************************************************/
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <string.h>
// #include <dlfcn.h>
#include "linear-algebra.h"
#include "matrix.h"
#include "array.h"
#include "class-array.h"
#include "utils.h"
#include "cmdline.h"
#include "rdb-matrix.h"
#include "compute-kernel.h"

#ifndef KERNEL_DEBUG
#define KERNEL_DEBUG 0
#endif

#ifndef RADIAL_DEBUG
#define RADIAL_DEBUG 0
#endif

#ifndef NORMALIZE_DEBUG
#define NORMALIZE_DEBUG 0
#endif

#ifndef DIAG_DEBUG
#define DIAG_DEBUG 0
#endif

// Enumerated type of kernel transformations.
typedef enum {INVALID_TRANSFORM,
	      CONSTANT_TRANSFORM,
	      COEFFICIENT_TRANSFORM,
	      POWER_TRANSFORM,
	      RADIAL_TRANSFORM,
	      DIAGONAL_TRANSFORM,
	      NORMALIZE_TRANSFORM,
	      CENTER_TRANSFORM,
	      DIFFUSION_TRANSFORM
} TRANSFORM_T;

/*****************************************************************************
 * Compute a jackknife dot product.
 *
 * Let d_i(x,y) be the dot-product of x and y with the ith element
 * deleted.  Then the jackknife dot-product for two vectors of length
 * N is min{d_1(x,y), d_2(x,y), ..., d_N(x,y), d(x,y)}.
 *****************************************************************************/
double jackknife_dot_product
(ARRAY_T* example1,
 ARRAY_T* example2)
{
  int num_values;
  int i_value;
  double value1;
  double value2;
  double minimum_dot_product;
  double new_dot_product;

  minimum_dot_product = dot_product(example1, example2);
  num_values = get_array_length(example1);
  for (i_value = 0; i_value < num_values; i_value++) {

    /* Store the values. */
    value1 = get_array_item(i_value, example1);
    value2 = get_array_item(i_value, example2);

    /* Replace them with zero. */
    set_array_item(i_value, 0.0, example1);
    set_array_item(i_value, 0.0, example2);

    /* Compute the new dot product. */
    new_dot_product = dot_product(example1, example2);

    /* Have we found a new minimum? */
    if (new_dot_product < minimum_dot_product) {
      minimum_dot_product = new_dot_product;
    }

    /* Put the original values back in. */
    set_array_item(i_value, value1, example1);
    set_array_item(i_value, value2, example2);
  }
  return(minimum_dot_product);
}


/*****************************************************************************
 * Normalize a kernel matrix.
 *
 * Let x be a vector of unnormalized features. Let \tilde{x} be the
 * vector of normalized features. 
 * 
 *         \tilde{x}_i = x_i/||x||
 * 
 * where ||x|| is the norm of x. This means that ||\tilde{x}|| = 1 for
 * all x. What is happening when you normalize in this way is that all
 * feature vectors x are getting their lengths scaled so that they lie
 * on the surface of the unit sphere in the input space.
 * 
 * It turns out that this kind of normalization can be done in general
 * for any kernel. If K(x,y) is any kernel supplied by a user, then
 * you can normalize it by defining
 * 
 *     \tilde{K}(x,y) = K(x,y)/\sqrt{K(x,x)}\sqrt{K(y,y)}
 * 
 * It turns out that \sqrt{K(x,x)} is the norm of x in the feature
 * space.  Hence, for any x, the norm of x in the feature space
 * defined by \tilde{K} is
 * 
 *           \sqrt{\tilde{K}(x,x)}  = 1
 * 
 * When you normalize the kernel this way, it ensures that all points
 * are mapped to the surface of the unit ball in some (possibly infinite
 * dimensional) feature space. 
 *****************************************************************************/
#define SMALL_VALUE  1E-10
static void normalize_kernel_matrix
  (ARRAY_T*  row_self_values,
   ARRAY_T*  col_self_values,
   MATRIX_T* kernel_matrix)
{
  int num_rows;
  int i_row;
  int num_cols;
  int i_col;
  double this_cell;
  double product;
  static BOOLEAN_T warning = FALSE;

  // Get the dimensions of the matrix.
  num_rows = get_num_rows(kernel_matrix);
  myassert(1, num_rows == get_array_length(row_self_values),
	   "Number of kernel matrix rows (%d) does not equal number of self values (%d)",
	   num_rows, get_array_length(row_self_values));
  num_cols = get_num_cols(kernel_matrix);
  myassert(1, num_cols == get_array_length(col_self_values),
	   "Number of kernel matrix columns (%d) does not equal number of self values (%d)",
	   num_cols, get_array_length(col_self_values));

  for (i_row = 0; i_row < num_rows; i_row++) {
    double row_diag = get_array_item(i_row, row_self_values);
    double row_diag_sqrt = sqrt(row_diag);

    for (i_col = 0; i_col < num_cols; i_col++) {
      double col_diag = get_array_item(i_col, col_self_values);

      // Get this cell's value.
      this_cell = get_matrix_cell(i_row, i_col, kernel_matrix);

      // If both diagonals are zero, then both points are at the origin.
      if ((row_diag == 0) && (col_diag == 0)) {
	this_cell = 1;
      }

      // If one diagonal is zero, then set kernel value to zero.
      else if ((row_diag == 0) || (col_diag == 0)) {
	this_cell = 0;
      }

      // Otherwise, just do the calculation.
      else {
	product = row_diag_sqrt * sqrt(col_diag);
	this_cell /= product;

	// Check for NaN.
	if (isnan(this_cell) || (this_cell > 1.0)) {

	  // Tell the user what's happening.
	  if (warning == FALSE) {
	    fprintf(stderr, "Warning: The normalization routine ran into ");
	    fprintf(stderr, "trouble at [%d,%d].\n", i_row, i_col);
	    fprintf(stderr, "The unnormalized kernel value is %g.\n",
		    get_matrix_cell(i_row, i_col, kernel_matrix));
	    fprintf(stderr, "The corresponding diagonals are %g and %g.\n",
		    row_diag, col_diag);
	    fprintf(stderr, "The product of their sqrts is %g * %g = %g.\n",
		    row_diag_sqrt, sqrt(col_diag), product);
	    fprintf(stderr, "The resulting normalized value is %g.\n", 
		    this_cell);
	  }

	  // Set this value manually.
	  if (i_row == i_col) {
	    this_cell = 1.0;
	  } else {
	    this_cell = 0.0;
	  }

	  // Finish the warning.
	  if (warning == FALSE) {
	    fprintf(stderr, "Setting k(%d,%d)=%g.\n", i_row, i_col, this_cell);
	    fprintf(stderr, "This warning will not be repeated if this ");
	    fprintf(stderr, "problem arises again.\n");
	    warning = TRUE;
	  }

	}
      }

      // Store the normalized result.
      set_matrix_cell(i_row, i_col, this_cell, kernel_matrix);
    }
  }

  // Set the self-kernel values to 1.0.
  init_array(1.0, row_self_values);
  init_array(1.0, col_self_values);
}

/*****************************************************************************
 * Given three parameters, A, B and C, compute (B(X + C))^A.
 *****************************************************************************/
static double polynomialize
  (double    power,
   double    coefficient,
   double    constant,
   double    value)
{
  /* Add the given constant. */
  value += constant;

  /* Multiply by the given coefficient. */
  value *= coefficient;

  /* Raise to the given power. */
  value = pow(value, power);

  return(value);
}

/*****************************************************************************
 * Given three parameters, A, B and C, replace each element X in a
 * given array by the value (B(X + C))^A.
 *****************************************************************************/
static void polynomialize_array
  (double    power,
   double    coefficient,
   double    constant,
   ARRAY_T*  array)
{
  int num_items;
  int i_item;
  double value;

  num_items = get_array_length(array);
  for (i_item = 0; i_item < num_items; i_item++) {

    /* Retrieve the value from the kernel matrix. */
    value = get_array_item(i_item, array);

    /* Polynomialize the kernel value. */
    value = polynomialize(power, coefficient, constant, value);

#ifdef DEBUG
    /* Check for NaN. */
    if (isnan(value)) {
      die("Failure to compute polynomial at [?,%d] (%g).\n",
	  i_item, value);
    }
#endif

    /* Store the new value. */
    set_array_item(i_item, value, array);
  }
}


/*****************************************************************************
 * Given three parameters, A, B and C, replace each element X in a
 * given matrix by the value (B(X + C))^A.  Also perform the same
 * operation on two given arrays.
 *****************************************************************************/
static void polynomialize_matrix
  (double    power,
   double    coefficient,
   double    constant,
   ARRAY_T*  row_self_values,
   ARRAY_T*  col_self_values,
   MATRIX_T* matrix)
{
  int num_rows;
  int i_row;

  /* Polynomialize the matrix. */
  num_rows = get_num_rows(matrix);
  for (i_row = 0; i_row < num_rows; i_row++) {
    polynomialize_array(power, coefficient, constant, 
			get_matrix_row(i_row, matrix));
  }

  /* Polynomialize the self-kernel values. */
  if (row_self_values != NULL) {
    polynomialize_array(power, coefficient, constant, row_self_values);
  }
  if (col_self_values != NULL) {
    polynomialize_array(power, coefficient, constant, col_self_values);
  }
}

/****************************************************************************
 * Define the squared distance as follows:
 *
 *     d^2(x,y) = K(x,x) - 2 K(x,y) + K(y,y)
 *
 * This is the squared Euclidean distance between points in feature
 * space. Here is the derivation:
 *
 * ||x-y|| = \sqrt((x_1 - y_1)^2 + (x_2 - y_2)^2 + ... + (x_n - y_n)^2)
 *         = \sqrt((x_1^2 - 2x_1y_1 + y_1^2) + (x_2^2 - 2x_2y_2 + y_2^2)
 *                  + ... + (x_n^2 - 2x_ny_n + y_n^2))
 *         = \sqrt((x_1^2 + x_2^2 + ... + x_n^2) - 2(x_1y_1 + x_2y_2
 *                  + ... + x_ny_n) + (y_1^2 + y_2^2 + ... + y_n^2))
 *         = \sqrt(K(x,x) - 2 K(x,y) + K(y,y))
 *
 ****************************************************************************/
static double compute_squared_distance
  (double Kxx,
   double Kxy,
   double Kyy)
{
  return(Kxx - (2 * Kxy) + Kyy);
}

/****************************************************************************
 * Set the width of a radial basis kernel.
 *
 * If kernel->width is zero, then the kernel->width_factor is used.
 * The width is set to be the median of the distances from each
 * positive example to the nearest negative example, multipled by the
 * width factor.  If no classifications are given, only the median of
 * the distance from each example to the nearest other example is
 * used.
 *
 * If the width is non-zero, then that width is used directly, but the
 * corresponding width_factor is computed and stored in
 * kernel->width_factor.
 * 
 * This function is only called during training, so we know the kernel
 * matrix is square.
 ****************************************************************************/
static double compute_two_squared_width
  (KERNEL_T* kernel,
   CLASS_ARRAY_T* classifications)
{
  int      num_positives;      /* Total number of positive examples. */
  int      num_examples;       /* Number of training examples. */
  int      i_example;          /* Index of this training example. */
  ARRAY_T* nearest_negatives;  /* Distances to nearest negative example. */
  int      i_positive;
  int      i_example1;
  int      i_example2;
  double   squared_distance;   /* Squared distance between two examples. */
  MATRIX_T*      kernel_matrix;
  double   median_distance;
  double   return_value;

  kernel_matrix = kernel->kernel_matrix;
  num_examples = get_num_rows(kernel_matrix);

  if (classifications == NULL) {
    /* If we ain't got no classifications, treat them all as positives. */
    num_positives = num_examples;
  } else {

    /* Count the number of positive examples. */
    num_positives = 0;
    for (i_example = 0; i_example < num_examples; i_example++) {
      if (get_class(i_example, classifications)) {
	num_positives++;
      }
    }
  }

  /* Allocate the array of distances to nearest negatives. */
  nearest_negatives = allocate_array(num_positives);

  /* Find the nearest negative example for each positive. */
  i_positive = -1;
  for (i_example1 = 0; i_example1 < num_examples; i_example1++) {

    /* Consider only positive examples. */
    if ((classifications == NULL) || 
	(get_class(i_example1, classifications))) {
      i_positive++;

      /* Initialize this position to something very large. */
      set_array_item(i_positive, HUGE_VAL, nearest_negatives);

      for (i_example2 = 0; i_example2 < num_examples; i_example2++) {

	/* Don't look at yourself. */
	if (i_example1 == i_example2) {
	  continue;
	}

	/* Consider only negative examples. */
	if ((classifications == NULL) || 
	    (!get_class(i_example2, classifications))) {

	  /* Compute the distance between these examples. */ 
	  squared_distance
	    = compute_squared_distance(get_matrix_cell(i_example1, 
						       i_example1,
						       kernel_matrix),
				       get_matrix_cell(i_example1,
						       i_example2,
						       kernel_matrix),
				       get_matrix_cell(i_example2, 
						       i_example2,
						       kernel_matrix));

	  /* Store the minimum, ignoring distances of zero. */
	  if ((get_array_item(i_positive, nearest_negatives)
	      > squared_distance) && (squared_distance != 0.0)) {
	    set_array_item(i_positive, squared_distance, nearest_negatives);
	  }
	}
      }
    }
  }

  if (RADIAL_DEBUG) {
    fprintf(stderr, "Distances to negatives: ");
    print_array(nearest_negatives, 5, 2, TRUE, stderr);
  }

  // Find the median distance.
  median_distance = compute_median(nearest_negatives);
  free_array(nearest_negatives);

  // Multiply in the given width factor.
  if (kernel->width == 0) {
    return_value = 2 * median_distance * 
      kernel->width_factor * kernel->width_factor;
    kernel->width = sqrt(return_value / 2);
  } 

  // Compute the corresponding width factor.
  else {
    return_value = (kernel->width * kernel->width) * 2;
    kernel->width_factor = sqrt(return_value / (2 * median_distance));
  }

  // If the median minimum distance is zero, fix it.
  if (return_value == 0.0) {
    return_value = SMALL_VALUE;
  }
  
  // Tell the user what's up.
  if (verbosity > NORMAL_VERBOSE) {
    fprintf(stderr, "median distance=%g width=%g width factor=%g\n",
	    median_distance, kernel->width, kernel->width_factor);
  }

  // Return the result.
  return(return_value);
}

/****************************************************************************
 * Compute a radial basis function kernel, defined by
 *
 *     K(x,y) = exp{-d^2(x,y)/2\sigma^2}
 *
 * for some constant \sigma (the width). 
 ****************************************************************************/
static double radial_kernel
  (double    two_squared_width,
   double    Kxx,
   double    Kxy,
   double    Kyy)
{
  double return_value;
  static BOOLEAN_T seen_underflow = FALSE;

  /* Compute the squared distances between the examples. */
  return_value = compute_squared_distance(Kxx, Kxy, Kyy);

  if (RADIAL_DEBUG) {
    fprintf(stderr, "squared distance=%g\n", return_value);
  }

  /* Divide by twice sigma squared. */
  return_value /= two_squared_width;

  if (RADIAL_DEBUG) {
    fprintf(stderr, "unexponentiated=%g\n", return_value);
  }

  /* Exponentiate the opposite. */
  return_value = exp(-return_value);

  /* Make sure we didn't hit zero. */
  if ((return_value == 0.0) && (!seen_underflow) && verbosity > NORMAL_VERBOSE) {
    fprintf(stderr, "Warning: possible underflow in radial_kernel\n");
    seen_underflow = TRUE;
  }

  if (RADIAL_DEBUG) {
    fprintf(stderr, "return_value=%g\n", return_value);
  }
  return(return_value);
}

/*****************************************************************************
 * Convert each value in a given kernel matrix to a radial basis
 * version.
 *****************************************************************************/
static void radialize_matrix
  (double         constant,
   ARRAY_T*       self_rows,
   ARRAY_T*       self_cols,
   double         two_squared_width,
   MATRIX_T*      kernel_matrix)
{
  int    num_rows;
  int    num_cols;
  int    i_row;
  int    i_col;
  double radial_value;

  /* Get the dimensions of the matrix. */
  num_rows = get_num_rows(kernel_matrix);
  assert(num_rows == get_array_length(self_rows));
  num_cols = get_num_cols(kernel_matrix);
  assert(num_cols == get_array_length(self_cols));

  /* Radialize each row of the matrix. */
  for (i_row = 0; i_row < num_rows; i_row++) {
    for (i_col = 0; i_col < num_cols; i_col++) {

      /* Compute the new value. */
      radial_value
	= radial_kernel(two_squared_width,
			get_array_item(i_row, self_rows),
			get_matrix_cell(i_row, i_col, kernel_matrix),
			get_array_item(i_col, self_cols));

      /* Add the constant back in. */
      radial_value += constant;

#ifdef DEBUG
      /* Check for NaN. */
      if (isnan(radial_value)) {
	die("Failure to radialize at [%d,%d] (%g).\n",
	    i_row, i_col, radial_value);
      }
#endif

      /* Store the new value. */
      set_matrix_cell(i_row, i_col, radial_value, kernel_matrix);
    }
  }
}

/***************************************************************************
 * One good way to set the constants to add to the diagonal is
 *
 *   n+ = number of positive examples
 *   n- = number of negative examples
 *   N  = total number of examples
 *   k  = some constant (given by diagonal_factor)
 *   m  = median value of diagonal
 *
 * Then set
 *
 *   positive_diagonal = (n+/N) * k * m
 *   negative_diagonal = (n-/N) * k * m
 *
 ***************************************************************************/
static void set_diagonal_constants
  (CLASS_ARRAY_T* classes, 
   double         diagonal_factor,
   MATRIX_T*      kernel_matrix,
   double*        positive_diagonal, 
   double*        negative_diagonal)
{
  ARRAY_T* self_kernel_values;
  int num_examples;
  int i_example;
  int num_positive;
  int num_negative;
  double median_diagonal;

  /* If the diagonal factor is zero, don't do nothing. */
  if (diagonal_factor == 0) {
    return;
  }

  /* Find the median self-kernel value. */
  self_kernel_values = extract_diagonal(kernel_matrix);
  median_diagonal = compute_median(self_kernel_values);

  /* Count the number of positives and negatives. */
  num_positive = 0;
  num_negative = 0;
  num_examples = get_num_items(classes);
  for (i_example = 0; i_example < num_examples; i_example++) {
    if (get_class(i_example, classes)) {
      num_positive++;
    } else {
      num_negative++;
    }
  }
  assert(num_positive + num_negative == get_array_length(self_kernel_values));
  free_array(self_kernel_values);

  /* Set the positive and negative diagonal constants. */
  *positive_diagonal = ((double)num_positive / (double)num_examples)
    * diagonal_factor * median_diagonal;
  *negative_diagonal = ((double)num_negative / (double)num_examples)
    * diagonal_factor * median_diagonal;

  if (DIAG_DEBUG) {
    fprintf(stderr, "num_positive=%d num_negative=%d ",
	    num_positive, num_negative);
    fprintf(stderr, "diagonal_factor=%g median_diagonal=%g\n",
	    diagonal_factor, median_diagonal);
    fprintf(stderr, "positive_diagonal=%g negative_diagonal=%g\n",
	    *positive_diagonal, *negative_diagonal);
  }
}

/***************************************************************************
 * Add a constant to the diagonal of the kernel matrix.
 *
 * This can be used to accomplish two things:
 * 
 * (1) If the kernel is not positive definite, then adding a
 *     sufficiently large constant to the diagonal will make it so.
 *
 * (2) Adding to the diagonal also effectively scales the weights.  A
 *     larger constant makes the weights smaller.  Adding different
 *     constants for the positives and negatives has the same effect
 *     as placing different constraint ceilings.
 *
 ***************************************************************************/
static void add_to_kernel_diagonal
  (CLASS_ARRAY_T* classes, 
   double*        positive_diagonal, 
   double*        negative_diagonal,
   double         diagonal_factor,
   MATRIX_T*      kernel_matrix)
{
  int num_rows;
  int i_row;

  /* Make sure we've got something to do. */
  if (diagonal_factor == 0.0) {
    return;
  }

  /* Use the diagonal factor to compute the positive and negative
     diagonal constants. */
  set_diagonal_constants(classes, diagonal_factor, kernel_matrix,
			 positive_diagonal, negative_diagonal);

  num_rows = get_num_rows(kernel_matrix);
  for (i_row = 0; i_row < num_rows; i_row++) {
    if (get_class(i_row, classes)) {
      incr_matrix_cell(i_row, i_row, *positive_diagonal, kernel_matrix);
    } else {
      incr_matrix_cell(i_row, i_row, *negative_diagonal, kernel_matrix);
    }
  }
}

/***************************************************************************
 * Do/Undo the kernel diagonal additon. We do this prior to classifying
 * the training examples.
 ***************************************************************************/
void revert_kernel_diagonal
  (CLASS_ARRAY_T* classes, 
   KERNEL_T*      kernel)
{
  int num_rows;
  int i_row;
  double positive_diagonal;
  double negative_diagonal;

  /* Make sure we've got something to do. */
  if (kernel->diagonal_factor == 0.0) {
    return;
  }
  
  /* If the kernel is already modified, undmodify it. If it has not,
     then remodify it */

  positive_diagonal = kernel->positive_diagonal;
  negative_diagonal = kernel->negative_diagonal;

  if (kernel->diagonal_isModified) {

    if (verbosity > NORMAL_VERBOSE)
      fprintf(stderr, "De-modifying the kernel diagonal\n");

    positive_diagonal = -positive_diagonal;
    negative_diagonal = -negative_diagonal;

  } else {

    if (verbosity > NORMAL_VERBOSE)
      fprintf(stderr, "Re-modifying the kernel diagonal\n");

  }

  num_rows = get_num_rows(kernel->kernel_matrix);
  for (i_row = 0; i_row < num_rows; i_row++) {
    if (get_class(i_row, classes)) {
      incr_matrix_cell(i_row, i_row, positive_diagonal, kernel->kernel_matrix);
    } else {
      incr_matrix_cell(i_row, i_row, negative_diagonal, kernel->kernel_matrix);
    }
  }

  kernel->diagonal_isModified = !kernel->diagonal_isModified;

}

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
  (BOOLEAN_T      square, // Is this a square kernel matrix for one data set.
   KERNEL_T*      kernel,
   CLASS_ARRAY_T* classes,
   ARRAY_T*       self_train,
   ARRAY_T*       self_test,
   MATRIX_T*      kernel_matrix)
{
  if (verbosity >= NORMAL_VERBOSE) {
    fprintf(stderr, "Transforming the kernel matrix.\n");
  }

  // Add to the diagonal, if requested. Note that this is not the soft margin
  if ( kernel->add_diag != 0.0 && square) {

    if (verbosity > NORMAL_VERBOSE) {
      fprintf(stderr, "Adding %g to kernel diagonal.\n", kernel->add_diag);
    }

    add_to_diagonal(kernel->add_diag, kernel_matrix);

    scalar_add(kernel->add_diag, self_train);
    scalar_add(kernel->add_diag, self_test);
  }

  /* Normalize the length of each vector. */
  if (kernel->normalize) {
    if (verbosity > NORMAL_VERBOSE) {
      fprintf(stderr, "Normalizing the kernel matrix.\n");
    }
    normalize_kernel_matrix(self_test, self_train, kernel_matrix);
  }

  /* Compute polynomial version of the matrix. */
  if (verbosity > NORMAL_VERBOSE) {
    fprintf(stderr, "Polynomializing the kernel matrix.\n");
  }
  polynomialize_matrix(kernel->power, 
		       kernel->coefficient, 
		       kernel->constant, 
		       self_train, 
		       self_test,
		       kernel_matrix);

  /* Convert to a radial basis kernel. */
  if (kernel->radial) {
    if (verbosity > NORMAL_VERBOSE) {
      fprintf(stderr, "Radializing the kernel matrix.\n");
    }

    /* Compute the kernel width, if it wasn't given. */
    if (kernel->two_squared_width == 0.0) {
      kernel->two_squared_width 
	= compute_two_squared_width(kernel, classes);
    }
    radialize_matrix(kernel->constant, 
		     self_test, 
		     self_train, 
		     kernel->two_squared_width,
		     kernel_matrix);
  }


  /* Add constants to the diagonal. (compute-weights only). 
     Don't do it if we are doing kernelout */
  if ((kernel->diagonal_factor != 0.0) && square  && !kernel->kernel_out) {

    if (verbosity > NORMAL_VERBOSE) {
      fprintf(stderr, "Adding to kernel diagonal.\n");
    }

    add_to_kernel_diagonal(classes, 
			   &(kernel->positive_diagonal),
			   &(kernel->negative_diagonal),
			   kernel->diagonal_factor,
			   kernel_matrix);

    kernel->diagonal_isModified = TRUE;
  }

  /* Subtract the mean from each feature. (Kernel PCA only) */
  if (kernel->zero_mean_col) {
    if (verbosity > NORMAL_VERBOSE) {
      fprintf(stderr, "Subtracting the mean from each feature.\n");
    }
    center_matrix(kernel_matrix, kernel_matrix);
  }
}

/*****************************************************************************
 * Given two input matrices, compute a matrix of pairwise dot products.
 *****************************************************************************/
void compute_base_kernel_matrix
(MATRIX_T* matrix1, // test if in classify 
 MATRIX_T* matrix2, // train
   MATRIX_T* kernel_matrix)
{
  ARRAY_T* example1;
  ARRAY_T* example2;
  int num_rows;
  int num_cols;
  int i_row;
  int i_col;
  double kernel;

  if (verbosity >= NORMAL_VERBOSE) {
    fprintf(stderr, "Computing base kernel matrix.\n");
  }

  /* Compute the size of the output matrix. */
  num_rows = get_num_rows(matrix1);
  assert(get_num_rows(kernel_matrix) == num_rows);
  num_cols = get_num_rows(matrix2);
  assert(get_num_cols(kernel_matrix) == num_cols);
  assert(get_num_cols(matrix1) == get_num_cols(matrix2));

  /* Compute the kernel for every pair of examples. */
  for (i_row = 0; i_row < num_rows; i_row++) {
    example1 = get_matrix_row(i_row, matrix1);

    if (verbosity >= DUMP_VERBOSE) {
      fprintf(stderr, "Computing kernels for example %d.\n", i_row);
    }

    for (i_col = 0; i_col < num_cols; i_col++) {
      example2 = get_matrix_row(i_col, matrix2);

      /* Compute the dot product. */
      kernel = dot_product(example1, example2);

#ifdef DEBUG
      /* Check for NaN. */
      if (isnan(kernel)) {
	die("Failure to compute base kernel at [%d,%d] (%g).\n",
	    i_row, i_col, kernel);
      }
#endif

      /* Store the resulting value in the matrix. */
      set_matrix_cell(i_row, i_col, kernel, kernel_matrix);

      if (KERNEL_DEBUG) {
	fprintf(stderr, "K(%d,%d)=%g\n", i_row, i_col, kernel);
      }
    }
  }
}

/*****************************************************************************
 * Read a kernel operation from an open file.
 *
 * Return NULL on EOF.
 *****************************************************************************/
#define MAX_LINE 100
TRANSFORM_T get_next_transformation
(FILE* kernel_file)
{
  char transform_string[MAX_LINE];
  int  num_scanned;

  num_scanned = fscanf(kernel_file, "%s", transform_string);
  if (num_scanned == EOF) {
    return(INVALID_TRANSFORM);
  } else if (num_scanned == 0) {
    die("Error reading kernel transformations file.\n");
  } else if (strcmp(transform_string, "constant") == 0) {
    return(CONSTANT_TRANSFORM);
  } else if (strcmp(transform_string, "coefficient") == 0) {
    return(COEFFICIENT_TRANSFORM);
  } else if (strcmp(transform_string, "power") == 0) {
    return(POWER_TRANSFORM);
  } else if (strcmp(transform_string, "radial") == 0) {
    return(RADIAL_TRANSFORM);
  } else if (strcmp(transform_string, "diagonal") == 0) {
    return(DIAGONAL_TRANSFORM);
  } else if (strcmp(transform_string, "normalize") == 0) {
    return(NORMALIZE_TRANSFORM);
  } else if (strcmp(transform_string, "center") == 0) {
    return(CENTER_TRANSFORM);
  } else if (strcmp(transform_string, "diffusion") == 0) {
    return(DIFFUSION_TRANSFORM);
  } else {
    die("Invalid transform (%s).\n", transform_string);
  }
  return(0);
}

/*****************************************************************************
 * Read a float from a file.
 *
 * Die on failure.
 *****************************************************************************/
float get_transform_float
(FILE* kernel_file)
{
  float transform_float;
  int   num_scanned;

  num_scanned = fscanf(kernel_file, "%f", &transform_float);
  if ((num_scanned == EOF) || (num_scanned == 0)) {
    die("Error reading kernel parameter.\n");
  }
  return(transform_float);
}


/*****************************************************************************
 * Compute all kernel transformations requested in a given file.
 *****************************************************************************/
void compute_kernel_transformations
(char*     kernel_filename,
 MATRIX_T* kernel_matrix)
{
  FILE* kernel_file = NULL;
  TRANSFORM_T transformation;

  // Allow no kernel operations.
  if (kernel_filename == NULL) { 
    return;
  }

  // Open the kernel transformations file.
  if (!open_file(kernel_filename, "r", TRUE, "kernel operations",
		 "kernel operations", &kernel_file)) {
    exit(1);
  }

  // Do each transform.
  transformation = get_next_transformation(kernel_file);
  while (transformation != INVALID_TRANSFORM) {
    float parameter;

    if (transformation == CONSTANT_TRANSFORM) {
      parameter = get_transform_float(kernel_file);
      if (verbosity >= NORMAL_VERBOSE) {
	fprintf(stderr, "Adding %g to kernel.\n", parameter);
      }
      scalar_add_matrix(parameter, kernel_matrix);
    }

    else if (transformation == COEFFICIENT_TRANSFORM) {
      parameter = get_transform_float(kernel_file);
      if (verbosity >= NORMAL_VERBOSE) {
	fprintf(stderr, "Multiply kernel by %g.\n", parameter);
      }
      scalar_mult_matrix(parameter, kernel_matrix);
    }

    else if (transformation == POWER_TRANSFORM) {
      parameter = get_transform_float(kernel_file);
      if (verbosity >= NORMAL_VERBOSE) {
	fprintf(stderr, "Raising kernel to %g power.\n", parameter);
      }
      power_matrix(parameter, kernel_matrix);
    }

    else if (transformation == RADIAL_TRANSFORM) {
      ARRAY_T* diagonal = extract_diagonal(kernel_matrix);
      parameter = get_transform_float(kernel_file);
      if (verbosity >= NORMAL_VERBOSE) {
	fprintf(stderr, "Radializing kernel with width %g.\n", parameter);
      }
      radialize_matrix(0, diagonal, diagonal, 2 * parameter * parameter,
		       kernel_matrix);
      free_array(diagonal);
    }

    else if (transformation == DIAGONAL_TRANSFORM) {
      parameter = get_transform_float(kernel_file);
      if (verbosity >= NORMAL_VERBOSE) {
	fprintf(stderr, "Adding %g to kernel diagonal.\n", parameter);
      }
      add_to_diagonal(parameter, kernel_matrix);
    }

    else if (transformation == NORMALIZE_TRANSFORM) {
      if (verbosity >= NORMAL_VERBOSE) {
	fprintf(stderr, "Normalizing kernel.\n");
      }
      ARRAY_T* diagonal = extract_diagonal(kernel_matrix);
      normalize_kernel_matrix(diagonal, diagonal, kernel_matrix);
      free_array(diagonal);
    }

    else if (transformation == CENTER_TRANSFORM) {
      if (verbosity >= NORMAL_VERBOSE) {
	fprintf(stderr, "Centering kernel.\n");
      }
      center_matrix(kernel_matrix, kernel_matrix);
    } else if (transformation == DIFFUSION_TRANSFORM) {
      parameter = get_transform_float(kernel_file);

      if (verbosity >= NORMAL_VERBOSE) {
	fprintf(stderr, "Converting kernel to distance matrix.\n");
      }

      // Convert from a kernel to a distance.
      MATRIX_T* distance_matrix = kernel_to_distance(kernel_matrix);

      if (verbosity >= NORMAL_VERBOSE) {
	fprintf(stderr, "Computing diffusion kernel with constant %g.\n",
		parameter);
      }
      // Do the diffusion.
      perform_diffusion(parameter, distance_matrix);

      // Copy the resulting kernel matrix.
      copy_matrix(distance_matrix, kernel_matrix);
      free_matrix(distance_matrix);
      
    } else {
      die("Invalid transformation.\n");
    }

    // Make sure the kernel is still OK.
    if (1) {
      // FIXME: Make this a function.
      int num_rows = get_num_rows(kernel_matrix);
      int num_cols = get_num_cols(kernel_matrix);
      int i_row;
      int i_col;

      for (i_row = 0; i_row < num_rows; i_row++) {
	for (i_col = 0; i_col < num_cols; i_col++) {
	  MTYPE this_cell = get_matrix_cell(i_row, i_col, kernel_matrix);
	  if (isnan(this_cell)) {
	    die("NaN at [%d,%d].\n", i_row, i_col);
	  }
	}
      }
    }

    transformation = get_next_transformation(kernel_file);
  }
}

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */


