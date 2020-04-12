/**************************************************************************
 * FILE: linear-algebra.c
 * AUTHOR: William Stafford Noble
 * CREATE DATE: 2/20/1999
 * PROJECT: SVM
 * COPYRIGHT: 1999-2001, Columbia University (see ../doc/copyright.html)
 * VERSION: $Revision: 1.7 $
 * DESCRIPTION: Linear algebra routines for matrices.
 **************************************************************************/
#include "linear-algebra.h"
#include "matrix.h"
#include "array.h"
#include "utils.h"
#include "eigens.h"
#include <math.h>
#include <assert.h>

#ifndef EIGEN_DEBUG
#define EIGEN_DEBUG 0
#endif

/***********************************************************************
 * Turn a given square matrix into an identity matrix.
 ***********************************************************************/
void make_identity_matrix
  (MATRIX_T* matrix)
{
  int num_rows;
  int num_cols;
  int i_row;
  int i_col;

  num_rows = get_num_rows(matrix);
  num_cols = get_num_cols(matrix);

  /* Make sure the matrix is square. */
  if (num_rows != num_cols) {
    die("Attempted to invert a non-square matrix (%d != %d).", 
	num_rows, num_cols);
  }

  /* Put 1s on the diagonal and 0s everywhere else. */
  for (i_row = 0; i_row < num_rows; i_row++) {
    for (i_col = 0; i_col < num_cols; i_col++) {
      if (i_row == i_col) {
	set_matrix_cell(i_row, i_col, (MTYPE)1, matrix);
      } else {
	set_matrix_cell(i_row, i_col, (MTYPE)0, matrix);
      }
    }
  }
}

/***********************************************************************
 * Transpose a matrix.
 ***********************************************************************/
MATRIX_T* transpose_matrix
  (MATRIX_T* matrix)
{
  MATRIX_T* new_matrix;
  int num_rows = get_num_rows(matrix);
  int num_cols = get_num_cols(matrix);
  int i_row;
  int i_col;

  /* Allocate the new matrix. */
  new_matrix = allocate_matrix(num_cols, num_rows);

  /* Copy elements from the old to the new matrix. */
  for (i_row = 0; i_row < num_rows; i_row++) {
    for (i_col = 0; i_col < num_cols; i_col++) {
      set_matrix_cell(i_col, i_row, get_matrix_cell(i_row, i_col, matrix),
		      new_matrix);
    }
  }

  /* Return the new matrix. */
  return(new_matrix);
}


/***********************************************************************
 * Helper function for the Gauss Jordan method of matrix inversion.
 *
 * Multiply a given array by a scalar and add the resulting array to a
 * given array.
 *
 * Assumes that the arrays are of the same length.
 ***********************************************************************/
static void add_array_product
  (ARRAY_T* source_array,
   MTYPE    multiplier,
   ARRAY_T* target_array)
{
  int i_item;
  int num_items = get_array_length(source_array);

  for (i_item = 0; i_item < num_items; i_item++) {
    set_array_item(i_item, get_array_item(i_item, target_array) + 
		   (multiplier * get_array_item(i_item, source_array)),
		   target_array);
  }
}


/***********************************************************************
 * Select from a given matrix a row that has a non-zero element in a
 * given column and a given range of rows.
 *
 * The parameters i_start and i_end tell which row to begin at and
 * which row to end at.  They must be non-equal, but it need not be
 * the case that i_start < i_end.
 * 
 * Returns the index of the selected row or -1 if there is no such
 * row.
 ***********************************************************************/
int select_other_row
  (int       i_start,
   int       i_end,
   int       i_col,
   MATRIX_T* matrix)
{
  int i_row;

  i_row = i_start;
  while (i_row != i_end) {
    if (!almost_equal(get_matrix_cell(i_row, i_col, matrix), 0.0, 
		      MATRIX_PRECISION)) {
      return(i_row);
    }

    /* We may be moving either up or down the matrix. */
    if (i_row < i_end) {
      i_row++;
    } else {
      i_row--;
    }
  }
  return(-1);
}

/***********************************************************************
 * Helper function for the Gauss Jordan method of matrix inversion.
 *
 * Returns a Boolean indicating success.
 ***********************************************************************/
static BOOLEAN_T make_lower_triangle_zeroes
  (MATRIX_T* matrix,
   MATRIX_T* inverse_matrix)
{
  int num_cols = get_num_cols(matrix);
  int i_col;
  int i_row;
  int i_prev_row;
  MTYPE prev_row_element;
  MTYPE this_row_element;
  MTYPE multiplier;

  for (i_col = 0; i_col < num_cols - 1; i_col++) {
    for (i_row = i_col + 1; i_row < num_cols; i_row++) {

      /* Select a row above this one that contains a non-zero element
         in this column. */
      i_prev_row = select_other_row(i_col, i_row, i_col, matrix);

      /* Punt if the matrix is not invertible. */
      if (i_prev_row == -1) {
	return(FALSE);
      }

      /* Get the two elements that are to be compared. */
      prev_row_element = get_matrix_cell(i_prev_row, i_col, matrix);
      this_row_element = get_matrix_cell(i_row, i_col, matrix);

      /* What do we need to multiply the previous row by to get a zero
         in this row? */
      multiplier = -this_row_element / prev_row_element;

      /* Multiply the previous row and add it to this row. */
      add_array_product(get_matrix_row(i_prev_row, matrix), multiplier, 
			get_matrix_row(i_row, matrix));

      /* Repeat this procedure for the inverse matrix. */
      add_array_product(get_matrix_row(i_prev_row, inverse_matrix), 
			multiplier, get_matrix_row(i_row, inverse_matrix));
    }
  }
  return(TRUE);
}


/***********************************************************************
 * Helper function for the Gauss Jordan method of matrix inversion.
 *
 * Returns a Boolean indicating success.
 ***********************************************************************/
static BOOLEAN_T make_upper_triangle_zeroes
  (MATRIX_T* matrix,
   MATRIX_T* inverse_matrix)
{
  int num_cols = get_num_cols(matrix);
  int i_col;
  int i_row;
  int i_next_row;
  MTYPE next_row_element;
  MTYPE this_row_element;
  MTYPE multiplier;

  for (i_col = num_cols - 1; i_col > 0; i_col--) {
    for (i_row = 0; i_row < i_col; i_row++) {

      /* Select a row below this one that contains a non-zero element
         in this column. */
      i_next_row = select_other_row(i_col, i_row, i_col, matrix);

      /* Punt if the matrix is not invertible. */
      if (i_next_row == -1) {
	return(FALSE);
      }

      /* Get the two elements that are to be compared. */
      next_row_element = get_matrix_cell(i_next_row, i_col, matrix);
      this_row_element = get_matrix_cell(i_row, i_col, matrix);

      /* What do we need to multiply the first row by to get a zero in
         this row? */
      multiplier = -this_row_element / next_row_element;

      /* Multiply the next row and add it to this row. */
      add_array_product(get_matrix_row(i_next_row, matrix), multiplier, 
			get_matrix_row(i_row, matrix));

      /* Repeat this procedure for the inverse matrix. */
      add_array_product(get_matrix_row(i_next_row, inverse_matrix), 
			multiplier, get_matrix_row(i_row, inverse_matrix));
    }
  }
  return(TRUE);
}


/***********************************************************************
 * Invert a matrix.
 ***********************************************************************/
MATRIX_T* invert_matrix
  (MATRIX_T* matrix)
{
  int num_rows = get_num_rows(matrix);
  int num_cols = get_num_cols(matrix);
  int i_row;
  MTYPE diagonal_element;
  MATRIX_T* matrix_copy;
  MATRIX_T* inverse_matrix;

  /* Make sure the matrix is square. */
  if (num_rows != num_cols) {
    die("Attempted to invert a non-square matrix (%d != %d).", 
	num_rows, num_cols);
  }

  /* Make a local copy of the matrix. */
  matrix_copy = allocate_matrix(num_rows, num_cols);
  copy_matrix(matrix, matrix_copy);

  /* Allocate the inverse matrix. */
  inverse_matrix = allocate_matrix(num_rows, num_cols);

  /* Start the new matrix as the identity matrix. */
  make_identity_matrix(inverse_matrix);

  /* Make the lower triangle of zeroes. */
  if (!make_lower_triangle_zeroes(matrix_copy, inverse_matrix)) {
    free_matrix(inverse_matrix);
    if (verbosity >= NORMAL_VERBOSE) {
      fprintf(stderr, "Matrix is not invertible.\n");
    }
    return(NULL);
  }

  /* Make the upper triangle of zeroes. */
  if (!make_upper_triangle_zeroes(matrix_copy, inverse_matrix)) {
    free_matrix(inverse_matrix);
    if (verbosity >= NORMAL_VERBOSE) {
      fprintf(stderr, "Matrix is not invertible.\n");
    }
    return(NULL);
  }

  /* Divide each row of the inverted matrix by the non-zero element of
     the corresponding row in the original matrix. */
  for (i_row = 0; i_row < num_rows; i_row++) {
    diagonal_element = get_matrix_cell(i_row, i_row, matrix_copy);
    scalar_mult(1.0 / diagonal_element, get_matrix_row(i_row, inverse_matrix));
  }

  free_matrix(matrix_copy);
  return(inverse_matrix);
}

/***********************************************************************
 * Compute the covariance matrix of a given matrix.
 ***********************************************************************/
MATRIX_T* compute_covariance
  (MATRIX_T* matrix) 
{
  int       i_col;
  int       num_cols;
  ARRAY_T*  this_column;
  MTYPE     factor;
  MATRIX_T* transposed;
  MATRIX_T* return_matrix;

  /* Subtract the mean from each column in the matrix. */
  num_cols = get_num_cols(matrix);
  for (i_col = 0; i_col < num_cols; i_col++) {
    this_column = get_matrix_column(i_col, matrix);
    sum_to_zero(this_column);
    set_matrix_column(this_column, i_col, matrix);
    free_array(this_column);
  }

  /* Get the transpose of the matrix. */
  transposed = transpose_matrix(matrix);

  /* Multiply the given matrix by its transpose. */
  return_matrix = matrix_multiply(transposed, matrix, NULL);
  free_matrix(transposed);

  /* Divide by the number of rows. */
  factor = (MTYPE)get_num_rows(matrix) - 1.0;
  factor = 1.0 / factor;
  scalar_mult_matrix(factor, return_matrix);

  /* Return the computed matrix. */
  return(return_matrix);
}


/***********************************************************************
 * Compute all correlations between pairs of rows or columns in two 
 * matrices.  Allocates and returns a new matrix.
 ***********************************************************************/
MATRIX_T* compute_pairwise_correlations
(BOOLEAN_T use_rows, // else use columns
 MATRIX_T* matrix1,
 MATRIX_T* matrix2) 
{
  int       num_target_rows;
  int       num_target_cols;
  int       i_target_row;
  int       i_target_col;
  ARRAY_T*  this_row;
  ARRAY_T*  this_col;
  MATRIX_T* return_matrix;

  // Allocate a matrix of the proper dimensionality.
  if (use_rows) {
    num_target_rows = get_num_rows(matrix1);
    num_target_cols = get_num_rows(matrix2);
    if (get_num_cols(matrix1) != get_num_cols(matrix2)) {
      die("Trying to compute row correlations for matrices with different numbers of columns (%d != %d).\n",
	  get_num_cols(matrix1), get_num_cols(matrix2));
    }
  } else {
    num_target_rows = get_num_cols(matrix1);
    num_target_cols = get_num_cols(matrix2);
    if (get_num_rows(matrix1) != get_num_rows(matrix2)) {
      die("Trying to compute column correlations for matrices with different numbers of rows (%d != %d).\n",
	  get_num_rows(matrix1), get_num_rows(matrix2));
    }
  }
  return_matrix = allocate_matrix(num_target_rows, num_target_cols);

  // Traverse the entire output matrix.
  for (i_target_row = 0; i_target_row < num_target_rows; i_target_row++) {

    if (use_rows) {
      this_row = get_matrix_row(i_target_row, matrix1);
    } else {
      this_row = get_matrix_column(i_target_row, matrix1);
    }
      
    for (i_target_col = 0; i_target_col < num_target_cols; i_target_col++) {

      if (use_rows) {
	this_col = get_matrix_row(i_target_col, matrix2);
      } else {
	this_col = get_matrix_column(i_target_col, matrix2);
      }

      // Compute and store the correlation.
      set_matrix_cell(i_target_row, i_target_col, 
		      correlation(this_row, this_col), return_matrix);

      // Free memory if this is a column.
      if (!use_rows) {
	free_array(this_col);
      }
    }

    // Free memory if this is a column.
    if (!use_rows) {
      free_array(this_row);
    }    
  }

  // Return the computed matrix.
  return(return_matrix);
}

/***********************************************************************
 * Compute Jacard similarity between pairs of rows in a matrix.
 ***********************************************************************/
MATRIX_T* compute_pairwise_jacard
  (MATRIX_T* matrix) 
{
  int       num_rows;
  int       i_row;
  int       i_col;
  ARRAY_T*  this_row;
  ARRAY_T*  this_col;
  MATRIX_T* return_matrix;

  /* Allocate a matrix of the proper dimensionality. */
  num_rows = get_num_rows(matrix);
  return_matrix = allocate_matrix(num_rows, num_rows);

  for (i_row = 0; i_row < num_rows; i_row++) {
    this_row = get_matrix_row(i_row, matrix);
    for (i_col = 0; i_col < num_rows; i_col++) {
      this_col = get_matrix_row(i_col, matrix);

      /* Compute and store the correlation. */
      set_matrix_cell(i_row, i_col, jacard(this_row, this_col), 
		      return_matrix);
    }
  }

  /* Return the computed matrix. */
  return(return_matrix);
}


/***********************************************************************
 * Multiply two matrices to get a third.
 *
 * If matrix3 is non-NULL, then it is used to store the output.
 ***********************************************************************/
MATRIX_T* matrix_multiply
  (MATRIX_T* matrix1,
   MATRIX_T* matrix2,
   MATRIX_T* matrix3)
{
  int i_row1;
  int i_col2;
  int i_common;
  MATRIX_T* new_matrix;
  MTYPE this_cell;
  MTYPE first_value;
  MTYPE second_value;

  /* Get the dimensions of both matrices. */
  int num_rows1 = get_num_rows(matrix1);
  int num_cols1 = get_num_cols(matrix1);
  int num_rows2 = get_num_rows(matrix2);
  int num_cols2 = get_num_cols(matrix2);

  /* Make sure the dimensions are OK. */
  if (num_cols1 != num_rows2) {
    die("Attempted to multiply matrices with incompatible dimensions (%d != %d).\n",
	num_cols1, num_rows2);
  }

  /* Tell the user what we're doing. */
  if ((verbosity > NORMAL_VERBOSE) &&
      (num_rows1 * num_cols1 * num_cols2 > 10000000)) {
    fprintf(stderr, "Begin matrix multiplication (%d x %d times %d x %d).\n",
	    num_rows1, num_cols1, num_rows2, num_cols2);
  }
	    
  /* Allocate the new matrix. */
  if (matrix3 == NULL) {
    new_matrix = allocate_matrix(num_rows1, num_cols2);
  } else {
    new_matrix = matrix3;
    if ((get_num_rows(matrix3) != num_rows1) ||
	(get_num_cols(matrix3) != num_cols2)) {
      die("Expected %d by %d matrix for matrix_multiply, got %d by %d.\n",
	  num_rows1, num_cols2, get_num_rows(matrix3),
	  get_num_cols(matrix3));
    }
  }

  /* Fill in the new matrix cell by cell. */
  for (i_row1 = 0; i_row1 < num_rows1; i_row1++) {
    if ((verbosity > HIGH_VERBOSE) &&
	(num_rows1 * num_cols1 * num_cols2 > 10000000)) {
      fprintf(stderr, "Matrix multiply: row %d.\n", i_row1);
    }

    for (i_col2 = 0; i_col2 < num_cols2; i_col2++) {
      
      /* Each cell is the sum of products. */
      this_cell = 0.0;
      for (i_common = 0; i_common < num_rows2; i_common++) {
	first_value = get_matrix_cell(i_row1, i_common, matrix1);
	second_value = get_matrix_cell(i_common, i_col2, matrix2);
	this_cell += first_value * second_value;
      }
      set_matrix_cell(i_row1, i_col2, this_cell, new_matrix);
    }
  }
  return(new_matrix);
}

#ifdef VERIFY_EIGENS
/****************************************************************************
 * If these are the correct eigenvectors and eigenvalues, then
 *
 *     matrix x eigenvectors = eigenvalues x eigenvectors.
 *****************************************************************************/
static void verify_eigens
  (MATRIX_T* matrix,
   MATRIX_T* eigenvectors,
   ARRAY_T*  eigenvalues)
{
  int num_rows;
  int i_row;
  MATRIX_T* eigenvalue_matrix;
  MATRIX_T* left_side;
  MATRIX_T* right_side;

  /* Put the eigenvalues along the diagonal of a matrix. */
  num_rows = get_num_rows(matrix);
  eigenvalue_matrix = allocate_matrix(num_rows, num_rows);
  init_matrix(0.0, eigenvalue_matrix);
  for (i_row = 0; i_row < num_rows; i_row++) {
    set_matrix_cell(i_row, i_row, get_array_item(i_row, eigenvalues),
		    eigenvalue_matrix);
  }

  /* Do the two matrix multiplies. */
  left_side = matrix_multiply(matrix, eigenvectors, NULL);
  right_side = matrix_multiply(eigenvectors, eigenvalue_matrix, NULL);

  /* Compare them. */
  if (!equal_matrices(MATRIX_PRECISION * 100.0, left_side, right_side)) {

    /* Print them. */
    fprintf(stderr, "Eigenvectors:\n");
    print_matrix(eigenvectors, 6, 3, FALSE, stderr);
    fprintf(stderr, "Eigenvalues:\n");
    print_array(eigenvalues, 6, 3, TRUE, stderr);
    fprintf(stderr, "Left side:\n");
    print_matrix(left_side, 6, 3, FALSE, stderr);
    fprintf(stderr, "Right side:\n");
    print_matrix(right_side, 6, 3, FALSE, stderr);

    die("Incorrect eigenvalue calculation.");
  }

  free_matrix(eigenvalue_matrix);
}
#endif

/**************************************************************************
 * Convert a symmetric matrix to an array in lower triangle form.
 * This is required by the find_symmetric_eigenvectors routine.
 **************************************************************************/
static void  matrix_to_lower_triangle
  (MATRIX_T* matrix,
   double*   lower_triangle)
{
  int i_row;
  int i_col;
  int matrix_size;
  int array_index;

  matrix_size = get_num_rows(matrix);
  array_index = 0;
  for (i_row = 0; i_row < matrix_size; i_row++) {
    for (i_col = 0; i_col <= i_row; i_col++) {
      lower_triangle[array_index] = get_matrix_cell(i_row, i_col, matrix);
      array_index++;
    }
  }
  assert(2 * array_index == matrix_size * (matrix_size + 1));
}

/**************************************************************************
 * Compute the eigenvectors and the eigenvalues of a symmetric matrix.
 *
 * The eigenvectors are stored as columns, and sorted by decreasing
 * eigenvalue.
 **************************************************************************/
void find_symmetric_eigenvectors
  (MATRIX_T*  matrix,
   ARRAY_T**  eigenvalues,
   MATRIX_T** eigenvectors)
{
  int      matrix_size;
  int      triangle_size;
  double*  lower_triangle;
  double*  raw_eigenvalues;
  double*  raw_eigenvectors;

  /* Verify that the matrix is symmetric. */
  if (!is_symmetric(TRUE, MATRIX_PRECISION, matrix)) {
    fprintf(stderr, "Warning: Finding the eigenvectors of a ");
    fprintf(stderr, "non-symmetric matrix.\n");
  }

  /* Convert the lower triangle of the given matrix to an array. */
  matrix_size = get_num_rows(matrix);
  triangle_size = (int)((double)(matrix_size * (matrix_size + 1)) / 2.0);
  lower_triangle = (double*)calloc(triangle_size, sizeof(double));
  matrix_to_lower_triangle(matrix, lower_triangle);

  /* Allocate memory for the raw eigenvalues and eigenvectors. */
  raw_eigenvalues = (double*)calloc(matrix_size, sizeof(double));
  raw_eigenvectors = (double*)calloc(matrix_size * matrix_size,
				     sizeof(double));

  /* Call the eigenvector routine. */
  eigens(lower_triangle, raw_eigenvectors, raw_eigenvalues, matrix_size);
  myfree(lower_triangle);
  
  /* Object-ify the raw values. */
  *eigenvalues = allocate_array(matrix_size);
  fill_array(raw_eigenvalues, *eigenvalues);
  myfree(raw_eigenvalues);

  *eigenvectors = allocate_matrix(matrix_size, matrix_size);
  fill_matrix(raw_eigenvectors, *eigenvectors);
  myfree(raw_eigenvectors);

  /* Sort the eigenvectors and the eigenvalues. */
  if (verbosity > NORMAL_VERBOSE) {
    fprintf(stderr, "Sorting eigenvectors.\n");
    print_array(*eigenvalues, 6, 3, TRUE, stderr);
  }
  sort_matrix_cols(TRUE, *eigenvalues, *eigenvectors);
  if (verbosity > NORMAL_VERBOSE) {
    print_array(*eigenvalues, 6, 3, TRUE, stderr);
  }

#ifdef VERIFY_EIGENS
  verify_eigens(matrix, *eigenvectors, *eigenvalues);
#endif
}

/**************************************************************************
 * Compute the first eigenvector of a given matrix.
 **************************************************************************/
#define CONVERGENCE 1E-10
#define MAXITER 1E6

ARRAY_T* find_first_eigenvector
  (MATRIX_T*  matrix)
{
  double    delta;
  int       iter;
  double    scale_factor;
  MATRIX_T* eigenvector;
  MATRIX_T* new_eigenvector = NULL;
  ARRAY_T*  return_array;    /* The eigenvector stored as an array. */

  /* Randomly guess an eigenvector. */
  eigenvector = allocate_matrix(get_num_cols(matrix), 1);
  randomize_matrix(1.0, eigenvector);

  if (EIGEN_DEBUG) {
    fprintf(stderr, "Initial eigenvector:\n");
    print_matrix(eigenvector, 5, 3, FALSE, stderr);
  }

  delta = CONVERGENCE;
  iter = 0;
  while (delta >= CONVERGENCE) {

    /* Multiply the matrix by the guess. */
    new_eigenvector = matrix_multiply(matrix, eigenvector, NULL);

    if (EIGEN_DEBUG) {
      fprintf(stderr, "Re-estimated eigenvector:\n");
      print_matrix(new_eigenvector, 5, 3, FALSE, stderr);
    }

    /* Scale so that the vector has length 1. */
    scale_factor = sqrt(sum_of_squares_matrix(new_eigenvector));
    assert(scale_factor != 0.0);    
    scalar_mult_matrix(1.0 / scale_factor, new_eigenvector);

    if (EIGEN_DEBUG) {
      fprintf(stderr, "Re-estimated eigenvector scaled by %g:\n", 
	      scale_factor);
      print_matrix(new_eigenvector, 5, 3, FALSE, stderr);
    }

    /* Compute the change. */
    delta = sum_of_square_diff_matrices(new_eigenvector, eigenvector);

    /* Replace the guess with the result. */
    free_matrix(eigenvector);
    eigenvector = new_eigenvector;

    /* Die if we've reached the maximum number of iterations. */
    if (iter > MAXITER) {
      die("Can't find the eigenvector in %d iterations.  Perhaps the %s",
	  MAXITER, "ratio of the largest and smallest eigenvalues is 1.0.\n");
    }
    iter++;
  }

  return_array = get_matrix_column(0, new_eigenvector);
  free_matrix(new_eigenvector);
  return(return_array);
}

/**************************************************************************
 * Given a matrix and an eigenvector, find the corresponding eigenvalue.
 *
 * This could be made more precise by repeating the operation across
 * the whole eigenvector and then averaging.
 **************************************************************************/
double find_corresponding_eigenvalue
  (ARRAY_T* eigenvector,
   MATRIX_T* matrix)
{
  ARRAY_T* first_row;
  double   product;
  double   return_value;

  /* Get the first row of the matrix. */
  first_row = get_matrix_row(0, matrix);

  /* Multiply with the eigenvector. */
  product = dot_product(first_row, eigenvector);

  /* Divide by the first element of the eigenvector. */
  return_value = product / get_array_item(0, eigenvector);

  return(return_value);
}

/**************************************************************************
 * Compute the last eigenvalue of a given matrix.
 **************************************************************************/
double find_last_eigenvalue
  (MATRIX_T*  matrix)
{
  MATRIX_T* eigenvectors;
  ARRAY_T*  eigenvalues;
  int       num_eigenvalues;
  double    last_eigenvalue;

  /* Compute all the eigenvalues. */
  find_symmetric_eigenvectors(matrix, &eigenvalues, &eigenvectors);

  /* Return the last one. */
  num_eigenvalues = get_array_length(eigenvalues);
  last_eigenvalue = get_array_item(num_eigenvalues - 1, eigenvalues);

  free_array(eigenvalues);
  free_matrix(eigenvectors);
  return(last_eigenvalue);

#ifdef DOES_NOT_WORK
  MATRIX_T* inverse_matrix;
  ARRAY_T*  inverse_eigenvector;
  double    inverse_eigenvalue;
  double    last_eigenvalue;

  /* Compute the inverse of the matrix. */
  inverse_matrix = invert_matrix(matrix);

  /* If the matrix is singular, return 0. */
  if (inverse_matrix == NULL) {
    fprintf(stderr, "Warning: Finding eigenvalue of a singular matrix.\n");
    return(0.0);
  }

  /* Find the first eigenvector of the inverted matrix. */
  inverse_eigenvector = find_first_eigenvector(inverse_matrix);

  /* Find the corresponding eigenvalue. */
  inverse_eigenvalue = find_corresponding_eigenvalue(inverse_eigenvector,
						     inverse_matrix);

  /* The first eigenvalue of the inverse matrix is the reciprocal
     of the smallest eigenvalue of the original matrix. */
  last_eigenvalue = 1.0 / inverse_eigenvalue;

  free_matrix(inverse_matrix);
  free_array(inverse_eigenvector);
  return(last_eigenvalue);
#endif
}

/****************************************************************************
 * Force a square matrix to be positive definite by adding to the diagonal.
 ****************************************************************************/
#define EXTRA_LITTLE_BIT 0.000001
void make_positive_definite
  (MATRIX_T* matrix)
{
  double value;

  /* Find the last eigenvalue. */
  value = find_last_eigenvalue(matrix);
  fprintf(stderr, "last eigenvalue=%g\n", value);
  

  /* If the smallest eigenvalue value is negative, add its opposite
     to the diagonal. */
  if (value < 0.0) {
    add_to_diagonal(-value * (1.0 + EXTRA_LITTLE_BIT), matrix);
  }
}

/*****************************************************************************
 * Convert a kernel matrix to a Euclidean distance matrix, or vice versa
 *****************************************************************************/
MATRIX_T* kernel_to_distance
(MATRIX_T* matrix)
{
  MATRIX_T* new_matrix;
  int num_rows = get_num_rows(matrix);
  int num_cols = get_num_cols(matrix);
  int i_row;
  int i_col;

  // Make sure the matrix is square.
  if (num_rows != num_cols) {
    die("Attempted to convert a non-square matrix (%d != %d).", 
	num_rows, num_cols);
  }

  // Allocate the new matrix.
  new_matrix = allocate_matrix(num_cols, num_rows);

  // Compute each distance.
  for (i_row = 0; i_row < num_rows; i_row++) {
    MTYPE row_diag = get_matrix_cell(i_row, i_row, matrix);
    for (i_col = 0; i_col < num_cols; i_col++) {
      MTYPE col_diag = get_matrix_cell(i_col, i_col, matrix);
      MTYPE kernel = get_matrix_cell(i_row, i_col, matrix);
      MTYPE distance = sqrt(row_diag + col_diag - (2 * kernel));
      set_matrix_cell(i_row, i_col, distance, new_matrix);
    }
  }
  
  return(new_matrix);
}

/*****************************************************************************
 * Compute the intrinsic dimensionality of a Euclidean distance
 * matrix, which is defined as the mean squared inter-object distance,
 * divided by the corresponding variance.
 *****************************************************************************/
double get_intrinsic_dimensionality
(MATRIX_T* matrix)
{
  int i_row;
  int i_col;

  // Make sure the matrix is square.
  int num_rows = get_num_rows(matrix);
  int num_cols = get_num_cols(matrix);
  if (num_rows != num_cols) {
    die("Attempted to convert a non-square matrix (%d != %d).", 
	num_rows, num_cols);
  }

  // Compute the average value.
  MTYPE average = 0.0;
  for (i_row = 0; i_row < num_rows; i_row++) {
    for (i_col = 0; i_col < num_cols; i_col++) {
      average += get_matrix_cell(i_row, i_col, matrix);
    }
  }
  average /= (MTYPE)(num_rows * num_cols);
  //fprintf(stderr, "Mean inter-object distance = %g\n", average);

  // Compute the total error.
  MTYPE total_error = 0.0;
  for (i_row = 0; i_row < num_rows; i_row++) {
    for (i_col = 0; i_col < num_cols; i_col++) {
      MTYPE error = get_matrix_cell(i_row, i_col, matrix) - average;
      total_error += error * error;
    }
  }

  // Compute the variance.
  MTYPE variance = total_error / (MTYPE)((num_rows * num_cols) - 1);
  //fprintf(stderr, "Variance on inter-object distance = %g\n", variance);

  // Compute the intrinsic dimensionality.
  return((average * average) / variance);
}

/*****************************************************************************
 * Perform a diffusion on a given symmetrix matrix.
 *****************************************************************************/
void perform_diffusion
(float     diffusion_constant,
 MATRIX_T* matrix)
{
  // Compute the laplacian of the kernel.
  laplacian(matrix);

  // Multiply by the diffusion constant.
  scalar_mult_matrix(diffusion_constant, matrix);

  // Matrix exponentiate.
  matrix_exponential(matrix);
}

/****************************************************************************
 * Sort a matrix by column, according to a given set of sort keys.
 * This is in here, rather than matrix.c, because it requires the
 * transpose function.
 ****************************************************************************/
void sort_matrix_cols
  (BOOLEAN_T reverse_sort,
   ARRAY_T*  keys,
   MATRIX_T* matrix)
{
  MATRIX_T* transposed_matrix;
  MATRIX_T* sorted_matrix;

  /* This function is an ugly, ugly memory hog. */
  transposed_matrix = transpose_matrix(matrix);
  sort_matrix_rows(reverse_sort, keys, transposed_matrix);
  sorted_matrix = transpose_matrix(transposed_matrix);
  free_matrix(transposed_matrix);
  copy_matrix(sorted_matrix, matrix);
  free_matrix(sorted_matrix);
}

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
(MATRIX_T* matrix)
{
  ARRAY_T*  eigenvalues;
  MATRIX_T* eigenvectors;
  MATRIX_T* transposed_eigenvectors;

  // Compute the eigenvalues.
  find_symmetric_eigenvectors(matrix,
			      &eigenvalues,
			      &eigenvectors);

  // Exponentiate the eigenvalues.
  exp_array(eigenvalues);

  // Transpose the eigenvectors.
  transposed_eigenvectors = transpose_matrix(eigenvectors);

  // Multiply the original eigenvectors by the eigenvalues.
  int num_rows = get_num_rows(eigenvectors);
  int i_row;
  for (i_row = 0; i_row < num_rows; i_row++) {
    scalar_mult(get_array_item(i_row, eigenvalues),
		get_matrix_row(i_row, transposed_eigenvectors));
  }
  free_array(eigenvalues);

  // Multiply the transposed by the original eigenvalues.
  matrix_multiply(eigenvectors, transposed_eigenvectors, matrix);
  free_matrix(eigenvectors);
  free_matrix(transposed_eigenvectors);
}

/****************************************************************************
 * Compute a diffusion on the given matrix.
 *
 * This function requires a square matrix.
 ****************************************************************************/
MATRIX_T* diffuse_matrix
(MTYPE     diffusion_constant,
 int       iterations,
 MATRIX_T* matrix)
{
  ARRAY_T* current_activations;
  ARRAY_T* next_activations;
  MATRIX_T* diffused_matrix;

  // Make sure the matrix is square.
  int matrix_size = get_num_rows(matrix);
  if (matrix_size != get_num_cols(matrix)) {
    die("Attempted to convert a non-square matrix (%d != %d).", 
	matrix_size, get_num_cols(matrix));
  }

  // Normalize each row to sum to 1.
  normalize_rows(0, matrix);
 
  // Allocate space.
  current_activations = allocate_array(matrix_size);
  next_activations = allocate_array(matrix_size);
  diffused_matrix = allocate_matrix(matrix_size, matrix_size);

  // Initialize the diffusion matrix.
  init_matrix(0, diffused_matrix);

  // Perform the diffusion from each node.
  int start_node;
  for (start_node = 0; start_node < matrix_size; start_node++) {

    // Initialize the current activations.
    init_array(0, current_activations);
    set_array_item(start_node, 1, current_activations);

    // Run the diffusion for the specified number of iterations.
    int iteration;
    for (iteration = 0; iteration < iterations; iteration++) {

      // Initialize the next activations to 0.
      init_array(0, next_activations);

      // Perform the diffusion from each node.
      int i_row;
      for (i_row = 0; i_row < matrix_size; i_row++) {
	MTYPE current_activation = get_array_item(i_row, current_activations);

	// Don't bother diffusing if there is no activation here.
	if (current_activation != 0.0) {
	  int i_col;
	  for (i_col = 0; i_col < matrix_size; i_col++) {
	    MTYPE edge_weight = get_matrix_cell(i_row, i_col, matrix);
	    MTYPE new_activation 
	      = current_activation * edge_weight * diffusion_constant;
	    incr_array_item(i_col, new_activation, next_activations);
	  }
	}
      }

      // Copy over the new activations.
      copy_array(next_activations, current_activations);
    }

    // Put the resulting activations into the diffusion matrix.
    set_matrix_row(start_node, current_activations, diffused_matrix);
  }

  // Free local space.
  free_array(next_activations);
  free_array(current_activations);

  // Return the diffused matrix.
  return(diffused_matrix);
}

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
