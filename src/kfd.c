/*****************************************************************************
 * FILE: kfd.c
 * AUTHOR: William Stafford Noble
 * CREATE DATE: 12/14/99
 * PROJECT: SVM
 * COPYRIGHT: 1999-2001, Columbia University (see ../doc/copyright.html)
 * VERSION: $Revision: 1.2 $
 * DESCRIPTION: Find the kernel Fisher discriminant.
 * 
 * This file follows the description given in "Fisher Discriminant
 * Analysis with Kernels" by Mika, R\"{a}tsch, Weston, Sch\"{o}lkopf
 * and M\"{u}ller.
 *
 * The implementation is currently absurdly space-inefficient.
 *****************************************************************************/
#include "kfd.h"
#include "matrix.h"
#include "class-array.h"
#include "linear-algebra.h"
#include <assert.h>

/*****************************************************************************
 * Given a set of data and class i, compute the following array M_i of
 * length N:
 *
 *      (M_i)_j := (1/l_i) \sum_{k=1}^{l_i} k(x_j, x_k^i)
 * 
 * where N is the total number of training examples, l_i is the number
 * of elements in class i, x_j is the jth training example, x_k^i is
 * the kth element in class i, and k is the kernel matrix.
 *****************************************************************************/
static ARRAY_T* compute_kernel_means
  (BOOLEAN_T      class,
   CLASS_ARRAY_T* classes,
   MATRIX_T*      kernel_matrix)
{
  int num_rows;
  ARRAY_T* kernel_means;
  int i_row;
  int i_col;
  int i_class;
  double class_total;

  /* Allocate the return array. */
  num_rows = get_num_rows(kernel_matrix);
  kernel_means = allocate_array(num_rows);

  for (i_row = 0; i_row < num_rows; i_row++) {

    /* Add up the kernel matrix entries for this class. */
    i_class = 0;
    class_total = 0.0;
    for (i_col = 0; i_col < num_rows; i_col++) {
      if (get_class(i_col, classes) == class) {
	class_total += get_matrix_cell(i_row, i_col, kernel_matrix);
	i_class++;
      }
    }

    /* Divide by the number of items in the class, and store. */
    set_array_item(i_row, class_total / i_class, kernel_means);
  }

  return(kernel_means);
}


/*****************************************************************************
 * Multiply a given array by its transpose.  The result, from an array
 * of length N, is an NxN matrix.
 *****************************************************************************/
MATRIX_T* multiply_by_transpose
  (ARRAY_T* array)
{
  int num_elements;
  MATRIX_T* matrix;
  int i_row; 
  int i_col;
  double this_row;
  double this_col;

  /* Allocate the matrix. */
  num_elements = get_array_length(array);
  matrix = allocate_matrix(num_elements, num_elements);

  /* Fill in the matrix. */
  for (i_row = 0; i_row < num_elements; i_row++) {
    this_row = get_array_item(i_row, array);
    for (i_col = 0; i_col < num_elements; i_col++) {
      this_col = get_array_item(i_col, array);
      set_matrix_cell(i_row, i_col, this_row * this_col, matrix);
    }
  }
  return(matrix);
}

/*****************************************************************************
 * Extract from a kernel matrix a matrix containing only the members
 * of a given class.
 *****************************************************************************/
MATRIX_T* extract_class_matrix
  (BOOLEAN_T      class,
   CLASS_ARRAY_T* classes,
   MATRIX_T*      kernel_matrix)
{
  int num_data;
  int num_class;
  MATRIX_T* class_matrix;
  int i_col;
  int i_class;
  ARRAY_T*  this_column;

  /* Allocate the return matrix. */
  num_data = get_num_rows(kernel_matrix);
  num_class = get_class_size(class, classes);
  assert(num_class != 0);
  class_matrix = allocate_matrix(num_data, num_class);

  /* Copy columns from one matrix to the other. */
  i_class = 0;
  for (i_col = 0; i_col < num_data; i_col++) {
    if (get_class(i_col, classes) == class) {
      this_column = get_matrix_column(i_col, kernel_matrix);
      set_matrix_column(this_column, i_class, class_matrix);
      free_array(this_column);
      i_class++;
    }
  }
  return(class_matrix);
}
/*****************************************************************************
 * I don't have a good name for this function because I don't really
 * know why it's doing what it's doing.
 *****************************************************************************/
static MATRIX_T*  make_n_matrix
  (MATRIX_T* matrix)
{
  MATRIX_T* identity_matrix;
  MATRIX_T* intermediate_matrix;
  MATRIX_T* transposed_matrix;
  MATRIX_T* n_matrix;
  int num_rows;
  int num_cols;
  int i_col;
  double value;

  /* Get the dimensions of the matrix. */
  num_rows = get_num_rows(matrix);
  num_cols = get_num_cols(matrix);

  /* Create an identity matrix. */
  identity_matrix = allocate_matrix(num_cols, num_cols);
  init_matrix(0.0, identity_matrix);
  for (i_col = 0; i_col < num_cols; i_col++) {
    set_matrix_cell(i_col, i_col, 1.0, identity_matrix);
  }

  /* Subtract off the requisite value. */
  value = -1.0 / num_cols;
  scalar_add_matrix(value, identity_matrix);

  /* Multiply them. */
  intermediate_matrix = matrix_multiply(matrix, identity_matrix, NULL);
  free_matrix(identity_matrix);

  /* Transpose the matrix and multiply again. */
  transposed_matrix = transpose_matrix(matrix);
  n_matrix = matrix_multiply(intermediate_matrix, transposed_matrix, NULL);
  free_matrix(intermediate_matrix);
  free_matrix(transposed_matrix);
  
  return(n_matrix);
}


/*****************************************************************************
 * Find the kernel Fisher discriminant.
 *****************************************************************************/
MATRIX_T* find_kfd
  (double         regularizer,
   CLASS_ARRAY_T* classes,
   MATRIX_T*      kernel_matrix)
{
  ARRAY_T*  positive_kernel_means;
  ARRAY_T*  negative_kernel_means;
  ARRAY_T*  diff_kernel_means;
  MATRIX_T* mean_matrix;
  MATRIX_T* positive_kernel_matrix;
  MATRIX_T* negative_kernel_matrix;
  MATRIX_T* positive_n_matrix;
  MATRIX_T* negative_n_matrix;
  MATRIX_T* n_matrix;
  MATRIX_T* n_invert;
  MATRIX_T* target_matrix;
  ARRAY_T*  kfd_array;
  MATRIX_T* kfd;

  /* Compute the mean of the kernel values in each class. */
  positive_kernel_means = compute_kernel_means(TRUE, classes, kernel_matrix);
  negative_kernel_means = compute_kernel_means(FALSE, classes, kernel_matrix);
  
  /* Subtract the two sets of means. */
  scalar_mult(-1.0, negative_kernel_means);
  sum_array(positive_kernel_means, negative_kernel_means);
  diff_kernel_means = negative_kernel_means;
  free_array(positive_kernel_means);

  /* Multiply this matrix by its transpose. */
  mean_matrix = multiply_by_transpose(diff_kernel_means);
  free_array(diff_kernel_means);

  /* Compute the class kernel matrices. */
  positive_kernel_matrix = extract_class_matrix(TRUE, classes, kernel_matrix);
  negative_kernel_matrix = extract_class_matrix(FALSE, classes, kernel_matrix);

  /* Combine the class matrices. */
  positive_n_matrix = make_n_matrix(positive_kernel_matrix);
  negative_n_matrix = make_n_matrix(negative_kernel_matrix);
  free_matrix(positive_kernel_matrix);
  free_matrix(negative_kernel_matrix);

  /* Combine the two matrices. */
  sum_matrices(positive_n_matrix, negative_n_matrix);
  free_matrix(positive_n_matrix);
  n_matrix = negative_n_matrix;

  /* Regularize. */
  add_to_diagonal(regularizer, n_matrix);

  /* Invert the N matrix. */
  n_invert = invert_matrix(n_matrix);
  if (n_invert == NULL) {
    die("Increase regularization (%g).", regularizer);
  }
  free_matrix(n_matrix);

  /* Multiply the inverted N matrix by the mean matrix. */
  target_matrix = matrix_multiply(n_invert, mean_matrix, NULL);
  free_matrix(n_invert);
  free_matrix(mean_matrix);

  /* Find the first eigenvector of the result. */
  kfd_array = find_first_eigenvector(target_matrix);
  free_matrix(target_matrix);
  
  /* Convert it to a matrix. */
  kfd = allocate_matrix(get_array_length(kfd_array), 1);
  set_matrix_column(kfd_array, 0, kfd);
  free_array(kfd_array);
  
  return(kfd);
}


/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
