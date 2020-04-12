/**************************************************************************
 * FILE: matrix.h
 * AUTHOR: William Stafford Noble
 * CREATE DATE: 2-4-97
 * PROJECT: shared
 * COPYRIGHT: 1999-2001, Columbia University
 * VERSION: $Revision: 1.4 $
 * DESCRIPTION: Some simple matrix manipulation routines.
 **************************************************************************/
#ifndef MATRIX_H
#define MATRIX_H

#ifdef ARRAY_H
#warning "array.h included before matrix.h"
#warning "Matrix type may be ignored"
#endif

#ifdef IMATRIX
#define IARRAY
#else
#ifdef SMATRIX
#define SARRAY
#else
#ifdef LMATRIX
#define LARRAY
#else
#define NOT_INT
#endif
#endif
#endif

#include "array.h"
#include <stdio.h>

/**************************************************************************
 * Uses floats by default.  Defining IMATRIX, SMATRIX, LMATRIX or
 * DMATRIX changes the type.
 **************************************************************************/

#define MTYPE ATYPE
#define MSCAN ASCAN

/***************************************************************************
 * Define a matrix type.
 ***************************************************************************/
typedef struct matrix_t {
  int       num_rows;
  int       num_cols;
  ARRAY_T** rows;
}MATRIX_T;

/**************************************************************************
 * Allocate a matrix.
 **************************************************************************/
MATRIX_T* allocate_matrix
  (int num_rows,
   int num_columns);

/**************************************************************************
 * Allocate a matrix and fill it with ones on the diagonal and zeroes
 * off the diagonal.
 **************************************************************************/
MATRIX_T* make_diagonal_matrix
  (int matrix_size);

/**************************************************************************
 * Grow a matrix by adding one row to it.  Signal an error if the row
 * does not have the same number of columns as the given matrix.
 **************************************************************************/
void grow_matrix
  (ARRAY_T*   one_row,
   MATRIX_T*  matrix);

/**************************************************************************
 * Basic access routines.
 **************************************************************************/
int get_num_rows
  (MATRIX_T* matrix);

int get_num_cols
  (MATRIX_T* matrix);

ARRAY_T* get_matrix_row
 (int       row,
  MATRIX_T* matrix);

void set_matrix_row
 (int       row,
  ARRAY_T* one_row,
  MATRIX_T* matrix);

#ifdef BOUNDS_CHECK
#define get_matrix_cell(row,col,matrix) \
   get_matrix_cell_defcheck(row,col,matrix)
#define set_matrix_cell(row,col,value,matrix) \
   set_matrix_cell_defcheck(row,col,value,matrix)
#define incr_matrix_cell(row,col,value,matrix) \
   incr_matrix_cell_defcheck(row,col,value,matrix)
#else
#define get_matrix_cell(row,col,matrix) \
   get_array_item(col, ((MATRIX_T*)matrix)->rows[row])
#define set_matrix_cell(row,col,value,matrix) \
   set_array_item(col, value, ((MATRIX_T *)matrix)->rows[row])
#define incr_matrix_cell(row,col,value,matrix) \
   set_array_item(col, get_array_item(col, ((MATRIX_T*)matrix)->rows[row]) + \
                                      value,((MATRIX_T*)matrix)->rows[row])
#endif

MTYPE get_matrix_cell_defcheck
  (int       row,
   int       col,
   MATRIX_T* matrix);

void set_matrix_cell_defcheck
  (int       row,
   int       col,
   MTYPE     value,
   MATRIX_T* matrix);

void incr_matrix_cell_defcheck
  (int       row,
   int       col,
   MTYPE     value,
   MATRIX_T* matrix);

/***********************************************************************
 * Get a column from a matrix.
 *
 * Returns a newly allocated copy of the requested column.
 ***********************************************************************/
ARRAY_T* get_matrix_column
  (int       i_col,
   MATRIX_T* matrix);

/***********************************************************************
 * Set a column in a matrix. 
 ***********************************************************************/
void set_matrix_column
  (ARRAY_T*  column,
   int       i_col,
   MATRIX_T* matrix);

/***********************************************************************
 * Turn an array into a matrix.
 ***********************************************************************/
MATRIX_T* array_to_matrix
  (BOOLEAN_T one_row, /* Put the array in one row, or in many. */
   ARRAY_T*  array);

/**************************************************************************
 * Copy a matrix.
 **************************************************************************/
void copy_matrix
  (MATRIX_T* source_matrix,
   MATRIX_T* target_matrix);

/**************************************************************************
 * Initialize all cells of a given matrix to a given value.
 **************************************************************************/
void init_matrix
  (MTYPE      value,
   MATRIX_T*  matrix);

/**************************************************************************
 * Convert all entries in a given matrix to -1/1s, using a given
 * threshold.
 **************************************************************************/
void binarize_matrix
  (MTYPE      value,
   MATRIX_T*  matrix);

/**************************************************************************
 * Fill a matrix with a given raw matrix of values.
 **************************************************************************/
void fill_matrix
  (MTYPE*     raw_matrix,
   MATRIX_T*  matrix);

/**************************************************************************
 * Compute the sum of two matrices, assuming they have the same
 * dimension.  The sum is stored in the second matrix.
 **************************************************************************/
void sum_matrices
  (MATRIX_T* matrix1,
   MATRIX_T* matrix2);

/**************************************************************************
 * Extract the diagonal from a square matrix and return it in a
 * newly-allocated array.
 **************************************************************************/
ARRAY_T* extract_diagonal
  (MATRIX_T* matrix);

/**************************************************************************
 * Determine whether a given matrix is symmetric.
 **************************************************************************/
BOOLEAN_T is_symmetric
  (BOOLEAN_T verbose,
   MTYPE     slop,
   MATRIX_T* matrix);

/**************************************************************************
 * Create a matrix M such that M[x,y] = (M1[x,y] + M2[y,x]) / 2.
 **************************************************************************/
MATRIX_T* average_across_diagonal
  (MATRIX_T* matrix1,
   MATRIX_T* matrix2);

/**************************************************************************
 * Compute the mean off-diagonal absolute value in a square matrix.
 **************************************************************************/
double mean_off_diagonal_absolute_value
  (MATRIX_T* matrix);

/**************************************************************************
 * Add a given value to each element in the diagonal of a square matrix.
 **************************************************************************/
void add_to_diagonal
  (MTYPE     value,
   MATRIX_T* matrix);

/***********************************************************************
 * Determine whether two matrices are equal, within a given bound.
 ***********************************************************************/
BOOLEAN_T equal_matrices
  (ATYPE    close_enough,
   MATRIX_T* matrix1,
   MATRIX_T* matrix2);

/**************************************************************************
 * Multiply all items in a matrix by a given scalar.
 **************************************************************************/
void scalar_mult_matrix
  (MTYPE      value,
   MATRIX_T*  matrix);

/**************************************************************************
 * Add a scalar to all items in a matrix.
 **************************************************************************/
void scalar_add_matrix
  (MTYPE      value,
   MATRIX_T*  matrix);

/**************************************************************************
 * Raise all entries in a given matrix to a given power.
 **************************************************************************/
void power_matrix
  (MTYPE      value,
   MATRIX_T*  matrix);

/**************************************************************************
 * Multiply together corresponding values in two matrices of equal
 * dimension.  Store the result in the second matrix.
 **************************************************************************/
void mult_matrix
  (MATRIX_T* matrix1,
   MATRIX_T* matrix2);

/**************************************************************************
 * Mix two matrices in log space.
 **************************************************************************/
void mix_log_matrices
  (float     mixing,      /* Percent of matrix2 to be retained. */
   MATRIX_T* matrix1,
   MATRIX_T* matrix2);

/**************************************************************************
 * Read a matrix from a file.
 *
 * Each row of the matrix must appear on its own line.
 **************************************************************************/
MATRIX_T* read_matrix
  (FILE * infile);

/**************************************************************************
 * Read a matrix from a file, when the dimensions are known in advance.
 **************************************************************************/
MATRIX_T* read_known_matrix
  (int    num_rows,
   int    num_cols,
   FILE * infile);

/**************************************************************************
 * Print the matrix, optionally with row and column indices.
 **************************************************************************/
void print_matrix
  (MATRIX_T* matrix,        /* The matrix to be printed. */
   int       width,         /* Width of each cell. */
   int       precision,     /* Precision of each cell. */
   BOOLEAN_T print_titles,  /* Include row and column indices? */
   FILE*     outfile);      /* File to which to write. */

/**************************************************************************
 * void free_matrix
 **************************************************************************/
void free_matrix
  (MATRIX_T* matrix);


/***********************************************************************
 * Fill a matrix with random values between 0 and a given number.
 *
 * Assumes that the random number generator is initialized.
 ***********************************************************************/
void randomize_matrix
  (MTYPE     max_value,
   MATRIX_T* matrix);

/***********************************************************************
 * Compute the sum of the elements in a matrix.
 ***********************************************************************/
MTYPE sum_of_matrix
  (MATRIX_T* matrix);

/***********************************************************************
 * Compute the sum of the squares of a matrix.
 ***********************************************************************/
MTYPE sum_of_squares_matrix
  (MATRIX_T* matrix);

/***********************************************************************
 * Compute the sum of the squares of the differences between two
 * matrices.
 ***********************************************************************/
MTYPE sum_of_square_diff_matrices
  (MATRIX_T* matrix1,
   MATRIX_T* matrix2);

/***********************************************************************
 * Subtract the mean from each row or column of a matrix.
 ***********************************************************************/
void zero_mean_matrix_rows
  (MATRIX_T* matrix);

void zero_mean_matrix_cols
  (MATRIX_T* matrix);

/***********************************************************************
 * Linearly rescale each row or column of a matrix.
 ***********************************************************************/
void rescale_matrix_rows
  (MTYPE     minimum,
   MTYPE     maximum,
   MATRIX_T* matrix);

void rescale_matrix_cols
  (MTYPE     minimum,
   MTYPE     maximum,
   MATRIX_T* matrix);

/***********************************************************************
 * Divide each matrix row by its standard deviation.
 ***********************************************************************/
void variance_one_matrix_rows
  (MATRIX_T* matrix);

/***********************************************************************
 * Rank transform each matrix column (Paul Pavlidis)
 ***********************************************************************/
void rank_transform_matrix_cols
(MATRIX_T* matrix);

/***********************************************************************
 * Normalize the rows of a matrix.
 ***********************************************************************/
void normalize_rows
  (float tolerance,
   MATRIX_T*   matrix);

/***********************************************************************
 * void normalize_matrix
 *
 * Rangarajan et al. (Neural Comp. (8) 1041-1060) mention a proof by
 * Sinkhorn (1964) that a square matrix can be converted to a double
 * stochastic matrix, in which each row and each column sums to 1.0,
 * via an iterative procedure in which rows and columns are normalized
 * alternately.  This program implements that procedure.
 ***********************************************************************/
void normalize_matrix
  (float tolerance,
   MATRIX_T*   matrix);

/*****************************************************************************
 * Extract one margin of a matrix and return it as an array.
 *****************************************************************************/
ARRAY_T* get_matrix_row_sums
  (MATRIX_T* matrix);
ARRAY_T* get_matrix_col_sums
  (MATRIX_T* matrix);

/*****************************************************************************
 * Sort a given matrix by row, according to a given set of sort keys.
 *****************************************************************************/
void sort_matrix_rows
  (BOOLEAN_T reverse_sort,
   ARRAY_T*  keys,
   MATRIX_T* matrix);

/*****************************************************************************
 * Randomly shuffle the order of rows in a given matrix
 *****************************************************************************/
void shuffle_matrix_row_order
  (MATRIX_T* matrix);

/*****************************************************************************
 * Randomly shuffle the entries of a given matrix.
 *****************************************************************************/
void shuffle_matrix
  (MATRIX_T* matrix);

/*****************************************************************************
 * Given a label vector (in single-column matrix form) consisting of
 * 1s and -1s, compute a square matrix in which each entry is the
 * product of the corresponding labels.
 *****************************************************************************/
MATRIX_T* create_label_matrix
  (MATRIX_T* labels);

/*****************************************************************************
 * Compute the alignment between two Gram matrices, a la Cristianini
 * et al. (NIPS 2002).  The formula is
 *
 *                 <K1, K2>
 *         ------------------------
 *         \sqrt{<K1, K1> <K2, K2>}
 *
 * where <,> is the inner product between matrices.
 *****************************************************************************/
MTYPE alignment
  (MATRIX_T* matrix1,
   MATRIX_T* matrix2);

/*****************************************************************************
 * Replace the diagonal of a matrix with the row sum minus the
 * diagonal value.
 *
 * This function is equivalent to step 2 of Ng's spectral clustering
 * algorithm.
 *****************************************************************************/
void replace_diagonal
  (MATRIX_T* matrix);

/*****************************************************************************
 * Compute a square matrix of Euclidean distances from a square kernel matrix.
 *
 *     d(x,y) = \sqrt( K(x,x) - 2 K(x,y) + K(y,y) )
 *
 * Here is the derivation:
 *
 * ||x-y|| = \sqrt((x_1 - y_1)^2 + (x_2 - y_2)^2 + ... + (x_n - y_n)^2)
 *         = \sqrt((x_1^2 - 2x_1y_1 + y_1^2) + (x_2^2 - 2x_2y_2 + y_2^2)
 *                  + ... + (x_n^2 - 2x_ny_n + y_n^2))
 *         = \sqrt((x_1^2 + x_2^2 + ... + x_n^2) - 2(x_1y_1 + x_2y_2
 *                  + ... + x_ny_n) + (y_1^2 + y_2^2 + ... + y_n^2))
 *         = \sqrt(K(x,x) - 2 K(x,y) + K(y,y))
 *
 *****************************************************************************/
MATRIX_T* euclidean
  (MATRIX_T* matrix);

/****************************************************************************
 * Center a possibly non-square kernel matrix with respect to a
 * second, square matrix.
 * 
 * After centering, each feature must have zero mean.  We can't do the
 * centering in the feature space because we don't have access to the
 * feature vectors.  We can, however, express the centered kernel
 * matrix in terms of the non-centered version.
 *
 * Let \phi translate from the input space to the kernel space.  Then
 * we want to find \~K, the centered kernel matrix.  \~K can be
 * expressed in terms of the original matrix K (of dimensionality M)
 * as follows:
 *
 * \~K_{ij} = [\phi(x_i) - 1/M \sum_{m=1}^M \phi(x_m)]
 *             . [\phi(x_j) - 1/M \sum_{n=1}^M \phi(x_n)]
 *
 *          = [\phi(x_i) \phi(x_j)] - [1/M \phi(x_i) \sum_{m=1}^M \phi(x_m)]
 *             - [1/M \phi(x_j) \sum_{n=1}^M \phi(x_n)]
 *             + [1/M^2 \sum_{m=1}^M \sum{n=1}^M \phi(x_m) \phi(x_n)]
 *
 *          = [K_{ij}] - [1/M \sum_{m=1}^M K_{mj}] - [1/M \sum_{n=1}^M K_{in}]
 *             + [1/M^2 \sum_{m=1}^M \sum{n=1}^M K_{mn}]
 *
 *
 * 
 * This last expression can be generalized if we assume that we have
 * an N x N train matrix and an M x N test matrix.  We want to
 * subtract from each example in the test matrix the feature means
 * computed from the train matrix.  Assume that example i comes from
 * the test matrix.  Then the second and third terms above corresond
 * to the row average for the ith row in the test matrix and the jth
 * row in the train matrix.  The fourth term is still the mean in the
 * train matrix.
 ****************************************************************************/
void center_matrix
  (MATRIX_T* train_matrix,
   MATRIX_T* test_matrix);

/*****************************************************************************
 * Compute the Laplacian of a square matrix.  This is done by
 * subtracting the row margins from the diagonal.
 *****************************************************************************/
void laplacian
(MATRIX_T* matrix);

#endif

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
