/*****************************************************************************
 * FILE: project.c
 * AUTHOR: William Stafford Noble
 * CREATE DATE: 1/26/2000
 * PROJECT: SVM
 * COPYRIGHT: 1999-2001, Columbia University (see ../doc/copyright.html)
 * VERSION: $Revision: 1.1 $
 * DESCRIPTION: Project a given set of data onto a set of eigenvectors.
 *****************************************************************************/
#include "project.h"
#include "matrix.h"
#include "array.h"

/*****************************************************************************
 * Compute the projections of a given set of data onto a given set of
 * eigenvectors, using a given kernel matrix.
 *
 * "For the extraction of principal components of a test point x,
 * compute the projections
 *
 *    a_k = \sum_{j=1}^N \alpha_{k,j} K(x_j, x),    k = 1, 2, ..., p
 *
 * where \alpha_{k,j} is the jth element of eigenvector \alpha_k"
 * (Haykin 1999 p. 437).
 *****************************************************************************/
void project_data
  (MATRIX_T* kernel_matrix,
   MATRIX_T* eigenvectors,
   MATRIX_T* output_matrix)
{
  int      i_output_row;
  int      num_output_rows;
  int      i_output_col;
  int      num_output_cols;
  ARRAY_T* eigenvector;
  ARRAY_T* kernel_row;
  double   projection;

  /* Get the output matrix dimensions. */
  num_output_rows = get_num_rows(output_matrix);
  num_output_cols = get_num_cols(output_matrix);

  /* Fill in the output matrix row by column. */
  for (i_output_row = 0; i_output_row < num_output_rows; i_output_row++) {

    if (verbosity >= DUMP_VERBOSE) {
      fprintf(stderr, "Projecting row %d.\n", i_output_row);
    }

    for (i_output_col = 0; i_output_col < num_output_cols; i_output_col++) {

      /* Get the corresponding eigenvector and kernel row. */
      eigenvector = get_matrix_column(i_output_col, eigenvectors);
      kernel_row = get_matrix_row(i_output_row, kernel_matrix);

      /* Compute the projection. */
      projection = dot_product(eigenvector, kernel_row);
      set_matrix_cell(i_output_row, i_output_col, projection, output_matrix);

      free_array(eigenvector);
    }
  }
}


/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
