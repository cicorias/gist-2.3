/*****************************************************************************
 * FILE: kernel-pca.c
 * AUTHOR: William Stafford Noble
 * CREATE DATE: 12/14/99
 * PROJECT: SVM
 * COPYRIGHT: 1999-2001, Columbia University (see ../doc/copyright.html)
 * VERSION: $Revision: 1.1 $
 * DESCRIPTION: Given a set of training examples in RDB format,
 * compute kernel-based eigenvectors.  This program follows the
 * outline given in "Neural Networks" by Haykin (1999), p. 435-436.
 *****************************************************************************/
#include "kernel-pca.h"
#include "compute-kernel.h"
#include "linear-algebra.h"
#include "matrix.h"
#include "array.h"
#include "utils.h"
#include <time.h>
#include <string.h>
#include <stdlib.h>

/*****************************************************************************
 * Normalize a given set of eigenvectors so that the dot product of
 * the eigenvector with itself equals the reciprocal of the
 * corresponding eigenvalue.
 *
 * This function assumes that the eigenvectors with zero eigenvalues
 * have already been removed.
 *****************************************************************************/
static void normalize_eigenvectors
  (ARRAY_T*  eigenvalues,
   MATRIX_T* eigenvectors)
{
  int num_eigens;
  int i_eigen;
  double   eigenvalue;
  ARRAY_T* eigenvector;
  double   eigen_product;
  double   scale_factor;

  num_eigens = get_num_cols(eigenvectors);
  for (i_eigen = 0; i_eigen < num_eigens; i_eigen++) {

    /* Get this eigenvector and self-multiply. */
    eigenvector = get_matrix_column(i_eigen, eigenvectors);
    eigen_product = dot_product(eigenvector, eigenvector);

    /* Get the corresponding eigenvalue. */
    eigenvalue = get_array_item(i_eigen, eigenvalues);

    /* Scale the eigenvector appropriately. */
    scale_factor = sqrt(eigenvalue / eigen_product), 
    scalar_mult(scale_factor, eigenvector);

#ifdef VERIFY_EIGENS
    if (!almost_equal(dot_product(eigenvector, eigenvector), eigenvalue,
		      0.000001)) {
      fprintf(stderr, "eigenvalue=%g eigen_product=%g\n", eigenvalue, 
	      eigen_product);
      die("Normalization failed (%g != %g).\n", 
	  dot_product(eigenvector, eigenvector), eigenvalue);
    } else {
      fprintf(stderr, "Verified eigenvalue=%g eigen_product=%g\n", eigenvalue, eigen_product);
    }
#endif

    /* Put the scaled eigenvector back in. */
    set_matrix_column(eigenvector, i_eigen, eigenvectors);
    free_array(eigenvector);
  }
}

/************************************************************************
 * Reduce the size of the eigenvector matrix by selecting only the
 * largest ones.
 ************************************************************************/
void select_eigenvectors
  (int        num_eigens,
   double     eigen_threshold,
   ARRAY_T**   eigenvalues,
   MATRIX_T** eigenvectors)
{
  MATRIX_T* new_eigenvectors;
  ARRAY_T*  new_eigenvalues;
  ARRAY_T*  eigenvector;
  int       i_eigen;
  double    eigen_sum;
  double    percent_variance;
  
  /* If zero were selected, then select all. */
  if (num_eigens == 0) {
    num_eigens = get_num_cols(*eigenvectors);
  }

  /* Store the sum of the eigenvalues. */
  eigen_sum = array_total(*eigenvalues);

  /* Make sure all the eigenvalues are above the given threshold. */
  for (i_eigen = 0; i_eigen < num_eigens; i_eigen++) {
    percent_variance = get_array_item(i_eigen, *eigenvalues) / eigen_sum;
    
    if (percent_variance < eigen_threshold) {
      if (verbosity >= NORMAL_VERBOSE) {
	fprintf(stderr, "Retaining %d eigenvalues.\n", i_eigen);
      }
      break;
    }
  }
  num_eigens = i_eigen;

  /* Allocate the new, smaller set of eigenvectors and eigenvalues. */
  new_eigenvalues = allocate_array(num_eigens);
  new_eigenvectors = allocate_matrix(get_num_rows(*eigenvectors),
				     num_eigens);

  /* Copy the eigenvectors and the eigenvalues. */
  for (i_eigen = 0; i_eigen < num_eigens; i_eigen++) {
    set_array_item(i_eigen, get_array_item(i_eigen, *eigenvalues), 
		   new_eigenvalues);
    eigenvector = get_matrix_column(i_eigen, *eigenvectors);
    set_matrix_column(eigenvector, i_eigen, new_eigenvectors);
    free_array(eigenvector);
  }

  /* Replace the old with the new. */
  free_array(*eigenvalues);
  *eigenvalues = new_eigenvalues;
  free_matrix(*eigenvectors);
  *eigenvectors = new_eigenvectors;
}


/*****************************************************************************
 * Find eigenvectors of a given matrix, normalized so that the dot
 * product of the eigenvector with itself equals the reciprocal of the
 * corresponding eigenvalue.  In the output, the eigenvectors are
 * sorted by increasing magnitude.
 *****************************************************************************/
MATRIX_T* find_normalized_eigenvectors
  (int       num_eigens,          /* Maximum number of eigenvectors to find. */
   double    eigen_threshold,     /* Only print eigenvectors whose 
				     corresponding eigenvalues account for at
				     least this much of the variance. */
   char*     eigenvalue_filename, /* Print the eigenvalues to this file. */
   MATRIX_T* matrix)              /* Matrix to analyze. */
{
  ARRAY_T*  eigenvalues;  /* The raw eigenvalues. */
  MATRIX_T* eigenvectors; /* The raw eigenvectors. */
  FILE*     eigenvalue_file;

  /* Compute the eigenvectors. */
  find_symmetric_eigenvectors(matrix, &eigenvalues, &eigenvectors);

  /* Select the biggest ones. */
  select_eigenvectors(num_eigens, eigen_threshold, &eigenvalues,
		      &eigenvectors);

  /* Print the eigenvalues. */
  if (verbosity > NORMAL_VERBOSE) {
    fprintf(stderr, "Eigenvalues: ");
    print_array(eigenvalues, 6, 3, TRUE, stderr);
  }

  /* Normalize the eigenvectors. */
  normalize_eigenvectors(eigenvalues, eigenvectors);

  /* Print the eigenvalues, if requested. */
  if (eigenvalue_filename != NULL) {
    if (verbosity >= NORMAL_VERBOSE) {
      fprintf(stderr, "Storing eigenvalues in %s.\n", eigenvalue_filename);
    }
    if (open_file(eigenvalue_filename, "w", TRUE, "eigenvalue", "eigenvalues",
		  &eigenvalue_file) == 0) {
      exit(1);
    }
    print_array(eigenvalues, 10, 5, TRUE, eigenvalue_file);
    fclose(eigenvalue_file);
  }

  /* Return the eigenvector matrix. */
  free_array(eigenvalues);
  return(eigenvectors);
}

/*****************************************************************************
 * A goofy helper function that assigns names to the columns in the
 * eigenvector matrix.
 *****************************************************************************/
void add_eigen_col_names
  (RDB_MATRIX_T* rdb_eigenvectors)
{
  int num_cols;
  int i_col;
  char col_name[100];

  num_cols = get_num_cols(get_raw_matrix(rdb_eigenvectors));
  for (i_col = 0; i_col < num_cols; i_col++) {
    sprintf(col_name, "eigen%d", i_col);
    add_col_name(col_name, rdb_eigenvectors);
  }
}

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
