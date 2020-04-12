/*****************************************************************************
 * FILE: gist-kernel.c
 * AUTHOR: William Stafford Noble
 * CREATE DATE: 6/24/04
 * PROJECT: SVM
 * VERSION: $Revision: 1.1 $
 * DESCRIPTION: Compute kernels
 *****************************************************************************/
#include "kernel.h"
#include "rdb-matrix.h"
#include "compute-kernel.h"
#include "linear-algebra.h"
#include "matrix.h"
#include "array.h"
#include "utils.h"
#include "cmdline.h"
#include <time.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>


#ifndef VERBOSITY
#define VERBOSITY
VERBOSE_T verbosity;
#endif


int main(int argc, char *argv[])
{

  // Required command line parameters.
  char*     data_filename = NULL;
  char*     kernel_filename = NULL;

  // Data structures.
  RDB_MATRIX_T* data_matrix = NULL;
  RDB_MATRIX_T* kernel_matrix = NULL;

  // Local variables.
  FILE*     data_file = NULL;
  BOOLEAN_T is_kernel = FALSE;            // Compute the scalar product.
  BOOLEAN_T format_line = FALSE;          // Include RDB format line.
  int       output_precision = 4;         // Number of digits in output.
  verbosity = NORMAL_VERBOSE;             // Set verbosity to stderr.

  // Parse the command line.
  DO_STANDARD_COMMAND_LINE
    (2,
     DATA_OPTN(1, data, <file>, data_filename = _OPTION_);
     DATA_OPTN(1, kernel, <file>, kernel_filename = _OPTION_);
     FLAG_OPTN(1, iskernel, is_kernel = TRUE);
     FLAG_OPTN(1, rdb, format_line = TRUE);
     DATA_OPTN(1, precision, <value> (default=4), 
	       output_precision = atoi(_OPTION_));
     DATA_OPTN(1, verbose, 1|2|3|4|5 (default=2),
	       verbosity = (VERBOSE_T)atoi(_OPTION_));
     );

  // Open both files.
  if (!open_file(data_filename, "r", TRUE, "data", "data", &data_file)) {
    exit(1);
  }

  // Read the data matrix.
  data_matrix = read_rdb_matrix(format_line, NULL, data_file);

  // Compute the base kernel, if required.
  if (is_kernel) {
    kernel_matrix = data_matrix;
  } else {
    kernel_matrix 
      = allocate_rdb_matrix(get_num_rows(get_raw_matrix(data_matrix)),
			    get_num_rows(get_raw_matrix(data_matrix)),
			    NULL);
    compute_base_kernel_matrix(get_raw_matrix(data_matrix),
			       get_raw_matrix(data_matrix),
			       get_raw_matrix(kernel_matrix));

    // Transfer the labels to the kernel matrix.
    set_corner_string(get_corner_string(data_matrix), kernel_matrix);    
    set_row_names(get_row_names(data_matrix), kernel_matrix);
    set_col_names(get_row_names(data_matrix), kernel_matrix);

  }

  // Do all the kernel operations.
  compute_kernel_transformations(kernel_filename,
				 get_raw_matrix(kernel_matrix));

  // Print the resulting matrix to stdout.
  print_rdb_matrix(kernel_matrix,
		   format_line,
		   output_precision + 2, 
		   output_precision,
		   stdout);

  return(0);
}


/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
