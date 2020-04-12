/****************************************************************************
 * FILE: fast-classify.c
 * AUTHOR: William Stafford Noble
 * CREATE DATE: 1/25/2000
 * PROJECT: SVM
 * COPYRIGHT: 1999-2001, Columbia University (see ../doc/copyright.html)
 * VERSION: $Revision: 1.3 $
 * DESCRIPTION: Main procedure for classifying or projecting.
 ****************************************************************************/
#include "test-main.h"
#include "classify.h"
#include "rdb-matrix.h"
#include "string-list.h"
#include "matrix.h"
#include "array.h"
#include "utils.h"
#include <stdio.h>
#include <string.h>
#include <unistd.h>

#ifndef VERBOSITY
#define VERBOSITY
VERBOSE_T verbosity;
#endif

/****************************************************************************
 * Main procedure
 ****************************************************************************/
#define MAX_ID 10000 // Maximum ID length.
#include "cmdline.h"
#include <time.h>

int main(int argc, char *argv[])
{
  // Command line options.
  char*     hyperplane_filename = NULL; // Hyperplane coordinates.
  FILE*     hyperplane_file = NULL;
  char*     test_filename = NULL;       // Test data to be classified.
  FILE*     test_file = NULL;
  BOOLEAN_T format_line = FALSE;     // Use format line in RDB files?
  BOOLEAN_T print_time = TRUE;       // Print timing data to output header.
  int       output_precision = 0;    // Number of digits in output.

  // Data structures.
  KERNEL_T*      kernel = NULL;      // Just used for storing kernel params.
  double         bias = 0.0;         // Bias term for discriminant.
  double         sum_of_weights = 0.0; // Sum of learned weights.
  ARRAY_T*       hyperplane = NULL;  // Hyperplane coordinates.
  char*          this_line = NULL;   // Storage for one input line.
  char*          row_name = NULL;    // Storage for row ID.
  ARRAY_T*       example = NULL;     // The current example.
  double         discriminant;       // Discriminant for this example.

  verbosity = NORMAL_VERBOSE;
  this_line = (char*)mymalloc(sizeof(char)*MAX_ROW);
  row_name = (char*)mymalloc(sizeof(char)*MAX_ROW);
  kernel = allocate_kernel(FALSE);

  // Get the command line arguments.
  DO_STANDARD_COMMAND_LINE
    (3,
     DATA_OPTN(1, hyperplane, <filename> (required),
	       hyperplane_filename = _OPTION_);
     DATA_OPTN(1, test, <filename> (required),
	       test_filename = _OPTION_);
     FLAG_OPTN(1, rdb, format_line = TRUE);
     FLAG_OPTN(1, notime, print_time = FALSE);
     DATA_OPTN(1,
	       precision, <value> (default = 6),
	       output_precision = atoi(_OPTION_));
     DATA_OPTN(1, verbose, 1|2|3|4|5 (default=2), 
	       verbosity = (VERBOSE_T)atoi(_OPTION_));
     );

  // Unless user-specified, set the output precision.
  if (output_precision == 0) {
    output_precision = DEFAULT_PRECISION;
  }

  // Make sure we got the required files.
  if (hyperplane_filename == NULL) {
    die("No hyperplane file given.");
  }
  if (test_filename == NULL) {
    die("No test set given.\n");
  }

  // Read the kernel parameters from the header of the learned file.
  // N.B. If we open the file once, and pass the file handle, then for
  // some reason 'project' no longer works.  -- WSN 4/4/03
  read_kernel_parameters(hyperplane_filename, &bias, &sum_of_weights, kernel);

  // Read the hyperplane coordinates.
  if (open_file(hyperplane_filename, "r", TRUE, "hyperplane", 
		"hyperplane coordinates", &hyperplane_file) == 0) {
    exit(1);
  }
  RDB_MATRIX_T* rdb = read_rdb_matrix(format_line, "NaN", hyperplane_file);
  hyperplane = get_matrix_column(0, get_raw_matrix(rdb));
  free_rdb_matrix(rdb);

  // Allocate the example.
  example = allocate_array(get_array_length(hyperplane));

  // Print the header.
  print_output_header(hyperplane_filename,
		      kernel,
		      print_time,
		      stdout);
  printf("gene\tclassification\tdiscriminant\n");

  // Open the test data file.
  if (open_file(test_filename, "r", TRUE, "test", "test data", 
		&test_file) == 0) {
    exit(1);
  }

  // Read the first row.
  read_one_row(test_file, MAX_ROW, this_line);

  // Keep reading till we get past the comments.
  while (is_comment(this_line)) {
    read_one_row(test_file, MAX_ROW, this_line);
  }

  // Read the next line, stopping if it's empty.
  int i_row = 1;
  while (fgets(this_line, MAX_ROW, test_file)) {

    // Parse the row.
    parse_rdb_row(this_line,
		  NULL,
		  NULL,
		  i_row, 
		  get_array_length(hyperplane),
		  row_name, 
		  example);

    // Classify this example.
    discriminant = fast_classify(example, hyperplane, sum_of_weights, kernel);
    
    // Print the prediction.
    if (discriminant >= 0.0) {
      printf("%s\t1.0\t%g\n", row_name, discriminant);
    } else {
      printf("%s\t-1.0\t%g\n", row_name, discriminant);
    }
  
    i_row++;
  }
  fclose(test_file);
  
  return(0);
}


/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
 
