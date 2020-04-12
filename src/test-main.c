/****************************************************************************
 * FILE: test-main.c
 * AUTHOR: William Stafford Noble
 * CREATE DATE: 1/25/2000
 * PROJECT: SVM
 * COPYRIGHT: 1999-2001, Columbia University (see ../doc/copyright.html)
 * VERSION: $Revision: 1.5 $
 * DESCRIPTION: Main procedure for classifying or projecting.
 ****************************************************************************/
#include "test-main.h"
#include "kernel.h"
#include "classify.h"
#include "kernel-pca.h"
#include "project.h"
#include "compute-kernel.h"
#include "rdb-matrix.h"
#include "string-list.h"
#include "matrix.h"
#include "array.h"
#include "utils.h"
#include <stdio.h>
#include <string.h>
#include <unistd.h>

#ifndef CLASSIFY_DEBUG
#define CLASSIFY_DEBUG 0
#endif

#ifndef CLASSIFY
#define CLASSIFY 0
#endif

#ifndef PROJECT
#define PROJECT 0
#endif


#ifndef VERBOSITY
#define VERBOSITY
VERBOSE_T verbosity;
#endif


/****************************************************************************
 * Read the kernel name and polynomial parameters from the header of a
 * given weights file.
 ****************************************************************************/
#define MAX_LINE 500
void read_kernel_parameters
  (char*      learned_filename,
   double*    bias,
   double*    sum_of_weights,
   KERNEL_T*  kernel)
{
  char  one_line[MAX_LINE];
  FILE* learned_file;
  int   num_scanned;
  char  matrix_from_file_string[MAX_LINE];
  char  zero_mean_string[MAX_LINE];
  char  variance_one_string[MAX_LINE];
  char  normalize_string[MAX_LINE];
  char  radial_string[MAX_LINE];
  char  feature_select_string[MAX_LINE];

  /* Open the weight file for reading. */
  if (!open_file(learned_filename, "r", TRUE, "weights", "weights",
		 &learned_file)) {
    exit(1);
  }

  // Read until we find the first line of parameters.
  num_scanned = 0;
  while (num_scanned != 4) {

    if (fgets(one_line, MAX_LINE, learned_file) == NULL) {
      die("Error parsing weight file header.\n%s", one_line);
    }

    num_scanned = 
      sscanf(one_line, "# matrix_from_file=%s zero_mean=%s variance_one=%s normalize=%s",
	       matrix_from_file_string, zero_mean_string, variance_one_string,
	       normalize_string);

  }
  kernel->matrix_from_file = boolean_from_string(matrix_from_file_string);
  kernel->zero_mean_row = boolean_from_string(zero_mean_string);
  kernel->variance_one = boolean_from_string(variance_one_string);
  kernel->normalize = boolean_from_string(normalize_string);

  DEBUG_CODE(1, if (kernel->matrix_from_file) fprintf(stderr, "Reading matrix from file\n"););


  /* Parse the fifth line. */
  fgets(one_line, MAX_LINE, learned_file);
  num_scanned 
    = sscanf(one_line, 
	     "# constant=%lf coefficient=%lf power=%lf bias=%lf",
	     &(kernel->constant), &(kernel->coefficient), &(kernel->power),
	     bias);
  if (num_scanned < 3) {
    die("Error parsing weight file header.\n%s", one_line);
  } else if (num_scanned == 3) {
    // Two-class SVM is zero-bias by default.
    *bias = 0.0;
  }

  /* Parse the sixth line. */
  fgets(one_line, MAX_LINE, learned_file);
  num_scanned 
    = sscanf(one_line, 
	     "# radial=%s width_factor=%lf two_squared_width=%lf",
	     radial_string, &(kernel->width_factor), 
	     &(kernel->two_squared_width));
  if (num_scanned < 2) {
    die("Error parsing weight file header.\n%s", one_line);
  }
  kernel->radial = boolean_from_string(radial_string);

  // Don't add to the diagonal during classification.
  kernel->add_diag = 0.0;
  kernel->diagonal_factor = 0.0;

  // Read feature selection parameters from the seventh line.
  fgets(one_line, MAX_LINE, learned_file);
  num_scanned = sscanf(one_line, "# feature_select=%s", feature_select_string);
  if (num_scanned == 1) {
    if (strcmp(feature_select_string, "none") != 0) {
      die("Cannot perform feature selection during classification.");
    }
  }

  // If requested, read the sum of the weights from the eighth line.
  if (sum_of_weights != NULL) {
    fgets(one_line, MAX_LINE, learned_file);
    num_scanned = sscanf(one_line, "# sum_of_weights=%lf", sum_of_weights);
    if (num_scanned < 1) {
      die("Error parsing weight file header.\n%s", one_line);
    }
  }

  // Close the file.
  fclose(learned_file);
}

/****************************************************************************
 * Put an RDB header on the output file.
 ****************************************************************************/
void print_output_header
  (char*     learned_filename,
   KERNEL_T* kernel,
   BOOLEAN_T print_time,
   FILE*     outfile)
{
  if (CLASSIFY) {
    fprintf(outfile, "# Generated by classify\n");
  } else if (PROJECT) {
    fprintf(outfile, "# Generated by project\n");
  } else {
    fprintf(outfile, "# Generated by fast-classify\n");
  }

  // Print release, web site and citation info.
  fprintf(outfile, "# Gist, version %s\n", VERSION);
  fprintf(outfile, "# For more information, go to http://svm.sdsc.edu\n# \n");
  fprintf(outfile, "# If you use this software in your research, please cite:\n");
  fprintf(outfile, "#   Paul Pavlidis, Ilan Wapinski and William Stafford Noble.\n");
  fprintf(outfile, "#   \"Support vector machine classification on the web.\"\n");
  fprintf(outfile, "#   Bioinformatics. 20(4):586-587, 2004.\n# \n");

  /* Print file names. */
  if (kernel->train_filename != NULL) {
    fprintf(outfile, "# train_file=%s\n", kernel->train_filename);
  }
  if (learned_filename != NULL) {
    fprintf(outfile, "# learned_file=%s\n", learned_filename);
  }
  if (kernel->test_filename != NULL) {
    fprintf(outfile, "# test_file=%s\n", kernel->test_filename);
  }
  if (kernel->self_train_filename != NULL) {
    fprintf(outfile, "# self_train_file=%s\n", kernel->self_train_filename);
  }
  if (kernel->self_test_filename != NULL) {
    fprintf(outfile, "# self_test_file=%s\n", kernel->self_test_filename);
  }

  /* Print the properties of the kernel. */
  fprintf(outfile, "# matrix_from_file=%s", 
	  boolean_to_string(kernel->matrix_from_file));
  fprintf(outfile, " zero_mean=%s", 
	  boolean_to_string(kernel->zero_mean_row));
  fprintf(outfile, " variance_one=%s",
	  boolean_to_string(kernel->variance_one));
  fprintf(outfile, " normalize=%s\n",
	  boolean_to_string(kernel->normalize));
  fprintf(outfile, "# constant=%g", kernel->constant);
  fprintf(outfile, " coefficient=%g", kernel->coefficient);
  fprintf(outfile, " power=%g\n", kernel->power);
  fprintf(outfile, "# radial=%s", boolean_to_string(kernel->radial));
  fprintf(outfile, " width_factor=%g", kernel->width_factor);
  fprintf(outfile, " two_squared_width=%g\n", kernel->two_squared_width);

  if (print_time) {
    fprintf(outfile, "# host=%s", hostname());
    fprintf(outfile, " date=%s\n", date_and_time());
  }
}

/****************************************************************************
 * Main procedure
 ****************************************************************************/



#ifdef TEST_MAIN
#include "cmdline.h"
#include <time.h>

int main(int argc, char *argv[])
{
  // Command line options.
  char*     learned_filename = NULL; // Training set weights or eigenvectors.
  FILE*     learned_file = NULL;
  BOOLEAN_T format_line = FALSE;     // Use format line in RDB files?
  BOOLEAN_T kernel_output = FALSE;   // Output kernel matrix and stop?
  BOOLEAN_T print_time = TRUE;       // Print timing data to output header.
  int       output_precision = 0;    // Number of digits in output.

  // Data structures.
  double         bias = 0.0;             // Bias term for discriminant.
  KERNEL_T*      kernel = allocate_kernel(FALSE);  // Kernel function.
  RDB_MATRIX_T*  learned_matrix = NULL;  // Weights / eigenvectors of training
  ARRAY_T*       classifications = NULL; // Predicted classifications.
  ARRAY_T*       discriminants = NULL;   // Discriminants of test set.
  RDB_MATRIX_T*  output_matrix = NULL;   // All output.
  verbosity = NORMAL_VERBOSE;

  // Get the command line arguments.
  DO_STANDARD_COMMAND_LINE
    (3,
     DATA_OPTN(1, train, <filename> (required unless test set is a kernel matrix), 
	       kernel->train_filename = _OPTION_);
     DATA_OPTN(1, learned, <filename> (required), learned_filename = _OPTION_);
     DATA_OPTN(1, test, <filename> (required),
	       kernel->test_filename = _OPTION_);
     DATA_OPTN(1, selftrain, <filename>, 
	       kernel->self_train_filename = _OPTION_);
     DATA_OPTN(1, selftest, <filename>, 
	       kernel->self_test_filename = _OPTION_);
     FLAG_OPTN(1, rdb, format_line = TRUE);
     FLAG_OPTN(1, kernelout, kernel_output = TRUE);
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
  if (kernel->test_filename == NULL) {
    die("No test file given.");
  }
  if (learned_filename == NULL) {
    if (CLASSIFY) {
      die("No training set weights given.\n");
    } else if (PROJECT) {
      die("No training set eigenvectors given.\n");
    } else {
      die("No hyperplane given.\n");
    }
  }

  // Read the kernel parameters from the header of the learned file.
  // N.B. If we open the file once, and pass the file handle, then for
  // some reason 'project' no longer works.  -- WSN 4/4/03
  read_kernel_parameters(learned_filename, &bias, NULL, kernel);

  // Check for more files.
  if ((kernel->train_filename == NULL) && (!kernel->matrix_from_file)) {
    die("No training file given.");
  }

  if (kernel->matrix_from_file) {
    if (kernel->normalize || kernel->radial) {
      if (kernel->self_train_filename == NULL) {
	die("No self-train file given.");
      }
      if (kernel->self_test_filename == NULL) {
	die("No self-test file given.");
      }
    } else {
      if (kernel->train_filename) {
	die("Unneeded training file entered.");
      }
    }
  }

  // Read the associated weights or eigenvectors.
  if (open_file(learned_filename, "r", TRUE, "learned", "learned values",
		&learned_file) == 0)
    exit(1);
  learned_matrix = read_rdb_matrix(format_line, "NaN", learned_file);

  // Read the data files.
  read_data_files(format_line, kernel);

  // Compute the kernel;
  compute_base_kernel(kernel, NULL, NULL);

  // Turbo-ize the kernel matrix.
  transform_kernel_matrix(FALSE,
			  kernel, 
			  NULL,
			  kernel->self_train,
			  kernel->self_test,
			  kernel->kernel_matrix);

  // Stop here if requested.
  if (kernel_output) {
    print_rdb_matrix(kernel->kernel_matrix_rdb, format_line, 
		     output_precision + 2, output_precision, stdout);
    exit(0);
  }

  if (CLASSIFY) {
    ARRAY_T* matrix_column;

    // Classify the test data.
    if (verbosity >= NORMAL_VERBOSE) {
      fprintf(stderr, "Classifying the test data.\n");
    }
    matrix_column = get_matrix_column(1, get_raw_matrix(learned_matrix));
    classify_list(get_num_rows(kernel->test_matrix),
		  bias,
		  matrix_column,
		  kernel->kernel_matrix,
		  &classifications,
		  &discriminants);
    free_array(matrix_column);

    // Create the output matrix.
    output_matrix = allocate_rdb_matrix(get_num_rows(kernel->test_matrix),
					2, NULL);
    set_corner_string(get_corner_string(kernel->kernel_matrix_rdb),
		      output_matrix);
    set_row_names(get_row_names(kernel->kernel_matrix_rdb), output_matrix);
    //    set_row_names(get_row_names(kernel->test_matrix_rdb), output_matrix);
    add_rdb_column("classification", classifications, output_matrix);
    add_rdb_column("discriminant", discriminants, output_matrix);

  } else {

    // Allocate space for the projected data.
    output_matrix
      = allocate_rdb_matrix(get_num_rows(kernel->test_matrix),
			    get_num_cols(get_raw_matrix(learned_matrix)),
			    NULL);
    set_corner_string(get_corner_string(kernel->kernel_matrix_rdb),
		      output_matrix);
    set_row_names(get_row_names(kernel->test_matrix_rdb), output_matrix);
    set_col_names(get_col_names(learned_matrix), output_matrix);

    // Project the data onto the eigenvectors.
    if (verbosity >= NORMAL_VERBOSE) {
      fprintf(stderr, "Projecting the test data.\n");
    }
    project_data(kernel->kernel_matrix,
		 get_raw_matrix(learned_matrix),
		 get_raw_matrix(output_matrix));

  }

  // Print the results to stdout.
  print_output_header(learned_filename,
		      kernel,
		      print_time,
		      stdout);
  //  print_rdb_matrix(output_matrix, format_line, 6, 4, stdout);
  print_rdb_matrix(output_matrix, format_line, output_precision + 2, 
		   output_precision, stdout);
  
#ifdef DEBUG
  // Free dynamic memory.
  free_kernel(kernel);
  free_rdb_matrix(learned_matrix);
  free_array(classifications);
  free_array(discriminants);
  free_rdb_matrix(output_matrix);
#endif
  return(0);
}
#endif


/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
 
