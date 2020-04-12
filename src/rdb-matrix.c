/*****************************************************************************
 * FILE: rdb-matrix.c
 * AUTHOR: William Stafford Noble
 * CREATE DATE: 4/19/99
 * PROJECT: shared
 * COPYRIGHT: 1999-2001, Columbia University
 * VERSION: $Revision: 1.9 $
 * DESCRIPTION: A matrix data structure with labels on rows and columns.
 *****************************************************************************/
#include "rdb-matrix.h"
#include "string-list.h"
#include "matrix.h"
#include "array.h"
#include "utils.h"
#include <string.h>
#include <stdio.h>
#include <assert.h>


#ifndef VERBOSITY
#define VERBOSITY
VERBOSE_T verbosity;
#endif


/***************************************************************************
 * Define the RDB matrix type.
 ***************************************************************************/
struct rdb_matrix_t {
  char*          corner_string; /* Upper left corner. */
  STRING_LIST_T* row_names;
  STRING_LIST_T* col_names;
  MATRIX_T*      matrix;
};

/***********************************************************************
 * Allocate memory for an RDB matrix.
 *
 * Optional third argument is a raw matrix to be stored in the new RDB
 * matrix.
 ***********************************************************************/
RDB_MATRIX_T* allocate_rdb_matrix
  (int num_rows,
   int num_cols,
   MATRIX_T* matrix)
{
  RDB_MATRIX_T* return_value;
  assert(num_rows != 0 && num_cols != 0);

  /* Allocate the new RDB matrix. */
  return_value = (RDB_MATRIX_T*)mymalloc(sizeof(RDB_MATRIX_T));
  return_value->row_names = new_string_list();
  return_value->col_names = new_string_list();
  if (matrix == NULL) {
    return_value->matrix = allocate_matrix(num_rows, num_cols);
  } else {
    return_value->matrix = matrix;
  }

  return(return_value);
}

/***********************************************************************
 * Convert a given matrix to an RDB matrix.
 ***********************************************************************/
RDB_MATRIX_T* rdbize_matrix
  (char*          corner_string,
   STRING_LIST_T* row_names,
   STRING_LIST_T* col_names,
   MATRIX_T*      matrix)
{
  RDB_MATRIX_T* new_matrix;
  assert(get_num_rows(matrix) == get_num_strings(row_names));
  assert(get_num_cols(matrix) == get_num_strings(col_names));

  new_matrix = allocate_rdb_matrix(get_num_rows(matrix),
				   get_num_cols(matrix),
				   matrix);

  set_corner_string(corner_string, new_matrix);
  set_row_names(row_names, new_matrix);
  set_col_names(col_names, new_matrix);
  
  return(new_matrix);
}

/***********************************************************************
 * Get or set the various pieces of the RDB matrix.
 ***********************************************************************/
STRING_LIST_T* get_row_names
  (RDB_MATRIX_T* rdb_matrix)
{
  if (rdb_matrix == NULL) {
    die("Attempted to access null matrix.");
  }
  return(rdb_matrix->row_names);
}

STRING_LIST_T* get_col_names
  (RDB_MATRIX_T* rdb_matrix)
{
  if (rdb_matrix == NULL) {
    die("Attempted to access null matrix.");
  }
  return(rdb_matrix->col_names);
}

void set_row_names
  (STRING_LIST_T* row_names,
   RDB_MATRIX_T* rdb_matrix)
{
  int num_rows;
  int i_row;

  if (rdb_matrix == NULL) {
    die("Attempted to access null matrix.");
  }

  num_rows = get_num_rows(get_raw_matrix(rdb_matrix));
  if (get_num_strings(row_names) != num_rows) {
    die("Adding %d row names to a matrix with %d rows.",
	get_num_strings(row_names), num_rows);
  }

  /* Add the strings. */
  for (i_row = 0; i_row < num_rows; i_row++) {
    add_string(get_nth_string(i_row, row_names), rdb_matrix->row_names);
  }
}

void set_col_names
  (STRING_LIST_T* col_names,
   RDB_MATRIX_T*  rdb_matrix)
{
  int i_col;
  int num_cols;

  if (rdb_matrix == NULL) {
    die("Attempted to access null matrix.");
  }

  num_cols = get_num_cols(get_raw_matrix(rdb_matrix));
  if (get_num_strings(col_names) != num_cols) {
    die("Adding %d column names to a matrix with %d columns.",
	get_num_strings(col_names), get_num_cols(get_raw_matrix(rdb_matrix)));
  }

  /* Add the strings. */
  for (i_col = 0; i_col < num_cols; i_col++) {
    add_string(get_nth_string(i_col, col_names), rdb_matrix->col_names);
  }
}

char* get_corner_string
  (RDB_MATRIX_T* rdb_matrix)

{
  if (rdb_matrix == NULL) {
    die("Attempted to access null matrix.");
  }
  return(rdb_matrix->corner_string);
}

void set_corner_string
  (char* corner_string,
   RDB_MATRIX_T* rdb_matrix)
{
  if (rdb_matrix == NULL) {
    die("Attempted to access null matrix.");
  }
  copy_string(&(rdb_matrix->corner_string), corner_string);
}


MATRIX_T* get_raw_matrix
  (RDB_MATRIX_T* rdb_matrix)
{
  if (rdb_matrix == NULL) {
    die("Attempted to access null matrix.");
  }
  return(rdb_matrix->matrix);
}

void set_raw_matrix
  (MATRIX_T*     matrix,
   RDB_MATRIX_T* rdb_matrix)
{
  if (rdb_matrix == NULL) {
    die("Attempted to access null matrix.");
  }
  free_matrix(rdb_matrix->matrix);
  rdb_matrix->matrix = matrix;
}

/***********************************************************************
 * Add a column to a matrix.
 ***********************************************************************/
void add_rdb_column
  (char*         col_name,
   ARRAY_T*      new_column,
   RDB_MATRIX_T* rdb_matrix)
{
  int i_col;
  int i_row;
  MATRIX_T* raw_matrix;
  int num_rows;
  MTYPE this_value;

  if (rdb_matrix == NULL) {
    die("Attempted to access null matrix.");
  }

  /* Add the new name. */
  add_string(col_name, rdb_matrix->col_names);

  /* Find out how many columns we have. */
  i_col = get_num_strings(rdb_matrix->col_names) - 1;

  /* Get the raw matrix. */
  raw_matrix = get_raw_matrix(rdb_matrix);
  assert(i_col < get_num_cols(raw_matrix));

  /* Copy from the array to the matrix. */
  num_rows = get_array_length(new_column);
  assert(num_rows == get_num_rows(raw_matrix));
  for (i_row = 0; i_row < num_rows; i_row++) {
    this_value = get_array_item(i_row, new_column);
    set_matrix_cell(i_row, i_col, this_value, raw_matrix);
  }
}

/***********************************************************************
 * Add a row or column name to a matrix.
 ***********************************************************************/
void add_row_name
  (char*         row_name,
   RDB_MATRIX_T* rdb_matrix)
{
  if (rdb_matrix == NULL) {
    die("Attempted to access null matrix.");
  }
  add_string(row_name, rdb_matrix->row_names);
}

void add_col_name
  (char*         col_name,
   RDB_MATRIX_T* rdb_matrix)
{
  if (rdb_matrix == NULL) {
    die("Attempted to access null matrix.");
  }
  add_string(col_name, rdb_matrix->col_names);
}

/***********************************************************************
 * Is the given line a comment?
 ***********************************************************************/
BOOLEAN_T is_comment
  (char* one_line)
{
  int  i_char;
  char this_char;

  /* Find the first non-space character. */
  i_char = 0;
  this_char = one_line[i_char];
  while (this_char == ' ') {
    this_char = one_line[i_char++];
  }

  /* Is it a hash mark? */
  return(this_char == '#');
}

/***********************************************************************
 * Read one row from a file, checking for long lines.
 ***********************************************************************/
void read_one_row
  (FILE* infile,
   int   length,
   char* this_row)
{
  char*     fgets_result;       /* Error indicator for 'fgets' function. */

  fgets_result = fgets(this_row, length, infile);

  if (this_row[strlen(this_row) - 1] != '\n') {
    if ((int)strlen(this_row) <= length) { // this occurs if there are backspace characters, for example.
      die("Error encountered while reading file, possible illegal character: row reads %s", this_row);
    } else {
      die("Matrix lines too long (more than  %d), or illegal characters encountered. ", length);
    }
  }
}

/***********************************************************************
 * Read one row from an RDB file.
 ***********************************************************************/
void parse_rdb_row
(char*    this_row, 
 char*    missing,
 char*    missing_marker,
 int      i_row,
 int      num_cols,
 char*    row_name, 
 ARRAY_T* row_array)
{
  char*     string_ptr = NULL;  // Storage for strtok function.
  int       i_column;           // Index of the current column.
  int       num_scanned;        // Error indicator for fscanf function.
  MTYPE     one_value;          // One value read from the file.

  /* Read the row name and store it. */
  string_ptr = strtok(this_row, "\t");

  /* clean up any stray dos linefeeds */
  if (string_ptr[strlen(string_ptr) - 1] == '\r') {
    string_ptr[strlen(string_ptr) - 1] = '\0';
  }

  // Store the name of this row.
  strcpy(row_name, string_ptr);

  /* Read the row. */
  // the +1 is so we read past the supposed end as a check.
  for (i_column = 0; i_column < num_cols+1; i_column++) { 
      
    /* This is an important check: if the number of
       headings is too few, we won't read in all of the data, which
       is almost not what is wanted. So check to see that there is
       no data there. Typically this is caused by a missing corner string.
       (PP) */
    if (i_column == num_cols) {
      string_ptr = strtok(NULL, "\t");
      if (string_ptr != NULL) {
	die("Data unexpectedly found past end of headings: (%d,%d). Check that all columns in the file have a heading, including the 'corner'.\n", 
	    i_row, i_column);
      } else {
	// good, we're done
	break;
      }
    }

    /* Read the first value. */
    string_ptr = strtok(NULL, "\t");

    /* Make sure we got a valid string. */
    if (string_ptr == NULL) {
      die("No entry found at (%d,%d).\n", i_row, i_column);
    }

    /* Extract a number */
    num_scanned = sscanf(string_ptr, MSCAN, &one_value);

    if (num_scanned == EOF) {
      die("Premature end-of-file at (%d,%d).\n", i_row, i_column);
    }

    /* If we didn't get a number, deal with missing values. */
    if (num_scanned == 0) {
      num_scanned = sscanf(string_ptr, "%s", missing);
      if ((missing_marker == NULL) ||
	  (num_scanned == 0) ||
	  (strcmp(missing, missing_marker) != 0)) {
	die("Invalid missing value (%s) at (%d,%d).\n", missing, i_row,
	    i_column);
      }
      if (strcmp(ATYPENAME, "double") != 0) {
	die("Missing values are only permitted in arrays of type double.\n");
      }
      one_value = NaN();
    }

    /* Store the value. */
    set_array_item(i_column, one_value, row_array);
  }
}


/***********************************************************************
 * Read an RDB file into a matrix.
 ***********************************************************************/
RDB_MATRIX_T* read_rdb_matrix
  (BOOLEAN_T format_line,
   char* missing_marker,
   FILE* infile)
{
  MATRIX_T* matrix;             /* The matrix to be read. */
  char*     corner_string;      /* Upper left corner of matrix. */
  STRING_LIST_T* row_names;     /* Row names in the matrix. */
  STRING_LIST_T* col_names;     /* Column names in the matrix. */
  char*     this_row = NULL;    /* One row of the matrix. */
  int       i_row;              /* Index of the current row. */
  int       num_rows;           /* Total number of rows. */
  int       num_cols;           /* Total number of rows. */
  char*     string_ptr = NULL;  /* Storage for strtok function. */
  ARRAY_T*  row_array = NULL;   /* One row of the matrix. */
  RDB_MATRIX_T* return_value;   /* The RDB matrix being created. */
  char*     missing = NULL;     /* Storage space for missing values. */
  char*     row_name = NULL;    /* ID associated with this row. */

  this_row = (char*)malloc(sizeof(char)*MAX_ROW);
  missing = (char*)malloc(sizeof(char)*MAX_ROW);
  row_name = (char*)malloc(sizeof(char)*MAX_ROW);

  /* Make sure the file is valid. */
  if (infile == NULL) {
    die("Attempted to read matrix from null file.");
  }

  /* Create the row names and column names lists. */
  row_names = new_string_list();
  col_names = new_string_list();

  /* Read the first row. */
  read_one_row(infile, MAX_ROW, this_row);

  /* Keep reading till we get past the comments. */
  while (is_comment(this_row)) {
    /* fprintf(stderr, "Skipping: %s", this_row); */
    read_one_row(infile, MAX_ROW, this_row);
  }

  /* Store the name of the first column. (if the corner string is
     missing, no error can be reported here, the next token is
     taken! This error is caught later -- PP) */
  string_ptr = strtok(this_row, "\t");
  copy_string(&corner_string, string_ptr); 

  /* Store the names of the columns. */
  for (string_ptr = strtok(NULL, "\t"); string_ptr != NULL;
       string_ptr = strtok(NULL, "\t")) {
 
    /* Remove EOL. */
    if (string_ptr[strlen(string_ptr) - 1] == '\n') {
      string_ptr[strlen(string_ptr) - 1] = '\0';
    }

    /* clean up any stray dos linefeed */
    if (string_ptr[strlen(string_ptr) - 1] == '\r') {
      string_ptr[strlen(string_ptr) - 1] = '\0';
    }
    
    /* Store string if it is non-empty. */
    if (strcmp(string_ptr, "") != 0) {
      add_string(string_ptr, col_names);
    }
  }
  num_cols = get_num_strings(col_names);

  /* Allocate the matrix. */
  matrix = allocate_matrix(0, num_cols);

  /* Allocate one row. */
  row_array = allocate_array(num_cols);

  /* Skip the format line, if necessary. */
  if (format_line) {
    read_one_row(infile, MAX_ROW, this_row);
  }

  /* Read the matrix. */
  for (i_row = 0; ; i_row++) {

    /* Read the next line, stopping if it's empty. */
    if (fgets(this_row, MAX_ROW, infile) == NULL) {
      break;
    }

    // Convert this line into a row.
    parse_rdb_row(this_row,
		  missing,
		  missing_marker, 
		  i_row, 
		  num_cols, 
		  row_name, 
		  row_array);

    // Store the name of the row.
    add_string(row_name, row_names);

    // Add this row to the matrix. */
    grow_matrix(row_array, matrix);
  
  }
  //  num_rows = i_row - 1; Modified by WSN 4/29/04
  num_rows = i_row;

  /* Assemble it all into an RDB matrix. */
  return_value = allocate_rdb_matrix(num_rows, num_cols, matrix);
  set_corner_string(corner_string, return_value);
  set_row_names(row_names, return_value);
  set_col_names(col_names, return_value);

  /* Free local dynamic memory. */
  myfree(corner_string);
  free_array(row_array);
  free_string_list(row_names);
  free_string_list(col_names);
  free(missing);
  free(this_row);
  return(return_value);
}

/***********************************************************************
 * Write a named array to a file.
 ***********************************************************************/
void print_rdb_array
  (STRING_LIST_T* names,         /* Structure storing names */ 
   char*          label1,        /* Label for column 1. */
   char*          label2,        /* Label for column 2. */
   ARRAY_T*       array,         /* The array to be printed. */
   int            width,         /* Width of each cell. */
   int            precision,     /* Precision of each cell. */
   FILE*          outfile)       /* File to which to write. */
{
  
  int   i_item;
  int   num_items;
  double item;

  fprintf(outfile, "%s\t%s\n", label1, label2);
  fprintf(outfile, "5S\t5N\n");

  num_items = get_array_length(array);
  for (i_item = 0; i_item < num_items; i_item++) {
    item = get_array_item(i_item, array);
    fprintf(outfile, "%s\t%*.*f\n", get_nth_string(i_item, names), width, 
            precision, item);
  }
}

/***********************************************************************
 * Write a labeled matrix in RDB format.
 ***********************************************************************/
void print_rdb_matrix
  (RDB_MATRIX_T*  rdb_matrix,    /* Matrix to be printed. */
   BOOLEAN_T      format_line,   /* Print format line? */
   int            width,         /* Width of each cell. */
   int            precision,     /* Precision of each cell. */
   FILE*          outfile)       /* File to which to write. */
{
  int   num_cols;
  int   num_rows;
  int   i_col;
  int   i_row;
  MTYPE item;

  /* Print the column names. */
  fprintf(outfile, "%s", rdb_matrix->corner_string);
  num_cols = get_num_cols(rdb_matrix->matrix);
  for (i_col = 0; i_col < num_cols; i_col++) {
    fprintf(outfile, "\t%s", get_nth_string(i_col, rdb_matrix->col_names));
  }
  fprintf(outfile, "\n");

  /* If requested, print the column widths. */
  if (format_line) {
    fprintf(outfile, "10S");
    for (i_col = 0; i_col < num_cols; i_col++) {
      fprintf(outfile, "\t%dN", width);
    }
    fprintf(outfile, "\n");
  }

  /* Print the matrix. */
  num_rows = get_num_rows(rdb_matrix->matrix);
  for (i_row = 0; i_row < num_rows; i_row++) {
    fprintf(outfile, "%s", get_nth_string(i_row, rdb_matrix->row_names));
    for (i_col = 0; i_col < num_cols; i_col++) {
      item = get_matrix_cell(i_row, i_col, rdb_matrix->matrix);
      // Make sure NaN is the same case on all architectures.
      if (isnan(item)) {
	fprintf(outfile, "\t%*s", width, "NaN");
      } else {
	fprintf(outfile, "\t%*.*g", width, precision, item);
      }
    }
    fprintf(outfile, "\n");    
  }
}


/***********************************************************************
 * Free an RDB matrix.
 ***********************************************************************/
void free_rdb_matrix
  (RDB_MATRIX_T* rdb_matrix)
{
  if (rdb_matrix != NULL) {
    if (rdb_matrix->corner_string != NULL) {
      myfree(rdb_matrix->corner_string);
    }
    free_string_list(rdb_matrix->row_names);
    free_string_list(rdb_matrix->col_names);
    free_matrix(rdb_matrix->matrix);
    myfree(rdb_matrix);
  }
}


#ifdef RDB_MATRIX_MAIN
#include "cmdline.h"
#include "linear-algebra.h"
#include <time.h>

/************************************************************************
 * Set the value on the diagonal of a given square matrix.
 ************************************************************************/
static void set_diagonal_matrix
(double    value,
 MATRIX_T* matrix)
{
  int num_cols;
  int i_col;

  num_cols = get_num_cols(matrix);

  // Verify that the matrix is square.
  if (num_cols != get_num_rows(matrix)) {
    die("Setting diagonal in non-square matrix (%d by %d).",
	get_num_rows(matrix), num_cols);
  }

  for (i_col = 0; i_col < num_cols; i_col++) {
    set_matrix_cell(i_col, i_col, value, matrix);
  }
}


/************************************************************************
 * Compute the variance of each column in an n x m matrix and stores
 * them in a 1 x m matrix.
 ************************************************************************/
static void compute_column_variances
  (MATRIX_T* matrix,
   MATRIX_T* variance_matrix)
{
  int      num_cols;
  int      i_col;
  ARRAY_T* this_column;
  double   variance;

  num_cols = get_num_cols(matrix);
  for (i_col = 0; i_col < num_cols; i_col++) {
    this_column = get_matrix_column(i_col, matrix);
    variance = array_variance(this_column);
    set_matrix_cell(0, i_col, variance, variance_matrix);
    free_array(this_column);
  }
}


/************************************************************************
 * Divide each row or column in a matrix by the corresponding element
 * from a given array.
 ************************************************************************/
static void vector_divide
  (BOOLEAN_T row,
   ARRAY_T* array,
   MATRIX_T* matrix)
{
  int num_rows;
  int i_row;
  int num_cols;
  int i_col;
  double denominator;

  num_rows = get_num_rows(matrix);
  num_cols = get_num_cols(matrix);
  if (row) {
    assert(get_array_length(array) == num_rows);
  } else {
    assert(get_array_length(array) == num_cols);
  }

  for (i_col = 0; i_col < num_cols; i_col++) {
    for (i_row = 0; i_row < num_rows; i_row++) {

      if (row) {
	denominator = get_array_item(i_row, array);
      } else {
	denominator = get_array_item(i_col, array);
      }

      // If there's nothing in the array, skip it.
      if (denominator != 0.0) {
	set_matrix_cell(i_row, i_col, 
			get_matrix_cell(i_row, i_col, matrix) / denominator,
			matrix);
      }
    }
  }
}

#define DEFAULT_WIDTH 8      // Total number of digits in output values.
#define DEFAULT_PRECISION 5  // Number of digits to right of zero.
#define DIFFUSION_ITERATIONS 20 // Number of iterations of diffusion.

int main(int argc, char *argv[])
{
  char*         matrix1_filename; /* Name of the file containing the matrix. */
  char*         matrix2_filename; /* Name of the file containing the matrix. */
  FILE*         matrix1_file;     /* Pointer to the input file. */
  FILE*         matrix2_file;     /* Pointer to the input file. */
  char*         operation;        /* Operation to perform on the matrix. */
  double        value;            /* Value to be passed to operation. */
  int           width;            /* Width of output values. */
  int           precision;        /* Precision of output values. */
  BOOLEAN_T     rdb_format;       /* Is the file in RDB format? */
  BOOLEAN_T     format_line;      /* Include format line in RDB? */
  long          seed;             /* Seed for random number generator. */
  RDB_MATRIX_T* rdb_matrix1;      /* The first matrix. */
  RDB_MATRIX_T* rdb_matrix2;      /* The second matrix. */
  RDB_MATRIX_T* rdb_new_matrix;   /* Output matrix. */
  MATRIX_T*     matrix1;          /* Non-rdb version of the above. */
  MATRIX_T*     matrix2;          /* Ditto. */
  MATRIX_T*     new_matrix;       /* Ditto. */
  ARRAY_T*      new_array;        /* Array used for temporary storage. */

  /* Initialize data structures. */
  rdb_matrix1 = NULL;
  rdb_matrix2 = NULL;
  rdb_new_matrix = NULL;
  matrix1 = NULL;
  matrix2 = NULL;
  new_matrix = NULL;

  /***** COMMAND LINE PROCESSING *****/

  /* Set the defaults. */
  matrix1_filename = NULL;
  matrix2_filename = NULL;
  operation = NULL;
  value = 0.000001;
  width = DEFAULT_WIDTH;
  precision = DEFAULT_PRECISION;
  rdb_format = TRUE;
  format_line = FALSE;
  seed = time(0);
  verbosity = NORMAL_VERBOSE;

  /* Get the command line arguments. */
  DO_STANDARD_COMMAND_LINE
    (1,
     DATA_OPTN(1, matrix1, <filename> (required), matrix1_filename = _OPTION_);
     DATA_OPTN(1, matrix2, <filename>, matrix2_filename = _OPTION_);
     DATA_OPTN(1, operation, 
	       none|size|randomize|covariance|square|invert|binarize|\n
               transpose|eigenvectors|eigenvalues|eigenvector1|symdiag|\n
               adddiag|getrowsums|getcolsums|getdiag|add|normalizerow|\n
               normalizecol|rescalerow|rescalecol|normalize|zeromeanrow|\n
               zeromeancol|varone|colvars|scalarmult|dotmult|posdef|\n
               row-correlation|col-correlation|jacard|shuffle|alignment|\n
               diagonal|label-align|scalediag|euclidean|center|setdiag|\n
               diffusion|euclidean|int-dim|simple-diffusion|modav\n
              (default=none),
	       operation = _OPTION_);
     DATA_OPTN(1, width, <value> (default=8), width = atoi(_OPTION_));
     DATA_OPTN(1, precision, <value> (default=5), precision = atoi(_OPTION_));
     DATA_OPTN(1, value, <value> (default=0.000001), value = atof(_OPTION_));
     DATA_OPTN(1, seed, <value> (default from clock), seed = atoi(_OPTION_));
     FLAG_OPTN(1, raw, rdb_format = FALSE);
     FLAG_OPTN(1, rdb, format_line = TRUE);
     DATA_OPTN(1, verbose, 1|2|3|4|5 (default=2), 
	       verbosity = (VERBOSE_T)atoi(_OPTION_));
     );

  /* Make sure we got the first file. */
  if (matrix1_filename == NULL) {
    die("No matrix filename given.\n");
  }

  /***** READING INPUT FILES *****/

  /* Open the first input file. */
  if (open_file(matrix1_filename, "r", TRUE, "matrix1", "matrix1",
		&matrix1_file) == 0)
    exit(1);

  /* Read the matrix from the file. */
  if (rdb_format) {
    rdb_matrix1 = read_rdb_matrix(format_line, NULL, matrix1_file);
    matrix1 = get_raw_matrix(rdb_matrix1);
  } else {
    matrix1 = read_matrix(matrix1_file);
  }
  fclose(matrix1_file);

  /* Do it again for the second matrix. */
  if (matrix2_filename != NULL) {
    if (open_file(matrix2_filename, "r", TRUE, "matrix", "matrix",
		  &matrix2_file) == 0)
      exit(1);
    if (rdb_format) {
      rdb_matrix2 = read_rdb_matrix(format_line, NULL, matrix2_file);
      matrix2 = get_raw_matrix(rdb_matrix2);
    } else {
      matrix2 = read_matrix(matrix2_file);
    }
    fclose(matrix2_file);
  }

  /****** PERFORM OPERATIONS *****/
  
  if ((operation == NULL) || (strcmp(operation, "none") == 0)) {

    /* Don't do anything. */
    new_matrix = allocate_matrix(get_num_rows(matrix1), get_num_cols(matrix1));
    copy_matrix(matrix1, new_matrix);

  } else if (strcmp(operation, "size") == 0) {
    printf("%d by %d\n", get_num_rows(matrix1), get_num_cols(matrix1));
    free_matrix(matrix1);
    return(0);

  } else if (strcmp(operation, "randomize") == 0) {

    /* Randomize the matrix. */
    my_srand(seed);
    new_matrix = allocate_matrix(get_num_rows(matrix1), get_num_cols(matrix1));
    copy_matrix(matrix1, new_matrix);
    randomize_matrix(1.0, new_matrix);

  } else if (strcmp(operation, "covariance") == 0) {

    /* Compute the covariance matrix. */
    new_matrix = compute_covariance(matrix1);

  } else if (strcmp(operation, "square") == 0) {

    /* Square the matrix. */
    new_matrix = matrix_multiply(matrix1, matrix1, NULL);

  } else if (strcmp(operation, "invert") == 0) {

    /* Invert the matrix. */
    new_matrix = invert_matrix(matrix1);
    if (new_matrix == NULL) {
      exit(0);
    }

  } else if (strcmp(operation, "binarize") == 0) {

    // Convert matrix to binary, using a given threshold.
    new_matrix = allocate_matrix(get_num_rows(matrix1), get_num_cols(matrix1));
    copy_matrix(matrix1, new_matrix);
    binarize_matrix(value, new_matrix);

  } else if (strcmp(operation, "transpose") == 0) {

    /* Transpose the matrix. */
    new_matrix = transpose_matrix(matrix1);

  } else if (strcmp(operation, "eigenvectors") == 0) {

    /* Compute the matrix of eigenvectors. */
    find_symmetric_eigenvectors(matrix1, &new_array, &new_matrix);
    free_array(new_array);

  } else if (strcmp(operation, "eigenvalues") == 0) {

    /* Compute the array of eigenvalues. */
    find_symmetric_eigenvectors(matrix1, &new_array, &new_matrix);
    free_matrix(new_matrix);
    new_matrix = array_to_matrix(TRUE, new_array);

  } else if (strcmp(operation, "eigenvector1") == 0) {

    // Find the first eigenvector.
    new_array = find_first_eigenvector(matrix1);
    new_matrix = array_to_matrix(FALSE, new_array);
    free_array(new_array);

  } else if (strcmp(operation, "symdiag") == 0) {

    /* Deal with optional second argument. */
    if (matrix2 == NULL) {
      matrix2 = allocate_matrix(get_num_rows(matrix1),
				get_num_cols(matrix1));
      copy_matrix(matrix1, matrix2);
    }

    /* Average across the diagonal. */
    new_matrix = average_across_diagonal(matrix1, matrix2);

  } else if (strcmp(operation, "adddiag") == 0) {

    /* Add to the diagonal. */
    new_matrix = allocate_matrix(get_num_rows(matrix1), get_num_cols(matrix1));
    copy_matrix(matrix1, new_matrix);
    add_to_diagonal(value, new_matrix);

  } else if (strcmp(operation, "getrowsums") == 0) {

    // Extract the sum of the rows.
    new_array = get_matrix_row_sums(matrix1);
    new_matrix = array_to_matrix(FALSE, new_array);
    free_array(new_array);

  } else if (strcmp(operation, "getcolsums") == 0) {

    // Extract the sum of the columns.
    new_array = get_matrix_col_sums(matrix1);
    new_matrix = array_to_matrix(FALSE, new_array);
    free_array(new_array);

  } else if (strcmp(operation, "getdiag") == 0) {

    /* Extract the diagonal. */
    new_array = extract_diagonal(matrix1);
    new_matrix = array_to_matrix(FALSE, new_array);
    free_array(new_array);
    
  } else if (strcmp(operation, "add") == 0) {
    
    /* Add two matrices. */
    new_matrix = allocate_matrix(get_num_rows(matrix1), get_num_cols(matrix1));
    copy_matrix(matrix1, new_matrix);
    sum_matrices(matrix2, new_matrix);

  } else if ((strcmp(operation, "normalizerow") == 0) ||
	     (strcmp(operation, "normalizecol") == 0)) {

    // If matrix2 was not given, then use matrix1.
    if (matrix2 == NULL) {
      matrix2 = matrix1;
      rdb_matrix2 = rdb_matrix1;
    }
    
    // Normalize the rows or columns.
    if (strcmp(operation, "normalizerow") == 0) {
      vector_divide(TRUE, get_matrix_row_sums(matrix1), matrix2);
    } else {
      vector_divide(FALSE, get_matrix_col_sums(matrix1), matrix2);
    }
    new_matrix = allocate_matrix(get_num_rows(matrix2), get_num_cols(matrix2));
    copy_matrix(matrix2, new_matrix);

  } else if ((strcmp(operation, "rescalerow") == 0) ||
	     (strcmp(operation, "rescalecol") == 0)) {
    
    // Rescale the row or column so that it ranges from -1 to 1.
    if (strcmp(operation, "rescalerow") == 0) {
      rescale_matrix_rows(-1.0, 1.0, matrix1);
    } else {
      rescale_matrix_cols(-1.0, 1.0, matrix1);
    }
    new_matrix = allocate_matrix(get_num_rows(matrix1), get_num_cols(matrix1));
    copy_matrix(matrix1, new_matrix);


  } else if (strcmp(operation, "normalize") == 0) {

    /* Normalize the matrix. */
    new_matrix = allocate_matrix(get_num_rows(matrix1), get_num_cols(matrix1));
    copy_matrix(matrix1, new_matrix);
    normalize_matrix(value, new_matrix);

  } else if (strcmp(operation, "zeromeanrow") == 0) {

    /* Set the mean of the matrix rows to zero. */
    new_matrix = allocate_matrix(get_num_rows(matrix1), get_num_cols(matrix1));
    copy_matrix(matrix1, new_matrix);
    zero_mean_matrix_rows(new_matrix);

  } else if (strcmp(operation, "zeromeancol") == 0) {

    /* Set the mean of the matrix columns to zero. */
    new_matrix = allocate_matrix(get_num_rows(matrix1), get_num_cols(matrix1));
    copy_matrix(matrix1, new_matrix);
    zero_mean_matrix_cols(new_matrix);

  } else if (strcmp(operation, "varone") == 0) {

    /* Set the variance of the matrix to one. */
    new_matrix = allocate_matrix(get_num_rows(matrix1), get_num_cols(matrix1));
    copy_matrix(matrix1, new_matrix);
    variance_one_matrix_rows(new_matrix);

  } else if (strcmp(operation, "colvars") == 0) {
    
    /* Compute the variance of each column. */
    new_matrix = allocate_matrix(1, get_num_cols(matrix1));
    compute_column_variances(matrix1, new_matrix);

  } else if (strcmp(operation, "scalarmult") == 0) {
    
    /* Multiply the matrix by a scalar value. */
    new_matrix = allocate_matrix(get_num_rows(matrix1), get_num_cols(matrix1));
    copy_matrix(matrix1, new_matrix);
    scalar_mult_matrix(value, new_matrix);

  } else if (strcmp(operation, "dotmult") == 0) {
    
    /* Multiply corresponding matrix values. */
    new_matrix = allocate_matrix(get_num_rows(matrix1), get_num_cols(matrix1));
    copy_matrix(matrix1, new_matrix);
    mult_matrix(matrix2, new_matrix);

  } else if (strcmp(operation, "posdef") == 0) {

    /* Make the matrix symmetric. */
    new_matrix = average_across_diagonal(matrix1, matrix1);
    
    /* Make the matrix positive definite. */
    make_positive_definite(new_matrix);

  } else if (strcmp(operation, "row-correlation") == 0) {

    // Deal with optional second argument.
    if (matrix2 == NULL) {
      matrix2 = matrix1;
      rdb_matrix2 = rdb_matrix1;
    }

    /* Compute all pairwise correlations. */
    new_matrix = compute_pairwise_correlations(TRUE, matrix1, matrix2);

  } else if (strcmp(operation, "col-correlation") == 0) {

    // Deal with optional second argument.
    if (matrix2 == NULL) {
      matrix2 = matrix1;
      rdb_matrix2 = rdb_matrix1;
    }

    /* Compute all pairwise correlations. */
    new_matrix = compute_pairwise_correlations(FALSE, matrix1, matrix2);

  } else if (strcmp(operation, "jacard") == 0) {

    // Compute the Jacard similarity between all pairs of rows.
    new_matrix = compute_pairwise_jacard(matrix1);

  } else if (strcmp(operation, "shuffle") == 0) {

    /* Randomly shuffle the matrix. */
    my_srand(seed);
    new_matrix = allocate_matrix(get_num_rows(matrix1), get_num_cols(matrix1));
    copy_matrix(matrix1, new_matrix);
    shuffle_matrix(new_matrix);

  } else if (strcmp(operation, "alignment") == 0) {

    // Compute the alignment between two matrices.
    printf("%*.*g\n", width, precision, alignment(matrix1, matrix2));
    free_matrix(matrix1);
    free_matrix(matrix2);
    return(0);

  } else if (strcmp(operation, "diagonal") == 0) {

    // Create a diagonal matrix.
    matrix2 = make_diagonal_matrix(get_num_rows(matrix1));

    // Align the given matrix to the diagonal one.
    printf("%*.*g\n", width, precision, alignment(matrix1, matrix2));
    free_matrix(matrix1);
    free_matrix(matrix2);
    return(0);

  } else if (strcmp(operation, "label-align") == 0) {
    
    // Create a YY' matrix from the labels.
    new_matrix = create_label_matrix(matrix2);

    // Align the given matrix to the label matrix.
    printf("%*.*g\n", width, precision, alignment(matrix1, new_matrix));
    free_matrix(new_matrix);
    free_matrix(matrix1);
    free_matrix(matrix2);
    return(0);

  } else if (strcmp(operation, "scalediag") == 0) {
    
    // Scale the diagonal of the given matrix.
    new_matrix = allocate_matrix(get_num_rows(matrix1), get_num_cols(matrix1));
    copy_matrix(matrix1, new_matrix);
    replace_diagonal(new_matrix);

  } else if (strcmp(operation, "euclidean") == 0) {
    
    // Make sure the matrix is square.
    new_matrix = euclidean(matrix1);

  } else if (strcmp(operation, "center") == 0) {

    // Make a copy and then center it.
    new_matrix = allocate_matrix(get_num_rows(matrix2), get_num_cols(matrix2));
    copy_matrix(matrix2, new_matrix);
    center_matrix(matrix1, new_matrix);

  } else if (strcmp(operation, "setdiag") == 0) {

    // Make a copy and then set the diagonal value.
    new_matrix = allocate_matrix(get_num_rows(matrix1), get_num_cols(matrix1));
    copy_matrix(matrix1, new_matrix);
    set_diagonal_matrix(value, new_matrix);

  } else if (strcmp(operation, "diffusion") == 0) {

    // Make a copy and then perform diffusion.
    new_matrix = allocate_matrix(get_num_rows(matrix1), get_num_cols(matrix1));
    copy_matrix(matrix1, new_matrix);
    perform_diffusion(value, new_matrix);

  } else if (strcmp(operation, "euclidean") == 0) {

    // Convert a kernel matrix to Euclidean distances.
    new_matrix = kernel_to_distance(matrix1);

  } else if (strcmp(operation, "int-dim") == 0) {

    // Compute the intrinsic dimension.
    new_matrix = kernel_to_distance(matrix1);
    value = get_intrinsic_dimensionality(new_matrix);

  } else if (strcmp(operation, "simple-diffusion") == 0) {

    // Compute a diffusion on the network.
    new_matrix = diffuse_matrix(value, DIFFUSION_ITERATIONS, matrix1);

  } else if (strcmp(operation, "modav") == 0) {

    // Compute the mean off-diagonal value.
    new_matrix = matrix1; // Avoid seg fault later.
    value = mean_off_diagonal_absolute_value(matrix1);

  } else {
    die("Unknown operation (%s).", operation);
  }


  /****** CREATE RDB VERSION OF THE NEW MATRIX. *****/

  if (rdb_format) {

    /* Allocate the RDB matrix. */
    rdb_new_matrix = allocate_rdb_matrix(get_num_rows(new_matrix),
					 get_num_cols(new_matrix), new_matrix);

    /* Transfer the corner string. */
    set_corner_string(get_corner_string(rdb_matrix1), rdb_new_matrix);

    /* Copy the row and column names. */
    if ((operation == NULL) ||
	(strcmp(operation, "none") == 0) ||
	(strcmp(operation, "randomize") == 0) ||
	(strcmp(operation, "square") == 0) ||
	(strcmp(operation, "invert") == 0) ||
	(strcmp(operation, "binarize") == 0) ||
	(strcmp(operation, "symdiag") == 0) ||
	(strcmp(operation, "adddiag") == 0) ||
	(strcmp(operation, "add") == 0) ||
	(strcmp(operation, "rescalerow") == 0) ||
	(strcmp(operation, "rescalecol") == 0) ||
	(strcmp(operation, "normalize") == 0) ||
	(strcmp(operation, "zeromeanrow") == 0) ||
	(strcmp(operation, "zeromeancol") == 0) ||
	(strcmp(operation, "varone") == 0) ||
	(strcmp(operation, "scalarmult") == 0) ||
	(strcmp(operation, "dotmult") == 0) ||
	(strcmp(operation, "posdef") == 0) ||
	(strcmp(operation, "shuffle") == 0) ||
	(strcmp(operation, "scalediag") == 0) ||
	(strcmp(operation, "euclidean") == 0) ||
	(strcmp(operation, "setdiag") == 0) ||
	(strcmp(operation, "diffusion") == 0) ||
	(strcmp(operation, "euclidean") == 0) || 
	(strcmp(operation, "simple-diffusion") == 0)) {
      
      set_row_names(get_row_names(rdb_matrix1), rdb_new_matrix);
      set_col_names(get_col_names(rdb_matrix1), rdb_new_matrix);

    } else if ((strcmp(operation, "normalizerow") == 0) ||
	       (strcmp(operation, "normalizecol") == 0) ||
	       (strcmp(operation, "center") == 0)) {
	
      set_row_names(get_row_names(rdb_matrix2), rdb_new_matrix);
      set_col_names(get_col_names(rdb_matrix2), rdb_new_matrix);

    } else if ((strcmp(operation, "covariance") == 0) ||
	       (strcmp(operation, "jacard") == 0)) {

      set_row_names(get_row_names(rdb_matrix1), rdb_new_matrix);
      set_col_names(get_row_names(rdb_matrix1), rdb_new_matrix);

    } else if (strcmp(operation, "row-correlation") == 0) {

      set_row_names(get_row_names(rdb_matrix1), rdb_new_matrix);
      set_col_names(get_row_names(rdb_matrix2), rdb_new_matrix);

    } else if (strcmp(operation, "col-correlation") == 0) {

      set_row_names(get_col_names(rdb_matrix1), rdb_new_matrix);
      set_col_names(get_col_names(rdb_matrix2), rdb_new_matrix);

    } else if (strcmp(operation, "transpose") == 0) {

      set_row_names(get_col_names(rdb_matrix1), rdb_new_matrix);
      set_col_names(get_row_names(rdb_matrix1), rdb_new_matrix);

    } else if (strcmp(operation, "eigenvectors") == 0) {

      set_row_names(get_row_names(rdb_matrix1), rdb_new_matrix);
      set_col_names(get_row_names(rdb_matrix1), rdb_new_matrix);

    } else if (strcmp(operation, "eigenvalues") == 0) {

      add_row_name("eigenvalues", rdb_new_matrix);
      set_col_names(get_row_names(rdb_matrix1), rdb_new_matrix);

    } else if (strcmp(operation, "eigenvector1") == 0) {

      add_col_name("eigenvector1", rdb_new_matrix);
      set_row_names(get_row_names(rdb_matrix1), rdb_new_matrix);

    } else if (strcmp(operation, "getcolsums") == 0) {

      set_row_names(get_col_names(rdb_matrix1), rdb_new_matrix);
      add_col_name("colsums", rdb_new_matrix);

    } else if (strcmp(operation, "getrowsums") == 0) {

      set_row_names(get_row_names(rdb_matrix1), rdb_new_matrix);
      add_col_name("rowsums", rdb_new_matrix);

    } else if (strcmp(operation, "getdiag") == 0) {

      set_row_names(get_row_names(rdb_matrix1), rdb_new_matrix);
      add_col_name("diagonal", rdb_new_matrix);

    } else if (strcmp(operation, "colvars") == 0) {

      add_row_name("variance", rdb_new_matrix);
      set_col_names(get_col_names(rdb_matrix1), rdb_new_matrix);

    }
  }

  /***** PRINT OUTPUT *****/

  // Special case: a single value.
  if ((strcmp(operation, "int-dim") == 0) ||
      (strcmp(operation, "modav") == 0)) {
    printf("%g\n", value);
  }

  // Print the resulting matrix.
  else if (rdb_format) {
    print_rdb_matrix(rdb_new_matrix, format_line, width, precision, stdout);
    if (rdb_matrix1 != rdb_matrix2) {
      free_rdb_matrix(rdb_matrix2);
    }
    free_rdb_matrix(rdb_matrix1);
    free_rdb_matrix(rdb_new_matrix);
  } else {
    print_matrix(new_matrix, width, precision, FALSE, stdout);
    free_matrix(matrix1);
    free_matrix(matrix2);
    free_matrix(new_matrix);
  }
  
  return(0);
}
#endif //  RDBMATRIX_MAIN

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
