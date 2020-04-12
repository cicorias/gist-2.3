/*****************************************************************************
 * FILE: class-array.c
 * AUTHOR: William Stafford Noble
 * CREATE DATE: 2/11/99
 * PROJECT: SVM
 * COPYRIGHT: 1999-2001, Columbia University (see ../doc/copyright.html)
 * VERSION: $Revision: 1.1 $
 * DESCRIPTION: Data class for an array of Boolean classifications.
 *****************************************************************************/
#include "matrix.h"
#include "class-array.h"
#include "rdb-matrix.h"
#include "array.h"
#include "utils.h"
#include <assert.h>
#include <stdio.h>

/*****************************************************************************
 * Type definition.
 *****************************************************************************/
struct class_array_t {
  int           num_items;
  int           num_positive;
  BOOLEAN_T*    classifications;
};

/*****************************************************************************
 * Allocate an empty classification array.
 *****************************************************************************/
CLASS_ARRAY_T* allocate_class_array
  (int num_items)
{
  CLASS_ARRAY_T* new_classes;

  new_classes = (CLASS_ARRAY_T*)mymalloc(sizeof(CLASS_ARRAY_T));
  new_classes->num_items = num_items;
  new_classes->num_positive = 0;
  new_classes->classifications 
    = (BOOLEAN_T*)mycalloc(num_items, sizeof(BOOLEAN_T));
  return(new_classes);
}

/*****************************************************************************
 * Find out how many items the array contains.
 *****************************************************************************/
int get_num_items
  (CLASS_ARRAY_T* classes)
{
  if (classes == NULL) {
    die("Attempted to access null class array.\n");
  }

  return(classes->num_items);
}

/*****************************************************************************
 * Get the classification of one element.  1 is TRUE; -1 is FALSE.
 *****************************************************************************/
BOOLEAN_T get_class
  (int             item,
   CLASS_ARRAY_T * classes)
{
  if (classes == NULL) {
    die("Attempted to access null class array.\n");
  }

  assert(item >= 0);
  assert(item < classes->num_items);

  return(classes->classifications[item]);
}

double get_class_sign
  (int             item,
   CLASS_ARRAY_T * classes)
{
  if (get_class(item, classes)) {
    return(1.0);
  } else {
    return(-1.0);
  }
}

/*****************************************************************************
 * Set the classification of one element in the matrix.
 *****************************************************************************/
void set_class
  (int             item,
   BOOLEAN_T       value,
   CLASS_ARRAY_T * classes)
{
  if (classes == NULL) {
    die("Attempted to access null class array.\n");
  }

  assert(item >= 0);
  assert(item < classes->num_items);

  /* Keep track of the class size. */
  if (classes->classifications[item] != value) {
    if (value) {
      classes->num_positive++;
    } else {
      classes->num_positive--;
    }
  }

  classes->classifications[item] = value;
}

/*****************************************************************************
 * Get the number of members in one class.
 *****************************************************************************/
int get_class_size
  (BOOLEAN_T       class,
   CLASS_ARRAY_T * classes)
{
  if (classes == NULL) {
    die("Attempted to access null class array.\n");
  }

  if (class) {
    return(classes->num_positive);
  } /* else */
  return(classes->num_items - classes->num_positive);
}

/*****************************************************************************
 * Return the classifications as an array of doubles.
 *****************************************************************************/
ARRAY_T* get_class_array
  (CLASS_ARRAY_T * classes)
{
  int num_items;
  int i_item;
  ARRAY_T* return_value;

  /* Allocate the output array. */
  num_items = get_num_items(classes);
  return_value = allocate_array(num_items);

  /* Fill the array with 1s and -1s. */
  for (i_item = 0; i_item < num_items; i_item++) {
    set_array_item(i_item, get_class_sign(i_item, classes), return_value);
  }

  /* Return the allocated array. */
  return(return_value);
}


/*****************************************************************************
 * Read classifications from an external file.
 *****************************************************************************/
CLASS_ARRAY_T* read_classifications
  (BOOLEAN_T format_line,
   FILE*     class_file,
   int  class_number)
{
  RDB_MATRIX_T* rdb_classes;
  MATRIX_T*     raw_classes;
  int num_items;
  int i_item;
  int class_index;
  double classification;
  CLASS_ARRAY_T* return_value;

  /* Read the class info into a matrix. */
  rdb_classes = read_rdb_matrix(format_line, NULL, class_file);
  raw_classes = get_raw_matrix(rdb_classes);
  num_items = get_num_rows(raw_classes);
  myassert(1, num_items > 0, "No items in class file");

  // if we are using a 'multiclass' file
  myassert(1, get_num_cols(raw_classes) >= class_number, 
	   "Invalid classnumber %d; there is/are only %d class(es) available in the class file.", 
	   class_number, get_num_cols(raw_classes));
  class_index = class_number - 1;
  if (verbosity > NORMAL_VERBOSE) 
    fprintf(stderr, "Using class number %d\n", class_number);

  /* Allocate the class. */
  return_value = allocate_class_array(num_items);

  /* Store the classifications. */
  return_value->num_positive = 0;

  for (i_item = 0; i_item < num_items; i_item++) {
    
    classification = get_matrix_cell(i_item, class_index, raw_classes);

    if (classification == 1.0) {
      return_value->classifications[i_item] = TRUE;
      return_value->num_positive++;
    } else if (classification == -1.0 || classification == 0.0) {
      return_value->classifications[i_item] = FALSE;
    } else {
      die("Illegal classification %d (%d).\n", i_item, classification);
    }
  }

  free_rdb_matrix(rdb_classes);
  return(return_value);
}

/*****************************************************************************
 * Free dynamic memory used by a class array.
 *****************************************************************************/
void free_class_array
  (CLASS_ARRAY_T * classes)
{
  if (classes == NULL) {
    return;
  }
  myfree(classes->classifications);
  myfree(classes);
}


/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
