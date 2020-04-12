/*****************************************************************************
 * FILE: class-array.h
 * AUTHOR: William Stafford Noble
 * CREATE DATE: 2/11/99
 * PROJECT: SVM
 * COPYRIGHT: 1999-2001, Columbia University (see ../doc/copyright.html)
 * VERSION: $Revision: 1.1 $
 * DESCRIPTION: Data class for an array of Boolean classifications.
 *****************************************************************************/
#ifndef CLASS_ARRAY_H
#define CLASS_ARRAY_H

#include "utils.h"
#include "array.h"
#include <stdio.h>

/*****************************************************************************
 * Type definition.
 *****************************************************************************/
typedef struct class_array_t CLASS_ARRAY_T;

/*****************************************************************************
 * Allocate an empty class array
 *****************************************************************************/
CLASS_ARRAY_T* allocate_class_array
  (int num_items);

/*****************************************************************************
 * Find out how many items the array contains.
 *****************************************************************************/
int get_num_items
  (CLASS_ARRAY_T* classes);

/*****************************************************************************
 * Get the classification of one element in the matrix.  1 is TRUE; -1
 * is FALSE.
 *****************************************************************************/
BOOLEAN_T get_class
  (int             item,
   CLASS_ARRAY_T * classes);

double get_class_sign
  (int             item,
   CLASS_ARRAY_T * classes);


/*****************************************************************************
 * Set the classification of one element in the matrix.
 *****************************************************************************/
void set_class
  (int             item,
   BOOLEAN_T       value,
   CLASS_ARRAY_T * classes);

/*****************************************************************************
 * Get the number of members in one class.
 *****************************************************************************/
int get_class_size
  (BOOLEAN_T       class,
   CLASS_ARRAY_T * classes);

/*****************************************************************************
 * Return the classifications as an array of doubles.
 *****************************************************************************/
ARRAY_T* get_class_array
  (CLASS_ARRAY_T * classes);

/*****************************************************************************
 * Read classifications from an external file.
 *****************************************************************************/
CLASS_ARRAY_T* read_classifications
  (BOOLEAN_T format_line,
   FILE*     class_file,
   int  class_number);

/*****************************************************************************
 * Free dynamic memory used by a class array.
 *****************************************************************************/
void free_class_array
  (CLASS_ARRAY_T* classes);

#endif


/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
