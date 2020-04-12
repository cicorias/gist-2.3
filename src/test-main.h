/****************************************************************************
 * FILE: test-main.h
 * AUTHOR: William Stafford Noble
 * CREATE DATE: 3/29/2005
 * PROJECT: SVM
 * COPYRIGHT: 2005, University of Washington
 * VERSION: $Revision: 1.2 $
 * DESCRIPTION: Main procedure for classifying or projecting.
 ****************************************************************************/
#ifndef TEST_MAIN_H
#define TEST_MAIN_H
#endif

#include "kernel.h"
#include "utils.h"

/****************************************************************************
 * Read the kernel name and polynomial parameters from the header of a
 * given weights file.
 ****************************************************************************/
void read_kernel_parameters
  (char*      learned_filename,
   double*    bias,
   double*    sum_of_weights,
   KERNEL_T*  kernel);

/****************************************************************************
 * Put an RDB header on the output file.
 ****************************************************************************/
void print_output_header
  (char*     learned_filename,
   KERNEL_T* kernel,
   BOOLEAN_T print_time,
   FILE*     outfile);

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
 
