/*****************************************************************************
 * FILE: kfd.h
 * AUTHOR: William Stafford Noble
 * CREATE DATE: 2/28/2000
 * PROJECT: SVM
 * COPYRIGHT: 1999-2001, Columbia University (see ../doc/copyright.html)
 * VERSION: $Revision: 1.1 $
 * DESCRIPTION: 
 *****************************************************************************/
#ifndef KFD_H
#define KFD_H
#include "matrix.h"
#include "class-array.h"

/*****************************************************************************
 * Find the kernel Fisher discriminant.
 *****************************************************************************/
MATRIX_T* find_kfd
  (double         regularizer,
   CLASS_ARRAY_T* classes,
   MATRIX_T*      kernel_matrix);

#endif
/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
