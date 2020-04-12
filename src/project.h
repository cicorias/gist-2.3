/*****************************************************************************
 * FILE: project.h
 * AUTHOR: William Stafford Noble
 * CREATE DATE: 1/26/2000
 * PROJECT: SVM
 * COPYRIGHT: 1999-2001, Columbia University (see ../doc/copyright.html)
 * VERSION: $Revision: 1.1 $
 * DESCRIPTION: Project a given set of data onto a set of eigenvectors.
 *****************************************************************************/
#ifndef PROJECT_H
#define PROJECT_H

#include "matrix.h"

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
   MATRIX_T* output_matrix);

#endif

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
