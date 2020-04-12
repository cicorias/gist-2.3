/*****************************************************************************
 * FILE: compute-weights.c
 * AUTHOR: William Stafford Noble
 * CREATE DATE: 2/11/99
 * PROJECT: SVM
 * COPYRIGHT: 1999-2001, Columbia University (see ../doc/copyright.html)
 * VERSION: $Revision: 1.1 $
 * DESCRIPTION: Given a matrix of pairwise distances for positive and
 * negative examples, compute weights that maximize a quadratic
 * objective function.
 *
 * This program is based upon
 *
 *   "A discriminative framework for detecting remote protein
 *   homologies."  T. Jaakkola, M. Diekhans, D. Haussler.  1998.
 *
 *****************************************************************************/
#include "compute-weights.h"
#include "compute-kernel.h"
#include "class-array.h"
#include "linear-algebra.h"
#include "matrix.h"
#include "array.h"
#include "utils.h"
#include <time.h>
#include <string.h>
#include <stdlib.h>

/*****************************************************************************
 * Shuffle an array of integers using Knuth's algorithm.
 *****************************************************************************/
static void shuffle
  (int   num_items,
   int*  data)
{
  int i_item;
  int i_rand;
  int temp_swap;

  // Shuffle the array.
  for (i_item = 0; i_item < num_items; i_item++) {  

    // Select a random position to the right of the current position.
    i_rand = (int)(my_drand() * (double)(num_items - i_item)) + i_item; 

    // Swap 'em.
    temp_swap = data[i_item];
    data[i_item] = data[i_rand];
    data[i_rand] = temp_swap;
  }

}


/*****************************************************************************
 * Randomly return an index between 0 and a given value.
 *
 * On the first call, current_item should be -1.
 *
 * The function allocates the data array and randomly fill it with
 * ascending integers.  On subsequent calls, iterate through the
 * array, returning the next index.  Once we iterate through the whole
 * thing, re-shuffle and start over.
 *****************************************************************************/
int select_random_item
  (int   num_items,
   int*  current_item,
   int** data)
{
  int i_item;
  int return_value;
  
  if (*current_item == -1) {
    *current_item = 0;

    /* Allocate memory for the array. */
    myfree(*data);
    *data = (int*)mymalloc(sizeof(int) * num_items);
    
    /* Fill the array with ascending integers. */
    for (i_item = 0; i_item < num_items; i_item++) {
      (*data)[i_item] = i_item;
    }
  }

  /* Shuffle if we're at the beginning of the array. */
  if (*current_item == 0) {
    shuffle(num_items, *data);
  }

  /* Store the return value. */
  return_value = (*data)[*current_item];

  /* Move to the next item. */
  (*current_item)++;
  if (*current_item >= num_items) {
    *current_item = 0;
  }

  /* Return. */
  return(return_value);
}

/*****************************************************************************
 * The discriminant function is used to determine whether a given
 * example is classified positively or negatively.
 *
 * This function implements equation (4) from the paper cited above.
 *****************************************************************************/
static double compute_discriminant
  (int            this_item,
   ARRAY_T*       weights,
   CLASS_ARRAY_T* classes,
   MATRIX_T*      kernel_matrix)
{
  double return_value;
  double this_value;
  int i_item;
  int num_items = get_num_items(classes);

  return_value = 0.0;
  for (i_item = 0; i_item < num_items; i_item++) {

    /* Weight the distance appropriately. */
    this_value = get_array_item(i_item, weights) *
      get_matrix_cell(this_item, i_item, kernel_matrix);

    /* Add or subtract, depending upon whether this is a positive or
       negative example. */
    return_value += this_value * get_class_sign(i_item, classes);
  }
  return(return_value);
}

/*****************************************************************************
 * Compute the objective function, equation (7) (but divide by the
 * number of items).
 *****************************************************************************/
double compute_objective
  (ARRAY_T*       weights,
   CLASS_ARRAY_T* classes,
   MATRIX_T*      kernel_matrix)
{
  int    num_items;
  int    i_item;
  double this_discriminant;
  double this_weight;
  double sum;

  sum = 0.0;
  num_items = get_num_items(classes);
  for (i_item = 0; i_item < num_items; i_item++) {
    this_weight = get_array_item(i_item, weights);
    this_discriminant = compute_discriminant(i_item, weights, classes, 
					     kernel_matrix);
    
    this_discriminant *= get_class_sign(i_item, classes);

    sum += this_weight * (2.0 - this_discriminant);
  }
  return(sum / num_items);
}

/*****************************************************************************
 * Keep a local copy of the weights array and signal if they've
 * stopped changing.
 *
 * Convergence is reached when the delta is below the convergence
 * threshold.
 *****************************************************************************/
BOOLEAN_T converged
  (BOOLEAN_T      reset,  /* Reset everything for a new optimization. */
   double         convergence_threshold,
   ARRAY_T*       weights,
   CLASS_ARRAY_T* classes,
   MATRIX_T*      kernel_matrix)
{
  double        objective;       /* Current value of the objective. */
  static double prev_objective = 0.0;  /* Prev objective for computing delta. */
  double        delta = 0.0;     /* Change in objective. */
  static int    iteration = 0;   /* Total number of iterations. */
  static BOOLEAN_T warned_once = FALSE;

  if (reset) {
    iteration = 0;
    delta = 0.0;
    prev_objective = 0.0;
    return(FALSE);
  }

  /* Compute the new objective. */
  objective = compute_objective(weights, classes, kernel_matrix);

  /* Compute the change in objective. */
  delta = objective - prev_objective;

  /* Store this objective for next time. */
  prev_objective = objective;

  /* Increment the iteration counter. */
  iteration++;

  /* Tell the user what's up. */
  if (verbosity >= DUMP_VERBOSE) {
    fprintf(stderr, "Weights: ");
    print_array(weights, 6, 4, TRUE, stderr);
  }
  if (verbosity > NORMAL_VERBOSE) {
    fprintf(stderr, "Iteration %d: objective=%g", iteration,
	    objective);
    if (iteration > 1) {
      fprintf(stderr, " delta=%g", delta);
    }
    fprintf(stderr, "\n");
  }

  /* Complain if the delta is negative. */
  if ((delta < 0.0) && 
      (fabs(delta) > convergence_threshold) && 
      (iteration != 1) && 
      (warned_once == FALSE)) {
    fprintf(stderr, "Negative delta.\n");
    warned_once = TRUE;
  }

  return(fabs(delta) < convergence_threshold);
}
      
/*****************************************************************************
 * Update one item's weight.  This update rule maximizes the
 * constrained maximization of J(\lambda).  This function implements
 * equations (9) and (10) in Jaakkola et al.
 *****************************************************************************/
static double update_weight
  (BOOLEAN_T      constrain_weights, /* Prevent weights from exceeding 1? */
   double         constraint,
   int            this_item,         /* Index of the current item. */
   ARRAY_T*       weights,           /* The weights being updated. */
   CLASS_ARRAY_T* classes,           /* Classifications of training set. */
   MATRIX_T*      kernel_matrix)      /* Matrix of kernel distances. */
{
  double this_discriminant;
  double self_distance;
  double this_weight;
  double new_weight;
  double class;
  
  this_discriminant = compute_discriminant(this_item, weights, classes,
					   kernel_matrix);
  self_distance = get_matrix_cell(this_item, this_item, kernel_matrix);
  this_weight = get_array_item(this_item, weights);

  /* Weight negative examples oppositely. */
  class = get_class_sign(this_item, classes);

  /* This is equation (8). */
  new_weight = 1.0 - (class * this_discriminant) 
    + (this_weight * self_distance);

  /* Divide by k(x,x), checking for divide-by-zero. */
  if (self_distance == 0.0) {
    new_weight /= new_weight;
  } else {
    new_weight /= self_distance;
  }

  if (verbosity >= DUMP_VERBOSE) {
    fprintf(stderr, "Class=%d Discriminant=%g\n", (int)class, 
	    this_discriminant);
  }

#ifdef DEBUG
  /* If we update the weights without the constraints, the
     discriminant of this item should be 1.0. */
  if (self_distance != 0.0) {
    this_weight = get_array_item(this_item, weights);
    set_array_item(this_item, new_weight, weights);
    this_discriminant = compute_discriminant(this_item, weights, classes,
					     kernel_matrix);
    if (!almost_equal(this_discriminant, class, 0.0001)) {
      die("discriminant(%d)=%g, should be %g\n", this_item, this_discriminant, class);
    }
    set_array_item(this_item, this_weight, weights);
  }
#endif

  /* Constrain the weight. */
  if ((constrain_weights) && (new_weight > constraint)) {
    new_weight = constraint;
  } else if (new_weight < 0.0) {
    new_weight = 0.0;
  }

  return(new_weight);
}

/**************************************************************************
 * From a given set of class labels, kernel matrix and weights,
 * produce corresponding smaller versions of each that contain only
 * the support vectors.  If the first (Boolean) parameter, is false,
 * do the reverse operation.
 **************************************************************************/
#define SLOP 0.0000001 
static void extract_support_vectors
  (BOOLEAN_T       extract,
   CLASS_ARRAY_T*  classes,
   MATRIX_T*       kernel_matrix,
   ARRAY_T*        weights,
   CLASS_ARRAY_T** sv_classes,
   MATRIX_T**      sv_kernel_matrix,
   ARRAY_T**       sv_weights)
{
  static int* sv_list;  /* Indices of the support vectors. */
  static int num_svs;
  int i_item;
  int num_items = get_num_items(classes);
  double weight;
  int i_sv_row;
  int i_sv_col;
  int i_row;
  int i_col;
  BOOLEAN_T class;
  double kernel_value;
  
  if (extract) {

    /* Re-set the number of support vectors. */
    num_svs = 0;

    /* Make a list of the support vector indices. */
    for (i_item = 0; i_item < num_items; i_item++) {
      weight = get_array_item(i_item, weights);
      if (!almost_equal(weight, 0.0, SLOP)) {
	sv_list = (int*)realloc(sv_list, (num_svs + 1) * sizeof(int));
	sv_list[num_svs] = i_item;
	num_svs++;
      }
    }

    /* Allocate space for the new data structures. */
    *sv_classes = allocate_class_array(num_svs);
    *sv_kernel_matrix = allocate_matrix(num_svs, num_svs);
    *sv_weights = allocate_array(num_svs);
  }

  /* Consider each support vector. */
  for (i_sv_row = 0; i_sv_row < num_svs; i_sv_row++) {
    i_row = sv_list[i_sv_row];

    /* Copy the weight. */
    if (extract) {
      weight = get_array_item(i_row, weights);
      set_array_item(i_sv_row, weight, *sv_weights);
    } else {
      weight = get_array_item(i_sv_row, *sv_weights);
      set_array_item(i_row, weight, weights);
    }

    /* Copy the class. */
    if (extract) {
      class = get_class(i_row, classes);
      set_class(i_sv_row, class, *sv_classes);
    }

    /* Copy the matrix row. */
    for (i_sv_col = 0; i_sv_col < num_svs; i_sv_col++) {
      i_col = sv_list[i_sv_col];

      if (extract) {
	kernel_value = get_matrix_cell(i_row, i_col, kernel_matrix);
	set_matrix_cell(i_sv_row, i_sv_col, kernel_value, *sv_kernel_matrix);
      } else {
	kernel_value = get_matrix_cell(i_sv_row, i_sv_col, *sv_kernel_matrix);
	set_matrix_cell(i_row, i_col, kernel_value, kernel_matrix);
      }
      
    }
  }
}

/*****************************************************************************
 * Optimize the weights so that the discriminant function puts the
 * negatives close to -1 and the positives close to +1.
 *****************************************************************************/

#define ITEREXIT 2
#define TIMEOUT 3
static long int optimize_weights
  (double         convergence_threshold,
   long int       maxiter,
   double         maxtime,
   double         positive_constraint,
   double         negative_constraint,
   CLASS_ARRAY_T* classes,
   MATRIX_T*      kernel_matrix,
   ARRAY_T*       weights)
{
  int      i_item;
  int      rand_item;
  int      num_items = get_num_items(classes);
  long int iter = 0;
  double   start_time;
  double   elapsed_time;
  double   new_weight;
  double   constraint;
  BOOLEAN_T   constrain_weights;
  static int  current_item; // Used by select_random_item.
  static int* data;         // Used by select_random_item.

  // Record the starting time.
  start_time = myclock() / 1.0E6;

  // Constrain weights if both are non-zero.
  constrain_weights = ((positive_constraint != 0.0) ||
		       (negative_constraint != 0.0));

  /* Iteratively improve the weights until convergence. */
  current_item = -1;
  while (!converged(FALSE, convergence_threshold, weights, classes,
		    kernel_matrix)) {

    for (i_item = 0; i_item < num_items; i_item++) {

      /* Randomly select a weight to update. */
      rand_item = select_random_item(num_items, &current_item, &data);

      /* Set the constraint, based upon the class of this item. */
      if (get_class(rand_item, classes)) {
	constraint = positive_constraint;
      } else {
	constraint = negative_constraint;
      }

      /* Calculate the new weight. */
      new_weight = update_weight(constrain_weights, constraint, rand_item,
				 weights, classes, kernel_matrix);

      /* Update the weight array accordingly. */
      if (verbosity >= DUMP_VERBOSE) {
	if (get_array_item(rand_item, weights) != new_weight) {
	  fprintf(stderr, "Changing weight %d from %g to %g.\n", rand_item,
		  get_array_item(rand_item, weights), new_weight);
	}
      }
      set_array_item(rand_item, new_weight, weights);
    }
    iter++;

    /* Check to be sure we haven't reached the iteration limit. */
    if ((maxiter != 0) &&  (iter >= maxiter)) {
      fprintf(stderr, "Warning: Terminating after failure to converge ");
      fprintf(stderr, "in %ld iterations.\n", maxiter);
      exit(ITEREXIT);
    }

    // Check to be sure we haven't reached the time limit.
    elapsed_time = (myclock() / 1.0E6) - start_time;
    if ((maxtime != 0.0) && (elapsed_time >= maxtime)) {
      fprintf(stderr, "Warning: Terminating after failure to converge ");
      fprintf(stderr, "in %g seconds.\n", maxtime);
      exit(TIMEOUT);
    }
  }
  return(iter + 1);
}

/*****************************************************************************
 * Optimize the weights on the entire training set, then extract the
 * support vectors and optimize them some more.
 *****************************************************************************/
long int double_optimize_weights
  (double         convergence_threshold,
   long int       maxiter,
   double         maxtime,
   double         positive_constraint,
   double         negative_constraint,
   CLASS_ARRAY_T* classes,
   MATRIX_T*      kernel_matrix,
   ARRAY_T*       weights)
{
  long int       iter;
  CLASS_ARRAY_T* sv_classes;
  MATRIX_T*      sv_kernel_matrix;
  ARRAY_T*       sv_weights;

  /* Make sure the dimensions all match. */
  myassert(1, get_num_items(classes) == get_num_rows(kernel_matrix), "Number of class labels (%d) does not equal number of rows (%d) in kernel matrix", get_num_items(classes), get_num_rows(kernel_matrix));
  myassert(1, get_num_items(classes) == get_num_cols(kernel_matrix), "Number of class labels (%d) does not equal number of columns (%d) in kernel matrix", get_num_items(classes), get_num_cols(kernel_matrix));
  myassert(1, get_num_items(classes) == get_array_length(weights),  "Number of class labels (%d) does not equal number of weights values (%d)", get_num_items(classes), get_array_length(weights));

  /* Optimize the weights. */
  if (verbosity > LOW_VERBOSE) {
    fprintf(stderr, "Optimizing weights on %d training examples.\n",
	    get_array_length(weights));
  }
  iter = optimize_weights(convergence_threshold, 
			  maxiter,
			  maxtime,
			  positive_constraint,
			  negative_constraint,
			  classes,
			  kernel_matrix,
			  weights);

  /* Pull out the support vectors. */
  extract_support_vectors(TRUE, classes, kernel_matrix, weights,
			  &sv_classes, &sv_kernel_matrix, &sv_weights);

  /* Optimize only the support vectors. */
  if ((maxiter == 0 || iter < maxiter) && (get_array_length(sv_weights) > 0)) {
    if (verbosity > LOW_VERBOSE) {
      fprintf(stderr, "Optimizing weights on %d support vectors.\n",
	      get_array_length(sv_weights));
    }
    iter += optimize_weights(convergence_threshold, 
			     maxiter,
			     maxtime,
			     positive_constraint,
			     negative_constraint,
			     sv_classes,
			     sv_kernel_matrix,
			     sv_weights);
  }

  /* Put the support vectors back in. */
  extract_support_vectors(FALSE, classes, kernel_matrix, weights,
			  &sv_classes, &sv_kernel_matrix, &sv_weights);

  /* Free local dynamic memory. */
  free_class_array(sv_classes);
  free_matrix(sv_kernel_matrix);
  free_array(sv_weights);

  return(iter);
}

/****************************************************************************
 * Encode the classifications in the weights by multiplying the negative
 * examples by -1.
 ****************************************************************************/
void sign_weights
  (CLASS_ARRAY_T* classes,
   ARRAY_T*       weights)
{
  int     i_item;
  double  item;
  int     num_items = get_num_items(classes);

  for (i_item = 0; i_item < num_items; i_item++) {
    item = get_array_item(i_item, weights);
    item *= get_class_sign(i_item, classes);
    set_array_item(i_item, item, weights);
  }
}




/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
