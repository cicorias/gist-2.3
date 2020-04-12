/*****************************************************************************
 * FILE: fselect.c
 * AUTHOR: William Stafford Noble, modifications and additions by Paul Pavlidis
 * CREATE DATE: 10/17/01
 * PROJECT: SVM
 * COPYRIGHT: 2001, Columbia University
 * VERSION: $Revision: 1.1 $
 * DESCRIPTION: Feature selection by Fisher criterion score.
 *****************************************************************************/
#include "fselect.h"
#include "cmdline.h"
#include "utils.h"
#include "matrix.h"
#include "array.h"
#include "linear-algebra.h"
#include "rdb-matrix.h"
#include "class-array.h"
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <time.h>

#ifndef FSELECT_DEBUG
#define FSELECT_DEBUG 0
#endif

/*****************************************************************************
 * Used by betai: Evaluates continued fraction for incomplete beta
 * function by modified Lentz's method.
 *
 * From Numerical Recipes in C.
 *****************************************************************************/
#define MAXIT 200
#define EPS 3.0e-7
#define FPMIN 1.0e-30
static double betacf
  (double a,
   double b,
   double x)
{
  int m, m2;
  float aa, c, d, del, h, qab, qam, qap;

  // These q's will be used in factors that occur in the coefficients.
  qab = a + b;
  qap = a + 1.0;
  qam = a - 1.0;

  // First step of Lentz's method.
  c = 1.0;
  d = 1.0 - (qab * (x / qap));

  if (fabs(d) < FPMIN) d = FPMIN;
  d = 1.0 / d;
  h = d;

  for (m = 1; m <= MAXIT; m++) {
    m2 = 2*m;
    aa= m * (b - m) * x / ((qam + m2) * (a + m2));

    // One step (the even one) of the recurrence.
    d = 1.0 + aa * d;
    if (fabs(d) < FPMIN) d = FPMIN;
    c = 1.0 + aa / c;
    if (fabs(c) < FPMIN) c = FPMIN;
    d = 1.0 / d;
    h *= d * c;
    aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2));

    // Next step of the recurrence (the odd one).
    d = 1.0 + aa * d;
    if (fabs(d) < FPMIN) d = FPMIN;
    c = 1.0 + aa / c;
    if (fabs(c) < FPMIN) c = FPMIN;
    d = 1.0 / d;
    del = d * c;
    h *= del;

    // Are we done?
    if (fabs(del - 1.0) < EPS) break;
  }
  if (m > MAXIT) {
    die("a (%g) or b (%g) too big, or MAXIT (%d) too small in betacf",
	a, b, MAXIT);
  }
  return(h);
}


/*****************************************************************************
 * Returns the incomplete beta function I_X(a,b).
 *
 * From Numerical Recipes in C.
 *****************************************************************************/
static double betai
  (double a,
   double b,
   double x)
{
  float bt;

  if (x < 0.0 || x > 1.0) {
    die("Bad x (%g) in routine betai.", x);
  }

  if (x == 0.0 || x == 1.0) {
    bt = 0.0;
  } else {
    bt = exp(lgamma(a + b) 
	     - lgamma(a)
	     - lgamma(b) 
	     + (a * log(x))
	     + (b * log(1.0 - x)));
  }

  if (x < (a + 1.0)/(a + b + 2.0)) {
    return(bt * betacf(a, b, x) / a);
  } // else
  return(1.0 - (bt * betacf(b, a, 1.0 - x) / b));
}


/***************************************************************************** 
 * ln(n!)
****************************************************************************/
static float lnfact
(int n) 
{
  if (n<0) {
    die("Attempt to measure lnfact of a negative value.");
    exit(1);
  }
  if (n<=1)
    return 0.0;
  else
    return lgamma(n+1.0);
}


/***************************************************************************** 
 * Binomial coefficient
****************************************************************************/
static double binomial_coeff
(int n, 
int k) 
{
  return floor(0.5+exp(lnfact(n)-lnfact(k)-lnfact(n-k)));
}


/***************************************************************************** 
 * Modified from Elisabetta Manduci's tpwy. I don't know where she
 * found, or if she is the originator of, this algorithm, it's not in
 * NRC afaik. It is used for the Mann-whitney u-test.
 *
 * Given a list a of k numbers a_0<a_1<...<a_(k-1) with a_i in {1,...,
 * n}, returns the list b: b_0<b_1<...<b_(k-1) which immediately
 * follows a in lexicographic order. It can be used to generate all
 * subsets of size k in {1,..., n}.
*****************************************************************************/
static void next_lex(int *ptr, int n, int k) {
  int h=k-1, s=n, count=0;
  
  while (ptr[h]==s) {
    count++;
    h--;
    s--;
  }
  h = k-1;
  s = n;
  while (ptr[h]==s) {
    ptr[h] = ptr[h-count]+count+1; 
    h--;
    s--;
    count--;
  }
  if (h>=0) 
    ptr[h] += 1;
}

/***************************************************************************
 * Calculate the U score for a feature. The data must already be rank
 * transformed. This could be simplified; this version does all the
 * checking for consistency. (Paul Pavlidis)
 ***************************************************************************/
static int calc_mwu_score
  (ARRAY_T* first_array,
   ARRAY_T* second_array)
{
  int num1 = get_array_length(first_array);
  int num2 = get_array_length(second_array);
  int R1, R2, u1, u2;
  R1 = array_total(first_array);
  R2 = array_total(second_array);

  // FIXME: This was off by one, so i commented it out. (WSN 1/19/02)
  // assert(R1+R2 == (num1+num2)*(num1+num2+1)/2);

  u1 =  num1*num2 + num1*(num1+1)/2 - R1;
  u2 =  num1*num2 + num2*(num2+1)/2 - R2;
  
  if (u1 > u2) {
    return u1;
  } else {
    return u2;
  }
}


/*****************************************************************************
 * Initialize the data structures to store the mann-whitney
 * distributions. 
 *
 * TODO: free this space at the end. Annoying because it is static and
 * freeing it at the right place will require either passing up a
 * pointer to it from compute_quality_scores or making it global.
 *****************************************************************************/
#ifndef MAXN
#define MAXN 100
#endif
static void initialize_mw_dist
  (double**** distributions)
{
  int i,j;
  *distributions = (double***)mymalloc((MAXN+1)*sizeof(double**));
  if (verbosity > NORMAL_VERBOSE) {
    fprintf(stderr, "Initializing MW distribution data structure with N=%d\n", MAXN);
  }
  for (i=0; i<=MAXN; i++) {
     (*distributions)[i] = (double**)mymalloc((MAXN+1)*sizeof(double*));
    for (j=0; j<=MAXN; j++) {
      /* We allocate the actual distributions later. These
	  placeholders let us tell when they are done already. */
      (*distributions)[i][j] = (double*)NULL;
    }
  }
}


/*****************************************************************************
 * Free the data structures used to store the mann-whitney
 * distributions. Not used and not tested yet!!!!
 *****************************************************************************/
void free_mw_dist
  (double**** distributions)
{
  int i,j;
  for (i=0; i<=MAXN; i++) {
    for (j=0; j<=MAXN; j++) {
      free((*distributions)[i][j]);
    }
    free((*distributions)[i]);
  }
  free(*distributions);
}


/*****************************************************************************
 * Calculate the mann-whitney u distribution. This must be called such
 * that num1 <= num2.
 *****************************************************************************/
static void calc_mwu_dist
  (int num1,
   int num2,
   double*** distribution)
{
  int k, m, u;
  double sum;
  int* klist = NULL;
  int groupmaxrank = num1*num2 + num1*(num1 + 1)/2;
  int samples = num1 + num2;
  int groupsumrank;

  assert(num1 <= num2 && num1 > 0);

  /* prepare space for the distribution */
  distribution[num1][num2] = (double*)mymalloc(groupmaxrank*sizeof(double));
  for (k=0; k<groupmaxrank; k++) 
    distribution[num1][num2][k] = 0.0;
  
  /* Initialize "group1". This is the starting subsample taken as group
   *  1 (1,2,3,...numgroup1) - the lowest possible set of ranks. */
  klist = (int*)mymalloc(num1*sizeof(int));
  for (k=0; k<num1; k++)
    klist[k] = k+1;

  /* for each possible increasing arrangement of values in the group,
   *   calculate u and add to the distribution
   */
  u = 1;
  m = 0;
  /* could equivalently use the binomial coeff for (num1, samples) */
  while (u>0) { 
    if (m>0)
      next_lex(klist, samples, num1);

    groupsumrank = 0;
    for (k=0; k<num1; k++) {
      groupsumrank += klist[k];
    }
    u = groupmaxrank - groupsumrank; 
    distribution[num1][num2][u]++;
    m++;
  }
  
  /* convert the counts to 1 - cumulative probabilities. Do from the
     right to avoid roundoff errors */
  sum = 0.0;
  for (k=groupmaxrank-1; k>=0; k--) { 
    sum+= distribution[num1][num2][k]/m;
    distribution[num1][num2][k] =  sum;
  }
}

/*****************************************************************************
 * Normal approximation to the mann-whitney test. Invoked when N is
 * too large. Single tailed.
 *****************************************************************************/
static double normal_approx_mann_whitney
(int u,
 int num1,
 int num2)
{
  double mean, stdev, z;
  
  mean = num1*num2/2;
  stdev = sqrt(num1*num2*(num1+num2+1)/12); /* Zar, Biostatistics, pg 151 */
  z = fabs(u - mean) / stdev; /* note we do not use the continuity correction (only needed if p>0.05 and we don't care) */
  if (z == 0.0) { /* this can definitely happen */
    return 0.5;
  } else {
    return 1.0 - 0.5*erfc(-z/sqrt(2.0)); 
  }
}

/*****************************************************************************
 * Calculate the two-sided pvalue for a mann-whitney u test
 * score. Call so num1 <= num2 (Paul Pavlidis)
 *****************************************************************************/

/* maxperms: the number of permutations before we switch to the normal
   approximation. 177 million permuatations takes about one minute on
   a 1 GHz intel processor */
#ifndef MAXPERMS
#define MAXPERMS 20e6
#endif
static double calc_mwu_pvalue
  (int u,
   int num1,
   int num2
   ) 
{
  static double*** distributions = NULL;
  static BOOLEAN_T init = FALSE;
  double numperms;

  assert(num1 <= num2);
  // If this is too large, use the normal approximation.
  numperms = binomial_coeff(num1 + num2, num1); 
  if ((num1 <= MAXN) && (num2 <= MAXN) && (numperms < MAXPERMS)) {
    if (verbosity >= DUMP_VERBOSE) {
      fprintf(stderr, "Using the exact Mann Whitney U test: ");
      fprintf(stderr, "there are %.0f permutations for this class\n",
	      numperms);
    }
    if (!init) {
      initialize_mw_dist(&distributions);
      init = TRUE;
    }
    /* Lazy calculation of distributions. */
    assert(distributions != NULL);
    if (distributions[num1][num2] == NULL) {
      calc_mwu_dist(num1, num2, distributions);
    }
    return 2 * distributions[num1][num2][u]; /* 2 for two tailed */
  } else {
    if (verbosity >= DUMP_VERBOSE) {
      fprintf(stderr, "Using the normal approximation to the Mann Whitney");
      fprintf(stderr, " U test: there are %g permutations", numperms);
      fprintf(stderr, " for this class\n");
    }
    return 2*normal_approx_mann_whitney(u, num1, num2);
  }
}


/*****************************************************************************
 * Compute the Fisher criterion score:
 *
 *        (mean1 - mean2)^2
 *        -----------------
 *           var1 + var2
 *
 *****************************************************************************/
double fisher_score
  (double mean1,
   double mean2,
   double variance1,
   double variance2)
{
  double variance_sum  = variance1 + variance2;
  double mean_diff = mean1 - mean2;

  // Avoid divide by zero.
  if (variance_sum == 0.0) {
    return(0.0);
  }

  if (isnan(variance_sum)) {
    return(0.0);
  }

  return((mean_diff * mean_diff) / variance_sum);
}

/*****************************************************************************
 * Given two arrays, compute the necessary statistics and pass them to
 * the Fisher score calculator.
 *****************************************************************************/
static double fisher_score_from_data
  (ARRAY_T* positive_data,
   ARRAY_T* negative_data)
{
  double mean1 = ave_array(positive_data);
  double mean2 = ave_array(negative_data);
  double variance1 = array_variance(positive_data);
  double variance2 = array_variance(negative_data);

  if (FSELECT_DEBUG) {
    fprintf(stderr, "mean1=%g mean2=%g var1=%g var2=%g\n",
	    mean1, mean2, variance1, variance2);
  }

  return(fisher_score(mean1, mean2, variance1, variance2));
}

/*****************************************************************************
 * Compute the p-value associated with a given t statistic.  Returns
 * the negative log of the p-value.
 *
 * Based upon p. 616 of Numerical Recipes in C.
 *****************************************************************************/
static double ttest_score_to_log_pvalue
  (int    degrees,  // Number of degrees of freedom.
   double tstat)    // t statistic
{
  return -my_log10(betai(0.5 * degrees, 0.5, degrees / 
			 (degrees + tstat * tstat)));
}

/*****************************************************************************
 * Compute the standard t-test score, or Welch's approximation
 * thereof.
 *
 * The standard t-test score is
 *
 *            |mean1 - mean2| 
 *         -----------------------
 *         sqrt(var/n1 + var/n2)
 *
 * where mean1 and mean2 are the means of the two sets of data, n1 and
 * n2 are the numbers of examples in each data set, and var is the
 * pooled variance, which is defined as follows:
 *
 *           sum-of-squares1 + sum-of-squares2
 *    var = ----------------------------------
 *                    n1 + n2 - 2
 *
 * Welch's approximation does not assume that the underlying variances
 * are equal.  The Welch's t-test score is:
 *
 *             |mean1 - mean2| 
 *         -----------------------
 *         sqrt(var1/n1 + var2/n2)
 *
 * This function computes either t-test score and returns the negative
 * log p-value (in base 10).  The p-value computation requires the
 * degrees of freedom, dof.  For the standard t-test,
 *
 *     dof = n1 + n2 - 2.
 *
 * For the Welch version, it is more complex. Define Ai =
 * (vari)/ni, then the degrees of freedom is
 *
 *                            (A1 + A2)^2
 *     dof = floor ---------------------------------
 *                (A1^2)/(n1 - 1) + (A2^2)/(n_2 - 1)
 *
 *****************************************************************************/
double ttest
  (BOOLEAN_T welch,
   double mean1,
   double mean2,
   double var1,
   double var2,
   double num1,
   double num2)
{
  double degrees_of_freedom = 0.0;
  double variance_sum = 0.0;
  double score = 0.0;

  // Standard t-test.
  if (!welch) {
    double var;

    // Compute the degrees of freedom.
    degrees_of_freedom = num1 + num2 - 2.0;

    // Compute the pooled variance.
    if (degrees_of_freedom == 0.0) {
      return(0.0);
    }
    var = (((num1 - 1.0) * var1) + ((num2 - 1.0) * var2)) / degrees_of_freedom;

    // Compute the denominator of the t-test score.
    if ((num1 == 0.0) || (num2 == 0.0)) {
      return(0.0);
    }
    variance_sum = (var/num1) + (var/num2);
  }

  // Welch's approximation.
  else {
    double A1;
    double A2;
    double numerator;
    double denominator;

    // Compute the Ai values.
    if ((num1 == 1.0) || (num2 == 1.0)) {
      return(0.0);
    }
    A1 = var1 / num1;
    A2 = var2 / num2;

    // Compute the degrees of freedom.
    numerator = (A1 + A2) * (A1 + A2);
    denominator = ((A1 * A1) / (num1 - 1.0)) + ((A2 * A2) / (num2 - 1.0));
    degrees_of_freedom = numerator / denominator;

    // Compute the denominator of the t-test score.
    if ((num1 == 0.0) || (num2 == 0.0)) {
      return(0.0);
    }
    variance_sum = (var1/num1) + (var2/num2);
  }

  // Compute the t-test score.
  if (variance_sum == 0.0) {
    return(0.0);
  }
  score = fabs(mean1 - mean2) / sqrt(variance_sum);

  if (FSELECT_DEBUG) {
    //    fprintf(stderr, "n1=%g n2=%g mean1=%g mean2=%g var1=%g var2=%g\n",
    //	    num1, num2, mean1, mean2, var1, var2);
    fprintf(stderr, "ttest score=|%g - %g| / sqrt(%g) = %g, dof=%g, p=%g\n",
	    mean1, mean2, variance_sum, score, degrees_of_freedom, ttest_score_to_log_pvalue(degrees_of_freedom, score));
  }

  // Return the negative log p-value.
  return ttest_score_to_log_pvalue(degrees_of_freedom, score);
}

/*****************************************************************************
 * Given two arrays, compute the necessary statistics and pass them to
 * the t-test calculator.
 *****************************************************************************/
static double ttest_from_data
  (BOOLEAN_T welch,
   ARRAY_T* positive_data,
   ARRAY_T* negative_data)
{
  double mean1 = ave_array(positive_data);
  double mean2 = ave_array(negative_data);
  double var1 = array_variance(positive_data);
  double var2 = array_variance(negative_data);
  double num1 = (double)get_array_length(positive_data);
  double num2 = (double)get_array_length(negative_data);
  return(ttest(welch, mean1, mean2, var1, var2, num1, num2));
}


/*****************************************************************************
 * Mann-whitney u-test (Paul Pavlidis)
 * 
 * The input is rank-transformed data!
 * 
 * The statistic is:
 *           
 *              U = n1*n2 + n1*(n1+1)/2 - R1
 *
 * Where n1 and n2 are the sizes of the two groups, and R1 is the sum
 * of the ranks in group1. This value is used to look up the p value,
 * assuming n1 is the smaller.  Otherwise, swap n1 and n2 in the
 * equation. Note that in this implementation, the ranks start at 1,
 * not zero.
 *****************************************************************************/
double mannwhitney
  (ARRAY_T* positive_data,
   ARRAY_T* negative_data)
{
  int num1 = get_array_length(positive_data);
  int num2 = get_array_length(negative_data);
  int u;
  double score = 0.0;

  // calculate the u score
  u = calc_mwu_score(positive_data, negative_data);
  assert(u > 0);

  // get the pvalue
  if (num1 < num2) {
    score = calc_mwu_pvalue(u, num1, num2);
  } else {
    score = calc_mwu_pvalue(u, num2, num1);
  }
  // Return the negative log p-value.
  return -my_log10(score);
}


/****************************************************************************
 * TNoM feature selection, implemented according to the algorithm
 * described in "Tissue Classification with Gene Expression Profiles",
 * Ben-Dor et al., JCompBio, vol7,pp559-583.
 *
 * The rank transformed data is used and assumed to be the input. All
 * possible rank thresholds are tested to see which yields the lowest
 * classification error. The error rate at that threshold is the TNoM
 * score. Given n samples, there are n thresholds (i ranging from 1 to
 * n) to be tested. For each i, count how many positives are
 * below(above) and negatives are above(below) the threshold.
 * (Paul)
 ***************************************************************************/
double tnom
  (ARRAY_T* positive_data,
   ARRAY_T* negative_data)
{
  int numpos, numneg;
  int bestscore;
  int numsamples;
  int i,j, score_one, score_two, sample;
  int numpos_below, numneg_below, numpos_above, numneg_above;

  numpos = get_array_length(positive_data);
  numneg = get_array_length(negative_data);
  numsamples = numpos + numneg;
  bestscore = numsamples; /* the worst we could do is to get all errors */

  /* test each threshold */
  for (i=0; i<numsamples; i++) {
    numpos_below = numneg_below = numpos_above = numneg_above = 0;
    /* count how many positives are above or below the threshold */
    for (j=0; j<numpos; j++) {
      sample = get_array_item(j, positive_data);
      if (sample < i) {
	numpos_below++;
      } else if (sample >= i) {
	numpos_above++;
      } /* not sure how to resolve ties, ben-dor doesn't say? */
    }

    /* count how many negatives are above or below that threshold */
    for (j=0; j<numneg; j++) {
      sample = get_array_item(j, negative_data);
      if (sample < i) {
	numneg_below++;
      } else if (sample >= i) {
	numneg_above++;
      }
    }

    score_one = numsamples - (numpos_below + numneg_above); /* number of errors */
    score_two = numsamples - (numneg_below + numpos_above); 

    if (score_one < score_two && score_one < bestscore) {
      bestscore = score_one;
    } else if (score_two < score_one && score_two < bestscore) {
      bestscore = score_two;
    }
  }
  return (numsamples - bestscore); /* make big numbers better */
} /* tnom */

/*****************************************************************************
 * Extract from a given array all the values for which the
 * corresponding value in a second array is equal to a given value.
 *****************************************************************************/
static ARRAY_T* extract_one_class
  (double   target_value,
   ARRAY_T* data_array,
   ARRAY_T* class_array)
{
  int num_data;
  int i_data;
  int num_selected;
  int i_selected;
  ARRAY_T* extracted_data;

  // Determine how many items are selected.
  num_data = get_array_length(data_array);
  assert(get_array_length(class_array) == num_data);
  i_selected = 0;
  for (i_data = 0; i_data < num_data; i_data++) {
    if (get_array_item(i_data, class_array) == target_value) {
      i_selected++;
    }
  }
  num_selected = i_selected;

  // Allocate the array.
  extracted_data = allocate_array(num_selected);

  // Fill the array.
  i_selected = 0;
  for (i_data = 0; i_data < num_data; i_data++) {
    if (get_array_item(i_data, class_array) == target_value) {
      set_array_item(i_selected, get_array_item(i_data, data_array), 
		     extracted_data);
      i_selected++;
    }
  }

  return(extracted_data);
}

/*****************************************************************************
 * SAM gene-specific scatter measurement for one gene (this is the
 * same as the ttest)
 *****************************************************************************/
static double calculate_sam_scatter
  (ARRAY_T* positive_data,
   ARRAY_T* negative_data)
{
  double a, ssqpos, ssqneg;
  double npos, nneg;

  npos = (double)get_array_length(positive_data);
  nneg = (double)get_array_length(negative_data);
  ssqpos = array_sumsquarederror(positive_data);
  ssqneg = array_sumsquarederror(negative_data);
  a = (ssqpos + ssqneg)/(npos + nneg - 2.0);
  a = a/npos + a/nneg;
  return sqrt(a);
} /* calculate_sam_scatter */


/*****************************************************************************
 * SAM relative difference for one gene
 *****************************************************************************/
static double calculate_sam_relative_difference
  (ARRAY_T* positive_data,
   ARRAY_T* negative_data,
   double szero)
{
  double mean1, mean2, s;

  mean1 = ave_array(positive_data);
  mean2 = ave_array(negative_data);
  s= calculate_sam_scatter(positive_data, negative_data);
  return fabs((mean1 - mean2)/(s + szero));
} /* calculate_sam_relative_difference */


/*****************************************************************************
 * For a gene+class labels calculate the sam scatter score.
 *****************************************************************************/
static ATYPE compute_pooled_variance
(ARRAY_T* data_array,
 ARRAY_T* class_array)
{
  ARRAY_T* positive_data = NULL;
  ARRAY_T* negative_data = NULL;
  ATYPE returnvalue;
  //  double degrees_of_freedom = 0.0;
  //  double var1, var2;
  //  double num1, num2;

  positive_data = extract_one_class(1, data_array, class_array);
  negative_data = extract_one_class(-1, data_array, class_array);
  returnvalue =  calculate_sam_scatter(positive_data, negative_data);
  free_array(positive_data);
  free_array(negative_data);
  return returnvalue;
}
			      

/*****************************************************************************
 * Determine quantiles of an array. Array must be be preallocated with
 * 100 spaces (move this to array.c?)
 *****************************************************************************/
static void quantiles
(ARRAY_T* data,
 ARRAY_T* quantiles)
{
  ARRAY_T* sorted_data;
  int i, qindex;
  int n = get_array_length(data);

  sorted_data = allocate_array(n);
  copy_array(data, sorted_data);
  sort_array(FALSE, sorted_data);
  if (n < 100) {
    die("Computing quantiles on a data set with too few rows.");
  }
  if (get_array_length(quantiles) < 100) {
    die("Quantile array must be of size 100 at least.");
  }
  
  for (i=0; i<100; i++) {
    qindex = (int)((double)i/100.00 * (double)n);
    set_array_item(i, get_array_item(qindex, sorted_data), quantiles);
    //    fprintf(stderr, "%d %.4f\n", qindex,  get_array_item(qindex,sorted_data));
  }
  //  fprintf(stderr, "\n");
}

/*****************************************************************************
 * SAM szero calculation. From description in SAM user's manual.  The
 * hard-coded values of 100 are the number of quantiles of the data,
 * and is fixed at 100 by definition.
 *****************************************************************************/
static double compute_szero
(MATRIX_T* data_matrix,
 ARRAY_T* class_array)
{
  ARRAY_T* var_quantiles;
  int num_features, i_feature;
  int i,j,k;
  int STEPSIZE = 5; // numbins/stepsize = number of candidate s0 values to generate.
  int BIG = 1e8; // a big number;

  ARRAY_T* variances = NULL;
  ARRAY_T* positive_data = NULL;
  ARRAY_T* negative_data = NULL;
  ARRAY_T* cvvs = allocate_array(100); // we will generate 20 candidate s0 values.
  double szero = 0.0;
  double median_d = 0.0;
  double mad_d = 0.0;
  double meanv, varv, mincv;
  int mincvindex;

  /* compute the variance values for all the data */
  num_features = get_num_cols(data_matrix);

  variances = allocate_array(num_features);
  for (i_feature = 0; i_feature < num_features; i_feature++) {
    ARRAY_T* feature = get_matrix_column(i_feature, data_matrix);
    set_array_item(i_feature, 
		   compute_pooled_variance(feature, class_array),
		   variances);
    free_array(feature);
  }

  /* compute the quantiles of the variance */
  var_quantiles = allocate_array(100);
  quantiles(variances, var_quantiles);

  /* for alpha in 0,0.5,0.1....1.0 compute v_j = mad(d_j^alpha|s_i in
     [q_j,q_j+1]), j = 1,2,...100, where mad is the median absolute
     deviation from the median divided by 0.64. */
  for (i=0; i<100; i+=STEPSIZE) {
    int m;
    double salpha= get_array_item(i, var_quantiles);
    ARRAY_T* vscores = allocate_array(100);

    for (k=0; k<99; k++) {
      ARRAY_T* window_indices; 
      ARRAY_T* window_dvalues;
      //      window_indices = allocate_array((int)((num_features/100)+1)); // plus one to make sure we don't run into roundoff problems. 
      window_indices = allocate_array((int)((num_features/10)+1)); // plus one to make sure we don't run into roundoff problems.  // I don't know why, but /100 was not working for some data sets!!
      window_dvalues = allocate_array(get_array_length(window_indices));

      // need the genes in the window k, k+1.  the brute force way is
      // to find the genes where the variance falls in the window.
      m = 0;
      for (j=0; j<num_features;j++) {
	if (get_array_item(j, variances) > get_array_item(k, var_quantiles) && get_array_item(j, variances) <= get_array_item(k+1, var_quantiles)) {
	  set_array_item(m, j, window_indices);
	  m++;
	}
      }

      // calculate d for all the genes in window_indices. Note that
      // salpha we're testing changes with each iteration of the outer loop.
      for (j=0; j<m; j++) {
	ARRAY_T* feature = get_matrix_column(get_array_item(j, window_indices), data_matrix);
	positive_data = extract_one_class(1, feature, class_array);
	negative_data = extract_one_class(-1, feature, class_array);
	set_array_item(j, calculate_sam_relative_difference(positive_data, negative_data, salpha), window_dvalues);
	free_array(feature);
	free_array(positive_data);
	free_array(negative_data);
      }

      // calculate the median of the d values
      median_d = compute_median(window_dvalues);
      // calculate the median absolute deviation from the mean.
      for (j=0; j<m; j++) {
	set_array_item(j, fabs(get_array_item(j, window_dvalues) - median_d), window_dvalues);
      }
      mad_d = compute_median(window_dvalues)/0.64; // 0.64 is a magic number...
      set_array_item(k, mad_d, vscores);
      free_array(window_indices);
      free_array(window_dvalues);
    } // loop over 100 windows

    /* calculate coefficient of variation of v_j values */
    meanv = ave_array(vscores);
    varv = array_variance(vscores);
    set_array_item(i, (varv*varv)/meanv, cvvs);
    free_array(vscores);
  } // loop over 20 alpha values.

  /* alphahat =  argmin(cv(alpha)) */
  mincv = BIG; 
  mincvindex = 0;
  for (j=0; j<100; j+=STEPSIZE) {
    if (get_array_item(j, cvvs) < mincv) {
      mincv = get_array_item(j, cvvs);
      mincvindex = j;
    }
  }

  /* compute and return s_ohat = salphahat. */
  //  fprintf(stderr, "%g %d %g\n", mincv, mincvindex, get_array_item(mincvindex, var_quantiles) );
  szero = get_array_item(mincvindex, var_quantiles);

  free_array(var_quantiles);
  free_array(variances);
  free_array(cvvs);
  return(szero);
}


/*****************************************************************************
 * A much faster variant of the SAM szero calculation. This uses just
 * a percentile of the sample standard deviations. Based on the
 * description in the SAM tech report "Microarrays and their use in a
 * comparative experiment", where the value is referred to as a0.  The
 * tech report actually says 90th percentile, but this yields entirely
 * different values for szero than the slow version. 

 * Here I choose a low percentile as a heuristic, because
 * experience shows that the slow calculation always chooses a very
 * low value for s0.
 * 
 * The hard-coded value of 100 is the number of quantiles of the
 * data, and is fixed at 100 by definition.
 *
 * A simple test suggests that this method is about 10-20x faster than
 * the slow method, and gives nearly identical results.
 *****************************************************************************/
static double compute_szero_fast
(MATRIX_T* data_matrix,
 ARRAY_T* class_array)
{
  ARRAY_T* variances = NULL;
  ARRAY_T* var_quantiles;  
  int num_features, i_feature;
  double szero;
  int QUANTILE_SELECTED = 0; // a heuristic for selecting s0
  num_features = get_num_cols(data_matrix);

  variances = allocate_array(num_features);
  for (i_feature = 0; i_feature < num_features; i_feature++) {
    ARRAY_T* feature = get_matrix_column(i_feature, data_matrix);
    set_array_item(i_feature, 
		   compute_pooled_variance(feature, class_array),
		   variances);
    free_array(feature);
  }
  var_quantiles = allocate_array(100);
  quantiles(variances, var_quantiles);
  szero = get_array_item(QUANTILE_SELECTED, var_quantiles);
  free_array(variances);
  free_array(var_quantiles);
  return szero;
}


/*****************************************************************************
 * Compute the quality of a single feature.
 *****************************************************************************/
static double compute_quality
  (FSELECT_T feature_select,
   ARRAY_T*  one_feature,
   ARRAY_T*  class_array,
   double szero)
{
  ARRAY_T* positive_data = NULL;
  ARRAY_T* negative_data = NULL;
  ARRAY_T* ranked_one_feature  = NULL;
  double quality = 0.0;

  if (feature_select == MANNWHITNEY_FSELECT || feature_select == TNOM_FSELECT) {
    /* Operate on a rank transformed copy of the data. The rank
       transformation is done here because we need to split the data
       into the two classes first (P.P.)  */
    
    rank_transform_array(one_feature);
    ranked_one_feature = allocate_array(get_array_length(one_feature));
    copy_array(one_feature, ranked_one_feature);
    
    // Split the data into positive and negative subsets.
    positive_data = extract_one_class(1, ranked_one_feature, class_array);
    negative_data = extract_one_class(-1, ranked_one_feature, class_array);

    assert(get_array_length(positive_data) +
	 get_array_length(negative_data) == 
	 get_array_length(ranked_one_feature));
  
  } else if ((feature_select == FISHER_FSELECT) || 
	     (feature_select == TTEST_FSELECT) ||
	     (feature_select == WELCH_FSELECT)) {
    
    // Split the data into positive and negative subsets.
    positive_data = extract_one_class(1, one_feature, class_array);
    negative_data = extract_one_class(-1, one_feature, class_array);

    assert(get_array_length(positive_data) + get_array_length(negative_data)
	   == get_array_length(one_feature));
  } else {
    // Hopefully, this is unreachable.
    die("Invalid or no feature selection metric specified.");
  }
  
  if (get_array_length(positive_data) <= 1 || get_array_length(negative_data) <= 1) {
    die("Feature selection is not possible on this data set:\nYou must have at least 2 positive and 2 negative examples.\nThis data set has %d positive and %d negative examples.", 
	get_array_length(positive_data), get_array_length(negative_data));
  }
  
  // Compute the quality score.
  if (feature_select == MANNWHITNEY_FSELECT) {
    quality = mannwhitney(positive_data, negative_data);
  } else if (feature_select == TNOM_FSELECT) {
    quality = tnom(positive_data, negative_data);
  } else if (feature_select == FISHER_FSELECT) {
    quality = fisher_score_from_data(positive_data, negative_data);
  } else if (feature_select == TTEST_FSELECT) {
    quality = ttest_from_data(FALSE, positive_data, negative_data);
  } else if (feature_select == WELCH_FSELECT) {
    quality = ttest_from_data(TRUE, positive_data, negative_data);
  } else if (feature_select == SAM_FSELECT || feature_select == SAM_FSELECT_FAST) {
    quality = calculate_sam_relative_difference(positive_data, negative_data, szero);
  }

  if (ranked_one_feature != NULL) {
    free_array(ranked_one_feature);
  }
  free_array(negative_data);
  free_array(positive_data);

  if (FSELECT_DEBUG) {
    fprintf(stderr, "quality=%g\n", quality);
  }

  return(quality);
}

/*****************************************************************************
 * Given a data matrix and a set of binary labels, compute the quality
 * of each feature (column) in the matrix.
 *****************************************************************************/
static ARRAY_T* compute_quality_scores
  (FSELECT_T feature_select,
   MATRIX_T* data_matrix,
   ARRAY_T*  class_array) 
{
  ARRAY_T* scores;
  int num_features;
  int i_feature;
  double szero = 0.0;

  if (feature_select == SAM_FSELECT) {
    szero = compute_szero(data_matrix, class_array);
    if (verbosity > NORMAL_VERBOSE)
      fprintf(stderr, "szero: %.3f\n", szero);
  } else if (feature_select == SAM_FSELECT_FAST) {
    szero = compute_szero_fast(data_matrix, class_array);
    if (verbosity > NORMAL_VERBOSE)
      fprintf(stderr, "szero (fast calc.): %.3f\n", szero);
  }

  // The score array has as many entries as columns in the data.
  num_features = get_num_cols(data_matrix);
  scores = allocate_array(num_features);

  // Compute each score.
  if (verbosity > NORMAL_VERBOSE) {
    fprintf(stderr, "Computing quality scores.\n");
  }
  for (i_feature = 0; i_feature < num_features; i_feature++) {
    ARRAY_T* feature = get_matrix_column(i_feature, data_matrix);
    set_array_item(i_feature, 
		   compute_quality(feature_select, feature, class_array, szero),
		   scores);
    free_array(feature);
  }
  return(scores);
}


/*****************************************************************************
 * Choose the score cutoff below which features are discarded.
 *****************************************************************************/
static double set_selection_cutoff
  (double   threshold,
   THRESH_T thresh_type,
   ARRAY_T* scores)
{
  double return_value;
  ARRAY_T* sorted_scores;
  int ith;

  return_value = 0.0;
  if (thresh_type == VALUE_THRESH) {
    return_value = threshold;

  } else {

    // Make a sorted copy of the list of scores.
    sorted_scores = allocate_array(get_array_length(scores));
    copy_array(scores, sorted_scores);
    sort_array(TRUE, sorted_scores);

    // Locate the cutoff item in the sorted array.
    if (thresh_type == PERCENT_THRESH) {
      ith = (int)((threshold / 100.0) * get_array_length(scores)) - 1;
      if (ith < 0) ith = 0.0;
      assert(ith < get_array_length(scores));
      return_value = get_array_item(ith, sorted_scores);
    } else if (thresh_type == NUMBER_THRESH) {
      // Don't select more than the given number of features.
      if (threshold > get_array_length(sorted_scores)) {
	fprintf(stderr, 
		"Warning: Requested more features than are available.\n");
	threshold = get_array_length(sorted_scores);
      }
      return_value = get_array_item(((int)threshold) - 1, sorted_scores);
    }
    free_array(sorted_scores);
  }
  return(return_value);
}

/***************************************************************************
 * Reduce the number of columns in an RDB matrix according to a given
 * set of quality scores and selection criteria.
 *
 * We build the transpose of the matrix, because it's easier to grow
 * matrices that way.
 *
 * If STRICT_THRESH is true, when using the NUMBER threshold, we only
 * select exactly that number of features, instead of taking all tied
 * features.
 ***************************************************************************/
typedef struct feat_score {
  int index;
  double score;
  double secondary_score;
} FEAT_SCORE_T;

static int sort_feat_score_compare
  (const void* elem1,
   const void* elem2)
{
  double num1 = (*(FEAT_SCORE_T**)elem1)->score;
  double num2 = (*(FEAT_SCORE_T**)elem2)->score;
  double num3 = (*(FEAT_SCORE_T**)elem1)->secondary_score;
  double num4 = (*(FEAT_SCORE_T**)elem2)->secondary_score;
  //  int k;
  //  long j;
  static int first_time = 1;
  if (first_time) {
    first_time = 0;
    srand( (long)time( NULL ));
  }
  
  // Note that this does a reverse sort - bigger values are first.
  if (num1 > num2) {
    return(-1);
  } else if (num1 < num2) {
    return(1);
  } else { ///* Primary score is tied, try to break tie with secondary score */
    if (num3 > num4) {
      return(-1);
    } else if (num3 < num4) {
      return(1);
    } else {
      return(0); /* They are tied for both scores */
    }
    //    j = 6;
    //    k = j & 1;
    //    fprintf(stderr, "%d : %d\n", j, k);
    //    if (1) {
    //   return(-1);
    //} else {
    //  return(1);
    // }
    //    return(0);
  }
  return(0);
}

static MATRIX_T* throw_out_features
  (MATRIX_T*      data_matrix,
   STRING_LIST_T* feature_names,
   ARRAY_T*       scores,
   //   ARRAY_T*       secondary_scores, 
   THRESH_T       thresh_type, 
   double         threshold,
   STRING_LIST_T* selected_feature_names)
{
  MATRIX_T*      reduced_matrix;
  MATRIX_T*      transposed_matrix;
  ARRAY_T*       one_feature;
  int            num_features;
  int            i_feature;
  int            num_selected;
  double         score;
  double         actual_threshold;
  int            num_finally_selected;
  FEAT_SCORE_T** features_taken_first_pass;

  // Choose value below which features are discarded.
  actual_threshold = set_selection_cutoff(threshold, thresh_type, scores);
  if (verbosity > NORMAL_VERBOSE) {
    if (thresh_type == NUMBER_THRESH) {
      fprintf(stderr, "Retaining features with scores better than %g, ",
	      actual_threshold);
      fprintf(stderr, "up to %d features.\n", (int)threshold);
    } else {
      fprintf(stderr, "Retaining features with scores better than %g.\n",
	      actual_threshold);
    }
  }

  num_features = get_num_cols(data_matrix);
  features_taken_first_pass = 
    (FEAT_SCORE_T**)mymalloc(num_features * sizeof(FEAT_SCORE_T*));
  num_selected = 0;

  // First pass: select features that have scores below or equal to threshold.
  srand(time(NULL));
  for (i_feature = 0; i_feature < num_features; i_feature++) {
    score = get_array_item(i_feature, scores);
    if (actual_threshold <= score) {
      
      if (verbosity > NORMAL_VERBOSE) {
	fprintf(stderr, "Retaining feature %d with score %g.\n",
		i_feature, get_array_item(i_feature, scores));
      }
      
      features_taken_first_pass[num_selected]
	= (FEAT_SCORE_T*)mymalloc(sizeof(FEAT_SCORE_T));
      features_taken_first_pass[num_selected]->score = score;

      /* use a random number to break ties */
      features_taken_first_pass[num_selected]->secondary_score 
      	= rand();

      features_taken_first_pass[num_selected]->index = i_feature;
      num_selected++;
    }
  }
  
  // Second pass: if need be, set up so we can whittle down to the
  // actual number of features needed.
  if ((thresh_type == NUMBER_THRESH) && (num_selected > threshold)) {
    if (verbosity > NORMAL_VERBOSE) {
      fprintf(stderr, "Selected %d, so trimming %d features due to tied scores.\n",
	      (int)(num_selected), (int)(num_selected) - (int)(threshold));
    }
    qsort(features_taken_first_pass, num_selected, sizeof(FEAT_SCORE_T*),
	  sort_feat_score_compare);
    num_selected = threshold;
  }
  if (verbosity > NORMAL_VERBOSE) {
    fprintf(stderr, "Selected %d features\n", (int)(num_selected));
  }

  // Do actual selection of the data
  reduced_matrix = allocate_matrix(0, get_num_rows(data_matrix));
  num_finally_selected = 0;
  for (i_feature = 0; i_feature < num_selected; i_feature++) {

    // Extract the column.
    one_feature
      = get_matrix_column(features_taken_first_pass[i_feature]->index, 
			  data_matrix);
    
    // Put it into the reduced matrix.
    grow_matrix(one_feature, reduced_matrix);
    if (selected_feature_names != NULL) {
      add_string(get_nth_string(features_taken_first_pass[i_feature]->index,
				feature_names),
		 selected_feature_names);
    }
    if (verbosity >= DUMP_VERBOSE) {
      fprintf(stderr, "Selecting feature %d with score %f.\n", 
	      features_taken_first_pass[i_feature]->index,
	      features_taken_first_pass[i_feature]->score);
    }

    num_finally_selected++;
    
    // Free the column.
    free_array(one_feature);
    
    // free the structs as we go.
    free(features_taken_first_pass[i_feature]);
  }
  free(features_taken_first_pass);
  assert(num_finally_selected == get_num_rows(reduced_matrix));

  if (num_selected == 0) {
    fprintf(stderr, "Warning: No features meet selection criterion.\n");
  }
  if (verbosity > NORMAL_VERBOSE) {
    fprintf(stderr, "Retaining %d of %d features.\n", num_selected,
	    get_array_length(scores));
  }

  // Transpose the matrix.
  transposed_matrix = transpose_matrix(reduced_matrix);
  free_matrix(reduced_matrix);

  // Put labels on the matrix.
  return transposed_matrix;
}

/***************************************************************************
 * Reduce a given data set using feature selection.
 ***************************************************************************/
void select_features
  (KERNEL_T*      kernel,        // Kernel parameters.
   CLASS_ARRAY_T* classes,       // Classes.
   STRING_LIST_T* feature_names, // Names of all features.
   MATRIX_T*      train_matrix,  // Data from which to score features.
   MATRIX_T*      test_matrix,   // Second data to remove features from.
   ARRAY_T**      score_array,   // Score array.
   STRING_LIST_T* selected_feature_names, // Names of selected features.
   MATRIX_T**     reduced_train, // Training set with fewer features.
   MATRIX_T**     reduced_test)  // Test set with fewer features.
{
  // Compute the quality scores.
  ARRAY_T* class_array = get_class_array(classes);

  myassert(1, get_array_length(class_array) == get_num_rows(train_matrix), 
	   "Number of training examples (%d) does not equal the number of class labels (%d)",  
	   get_array_length(class_array), 
	   get_num_rows(train_matrix));

  // If no feature selection was requested, don't bother.
  if (kernel->feature_select == NONE_FSELECT) {
    return;
  }
  
  // Compute the quality score asscoiated with each feature.
  *score_array = compute_quality_scores(kernel->feature_select,
					train_matrix,
					class_array);

  free_array(class_array);

  // Reduce the train and test matrices according to the scores.
  *reduced_train = throw_out_features(train_matrix,
				      feature_names,
				      *score_array, 
				      kernel->thresh_type,
				      kernel->fthreshold,
				      selected_feature_names);
  if (train_matrix == test_matrix) {
    *reduced_test = *reduced_train;
  } else {
    *reduced_test = throw_out_features(test_matrix,
				       NULL,
				       *score_array, 
				       kernel->thresh_type,
				       kernel->fthreshold,
				       NULL);
  }

}


#ifdef FSELECT_MAIN

/*****************************************************************************
 * A goofy helper function required because the command line macro
 * takes single lines as input.
 *****************************************************************************/
static void set_required_files
 (char** class_filename,
  char** primary_filename,
  char** secondary_filename,
  char*  given_filename) {
  if (*class_filename == NULL) {
    *class_filename = given_filename;
  } else if (*primary_filename == NULL) {
    *primary_filename = given_filename;
  } else if (*secondary_filename == NULL) {
    *secondary_filename = given_filename;
  } else {
    die("Extra unlabeled option: %s", given_filename);
  }
}

/***************************************************************************
 * Main procedure
 ***************************************************************************/

int main(int argc, char *argv[]) {

  // Command-line options.
  char*     class_filename = NULL;  // Input label file.
  FILE*     class_file = NULL;
  int       class_number = 1;       // In a multiclass file, which class to use. First one is 1
  char*     score_filename = NULL;  // Output score file.
  BOOLEAN_T format_line = FALSE;    // Read and write RDB format line?

  // Data structures.
  KERNEL_T*      kernel = allocate_kernel(FALSE); // Store options.
  CLASS_ARRAY_T* classes = NULL;          // The training set classifications.
  MATRIX_T*      reduced_train = NULL;    // Data with features removed.
  MATRIX_T*      reduced_test = NULL;   
  RDB_MATRIX_T*  reduced_test_rdb = NULL; 
  STRING_LIST_T* selected_feature_names = NULL; // Names of selected features.
  ARRAY_T*       score_array = NULL;
  //  ARRAY_T*       secondary_score_array = NULL;

  // Set the default feature selection method.
  kernel->feature_select = FISHER_FSELECT;

  // Parse the command line.
  DO_STANDARD_COMMAND_LINE
    (2,
     DATA_OPTN(1, scores, <file>, score_filename = _OPTION_);
     DATA_OPTN(1, metric, fisher|ttest|welch|mannwhitney|sam|tnom, // hide the samf option for now
	       kernel->feature_select = convert_enum_type_str(_OPTION_,
							      INVALID_FSELECT,
 							      FSELECT_STRS,
	 						      NUM_FSELECT_T));
     DATA_OPTN(1, threshtype, percent|number|value,
	       kernel->thresh_type = convert_enum_type_str(_OPTION_,
							   INVALID_THRESH,
							   THRESH_STRS,
							   NUM_THRESH_T));
     DATA_OPTN(1, fthreshold, <value>, set_fthreshold(kernel, atof(_OPTION_)));
     SIMPLE_FLAG_OPTN(1, rdb, format_line);
     DATA_OPTN(1, verbose, 1|2|3|4|5 (default=2),
	       verbosity = (VERBOSE_T)atoi(_OPTION_));
     DATA_OPTN(1, useclassnumber, <value>, class_number = atoi(_OPTION_));
     NON_SWITCH(1, <labels> <primary> <secondary>,
		set_required_files(&class_filename,
				   &(kernel->train_filename), 
				   &(kernel->test_filename),
				   _OPTION_));
     );

  // Set the default threshold.
  set_default_feature_selection(kernel);

  // Make sure we got the required files.
  if (class_filename == NULL) {
    die("No classification file given.");
  }
  if (kernel->train_filename == NULL) {
    die("No primary data file given.");
  }
  if (kernel->test_filename == NULL) {
    die("No secondary data file given.");
  }

  // Read both data sets.
  if (open_file(kernel->train_filename, "r", TRUE, "primary data",
		"primary data set", &(kernel->train_file)) == 0)
    exit(1);
  set_train_matrix(read_rdb_matrix(format_line, NULL, kernel->train_file),
		   kernel);
  fclose(kernel->train_file);

  if (open_file(kernel->test_filename, "r", TRUE, "secondary data", 
		"secondary data set", &kernel->test_file) == 0)
    exit(1);
  set_test_matrix(read_rdb_matrix(format_line, NULL, kernel->test_file),
		  kernel);
  fclose(kernel->test_file);

  // Read the classifications.
  if (open_file(class_filename, "r", TRUE, "class", "classes",
		&class_file) == 0)
    exit(1);
  classes = read_classifications(format_line, class_file, class_number);
  fclose(class_file);

  // Create list of selected features names.
  selected_feature_names = new_string_list();

  // Select features.
  select_features(kernel,
		  classes,
		  get_col_names(kernel->train_matrix_rdb),
		  kernel->train_matrix,
		  kernel->test_matrix,
		  &score_array,
		  selected_feature_names,
		  &reduced_train,
		  &reduced_test);

  // Print the score file.
  if (score_filename != NULL) {
    FILE* score_file = NULL;
    STRING_LIST_T* col_names = NULL;
    MATRIX_T* score_matrix = NULL;
    RDB_MATRIX_T* score_matrix_rdb = NULL;

    // Convert the scores to a matrix.
    score_matrix = array_to_matrix(FALSE, score_array);

    // Add labels to the score matrix.
    col_names = new_string_list();
    add_string(convert_enum_type(kernel->feature_select, FSELECT_STRS,
				 NUM_FSELECT_T), col_names);
    score_matrix_rdb = rdbize_matrix("fselect",
				     get_col_names(kernel->train_matrix_rdb),
				     col_names,
				     score_matrix);

    // Open the file.
    if (open_file(score_filename, "w", FALSE, "score", "scores",
		  &score_file) == 0)
      exit(1);

    // Print the scores.
    print_rdb_matrix(score_matrix_rdb, format_line, 8, 6, score_file);
    fclose(score_file);
    free_rdb_matrix(score_matrix_rdb);
  }

  // Add labels to the output matrix.
  reduced_test_rdb = rdbize_matrix(get_corner_string(kernel->test_matrix_rdb), 
				   get_row_names(kernel->test_matrix_rdb), 
				   selected_feature_names, 
				   reduced_test);

  // Print the output.
  print_rdb_matrix(reduced_test_rdb, format_line, 8, 6, stdout);

  // Free dynamic memory.
  free_kernel(kernel);
  free_class_array(classes);
  free_matrix(reduced_train);
  free_rdb_matrix(reduced_test_rdb);
  free_string_list(selected_feature_names);
  free_array(score_array);
  //  free_array(secondary_score_array);

  return(0);
}
#endif

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
