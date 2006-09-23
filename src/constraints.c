/*##############################################################################
#
# Copyright © 2005 Michel Grabisch and Ivan Kojadinovic   
#
# Ivan.Kojadinovic@polytech.univ-nantes.fr
#
# This software is a package for the statistical system GNU R:
# http://www.r-project.org 
#
# This software is governed by the CeCILL license under French law and
# abiding by the rules of distribution of free software.  You can  use, 
# modify and/ or redistribute the software under the terms of the CeCILL
# license as circulated by CEA, CNRS and INRIA at the following URL
# "http://www.cecill.info". 
#
# As a counterpart to the access to the source code and  rights to copy,
# modify and redistribute granted by the license, users are provided only
# with a limited warranty  and the software's author,  the holder of the
# economic rights,  and the successive licensors  have only  limited
# liability. 
#
# In this respect, the user's attention is drawn to the risks associated
# with loading,  using,  modifying and/or developing or reproducing the
# software by the user in light of its specific status of free software,
# that may mean  that it is complicated to manipulate,  and  that  also
# therefore means  that it is reserved for developers  and  experienced
# professionals having in-depth computer knowledge. Users are therefore
# encouraged to load and test the software's suitability as regards their
# requirements in conditions enabling the security of their systems and/or 
# data to be ensured and,  more generally, to use and operate it in the 
# same conditions as regards security. 
#
# The fact that you are presently reading this means that you have had
# knowledge of the CeCILL license and that you accept its terms.
#
##############################################################################
*/


/*****************************************************************************

  Minimum variance capacity identification
  Ivan Kojadinovic, 01/03/2005

*****************************************************************************/

#include <R.h>
#include "core.h"
#include "constraints.h"


/*****************************************************************************

   Monotonicity constraints

*****************************************************************************/

void monotonicity_constraints(int *n, int *k, int *subset, int *A) {
	
  int i,j,l,m;
  int pow = 1<<*n;
  int sb = (int)sum_binom(*n,*k);
  
  m = 0;
  for (i=0; i<*n; i++)
    for (j=1; j<pow; j++)
      if(j & 1<<i)
	for(l=1;l<sb;l++)
	  if (((subset[l] & j) == subset[l]) 
	      && (subset[l] & 1<<i))
	    A[m++] = 1;
	  else
	    A[m++] = 0;
}


/*****************************************************************************

  Constraint C(a) >= C(b)

*****************************************************************************/

void Choquet_preorder_constraint(int *n, int *k, int *subset, double *a, 
				 double *b, double *A) {

  int i,j,l;
  int sb = (int)sum_binom(*n,*k);
  double min_a, min_b;

  for (i=1;i<sb;i++) {
		
    for (j=0;j<*n;j++)
      if (subset[i] & 1<<j) {

	min_a = a[j];
	min_b = b[j];
	break;
      }
    for (l=j+1;l<*n;l++)
      if (subset[i] & 1<<l) {
				
	min_a = INF(min_a, a[l]);
	min_b = INF(min_b, b[l]);
      }

    A[i-1] = min_a - min_b;
  }

} 

/*****************************************************************************

  Constraint phi(a) >= phi(b)

*****************************************************************************/

void Shapley_preorder_constraint(int *n, int *k, int *subset, int *a, int *b, 
			double *A) {
  int j; 
  int sb = (int)sum_binom(*n,*k);
  for (j=1; j<sb; j++) {

    A[j-1] = 0.0;
    
    if (subset[j] & 1<<*a)
      A[j-1] = 1 / (double)cardinal(subset[j]);
    if (subset[j] & 1<<*b)
      A[j-1] -= 1 / (double)cardinal(subset[j]);
      
  }
}

/*****************************************************************************

  Constraint min <= phi(a) <= max

*****************************************************************************/

void Shapley_interval_constraint(int *n, int *k, int *subset, int *a, 
				 double *A) {

  int j; 
  int sb = (int)sum_binom(*n,*k);

  for (j=1; j<sb; j++) {

    if (subset[j] & 1<<*a)
      A[j-1] = 1 / (double)cardinal(subset[j]);
    else
      A[j-1] = 0.0;
   
  }
}

/*****************************************************************************

  Constraint I(ab) >= I(cd)

*****************************************************************************/

void interaction_preorder_constraint(int *n, int *k, int *subset, int *a, 
				     int *b, int *c, int *d, double *A) {

  int l;
  int sb = (int)sum_binom(*n,*k);
  
  for(l=1; l<sb; l++) {

    A[l-1] = 0.0;
    if((subset[l] & 1<<*a) && (subset[l] & 1<<*b))
      A[l-1] = 1/(double)(cardinal(subset[l])-1);
    if((subset[l] & 1<<*c) && (subset[l] & 1<<*d))
      A[l-1] -= 1/(double)(cardinal(subset[l])-1);    
  }
}

/*****************************************************************************

  Constraint min <= I(ab) <= max

*****************************************************************************/

void interaction_interval_constraint(int *n, int *k, int *subset, int *a, 
				     int *b, double *A) {

  int l;
  int sb = (int)sum_binom(*n,*k);
  
  for(l=1; l<sb; l++) {

    if((subset[l] & 1<<*a) && (subset[l] & 1<<*b))
      A[l-1] = 1/(double)(cardinal(subset[l])-1);
    else
      A[l-1] = 0.0;
  }
}

/*****************************************************************************

  Partition {A1,...,Ap}
  Constraint: Mobius(S) = 0 if S not included in one Ai
  Partition given under the form of an index vector

*****************************************************************************/

void inter_additive_constraint(int *n, int *k, int *subset, int *partition, 
			       int *p, double *A) {

  int i,j,l,c;
  int sb = (int)sum_binom(*n,*k);
  /* binary coding of the partition */
  int *part_binary = (int *) R_alloc(*p, sizeof(int));
  /* binary coding of one element of the parition */
  int *set = (int *) R_alloc(*n, sizeof(int));
  int included,  max_card_subset, max_card = 0;

  /* translate the partition to binary */
  for (i=0;i<*p;i++) {

    l = 0;
    for (j=0;j<*n;j++)
      if (partition[j] == (i+1))
	set[l++] = j;
    
    part_binary[i] = subset2binary(set,l);

    /* find max card subset */
    c = cardinal(part_binary[i]);
    if (c > max_card) {
      
      max_card = c;
      max_card_subset = i;
    }
  }
  
  /* subsets whose Mobius has to be zero */
  for(l=1; l<sb; l++) {
    
    included = 0; /* not included */
    A[l-1] = 0;

    if (cardinal(subset[l]) <= max_card) 
      for (i=0;i<*p;i++)
	if((subset[l] & part_binary[i]) == subset[l]) { 
	  
	  included = 1;
	  break;
	}

    if (!included) 
      A[l-1] = 1;
  }
}
