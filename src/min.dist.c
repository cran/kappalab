/*##############################################################################
#
# Copyright 2005, 2006, 2007 Michel Grabisch and Ivan Kojadinovic   
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

  Minimum distance capacity identification
  Ivan Kojadinovic, 11/2005

*****************************************************************************/

#include <R.h>
#include "core.h"
#include "min.dist.h"


/*****************************************************************************

  Diagonal matrix partially defining the objective function for distance d2

*****************************************************************************/

void objective_function_Choquet_coefficients(int *n, double *D) {

  int i, j, m;
  int pow = 1<<*n;

  m = 0;
  for (i=0; i<*n; i++)
    for (j=0; j<pow; j++)
      if(j & 1<<i)
	D[m++] = gamm(cardinal(j)-1, *n);
      
}

/*****************************************************************************

  Matrix B partially defining the objective function for distance d1

*****************************************************************************/

void objective_function_binary_alternatives(int *n, int *k, int *subset, int *B) {

  int i, j, m;
  int sb = (int)sum_binom(*n,*k);

  m = 0; 
  for (i=1; i<(1<<*n); i++) 
    for (j=1; j<sb; j++) {

      if ((subset[j] & i) == subset[j])
	B[m++] = 1;
      else
	B[m++] = 0;
      /*if (subset[j] == i)
	break;*/
    }
}

/*****************************************************************************

  Matrix Q partially defining the objective function for distance d3

*****************************************************************************/

void objective_function_global_scores(int *n, int *k1, int*k2, int *subset, double *Q) {

  int i, j, m;
  double s;
  int sb1 = (int)sum_binom(*n,*k1);
  int sb2 = (int)sum_binom(*n,*k2);

  m = 0; 
  for (i=1; i<sb1; i++) {
    s = 1.0/(double)(cardinal(subset[i])+1);
    for (j=1; j<sb2; j++) 
      Q[m++] = (s + 1.0/ (double)(cardinal(subset[j])+1)) / (double)(cardinal(subset[i] | subset[j])+2);
    
  }
}

/*****************************************************************************/
