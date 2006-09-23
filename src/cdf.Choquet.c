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

  Computation of the c.d.f. of the Choquet integral 
  in the uniform case
  
  Ivan Kojadinovic, 04/2006

*****************************************************************************/

#include <R.h>
#include "core.h"

/*****************************************************************************

  Computation of the divided difference [(x-y)^n_{-} : a_0,...,a_n]

*****************************************************************************/

double div_diff_xn_y_minus(int n, double *a, double y) {

  int i, j, k, r=0, s=0;
  double *dd, *b, *c;
  for (i=0;i<=n;i++) {
    if (a[i] < y)
      r++;
    else
      s++;
  }

  b = (double *)R_alloc(r,sizeof(double));
  c = (double *)R_alloc(s,sizeof(double));
  
  j=0;
  k=0;
  for (i=0;i<=n;i++)
    if (a[i] < y)
      b[j++]=a[i]-y;
    else
      c[k++]=a[i]-y;

  dd = (double *)R_alloc(s+1,sizeof(double));
  
  /* Initialize dd */
  dd[0] = 1.0;
  for (j=1; j<=s; j++)
    dd[j] = 0.0;

  /* Computation by induction */
  for (i=1; i<=r; i++) 
    for (j=1; j<=s; j++) 
      dd[j] = (c[j-1] * dd[j] - b[i-1] * dd[j-1])/(c[j-1] - b[i-1]);

  return dd[s];
}

/*****************************************************************************

  Computation of the divided difference [(x-y)^(n-1)_{+} : a_0,...,a_n]

*****************************************************************************/

double div_diff_xn_1_y_plus(int n, double *a, double y) {

  int i, j, k, r=0, s=0;
  double *dd, *b, *c;
  for (i=0;i<=n;i++) {
    if (a[i] < y)
      r++;
    else
      s++;
  }

  if (r==0 || s==0)
    return 0.0;

  b = (double *)R_alloc(r,sizeof(double));
  c = (double *)R_alloc(s,sizeof(double));
  
  j=0;
  k=0;
  for (i=0;i<=n;i++)
    if (a[i] < y)
      b[j++]=a[i]-y;
    else
      c[k++]=a[i]-y;

  dd = (double *)R_alloc(s+1,sizeof(double));
  
  /* Initialize dd */
  dd[0] = 0.0;
  dd[1] = 1.0/(c[0] - b[0]);

  /* Computation by induction */
  for (j=2; j<=s; j++) 
    dd[j] = - b[0] * dd[j-1]/(c[j-1] - b[0]);
   
  for (i=2; i<=r; i++) 
    for (j=1; j<=s; j++)
      dd[j] = (c[j-1] * dd[j] - b[i-1] * dd[j-1])/(c[j-1] - b[i-1]);
   
  return dd[s];
}

/*****************************************************************************

  p-th lexicographic permutation

*****************************************************************************/

void lex_permut(int n, int p, int *x, int *res) {

  int i,j,k,q,l,ifact;
  for (i=n-1;i>=0;i--) {
    ifact = (int)fact(i);
    p = p % ((i+1)*ifact);
    q = (int)(p / ifact);
    k = x[q];
    
    for (j=0;j<=i;j++)
      if (x[j] == k) {
	l = j;
	break;
      }
    for (j=l;j<i;j++)
      x[j] = x[j+1];
    
    res[n-1-i]=k;
  }
}

/*****************************************************************************

  c.d.f. of the Choquet integral in the uniform case

*****************************************************************************/

void cdf_Choquet(int *n, double *mu, double *y, double *Fy) {

  int i,j;
  double f = fact(*n), cdf = 0.0;
  int *id = (int *)R_alloc(*n,sizeof(int));
  int *sigma = (int *)R_alloc(*n,sizeof(int));
  double *a = (double *)R_alloc(*n+1,sizeof(double));

  for (i=0;i<(int)f;i++) {
 
    for (j=0;j<*n;j++)
      id[j]=j;

    lex_permut(*n, i, id, sigma);

    a[0] = mu[(1<<*n) - 1];
    for(j=1; j<*n; j++)
      a[j] = mu[subset2binary(sigma + j, *n - j)];
    a[*n] = 0.0; 

    cdf += div_diff_xn_y_minus(*n,a,*y);
  }
  
  *Fy = cdf/f;
}

/*****************************************************************************

  p.d.f. of the Choquet integral in the uniform case

*****************************************************************************/

void pdf_Choquet(int *n, double *mu, double *y, double *py) {

  int i,j;
  double f = fact(*n), pdf = 0.0;
  int *id = (int *)R_alloc(*n,sizeof(int));
  int *sigma = (int *)R_alloc(*n,sizeof(int));
  double *a = (double *)R_alloc(*n+1,sizeof(double));


  for (i=0;i<f;i++) {
 
    for (j=0;j<*n;j++)
      id[j]=j;

    lex_permut(*n, i, id, sigma);

    a[0] = mu[(1<<*n) - 1];
    for(j=1; j<*n; j++)
      a[j] = mu[subset2binary(sigma + j, *n - j)];
    a[*n] = 0.0; 

    pdf += div_diff_xn_1_y_plus(*n,a,*y);
  }

  *py = pdf * (double)(*n) / f;
}


/*****************************************************************************/
