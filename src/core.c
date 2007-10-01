/*##############################################################################
#
# Copyright © 2005, 2006, 2007 Michel Grabisch and Ivan Kojadinovic   
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

  Non-additive measure and integral manipulation routines.
  C functions source file. Binary coding of sets.
  Univers X, X={0,1,2,3,...,n-1} with n=32 maximum.

  Initial work: Michel Grabisch, 21/12/94
  Contributions: Ivan Kojadinovic, 01/03/2005

*****************************************************************************/


#include <string.h>
#include <R.h>
#include "core.h"

/*****************************************************************************

  Conversion list of elements -> binary
  x : set to be converted
  nx : number of elements of x

*****************************************************************************/

int subset2binary(int *x, int nx) {

  int i, b;

  b=0;
  for(i=0; i<nx; i++)
    b += 1 << x[i];

  return(b);
}


/*****************************************************************************

  Conversion binary -> list of elements
  b : binary code to be converted
  x : result in an array (ordered liste of elements)

*****************************************************************************/

void binary2subset(int n, int b, int *x) {

  int i;

  for(i=0; i<n; i++)
    if(b & 1 << i) {

      *x = i;
      x++;
    }
}


/*****************************************************************************
 
   Pointer version of the binary2subset function, interfaced in kappalab.R
   l: the real length of x

*****************************************************************************/

void binary2subsetR(int *n, int *b, int *x, int *l) {

  int i;
  
  *l=0;

  for(i=0; i<*n; i++)
    if(*b & 1 << i) {

      *x = i+1; /* R arrays are 1-based */
      x++;
      (*l)++;
    }
}


/*****************************************************************************

   Sorting by tournament (M. Hebert)
   n: dimension of vec and tourn
   vec: vector to be sorted
   iv: indices of the sorted elements (out)

*****************************************************************************/

void tri(register int n, register int *tourn, register double *vec, 
	 register int *iv) {
  
  int i;
  register int k, kk;
  double x;
  int nb, nb1, nnb1, id;

  tourn--;
  vec--;
  iv--;

  /* adjonction dans le tournoi */
  for(i=1; i<=n; i++) {

    x = vec[i];
    k = i;

    while((k > 1) && (x < vec[tourn[(kk = k >> 1)]])) {
	      
      tourn[k] = tourn[kk];
      k = kk;
    }

    tourn[k] = i;
  }

  for(i=1; i<=n; i++) {

    nb = n - i + 1;
    iv[i] = tourn[1] - 1;

    /* extraction de la racine (min) */
    id = tourn[nb];
    k = 1;
    nb1 = nb + 1;
    nnb1 = nb1 >> 1;

    while(k < nnb1)
      if(vec[tourn[k<<1]] > vec[tourn[(k<<1)|1]]) {

	tourn[k] = tourn[(k << 1) | 1];
	k = (k << 1) | 1;
      }
      else {
		
	tourn[k] = tourn[k << 1];
	k = k << 1;
      }

    if(nb1 == ((k << 1) | 1)) {

      tourn[k] = tourn[k << 1];
      k = k << 1;
    }

    while((k > 1) && (vec[id] <= vec[tourn[(kk = k >> 1)]])) {

      tourn[k] = tourn[kk];
      k = kk;
    }

    tourn[k] = tourn[nb];
  }
}


/*****************************************************************************

  Counts the number of bits equal to 1

*****************************************************************************/

int cardinal(int n) {

  int i;

  for(i=0; n; n >>= 1)
    i += n & 1;

  return(i);
}

/**************************************************************************
    
  Conjugate transform of a set function
  mu: set function 
  The result is written in mu_out

***************************************************************************/

void setfunction2conjugate(double *mu, int *n, double *mu_out) {

  int i;
  int pow = 1<<*n;
  for (i=0;i<pow;i++)
    mu_out[i] = mu[pow-1] - mu[pow-1-i];
  
}

/**************************************************************************
    
  Quick computation of the lower cardinality transform
  c: factor used in the transformation
  mu: set function 
  The result is written in mu

***************************************************************************/

void fast_lower_cardinality_transform(double *mu, double c, int n) {

  int i;
  int j,k;
 
  for (i=1;i<=n;i++)
    for (j=1;j<(1<<i);j+=2)
      for (k=0;k<(1<<(n-i));k++) 
	mu[j*(1<<(n-i))+k] += c*mu[(j-1)*(1<<(n-i))+k];
}

/****************************************************************************

   Computation of the Mobius transform using Robert Kennes algorithm (1990)
   The result is written in mu.

***************************************************************************/

void setfunction2Mobius(double *mu, int *n) {

  fast_lower_cardinality_transform(mu,-1.0,*n);
}

/****************************************************************************

   Computation of the inverse Mobius transform.

   Old version: Robert Kennes algorithm (1990)
   Pre: v is given in binary order.

   void Mobius2setfunction(double *v, int *n)
   {
     fast_lower_cardinality_transform(v,1.0,*n);
   }

   Current version: less efficient
   Pre: a is given in "natural" order

***************************************************************************/

void Mobius2setfunction(int *n, int *k, double *a, int *subset, double *mu) {

  int i,j;
  int sb = (int)sum_binom(*n,*k);
  for (i=0; i<(1<<*n); i++) {

    mu[i] = a[0];
    for (j=1; j<sb; j++) {

      if ((subset[j] & i) == subset[j])
	mu[i] += a[j];
      if (subset[j] == i)
	break;
    }
  }
}

/****************************************************************************

   Normalizes a Mobius transform ( mu(X) = 1 )

***************************************************************************/

void normalize_Mobius(int n, int k, double *a) {

  int i;
  int sb = (int)sum_binom(n,k);
  double s = 0.0;
  for (i=0; i<sb; i++)
    s += a[i];
  for (i=0; i<sb; i++)
    a[i] /= s;
}

/*****************************************************************************

  k-truncation of a set function mu. 
  Returns its truncated Mobius representation.
  Result in a (dim sum_binom(n,k)).

*****************************************************************************/

void k_truncation(int *n, int *k, double *mu, int *subset, double *a) {

  int i;
  setfunction2Mobius(mu, n);
  for (i=0;i<sum_binom(*n,*k);i++)
    a[i] = mu[subset[i]];
}

/*****************************************************************************

  Generates a cardinal set function from an array of size n+1 (csf)
  Post: sf given in binary order

*****************************************************************************/

void cardinal2setfunction(int *n, double *csf, double *sf) {

  int i;
  for (i=0; i<(1<<*n); i++) 
    sf[i] = csf[cardinal(i)];
}

/*****************************************************************************

  Generates an array of size n+1 (csf) from a cardinal set function
  Pre: sf given in "natural" order and may be truncated

*****************************************************************************/

void setfunction2cardinal(int *n, int *k, double *sf, double *csf) {

  int i;
  int sb = 0;
  for (i=0; i<=*k; i++) {

    csf[i] = sf[sb];
    sb += binom(*n,i);
  }
  for (i=*k+1;i<=*n;i++)
    csf[i] = 0.0;
}

/*****************************************************************************

  Calculation of the Choquet integral.
  mu : game (2^n coefficients)
  f : function (n coefficients)

*****************************************************************************/

void Choquet_integral_game(int *n, double *mu, double *f, double *resul) {

  int i;
  int tourn[NMAX];
  int *index = (int *) R_alloc(*n, sizeof(int));

  tri(*n, tourn, f, index);
  *resul = f[index[0]] * mu[(1<<*n)-1];

  for(i=1; i<*n; i++)
    *resul += (f[index[i]] - f[index[i - 1]]) 
      * mu[subset2binary(index + i, *n - i)];

  index = NULL;
}

/*****************************************************************************

  Calculation of the Choquet integral.
  a : Mobius representation of a game
  f : function (n coefficients)

*****************************************************************************/

void Choquet_integral_Mobius(int *n, int *k, double *a, int *subset, 
			     double *f, double *resul) {

  int i,j,l;
  int sb = (int)sum_binom(*n,*k);
  double min_f;

  *resul = 0.0;
  for (i=1;i<sb;i++)
    {
      for (j=0;j<*n;j++)
	if (subset[i] & 1<<j)
	  {
	    min_f = f[j];
	    break;
	  }
      for (l=j+1;l<*n;l++)
	if (subset[i] & 1<<l)
	  min_f = INF(min_f, f[l]);

      *resul += a[i] * min_f;
    }
}

/*****************************************************************************

  Calculation of the Sugeno integral.
  mu : game (2^n coefficients)
  f : function (n coefficients) 

*****************************************************************************/

void Sugeno_integral_game(int *n, double *mu, double *f, double *resul) {

  int i;
  int tourn[NMAX];
  int *index = (int *) R_alloc(*n, sizeof(int));

  tri(*n, tourn, f, index);

  *resul = INF(f[index[0]], mu[(1<<*n) -1]);

  for(i=1; i<*n; i++)
    *resul = SUP(*resul, INF(f[index[i]], mu[subset2binary(index+i, *n-i)]));

  index = NULL;
}

/*****************************************************************************

  Calculation of the Sugeno integral.
  a : Mobius representation of a game 
  f : function (n coefficients) 

  VERY VERY VERY SUBOPTIMAL!!!

*****************************************************************************/

void Sugeno_integral_Mobius(int *n, int *k, double *a, int *subset, 
			    double *f, double *resul) {

  int i,j,s;
  int sb = (int)sum_binom(*n,*k);
  int tourn[NMAX];
  int *index = (int *) R_alloc(*n, sizeof(int));
  double mu;

  tri(*n, tourn, f, index);

  mu = 0.0;
  for (i=1;i<sb;i++)
    mu += a[i];

  *resul = INF(f[index[0]], mu);

  for(i=1; i<*n; i++) {
    s = subset2binary(index+i, *n-i);
    mu = 0.0;
    for (j=1; j < sb; j++)
      if ((subset[j] & s) == subset[j]) 
	mu += a[j];
    *resul = SUP(*resul, INF(f[index[i]], mu));
  }

  index = NULL;
}

/*****************************************************************************
  
  Calculation of the Sipos integral.
  mu : game (2^n coefficients)
  f : function (n coefficients) 

*****************************************************************************/

void Sipos_integral_game(int *n, double *mu, double *f, double *resul) {

  int i, p;
  int tourn[NMAX];
  int *index = (int *) R_alloc(*n, sizeof(int));

  /* sorts f */
  tri(*n, tourn, f, index);

  /* Computes p such that
     f[index[0]] <= ... <= f[index[p-1]] < 0 <= f[index[p]] <= ... 
     <= f[index[n-1]]
  */
  for(p=0; (p<*n) && (f[index[p]]<0); p++);

  *resul=0.0;

  /* Computes the terms of the integral for i<=p */
  if(p>0) {

    for(i=0; i<p-1; i++)
      *resul += (f[index[i]] - f[index[i+1]]) 
	* mu[subset2binary(index, i+1)];

    *resul += f[index[p-1]] * mu[subset2binary(index, p)];
  }

  /* Computes the terms of the integral for i>p */
  if(p<*n) {

    *resul += f[index[p]] * mu[subset2binary(index+p, *n-p)];

    for(i=p+1; i<*n; i++)
      *resul += (f[index[i]] - f[index[i-1]]) 
	* mu[subset2binary(index+i, *n-i)];
  }
  
  index = NULL;
}

/*****************************************************************************

  Generation of the first k + 1 levels of the power set of X 
  in the "natural" order. Recursive function.

*****************************************************************************/

void k_power_set_rec(int n, int k, int last, int *power_set, int *b) {

  int i, istart;

  /* look for the leftmost 1 in b and start to fill blank cases 
     with 1 left from this position */
  istart = n;

  if(*b > 0)
    while(!(*b & 1<<(istart-1)))
      istart--;
  else
    istart = 0;

  for(i=istart; i<n; i++) {

    last++;
    *(power_set+last) = *b + (1<<i);
  }

  if(last != (int)sum_binom(n,k) - 1)
    k_power_set_rec(n, k, last, power_set, b+1);
    
}

/*****************************************************************************

  Generation of the first k + 1 levels of the power set of X 
  in the "natural" order. Wrapps the previous function. 

*****************************************************************************/

void k_power_set(int *n, int *k, int *power_set) {

  power_set[0] = 0;
  k_power_set_rec(*n, *k, 0, power_set, power_set);
}

/*****************************************************************************

  Converts the k power set of X in the "natural" order to char**.

*****************************************************************************/

void k_power_set_char(int *n, int *k, int *k_power_set, char **subset) {

  int i, j;
  int x[NMAX];
  char string[255];
  
  sprintf(subset[0],"{}");

  for(i=1; i<sum_binom(*n,*k); i++) {

    for(j=0; j<*n; j++)
      x[j]=0;

    binary2subset(*n,k_power_set[i],x);
      
    sprintf(subset[i],"{%d",x[0]+1);

    for(j=1; j<cardinal(k_power_set[i]); j++) {

      sprintf(string,",%d", x[j]+1);
      strcat(subset[i],string);
    }

    strcat(subset[i],"}");
  }
}

/*****************************************************************************

  Converts the power set of X in the binary order to char**.

*****************************************************************************/

void power_set_binary_char(int *n, char **power_set) {

  int i, j;
  int x[NMAX];
  char string[255];

  sprintf(power_set[0],"{}");

  for(i=1; i<(1<<*n); i++) {

    for(j=0; j<*n; j++)
      x[j]=0;

    binary2subset(*n,i,x);
      
    sprintf(power_set[i],"{%d",x[0]+1);

    for(j=1; j<cardinal(i); j++) {

      sprintf(string,",%d", x[j]+1);
      strcat(power_set[i],string);
    }

    strcat(power_set[i],"}");
  }
}

/*****************************************************************************

  Prints the set function in "natural" order in the R terminal

*****************************************************************************/

void Rprint_setfunction(int *n, int *k, double *mu, int *subset, int *mobius) {

  int	i,j;
  int	x[NMAX];

  Rprintf("{}\t\t%lf\n",mu[0]);

  for(i=1; i<sum_binom(*n,*k); i++) {

      for(j=0; j<*n; j++) 
	x[j]=0;
      binary2subset(*n,subset[i],x);

      Rprintf("{%d",x[0]+1);

      for(j=1; j<cardinal(subset[i]); j++)
	Rprintf(",%d",x[j]+1);
      
      if (*mobius)
	Rprintf("}\t\t%lf\n",mu[i]);
      else
	Rprintf("}\t\t%lf\n",mu[subset[i]]);
    }  
}

/*****************************************************************************

  Writing a set function given in binary order in the "natural" order
  Pre: power_set contains the power_set in "natural" order

*****************************************************************************/

void binary2natural(int *n, double *sf, int *power_set, double *sf_out) {

  int i;
  for(i=0; i<(1<<*n); i++)
    sf_out[i] =  sf[power_set[i]];
}

/*****************************************************************************
 
  Writing a set function given in "natural" order in the binary order
  Pre: power_set contains the power_set in "natural" order

*****************************************************************************/

void natural2binary(int *n, double *sf, int *power_set, double *sf_out) {

  int i;
  for(i=0; i<(1<<*n); i++)
    sf_out[power_set[i]] = sf[i];
}


/*****************************************************************************

  Is sf a k-cardinal set function? Used to test both the cardinality of 
  a set function or of its Mobius representation.
  Pre: sf is given in "natural" order.

*****************************************************************************/

void is_kcardinal(int *n, int *k, double *sf, int *flag) {

  int i, j, l;
  int cni;

  *flag = 0;
  l = 1;
  for(i=1; i<=INF(*k,*n-1); i++) {

    cni = (int)binom(*n,i);
    for (j=1; j < cni; j++) {

      if (sf[l] != sf[l+1]) {
 
	*flag = 1;
	return;
      }
      l++;
    }
    l++;
  }
}

/*****************************************************************************

  Is mu a cardinal set function? 
  Pre: mu given in binary order

*****************************************************************************/

void is_cardinal_setfunction(int *n, double *mu, int *power_set, int *flag) {

  double *mu_out = (double *) R_alloc((1<<*n), sizeof(double));
  binary2natural(n, mu, power_set, mu_out);

  is_kcardinal(n, n, mu_out, flag);

  mu = NULL;
}

/*****************************************************************************

  Is mu a k-additive set function?
  Pre: 1 <= k <= n
  kmax: the order of truncation of mu.
  a: Mobius representation of mu, given in the "natural" order 

*****************************************************************************/

void is_kadditive_Mobius(int *n, int *kmax, int *k, double *a, double *epsilon, 
			 int *flag) {

  int i;
  int sb = (int)sum_binom(*n,*k-1);
  int cnk = (int)binom(*n,*k);

  *flag = 1;
  for (i=0; i<cnk; i++)
    if (fabs(a[sb + i]) > *epsilon) {

      *flag = 0;
      break;
    }

  if (*flag == 1)
    return;

  sb += cnk;

  for (i=sb; i<(int)sum_binom(*n,*kmax); i++)
    if (fabs(a[i]) > *epsilon) {
 
      *flag = 1;
      return;
    }
}

/*****************************************************************************

  Is mu a k-additive set function?
  Pre: mu given in binary order 

*****************************************************************************/

void is_kadditive_setfunction(int *n, int *k, double *mu, int *power_set,
			      double *epsilon, int *flag) {

  double *mu_out = (double *) R_alloc((1<<*n), sizeof(double));
  
  setfunction2Mobius(mu,n);
  
  binary2natural(n, mu, power_set, mu_out);

  is_kadditive_Mobius(n, n, k, mu_out, epsilon, flag);

  mu = NULL;
}


/*****************************************************************************

  Adds a strict veto criterion to mu_init (n criteria) at position n. 
  Result in mu (n+1 criteria).

*****************************************************************************/

void add_veto_setfunction(int *n, double *mu_init, double *mu) {

  int i;
  int power2n = 1<<*n;

  for(i=0; i<power2n; i++)
    mu[i+power2n] = mu_init[i];
}

/*****************************************************************************

  Searchs for the upper neighbors 

*****************************************************************************/

void search_upper_neighbors(int n, int node, int *neighbors) {

  int i, count = 0;

  for(i=0; i<n; i++)
    if(!(1<<i & node)) {

      neighbors[count] = node + (1<<i);
      count++;
    }
}

/*****************************************************************************

  Searchs for the lower neighbors 

*****************************************************************************/

void search_lower_neighbors(int n, int node, int *neighbors) {

  int i, count = 0;

  for(i=0; i<n; i++)
    if(1<<i & node) {

      neighbors[count] = node - (1<<i);
      count++;
    }
}

/*****************************************************************************

  Is the set function mu monotone?
  Pre: mu is given in binary order

*****************************************************************************/

void is_monotone_setfunction(int *n, double *mu, int *print, double *epsilon,
			     int *flag) {

  int i, j, k, level, c;
  int neighbors[NMAX];
  int subset[NMAX];
  
  *flag = 0;

  for(i=0; i<(1<<*n)-1; i++) {

    search_upper_neighbors(*n, i, neighbors);
    level = cardinal(i);

    for(j=0; j<*n-level; j++) {

      if(mu[ i ] - mu[neighbors[j]] > *epsilon) {

	*flag = 1;

	if(*print) {

	  Rprintf("Violation of monotonicity constraint between {");
	  binary2subset(*n, i, subset);
	  c = cardinal(i);
		  
	  for(k=0; k<cardinal(i); k++)
	    Rprintf(" %d", subset[k] + 1);

	  Rprintf(" } and {");

	  binary2subset(*n, neighbors[j], subset);

	  for(k=0; k<c+1; k++)
	    Rprintf(" %d", subset[k] + 1);

	  Rprintf(" }. \n");
	}
	else
	  return;
      }
    }
  }
}

/*****************************************************************************

  Does the Mobius representation a correspond to a monotone set function?
  Pre: a given in "natural" order

*****************************************************************************/

void is_monotone_Mobius(int *n, int *k, double *a, int *subset, int *print, 
			double *epsilon, int *flag) {

  int i,j,l,m,c,r;
  int pow = 1<<*n;
  int sb = (int)sum_binom(*n,*k);
  double s; 
  int list[NMAX];

  *flag = 0;
  
  for (i=0; i<*n; i++)
    for (j=1; j<pow; j++)
      if(j & 1<<i) {

	s = 0.0;
	for(l=1;l<sb;l++)
	  if (((subset[l] & j) == subset[l]) && (subset[l] & 1<<i))
	    s += a[l];

	if (s < - *epsilon) {

	  *flag = 1;
	  if (*print) {
 
	    r = j ^ 1<<i;
	    Rprintf("Violation of monotonicity constraint between {");
	    binary2subset(*n, r, list);
	    c = cardinal(r);

	    for(m=0; m<c; m++)
	      Rprintf(" %d", list[m] + 1);

	    Rprintf(" } and {");

	    binary2subset(*n, j, list);

	    for(m=0; m<c+1; m++)
	      Rprintf(" %d", list[m] + 1);

	    Rprintf(" }. \n");
	  }
	  else
	    return;
	}
	    
      
      }
}

/*****************************************************************************

  Factorial n

*****************************************************************************/

double fact(int n) {
  int i;
  double f = 1;

  for(i=n; i>0; i--)
    f *= i;

  return(f);
}

/*****************************************************************************

  The gamma function

*****************************************************************************/

double gamm(int a, int n) {

  return(fact(n - a - 1) * fact(a) / fact(n));
}

/*****************************************************************************

  The zeta function
 
*****************************************************************************/

double zeta(int a, int n) {

  return(fact(n - a - 2) * fact(a) / fact(n - 1));
}

/*****************************************************************************

  Choose(n,k)

*****************************************************************************/

double binom(int n, int k) {

  /* return(fact(n) / (fact(k) * fact(n-k))); */

  if ((k==0) || (k==n))
    return 1.0;
  else
    return binom(n-1,k) + binom(n-1,k-1);

}

/*****************************************************************************

  sum_{i=0}^k Choose(n,i)

*****************************************************************************/

double sum_binom(int n, int k) {

  int i;
  double s = 1.0; 
  for (i=1;i<=k;i++)
    s += binom(n,i);
  return s;
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

  Calculation of the Shapley value of mu. 
  Result in phi (dimension n).

*****************************************************************************/

void Shapley_value_setfunction(int *n, double *mu, double *phi) {

  int i, k;
  int pow = 1<<*n;

  for(k=0; k<*n; k++) {

    phi[k] = 0.0;
    for(i=0; i<pow; i++)
      if(!(i & 1<<k))
	phi[k] += gamm(cardinal(i), *n) * (mu[i + (1<<k)] - mu[i]);
  }
}

/*****************************************************************************

  Calculation of the Shapley value of a set function using its 
  Mobius representation a. 
  Result in phi (dimension n).

*****************************************************************************/

void Shapley_value_Mobius(int *n, int *k, double *a, int *subset, 
			  double *phi) {
  int i,j; 
  int sb = (int)sum_binom(*n,*k);

  for (i=0;i<*n;i++) { 
    phi[i] = 0.0;
    for (j=1; j<sb; j++)
      if (subset[j] & 1<<i)
	phi[i] += a[j] / (double)cardinal(subset[j]);
  }
}

/*****************************************************************************

  Calculation of the Shapley interaction indices of mu. 
  Result in phi (dimension n*n).

*****************************************************************************/

void interaction_indices_setfunction(int *n, double *mu, double *phi) {

  int i, j, k;
  int pow = 1<<*n;

  for(k=0; k<*n; k++) {

    for(j=k+1; j<*n; j++) {

      phi[k*(*n)+j] = 0.0;
      for(i=0; i<pow; i++)
	if(!(i & 1<<k) && !(i & 1<<j))
	  phi[k*(*n)+j] += zeta(cardinal(i), *n) * 
	    (mu[i+(1<<k)+(1<<j)] - mu[i+(1<<k)] - mu[i+(1<<j)]+ mu[i]);
	       
    }
  }

  for(k=0; k<*n; k++)
    for(j=k+1; j<*n;j++)
      phi[j*(*n)+k] = phi[k*(*n)+j];
     

  for(k=0; k<*n; k++)
    phi[k*(*n)+k] = 0.0;
}

/*****************************************************************************

  Calculation of the Shapley interaction indices of a set function using its 
  Mobius representation a. 
  Result in phi (dimension n*n).

*****************************************************************************/

void interaction_indices_Mobius(int *n, int *k, double *a, int *subset, 
				double *phi) {
  int i, j, l;
  int sb = (int)sum_binom(*n,*k);

  for(i=0; i<*n; i++)
    for(j=i+1; j<*n; j++) {

      phi[i*(*n)+j] = 0.0;
      for(l=1; l<sb; l++)
	if((subset[l] & 1<<i) && (subset[l] & 1<<j))
	  phi[i*(*n)+j] += a[l]/(double)(cardinal(subset[l])-1);    
    }

  for(i=0; i<*n; i++)
    for(j=i+1; j<*n;j++)
      phi[j*(*n)+i] = phi[i*(*n)+j];
     

  for(i=0; i<*n; i++)
    phi[i*(*n)+i] = 0.0;

}

/*****************************************************************************

  Calculation of veto indices of a capacity mu 
  Result in v (dimension n)

*****************************************************************************/

void veto_capacity(int *n, double *mu, double *v) {
  int i, k;
  int pow = 1<<*n;

  for(k=0; k<*n; k++) {

    v[k] = 0.0;
    for(i=1; i<pow; i++)
      if(!(i & 1<<k))
	v[k] += mu[i] / binom(*n-1,cardinal(i));
	
    v[k] /= (double)(*n - 1) * mu[pow-1];
    v[k] = 1.0 - v[k];
  }
}

/*****************************************************************************

  Calculation of veto indices from the Mobius representation of a capacity
  Result in v (dimension n)

*****************************************************************************/

void veto_Mobius(int *n, int *k, double *a, int *subset, double *v) {

  int i,j;
  int sb = (int)sum_binom(*n,*k);

  normalize_Mobius(*n,*k,a);

  for(i=0; i<*n; i++) {

    v[i] = 0.0;
    for(j=1; j<sb; j++)
      if(!(subset[j] & 1<<i))
	v[i] += a[j] / (double)(cardinal(subset[j]) + 1);
	
    v[i] *= (double)(*n)/(double)(*n - 1);
    v[i] = 1.0 - v[i];
  }
}

/*****************************************************************************

  Calculation of favor indices of a capacity mu 
  Result in favor (dimension n)

*****************************************************************************/

void favor_capacity(int *n, double *mu, double *f) {

  int i, k;
  int pow = 1<<*n;

  for(k=0; k<*n; k++) {

    f[k] = 0.0;
    for(i=0; i<pow; i++)
      if(!(i & 1<<k))
	f[k] += mu[i + (1<<k)] / binom(*n-1,cardinal(i));
	
    f[k] /= (double)(*n - 1) * mu[pow-1];
    f[k] -= 1.0/(double)(*n - 1);
  }
}

/*****************************************************************************

  Calculation of favor indices from the Mobius representation of a capacity
  Result in f (dimension n)

*****************************************************************************/

void favor_Mobius(int *n, int *k, double *a, int *subset, double *f) {
  int i,j,l;
  int sb = (int)sum_binom(*n,*k); 

  normalize_Mobius(*n,*k,a);

  for(i=0; i<*n; i++) {

    f[i] = 0.0;
    for(j=0; j<sb; j++)
      if(!(subset[j] & 1<<i)) {

	for (l=j+1;l<sb;l++)
	  if ((subset[j] | 1<<i) == subset[l])
	    break;
	if (l < sb)
	  f[i] += (a[j] + a[l]) 
	    / (double)(cardinal(subset[j]) + 1);
	else
	  f[i] += a[j] / (double)(cardinal(subset[j]) + 1);
      }
	
    f[i] *= (double)(*n)/(double)(*n - 1);
    f[i] -= 1.0/(double)(*n - 1);
  }
}

/*****************************************************************************

  Calculation of the orness of a capacity mu 

*****************************************************************************/

void orness_capacity(int *n, double *mu, double *resul) {

  int i; 
  int pow = 1<<*n;

  *resul = 0.0;
  for(i=1; i<pow-1; i++) /* could be optimized */ 
    *resul += mu[i] / binom(*n, cardinal(i)); 

  *resul /= (double)(*n - 1) *  mu[pow-1];
}

/*****************************************************************************

  Calculation of the orness from the Mobius representation of a capacity 

*****************************************************************************/

void orness_Mobius(int *n, int *k, double *a, int *subset, double *resul) {

  int i,c; 
  int sb = (int)sum_binom(*n,*k); 

  normalize_Mobius(*n,*k,a);

  *resul = 0.0;
  for(i=1; i<sb; i++) {

    c = cardinal(subset[i]);
    *resul +=  (double)(*n - c) * a[i] / (double)(c + 1); 
  }

  *resul /= (double)(*n - 1);
}

/*****************************************************************************

  Calculation of the variance of a capacity mu

*****************************************************************************/

void variance_capacity(int *n, double *mu, double *resul) {

  int i, k;
  int pow = 1<<*n;

  *resul = 0.0;

  for(k=0; k<*n; k++)
    for(i=0; i<pow; i++)
      if(!(i & 1<<k))
	*resul += gamm(cardinal(i), *n) 
	  * SQR((mu[i + (1<<k)] - mu[i])
		/ mu[pow-1]); 

  *resul -= 1.0/(double)*n;

  *resul /= 1 - 1.0/(double)*n;
}

/*****************************************************************************

  Calculation of the variance from the Mobius representation of a capacity 

*****************************************************************************/

void variance_Mobius(int *n, int *k, double *a, int *subset, double *resul) {
 
  int i,j,l;
  int pow = 1<<*n;
  int sb = (int)sum_binom(*n,*k);
  double s; 

  normalize_Mobius(*n,*k,a);

  for (i=0; i<*n; i++)
    for (j=1; j<pow; j++)
      if(j & 1<<i) {

	s = 0.0;
	for(l=1;l<sb;l++)
	  if (((subset[l] & j) == subset[l]) 
	      && (subset[l] & 1<<i))
	    s += a[l];
	    
	*resul += gamm(cardinal(j)-1, *n) * SQR(s);
      }

  *resul -= 1.0/(double)*n;

  *resul /= 1 - 1.0/(double)*n;
}

/*****************************************************************************

  Calculation of the normalized Marichal entropy of a capacity mu 

*****************************************************************************/

void entropy_capacity(int *n, double *mu, double *resul) {

  int i, k;
  int pow = 1<<*n;

  *resul = 0.0;

  for(k=0; k<*n; k++)
    for(i=0; i<pow; i++)
      if(!(i & 1<<k))
	*resul += gamm(cardinal(i), *n) 
	  * H((mu[i + (1<<k)] - mu[i])
	      /mu[pow-1]);

  *resul /= log(*n); 
}

/*****************************************************************************

  Calculation of the normalized Marichal 
  from the Mobius representation of a capacity 

*****************************************************************************/

void entropy_Mobius(int *n, int *k, double *a, int *subset, double *resul) {
 
  int i,j,l;
  int pow = 1<<*n;
  int sb = (int)sum_binom(*n,*k);
  double s; 

  normalize_Mobius(*n,*k,a);
  
  for (i=0; i<*n; i++)
    for (j=1; j<pow; j++)
      if(j & 1<<i) {

	s = 0.0;
	for(l=1;l<sb;l++)
	  if (((subset[l] & j) == subset[l]) 
	      && (subset[l] & 1<<i))
	    s += a[l];
	    
	*resul += gamm(cardinal(j)-1, *n) * H(s);
      }

  *resul /= log(*n);
}

/****************************************************************************/
