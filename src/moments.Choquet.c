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

  Computation of the moments of the Choquet integral 
  
  Ivan Kojadinovic, 02/2007

*****************************************************************************/

#include <R.h>
#include <Rmath.h>
#include "core.h"


/*****************************************************************************

  First derivative of G in the standard normal case
  G is the inverse of the c.d.f.

*****************************************************************************/

double G1(double phiGp) {

  return 1.0/phiGp;
}

/*****************************************************************************

  Second derivative of G in the standard normal case

*****************************************************************************/

double G2(double Gp, double phiGp) {

  return Gp/(phiGp * phiGp);
}
 
/*****************************************************************************

  Third derivative of G in the standard normal case

*****************************************************************************/

double G3(double Gp, double phiGp) {

  return (1.0 + 2.0 * Gp * Gp)/(phiGp * phiGp * phiGp);
}


/*****************************************************************************

  Fourth derivative of G in the standard normal case

*****************************************************************************/

double G4(double Gp, double phiGp) {

  return Gp * (7.0 + 6.0 * Gp * Gp)/(phiGp * phiGp * phiGp * phiGp);
}

/*****************************************************************************

  Fifth derivative of G in the standard normal case

*****************************************************************************/

double G5(double Gp, double phiGp) {

  return ( 7.0 + 46.0 * Gp * Gp + 24.0 * Gp * Gp * Gp * Gp) 
    / (phiGp * phiGp * phiGp * phiGp * phiGp);
}

/*****************************************************************************

  Sixth derivative of G in the standard normal case

*****************************************************************************/

double G6(double Gp, double phiGp) {

  return ( 127.0 * Gp + 326.0 * Gp * Gp * Gp + 96.0 * Gp * Gp * Gp * Gp * Gp) 
    / (phiGp * phiGp * phiGp * phiGp * phiGp * phiGp);
}

/*****************************************************************************

  Approximate expectation of Xi:n in the standard normal case

*****************************************************************************/

double expectation_Xin(double i, double n) {

  double pi = i/(n+1.0);
  double qi = 1.0 - pi;
  double Gpi = qnorm(pi, 0.0, 1.0, 1, 0);
  double phiGpi = dnorm(Gpi, 0.0, 1.0, 0);

  return Gpi + pi * qi / (2.0 * (n + 2.0)) * G2(Gpi, phiGpi)
    + pi * qi  / SQR(n + 2.0) * ((qi - pi) / 3.0 * G3(Gpi, phiGpi)  + pi * qi / 8.0 * G4(Gpi, phiGpi)) 
    + pi * qi  / (SQR(n + 2.0) * (n + 2.0) ) 
    * ( - (qi - pi) / 3.0 * G3(Gpi, phiGpi) 
	+ (SQR(qi - pi) - pi * qi) / 4.0 * G4(Gpi, phiGpi) 
	+ pi * qi * (qi - pi) / 6.0 * G5(Gpi, phiGpi) 
	+ SQR(pi) * SQR(qi)/48.0 * G6(Gpi, phiGpi));

}

/*****************************************************************************

  Approximate covariance of Xi:n and  Xj:n in the standard normal case

*****************************************************************************/

double covariance_XinXjn(double i, double j, double n) {

  double pi = i/(n+1.0);
  double qi = 1.0 - pi;
  double pj = j/(n+1.0);
  double qj = 1.0 - pj;
  double Gpi = qnorm(pi, 0.0, 1.0, 1, 0);
  double Gpj = qnorm(pj, 0.0, 1.0, 1, 0);
  double phiGpi = dnorm(Gpi, 0.0, 1.0, 0);
  double phiGpj = dnorm(Gpj, 0.0, 1.0, 0);

  if (i > j) /* Probleme if i > j. Why? Don't know!!!! */
    return covariance_XinXjn(j,i,n);
  else
    return pi * qj / (n+2) * G1(phiGpi) * G1(phiGpj)
      + pi * qj / SQR(n+2) * ( (qi - pi) * G2(Gpi, phiGpi) * G1(phiGpj)
			       + (qj - pj) * G1(phiGpi) * G2(Gpj, phiGpj)
			       + pi * qi / 2.0 * G3(Gpi, phiGpi) * G1(phiGpj)
			       + pj * qj / 2.0 * G1(phiGpi) * G3(Gpj, phiGpj)
			       + pi * qj / 2.0 * G2(Gpi, phiGpi) * G2(Gpj, phiGpj) )
      + pi * qj / ( SQR(n+2) * (n+2) )  
      * ( - (qi - pi) * G2(Gpi, phiGpi) * G1(phiGpj)
	  - (qj - pj) *  G1(phiGpi) * G2(Gpj, phiGpj) 
	  + ( SQR(qi - pi) - pi * qi ) * G3(Gpi, phiGpi) * G1(phiGpj)
	  + ( SQR(qj - pj) - pj * qj) * G1(phiGpi) * G3(Gpj, phiGpj)
	  + ( 3.0 / 2.0 * (qi - pi) * (qj - pj) + pj * qi / 2.0 - 2.0 * pi * qj) * G2(Gpi, phiGpi) * G2(Gpj, phiGpj)
	  + 5.0 / 6.0 * pi * qi * (qi - pi) * G4(Gpi, phiGpi) * G1(phiGpj)
	  + 5.0 / 6.0 * pj * qj * (qj - pj) * G1(phiGpi) * G4(Gpj, phiGpj)
	  + ( pi * qj * (qi - pi) + pi * qi / 2.0 * (qj - pj)) * G3(Gpi, phiGpi) * G2(Gpj, phiGpj)
	  + ( pi * qj * (qj - pj) + pj * qj / 2.0 * (qi - pi)) * G2(Gpi, phiGpi) * G3(Gpj, phiGpj)
	  + SQR(pi) * SQR(qi) / 8.0 *  G5(Gpi, phiGpi) * G1(phiGpj)
	  + SQR(pj) * SQR(qj) / 8.0 * G1(phiGpi) * G5(Gpj, phiGpj)
	  + SQR(pi) * qi * qj / 4.0 * G4(Gpi, phiGpi) * G2(Gpj, phiGpj)
	  + pi * pj * SQR(qj) / 4.0 * G2(Gpi, phiGpi) * G4(Gpj, phiGpj)
	  + ( 2 * SQR(pi) * SQR(qj) + 3 * pi * pj * qi * qj) / 12.0 * G3(Gpi, phiGpi) * G3(Gpj, phiGpj));

}

/*****************************************************************************

  Approximate product-moment of Xi:n and  Xj:n in the standard normal case

*****************************************************************************/

double product_moment_XinXjn(double i, double j, double n) {
 
 return covariance_XinXjn(i,j,n) + expectation_Xin(i,n) * expectation_Xin(j,n);
}

/*****************************************************************************

  Approximate expectation of the Choquet integral in the standard normal case

*****************************************************************************/

void expectation_Choquet_norm_game(int *n, double *mu, double *E) {
  
  int i, k, c;
  int pow = 1<<*n;

  *E = 0.0;
  for(k=0; k<*n; k++) {

    for(i=0; i<pow; i++)
      if(!(i & 1<<k)) {

	c = cardinal(i);
	*E += gamm(c, *n) * (mu[i + (1<<k)] - mu[i]) 
	  * expectation_Xin((double)(*n -c),(double)(*n)); 
      }
  }
}

/*****************************************************************************

  Approximate expectation of the Choquet integral in the standard normal case

*****************************************************************************/

void expectation_Choquet_norm_Mobius(int *n, int *k, double *a, int *subset, double *E) {

  int i, sb = (int)sum_binom(*n,*k); 

  *E = 0.0;
  for(i=1; i<sb; i++) {

    *E +=  a[i] *  expectation_Xin(1.0,(double)cardinal(subset[i]));
  }
}

/*****************************************************************************

  Approximate standard deviation of the Choquet integral in the standard normal case
  Very suboptimal!!!

*****************************************************************************/

void sd_Choquet_norm(int *n, double *mu, double *SD) {

  int i,j,k;
  double f = fact(*n), EC = 0.0, EC2 = 0.0;
  int *id = (int *)R_alloc(*n,sizeof(int));
  int *sigma = (int *)R_alloc(*n,sizeof(int));
  double *p = (double *)R_alloc(*n,sizeof(double));

  for (i=0;i<(int)f;i++) {
 
    for (j=0;j<*n;j++)
      id[j]=j;

    lex_permut(*n, i, id, sigma);

    for(j=0; j<*n; j++) 
      p[j] = mu[subset2binary(sigma + j, *n - j)] - mu[subset2binary(sigma + j + 1, *n - j - 1)];
    
    for(j=0; j<*n; j++) {
    
      EC += p[j] * expectation_Xin((double)(j+1), (double)(*n)); 

      for(k=0; k<*n; k++) 
	EC2 += p[j] * p[k] *  product_moment_XinXjn((double)(j+1), (double)(k+1), (double)(*n));  
    }
  }
  
  *SD = sqrt(EC2/f - SQR(EC/f));
}

/*****************************************************************************

  Expectation of the Choquet integral in the uniform case

*****************************************************************************/

void expectation_Choquet_unif(int *n, double *mu, double *E) {

  double EC = 0.0;
  int i;
  for(i=0; i<(1<<*n); i++)
    EC +=  mu[i]/binom(*n,cardinal(i));

  *E = EC/(double)(*n+1) ;
}

/*****************************************************************************

  Standard deviation of the Choquet integral in the uniform case

*****************************************************************************/

void sd_Choquet_unif(int *n, double *mu, double *SD) {

  double EC = 0.0, EC2 = 0.0, s;
  int i,j,pow=1<<*n;

  for(i=0; i<pow; i++) {

    EC +=  mu[i]/binom(*n,cardinal(i));
    
    s = 0.0;
    for(j=0;j<pow;j++)
      if ((i & j) == j) 
	s += mu[j]/binom(cardinal(i),cardinal(j));
    
    EC2 +=  mu[i]/binom(*n,cardinal(i)) * s;
  }

  *SD = sqrt(EC2 * 2.0/((double)((*n +1)*(*n +2))) - SQR(EC/(double)(*n+1)));
}
