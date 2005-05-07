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

  Non-additive measure and integral manipulation routines.
  C functions source file. Binary coding of sets.
  Univers X, X={0,1,2,3,...,n-1} with n=32 maximum.

*****************************************************************************/

#ifndef KAPPALAB_CORE_H
#define KAPPALAB_CORE_H

/*****************************************************************************

	Parameters

*****************************************************************************/

#define		NMAX		32

/*****************************************************************************

	Macros 

*****************************************************************************/

#define	INF(x,y) ((x) > (y) ? (y) : (x))
#define	SUP(x,y) ((x) < (y) ? (y) : (x))
#define	ABS(x) SUP((x),-(x))
#define	SQR(x)	((x) * (x))
#define H(x) ((x) > 0.0 ? (-(x) * log(x)) : (0.0)) 

/*****************************************************************************

  	Non interfaced functions

*****************************************************************************/

int subset2binary(int *x, int nx);
void binary2subset(int n, int b, int *x);
int cardinal(int n);

void tri(register int n, register int *tourn, register double *vec, 
	 register int *iv);

void k_power_set_rec(int n, int k, int last, int *power_set, 
			       int *b);
void search_upper_neighbors(int n, int node, int *neighbors);
void search_lower_neighbors(int n, int node, int *neighbors);

double fact(int n);
double binom(int n, int k);
double sum_binom(int n, int k);
double gamm(int a, int n);
double zeta(int a, int n);

void fast_lower_cardinality_transform(double *mu, double c, int n);

void normalize_Mobius(int n, int k, double *a);

/*****************************************************************************

  	Interfaced functions

*****************************************************************************/

void binary2subsetR(int *n, int *b, int *x, int *l);

void setfunction2Mobius (double *mu, int *n);
void Mobius2setfunction(int *n, int *k, double *a, int *subset, double *mu);
void k_truncation(int *n, int *k, double *mu, int *subset, double *a);
void cardinal2setfunction(int *n, double *csf, double *sf);
void setfunction2cardinal(int *n, int *k, double *sf, double *csf);
void setfunction2conjugate(double *mu, int *n, double *mu_out);

void Choquet_integral_game(int *n, double *mu, double *f, double *resul);
void Choquet_integral_Mobius(int *n, int *k, double *a, int *subset, double *f, 
			     double *resul);
void Sugeno_integral_game(int *n, double *mu, double *f, double *resul);
void Sugeno_integral_Mobius(int *n, int *k, double *a, int *subset, double *f, 
			    double *resul);
void Sipos_integral_game(int *n, double *mu, double *f, double *resul);

void k_power_set(int *n, int *k, int *power_set);
void k_power_set_char(int *n, int *k, int *k_power_set, char **subset);
void power_set_binary_char(int *n, char **power_set);
void binary2natural(int *n, double *sf, int *power_set, double *sf_out);
void natural2binary(int *n, double *sf, int *power_set, double *sf_out);

void Rprint_setfunction(int *n, int *k, double *mu, int *subset, int *mobius);

void is_monotone_setfunction(int *n, double *mu, int *print, double *epsilon, 
			     int *flag);
void is_monotone_Mobius(int *n, int *k, double *a, int *subset, int *print, 
			double *epsilon, int *flag);
void is_kcardinal(int *n, int *k, double *sf, int *flag);
void is_cardinal_setfunction(int *n, double *mu, int *power_set, int *flag);
void is_kadditive_setfunction(int *n, int *k, double *mu, int *power_set,
			      double *epsilon, int *flag);
void is_kadditive_Mobius(int *n, int *kmax, int *k, double *a, 
			 double *epsilon, int *flag);

void add_veto_setfunction(int *n, double *mu_init, double *mu);

void Shapley_value_setfunction(int *n, double *mu, double *phi);
void Shapley_value_Mobius(int *n, int *k, double *a, int *subset, double *phi);
void interaction_indices_setfunction(int *n, double *mu, double *phi);
void interaction_indices_Mobius(int *n, int *k, double *a, int *subset, 
				double *phi);
void veto_capacity(int *n, double *mu, double *v);
void veto_Mobius(int *n, int *k, double *a, int *subset, double *v);
void favor_capacity(int *n, double *mu, double *f);
void favor_Mobius(int *n, int *k, double *a, int *subset, double *v);
void orness_capacity(int *n, double *mu, double *resul);
void orness_Mobius(int *n, int *k, double *a, int *subset, double *resul);
void variance_capacity(int *n, double *mu, double *resul);
void variance_Mobius(int *n, int *k, double *a, int *subset, double *resul);
void entropy_capacity(int *n, double *mu, double *resul);
void entropy_Mobius(int *n, int *k, double *a, int *subset, double *resul);

#endif /* ! KAPPALAB_CORE_H */
