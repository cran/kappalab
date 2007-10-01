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

*****************************************************************************/

#ifndef KAPPALAB_CONSTRAINTS_H
#define KAPPALAB_CONSTRAINTS_H

void monotonicity_constraints(int *n, int *k, int *subset, int *A);
void Choquet_preorder_constraint(int *n, int *k, int *subset, double *a, 
				 double *b, double *A);
void Shapley_preorder_constraint(int *n, int *k, int *subset, int *a, int *b, 
			double *A);
void Shapley_constant_constraint(int *n, int *k, int *subset, int *a, double *A);
void interaction_preorder_constraint(int *n, int *k, int *subset, int *a, 
				     int *b, int *c, int *d, double *A);
void interaction_constant_constraint(int *n, int *k, int *subset, int *a, 
				     int *b, double *A);
void inter_additive_constraint(int *n, int *k, int *subset, int *partition, 
			       int *p, double *A);

#endif /* ! KAPPALAB_CONSTRAINTS_H */
