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

  Minimum least squares capacity identification
  Michel Grabisch, 1995
  Contributions: Ivan Kojadinovic, 2005

*****************************************************************************/

#include <R.h>
#include "Mobius.bounds.h"
#include "core.h"
#include "least.squares.h"

/********************************************************************

       Calcul de la matrice objectif
       R: matrice objectif est R'R 
       X: matrice de scores sur les elements en une ligne
       xl: borne inf. des Möbius
       xu: borne sup. des Möbius
       Integral: 1: Choquet, autre: Sipos

       Fonction objectif: x'R'Rx - 2y'Rx

********************************************************************/

void k_additive_objectif(int *n, int *k, int *subset, int *Integral, double *X,
			 int *nX, double *R, double *xl, double *xu)
{
  int i,j,l,m,p,q;
  int sb = (int)sum_binom(*n,*k);
  double min_x, min_xplus, min_xminus;

  p = 0;
  q = 0;

  /* On parcourt les "lignes" de X */ 
  for (m=0;m<*nX;m++) {

    if (*Integral==1) {
      /* Choquet */
      for (i=1;i<sb;i++) {
	
	for (j=0;j<*n;j++)
	  if (subset[i] & 1<<j) {
	    
	    min_x = X[j+p];
	    break;
	  }
	for (l=j+1;l<*n;l++)
	  if (subset[i] & 1<<l) 
		min_x = INF(min_x, X[l+p]);
	
	R[i - 1 + q] = min_x;
      }
    }
    else {
      /* Sipos */
      for (i=1;i<sb;i++) {
	
	for (j=0;j<*n;j++)
	  if (subset[i] & 1<<j) {
	    
	    min_xplus = SUP(X[j+p],0);
	    min_xminus = SUP(-X[j+p],0);
	    break;
	  }
	for (l=j+1;l<*n;l++)
	  if (subset[i] & 1<<l) {

		min_xplus = INF(min_xplus, SUP(X[l+p],0));
		min_xminus = INF(min_xminus, SUP(-X[l+p],0));
	  }
	
	R[i - 1 + q] = min_xplus - min_xminus;
      }
    }
    p += *n;
    q += sb - 1;
  }

  /*  Initialise la borne inf. et sup. de la representation de Möbius */
  for (i=0;i<sb-1;i++) {

    xl[i]=(double)lower_bound(i+1, *n);
    xu[i]=(double)upper_bound(i+1, *n);
  }
  

}

/*****************************************************************************/
