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

  Michel Grabisch, 1995
  Contributions: Ivan Kojadinovic, 2005

*****************************************************************************/

#include <R.h>
#include "core.h"
#include "least.squares.h"


/*****************************************************************************

      Calcul de la cardinalite de l'ensemble difference entre deux
      ensembles I et J donnes par ses codes binaires.
      Il faut que I soit contenu dans J.

 ******************************************************************************/

int difference (int i,int j,int n)
{
 int k;
 int l,dif;
 
 /*  Si i=b000001100 et j=b101001101 (on a bien i inclu dans j),
     alors ~i=b111110011 et donc j&(~i)=b001000001.
     Donc k=j&(~i) est l'ensemble difference J\I.
     Attention: Dans ~i, l'operateur "~" est un TILDE et non un "-".  */
 k=j&(~i);
 
 dif=0;
 for (l=0;l<n;l++)
   if (k&(1<<l))
     dif++;
 return (dif);
}

/************************************************************************

      Lower bound for Mobius, Shapley et Banzhaf
      i: ensemble

 ************************************************************************/
int lower_bound (int i,int n)
{
  int l,p,diff;
  int bound;

  diff=difference(0,i,n);
  p=diff%4;
  if (p==0) l=diff/2-1;
  if (p==1) l=(diff+1)/2;
  if (p==2) l=diff/2;
  if (p==3) l=(diff-1)/2;
  /*  l est impair */
  bound=-binom(diff,l);
  return(bound);
}

/************************************************************************

      Upper bound for Mobius,Shapley et Banzhaf
      i :ensemble

 ************************************************************************/ 
int upper_bound (int i,int n)
{
  int l,p,diff;
  int bound;
 
  diff=difference(0,i,n);
  p=diff%4;
  if (p==0) l=diff/2;
  if (p==1) l=(diff-1)/2;
  if (p==2) l=diff/2-1;
  if (p==3) l=(diff+1)/2;
  /*  l est pair */
  bound=binom(diff,l);
  return(bound);
}

/****************************************************************************

    Si i s'ecrit en ordre binaire : i=b10110101, alors f[maximal(p,i,index)]
    doit valoir $M:=max_{j\in {0,2,4,5,7}\cap{0,...,p} } f[j]$. L'ordre
    de $f$ est donn\'e par "index":
        f[index[0]] \leq f[index[1]] \leq ... \leq f[index[n-1]] .
    On regarde donc si index[p]\in I:={0,2,4,5,7}. Si c'est le cas, c'est
    que $M=f[index[p]]$. Si ce n'est pas le cas, on regarde si
    index[p-1]\in I:={0,2,4,5,7}. Et ainsi de suite.

****************************************************************************/

int maximal(int p,int i,int *index)
{
  int j;
  
  j=p;
  while (!(i&(1<<index[j])))
    j--;
  return index[j];
}

/****************************************************************************

    Si i s'ecrit en ordre binaire : i=b10110101, alors f[minimal(p,i,index)]
    doit valoir $M:=min_{j\in {0,2,4,5,7}\cap{p,...,n-1} } f[j]$. L'ordre
    de $f$ est donn\'e par "index":
        f[index[0]] \leq f[index[1]] \leq ... \leq f[index[n-1]] .
    On regarde donc si index[p]\in I:={0,2,4,5,7}. Si c'est le cas, c'est
    que $M=f[index[p]]$. Si ce n'est pas le cas, on regarde si
    index[p+1]\in I:={0,2,4,5,7}. Et ainsi de suite.

****************************************************************************/

int minimal(int p,int i,int *index)
{
  int j;
  
  j=p;
  while (!(i&(1<<index[j])))
    j++;
  return index[j];
}

/********************************************************************

    Borne inf. de la representation de Möbius
    Utilisé dans la cadre de la prog. lin.


********************************************************************/

void Mobius_lower_bound(int *n, int *k, int *subset, double *lower)
{
  int i;
  int sb = (int)sum_binom(*n,*k);

   for (i=0;i<sb-1;i++) 
    lower[i]=(double)lower_bound(i+1, *n);
}

/*****************************************************************************/
