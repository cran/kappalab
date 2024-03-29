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

  Computation of divided differences of the minus function
  Ivan Kojadinovic, 04/2006

*****************************************************************************/

#ifndef KAPPALAB_CDFCHOQUET_H
#define KAPPALAB_CDFCHOQUET_H

void cdf_Choquet_unif(int *n, double *mu, double *y, double *Fy);
void pdf_Choquet_unif(int *n, double *mu, double *y, double *py);
void pdf_Choquet_exp(int *n, double *mu, double *y, double *py);

#endif /* ! KAPPALAB_CDFCHOQUET_H */
