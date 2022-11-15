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

#include <R.h>
#include "core.h"

/******************************************************************************

  Calculation of the Choquet integral.
  mu : game (2^n coefficients)
  f : function (n coefficients)

*****************************************************************************/

double Choquet_integral(int n, double *mu, double *f, int *index) {

  int i;
  double resul;
  int tourn[NMAX];

  tri(n, tourn, f, index);
  resul = f[index[0]] * mu[(1<<n)-1];

  for(i=1; i<n; i++)
    resul += (f[index[i]] - f[index[i - 1]]) 
      * mu[subset2binary(index + i, n - i)];

  return resul;
}

/*****************************************************************************
  
  Calculation of the Sipos integral.
  mu : game (2^n coefficients)
  f : function (n coefficients) 

*****************************************************************************/

double Sipos_integral(int n, double *mu, double *f, int *index) {

  int i, p;
  double resul;
  int tourn[NMAX];

  /* sorts f */
  tri(n, tourn, f, index);

  /* Computes p such that
     f[index[0]] <= ... <= f[index[p-1]] < 0 <= f[index[p]] <= ... 
     <= f[index[n-1]]
  */
  for(p=0; (p<n) && (f[index[p]]<0); p++);

  resul=0.0;

  /* Computes the terms of the integral for i<=p */
  if(p>0) {

    for(i=0; i<p-1; i++)
      resul += (f[index[i]] - f[index[i+1]]) 
	* mu[subset2binary(index, i+1)];

    resul += f[index[p-1]] * mu[subset2binary(index, p)];
  }

  /* Computes the terms of the integral for i>p */
  if(p<n) {

    resul += f[index[p]] * mu[subset2binary(index+p, n-p)];

    for(i=p+1; i<n; i++)
      resul += (f[index[i]] - f[index[i-1]]) 
	* mu[subset2binary(index+i, n-i)];
  }
  
  return resul;
}

/******************************************************************************

  Apprentissage de la mesure floue mu a 2**n coefficients par une methode
  heuristique iterative (modification de la precedente)
  mu doit etre initialisee en dehors de la routine, typiquement par une 
  mesure equidistribuee additive.
  L'algorithme utilise le codage binaire des mesures floues

  Resume de la modification :
  step 1.1, 1.2 : idem (update)
  step 1.3 : test avec TOUS les voisins inf ou sup (modifies ou pas)
  on effectue les steps 1.2 et 1.3 sur toutes les donnees (1 iteration)
  step 2.1 : inutile car par construction la mesure est monotone
  step 2.2 : idem avant
  retour au step 1.2 pour une nouvelle iteration	

  Michel GRABISCH, Univ. of Paris I, 2004

  n             nombre d'elements de X; dimension des vecteurs d'apprentissage 
  Integral      l'agregation se fait par l'integrale de Choquet si Integrale==1 
                et par l'integrale de Sipos sinon
                           si Integrale=choiceSipos.  
  itmax        nombre max d'iterations 
  mu           mesure floue (2**n coefficients) 
  n_data       nombre de vecteurs d'apprentissage 
  vectors      vecteurs d'apprentissage:  vectors[0]...vectors[n-1] = 1er vecteur
	       vectors[n]...vectors[2n-1] = 2eme vecteur
	       ...
	       macro : VECT(i,j) = jieme composante du ieme vecteur
			   
  out	       les p reponses desirees aux p vecteurs d'entree 
  alpha	       coeff. dans [0,1] pour step 1. Correction maximale pour alpha = 1 
  epsilon
  error2      somme des erreurs carrees (critere) 

******************************************************************************/

#define	VECT(i,j)	vectors[ (i)*(*n) + (j) ]

void hlms(int *n, int *Integral, int *itmax, double *mu, int *n_data, 
	   double *vectors, double *out, double *alpha, /* double *beta, */ 
	  double *epsilon, double *error2)
{
  int i, j, jj, iter,l,p;
  int *current_node_plus = (int *) R_alloc(*n, sizeof(int) );
  int *current_node_minus = (int *) R_alloc(*n, sizeof(int) );
  
  int *index = (int *) R_alloc(*n, sizeof(int) );
  double error, criterion = 0.0, prec_criterion;
  double *f;
  double *path_mu_plus = (double *) R_alloc( (*n+1), sizeof(double) );
  double *path_mu_minus = (double *) R_alloc( (*n+1), sizeof(double) );
  int *touched_mu = (int *) R_alloc( (1<<*n), sizeof(int) );  
                       /* touched_mu[i] indique si le coefficient 
                       mu[i] a ete modifie (=1) ou non (=0) */

  int *upper_neighbors = (int *) R_alloc(*n, sizeof(int) );
  int *lower_neighbors = (int *) R_alloc(*n, sizeof(int) );

  double min_upper_dist, min_lower_dist;
  double mean_upper_mu = 0.0, mean_lower_mu = 0.0;
  double fading = 1.0, grad,emax;
  double outmax, res, norm,normax;


  /*  On met touched_mu a 1 pour l'ensemble vide et X afin de ne pas
      modifier mu(0) et mu(X) dans STEP 2.  */
  touched_mu[0] = 1;
  touched_mu[ (1<<*n) - 1 ] = 1;
  for( i=1; i<(1<<*n)-1; i++ )
    touched_mu[i] = 0;

  /* Calcule le module de out */
  outmax=ABS(out[0]);
  for( i=1; i<*n_data; i++ )
    outmax = SUP( outmax , ABS(out[i]) );
  outmax = 1./outmax;
  
  if (*Integral==1)
    {
      /*         Cas de l'integrale de CHOQUET  */

      iter = 0;
      while( (iter<2) || ( (iter<*itmax) && (criterion - prec_criterion < 0.)
			   && (normax>*epsilon) )  )
	{      

	  /* STEP 1.2 : apply the gradient */

	  prec_criterion = criterion;
	  criterion = 0.0;
	  normax=0.0;
	  /* printf(" iter = %d : ", iter); */
	  for( i=0; i<*n_data; i++ )
	    {
	      /* current vector */
	      f = vectors+(i*(*n));
	      /* calcul de l'erreur de modele */
	      res = Choquet_integral(*n, mu, f, index);
	      
	      /*  Cacul de l'erreur L2 */
	      error = res - out[i];
	      criterion += SQR(error);
	      /*  Calcul de l'erreur max entre le modele et la mesure */
	      norm=outmax*ABS(error);
	      normax = SUP(normax,norm);
	      /* mise en memoire des coefficients de mu utilises */
	      path_mu_plus[0] = 0;
	      path_mu_plus[*n] = 1;
	      for( j=1; j<*n; j++ )
		{
		  current_node_plus[j] = subset2binary(index+(*n)-j, j);
		  path_mu_plus[j] =  mu[current_node_plus[j]];
		}
	      /*  Calcul la valeur de emax */
	      emax=0.0;
	      for( j=1; j<*n; j++ )
		emax += SQR( f[index[*n-j]]-f[index[*n-j-1]] );
	      /* modification des coefficients utilises dans le calcul precedent */
	      if( error > 0. )
		/* correction vers diminution : on commence par modifier 
		   path_mu_plus[1], path_mu_plus[2], ... */
		for( j=1; j<*n; j++)
		  {
		    if( f[index[*n-j]] - f[index[*n-j-1]] != 0. )
		      {
			path_mu_plus[j] -= *alpha * error / emax * ( f[index[*n-j]] -
								     f[index[*n-j-1]] );
			touched_mu[current_node_plus[j]] = 1;
              
			/* STEP 1.3 : verification des contraintes de monotonie 
			   avec les voisins inferieurs (il y en a j)*/
			search_lower_neighbors(*n, current_node_plus[j], 
					       lower_neighbors);
			for( l=0; l<j; l++ )
			  {
			    if( mu[ lower_neighbors[l] ] > path_mu_plus[ j ] )
			      path_mu_plus[ j ] = mu[ lower_neighbors[l] ];
			  }
			/* update */
			mu[current_node_plus[j]] = path_mu_plus[j];
		      }
		  }
	      else if( error < 0. )
		/* correction vers augmentation : on commence par modifier 
		   path_mu_plus[*n-1], path_mu_plus[*n-2], ... */
		for( j=*n-1; j>0; j--)
		  {
		    if( f[index[*n-j]] - f[index[*n-j-1]] != 0. )
		      {
			path_mu_plus[j] -= *alpha * error / emax * ( f[index[*n-j]] -
								     f[index[*n-j-1]] );
			touched_mu[current_node_plus[j]] = 1;
			/* STEP 1.3 : verification des contraintes de monotonie 
			   avec les voisins superieurs (il y en a *n-j)*/
			search_upper_neighbors(*n, current_node_plus[j], 
					       upper_neighbors);
			for( l=0; l<*n-j; l++ )
			  {
			    if( mu[ upper_neighbors[l] ] < path_mu_plus[ j ] )
			      path_mu_plus[ j ] = mu[ upper_neighbors[l] ];
			  }
			/* update */
			mu[current_node_plus[j]] = path_mu_plus[j];
		      }
		  }
        
	    } /* end of for(i=0 ... */
         
	  *alpha *= fading;

	  /* STEP 2 : equilibrate the fuzzy measure for untouched coefficients */
	  /* lower levels first */
	  for( jj=1; jj<*n; jj++ ) {
	    for ( i=1; i<(1<<*n)-1; i++ )
	      if( touched_mu[i] == 0 )
		if( cardinal( i ) == jj )
		  {
		    /* look for upper neighbors */
		    search_upper_neighbors(*n, i, upper_neighbors);
		    /* look for lower neighbors */
		    search_lower_neighbors(*n, i, lower_neighbors);
		    /* calcul de diverses constantes */
		    mean_upper_mu = 0;
		    for( l=0; l<*n-jj; l++)
		      mean_upper_mu += mu[ upper_neighbors[l] ];
		    mean_upper_mu /= *n-jj;
		    min_upper_dist = mu[ upper_neighbors[0] ] - mu[i];
		    for( l=1; l<*n-jj; l++ )
		      if( mu[ upper_neighbors[l] ] - mu[i] < min_upper_dist )
			min_upper_dist = mu[ upper_neighbors[l] ] - mu[i];
	      
		    mean_lower_mu = 0;
		    for( l=0; l<jj; l++)
		      mean_lower_mu += mu[ lower_neighbors[l] ];
		    mean_lower_mu /= jj;
		    min_lower_dist = mu[i] - mu[ lower_neighbors[0] ];
		    for( l=1; l<jj; l++ )
		      if( mu[i] - mu[ lower_neighbors[l] ] < min_lower_dist )
			min_lower_dist = mu[i] - mu[ lower_neighbors[l] ];
	    
		    /* update */
		    /*if( mean_upper_mu + mean_lower_mu - 2*mu[i] > 0. )*/
		      /* correction vers augmentation */
		      /*mu[i] += *beta * ( mean_upper_mu + mean_lower_mu - 2*mu[i] ) /
			( mean_upper_mu - mean_lower_mu ) * min_upper_dist / 2.;*/
		    /*else*/
		      /* correction vers diminution */
		      /*mu[i] += *beta * ( mean_upper_mu + mean_lower_mu - 2*mu[i] ) /
			( mean_upper_mu - mean_lower_mu ) * min_lower_dist / 2.;*/
		    mu[i] = mu[i] + (min_upper_dist - min_lower_dist)/2.;
		  }      
	  }	
	  /* *beta *= fading; */
	  iter++;
	}   /* end of while(iter... */
    }
  else
    {
      /* Cas de l'integrale de SIPOS */

      iter = 0;
      while( (iter<2) || ( (iter<*itmax) && (criterion - prec_criterion < 0.)
			   && (normax>*epsilon) )  )
	{ 

	  /* STEP 1.2 : apply the gradient */

	  prec_criterion = criterion;
	  criterion = 0.0;
	  normax=0.0;
	  /* printf(" iter = %d : ", iter); */
	  for( i=0; i<*n_data; i++ )
	    {
	      /* current vector */
	      f = vectors+(i*(*n));
	      /* calcul de l'erreur de modele */
	      res = Sipos_integral(*n, mu, f, index);
	      /*  Cacul de l'erreur L2 */
	      error = res - out[i];
	      criterion += SQR(error);
	      /*  Calcul de l'erreur max entre le modele et la mesure */
	      norm=outmax*ABS(error);
	      normax = SUP(normax,norm);
	      /*  Calcule l'indice p tel que
		  f[index[0]] <= ... <= f[index[p-1]] < 0 <= f[index[p]] <= ... 
		  f[index[n-1]] */
	      for (p=0; (p<*n) && (f[index[p]]<0); p++);
	      /* mise en memoire des coefficients de mu utilises */
	      /* Cas des valeurs positives */
	      path_mu_plus[0] = 0.0;
	      emax=0.0;
	      if (p<*n)
		for (j=1;j<=INF(*n-p,*n-1);j++)
		  {
		    if (j<*n-p)
		      emax += SQR( f[index[*n-j]] - f[index[*n-j-1]] );
		    else
		      emax += SQR( f[index[p]] );
		    current_node_plus[j]=subset2binary(index+*n-j,j);
		    path_mu_plus[j]=mu[current_node_plus[j]];
		  }
	      /* Cas des valeurs negatives */
	      path_mu_minus[0] = 0.0;
	      if (p>0)
		for (j=1;j<=INF(p,*n-1);j++)
		  {
		    if (j<p)
		      emax += SQR( f[index[j-1]] - f[index[j]] );
		    else
		      emax += SQR( f[index[p-1]] );
		    current_node_minus[j]=subset2binary(index,j);
		    path_mu_minus[j]=mu[current_node_minus[j]];
		  }
        
	      /* modification des coefficients utilises dans le calcul precedent */
	      if( error > 0. )
		{
		  /* Pour les coefficients de la mesure correspondant aux valeurs
		     positives, on commence par modifier path_mu_plus[1], 
		     path_mu_plus[2], ... */
		  if (p<*n)
		    for( j=1; j<=INF(*n-p,*n-1); j++)
		      {
			/* Calcul de la direction du gradient */
			if (j<*n-p)
			  grad=f[index[*n-j]] - f[index[*n-j-1]];
			else
			  grad=f[index[p]];
			/* Effectue la modification */
			if( grad != 0.0 )
			  {
			    path_mu_plus[j] -= *alpha * error / emax * grad;
			    touched_mu[current_node_plus[j]] = 1;
              
			    /* STEP 1.3 : verification des contraintes de 
			       monotonie avec les voisins inferieurs 
			       (il y en a j)*/
			    search_lower_neighbors(*n, current_node_plus[j], 
						   lower_neighbors);
			    for( l=0; l<j; l++ )
			      {
				if( mu[ lower_neighbors[l] ] > path_mu_plus[ j ] )
				  path_mu_plus[ j ] = mu[ lower_neighbors[l] ];
			      }
			    /* update */
			    mu[current_node_plus[j]] = path_mu_plus[j];
			  }
		      }
		  /* Pour les coefficients de la mesure correspondant aux valeurs 
		     negatives, on commence par modifier path_mu_minus[p], 
		     path_mu_minus[p-1], ... */
		  if (p>0)
		    for( j=INF(p,*n-1); j>0; j--)
		      {
			/* Calcul de la direction du gradient */
			if (j<p)
			  grad=f[index[j-1]] - f[index[j]];
			else
			  grad=f[index[p-1]];
			/* Effectue la modification */
			if( grad != 0. )
			  {
			    path_mu_minus[j] -= *alpha * error / emax * grad;
			    touched_mu[current_node_minus[j]] = 1;
                
			    /* STEP 1.3 : verification des contraintes de 
			       monotonie avec les voisins 
			       superieurs (il y en a n-j)*/
			    search_upper_neighbors(*n, current_node_minus[j], 
						   upper_neighbors);
			    for( l=0; l<*n-j; l++ )
			      {
				if( mu[ upper_neighbors[l] ] < path_mu_minus[ j ] )
				  path_mu_minus[ j ] = mu[ upper_neighbors[l] ];
			      }
			    /* update */
			    mu[current_node_minus[j]] = path_mu_minus[j];
			  }
		      }
		}
	      else if( error < 0. )
		{
		  /* Pour les coefficients de la mesure correspondant aux valeurs 
		     positives, on commence par modifier path_mu_plus[n-p], 
		     path_mu_plus[n-p-1], ... */
		  if (p<*n)
		    for( j=INF(*n-p,*n-1); j>0; j--)
		      {
			/* Calcul de la direction du gradient */
			if (j<*n-p)
			  grad=f[index[*n-j]] - f[index[*n-j-1]];
			else
			  grad=f[index[p]];
			/* Effectue la modification */
			if( grad != 0.0 )
			  {
			    path_mu_plus[j] -= *alpha * error / emax * grad;
			    touched_mu[current_node_plus[j]] = 1;
                
			    /* STEP 1.3 : verification des contraintes de 
			       monotonie avec les voisins 
			       superieurs (il y en a n-j)*/
			    search_upper_neighbors(*n, current_node_plus[j], 
						   upper_neighbors);
			    for( l=0; l<*n-j; l++ )
			      {
				if( mu[ upper_neighbors[l] ] < path_mu_plus[ j ] )
				  path_mu_plus[ j ] = mu[ upper_neighbors[l] ];
			      }
			    /* update */
			    mu[current_node_plus[j]] = path_mu_plus[j];
			  }
		      }
		  /* Pour les coefficients de la mesure correspondant aux valeurs 
		     negatives, on commence par modifier path_mu_minus[1], 
		     path_mu_minus[2], ... */
		  if (p>0)
		    for( j=1; j<=INF(p,*n-1); j++)
		      {
			/* Calcul de la direction du gradient */
			if (j<p)
			  grad=f[index[j-1]] - f[index[j]];
			else
			  grad=f[index[p-1]];
			/* Effectue la modification */
			if( grad != 0. )
			  {
			    path_mu_minus[j] -= *alpha * error / emax * grad;
			    touched_mu[current_node_minus[j]] = 1;

			    /* STEP 1.3 : verification des contraintes de 
			       monotonie avec les voisins 
			       inferieurs (il y en a j)*/
			    search_lower_neighbors(*n, current_node_minus[j], 
						   lower_neighbors);
			    for( l=0; l<j; l++ )
			      {
				if( mu[ lower_neighbors[l] ] > path_mu_minus[ j ] )
				  path_mu_minus[ j ] = mu[ lower_neighbors[l] ];
			      }
			    /* update */
			    mu[current_node_minus[j]] = path_mu_minus[j];
			  }
		      }
		}
	    } /* end of for(i=0 ... */

	  *alpha *= fading;

	  /* STEP 2 : equilibrate the fuzzy measure for untouched coefficients */

	  /* lower levels first */
	  for( jj=1; jj<*n; jj++ ) {
	    for ( i=1; i<(1<<*n)-1; i++ )
	      if( touched_mu[i] == 0 )
		if( cardinal( i ) == jj )
		  {
		    /* look for upper neighbors */
		    search_upper_neighbors(*n, i, upper_neighbors);
	    
		    /* look for lower neighbors */
		    search_lower_neighbors(*n, i, lower_neighbors);
		    /* calcul de diverses constantes */
		    mean_upper_mu = 0;
		    for( l=0; l<*n-jj; l++)
		      mean_upper_mu += mu[ upper_neighbors[l] ];
		    mean_upper_mu /= *n-jj;
		    min_upper_dist = mu[ upper_neighbors[0] ] - mu[i];
		    for( l=1; l<*n-jj; l++ )
		      if( mu[ upper_neighbors[l] ] - mu[i] < min_upper_dist )
			min_upper_dist = mu[ upper_neighbors[l] ] - mu[i];
	      
		    mean_lower_mu = 0;
		    for( l=0; l<jj; l++)
		      mean_lower_mu += mu[ lower_neighbors[l] ];
		    mean_lower_mu /= jj;
		    min_lower_dist = mu[i] - mu[ lower_neighbors[0] ];
		    for( l=1; l<jj; l++ )
		      if( mu[i] - mu[ lower_neighbors[l] ] < min_lower_dist )
			min_lower_dist = mu[i] - mu[ lower_neighbors[l] ];
	    
		    /* update */ 
		    /* if( mean_upper_mu + mean_lower_mu - 2*mu[i] > 0. ) */
		      /* correction vers augmentation */
		      /* mu[i] += *beta * ( mean_upper_mu + mean_lower_mu - 2*mu[i] )
			 / ( mean_upper_mu - mean_lower_mu ) 
			 * min_upper_dist / 2.; */
		    /* else */
		      /* correction vers diminution */
		      /* mu[i] += *beta * ( mean_upper_mu + mean_lower_mu - 2*mu[i] ) 
			/ ( mean_upper_mu - mean_lower_mu ) 
			* min_lower_dist / 2.; */

		    mu[i] = mu[i] + (min_upper_dist - min_lower_dist)/2.;
		  }      
	  }
	  /* *beta *= fading; */
	  iter++;
	}  /* end of while(iter... */
    }
  *error2 = sqrt( criterion/(double) *n_data );
  Rprintf("hlms: L2 ending value of criterion: %g at iteration: %d\n", *error2, --iter);
  *itmax = iter;

  path_mu_plus = NULL;
  path_mu_minus = NULL;
  touched_mu = NULL;
  index = NULL;
  current_node_plus = NULL;
  current_node_minus = NULL;
  touched_mu = NULL;
  upper_neighbors = NULL;
  lower_neighbors = NULL;
}
