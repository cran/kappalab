##############################################################################
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
##################################################################################

## Class game extends set.func

##############################################################################

## Constructor from numeric
game <- function(object) {
    
    sf <- set.func.internal(object)
    new("game", data = sf$data, subsets = sf$subsets, n = sf$n)
}

## Constructor from set.func
setMethod("as.game", signature(object = "set.func"),
          function(object,...) {
              
              new("game", data = object@data, subsets = object@subsets,
                  n = object@n)
          }
          )

## Constructor
setMethod("as.game", signature(object = "card.game"),
          function(object,...) {
              
              sf <- as.set.func.internal(object)
              new("game", data = sf$data, subsets = sf$subsets, n = sf$n)
          }
          )

## Constructor from Mobius.game
setMethod("zeta", signature(object = "Mobius.game"),
          function(object, ...) {
              
              sf <- zeta.internal(object)
              new("game", data = sf$data, subsets = sf$subsets, n = sf$n)
          }
          )

## Constructor from set.func
## The dual of a set function is a game
setMethod("dual", signature(object = "set.func"),
          function(object, ...) {
              
              sf <- dual.internal(object)
              new("game", data = sf$data, subsets = sf$subsets, n = sf$n)
          }
          )

##############################################################################

## Computes the Choquet integral of f w.r.t mu 
setMethod("Choquet.integral",signature(object = "game", f = "numeric"),
          function(object,f,...) {
              
              if (!(object@n == length(f)))
                  stop("wrong arguments")
          
              .C("Choquet_integral_game", 
                 as.integer(object@n), 
                 as.double(object@data), 
                 as.double(f),  
                 res = double(1),
                 PACKAGE="kappalab")$res[1]
          }
          )

## Computes the Sugeno integral of f w.r.t mu (f >= 0)
setMethod("Sugeno.integral",signature(object = "game", f = "numeric"),
          function(object,f,...) {
              
              if (!(is.positive(f) && object@n == length(f)))
                  stop("wrong arguments")

              if (is(object)[1] == "capacity" && sum(f <=
                    object@data[2^object@n]) != object@n)
                  stop("wrong arguments")
              
              .C("Sugeno_integral_game", 
                 as.integer(object@n), 
                 as.double(object@data), 
                 as.double(f),  
                 res = double(1),
                 PACKAGE="kappalab")$res[1]
          }
          )

## Computes the Sipos integral of f w.r.t mu 
setMethod("Sipos.integral",signature(object = "game", f = "numeric"),
          function(object,f,...) {
              
              if (!(is.double(f) && object@n == length(f)))
                  stop("wrong arguments")
              
              .C("Sipos_integral_game", 
                 as.integer(object@n), 
                 as.double(object@data), 
                 as.double(f), 
                 res = double(1),
                 PACKAGE="kappalab")$res[1]
          }
          )

##############################################################################
