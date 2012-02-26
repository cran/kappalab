pkgname <- "kappalab"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('kappalab')

assign(".oldSearch", search(), pos = 'CheckExEnv')
cleanEx()
nameEx("Choquet.integral-methods")
### * Choquet.integral-methods

flush(stderr()); flush(stdout())

### Name: Choquet.integral-methods
### Title: Choquet integral
### Aliases: Choquet.integral Choquet.integral-methods
###   Choquet.integral,Mobius.game,numeric-method
###   Choquet.integral,card.game,numeric-method
###   Choquet.integral,game,numeric-method
### Keywords: methods

### ** Examples

## a normalized capacity
mu <- capacity(c(0:13/13,1,1))

## and its Mobius transform
a <- Mobius(mu)

## a discrete positive function f
f <- c(0.1,0.9,0.3,0.8)

## the Choquet integral of f w.r.t mu
Choquet.integral(mu,f)
Choquet.integral(a,f)

## a similar example with a cardinal capacity
mu <- uniform.capacity(4)
Choquet.integral(mu,f)



cleanEx()
nameEx("Mobius-methods")
### * Mobius-methods

flush(stderr()); flush(stdout())

### Name: Mobius-methods
### Title: The Möbius transform
### Aliases: Mobius Mobius-methods Mobius,capacity-method
###   Mobius,card.set.func-method Mobius,game-method Mobius,set.func-method
### Keywords: methods

### ** Examples

## a capacity
mu <- capacity(0:15)
mu

## its Mobius transform
a <- Mobius(mu)
a

## its zeta transform
zeta(a)

## a similar example with a game object
mu <- game(c(0,-2:12))
mu
Mobius(mu)
zeta(Mobius(mu))

## a similar example with a set.func object
mu <- set.func(-7:8)
mu
Mobius(mu)
zeta(Mobius(mu))

## a similar example with a card.set.func object
mu <- card.set.func(-3:4)
mu
Mobius(mu)
zeta(Mobius(mu))



cleanEx()
nameEx("Mobius.capacity-class")
### * Mobius.capacity-class

flush(stderr()); flush(stdout())

### Name: Mobius.capacity-class
### Title: Class "Mobius.capacity"
### Aliases: Mobius.capacity-class
### Keywords: classes

### ** Examples

## a capacity
mu <- capacity(c(0,0,0:13))
## and its Mobius representation
a <- Mobius(mu)
a

# the attributes of object a
a@n
a@k
a@data
a@subsets

## a test
is.normalized(a)
## normalize it
normalize(a)

## a transformation
zeta(a)
## Let us check ...
Mobius(zeta(a))

## some summary indices
orness(a)
veto(a)
favor(a)
variance(a)
entropy(a)
## the same
summary(a)



cleanEx()
nameEx("Mobius.card.gen")
### * Mobius.card.gen

flush(stderr()); flush(stdout())

### Name: Mobius.card.set.func
### Title: Creates an object representing the Möbius transform of a
###   cardinal set function.
### Aliases: Mobius.card.set.func
### Keywords: math

### ** Examples

Mobius.card.set.func(4:-2)



cleanEx()
nameEx("Mobius.card.set.func-class")
### * Mobius.card.set.func-class

flush(stderr()); flush(stdout())

### Name: Mobius.card.set.func-class
### Title: Class "Mobius.card.set.func"
### Aliases: Mobius.card.set.func-class
### Keywords: classes

### ** Examples

## the Mobius representation of a cardinal set function
a <- Mobius.card.set.func(-3:2)

# the attributes of the object
a@n
a@data

## some transformations
as.set.func(a)
zeta(a)
## let us check ...
Mobius(zeta(a))



cleanEx()
nameEx("Mobius.game-class")
### * Mobius.game-class

flush(stderr()); flush(stdout())

### Name: Mobius.game-class
### Title: Class "Mobius.game"
### Aliases: Mobius.game-class
### Keywords: classes

### ** Examples

## a game (which is a capacity)
mu <- game(c(0,rep(1,15)))
## and its Mobius representation
a <- Mobius(mu)

# the attributes of object a
a@n
a@k
a@data
a@subsets

## a transformation
zeta(a)
## let us check ...
Mobius(zeta(a))

## integral calculations 
f <- c(0.2,0.3,0.1,0.7)
Choquet.integral(a,f)
Sugeno.integral(a,f)
f <- c(0.2,-0.3,0.1,-0.7)
Sipos.integral(a,f)



cleanEx()
nameEx("Mobius.gen")
### * Mobius.gen

flush(stderr()); flush(stdout())

### Name: Mobius.set.func
### Title: Create objects representing the Möbius transform of a set
###   function.
### Aliases: Mobius.set.func Mobius.game Mobius.capacity additive.capacity
### Keywords: math

### ** Examples

Mobius.set.func(c(0,1,1,2,1,2,1,2,2,3,2),4,2)
Mobius.game(c(0,1,1,2,1,2,1,2,2,3,2),4,2)
Mobius.capacity(c(0,1,1,2,1,2,1,2,2,3,2),4,2)
additive.capacity(c(1,2,3,4))



cleanEx()
nameEx("Mobius.set.func-class")
### * Mobius.set.func-class

flush(stderr()); flush(stdout())

### Name: Mobius.set.func-class
### Title: Class "Mobius.set.func"
### Aliases: Mobius.set.func-class
### Keywords: classes

### ** Examples

## the Mobius transform of a set function directly
a <- Mobius.set.func(1:16,4,4)

## the attributes of the object
a@n
a@k
a@data
a@subsets

## a set function
mu <- set.func(7:-8)
## and its Mobius transform
a <- Mobius(mu)

## some conversions that cannot work
## as.game(a)
## as.capacity(a)
## as.card.set.func(a)

## some tests
is.cardinal(a)
is.kadditive(a,2)
is.monotone(a)

## some transformations
zeta(a)
k.truncate.Mobius(a,2)

## summary 
Shapley.value(a)
interaction.indices(a)
# the same
summary(a)

## save the Mobius transform to a file
d <- to.data.frame(a)
write.table(d,"my.Mobius.set.func.csv",sep="\t")

# finally, some conversions that should work
mu <- set.func(c(0,1,1,1,2,2,2,3))
a <- Mobius(mu)
as.Mobius.game(a)
as.Mobius.capacity(a)
as.Mobius.card.set.func(a)



cleanEx()
nameEx("Shapley.value-methods")
### * Shapley.value-methods

flush(stderr()); flush(stdout())

### Name: Shapley.value-methods
### Title: The Shapley value
### Aliases: Shapley.value Shapley.value-methods
###   Shapley.value,Mobius.set.func-method
###   Shapley.value,card.set.func-method Shapley.value,set.func-method
### Keywords: methods

### ** Examples

## a set function
mu <- set.func(c(0:13/13,1,1))

## the Shapley value
Shapley.value(mu)

## the efficiency property should be satisfied
sum(Shapley.value(mu))

## a similar example using a Mobius.set.func object
a <- Mobius(mu)
Shapley.value(a)

## a similar example using a card.set.func object
mu <- upper.capacity(6)
Shapley.value(mu)
## the efficiency property should be satisfied
Shapley.value(mu)*6



cleanEx()
nameEx("Sipos.integral-methods")
### * Sipos.integral-methods

flush(stderr()); flush(stdout())

### Name: Sipos.integral-methods
### Title: Sipos integral
### Aliases: Sipos.integral Sipos.integral-methods
###   Sipos.integral,Mobius.game,numeric-method
###   Sipos.integral,card.game,numeric-method
###   Sipos.integral,game,numeric-method
### Keywords: methods

### ** Examples

## a normalized capacity
mu <- capacity(c(0:13/13,1,1))

## and its Mobius transform
a <- Mobius(mu)

## a discrete function f
f <- c(0.1,-0.9,-0.3,0.8)

## the Sugeno integral of f w.r.t mu
Sipos.integral(mu,f)
Sipos.integral(a,f)

## a similar example with a cardinal capacity
mu <- uniform.capacity(4)
Sipos.integral(mu,f)



cleanEx()
nameEx("Sugeno.integral-methods")
### * Sugeno.integral-methods

flush(stderr()); flush(stdout())

### Name: Sugeno.integral-methods
### Title: Sugeno integral
### Aliases: Sugeno.integral Sugeno.integral-methods
###   Sugeno.integral,Mobius.game,numeric-method
###   Sugeno.integral,card.game,numeric-method
###   Sugeno.integral,game,numeric-method
### Keywords: methods

### ** Examples

## a normalized capacity
mu <- capacity(c(0:13/13,1,1))

## and its Mobius transform
a <- Mobius(mu)

## a discrete function f
f <- c(0.1,0.9,0.3,0.8)

## the Sugeno integral of f w.r.t mu
Sugeno.integral(mu,f)
Sugeno.integral(a,f) 

## a similar example with a cardinal capacity
mu <- uniform.capacity(4)
Sugeno.integral(mu,f)



cleanEx()
nameEx("capacity-class")
### * capacity-class

flush(stderr()); flush(stdout())

### Name: capacity-class
### Title: Class "capacity"
### Aliases: capacity-class
### Keywords: classes

### ** Examples

## a capacity
mu <- capacity(c(0:13,13,13)/13)

## the attributes of the object
mu@n
mu@data
mu@subsets

## a test
is.normalized(mu)
normalize(mu)

## a conversion that should not work
## as.card.capacity(mu)

## some transformations
conjugate(mu)
Mobius(mu)
## let us check ...
zeta(Mobius(mu))

## some summary indices
orness(mu)
veto(mu)
favor(mu)
variance(mu)
entropy(mu)
## the same
summary(mu)



cleanEx()
nameEx("card.capacity-class")
### * card.capacity-class

flush(stderr()); flush(stdout())

### Name: card.capacity-class
### Title: Class "card.capacity"
### Aliases: card.capacity-class
### Keywords: classes

### ** Examples

## a capacity
mu <- card.capacity(0:6/6)
## the same
mu <- uniform.capacity(6)

# the attributes of the object
mu@n
mu@data

## a test
is.normalized(mu)
normalize(mu)

## a transformation
conjugate(mu)

## some summary indices
orness(mu)
veto(mu)
favor(mu)
variance(mu)
entropy(mu)
## the same
summary(mu)



cleanEx()
nameEx("card.game-class")
### * card.game-class

flush(stderr()); flush(stdout())

### Name: card.game-class
### Title: Class "card.game"
### Aliases: card.game-class
### Keywords: classes

### ** Examples

## a cardinal game (which is a capacity)
mu <- card.game(c(0,rep(1,4)))

# the attributes of the object
mu@n
mu@data

## a conversion
as.game(mu)

## integral calculations 
f <- c(0.2,0.3,0.1,0.7)
Choquet.integral(mu,f)
Sugeno.integral(mu,f)
f <- c(0.2,-0.3,0.1,-0.7)
Sipos.integral(mu,f)



cleanEx()
nameEx("card.set.func-class")
### * card.set.func-class

flush(stderr()); flush(stdout())

### Name: card.set.func-class
### Title: Class "card.set.func"
### Aliases: card.set.func-class
### Keywords: classes

### ** Examples

## a cardinal set function
mu <- card.set.func(-3:2)

# the attributes of the object
mu@n
mu@data

## some conversions that cannot work
## Not run: as.card.game(mu)
## Not run: as.card.capacityfunc(mu)

## the following should work
as.set.func(mu)

## some tests
is.cardinal(mu)
is.kadditive(mu,2)
is.monotone(mu)

## some transformations
conjugate(mu)
Mobius(mu)
## let us check ...
zeta(Mobius(mu))

## summary 
Shapley.value(mu)
interaction.indices(mu)
# the same
summary(mu)

## save the set function to a file
d <- to.data.frame(mu)
write.table(d,"my.card.set.func.csv",sep="\t")

# finally, some other conversions that should work
mu <- card.set.func(0:5)
as.card.game(mu)
as.card.capacity(mu)



cleanEx()
nameEx("card.set.func.gen")
### * card.set.func.gen

flush(stderr()); flush(stdout())

### Name: card.set.func
### Title: Create objects representing cardinal set functions.
### Aliases: card.set.func card.game card.capacity lower.capacity
###   upper.capacity uniform.capacity
### Keywords: math

### ** Examples

card.set.func(4:-2)
card.game(c(0,-2:2))
card.capacity(0:5)
lower.capacity(3)
upper.capacity(4)
uniform.capacity(5)



cleanEx()
nameEx("conjugate-methods")
### * conjugate-methods

flush(stderr()); flush(stdout())

### Name: conjugate-methods
### Title: The conjugate (or dual) transform
### Aliases: conjugate conjugate-methods conjugate,capacity-method
###   conjugate,card.capacity-method conjugate,card.set.func-method
###   conjugate,set.func-method
### Keywords: methods

### ** Examples

## a game
mu <- game(c(0,-7:7))
mu

## its conjugate 
conjugate(mu)
## and mu again 
conjugate(conjugate(mu))

## a similar example with the upper capacity
mu <- capacity(c(0,rep(1,15)))
mu
conjugate(mu)
conjugate(conjugate(mu))

## a similar example with an object of class card.capacity
mu <- upper.capacity(6)
mu
conjugate(mu)
conjugate(conjugate(mu))

## the conjugate of a set function is a game
mu <- set.func(-7:8)
mu
conjugate(mu)
mu <- card.set.func(-2:5)
conjugate(mu)



cleanEx()
nameEx("entropy-methods")
### * entropy-methods

flush(stderr()); flush(stdout())

### Name: entropy-methods
### Title: Normalized entropy of a capacity
### Aliases: entropy entropy-methods entropy,Mobius.capacity-method
###   entropy,capacity-method entropy,card.capacity-method
### Keywords: methods

### ** Examples

## a capacity
mu <- capacity(c(0,0,0:13))

## its Mobius transform
a <- Mobius(mu)

## their normalized entropy
entropy(mu)
entropy(a)

## similar examples with card.capacity objects
mu <- lower.capacity(4)
entropy(mu)
mu <- uniform.capacity(4)
entropy(mu)



cleanEx()
nameEx("entropy.capa.ident")
### * entropy.capa.ident

flush(stderr()); flush(stdout())

### Name: entropy.capa.ident
### Title: Unsupervised identification of a capacity from profiles
### Aliases: entropy.capa.ident
### Keywords: math

### ** Examples

## a set of randomly generated data
## for instance, marks on a [0,20] scale
p <- data.frame(matrix(runif(500,0,20),100,5))
names(p) <- c("Stat","Prob","Alg","Cal","Eng")

## discretization
p[p <= 5] <- 1  
p[p > 5 & p <= 10] <- 2 
p[p > 10 & p <= 15] <- 3 
p[p > 15] <- 4

d <- data.frame(factor(p[[1]]),
                factor(p[[2]]),
                factor(p[[3]]),
                factor(p[[4]]),
                factor(p[[5]]))

## associated unsupervised capacity
mu <- entropy.capa.ident(d)
mu



cleanEx()
nameEx("expect.Choquet.unif-methods")
### * expect.Choquet.unif-methods

flush(stderr()); flush(stdout())

### Name: expect.Choquet.unif-methods
### Title: Expectation and standard deviation of the Choquet integral in
###   the uniform and normal cases
### Aliases: expect.Choquet.unif sd.Choquet.unif
###   expect.Choquet.unif-methods sd.Choquet.unif-methods
###   expect.Choquet.unif,game-method sd.Choquet.unif,game-method
###   expect.Choquet.norm sd.Choquet.norm expect.Choquet.norm-methods
###   sd.Choquet.norm-methods expect.Choquet.norm,game-method
###   sd.Choquet.norm,game-method expect.Choquet.norm,Mobius.game-method
### Keywords: methods

### ** Examples


## a capacity
mu <- capacity(c(0,0.1,0.6,rep(0.9,4),1))

## the expectation and the standard deviation
## of the Choquet integral in the uniform case
expect.Choquet.unif(mu)
sd.Choquet.unif(mu)

## the same but empirically
m <- 10000
ch <- numeric(m)
for (i in 1:m) {
     f <- runif(3) 
     ch[i] <- Choquet.integral(mu,f)
}
mean(ch)
sd(ch)

## the expectation and the standard deviation
## of the Choquet integral in the normal case
expect.Choquet.norm(mu)
sd.Choquet.norm(mu)
expect.Choquet.norm(Mobius(mu))

## the same but empirically
for (i in 1:m) {
     f <- rnorm(3) 
     ch[i] <- Choquet.integral(mu,f)
}
mean(ch)
sd(ch)




cleanEx()
nameEx("favor-methods")
### * favor-methods

flush(stderr()); flush(stdout())

### Name: favor-methods
### Title: Favor indices
### Aliases: favor favor-methods favor,Mobius.capacity-method
###   favor,capacity-method favor,card.capacity-method
### Keywords: methods

### ** Examples

## a capacity
mu <- capacity(c(0:13,13,13))

## its Mobius transform
a <- Mobius(mu)

## their favor indices
favor(mu)
favor(a)

## the same with a card.capacity object
mu <- lower.capacity(4)
favor(mu)



cleanEx()
nameEx("game-class")
### * game-class

flush(stderr()); flush(stdout())

### Name: game-class
### Title: Class "game"
### Aliases: game-class
### Keywords: classes

### ** Examples

## a game (which is a capacity)
mu <- game(c(0,rep(1,15)))

## the attributes of the object
mu@n
mu@data
mu@subsets

## a conversion
as.card.game(mu)

## a transformation
Mobius(mu)
## let us check ...
zeta(Mobius(mu))

## integral calculations 
f <- c(0.2,0.3,0.1,0.7)
Choquet.integral(mu,f)
Sugeno.integral(mu,f)
f <- c(0.2,-0.3,0.1,-0.7)
Sipos.integral(mu,f)



cleanEx()
nameEx("heuristic.ls.capa.ident")
### * heuristic.ls.capa.ident

flush(stderr()); flush(stdout())

### Name: heuristic.ls.capa.ident
### Title: Heuristic least squares capacity identification
### Aliases: heuristic.ls.capa.ident
### Keywords: math

### ** Examples


## number of criteria
n <- 4

## the number of alternatives
n.a <- 1000

## a randomly generated 5-criteria matrix
C <- matrix(rnorm(n*n.a,10,2),n.a,n)

## the corresponding global scores
g <- numeric(n.a)

## generate a random capacity
x <- runif(2^n-1)
for (i in 2:(2^n-1))
    x[i] <- x[i] + x[i-1]
mu <- normalize(capacity(c(0,x)))
for (i in 1:n.a)
  g[i] <- Choquet.integral(mu,C[i,])

## the initial capacity
## here the uniform capacity
mu.in <- as.capacity(uniform.capacity(n))

## the solution 
hlsc <- heuristic.ls.capa.ident(n,mu.in,C,g)
mu.sol <- hlsc$solution

## the difference between mu and mu.sol
mu@data - mu.sol@data

hlsc



cleanEx()
nameEx("interaction.indices-methods")
### * interaction.indices-methods

flush(stderr()); flush(stdout())

### Name: interaction.indices-methods
### Title: The Shapley interaction indices
### Aliases: interaction.indices interaction.indices-methods
###   interaction.indices,Mobius.set.func-method
###   interaction.indices,card.set.func-method
###   interaction.indices,set.func-method
### Keywords: methods

### ** Examples

## a set function
mu <- set.func(c(-7:6,6,6))

## the associated interaction indices
interaction.indices(mu)


## a similar example using a Mobius.set.func object
a <- Mobius(mu)
interaction.indices(a)

## a similar example using a card.set.func object
mu <- upper.capacity(6)
interaction.indices(mu)



cleanEx()
nameEx("is.cardinal-methods")
### * is.cardinal-methods

flush(stderr()); flush(stdout())

### Name: is.cardinal-methods
### Title: Test method
### Aliases: is.cardinal is.cardinal-methods
###   is.cardinal,Mobius.set.func-method is.cardinal,card.set.func-method
###   is.cardinal,set.func-method
### Keywords: methods

### ** Examples

is.cardinal(set.func(-7:8))
is.cardinal(uniform.capacity(8))
is.cardinal(Mobius.game(0:10,4,2))



cleanEx()
nameEx("is.kadditive-methods")
### * is.kadditive-methods

flush(stderr()); flush(stdout())

### Name: is.kadditive-methods
### Title: Test method
### Aliases: is.kadditive is.kadditive-methods
###   is.kadditive,Mobius.set.func,numeric-method
###   is.kadditive,card.set.func,numeric-method
###   is.kadditive,set.func,numeric-method
### Keywords: methods

### ** Examples

## a set function
mu <- set.func(c(0,1,1,1,2,2,2,3))
mu
is.kadditive(mu,2)
is.kadditive(mu,1)

## the Mobius representation of a set function, 2-additive by construction 
a <- Mobius.set.func(c(0,1,2,1,3,1,2,1,2,3,1),4,2)
is.kadditive(a,2)
is.kadditive(a,1)



cleanEx()
nameEx("is.monotone-methods")
### * is.monotone-methods

flush(stderr()); flush(stdout())

### Name: is.monotone-methods
### Title: Test method
### Aliases: is.monotone is.monotone-methods
###   is.monotone,Mobius.set.func-method is.monotone,card.set.func-method
###   is.monotone,set.func-method
### Keywords: methods

### ** Examples

## a monotone set function
mu <- set.func(c(0,1,1,1,2,2,2,3))
mu
is.monotone(mu)

## the Mobius representation of a monotone set function
a <- Mobius.set.func(c(0,1,2,1,3,1,2,1,2,3,1),4,2)
is.monotone(a)

## non-monotone examples
mu <- set.func(c(0,-7:7))
is.monotone(mu,verbose=TRUE)
a <- Mobius(mu)
is.monotone(a,verbose=TRUE)



cleanEx()
nameEx("is.normalized-methods")
### * is.normalized-methods

flush(stderr()); flush(stdout())

### Name: is.normalized-methods
### Title: Test method
### Aliases: is.normalized is.normalized-methods
###   is.normalized,Mobius.capacity-method is.normalized,capacity-method
###   is.normalized,card.capacity-method
### Keywords: methods

### ** Examples

## a capacity
mu <- capacity(0:15)
is.normalized(mu)
normalize(mu)

## its Mobius transform
a <- Mobius(mu)
is.normalized(a)
normalize(a)

## a cardinal capacity
mu <- uniform.capacity(7)
is.normalized(mu)



cleanEx()
nameEx("k.truncate.Mobius-methods")
### * k.truncate.Mobius-methods

flush(stderr()); flush(stdout())

### Name: k.truncate.Mobius-methods
### Title: k-order truncation of the Möbius representation of a set
###   function.
### Aliases: k.truncate.Mobius k.truncate.Mobius-methods
###   k.truncate.Mobius,set.func,numeric-method
###   k.truncate.Mobius,Mobius.set.func,numeric-method
### Keywords: methods

### ** Examples

## a set function
mu <- set.func(c(0,1,1,1,2,2,2,3))
mu

## 2-truncate it
k.truncate.Mobius(mu,2)
## 2-truncate it
k.truncate.Mobius(Mobius(mu),2)



cleanEx()
nameEx("least.squares.capa.ident")
### * least.squares.capa.ident

flush(stderr()); flush(stdout())

### Name: least.squares.capa.ident
### Title: Least squares capacity identification
### Aliases: least.squares.capa.ident
### Keywords: math

### ** Examples


## the number of data
n.d <- 20

## a randomly generated 5-criteria matrix
C <- matrix(rnorm(5*n.d,10,2),n.d,5)


## the corresponding global scores
g <- numeric(n.d)
mu <- capacity(c(0:29,29,29)/29)
for (i in 1:n.d)
  g[i] <- Choquet.integral(mu,C[i,])

## Not run: 
##D ## the full solution 
##D lsc <- least.squares.capa.ident(5,5,C,g)
##D a <- lsc$solution
##D a
##D mu.sol <- zeta(a)
##D 
##D ## the difference between mu and mu.sol
##D mu@data - mu.sol@data
##D 
##D ## the residuals
##D lsc$residuals
##D 
##D ## the mean square error
##D mean(lsc$residuals^2)
##D 
##D ## a 3-additive solution 
##D lsc <- least.squares.capa.ident(5,3,C,g)
##D a <- lsc$solution
##D mu.sol <- zeta(a)
##D mu@data - mu.sol@data
##D lsc$residuals
## End(Not run)



## a similar example based on the Sipos integral

## a randomly generated 5-criteria matrix
C <- matrix(rnorm(5*n.d,0,2),n.d,5)

## the corresponding global scores
g <- numeric(n.d)
mu <- capacity(c(0:29,29,29)/29)
for (i in 1:n.d)
  g[i] <- Sipos.integral(mu,C[i,])

## Not run: 
##D ## the full solution 
##D lsc <- least.squares.capa.ident(5,5,C,g,Integral = "Sipos")
##D a <- lsc$solution
##D mu.sol <- zeta(a)
##D mu@data - mu.sol@data
##D lsc$residuals
##D 
##D ## a 3-additive solution 
##D lsc <- least.squares.capa.ident(5,3,C,g,Integral = "Sipos")
##D a <- lsc$solution
##D mu.sol <- zeta(a)
##D mu@data - mu.sol@data
##D lsc$residuals
## End(Not run)



## additional constraints

## a Shapley preorder constraint matrix
## Sh(1) - Sh(2) >= -delta.S
## Sh(2) - Sh(1) >= -delta.S
## Sh(3) - Sh(4) >= -delta.S
## Sh(4) - Sh(3) >= -delta.S
## i.e. criteria 1,2 and criteria 3,4
## should have the same global importances
delta.S <- 0.01    
Asp <- rbind(c(1,2,-delta.S),
             c(2,1,-delta.S),
             c(3,4,-delta.S),
             c(4,3,-delta.S)
            )

## a Shapley interval constraint matrix
## 0.3 <= Sh(1) <= 0.9 
Asi <- rbind(c(1,0.3,0.9))


## an interaction preorder constraint matrix
## such that I(12) = I(45)
delta.I <- 0.01
Aip <- rbind(c(1,2,4,5,-delta.I),
             c(4,5,1,2,-delta.I))

## an interaction interval constraint matrix
## i.e. -0.20 <= I(12) <= -0.15 
delta.I <- 0.01
Aii <- rbind(c(1,2,-0.2,-0.15))

## an inter-additive partition constraint
## criteria 1,2,3 and criteria 4,5 are independent 
Aiap <- c(1,1,1,2,2)





## a more constrained solution

lsc <- least.squares.capa.ident(5,5,C,g,Integral = "Sipos",
                                 A.Shapley.preorder = Asp,
                                 A.Shapley.interval = Asi,
                                 A.interaction.preorder = Aip,
                                 A.interaction.interval = Aii,
                                 A.inter.additive.partition = Aiap,
                                 sigf = 5)

a <- lsc$solution
mu.sol <- zeta(a)
mu@data - mu.sol@data
lsc$residuals
summary(a)



cleanEx()
nameEx("lin.prog.capa.ident")
### * lin.prog.capa.ident

flush(stderr()); flush(stdout())

### Name: lin.prog.capa.ident
### Title: Capacity identification based on linear programming
### Aliases: lin.prog.capa.ident
### Keywords: math

### ** Examples

## some alternatives
a <- c(18,11,18,11,11)
b <- c(18,18,11,11,11)
c <- c(11,11,18,18,11)
d <- c(18,11,11,11,18)
e <- c(11,11,18,11,18)
    
## preference threshold relative
## to the preorder of the alternatives
delta.C <- 1

## corresponding Choquet preorder constraint matrix 
Acp <- rbind(c(d,a,delta.C),
             c(a,e,delta.C),
             c(e,b,delta.C),
             c(b,c,delta.C)
            )

## a Shapley preorder constraint matrix
## Sh(1) - Sh(2) >= -delta.S
## Sh(2) - Sh(1) >= -delta.S
## Sh(3) - Sh(4) >= -delta.S
## Sh(4) - Sh(3) >= -delta.S
## i.e. criteria 1,2 and criteria 3,4
## should have the same global importances
delta.S <- 0.01    
Asp <- rbind(c(1,2,-delta.S),
             c(2,1,-delta.S),
             c(3,4,-delta.S),
             c(4,3,-delta.S)
            )

## a Shapley interval constraint matrix
## 0.3 <= Sh(1) <= 0.9 
Asi <- rbind(c(1,0.3,0.9))


## an interaction preorder constraint matrix
## such that I(12) = I(34)
delta.I <- 0.01
Aip <- rbind(c(1,2,3,4,-delta.I),
             c(3,4,1,2,-delta.I))

## an interaction interval constraint matrix
## i.e. -0.20 <= I(12) <= -0.15 
Aii <- rbind(c(1,2,-0.2,-0.15))


## Not run: 
##D ## a LP 2-additive solution
##D lin.prog <- lin.prog.capa.ident(5,2,A.Choquet.preorder = Acp)              
##D m <- lin.prog$solution
##D m
##D 
##D ## the resulting global evaluations
##D rbind(c(a,mean(a),Choquet.integral(m,a)),
##D       c(b,mean(b),Choquet.integral(m,b)),
##D       c(c,mean(c),Choquet.integral(m,c)),
##D       c(d,mean(d),Choquet.integral(m,d)),
##D       c(e,mean(e),Choquet.integral(m,e)))
##D 
##D ## the Shapley value
##D Shapley.value(m)
##D 
##D ## a LP 3-additive more constrained solution
##D lin.prog2 <- lin.prog.capa.ident(5,3,
##D                                    A.Choquet.preorder = Acp,
##D                                    A.Shapley.preorder = Asp)
##D m <- lin.prog2$solution
##D m
##D rbind(c(a,mean(a),Choquet.integral(m,a)),
##D       c(b,mean(b),Choquet.integral(m,b)),
##D       c(c,mean(c),Choquet.integral(m,c)),
##D       c(d,mean(d),Choquet.integral(m,d)),
##D       c(e,mean(e),Choquet.integral(m,e)))
##D Shapley.value(m)
##D 
##D ## a LP 5-additive more constrained solution
##D lin.prog3 <- lin.prog.capa.ident(5,5,
##D                                    A.Choquet.preorder = Acp,
##D                                    A.Shapley.preorder = Asp,
##D                                    A.Shapley.interval = Asi,
##D                                    A.interaction.preorder = Aip,
##D                                    A.interaction.interval = Aii)
##D 
##D m <- lin.prog3$solution
##D m
##D rbind(c(a,mean(a),Choquet.integral(m,a)),
##D       c(b,mean(b),Choquet.integral(m,b)),
##D       c(c,mean(c),Choquet.integral(m,c)),
##D       c(d,mean(d),Choquet.integral(m,d)),
##D       c(e,mean(e),Choquet.integral(m,e)))
##D summary(m)
## End(Not run)



cleanEx()
nameEx("ls.ranking.capa.ident")
### * ls.ranking.capa.ident

flush(stderr()); flush(stdout())

### Name: ls.ranking.capa.ident
### Title: Least squares capacity identification in the framework of a
###   ranking procedure
### Aliases: ls.ranking.capa.ident
### Keywords: math

### ** Examples


arthur <- c(1, 1, 0.75, 0.25)
lancelot <- c(0.75, 0.75, 0.75, 0.75)
yvain <- c(1, 0.625, 0.5, 1)
perceval <- c(0.25, 0.5, 0.75, 0.75)
erec <- c(0.375, 1, 0.5 , 0.75)

P <- rbind(arthur, lancelot, yvain, perceval, erec)

# lancelot > erec; yvain > erec, erec > perceval, erec > arthur
rk.proto <- rbind(c("lancelot","erec"), c("yvain","erec"), c("erec","perceval"), c("erec","arthur"))

n<-4
k<-2
d<-0.1

## search for a capacity which satisfies the constraints
lrc <- ls.ranking.capa.ident(n ,k, P, rk.proto, d)

lrc




cleanEx()
nameEx("ls.sorting.capa.ident")
### * ls.sorting.capa.ident

flush(stderr()); flush(stdout())

### Name: ls.sorting.capa.ident
### Title: Least squares capacity identification in the framework of a
###   sorting procedure
### Aliases: ls.sorting.capa.ident
### Keywords: math

### ** Examples


## generate a random problem with 10 prototypes and 4 criteria
n.proto <- 10 ## prototypes
n <- 4  ## criteria
k <- 4  
d <- 0.1
	
## generating random data for the prototypes
C <- matrix(runif(n.proto*n,0,1),n.proto,n)
cl <- numeric(n.proto)

## the corresponding global evaluations
glob.eval <- numeric(n.proto)
a <- capacity(c(0:(2^n-3),(2^n-3),(2^n-3))/(2^n-3))
for (i in 1:n.proto)
  glob.eval[i] <- Choquet.integral(a,C[i,])

## and the classes for the prototypes
cl[glob.eval <= 0.33] <- 1
cl[glob.eval > 0.33 & glob.eval <= 0.66] <-2
cl[glob.eval > 0.66] <- 3

cl

## Not run: 
##D # starting the calculations
##D # search for a capacity which satisfies the constraints
##D lsc <- ls.sorting.capa.ident(n ,k, C, cl, d)
##D 
##D ## output of the quadratic program (ipop, package kernlab)
##D lsc$how
##D 
##D ## the capacity satisfying the constraints
##D lsc$solution
##D summary(lsc$solution)
##D ## the global evaluations satisfying the constraints
##D lsc$glob.eval
## End(Not run)

## let us now add some constraints        

## a Shapley preorder constraint matrix
## Sh(1) > Sh(2)
## Sh(3) > Sh(4)
delta.S <-0.01
Asp <- rbind(c(1,2,delta.S), c(3,4,delta.S))

## a Shapley interval constraint matrix
## 0.1 <= Sh(1) <= 0.2 
Asi <- rbind(c(1,0.1,0.2))
        
## an interaction preorder constraint matrix
## such that I(12) > I(34)
delta.I <- 0.01
Aip <- rbind(c(1,2,3,4,delta.I))
        
## an interaction interval constraint matrix
## i.e. 0.2 <= I(12) <= 0.4 
## and 0 < I(34) <= 1
Aii <- rbind(c(1,2,0.2,0.4), c(3,4,delta.I,1))
        
## an inter-additive partition constraint
## criteria 1,2 and criteria 3,4 are independent 
Aiap <- c(1,1,2,2)


## starting the calculations
## search for a capacity which satisfies the constraints
lsc <- ls.sorting.capa.ident(n ,k, C, cl, d,
                                  A.Shapley.preorder = Asp,
                                  A.Shapley.interval = Asi,
                                  A.interaction.preorder = Aip,
                                  A.interaction.interval = Aii,
                                  A.inter.additive.partition = Aiap)

## output of ipop
lsc$how
## the capacity satisfying the constraints
lsc$solution
summary(lsc$solution)
## the global evaluations satisfying the constraints
lsc$glob.eval




cleanEx()
nameEx("ls.sorting.treatment")
### * ls.sorting.treatment

flush(stderr()); flush(stdout())

### Name: ls.sorting.treatment
### Title: Least squares capacity identification in the framework of a
###   sorting procedure: evaluation of the determined capacity
### Aliases: ls.sorting.treatment
### Keywords: math

### ** Examples


## generate a random problem with 10 prototypes and 4 criteria
n.proto <- 10 ## prototypes
n <- 4  ## criteria
P <- matrix(runif(n.proto*n,0,1),n.proto,n)

## the corresponding global scores, based on a randomly generated
## capacity a
glob.eval <- numeric(n.proto)
a <- capacity(c(0:(2^n-3),(2^n-3),(2^n-3))/(2^n-3))
for (i in 1:n.proto)
  glob.eval[i] <- Choquet.integral(a,P[i,])

## based on these global scores, let us create a classification (3 classes)
cl.proto<-numeric(n.proto)
cl.proto[glob.eval <= 0.33] <- 1
cl.proto[glob.eval > 0.33 & glob.eval<=0.66] <-2
cl.proto[glob.eval > 0.66] <- 3

## search for a capacity which satisfies the constraints
lsc <- ls.sorting.capa.ident(n ,4, P, cl.proto, 0.1)

## output of the QP
lsc$how

## analyse the quality of the model (classify the prototypes by the
## model and compare both assignments)
lst <- ls.sorting.treatment(P,cl.proto,lsc$solution,P,cl.proto)

## assignments of the prototypes
lst$class.A
## assignment types
lst$correct.A
## evaluation
lst$eval.correct

## generate a second set of random alternatives (A)
## their "correct" class is determined as beforehand with the
## randomly generated capacity a
## the goal is to see if we can reproduce this classification
## by the capacity learnt from the prototypes

## a randomly generated criteria matrix of 10 alternatives
A <- matrix(runif(10*n,0,1),10,n)
cl.orig.A <-numeric(10)
	
## the corresponding global scores
glob.eval.A <- numeric(10)
for (i in 1:10)
  glob.eval.A[i] <- Choquet.integral(a,A[i,])

## based on these global scores, let us determine a classification
cl.orig.A[glob.eval.A <= 0.33] <- 1
cl.orig.A[glob.eval.A>0.33 & glob.eval.A<=0.66] <-2
cl.orig.A[glob.eval.A > 0.66] <- 3

## let us now classify the alternatives of A according to the model
## built on P

lst <- ls.sorting.treatment(P,cl.proto,lsc$solution,A,cl.orig.A)

## assignment of the alternatives of A
lst$class.A
## type of assignments
lst$correct.A
## evaluation
lst$eval.correct

## show the learnt capacity
## x11()
## barplot(Shapley.value(lsc$solution), main="Learnt capacity", sub="Shapley")
## summary of the learnt capacity
lsc$solution
summary(lsc$solution)




cleanEx()
nameEx("mini.dist.capa.ident")
### * mini.dist.capa.ident

flush(stderr()); flush(stdout())

### Name: mini.dist.capa.ident
### Title: Minimum distance capacity identification
### Aliases: mini.dist.capa.ident
### Keywords: math

### ** Examples

## some alternatives
a <- c(18,11,18,11,11)
b <- c(18,18,11,11,11)
c <- c(11,11,18,18,11)
d <- c(18,11,11,11,18)
e <- c(11,11,18,11,18)
    
## preference threshold relative
## to the preorder of the alternatives
delta.C <- 1

## corresponding Choquet preorder constraint matrix 
Acp <- rbind(c(d,a,delta.C),
             c(a,e,delta.C),
             c(e,b,delta.C),
             c(b,c,delta.C)
            )

## a Shapley preorder constraint matrix
## Sh(1) - Sh(2) >= -delta.S
## Sh(2) - Sh(1) >= -delta.S
## Sh(3) - Sh(4) >= -delta.S
## Sh(4) - Sh(3) >= -delta.S
## i.e. criteria 1,2 and criteria 3,4
## should have the same global importances
delta.S <- 0.01    
Asp <- rbind(c(1,2,-delta.S),
             c(2,1,-delta.S),
             c(3,4,-delta.S),
             c(4,3,-delta.S)
            )

## a Shapley interval constraint matrix
## 0.3 <= Sh(1) <= 0.9 
Asi <- rbind(c(1,0.3,0.9))


## an interaction preorder constraint matrix
## such that I(12) = I(34)
delta.I <- 0.01
Aip <- rbind(c(1,2,3,4,-delta.I),
             c(3,4,1,2,-delta.I))

## an interaction interval constraint matrix
## i.e. -0.20 <= I(12) <= -0.15 
Aii <- rbind(c(1,2,-0.2,-0.15))

## the capacity that we want to approach
x <- runif(31)
for (i in 2:31)
    x[i] <- x[i] + x[i-1]
mu <- normalize(capacity(c(0,x)))
## and its Mobius transform
a.mu <- Mobius(mu)

## some basic checks
## Not run: 
##D mini.dist.capa.ident(a.mu,5)
##D mini.dist.capa.ident(a.mu,5,"binary.alternatives")
##D mini.dist.capa.ident(a.mu,5,"global.scores")
##D mini.dist.capa.ident(a.mu,3)
##D mini.dist.capa.ident(a.mu,3,"binary.alternatives")
##D mini.dist.capa.ident(a.mu,3,"global.scores")
## End(Not run)

## a minimum distance 2-additive solution
min.dist <- mini.dist.capa.ident(a.mu,2,"binary.alternatives",A.Choquet.preorder = Acp)              
m <- min.dist$solution
m

## a minimum distance 3-additive more constrained solution
min.dist2 <- mini.dist.capa.ident(a.mu,3,"global.scores",
                                   A.Choquet.preorder = Acp,
                                   A.Shapley.preorder = Asp)
m <- min.dist2$solution
m
rbind(c(a,mean(a),Choquet.integral(m,a)),
      c(b,mean(b),Choquet.integral(m,b)),
      c(c,mean(c),Choquet.integral(m,c)),
      c(d,mean(d),Choquet.integral(m,d)),
      c(e,mean(e),Choquet.integral(m,e)))
Shapley.value(m)

## Not run: 
##D ## a minimum distance 5-additive more constrained solution
##D min.dist3 <- mini.dist.capa.ident(a.mu,5,
##D                                    A.Choquet.preorder = Acp,
##D                                    A.Shapley.preorder = Asp,
##D                                    A.Shapley.interval = Asi,
##D                                    A.interaction.preorder = Aip,
##D                                    A.interaction.interval = Aii)
##D 
##D m <- min.dist3$solution
##D m
##D rbind(c(a,mean(a),Choquet.integral(m,a)),
##D       c(b,mean(b),Choquet.integral(m,b)),
##D       c(c,mean(c),Choquet.integral(m,c)),
##D       c(d,mean(d),Choquet.integral(m,d)),
##D       c(e,mean(e),Choquet.integral(m,e)))
##D summary(m)
## End(Not run)




cleanEx()
nameEx("mini.var.capa.ident")
### * mini.var.capa.ident

flush(stderr()); flush(stdout())

### Name: mini.var.capa.ident
### Title: Minimum variance capacity identification
### Aliases: mini.var.capa.ident
### Keywords: math

### ** Examples

## some alternatives
a <- c(18,11,18,11,11)
b <- c(18,18,11,11,11)
c <- c(11,11,18,18,11)
d <- c(18,11,11,11,18)
e <- c(11,11,18,11,18)
    
## preference threshold relative
## to the preorder of the alternatives
delta.C <- 1

## corresponding Choquet preorder constraint matrix 
Acp <- rbind(c(d,a,delta.C),
             c(a,e,delta.C),
             c(e,b,delta.C),
             c(b,c,delta.C)
            )

## a Shapley preorder constraint matrix
## Sh(1) - Sh(2) >= -delta.S
## Sh(2) - Sh(1) >= -delta.S
## Sh(3) - Sh(4) >= -delta.S
## Sh(4) - Sh(3) >= -delta.S
## i.e. criteria 1,2 and criteria 3,4
## should have the same global importances
delta.S <- 0.01    
Asp <- rbind(c(1,2,-delta.S),
             c(2,1,-delta.S),
             c(3,4,-delta.S),
             c(4,3,-delta.S)
            )

## a Shapley interval constraint matrix
## 0.3 <= Sh(1) <= 0.9 
Asi <- rbind(c(1,0.3,0.9))


## an interaction preorder constraint matrix
## such that I(12) = I(34)
delta.I <- 0.01
Aip <- rbind(c(1,2,3,4,-delta.I),
             c(3,4,1,2,-delta.I))

## an interaction interval constraint matrix
## i.e. -0.20 <= I(12) <= -0.15 
Aii <- rbind(c(1,2,-0.2,-0.15))


## a minimum variance 2-additive solution
min.var <- mini.var.capa.ident(5,2,A.Choquet.preorder = Acp)              
m <- min.var$solution
m

## the resulting global evaluations
rbind(c(a,mean(a),Choquet.integral(m,a)),
      c(b,mean(b),Choquet.integral(m,b)),
      c(c,mean(c),Choquet.integral(m,c)),
      c(d,mean(d),Choquet.integral(m,d)),
      c(e,mean(e),Choquet.integral(m,e)))

## the Shapley value
Shapley.value(m)

## a minimum variance 3-additive more constrained solution
min.var2 <- mini.var.capa.ident(5,3,
                                   A.Choquet.preorder = Acp,
                                   A.Shapley.preorder = Asp)
m <- min.var2$solution
m
rbind(c(a,mean(a),Choquet.integral(m,a)),
      c(b,mean(b),Choquet.integral(m,b)),
      c(c,mean(c),Choquet.integral(m,c)),
      c(d,mean(d),Choquet.integral(m,d)),
      c(e,mean(e),Choquet.integral(m,e)))
Shapley.value(m)

## a minimum variance 5-additive more constrained solution
min.var3 <- mini.var.capa.ident(5,5,
                                   A.Choquet.preorder = Acp,
                                   A.Shapley.preorder = Asp,
                                   A.Shapley.interval = Asi,
                                   A.interaction.preorder = Aip,
                                   A.interaction.interval = Aii)

m <- min.var3$solution
m
rbind(c(a,mean(a),Choquet.integral(m,a)),
      c(b,mean(b),Choquet.integral(m,b)),
      c(c,mean(c),Choquet.integral(m,c)),
      c(d,mean(d),Choquet.integral(m,d)),
      c(e,mean(e),Choquet.integral(m,e)))
summary(m)




cleanEx()
nameEx("normalize-methods")
### * normalize-methods

flush(stderr()); flush(stdout())

### Name: normalize-methods
### Title: Normalizes a capacity.
### Aliases: normalize normalize-methods normalize,Mobius.capacity-method
###   normalize,capacity-method normalize,card.capacity-method
### Keywords: methods

### ** Examples

## a capacity
mu <- capacity(0:15)

## normalize it
is.normalized(mu)
normalize(mu)

## a similar example with a Mobius.capacity object
a <- Mobius(mu)
is.normalized(a)
normalize(a)
zeta(normalize(a))

## a similar example with a card.capacity object
mu <- card.capacity(0:6)
is.normalized(mu)
normalize(mu)
Mobius(normalize(mu))




cleanEx()
nameEx("orness-methods")
### * orness-methods

flush(stderr()); flush(stdout())

### Name: orness-methods
### Title: Orness degree
### Aliases: orness orness-methods orness,Mobius.capacity-method
###   orness,capacity-method orness,card.capacity-method
### Keywords: methods

### ** Examples

## the upper capacity
mu <- capacity(c(0,rep(1,15)))

## the Choquet integral w.r.t mu behaves like the maximum operator
f <- c(0.1,0.1,0,0.9)
Choquet.integral(mu,f)

## its orness is 1
orness(mu)

## the same example with a Mobius.capacity object
a <- Mobius(mu)
Choquet.integral(a,f)
orness(a)

## the same example with a card.capacity object
mu <- upper.capacity(4)
Choquet.integral(mu,f)
orness(mu)



cleanEx()
nameEx("pdf.Choquet.unif-methods")
### * pdf.Choquet.unif-methods

flush(stderr()); flush(stdout())

### Name: pdf.Choquet.unif-methods
### Title: Distribution of the Choquet integral for evaluations uniformly
###   distributed on the unit hypercube
### Aliases: pdf.Choquet.unif cdf.Choquet.unif pdf.Choquet.unif-methods
###   cdf.Choquet.unif-methods pdf.Choquet.unif,game,numeric-method
###   cdf.Choquet.unif,game,numeric-method
### Keywords: methods

### ** Examples


## a capacity
mu <- capacity(c(0,0.1,0.6,rep(0.9,4),1))
## the cdf of the Choquet integral at 0.7
cdf.Choquet.unif(mu,0.7)

## the same but empirically
m <- 10000
ch <- numeric(m)
for (i in 1:m) {
     f <- runif(3) 
     ch[i] <- Choquet.integral(mu,f)
}
sum(ifelse(ch<=0.7,1,0))/m



cleanEx()
nameEx("rnd-methods")
### * rnd-methods

flush(stderr()); flush(stdout())

### Name: rnd-methods
### Title: Rounding of set function coefficients
### Aliases: rnd rnd-methods rnd,superclass.set.func-method
### Keywords: methods

### ** Examples

## a capacity
mu <- capacity(0:15/15)
mu
rnd(mu,2)

## a similar example with a Mobius.capacity object
a <- Mobius(mu)
a
rnd(a,1)

## a similar example with a card.capacity object
mu <- uniform.capacity(6)
mu
rnd(mu)



cleanEx()
nameEx("set.func-class")
### * set.func-class

flush(stderr()); flush(stdout())

### Name: set.func-class
### Title: Class "set.func"
### Aliases: set.func-class
### Keywords: classes

### ** Examples

## a set function
mu <- set.func(c(1:8,8:1)/8)

## the attributes of the object
mu@n
mu@data
mu@subsets

## some conversions that cannot work
## Not run: as.game(mu)
## Not run: as.capacity(mu)
## Not run: as.card.set.func(mu)

## some tests
is.cardinal(mu)
is.kadditive(mu,2)
is.monotone(mu)

## some transformations
conjugate(mu)
Mobius(mu)
k.truncate.Mobius(mu,2)

## summary 
Shapley.value(mu)
interaction.indices(mu)
# the same
summary(mu)

## save the set function to a file
d <- to.data.frame(mu)
write.table(d,"my.set.func.csv",sep="\t")

# finally, some conversions that should work
mu <- set.func(c(0,1,1,1,2,2,2,3))
as.game(mu)
as.capacity(mu)
as.card.set.func(mu)



cleanEx()
nameEx("set.func.gen")
### * set.func.gen

flush(stderr()); flush(stdout())

### Name: set.func
### Title: Create objects of class "set.func", "game", or "capacity".
### Aliases: set.func game capacity
### Keywords: math

### ** Examples

set.func(-12:3)
game(c(0,rep(-1,15)))
capacity(0:15)



cleanEx()
nameEx("show-methods")
### * show-methods

flush(stderr()); flush(stdout())

### Name: show-methods
### Title: Methods for Function show in Package 'kappalab'
### Aliases: show-methods show,superclass.set.func-method
###   show,set.func-method show,Mobius.set.func-method
###   show,summary.superclass.set.func-method
###   show,summary.superclass.capacity-method
### Keywords: methods

### ** Examples

## a set function
mu <- set.func(0:15/15)
show(mu)
## the same
mu

## a Mobius transform
a <- Mobius.set.func(0:10,4,2)
show(a)
a

## a cardinal capacity
mu <- lower.capacity(5)
show(mu)
mu



cleanEx()
nameEx("summary.superclass.capacity-class")
### * summary.superclass.capacity-class

flush(stderr()); flush(stdout())

### Name: summary.superclass.capacity-class
### Title: Class "summary.superclass.capacity"
### Aliases: summary.superclass.capacity-class
### Keywords: classes

### ** Examples

## a capacity
mu <- capacity(c(0:13,13,13)/13)
a <- Mobius(mu)

## its summary
summary(mu)
summary(a)



cleanEx()
nameEx("summary.superclass.set.func-class")
### * summary.superclass.set.func-class

flush(stderr()); flush(stdout())

### Name: summary.superclass.set.func-class
### Title: Class "summary.superclass.set.func"
### Aliases: summary.superclass.set.func-class
### Keywords: classes

### ** Examples

## a capacity
mu <- set.func(c(0:13,13,13)/13)
a <- Mobius(mu)

## its summary
summary(mu)
summary(a)



cleanEx()
nameEx("superclass.capacity-class")
### * superclass.capacity-class

flush(stderr()); flush(stdout())

### Name: superclass.capacity-class
### Title: Class "superclass.capacity"
### Aliases: superclass.capacity-class
### Keywords: classes

### ** Examples

## three capacities
mu1 <- uniform.capacity(5)
mu2 <- capacity(c(0,0,0:13))
a <- Mobius(mu2)

## compute indices to summarize them
summary(mu1)
summary(mu2)
summary(a)



cleanEx()
nameEx("superclass.set.func-class")
### * superclass.set.func-class

flush(stderr()); flush(stdout())

### Name: superclass.set.func-class
### Title: Class "superclass.set.func"
### Aliases: superclass.set.func-class
### Keywords: classes

### ** Examples

## three set functions
mu1 <- card.set.func(-2:4)
mu2 <- set.func(c(-2,-2,-2:11/11))
mu3 <- Mobius.set.func(c(-7:6,6,6),4,4)

## print mu1
show(mu1)
## the same
mu1
## the others
mu2
mu3

## round mu2
rnd(mu2,2)

## compute indices to summarize them
summary(mu1)
summary(mu2)
summary(mu3)



cleanEx()
nameEx("to.data.frame-methods")
### * to.data.frame-methods

flush(stderr()); flush(stdout())

### Name: to.data.frame-methods
### Title: Puts a set function under the form of a data.frame
### Aliases: to.data.frame to.data.frame-methods
###   to.data.frame,Mobius.card.set.func-method
###   to.data.frame,Mobius.set.func-method
###   to.data.frame,card.set.func-method to.data.frame,set.func-method
### Keywords: methods

### ** Examples

## the Mobius representation of set function
a <- Mobius.set.func(-1:-16,4,4)

## to data.frame
d <- to.data.frame(a)
write.table(d,"my.set.func.csv",sep="\t") 



cleanEx()
nameEx("variance-methods")
### * variance-methods

flush(stderr()); flush(stdout())

### Name: variance-methods
### Title: Normalized variance of a capacity
### Aliases: variance variance-methods variance,Mobius.capacity-method
###   variance,capacity-method variance,card.capacity-method
### Keywords: methods

### ** Examples

## a capacity
mu <- capacity(c(0,0,0,0:12)/12)

## its Mobius transform
a <- Mobius(mu)

## their normalized variance
variance(mu)
variance(a)

## similar examples with card.capacity objects
mu <- lower.capacity(4)
variance(mu)
mu <- uniform.capacity(4)
variance(mu)



cleanEx()
nameEx("veto-methods")
### * veto-methods

flush(stderr()); flush(stdout())

### Name: veto-methods
### Title: Veto indices
### Aliases: veto veto-methods veto,Mobius.capacity-method
###   veto,capacity-method veto,card.capacity-method
### Keywords: methods

### ** Examples

## a capacity
mu <- capacity(c(0:13,13,13)/13)

## its Mobius transform
a <- Mobius(mu)

## their veto indices
veto(mu)
veto(a)

## the same with a card.capacity object
mu <- lower.capacity(4)
veto(mu)



cleanEx()
nameEx("zeta-methods")
### * zeta-methods

flush(stderr()); flush(stdout())

### Name: zeta-methods
### Title: The zeta transform
### Aliases: zeta zeta-methods zeta,Mobius.capacity-method
###   zeta,Mobius.card.set.func-method zeta,Mobius.game-method
###   zeta,Mobius.set.func-method
### Keywords: methods

### ** Examples

## the Mobius transform of a capacity
a <- Mobius.capacity(c(rep(0,15),1),4,4)
a

## its zeta transform
zeta(a)

## let us check
Mobius(zeta(a))

## a similar example with a Mobius.card.set.func object
mu <- card.set.func(-3:4)
a <- Mobius(mu)
zeta(a)



### * <FOOTER>
###
cat("Time elapsed: ", proc.time() - get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
