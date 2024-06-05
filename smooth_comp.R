# Comparison of four types of smoothing splines

# Make example data
n <- 15
x <- runif( n, 0, 10 )
y <- sin(x) + rnorm( n, sd=.1 )

# Make grid to plot smooth curves
xg <- seq( 0, 10, length.out=1000 )

# Determine knot positions:
#   The knots should have the same range as the data x values
#   and the number of data points between knots should be similar
nknots <- 10
stopifnot( knots < n )
f <- approxfun( 1:length(x), sort(x) )
knots <- f( seq( 1, length(x), length.out=nknots ) )  

# Make plot to justify knot positions
plot( 1:length(x), sort(x) )
lines( seq( 0,length(x), length.out=1000), f(seq( 0,length(x), length.out=1000)), type="l" )
abline( v = seq( 1, length(x), length.out=nknots ), col="gray" )
abline( h = knots, col="gray" )

# Plot data and knot positions
plot( x, y, xlim=c(-1, 11), ylim=c(-1.3,1.3) )
abline( v = knots, col="gray" )
abline( v = range(x) )

# With bs
mysp <- function(x) bs( x, intercept=TRUE, 
    Boundary.knots=knots[c(1,length(knots))], knots=knots[2:(length(knots)-1)] )
fit <- lm.fit( mysp(x), y )
yg <- mysp(xg) %*% fit$coefficients
lines( xg, yg, col="pink" )

# With ns (like bs, but zero curvature at boundaries)
mysp <- function(x) ns( x, intercept=TRUE, 
    Boundary.knots=knots[c(1,length(knots))], knots=knots[2:(length(knots)-1)] )
fit <- lm.fit( mysp(x), y )
yg <- mysp(xg) %*% fit$coefficients
lines( xg, yg, col="darkgreen" )

# With smooth.spline (like bs, but with curvature penalty)
fit <- smooth.spline( x, y, df=df )
yg <- predict( fit, xg )$y
lines( xg, yg, col="blue" )

# With P-splines
m <- 50
mysp <- function(xx) bs( xx, intercept=TRUE, 
    Boundary.knots=range(x), knots=seq(min(x),max(x),length.out=m-4) )

pty <- diag(2,m)
pty[ row(pty) == col(pty)+1 ] <- -1
pty[ row(pty) == col(pty)-1 ] <- -1
pty[1,1] <- 1
pty[m,m] <- 1

lambda_scale <- sum(diag( t(mysp(x)) %*% mysp(x) )) / sum(diag(pty))
lambda <- 3 * lambda_scale   # <- needs to be optimized
beta <- solve( t(mysp(x)) %*% mysp(x) + lambda * pty , t(mysp(x)) %*% y )
lines( xg, mysp(xg) %*% beta, col="brown" )
