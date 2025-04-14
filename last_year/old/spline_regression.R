b <- function( t, i, p, knots ) {
  if( p==0 )
    as.numeric( t >= knots[i] & t < knots[i+1] )
  else
    ( t - knots[i] ) / ( knots[i+p]-knots[i] ) * b( t, i, p-1, knots ) +
    ( knots[i+p+1] - t ) / ( knots[i+p+1]-knots[i+1] ) * b( t, i+1, p-1, knots )
}    

knots <- seq( 0, 1, length.out=10 )
tg <- seq( 0, 1, length.out=1000 )


m <- sapply( 1:(length(knots)-1), function(i) b( tg, i, 0, knots ) )
matplot( tg, m, type="l" )

m <- sapply( 1:(length(knots)-2), function(i) b( tg, i, 1, knots ) )
matplot( tg, m, type="l" )
lines( tg, rowSums(m), lwd=3, col="lightblue" )
abline( v=knots, lty="dashed", col="gray80" )

m <- sapply( 1:(length(knots)-3), function(i) b( tg, i, 2, knots ) )
matplot( tg, m, type="l", ylim=c(0,1) )
lines( tg, rowSums(m), lwd=3, col="lightblue" )
abline( v=knots, lty="dashed", col="gray80" )

m <- sapply( 1:(length(knots)-4), function(i) b( tg, i, 3, knots ) )
matplot( tg, m, type="l", ylim=c(0,1) )
lines( tg, rowSums(m), lwd=3, col="lightblue" )
abline( v=knots, lty="dashed", col="gray80" )


bd <- function( t, i, p, knots ) {
  if( p==0 )
    0
  elsif( p==1 )
  else
    ( t - knots[i] ) / ( knots[i+p]-knots[i] ) * bd( t, i, p-1, knots ) -
    ( knots[i+p+1] - t ) / ( knots[i+p+1]-knots[i+1] ) * bd( t, i+1, p-1, knots )
}    

bdd <- function( t, i, p, knots ) {
  if( p==0 )
    0
  else
    ( t - knots[i] ) / ( knots[i+p]-knots[i] ) * bdd( t, i, p-1, knots ) -
    ( knots[i+p+1] - t ) / ( knots[i+p+1]-knots[i+1] ) * bdd( t, i+1, p-1, knots )
}    







t <- seq( -1, 5, length.out=1000 )

f1 <- function(t) case_when( 
   t < 0  ~  0,
   t < 1   ~  t,
   t < 2   ~  2-t,
   TRUE    ~ 0 )

plot( t, f1(t), type="l", col="blue" )
lines( t, f1(t-1), type="l", col="red" )

f2 <- function(t)
  t/3 * f1(t) + (3-t)/3 * f1(t-1)

plot( t, f2(t), type="l" )

plot( t, f2(t), type="l", col="blue" )
lines( t, f2(t-1), type="l", col="red" )

f3 <- function(t)
  t/4 * f2(t) + (4-t)/4 * f2(t-1)
plot( t, f3(t), type="l" )

plot( t, ( f3(t) - f3(lag(t)) ) / (t-lag(t)), type="l" )
plot( t, ( f3(lead(t)) - 2*f3(t) + f3(lag(t)) ) / (t-lag(t))^2, type="l" )

f3tt <- function(t) case_when(
  t < 0   ~  0,
  t < 1   ~  0.5*t,
  t < 2   ~  2 - 1.5*t,
  t < 3   ~  1.5*t - 4,
  t < 4   ~  2 - 0.5*t,
  TRUE    ~  0 )
  
plot( t, ( f3(lead(t)) - 2*f3(t) + f3(lag(t)) ) / (t-lag(t))^2, type="l" )
lines( t, f3tt(t), col="red", lty="dashed" )

ioanide

