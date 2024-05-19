t <- seq( 2*pi, 3.5*pi, length.out=1000 )
plot( t*cos(t), t*sin(t) )

t <- runif( 1000, 2*pi, 3.5*pi )
x <- t*cos(t) + rnorm(  1000, 0, 1 )
y <- t*sin(t) + rnorm(  1000, 0, 1 )
plot( x, y, asp=1 )

pc <- principal_curve( cbind( x, y ) )
str(pc)

tibble( x, y, lambda=pc$lambda ) %>%
mutate( idx=row_number() ) %>%
ggplot + geom_point( aes( x=x, y=y, col=lambda ) )

plot( x, y, asp=1, col="gray" )
points( pc$s, col="red" )

plot( pc$lambda, x )

library(splines)
# a grid of values, spanning from 0 to 1
tg <- seq(0,1,length.out=1000)
# the spline basis, with 6 cloumns
b <- bs( tg, intercept=TRUE, df=6 )
# the 6 coumns
plot( NULL, xlim=c(0,1), ylim=c(0,1) )
for( i in 1:6 )
   lines( tg, b[,i] )
# Making a new curve as linear combination of basis splines
plot( tg, b %*% c( 0, 1, .5, 1, 1.2, 0 ) )

plot( pc$lambda, x )
fit <- lm( x ~ bs( pc$lambda, intercept=TRUE, df=6 ) + 0 )
coef(fit)
points( pc$lambda, fitted.values(fit), col="red" )
