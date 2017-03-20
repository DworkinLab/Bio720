# Source file for BIO708, multivariate linear models tutorials.


ShapeScore <- function(Beta, Y) {
    ## This computes the "shape score" of Drake and Klingenberg 2008
    ## A simple projection of the shape data in the direction of the vector of regression coefficients (best fit line)
    ## Beta is the coefficients for each procrustes coordinate from a multivariate multiple regression. If entering from coef from lm (as a vector!) it needs to be transposed in the call.
    ## Y is the raw procrustes coordinates
    s = Y %*% t(Beta) %*% ((Beta %*% t(Beta))^-0.5)
    return(shapeScore=s)
}


# A method to compute a coefficient of determination for a multivariate model. See Pitchers et al. 2013 (Evolution)
shapeRsq <- function( model ){
    fitted.variance <- sum(diag(var(model$fitted)))
    total.variance	<- sum(diag(var(model$fitted + model$resid)))
    fitted.variance / total.variance
}

comment(shapeRsq) <- "this function takes a multivariate model object and returns the fitted variance / total variance: equivalent to an Rsquared"


shapePRsq <- function( model ){
    # Based on the derivation from page 269 of Kutner et. al. (Applied linear statistical models edition 5.)
    residual.variance <- var(model$resid)
    variables <- attr(terms(model), "term.labels")
    model.length <- length(variables)
    variable.name <- rep(NA, model.length )
    partial.Rsq <- rep(NA, model.length )
    
    for (i in 1:model.length){
        variable.name[i] <- variables[i]
        drop <- parse( text=variables[i] )
        new.formula <- as.formula( paste( ".~.-", variables[i], sep=""))
        new.model <- update(model, new.formula )
        partial.Rsq[i] <- (sum ( diag( var(new.model$resid))) - sum( diag( residual.variance)) ) / sum( diag( var(new.model$resid)))
    }
    R2 <- shapeRsq( model )
    list(Rsquared=R2, partials=data.frame( cbind( variable.name, partial.Rsq ))	)
}

# calculate tangent approximaton for tangent approximates Procrustes Distance (Euclidean Distance)
# This is just the magnitude of the vector!

PD <- function( x ) {
    sqrt( t( x ) %*% x )}
comment(PD) <- c("This just computes the Euclidean Distance (norm) for a vector")



####### When the vectors can be of arbitrary sign, use this which computes the magnitude of the vector correlation, and then computes the angle.

ang.vec.abs <- function(vec1, vec2){
    vec.cor <- abs((t(vec1) %*% vec2)/(PD(vec1)*PD(vec2)))
    vec.angle <- acos(vec.cor)*(180/pi)
    return(c(vector.cor=vec.cor, vec.angle=vec.angle))}
comment(ang.vec.abs) <- c(" This computes both the vector correlation, and angle, between two vectors.", " to compare to the Pearson correlation coefficient make sure to center and standardize vectors", "set it up to compute the absolute values of the vector correlation")

#################################
