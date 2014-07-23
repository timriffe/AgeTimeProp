library(compiler)

# help functions for diagram of lifeline
# clockwise quarter arc (90 degrees)
degrees2radians <- function(degrees){
    degrees * (pi / 180)
}

quarterArc <- function(x, y, radius = 1, fromDegrees = 180, ...){
    xx <- degrees2radians(seq(fromDegrees, fromDegrees + 90, by = .5))
    x <- cos(xx) * radius + x
    y <- sin(xx) * radius + y
    lines(x, y, ...)
}

curlyBrace1 <- function(xl, y, length = 5, radius1 = .5, radius2 = .25, top = TRUE, ...){  
    # i.e. the pointy part on top or on bottom?
    if (top){
        quarterArc(xl + radius1, y - radius1, radius = radius1, fromDegrees = 90, ...)
        quarterArc(xl + length - radius1, y - radius1 , radius = radius1, fromDegrees = 0, ...)
        quarterArc(xl + length / 2 - radius2, y + radius2, radius = radius2, fromDegrees = 270, ...)
        quarterArc(xl + length / 2 + radius2, y + radius2, radius = radius2, fromDegrees = 180, ...)
    } else {
        quarterArc(xl + radius1, y + radius1, radius = radius1, fromDegrees = 180, ...)
        quarterArc(xl + length - radius1, y + radius1 , radius = radius1, fromDegrees = 0 - 90, ...)
        quarterArc(xl + length / 2 - radius2, y - radius2, radius = radius2, fromDegrees = 270 + 90, ...)
        quarterArc(xl + length / 2 + radius2, y - radius2, radius = radius2, fromDegrees = 180 - 90, ...)        
    }
    segments(xl + radius1, y, xl + length / 2 - radius2, y, ...)
    segments(xl + length - radius1, y, xl + length / 2 + radius2, y, ...)   
}

# -----------------------------------------------------
# for the sake of matrix speculation:
Mna0 <- cmpfun(function(M){
            M[is.na(M)]  <- 0
            M[is.nan(M)] <- 0
            M
        })

Minf0 <- cmpfun(function(M){
            M[is.infinite(M)]  <- 0
            M
        })
MinfNA <- cmpfun(function(M){
            M[is.infinite(M)]  <- NA
            M
        })
Thano <- cmpfun(function(Px, dx, stagger = TRUE){
            Np <- length(Px)
            Nd <- length(dx)
            if (Np != Nd){
                N <- max(Np, Nd)
                Px <- c(Px, rep(0, N - Np))
                dx <- c(dx, rep(0, N - Nd))
            } else {
                N <- Np
            }
            
            ay      <- 1:N - 1
            
            dx      <- Mna0(dx)   # remove NAs if any       
            dx      <- c(dx, dx * 0) / sum(dx) # pad out with 0s
            EDx     <- matrix(dx[col(matrix(nrow = N, 
                                            ncol = N)) + ay], 
                    nrow = N, 
                    ncol = N, 
                    dimnames = list(Ex = ay, 
                            Age = ay)
            )
            if (stagger){
                EDx <- (EDx + cbind(EDx[, 2:ncol(EDx)], 0)) / 2
            }
            t(Px * Minf0(Mna0(EDx / rowSums(EDx))))
        })

# this construction is rough, as things need to be staggered still,
# as with the Leslie matrix.
ThanoProjMatrix <- cmpfun(function(Fy, da, lambda = .9){
            N       <- length(Fy)
            # discount for part of infant mortality not surviving until end of year
            da[1]   <- da[1] * (1-lambda)
            
            # NxN matrix
            # fertility component
            Y       <- outer(da, Fy, "*") + 
                    # add survival element-wise
                    rbind(cbind(0,diag(N - 1)), 0)
            
            # reduce e0 fertility by 1/2, as only exposed for part of year
            Y[, 1]  <- Y[, 1] / 2
            # do not allow for Inf or NA values: impute 0
            Y       <- Mna0(Minf0(Y))
            # return projection matrix
            Y
        })

Leslie <- function(Fx, Lx){
    NN     <- length(Lx)
    N      <- NN - 1
    lambda <- -diff(Lx)/2
    Fx     <- (Fx[1:N] + Fx[-1])/2 * (1-lambda) 
    Sx     <- Lx[2:NN] / Lx[1:N]
    Mna0(cbind(rbind(Fx, diag(Sx)), 0))
}

st <- function(x,tot = 1){
    (x / sum(x, na.rm = TRUE)) * tot
}

# modified from popbio, originally from James Holland Jones DemogR package.
# Changed to continue working if solve() returns a singularity error
# by using ginv() from MASS. Possibly bad, but I ran into computational singularity with
# some pretty standard Leslie matrices with contemporary data. Can be toggled off
eigen.analysis <- function(A, zero=TRUE, ginv.if.singular = TRUE){
    ev <- eigen(A)
    # R sorts eigenvalues in decreasing order, according to Mod(values)
    #  ususally dominant eigenvalue is first (ev$values[1]), except for imprimitive matrices with d eigenvalues of equal modulus
    # this should work for most cases
    lmax <- which.max(Re(ev$values))
    lambda <- Re(ev$values[lmax])
    ## Damping ratio. Use second eigenvalue
    # dr<- lambda/abs(ev$values[2])
    ## OR  second largest magnitude in case of ties using rle - round needed for imprimitive matrices
    dr<-rle(round(Mod(ev$values), 5  ))$values
    dr<-dr[1]/dr[2]
    
    W <- ev$vectors
    w <- abs(Re(W[, lmax]))
    ## check if matrix is singular-and output NAs rather than stop (better for loops and bootstrapping)
    V <- try(Conj(solve(W)), silent=TRUE)
    if (class(V) == "try-error" & ginv.if.singular) {
        cat("\n!!!!!!!!!!\nWarning: ginv() used because solve() returned error\nThis could give erroneous results!\n")
        # this is the difference
        V <- Conj(MASS::ginv(W))
        
        
    }
    v <- abs(Re(V[lmax, ]))
    s <- v %o% w
    if (zero) {
        s[A == 0] <- 0
    }
    e <- s * A/lambda
    
    x <- dimnames(A)
    dimnames(s) <- x
    names(w) <- x[[1]]
    names(v) <- x[[1]]
    eigen.analysis <- list(lambda1 = lambda, stable.stage = w/sum(w), 
            sensitivities = s, elasticities = e, repro.value = v/v[1], 
            damping.ratio = dr)
    
    eigen.analysis
}