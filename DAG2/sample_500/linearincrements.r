# This is a function for fitting the various linear increments methods.  It
# is designed to be for public distribution.  It is a tidied up/improved
# version of the linearincr() function in the file linearincr_function.r.
# It is also designed (unlike linearincr()) to deal with situation where
# data can be missing at time 1.


library(norm) # This is for em.norm and imp.norm functions

expit <- function(x) { return(exp(x) / (1+exp(x))) }
logit <- function(x) { return( log(x / (1-x)) ) }

# The function decimal2binary has been copied from the `GA' R library in CRAN.

decimal2binary <- function (x, length) 
{
  x <- as.integer(x)
  b <- if (missing(length)) 
    NULL
  else rep(0, length)
  i <- 0
  while (x >= 1) {
    i <- i + 1
    b[i] <- x%%2
    x <- x%/%2
  }
  return(rev(b))
}


linearincrements <- function(y, x=NULL, alive=NULL, autoreg=1, method="LI-LS", makemonotone=F, delete.partmiss=T, maxiter.em=1000, toler.emnr=1e-5, returnimps=F, printsigmas=F, verbose=F)
  {
# This function applies one of eight methods: "compensator" (estimating the compensator);
# "LI-LS" (LI-LS imputation); "LI-uMVN" (LI-uMVN imputation); "LI-aMVN" (LI-aMVN imputation);
# "LI-rMVN" (LI-rMVN imputation); "uMVN" (uMVN imputation); "aMVN" (aMVN imputation); and
# "rMVN" (rMVN imputation).
#
# autoreg determines the order of the autoregression used in the LI model.  If method equals
# LI-rMVN imputation or rMVN imputation, this argument is ignored and autoreg is set to zero.
# The LI-aMVN imputation or aMVN imputation methods have only been implemented with autoreg=1.
#   
# maxiter.em is maximum number of iterations in my EM/Newton Raphson algorithm used for
# LI-aMVN imputation (EM algorithm), aMVN imputation (EM algorithm), LI-rMVN imputation (NR
# algorithm) or rMVN imputation (NR algorithm).
#
# toler.em is tolerance parameter in my EM/NR algorithm.
#
# N is number of individuals
# M is number of timepoints
#
# If makemonotone==T, drop all measurements after the first dropout.  This is
# not the default.
#
# If delete.partmiss==T (which is the default) and autoreg>1, delete outcomes
# that are observed but that do not belong to any fully-observed vectors of
# outcomes, and hence are treated as missing.  By choosing delete.partmiss==F,
# this deletion of outcomes that are technically missing is not performed and
# these outcomes are retained.
#    
# If returnimps==T, the imputed values are returned
#
# If printsigmas==T, print the MLEs of the variance matrices of the MVN models

# Note that beta(2:(autoreg+1),] are the slopes for regression of (Y_t - Y_{t-1}) on
# Y_{t-autoreg}, ..., Y_{t-1}, rather than Y_t on Y_{t-autoreg}, ..., Y_{t-1}
    
    N <- nrow(y)
    M <- ncol(y)
    if (is.null(x))
      dimx <- 0 else
    {x <- as.matrix(x)
     dimx <- ncol(x)
   }

    if (delete.partmiss & autoreg>1)
      {# r.vec[i,j] will be an indicator that equals one if and only if
       # individual i's jth, (j-1)th, ... and (j-autoreg+1)th outcomes y are
       # all observed.
       r.vec <- !is.na(y)
       obs.before <- sum(!is.na(y))
       for (j in 2:M)
         {furthest <- max(j-autoreg+1, 1)
           for (k in furthest:(j-1))
             r.vec[is.na(y[,k]), j] <- F
        }
       # Now delete outcome y[i,j] if none of the indicators r.vec[i,j],
       # r.vec[i,j+1], ... , r.vec[i,j+autoreg-1] equal true.
       deletey <- matrix(T, N, M)
       deletey[, 1] <- F
       for (j in 2:M)
         {furthest <- min(j+autoreg-1, M)
           for (k in j:furthest)
             deletey[r.vec[,k], j] <- F
        }
       y[deletey==T] <- NA
       obs.after <- sum(!is.na(y))
       print(paste("Dropping", obs.before - obs.after, "observations.  Use delete.partmiss=F to prevent this."))
     }
    
    r <- !is.na(y)

    if (method=="LI-rMVN" | method=="rMVN")
      autoreg <- 0

    if (( method=="LI-aMVN" | method=="aMVN") & autoreg!=1)
      stop("ERROR: LI-aMVN and aMVN methods only currently implemented for autoreg=1")

    autoreg1 <- max(autoreg, 1)

    if (method!="compensator" & method!="LI-LS" & method!="LI-uMVN" & method!="LI-aMVN" & method!="LI-rMVN" & method!="uMVN" & method!="aMVN" & method!="rMVN")
      stop("ERROR: Invalid method.  Valid methods are compensator, LI-LS, LI-uMVN, LI-aMVN, LI-rMVN, uMVN, aMVN and rMVN")

    if (!is.null(x))
      if (sum(is.na(x))>0)
        stop("ERROR: This function cannot handle missing values in the covariates X")
    
    if (makemonotone)
      for (j in 2:M)
        r[ r[,j-1]==0, j ] <- 0

    iter <- 0
# iter will be used by the EM algorithm and Newton-Raphson algorithm only, but is returned
# by this function even if it is not used.
    
# Work out which method is being used to estimate the betas.
    if (method=="compensator" | method=="LI-LS")
      beta.method <- "LS"
    if (method=="LI-uMVN" | method=="uMVN")
      beta.method <- "uMVN"
    if (method=="LI-aMVN" | method=="aMVN")
      beta.method <- "aMVN"
    if (method=="LI-rMVN" | method=="rMVN")
      beta.method <- "rMVN"


# If beta.method==LS, estimate beta in the Aalen and Gunnes way.
    
    if (beta.method=="LS")
      {
# Calculate the missingness indicator r.ywhole that will be used by the Aalen
# and Gunnes method for estimating the parameters of the linear increments
# model.  It equals 1 only if the current outcome and the previous autoreg
# outcomes are all observed.
        r.ywhole <- r
        for (j in 2:M)
          {furthest <- max(j-autoreg, 1)
           for (k in furthest:(j-1))
             r.ywhole[r[,k]==0, j] <- 0
         }
        
        beta.ag <- matrix(0, 1+autoreg1+dimx, M)
# Each column of beta.ag is for one of the M timepoints.  The first row of
# beta.ag is for the intercept (called alpha in our article).  The next
# autoreg rows are for the coefficients of the previous autoreg outcomes
# (called beta in our article).  The final rows are for the coefficients of
# the covariates x (called gamma in our article).  Beta.uc and beta.ar are
# analogous.

# Estimate beta.ag for first timepoint.
        if (dimx>0)
          beta.ag[c(1, (1+autoreg1)+1:dimx), 1] <- lm(y[,1] ~ x, subset=r.ywhole[,1]==1)$coefficients
        else
          beta.ag[1, 1] <- mean(y[r.ywhole[,1]==1, 1])
        
# Estimate beta.ag for second to Mth timepoints.
        for (j in 2:M)
          {if (autoreg>0)
             {furthest <- max(j-autoreg, 1)
# furthest is the time of the earliest outcome on which the outcome at time j
# will be regressed.  It will be regressed on the oucomes at time furthest,
# furthest+1, ... , j-1 and the covariates x.
              covars <- cbind(y[, furthest:(j-1)], x)
              firstindex <- max(1, autoreg-j+2)
# firstindex equals 1 if the number of earlier outcomes on which the outcome
# at time j is regressed is the maximum number (i.e. autoreg).  In this case,
# all elements of this column of beta.ag will need to be filled.  If
# firstindex>1 (meaning that the outcome at time j is regressed on fewer than
# autoreg outcomes), not all elements of this column of beta.ag will need to be
# filled.
              beta.ag[c(1, (firstindex+1):(1+autoreg+dimx)), j] <- lm(y[,j] - y[,j-1] ~ covars, subset=r.ywhole[,j]==1)$coefficients
            } else
           {# next bit is for autoreg=0 case
             if (dimx>0)  
               {covars <- x
                beta.ag[c(1,3:(2+dimx)), j] <- lm(y[,j] - y[,j-1] ~ covars, subset=r[,j]==1 & r[,j-1]==1)$coefficients
              } else
             beta.ag[1, j] <- lm(y[,j] - y[,j-1] ~ 1, subset=r[,j]==1 & r[,j-1]==1)$coefficients
           } # that is end of autoreg=0 case
         } # end of for (j in 2:M) loop
      }
    
# That is the end of the calculation of beta.ag

# Now, if beta.method=="uMVN", calculate beta using the unconstrained MVN method.  The
# beginning of the calculation will also be done if beta.method=="aMVN", because the
# estimates from the unconstrained method are used as starting values in my EM algorithm
# for the autoregressive MVN method.
    
    if (beta.method=="uMVN" | beta.method=="aMVN")
      {yx <- cbind(y, x)
       prelim <- prelim.norm(yx)
       mleobj.uc <- em.norm(prelim, showits=F)
       mlemean.uc <- getparam.norm(prelim, mleobj.uc)$mu
       mlevar.uc <- getparam.norm(prelim, mleobj.uc)$sigma
     }

    if (beta.method=="uMVN")
      {beta.uc <- matrix(0, 1+autoreg1+dimx, M)
        if (dimx==0)
          {beta.uc[1, 1] <- mlemean.uc[1]
           for (j in 2:M)
            {furthest <- max(j-autoreg, 1)
             firstindex <- max(1, autoreg-j+2)
             beta.uc[(firstindex:autoreg)+1, j] <- solve(mlevar.uc[furthest:(j-1), furthest:(j-1)]) %*% mlevar.uc[furthest:(j-1), j]
             beta.uc[1, j] <- mlemean.uc[j] - t(beta.uc[(firstindex:autoreg)+1, j]) %*% mlemean.uc[furthest:(j-1)]
           }
         } else
        {invcovx <- solve(mlevar.uc[M+1:dimx, M+1:dimx])
         beta.uc[1+autoreg+1:dimx, 1] <- invcovx %*% mlevar.uc[M+1:dimx, 1]
         beta.uc[1, 1] <- mlemean.uc[1] - t(beta.uc[1+autoreg+1:dimx, 1]) %*% mlemean.uc[M+1:dimx]
         for (j in 2:M)
           {furthest <- max(j-autoreg, 1)
            firstindex <- max(1, autoreg-j+2)
            beta.uc[(firstindex:autoreg)+1, j] <- solve( mlevar.uc[furthest:(j-1), furthest:(j-1)] - mlevar.uc[furthest:(j-1), M+1:dimx] %*% invcovx %*% mlevar.uc[M+1:dimx, furthest:(j-1)] ) %*% ( mlevar.uc[furthest:(j-1), j] - mlevar.uc[furthest:(j-1), M+1:dimx] %*% invcovx %*% mlevar.uc[M+1:dimx, j] )
            if (furthest!=j-1)
              beta.uc[autoreg+1+1:dimx, j] <- invcovx %*% ( mlevar.uc[M+1:dimx, j] - mlevar.uc[M+1:dimx, furthest:(j-1)] %*% beta.uc[(firstindex:autoreg)+1, j] ) else
           beta.uc[autoreg+1+1:dimx, j] <- invcovx %*% ( mlevar.uc[M+1:dimx, j] - mlevar.uc[M+1:dimx, furthest] * beta.uc[(firstindex:autoreg)+1, j] )
          }
        for (j in 2:M)
          {furthest <- max(j-autoreg, 1)
           firstindex <- max(1, autoreg-j+2)
           beta.uc[1, j] <- mlemean.uc[j] - t(beta.uc[(firstindex:autoreg)+1, j]) %*% mlemean.uc[furthest:(j-1)] - mlemean.uc[M+1:dimx] %*% beta.uc[autoreg+1+1:dimx, j]
         }
       }
    
# But beta[2:(autoreg+1),] are the slopes for regression of Y_t on
# Y_{t-autoreg}, ..., Y_{t-1}, rather than (Y_t - Y_{t-1}) on
# Y_{t-autoreg}, ..., Y_{t-1}, so make it for (Y_t - Y_{t-1}) on
# Y_{t-autoreg}, ..., Y_{t-1} by subtracting 1
        beta.uc[autoreg+1, 2:M] <- beta.uc[autoreg+1, 2:M] - 1
      }

# If beta.method=="uMVN" or "aMVN", then beta.uc has now been calculated

# Now, if beta.method="aMVN", then estimate beta in the autoregressive MVN using my EM
# algorithm.

# THE AUTOREGRESSIVE METHOD IS NOT IMPLEMENTED FOR AUTOREG>1.

    if (beta.method=="aMVN")
      {mlemean.ar <- mlemean.uc
       mlevar.ar <- mlevar.uc

       flag <- 0
       iter <- 0
       
# Start of my EM algorithm for autoregressive model
       
       while (!flag & iter < maxiter.em)
         {iter <- iter+1
          currentval <- makeparam.norm(prelim, thetalist=list(mlemean.ar, mlevar.ar))
          mlemean.ar0 <- mlemean.ar
          mlevar.ar0 <- mlevar.ar
          mleobj.ar <- em.norm(prelim, maxits=1, start=currentval, showits=F)
          mlemean.ar <- getparam.norm(prelim, mleobj.ar)$mu
          mlevar.ar <- getparam.norm(prelim, mleobj.ar)$sigma
          beta.ar <- matrix(0, 2+dimx, M)
           
          if (dimx==0)
            {for (j in 2:M)
               beta.ar[2,j] <- mlevar.ar[j, j-1] / mlevar.ar[j-1, j-1]
             beta.ar[1, 1] <- mlemean.ar[1]
             beta.ar[1, 2:M] <- mlemean.ar[2:M] - beta.ar[2, 2:M] * mlemean.ar[1:(M-1)]
             for (j in 1:(M-1))
               for (k in (j+1):M)
                 mlevar.ar[j, k] <- prod(beta.ar[2, (j+1):k]) * mlevar.ar[j,j]
             for (j in 1:(M-1))
               for (k in (j+1):M)
                 mlevar.ar[k,j] <- mlevar.ar[j,k]
           } else
           {# Now to do these calculations  when there are covariates
             covx <- mlevar.ar[M+1:dimx, M+1:dimx]
             invcovx <- solve(covx)
             beta.ar[2+1:dimx, 1] <- invcovx %*% mlevar.ar[M+1:dimx, 1]
             beta.ar[1, 1] <- mlemean.ar[1] - beta.ar[2+1:dimx, 1] %*% mlemean.ar[M+1:dimx]
             for (j in 2:M)
               {beta.ar[2, j] <- ( mlevar.ar[j-1, j] - mlevar.ar[j-1, M+1:dimx] %*% invcovx %*% mlevar.ar[M+1:dimx, j] ) / ( mlevar.ar[j-1, j-1] - mlevar.ar[j-1, M+1:dimx] %*% invcovx %*% mlevar.ar[M+1:dimx, j-1] )
                beta.ar[2+1:dimx, j] <- invcovx %*% ( mlevar.ar[M+1:dimx, j] - mlevar.ar[M+1:dimx, j-1] * beta.ar[2,j] )
              }
             for (j in 2:M)
               beta.ar[1, j] <- mlemean.ar[j] - beta.ar[2, j] * mlemean.ar[j-1] - mlemean.ar[M+1:dimx] %*% beta.ar[2+1:dimx, j]
             
             for (j in 1:(M-1))
               for (k in (j+1):M)
                 mlevar.ar[j, k] <- mlevar.ar[k, M+1:dimx] %*% invcovx %*% mlevar.ar[M+1:dimx, j] + prod(beta.ar[2, (j+1):k]) * ( mlevar.ar[j,j] - mlevar.ar[j, M+1:dimx] %*% invcovx %*% mlevar.ar[M+1:dimx, j] )
             
# So far, I have only calculated the upper triangular part of the top-left
# M-by-M submatrix of mlevar.ar.  Now to do the lower triangular part,
# remembering that this is a symmetric matrix.
             for (j in 1:(M-1))
               for (k in (j+1):M)
                 mlevar.ar[k,j] <- mlevar.ar[j,k]
           } # End of calculation when there are covariates
          
          sumdiff <- sum(abs(mlemean.ar - mlemean.ar0)) + sum(abs(mlevar.ar - mlevar.ar0))
          if (sumdiff < toler.emnr)
            flag <- 1
          
# But beta[2,] are the slopes for regression of Y_t on Y_{t-1}, rather than
# (Y_t - Y_{t-1}) on Y_{t-1}, so make it for (Y_t - Y_{t-1}) on Y_{t-1} by
# subtracting 1
          beta.ar[2, 2:M] <- beta.ar[2, 2:M] - 1
# End of EM algorithm for autoregressive model (i.e. autoreg=T)
         }
      }

# If beta.method=="aMVN", then beta.ar has now been calculated using a LI model with
# autoregression of order 1.

# Now, if beta.method=="rMVN", calculate beta.ar using a LI model with a random-walk
# structure (autoregression of order zero).  This is done using a Newton-Raphson algorithm.
    
    if (beta.method=="rMVN")
      {temp <- mvn.nonauto(y, x, verbose=verbose)
       iter <- temp$iter.nr
       theta.nonauto <- temp$theta
      
# Extract estimates from theta.non.auto
       beta.ar <- matrix(0, 2+dimx, M)
# Note that the second column of beta.ar is not used.
       beta.ar[1,] <- theta.nonauto[1:M]
       if (dimx>0)
         for (l in 1:dimx)
           beta.ar[2+l,] <- theta.nonauto[M*l+1:M]
       sigma.ar <- theta.nonauto[M*(dimx+1)+1:M]
      
# Now calculate mlemean.ar
       mlemean.ar <- rep(0, M+dimx)
       if (dimx>0)
         {mlemean.ar[M+1:dimx] <- apply(x, 2, mean)
          mlemean.ar[1] <- beta.ar[1, 1] + beta.ar[2+1:dimx, 1] %*% mlemean.ar[M+1:dimx]
          for (j in 2:M)
            mlemean.ar[j] <- mlemean.ar[j-1] + beta.ar[1, j] + beta.ar[2+1:dimx, j] %*% mlemean.ar[M+1:dimx]
        } else
       {mlemean.ar[1] <- beta.ar[1, 1]
        for (j in 2:M)
          mlemean.ar[j] <- mlemean.ar[j-1] + beta.ar[1, j]
      }

# Now calculate mlevar.ar
       mlevar.ar <- matrix(0, M+dimx, M+dimx)
       if (dimx>0)
         {# First get variance-covariance matrix of covariates x
           mlevar.ar[M+1:dimx, M+1:dimx] <- var(x) * (N-1) / N
# multiplying by (N-1)/N gives the mle, rather than the unbiased estimate that the var() function gives (because it divides by N-1 rather than N).
# Second, get covariance between Yj and covariates x.
           mlevar.ar[1, M+1:dimx] <- beta.ar[2+1:dimx, 1] %*% mlevar.ar[M+1:dimx, M+1:dimx]
           for (j in 2:M)
             mlevar.ar[j, M+1:dimx] <- mlevar.ar[j-1, M+1:dimx] + beta.ar[2+1:dimx, j] %*% mlevar.ar[M+1:dimx, M+1:dimx]
           mlevar.ar[M+1:dimx, 1:M] <- t(mlevar.ar[1:M, M+1:dimx])
# Third, get variance of Yj.
           mlevar.ar[1, 1] <- t(beta.ar[2+1:dimx, 1]) %*% mlevar.ar[M+1:dimx, M+1:dimx] %*% beta.ar[2+1:dimx, 1] + sigma.ar[1]
           for (j in 2:M)
             mlevar.ar[j, j] <- mlevar.ar[j-1, j-1] + t(beta.ar[2+1:dimx, j]) %*% mlevar.ar[M+1:dimx, M+1:dimx] %*% beta.ar[2+1:dimx, j] + 2 * t(beta.ar[2+1:dimx, j]) %*% mlevar.ar[j-1, M+1:dimx] + sigma.ar[j]
         } else
       {# All that is simpler if there are no covariates.
         mlevar.ar[1, 1] <- sigma.ar[1]
         for (j in 2:M)
           mlevar.ar[j, j] <- mlevar.ar[j-1, j-1] + sigma.ar[j]
       }
# Fourth, get covariance between Yj and Yk.
       for (j in 1:(M-1))
         for (k in (j+1):M)
           {if (dimx>0)
              mlevar.ar[j, k] <- mlevar.ar[j, j] + t(mlevar.ar[j, M+1:dimx]) %*% solve(mlevar.ar[M+1:dimx, M+1:dimx]) %*% mlevar.ar[k, M+1:dimx] - t(mlevar.ar[j, M+1:dimx]) %*% solve(mlevar.ar[M+1:dimx, M+1:dimx]) %*% mlevar.ar[M+1:dimx, j] else
            mlevar.ar[j, k] <- mlevar.ar[j, j]
          }
       for (j in 1:(M-1))
         for (k in (j+1):M)
           mlevar.ar[k, j] <- mlevar.ar[j, k]
     }

# If beta.method=="rMVN", then beta.ar has now been calculated using a LI model with
# a random-walk structure.

# LIparams will now contain the betas, i.e. the estimated parameters of the LI model.
    if (beta.method=="compensator" | beta.method=="LS")
      LIparams <- beta.ag
    if (beta.method=="uMVN")
      LIparams <- beta.uc
    if (beta.method=="aMVN" | beta.method=="rMVN")
      LIparams <- beta.ar

    if (dimx==0)
      dimnames(LIparams) <- list(c("alpha", rep("beta", autoreg1)), paste("t", 1:M, sep="")) else
    dimnames(LIparams) <- list(c("alpha", rep("beta", autoreg1), paste("gamma", 1:dimx, sep="")), paste("t", 1:M, sep=""))

# Now, if method=="compensator", then estimate the compensator (this is done using the LS
# estimates of the parameters in the LI model).

    if (method=="compensator")
      {Delta <- matrix(0, N, M)
       Delta[,1] <- rep(LIparams[1,1], N)
       if (dimx>0)
         for (k in 1:dimx)
           Delta[,1] <- Delta[,1] + x[,k] * LIparams[1+autoreg1+k, 1]
       
       for (j in 2:M)
         {furthest <- max(j-autoreg, 1)
          firstindex <- max(1, autoreg-j+2)
          Delta[,j] <- LIparams[1, j] + Delta[, j-1]
          if (autoreg>0)
            for (k in 1:autoreg)
              if (k<=j-1)
                Delta[,j] <- Delta[,j] + LIparams[autoreg+2-k, j] * Delta[, j-k]
          if (dimx>0)
            for (k in 1:dimx)
              Delta[,j] <- Delta[,j] + x[,k] * LIparams[1+autoreg1+k, j]
        }
       meanDelta <- apply(Delta, 2, mean)
       names(meanDelta) <- paste("t", 1:M, sep="")
     }

# End of estimating the compensator.

# Now, if method=="LI-LS", "LI-uMVN", "LI-aMVN" or "LI-rMVN", perform LI imputation in the
# Aalen and Gunnes way using the corresponding estimate (LS, uMVN, aMVN or rMVN) of the
# parameters in the LI model.

    if (method=="LI-LS" | method=="LI-uMVN" | method=="LI-aMVN" | method=="LI-rMVN")
      {yimp <- y
       dimnames(yimp) <- list(1:N, paste("t", 1:M, sep=""))
       
       missing <- !r[,1]
       yimp[missing==T, 1] <- LIparams[1,1]
       if (dimx>0)
         for (k in 1:dimx)
           yimp[missing==T, 1] <- yimp[missing==T, 1] + x[missing==T, k] * LIparams[1+autoreg1+k, 1]

       for (j in 2:M)
         {missing <- !r[,j]
          furthest <- max(j-autoreg, 1)
          firstindex <- max(1, autoreg-j+2)
          yimp[missing==T, j] <- LIparams[1, j] + yimp[missing==T, j-1]
          if (autoreg>0)
            for (k in 1:autoreg)
              if (k<=j-1)
                yimp[missing==T, j] <- yimp[missing==T, j] + yimp[missing==T, j-k] * LIparams[autoreg1+2-k, j]
          if (dimx>0)
            yimp[missing==T, j] <- yimp[missing==T, j] + as.matrix(x[missing==T,]) %*% LIparams[1+autoreg1+1:dimx, j]
        }
       meanimpy <- apply(yimp, 2, mean)
       names(meanimpy) <- paste("t", 1:M, sep="")
       if (!is.null(alive))
         {meanimpy.alive <- apply(yimp*alive, 2, mean) / apply(alive, 2, mean)
          names(meanimpy.alive) <- paste("t", 1:M, sep="")
        }
     }

# That is the end of LI-LS, LI-uMVN, LI-aMVN or LI-rMVN imputation.


# Now, if method=="uMVN", "aMVN" or "rMVN", perform MVN imputation.
    
    if (method=="uMVN" | method=="aMVN" | method=="rMVN")
      {if (method=="uMVN")
         {mlemean <- mlemean.uc
          mlevar <- mlevar.uc
        } else
       {mlemean <- mlemean.ar
        mlevar <- mlevar.ar
      }

       yx <- cbind(y, x)
       yimp <- expect.cond.norm(yx, mlemean, mlevar)[, 1:M]
       meanimpy <- apply(yimp, 2, mean)
       names(meanimpy) <- paste("t", 1:M, sep="")
       if (!is.null(alive))
         {meanimpy.alive <- apply(yimp*alive, 2, mean) / apply(alive, 2, mean)
          names(meanimpy.alive) <- paste("t", 1:M, sep="")
        }
       
       if (printsigmas)
         {print("MLE of variance matrix for MVN model is")
          print(mlevar)
        }
     }
       
# Now return the results
    
    if (method=="compensator")
      return( list(LIparams=LIparams, meanDelta=meanDelta) ) else
    {if (returnimps & !is.null(alive))
       return( list(LIparams=LIparams, meanimpy=meanimpy, meanimpy.alive=meanimpy.alive, numiters.emnr=iter, yimp=yimp) )
     if (returnimps & is.null(alive))
       return( list(LIparams=LIparams, meanimpy=meanimpy, numiters.emnr=iter, yimp=yimp) )
     if (!returnimps & !is.null(alive))
       return( list(LIparams=LIparams, meanimpy=meanimpy, meanimpy.alive=meanimpy.alive, numiters.emnr=iter) )
     if (!returnimps & is.null(alive))
       return( list(LIparams=LIparams, meanimpy=meanimpy, numiters.emnr=iter) )
   }
    
  } # end of linearincrements function




























# The following function linearincrements.doall() has been (largely) tested.  However, it
# applies all methods and I want instead a function that only applies one method at a time.
# Hence this method has been superceded by linearincrements().

linearincrements.doall <- function(y, r, x=NULL, alive=matrix(1, nrow(y), ncol(y)), autoreg=1, makemonotone=F, maxiter.em=1000, toler.emnr=1e-5, returnimps=F, printsigmas=F, verbose=F)
  {# This function applies all the methods.
# If autoreg>0 then fit the autoregressive LI model of order autoreg (i.e.
# estimate the effect of Y_{t-1}, ..., Y_{t-autoreg} on Y_t; if autoreg=F then
# constrain beta to equal one for all methods.  For the constrained MVN
# method, this uses Newton-Raphson (NR) instead of EM.
# maxiter.em is maximum number of iterations in my EM/NR algorithm.
# toler.em is tolerance parameter in my EM/NR algorithm
# N is number of individuals
# M is number of timepoints
# If makemonotone==T, drop all measurements after the first dropout
# If returnimps==T, the imputed values are returned
# If printsigmas==T, print the MLEs of the variance matrices of the MVN models

# Note that beta(2:(autoreg+1),] are the slopes for regression of (Y_t - Y_{t-1}) on
# Y_{t-autoreg}, ..., Y_{t-1}, rather than Y_t on Y_{t-autoreg}, ..., Y_{t-1}
    
    N <- nrow(y)
    M <- ncol(y)
    if (is.null(x))
      dimx <- 0 else
    {x <- as.matrix(x)
     dimx <- ncol(x)
   }
    
    autoreg1 <- max(autoreg, 1)
    
    if (makemonotone)
      for (j in 2:M)
        r[ r[,j-1]==0, j ] <- 0

# Calculate the missingness indicator r.ywhole that will be used by the Aalen
# and Gunnes method for estimating the parameters of the linear increments
# model.  It equals 1 only if the current outcome and the previous autoreg
# outcomes are all observed.
    r.ywhole <- r
    for (j in 2:M)
      {furthest <- max(j-autoreg, 1)
       for (k in furthest:(j-1))
         r.ywhole[r[,k]==0, j] <- 0
     }
        
# Estimate beta in the Aalen and Gunnes way
    beta.ag <- matrix(0, 1+autoreg1+dimx, M)
# Each column of beta.ag is for one of the M timepoints.  The first row of
# beta.ag is for the intercept (called alpha in our article).  The next
# autoreg rows are for the coefficients of the previous autoreg outcomes
# (called beta in our article).  The final rows are for the coefficients of
# the covariates x (called gamma in our article).

# Estimate beta for first timepoint.
    if (dimx>0)
      beta.ag[c(1, (1+autoreg1)+1:dimx), 1] <- lm(y[,1] ~ x, subset=r.ywhole[,1]==1)$coefficients
    else
      beta.ag[1, 1] <- mean(y[r.ywhole[,1]==1, 1])
    
# Estimate beta for second to Mth timepoints.
    for (j in 2:M)
      {if (autoreg>0)
         {furthest <- max(j-autoreg, 1)
# furthest is the time of the earliest outcome on which the outcome at time j
# will be regressed.  It will be regressed on the oucomes at time furthest,
# furthest+1, ... , j-1 and the covariates x.
          covars <- cbind(y[, furthest:(j-1)], x)
          firstindex <- max(1, autoreg-j+2)
# firstindex equals 1 if the number of earlier outcomes on which the outcome
# at time j is regressed is the maximum number (i.e. autoreg).  In this case,
# all elements of this column of beta.ag will need to be filled.  If
# firstindex>1 (meaning that the outcome at time j is regressed on fewer than
# autoreg outcomes), not all elements of this column of beta.ag will need to be
# filled.
          beta.ag[c(1, (firstindex+1):(1+autoreg+dimx)), j] <- lm(y[,j] - y[,j-1] ~ covars, subset=r.ywhole[,j]==1)$coefficients
        } else
       {# next bit is for autoreg=0 case
         if (dimx>0)  
           {covars <- x
            beta.ag[c(1,3:(2+dimx)), j] <- lm(y[,j] - y[,j-1] ~ covars, subset=r[,j]==1 & r[,j-1]==1)$coefficients
          } else
         beta.ag[1, j] <- lm(y[,j] - y[,j-1] ~ 1, subset=r[,j]==1 & r[,j-1]==1)$coefficients
       } # that is end of autoreg=0 case
     } # end of for (j in 2:M) loop
    
# Estimate beta in the unconstrained MVN way
    yx <- cbind(y, x)
    prelim <- prelim.norm(yx)
    mleobj.uc <- em.norm(prelim, showits=F)
    mlemean.uc <- getparam.norm(prelim, mleobj.uc)$mu
    mlevar.uc <- getparam.norm(prelim, mleobj.uc)$sigma
    beta.uc <- matrix(0, 1+autoreg1+dimx, M)
    if (autoreg>0)
      {# If autoreg==0, don't calculate beta.uc at all.
        if (dimx==0)
          {beta.uc[1, 1] <- mlemean.uc[1]
           for (j in 2:M)
            {furthest <- max(j-autoreg, 1)
             firstindex <- max(1, autoreg-j+2)
             beta.uc[(firstindex:autoreg)+1, j] <- solve(mlevar.uc[furthest:(j-1), furthest:(j-1)]) %*% mlevar.uc[furthest:(j-1), j]
             beta.uc[1, j] <- mlemean.uc[j] - t(beta.uc[(firstindex:autoreg)+1, j]) %*% mlemean.uc[furthest:(j-1)]
        }
         } else
        {invcovx <- solve(mlevar.uc[M+1:dimx, M+1:dimx])
         beta.uc[1+autoreg+1:dimx, 1] <- invcovx %*% mlevar.uc[M+1:dimx, 1]
        beta.uc[1, 1] <- mlemean.uc[1] - t(beta.uc[1+autoreg+1:dimx, 1]) %*% mlemean.uc[M+1:dimx]
        for (j in 2:M)
          {furthest <- max(j-autoreg, 1)
           firstindex <- max(1, autoreg-j+2)
           beta.uc[(firstindex:autoreg)+1, j] <- solve( mlevar.uc[furthest:(j-1), furthest:(j-1)] - mlevar.uc[furthest:(j-1), M+1:dimx] %*% invcovx %*% mlevar.uc[M+1:dimx, furthest:(j-1)] ) %*% ( mlevar.uc[furthest:(j-1), j] - mlevar.uc[furthest:(j-1), M+1:dimx] %*% invcovx %*% mlevar.uc[M+1:dimx, j] )
           if (furthest!=j-1)
             beta.uc[autoreg+1+1:dimx, j] <- invcovx %*% ( mlevar.uc[M+1:dimx, j] - mlevar.uc[M+1:dimx, furthest:(j-1)] %*% beta.uc[(firstindex:autoreg)+1, j] ) else
           beta.uc[autoreg+1+1:dimx, j] <- invcovx %*% ( mlevar.uc[M+1:dimx, j] - mlevar.uc[M+1:dimx, furthest] * beta.uc[(firstindex:autoreg)+1, j] )
         }
        for (j in 2:M)
          {furthest <- max(j-autoreg, 1)
           firstindex <- max(1, autoreg-j+2)
           beta.uc[1, j] <- mlemean.uc[j] - t(beta.uc[(firstindex:autoreg)+1, j]) %*% mlemean.uc[furthest:(j-1)] - mlemean.uc[M+1:dimx] %*% beta.uc[autoreg+1+1:dimx, j]
         }
       }
    
# But beta[2:(autoreg+1),] are the slopes for regression of Y_t on
# Y_{t-autoreg}, ..., Y_{t-1}, rather than (Y_t - Y_{t-1}) on
# Y_{t-autoreg}, ..., Y_{t-1}, so make it for (Y_t - Y_{t-1}) on
# Y_{t-autoreg}, ..., Y_{t-1} by subtracting 1
        beta.uc[autoreg+1, 2:M] <- beta.uc[autoreg+1, 2:M] - 1
      } # end of if (autoreg>0) condition
    
# Estimate beta in the autoregressive/constrained MVN way using my EM algorithm
# (for autoreg=1) or NR (for autoreg=0).
    
# THE AUTOREGRESSIVE METHOD IS NOT IMPLEMENTED FOR AUTOREG>1.

    mlemean.ar <- mlemean.uc
    mlevar.ar <- mlevar.uc

    flag <- 0
    iter <- 0
    
    if (autoreg==1)
      {# Start of my EM algorithm for autoregressive model (i.e. autoreg=1)
        while (!flag & iter < maxiter.em)
          {iter <- iter+1
           currentval <- makeparam.norm(prelim, thetalist=list(mlemean.ar, mlevar.ar))
           mlemean.ar0 <- mlemean.ar
           mlevar.ar0 <- mlevar.ar
           mleobj.ar <- em.norm(prelim, maxits=1, start=currentval, showits=F)
           mlemean.ar <- getparam.norm(prelim, mleobj.ar)$mu
           mlevar.ar <- getparam.norm(prelim, mleobj.ar)$sigma
           beta.ar <- matrix(0, 2+dimx, M)
           
           if (dimx==0)
             {for (j in 2:M)
                beta.ar[2,j] <- mlevar.ar[j, j-1] / mlevar.ar[j-1, j-1]
              beta.ar[1, 1] <- mlemean.ar[1]
              beta.ar[1, 2:M] <- mlemean.ar[2:M] - beta.ar[2, 2:M] * mlemean.ar[1:(M-1)]
              for (j in 1:(M-1))
                for (k in (j+1):M)
                  mlevar.ar[j, k] <- prod(beta.ar[2, (j+1):k]) * mlevar.ar[j,j]
              for (j in 1:(M-1))
                for (k in (j+1):M)
                  mlevar.ar[k,j] <- mlevar.ar[j,k]
            } else
           {# Now to do these calculations  when there are covariates
             covx <- mlevar.ar[M+1:dimx, M+1:dimx]
             invcovx <- solve(covx)
             beta.ar[2+1:dimx, 1] <- invcovx %*% mlevar.ar[M+1:dimx, 1]
             beta.ar[1, 1] <- mlemean.ar[1] - beta.ar[2+1:dimx, 1] %*% mlemean.ar[M+1:dimx]
             for (j in 2:M)
               {beta.ar[2, j] <- ( mlevar.ar[j-1, j] - mlevar.ar[j-1, M+1:dimx] %*% invcovx %*% mlevar.ar[M+1:dimx, j] ) / ( mlevar.ar[j-1, j-1] - mlevar.ar[j-1, M+1:dimx] %*% invcovx %*% mlevar.ar[M+1:dimx, j-1] )
                beta.ar[2+1:dimx, j] <- invcovx %*% ( mlevar.ar[M+1:dimx, j] - mlevar.ar[M+1:dimx, j-1] * beta.ar[2,j] )
              }
             for (j in 2:M)
               beta.ar[1, j] <- mlemean.ar[j] - beta.ar[2, j] * mlemean.ar[j-1] - mlemean.ar[M+1:dimx] %*% beta.ar[2+1:dimx, j]
             
             for (j in 1:(M-1))
               for (k in (j+1):M)
                 mlevar.ar[j, k] <- mlevar.ar[k, M+1:dimx] %*% invcovx %*% mlevar.ar[M+1:dimx, j] + prod(beta.ar[2, (j+1):k]) * ( mlevar.ar[j,j] - mlevar.ar[j, M+1:dimx] %*% invcovx %*% mlevar.ar[M+1:dimx, j] )
             
# So far, I have only calculated the upper triangular part of the top-left
# M-by-M submatrix of mlevar.ar.  Now to do the lower triangular part,
# remembering that this is a symmetric matrix.
             for (j in 1:(M-1))
               for (k in (j+1):M)
                 mlevar.ar[k,j] <- mlevar.ar[j,k]
           } # End of calculation when there are covariates
          
           sumdiff <- sum(abs(mlemean.ar - mlemean.ar0)) + sum(abs(mlevar.ar - mlevar.ar0))
           if (sumdiff < toler.emnr)
             flag <- 1
          
# But beta[2,] are the slopes for regression of Y_t on Y_{t-1}, rather than
# (Y_t - Y_{t-1}) on Y_{t-1}, so make it for (Y_t - Y_{t-1}) on Y_{t-1} by
# subtracting 1
           beta.ar[2, 2:M] <- beta.ar[2, 2:M] - 1
# End of EM algorithm for autoregressive model (i.e. autoreg=T)
         }
      } # end of if (autoreg==1) condition
    
    if (autoreg==0)
      {# Use Newton-Raphson algorithm for non-autoregressive case
        temp <- mvn.nonauto(y, x, verbose=verbose)
        iter <- temp$iter.nr
        theta.nonauto <- temp$theta
      
# Extract estimates from theta.non.auto
        beta.ar <- matrix(0, 2+dimx, M)
# Note that the second column of beta.ar is not used.
        beta.ar[1,] <- theta.nonauto[1:M]
        if (dimx>0)
          for (l in 1:dimx)
            beta.ar[2+l,] <- theta.nonauto[M*l+1:M]
        sigma.ar <- theta.nonauto[M*(dimx+1)+1:M]
      
# Now calculate mlemean.ar
        if (dimx>0)
          {mlemean.ar[M+1:dimx] <- apply(x, 2, mean)
           mlemean.ar[1] <- beta.ar[1, 1] + beta.ar[2+1:dimx, 1] %*% mlemean.ar[M+1:dimx]
           for (j in 2:M)
             mlemean.ar[j] <- mlemean.ar[j-1] + beta.ar[1, j] + beta.ar[2+1:dimx, j] %*% mlemean.ar[M+1:dimx]
         } else
        {mlemean.ar[1] <- beta.ar[1, 1]
         for (j in 2:M)
           mlemean.ar[j] <- mlemean.ar[j-1] + beta.ar[1, j]
       }

# Now calculate mlevar.ar
        if (dimx>0)
          {# First get variance-covariance matrix of covariates x
            mlevar.ar[M+1:dimx, M+1:dimx] <- var(x) * (N-1) / N
# multiplying by (N-1)/N gives the mle, rather than the unbiased estimate that the var() function gives (because it divides by N-1 rather than N).
# Second, get covariance between Y1 and covariates x.
            mlevar.ar[1, M+1:dimx] <- beta.ar[2+1:dimx, 1] %*% mlevar.ar[M+1:dimx, M+1:dimx]
            mlevar.ar[M+1:dimx, 1] <- mlevar.ar[1, M+1:dimx]
# Third, get covariate between Yj (j>=2) and covariates x.
            for (j in 2:M)
             mlevar.ar[j, M+1:dimx] <- mlevar.ar[j-1, M+1:dimx] + beta.ar[2+1:dimx, j] %*% mlevar.ar[M+1:dimx, M+1:dimx]
           mlevar.ar[M+1:dimx, 1:M] <- mlevar.ar[1:M, M+1:dimx]
# Fourth, get variance of Y1.
            mlevar.ar[1, 1] <- t(beta.ar[2+1:dimx, 1]) %*% mlevar.ar[M+1:dimx, M+1:dimx] %*% beta.ar[2+1:dimx, 1] + sigma.ar[1]
# Fifth, get variance of Yj (j>=2).
            for (j in 2:M)
              mlevar.ar[j, j] <- mlevar.ar[j-1, j-1] + t(beta.ar[2+1:dimx, j]) %*% mlevar.ar[M+1:dimx, M+1:dimx] %*% beta.ar[2+1:dimx, j] + 2 * t(beta.ar[2+1:dimx, j]) %*% mlevar.ar[j-1, M+1:dimx] + sigma.ar[j]
          } else
        {# All that is simpler if there are no covariates.
          mlevar.ar[1, 1] <- sigma.ar[1]
          for (j in 2:M)
            mlevar.ar[j, j] <- mlevar.ar[j-1, j-1] + sigma.ar[j]
        }

# Sixth, get covariance between Yj and Yk.
        for (j in 1:(M-1))
          for (k in (j+1):M)
            {if (dimx>0)
               mlevar.ar[j, k] <- mlevar.ar[j, j] + t(mlevar.ar[j, M+1:dimx]) %*% solve(mlevar.ar[M+1:dimx, M+1:dimx]) %*% mlevar.ar[k, M+1:dimx] - t(mlevar.ar[j, M+1:dimx]) %*% solve(mlevar.ar[M+1:dimx, M+1:dimx]) %*% mlevar.ar[M+1:dimx, j] else
             mlevar.ar[j, k] <- mlevar.ar[j, j]
           }
        for (j in 1:(M-1))
          for (k in (j+1):M)
            mlevar.ar[k, j] <- mlevar.ar[j, k]
      } # End of Newton-Raphson algorithm for non-autoregressive case

    if (autoreg>1)
      {# Just set up an empty beta vector for the autoregressive/unconstrained method
        beta.ar <- matrix(0, 1+autoreg1+dimx, M)
      }

    betatab <- rbind(beta.ag, beta.uc, beta.ar)
    if (dimx==0)
      dimnames(betatab) <- list(c("alpha.ag", rep("beta.ag", autoreg1), "alpha.uc", rep("beta.uc", autoreg1), "alpha.ar", rep("beta.ar", autoreg1)), paste("t", 1:M, sep="")) else
    dimnames(betatab) <- list(c("alpha.ag", rep("beta.ag", autoreg1), paste("gamma", 1:dimx, ".ag", sep=""), "alpha.uc", rep("beta.uc", autoreg1), paste("gamma", 1:dimx, ".uc", sep=""), "alpha.ar", rep("beta.ar", autoreg1), paste("gamma", 1:dimx, ".ar", sep="")), paste("t", 1:M, sep=""))

# Estimating the compensator

    meanDelta <- matrix(0, 3, M)
    dimnames(meanDelta) <- list(c("ag", "uc", "ar"), paste("t", 1:M, sep=""))
    Delta <- matrix(0, N, M)

    for (method in 1:3)
      {beta <- betatab[1:(1+autoreg1+dimx) + (method-1)*(1+autoreg1+dimx),]
       
       Delta[,1] <- rep(beta[1,1], N)
          if (dimx>0)
            for (k in 1:dimx)
              Delta[,1] <- Delta[,1] + x[,k] * beta[1+autoreg1+k, 1]

       for (j in 2:M)
         {furthest <- max(j-autoreg, 1)
          firstindex <- max(1, autoreg-j+2)
           Delta[,j] <- beta[1, j] + Delta[, j-1]
          if (autoreg>0)
            for (k in 1:autoreg)
              if (k<=j-1)
                Delta[,j] <- Delta[,j] + beta[autoreg+2-k, j] * Delta[, j-k]
          if (dimx>0)
            for (k in 1:dimx)
              Delta[,j] <- Delta[,j] + x[,k] * beta[1+autoreg1+k, j]
        }
       meanDelta[method,] <- apply(Delta, 2, mean)
     }
   
# Perform the imputation in the Aalen and Gunnes way

    meany <- matrix(0, 3, M)
    meany.alive <- matrix(0, 3, M)
    meany[1,] <- apply(y, 2, mean)
    meany.alive[1,] <- apply(y*alive, 2, mean) / apply(alive, 2, mean)

    if (returnimps)
      {yimp.ag <- array(0, c(N, M, 3))
       dimnames(yimp.ag) <- list(1:N, paste("t", 1:M, sep=""), c("ag", "uc", "ar"))
     }

    for (method in 1:3)
      {beta <- betatab[1:(1+autoreg1+dimx) + (method-1)*(1+autoreg1+dimx),]
       yest <- yx[, 1:M]
       missing <- is.na(yest[,1])
       yest[missing==T, 1] <- rep(beta[1,1], sum(missing))
       if (dimx>0)
         for (k in 1:dimx)
           yest[missing==T, 1] <- yest[missing==T, 1] + x[missing==T, k] * beta[1+autoreg1+k, 1]

       for (j in 2:M)
         {missing <- is.na(yest[,j])
          furthest <- max(j-autoreg, 1)
          firstindex <- max(1, autoreg-j+2)
          yest[missing==T, j] <- beta[1, j] + yest[missing==T, j-1]
          if (autoreg>0)
            for (k in 1:autoreg)
              if (k<=j-1)
                yest[missing==T, j] <- yest[missing==T, j] + yest[missing==T, j-k] * beta[autoreg1+2-k, j]
          if (dimx>0)
            yest[missing==T, j] <- yest[missing==T, j] + as.matrix(x[missing==T,]) %*% beta[1+autoreg1+1:dimx, j]
        }
       meany[method,] <- apply(yest, 2, mean)
       meany.alive[method,] <- apply(yest*alive, 2, mean) / apply(alive, 2, mean)
       if (returnimps)
         yimp.ag[,, method] <- yest
     }
    
# Perform mean imputation using MVN with MAR    
    yimp.uc <- expect.cond.norm(yx, mlemean.uc, mlevar.uc)[, 1:M]
    mean.mvn.uc <- apply(yimp.uc, 2, mean)
    yimp.ar <- expect.cond.norm(yx, mlemean.ar, mlevar.ar)[, 1:M]
    mean.mvn.ar <- apply(yimp.ar, 2, mean)
    mean.mvn.alive.uc <- apply(yimp.uc * alive, 2, mean) / apply(alive, 2, mean)
    mean.mvn.alive.ar <- apply(yimp.ar * alive, 2, mean) / apply(alive, 2, mean)

    if (printsigmas)
      {print("MLE of variance matrix for unstructured model is")
       print(mlevar.uc)
       print("MLE of variance matrix for autoregressive model is")
       print(mlevar.ar)
     }
   
# Now return the results
    mlemeanmvn <- rbind(mlemean.uc, mlemean.ar)
    if (dimx==0)
      dimnames(mlemeanmvn) <- list(c("uc", "ar"), paste("t", 1:M, sep="")) else
    dimnames(mlemeanmvn) <- list(c("uc", "ar"), c(paste("t", 1:M, sep=""), paste("x", 1:dimx, sep="")))
    meanimpy <- rbind(meany, mean.mvn.uc, mean.mvn.ar)
    meanimpy.alive <- rbind(meany.alive, mean.mvn.alive.uc, mean.mvn.alive.ar)
    dimnames(meanimpy) <- list(c("ag", "uc", "ar", "mvn.uc", "mvn.ar"), paste("t", 1:M, sep=""))
    dimnames(meanimpy.alive) <- list(c("ag", "uc", "ar", "mvn.uc", "mvn.ar"), paste("t", 1:M, sep=""))

    if (!returnimps)
      return( list(betatab=betatab, meanDelta=meanDelta, mlemeanmvn=mlemeanmvn, meanimpy=meanimpy, meanimpy.alive=meanimpy.alive, numiters.emnr=iter) ) else
    return( list(betatab=betatab, meanDelta=meanDelta, mlemeanmvn=mlemeanmvn, meanimpy=meanimpy, meanimpy.alive=meanimpy.alive, numiters.emnr=iter, yimp.ag=yimp.ag, yimp.uc=yimp.uc, yimp.ar=yimp.ar) )
  } # end of linearincrements.doall() function



expect.cond.norm <- function(yx, mu, sigma)
  {# This function is like Schafer's imp.norm, but is deterministic, rather
# than stochastic.  It replaces missing values of YX by their conditional
# expectations given the observed values of YX and assuming that
# YX ~ MVN(mu, sigma).  In practice, YX is the matrix in which row i is the vector
# (Y_i1, ..., Y_iT, X_i).

    M <- ncol(yx)
    N <- nrow(yx)
    n.pattern <- 2^M
    expected.y <- yx
    r <- !is.na(yx)

    for (pat in 1:n.pattern)
      {pattern <- decimal2binary(pat-1, length=M)
       use <- rep(1, N)
       for (j in 1:M)
         use[r[,j]!=pattern[j]] <- 0
       if (sum(use)>0)
         {mu.mis <- mu[pattern==0]
          mu.obs <- mu[pattern==1]
          sigma.obs.obs <- sigma[pattern==1, pattern==1]
          sigma.mis.obs <- sigma[pattern==0, pattern==1]
          sigmaprod <- sigma.mis.obs %*% solve(sigma.obs.obs)
          deviations <- yx[use==T, pattern==1] - matrix(mu.obs, nrow=sum(use), ncol=sum(pattern==1), byrow=T)
         correction <- t( sigmaprod %*% t(deviations) )
          expected.y[use==T, pattern==0] <- matrix(mu.mis, nrow=sum(use), ncol=sum(pattern==0), byrow=T) + correction
        }
     }
    return(expected.y)
  }




mvn.nonauto <- function(y, x=NULL, max.it=10000, toler.nr=1e-5, verbose=F)
  {# This function fits the MVN model with the regression coefficient of y[t-1]
#constrained to equal one.  It uses the Newton-Raphson algoritm.

# To get initial values for the Newton-Raphson algorithm, use the Aalen and
# Gunnes method for estimating parameters with the regr coeff of y[t-1]
# constrained to equal 1.

    N <- nrow(y)
    M <- ncol(y)

    if (is.null(x))
      dimx <- 0 else
    {if (is.vector(x))
       x <- as.matrix(x, N, 1)
     dimx <- ncol(x)
   }

    r <- !is.na(y)
    
    alphagamma.initial <- matrix(0, 1+dimx, M)
    sigma2.initial <- rep(0, M)

    n.used <- sum(r[,1]==1)
    if (dimx>0)
      {temp <- lm(y[,1] ~ x, subset=r[,1]==1)
       alphagamma.initial[, 1] <- temp$coefficients
       sigma2.initial[1] <- summary(temp)$sigma^2 * (n.used-1-dimx) / n.used
     } else
      {temp <- lm(y[,1] ~ 1, subset=r[,1]==1)
       alphagamma.initial[, 1] <- mean(y[r[,1]==1, 1])
       sigma2.initial[1] <- var(y[r[,1]==1, 1]) * (n.used-1) / n.used
     }

    if (dimx>0)
      for (j in 2:M)
        {temp <- lm(y[,j] - y[,j-1] ~ x, subset=r[,j]==1 & r[,j-1]==1)
         alphagamma.initial[, j] <- temp$coefficients
         n.used <- sum(r[,j]==1 & r[,j-1]==1)
         sigma2.initial[j] <- summary(temp)$sigma^2 * (n.used-1-dimx) / n.used
       } else
    for (j in 2:M)
      {temp <- lm(y[,j] - y[,j-1] ~ 1, subset=r[,j]==1 & r[,j-1]==1)
       alphagamma.initial[1, j] <- temp$coefficients
       n.used <- sum(r[,j]==1 & r[,j-1]==1)
       sigma2.initial[j] <- summary(temp)$sigma^2 * (n.used-1) / n.used
     }
         
# Put all parameters into the vectors alphas, gammas and sigma2s
    alpha <- alphagamma.initial[1,]
    if (dimx>0)
      gamma <- matrix(alphagamma.initial[(1:dimx)+1,], nrow=dimx, ncol=M)
    sigma2 <- sigma2.initial

# Put all parameters in the vector theta (first alphas, then gammas and
# finally sigma2s)
    theta <- rep(0, (2+dimx)*M)
    theta[1:M] <- alpha
    if (dimx>0)
      for (j in 1:dimx)
        theta[1:M + j*M] <- gamma[j,]
    theta[1:M + (dimx+1)*M] <- sigma2

# For each person and each timepoint, find the last observation time strictly
# before that timepoint (a) and the first observation time at or after that
# timepoint (b).

    a <- matrix(0, N, M)
    b <- matrix(0, N, M)

    for (k in M:1)
      b[r[,k]==1, 1] <- k

    for (j in 2:M)
      {for (k in 1:(j-1))
         a[r[,k]==1,j] <- k
       for (k in M:j)
         b[r[,k]==1,j] <- k
     }

    flag <- T
    it <- 1
    max.it <- 100
    
    loglike <- -1e99

# Now begin the iterative algorithm.  flag==T until converge or reached
# maximum iteration number.
    
    while (flag)
      {# Extracting alpha, gamma and sigma2 from theta
        alpha <- theta[1:M]
        if (dimx>0)
          {gamma <- matrix(0, dimx, M)
           for (j in 1:dimx)
             gamma[j,] <- theta[M*j+1:M]
         }
        sigma2 <- theta[M*(dimx+1)+1:M]

        if (dimx==1)
          x <- matrix(x, nrow=N, ncol=1)

# Calculate the log likelihood
        loglike.old <- loglike
        loglike <- 0
        for (j in 1:M)
          for (k in (j-1):0)
            {use <- r[,j]==1 & a[,j]==k
             sigma2sum <- sum(sigma2[(k+1):j])
             loglike <- loglike - sum(use) * 0.5 * log(sigma2sum)
             if (k>0)
               ydiff <- y[,j] - y[,k]
             else
               ydiff <- y[,j]
             
             if (dimx>1 & j>k+1)
               xgamma <- x %*% apply(gamma[,(k+1):j], 1, sum)
             if (dimx>1 & j==k+1)
               xgamma <- x %*% gamma[,j]
             if (dimx==1)
               xgamma <- x * sum(gamma[1, (k+1):j])
             if (dimx==0)
               xgamma <- 0
             loglike <- loglike - 0.5 * sigma2sum^(-1) * sum( ( (ydiff - sum(alpha[(k+1):j]) - xgamma)^2)[use==T] )
           }
        
        dll.dsigma2 <- rep(0, M)
        dll.dalpha <- rep(0, M)
        if (dimx>0)
          dll.dgamma <- matrix(0, dimx, M)
    
        for (j in 1:M)
          for (kb in j:M)
            for (ka in 0:(j-1))
              {use <- a[,j]==ka & b[,j]==kb
               sigma2sum <- sum(sigma2[(ka+1):kb])
               if (ka>0)
                 ydiff <- y[,kb] - y[,ka]
               else
                 ydiff <- y[,kb]

               if (dimx>1 & kb>ka+1)
                 xgamma <- x %*% apply(gamma[,(ka+1):kb], 1, sum)
               if (dimx>1 & kb==ka+1)
                 xgamma <- x %*% gamma[, kb]
               if (dimx==1)
                 xgamma <- x * sum(gamma[1, (ka+1):kb])
               if (dimx==0)
                 xgamma <- 0
               
               alphasum <- sum(alpha[(ka+1):kb])
               
               dll.dsigma2[j] <- dll.dsigma2[j] - sum(use) * 0.5 * sigma2sum^(-1)
               dll.dsigma2[j] <- dll.dsigma2[j] + 0.5 * sigma2sum^(-2) * sum( ( (ydiff - alphasum - xgamma)^2 )[use==T] )

               dll.dalpha[j] <- dll.dalpha[j] + sigma2sum^(-1) * sum( (ydiff - alphasum - xgamma)[use==T] )
               
               if (dimx>0)
                 for (l in 1:dimx)
                   dll.dgamma[l, j] <- dll.dgamma[l, j] + sigma2sum^(-1) * sum( ( (ydiff - alphasum - xgamma) * x[,l] )[use==T] )
      }
        
        d2ll.dsigma2.dsigma2 <- matrix(0, M, M)
        d2ll.dalpha.dalpha <- matrix(0, M, M)
        d2ll.dsigma2.dalpha <- matrix(0, M, M)
        if (dimx>0)
          {d2ll.dgamma.dgamma <- array(0, c(M, M, dimx, dimx))
           d2ll.dsigma2.dgamma <- array(0, c(M, M, dimx))
           d2ll.dalpha.dgamma <- array(0, c(M, M, dimx))
         }
 
        for (j1 in 1:M)
          for (j2 in j1:M)
            for (ka in 0:(j1-1))
              for (kb in j2:M)
                {use <- a[,j1]==ka & a[,j2]==ka & b[,j1]==kb & b[,j2]==kb
                 sigma2sum <- sum(sigma2[(ka+1):kb])
                 if (ka>0)
                   ydiff <- y[,kb] - y[,ka]
                 else
                   ydiff <- y[,kb]

                 if (dimx>1 & kb>ka+1)
                   xgamma <- x %*% apply(gamma[,(ka+1):kb], 1, sum)
                 if (dimx>1 & kb==ka+1)
                   xgamma <- x %*% gamma[, kb]
                 if (dimx==1)
                   xgamma <- x * sum(gamma[1, (ka+1):kb])
                 if (dimx==0)
                   xgamma <- 0
                 
                 alphasum <- sum(alpha[(ka+1):kb])
                 
                 d2ll.dsigma2.dsigma2[j1, j2] <- d2ll.dsigma2.dsigma2[j1, j2] + sum(use) * 0.5 * sigma2sum^(-2)
                 d2ll.dsigma2.dsigma2[j1, j2] <- d2ll.dsigma2.dsigma2[j1, j2] - sigma2sum^(-3) * sum( ( (ydiff - alphasum - xgamma)^2 )[use==T] )

                 d2ll.dalpha.dalpha[j1, j2] <- d2ll.dalpha.dalpha[j1, j2] - sum(use) * sigma2sum^(-1)
                 
                 d2ll.dsigma2.dalpha[j1, j2] <- d2ll.dsigma2.dalpha[j1, j2] - sigma2sum^(-2) * sum( (ydiff - alphasum - xgamma)[use==T] )
                 
                 if (dimx>0)
                   {for (l1 in 1:dimx)
                      for (l2 in 1:dimx)
                        d2ll.dgamma.dgamma[j1, j2, l1, l2] <- d2ll.dgamma.dgamma[j1, j2, l1, l2] - sum(x[use, l1] * x[use, l2]) * sigma2sum^(-1)
                    for (l in 1:dimx)
                      {d2ll.dsigma2.dgamma[j1, j2, l] <- d2ll.dsigma2.dgamma[j1, j2, l] - sigma2sum^(-2) * sum( ( (ydiff - alphasum - xgamma) * x[,l] )[use==T] )
                       d2ll.dalpha.dgamma[j1, j2, l] <- d2ll.dalpha.dgamma[j1, j2, l] - sum(x[use, l]) * sigma2sum^(-1)
                     }
                  }
               }

# I have only calculated the diagonal and above diagonal so far, so now use
# fact that these matrices are symmetric
        if (M>2)
          {for (j1 in 1:(M-1))
             for (j2 in (j1+1):M)
               {d2ll.dsigma2.dsigma2[j2, j1] <- d2ll.dsigma2.dsigma2[j1, j2]
                d2ll.dalpha.dalpha[j2, j1] <- d2ll.dalpha.dalpha[j1, j2]
                d2ll.dsigma2.dalpha[j2, j1] <- d2ll.dsigma2.dalpha[j1, j2]
                if (dimx>0)
                  {for (l1 in 1:dimx)
                     for (l2 in 1:dimx)
                       d2ll.dgamma.dgamma[j2, j1, l1, l2] <- d2ll.dgamma.dgamma[j1, j2, l1, l2]
                   for (l in 1:dimx)
                     {d2ll.dsigma2.dgamma[j2, j1, l] <- d2ll.dsigma2.dgamma[j1, j2, l]
                      d2ll.dalpha.dgamma[j2, j1, l] <- d2ll.dalpha.dgamma[j1, j2, l]
                    }
                 }
              }
         } # end of if (M>2) condition
        
        if (dimx>0)
          {dll.dgamma.vec <- rep(0, M*dimx)
           for (l in 1:dimx)
             dll.dgamma.vec[(l-1)*M + 1:M] <- dll.dgamma[l,]
           dll.dtheta <- c(dll.dalpha, dll.dgamma.vec, dll.dsigma2)
         } else
        dll.dtheta <- c(dll.dalpha, dll.dsigma2)
        
        d2ll.dtheta.dtheta <- matrix(0, M*(2+dimx), M*(2+dimx))
        
        d2ll.dtheta.dtheta[1:M, 1:M] <- d2ll.dalpha.dalpha
        d2ll.dtheta.dtheta[M*(1+dimx)+1:M, M*(1+dimx)+1:M] <- d2ll.dsigma2.dsigma2
        d2ll.dtheta.dtheta[1:M, M*(1+dimx)+1:M] <- d2ll.dsigma2.dalpha

        if (dimx>0)
          {for (l1 in 1:dimx)
             for (l2 in 1:dimx)
               d2ll.dtheta.dtheta[M*l1+1:M, M*l2+1:M] <- d2ll.dgamma.dgamma[,, l1, l2]
           for (l in 1:dimx)
             {d2ll.dtheta.dtheta[1:M, M*l+1:M] <- d2ll.dalpha.dgamma[,,l]
              d2ll.dtheta.dtheta[M*l+1:M, M*(1+dimx)+1:M] <- d2ll.dsigma2.dgamma[,,l]
            }
         }

# I have only calculated the diagonal and above diagonal so far, so now use
# fact that this matrices is symmetric
        for (j1 in 1:(M*(2+dimx)-1))
          for (j2 in (j1+1):(M*(2+dimx)))
            d2ll.dtheta.dtheta[j2, j1] <- d2ll.dtheta.dtheta[j1, j2]

        theta.old <- theta
        theta <- theta - as.vector( dll.dtheta %*% solve(d2ll.dtheta.dtheta) )

        if (verbose==T)
          {print(paste("iteration", it))
           print("theta at beginning of iteration")
           print(theta.old)
           print(paste("sum of absolute changes in theta =", sum(abs(theta - theta.old))))
           print(paste("change in loglike =", loglike - loglike.old))
         }
        
        it <- it+1
        if (it>max.it)
          flag <- F
        if (sum(abs(theta - theta.old))<toler.nr & loglike.old - loglike < toler.nr)
          flag <- F
      }

    return(list(theta=theta, iter.nr=it))
  }
