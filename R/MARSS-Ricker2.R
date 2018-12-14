#MARSS DLM Ricker Model - Fixed Intercept + Time Varying PDO/NPGO effects
MARSS_Ricker2 <- function(dat.input, covars) {
  require(MARSS)

##########################
#RUN DLM MODEL

years <- dat.input$broodYr
TT <- length(years)

## get response data: logit(survival)
dat <- matrix(dat.input$ln.rps, nrow = 1)

## get regressor variable
pred <- t(covars) #PDO then NPGO

## number of regr params (slope + intercept)
m = dim(pred)[1]# + 1

##########################
#PROCESS EQUATION
B = diag(m) # 2x2; Identity
# B[1,1] <- 0
# U = matrix(0,nrow=m,ncol=1) # 2x1; both elements = 0 --- WE ONLY DO THIS IF WE HAVE DE-MEANED THE DATA
Q = matrix(list(0),m,m) # 2x2; all 0 for now
diag(Q) = c("q.pdo","q.npgo") # 2x2; diag = (q1,q2)

#DEFINE OBSERVATION MODEL
# for observation eqn
Z = array(NA, c(1,m,TT)) # NxMxT; empty for now
# Z[1,1,] = rep(1,TT) # Nx1; 1's for intercept
Z[1,1:2,] = pred # Nx1; regr variable
A = matrix(0) # 1x1; scalar = 0
R = matrix("r") # 1x1; scalar = r

#DEFINE STARTING VALUES
# only need starting values for regr parameters
#inits.list = list(x0=matrix(c(mean(resp.data), mean(resp.data/(pred.data+0.00001))), nrow=m))
inits.list <- list(x0=matrix(c(0,0),nrow=m))

# list of model matrices & vectors
mod.list = list(B=B, Q=Q, Z=Z, A=A, R=R, D='unconstrained', d=t(cbind(1,dat.input$spawn)))#d=matrix(1))

# mod.list = list(B=B, Q=Q, Z=Z, A=A, R=R, D='unconstrained', d=matrix(1))


#FIT THE MODEL
# fit univariate DLM
# foo <- MARSS(dat, inits=inits.list, model=mod.list, fit=FALSE)
dlm <- MARSS(dat, inits=inits.list, model=mod.list, control=list(maxit=1e4), silent=TRUE)

# names(dlm)

# #Plot the Fit
# plot(ln.rps~broodYr, data=dat.input, type='b')
# lines(x=dat.input$broodYr, y=dlm$y.se, col='red')
# 
# #Get Parameter states
# dlm$states
# 
# MARSSparamCIs(dlm)
  

  #Return Section
  out <- NULL
  out$dlm <- dlm
  out$dat.input <- dat.input
  # out$covars <- covars
  return(out)

}