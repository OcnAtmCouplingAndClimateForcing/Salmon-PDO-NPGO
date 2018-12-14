#Plot Other Time-varying Ricker Parameters

plot_other_pars <- function(fit) {
  require(BEST)
  require(rstan)
  #Extract Parameters
  pars <- rstan::extract(fit)
  #Plot
  par(mfrow=c(3,2), oma=rep(0,4), mar=c(4,3,1,3))
  BEST::plotPost(pars$alpha, xlab='alpha: log(max RpS)')
  BEST::plotPost(pars$beta, xlab='beta')
  BEST::plotPost(pars$sigma_oe, xlab='sigma oe')
  BEST::plotPost(pars$sigma_pe[,1], xlab='sigma pe: PDO')
  BEST::plotPost(pars$sigma_pe[,2], xlab='sigma pe: NPGO')
}