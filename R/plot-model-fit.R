#Plots Model Model Fit
plot_model_fit <- function(fit, dat.input) {
  require(rstan)
  #Extract pars
  pars <- rstan::extract(fit)
  #Extract Location Info
  temp.species <- unique(dat.input$species)
  temp.loc <- unique(dat.input$stock)
  #Calculate quantiles
  pred.mat <- apply(pars$pred, 2, quantile, probs=c(0.025,0.25,0.5,0.75,0.975))
  y.lim <- c(min(dat.input$ln.rps,pred.mat), max(dat.input$ln.rps, pred.mat))
  #Observed
  par(mfrow=c(1,1), oma=rep(0,4), mar=c(4,4,2,1))
  plot(ln.rps~broodYr, data=dat.input, type='l', ylim=y.lim, col='blue',
       main=paste0(temp.species,": ",temp.loc),
       xlab='Brood Year', ylab='ln(Recruits per Spawner)')
  points(ln.rps~broodYr, data=dat.input, pch=21, bg='blue')
  #Predicted
  polygon(x=c(dat.input$broodYr,rev(dat.input$broodYr)), 
          y=c(pred.mat[1,],rev(pred.mat[5,])), 
          col=rgb(1,0,0, alpha=0.25), border=FALSE)
  polygon(x=c(dat.input$broodYr,rev(dat.input$broodYr)), 
          y=c(pred.mat[2,],rev(pred.mat[4,])), 
          col=rgb(1,0,0, alpha=0.25), border=FALSE)
  lines(x=dat.input$broodYr, y=pred.mat[3,], lwd=2, col='red')
  #Return Section
  out <- NULL
  out$pred.mat <- pred.mat
  out$broodYr <- dat.input$broodYr
  out$ln.rps <- dat.input$ln.rps
  return(out)
}
