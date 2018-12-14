#Plot Time-varying Coefficients

plot_coef_ts <- function(fit, dat.input, temp.offset=2) {
  require(rstan)
  #Extract pars
  pars <- rstan::extract(fit)
  #Extract Location Info
  temp.species <- unique(dat.input$species)
  temp.loc <- unique(dat.input$stock)
  #Calculate quantiles
  coef.mat <- apply(pars$coef, c(2,3), quantile, probs=c(0.025,0.25,0.5,0.75,0.975))
  #Determine Years
  broodYr <- dat.input$broodYr
  years <- dat.input$broodYr + temp.offset
  
  #Observed
  y.lim <- c(min(coef.mat), max(coef.mat))
  par(mfrow=c(1,1), oma=rep(0,4), mar=c(4,4,2,1))
  
  plot(x=NULL, y=NULL, xlim=c(min(years),max(years)), ylim=y.lim,
         main=paste0(temp.species,": ",temp.loc),
         xlab='Calendar Year', ylab='Coefficient Value')
  abline(h=0, lty=2)
  #PDO - Blue
  polygon(x=c(years,rev(years)), 
          y=c(coef.mat[1,,1],rev(coef.mat[5,,1])), 
          col=rgb(0,0,1, alpha=0.25), border=FALSE)
  polygon(x=c(years,rev(years)), 
          y=c(coef.mat[2,,1],rev(coef.mat[4,,1])), 
          col=rgb(0,0,1, alpha=0.25), border=FALSE)
  lines(x=years, y=coef.mat[3,,1], lwd=2, col='blue')
  
  #NPGO - Red
  polygon(x=c(years,rev(years)), 
          y=c(coef.mat[1,,2],rev(coef.mat[5,,2])), 
          col=rgb(1,0,0, alpha=0.25), border=FALSE)
  polygon(x=c(years,rev(years)), 
          y=c(coef.mat[2,,2],rev(coef.mat[4,,2])), 
          col=rgb(1,0,0, alpha=0.25), border=FALSE)
  lines(x=years, y=coef.mat[3,,2], lwd=2, col='red')
  

  
  #Legend
  legend('bottomleft', legend=c('PDO','NPGO'), text.col=c('blue','red'), 
         bty='n')
  
  #Return Section
  out <- NULL
  out$coef.mat <- coef.mat
  out$broodYr <- broodYr
  out$years <- years
  return(out)
}