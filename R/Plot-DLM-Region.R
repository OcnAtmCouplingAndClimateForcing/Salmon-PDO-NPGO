  species <- 'Chinook'
  
  plot_DLM_Region <- function(dir.figs=file.path(getwd(),"Figs","DLM-Region"), 
                              dir.output=file.path(getwd(),"Output","DLM-Region"),
                              species="Sockeye") {
    
    ###TESTING###
    # species <- "Coho" #"Sockeye" "Chinook" "Pink"    "Chum"    "Coho"   
    # dir.figs <- file.path(getwd(),"Figs","DLM-Region")
    # dir.output <- file.path(getwd(),"Output","DLM-Region")
    
    #Load ======================
    require(rstan)
    require(tidyverse)
    require(dplyr)
    require(bayesplot)
    require(rstanarm)
    require(ggplot2)
    require(ggthemes)
    require(reshape2)
    
    #Load the STAN Fit & Extras ============================
    fit <- readRDS(file=file.path(dir.output,paste0(species,'-fit.rds')))
    extras <- readRDS(file=file.path(dir.output,paste0(species,'-extras.rds')))             
    
    #Extract parameters
    pars <- rstan::extract(fit)
    
    #Plot Convergence ==================================
    summary <- data.frame(summary(fit)$summary)
    conv.diag <- reshape2::melt(summary[,c(9:10)])
    g <- ggplot(conv.diag, aes(value, fill=variable)) +
      theme_linedraw() +
      scale_fill_colorblind() +
      geom_histogram(color='red') +
      facet_wrap(~variable, scales='free') + 
      ggtitle(species, subtitle='Convergence Diagnostics') +
      theme(legend.position = 'none')
    g
    ggsave(file.path(dir.figs,paste0(species,"-Diagnostics.png")), plot=g, height=4, width=4, units='in',
           dpi=500)
    #Caterpillar plot for alpha and beta coefficients =========================
    # posterior <- as.array(fit)
    # 
    # g <- bayesplot::mcmc_intervals(post$alpha)
    # g
    
    #Plot PDO/NPGO Trajectories ======================
    coefs.pdo <- apply(pars$coef_PDO, c(2,3), quantile, probs=c(0.025,0.25,0.5,0.75,0.975))
    coefs.npgo <- apply(pars$coef_NPGO, c(2,3), quantile, probs=c(0.025,0.25,0.5,0.75,0.975))
  
    dimnames(coefs.pdo)[[2]] <- extras$stock.regions
    dimnames(coefs.pdo)[[3]] <- extras$stock.years
  
    dimnames(coefs.npgo)[[2]] <- extras$stock.regions
    dimnames(coefs.npgo)[[3]] <- extras$stock.years
  
    list.coefs.pdo <- melt(coefs.pdo)
    list.coefs.npgo <- melt(coefs.npgo)
    names(list.coefs.pdo) <- c('quant','large.region','broodYr','value')
    names(list.coefs.npgo) <- c('quant','large.region','broodYr','value')
    
    #Combine 
    list.coefs.npgo$coef <- 'NPGO'
    list.coefs.pdo$coef <- 'PDO'
    list.coefs <- rbind(list.coefs.pdo, list.coefs.npgo)
    list.coefs$large.region <- factor(list.coefs$large.region, levels=c('BS','GOA','CA','WC'))
    
  
    g <- ggplot(filter(list.coefs, quant=='50%'), aes(x=broodYr, y=value, color=factor(large.region))) +
           theme_linedraw() +
           geom_hline(yintercept=0) +
           geom_vline(xintercept=1989, lty=2) +
           geom_line() +
           ggtitle(paste('Random Walk Coefficients'), subtitle=paste(species,"salmon")) +
           facet_wrap(~coef, ncol=1, scales='free_y') +
           xlab('Interaction Year')
    g
    
    ggsave(file.path(dir.figs,paste0(species,"-Coefficient Trends.png")), plot=g, height=6, width=8, units='in',
             dpi=500)
    
    #Wide format
    list.coefs.wide <- list.coefs %>% spread(key=quant, value=value)
    names(list.coefs.wide)[4:8] <- c('low.95','low.50','med','up.50','up.95')
    list.coefs.wide$large.region <- factor(list.coefs.wide$large.region, levels=c('BS','GOA','CA','WC'))
    
    g2 <- ggplot(list.coefs.wide, aes(x=broodYr, y=med, color=factor(large.region), fill=factor(large.region))) +
            theme_linedraw() +
            geom_hline(yintercept=0) +
            geom_vline(xintercept=1989, lty=2) +
            geom_line() +
            geom_ribbon(aes(ymin=low.50, ymax=up.50), alpha=0.25, lwd=1e-6) +
            ggtitle(paste('Random Walk Coefficients'), subtitle=paste(species,"salmon")) +
            facet_wrap(~coef, ncol=1, scales='free_y') +
            xlab('Interaction Year')
    g2
    ggsave(file.path(dir.figs,paste0(species,"-Coefficient Trends_50CI.png")), plot=g2, height=6, width=8, units='in',
           dpi=500)
    
    g3 <- ggplot(list.coefs.wide, aes(x=broodYr, y=med, color=factor(large.region), fill=factor(large.region))) +
      theme_linedraw() +
      geom_hline(yintercept=0) +
      geom_vline(xintercept=1989, lty=2) +
      geom_ribbon(aes(ymin=low.95, ymax=up.95), alpha=0.25, lwd=1e-6) +
      geom_ribbon(aes(ymin=low.50, ymax=up.50), alpha=0.25, lwd=1e-6) +
      geom_line() +
      ggtitle(paste('Random Walk Coefficients'), subtitle=paste(species,"salmon")) +
      facet_grid(coef~large.region, scales='free_y') +
      xlab('Interaction Year')
    g3
    ggsave(file.path(dir.figs,paste0(species,"-Coefficient Trends_pannel.png")), plot=g3, height=8, width=10, units='in',
           dpi=500)
    
    # Plot histogram of phi if present
    if("phi" %in% names(pars)) {
      pdf(file.path(dir.figs,paste0(species,"-Autocorr Plots.pdf")), height=8, width=8)
      # Plot Phi
      hist(pars$phi, col='lightblue')
      abline(v=0, col='red', lty=3)
      
      # str(pars$init_residual)
      # hist(pars$init_residual[,10])
      dev.off()
    }
    
  }
  
  # Single species
    # plot_DLM_Region(dir.figs=file.path(getwd(),"Figs","DLM-Region"),
    #                             dir.output=file.path(getwd(),"Output","DLM-Region"),
    #                             species="Chinook")
  
  # Plot ALL species ========================================
  for(species in c("Sockeye","Pink", "Chum", "Coho", "Chinook")) {
    print(species)
    plot_DLM_Region(dir.figs=file.path(getwd(),"Figs","DLM-Region"),
                                dir.output=file.path(getwd(),"Output","DLM-Region"),
                                species=species)
  }
