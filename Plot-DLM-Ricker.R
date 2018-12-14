#==================================================================================================
#Project Name: COAST-WIDE SALMON RECRUITMENT ANALYSIS - Plot DLM-Ricker
#Creator: Curry James Cunningham, NOAA/NMFS, ABL
#Date: 11.27.18
#
#Purpose: To plot and summarize results from a dynamic linear Ricker model with time-varying PDO/NPGO 
#
#
#==================================================================================================
#NOTES:
#
#==================================================================================================

require(tidyverse)
require(dplyr)
require(ggthemes)
require(cowplot)
require(viridis)
require(rstan)
require(rstantools)
require(bayesplot)
require(shinystan)
require(BEST)

#CONTROL SECTION ==========================================================
read.rds <- FALSE #Whether to read in data from stan output files

#Define Workflow Paths ====================================================
# *Assumes you are working from the Coastwide_Salmon_Analysis R project
wd <- getwd()
dir.data <- file.path(wd,"Data")
dir.figs <- file.path(wd,"Figs","DLM-Ricker")
dir.output <- file.path(wd,"Output","DLM-Ricker")
dir.R <- file.path(wd,"R")

#Create
dir.create(dir.figs, recursive=TRUE)
dir.create(dir.output, recursive=TRUE)

#Source Necessary Functions =====================================
source(file.path(dir.R, "plot-model-fit.R"))
source(file.path(dir.R, "plot-coef-ts.R"))
source(file.path(dir.R, "plot-other-pars.R"))

#Read SR Data ===================================================
dat <- read.csv(file.path(dir.data,"full.dat.4.csv"), header=TRUE, stringsAsFactors=FALSE)

#Add rps and ln.rps
dat$rps <- dat$rec/dat$spawn
dat$ln.rps <- log(dat$rps)

#Read PDO/NPGO Data =============================================
covar.dat <- read.csv(file.path(dir.data,"Covariates","covar.dat.csv"))
(covars <- unique(covar.dat$name))

#Extract Metadata ===================================================
species <- unique(dat$species)
n.species <- length(species)

#Specify Offsets for Covariates ================================
#From BroodYear
offset.table <- read.csv(file=file.path(dir.figs,"offset.table.csv"))
offset <- offset.table$offset

#Subset SR Data ==============================================
dat.2 <- dat %>% filter(!is.infinite(ln.rps), !is.na(ln.rps), broodYr>=1950, broodYr<=2010)
dim(dat.2)

#Create Output Data Object: Species, Stock, broodYr, covYr, covar, mean, med, low.95, low.50,high.50, high.95
effect.df.names <- c('species', 'region','sub.region','stock', 'broodYr', 'covYr', 'covar', 
                     'mean', 'low.95', 'low.50', 'med','high.50', 'high.95')
effect.df <- data.frame(array(dim=c(0,length(effect.df.names))))
names(effect.df) <- effect.df.names


#Parameter List
par.df.names <- c('species', 'region','sub.region','stock', 'par',
                    'mean', 'low.95', 'low.50', 'med','high.50', 'high.95')
par.df <- data.frame(array(dim=c(0,length(par.df.names))))
names(par.df) <- par.df.names

pars <- c('alpha','beta','beta.std','sigma_oe','sigma_pdo', 'sigma_npgo')

#Read STAN Models ==============================================
if(read.rds==TRUE) {
start <- date()
s <- 1
for(s in 1:n.species) {

  print(paste('###### s',s,'of',n.species))
  temp.species <- species[s]
  dat.3 <- dat.2 %>% filter(species==temp.species)
  stocks <- unique(dat.3$stock)
  n.stocks <- length(stocks)

  p <- 1
  for(p in 1:n.stocks) {
    print(paste('# p',p,'of',n.stocks))
    temp.stock <- stocks[p]
    dat.input <- dat.3 %>% filter(stock==temp.stock) %>% arrange(broodYr)
  
    #Determine Covariate Years
    covar.year <- dat.input$broodYr + offset[s]
    
    #Load Data
    if(file.exists(file.path(dir.output,paste0(temp.species,'-',temp.stock,'-fit.rds')))) {
      
      fit <- readRDS(file=file.path(dir.output,paste0(temp.species,'-',temp.stock,'-fit.rds')))
    
      #Extract coefficients
      ex <- rstan::extract(fit)
      
      #EFFECTS
      coefs <- ex$coef
    
      #PDO df
      temp.pdo <- data.frame(temp.species, dat.input$region, dat.input$sub.region, temp.stock, 
                             dat.input$broodYr, covar.year,
                            'PDO', apply(coefs[,,1], 2, mean),
                             t(apply(coefs[,,1],2, quantile, probs=c(0.025,0.25,0.5,0.75,0.975))))
      names(temp.pdo) <- effect.df.names
      #NPGO df
      temp.npgo <- data.frame(temp.species, dat.input$region, dat.input$sub.region, temp.stock, 
                             dat.input$broodYr, covar.year,
                             'NPGO', apply(coefs[,,2], 2, mean),
                             t(apply(coefs[,,2],2, quantile, probs=c(0.025,0.25,0.5,0.75,0.975))))
      names(temp.npgo) <- effect.df.names
      #Append
      effect.df <- rbind(effect.df, temp.pdo, temp.npgo)
      
      #PARAMETERS
      
      vals <- rbind(c(mean(ex$alpha), quantile(ex$alpha, probs=c(0.025,0.25,0.5,0.75,0.975))),
                    c(mean(ex$beta), quantile(ex$beta, probs=c(0.025,0.25,0.5,0.75,0.975))),
                    c(mean(ex$beta), quantile(ex$beta, probs=c(0.025,0.25,0.5,0.75,0.975)))*mean(dat.input$spawn),
                    c(mean(ex$sigma_oe), quantile(ex$sigma_oe, probs=c(0.025,0.25,0.5,0.75,0.975))),
                    c(mean(ex$sigma_pe[,1]), quantile(ex$sigma_pe[,1], 
                                                      probs=c(0.025,0.25,0.5,0.75,0.975))),
                    c(mean(ex$sigma_pe[,2]), quantile(ex$sigma_pe[,2], 
                                                      probs=c(0.025,0.25,0.5,0.75,0.975))))
      temp.par.df <- data.frame(pars, vals)
      temp.par.df <- data.frame(temp.stock, temp.par.df)
      temp.par.df <- data.frame(unique(dat.input$sub.region), temp.par.df)
      temp.par.df <- data.frame(unique(dat.input$region), temp.par.df)
      temp.par.df <- data.frame(temp.species, temp.par.df)
      #Rename
      names(temp.par.df) <- par.df.names
      #Append
      par.df <- rbind(par.df,temp.par.df)
    }#end if
  }#next p
}#next s

#Write Output
write.csv(effect.df, file.path(dir.output,"effect.df.csv"), row.names=FALSE)
write.csv(par.df, file.path(dir.output,"par.df.csv"), row.names=FALSE)
}else {
  effect.df <- read.csv(file=file.path(dir.output,"effect.df.csv"))
  par.df <- read.csv(file=file.path(dir.output,"par.df.csv"))
}




#FIGURES ===============================================

#Plot Effects Timeseries ===============================
pdf(file.path(dir.figs, "Effect Timeseries_stdx.pdf"), height=8, width=11)
s <- 1
for(s in 1:n.species) {
  #PDO
  covar <- 'PDO'
  g <- ggplot(effect.df[effect.df$species==species[s] & 
                          effect.df$covar==covar,], 
              aes(x=covYr, y=med, color=region, group=stock)) +
        theme_bw() +
        # scale_color_gdocs() +
        facet_wrap(~region, scales='free_y') +
        geom_line() +
        theme(legend.position='na')+
        geom_hline(yintercept=0) +
         ggtitle(paste(covar,'Effect'), subtitle=species[s]) +
         xlab('Covariate Year') +
         ylab('Median Effect')
  g
  plot(g)
  
  #NPGO
  covar <- 'NPGO'
  g <- ggplot(effect.df[effect.df$species==species[s] & 
                          effect.df$covar==covar,], 
              aes(x=covYr, y=med, color=region, group=stock)) +
    theme_bw() +
    facet_wrap(~region, scales='free_y') +
    geom_line() +
    theme(legend.position='na')+
    geom_hline(yintercept=0) +
    ggtitle(paste(covar,'Effect'), subtitle=species[s]) +
    xlab('Covariate Year') +
    ylab('Median Effect')
  g
  plot(g)

}#next s
dev.off()

#Geom_smooth
pdf(file.path(dir.figs, "Effect Timeseries_smooth.pdf"), height=8, width=11)
s <- 1
for(s in 1:n.species) {
  #PDO
  covar <- 'PDO'
  g <- ggplot(effect.df[effect.df$species==species[s] & 
                          effect.df$covar==covar,], 
              aes(x=covYr, y=med, color=region, group=stock)) +
    theme_bw() +
    # scale_color_gdocs() +
    facet_wrap(~region, scales='free_y') +
    geom_line(alpha=0.5) +
    theme(legend.position='na')+
    geom_hline(yintercept=0) +
    ggtitle(paste(covar,'Effect'), subtitle=species[s]) +
    xlab('Covariate Year') +
    ylab('Median Effect') +
    geom_smooth(data=effect.df[effect.df$species==species[s] & 
                            effect.df$covar==covar,], aes(x=covYr, y=med, group=region), color='black', lwd=0.5)
  g
  plot(g)
  
  #NPGO
  covar <- 'NPGO'
  g <- ggplot(effect.df[effect.df$species==species[s] & 
                          effect.df$covar==covar,], 
              aes(x=covYr, y=med, color=region, group=stock)) +
    theme_bw() +
    facet_wrap(~region, scales='free_y') +
    geom_line(alpha=0.5) +
    theme(legend.position='na')+
    geom_hline(yintercept=0) +
    ggtitle(paste(covar,'Effect'), subtitle=species[s]) +
    xlab('Covariate Year') +
    ylab('Median Effect') +
    geom_smooth(data=effect.df[effect.df$species==species[s] & 
                                 effect.df$covar==covar,], aes(x=covYr, y=med, group=region), color='black', lwd=0.5)
  g
  plot(g)
  
}#next s

dev.off()


#Plot Parameters ===============================
names(par.df)

#Pars-pannels, Region Colors, Species Pages
pdf(file.path(dir.figs,'Pars_1.pdf'), height=6, width=9)
s <- 1
for(s in 1:n.species) {

  g <- ggplot(par.df[par.df$species==species[s],], 
              aes(y=med, x=region, fill=region)) +
    theme_linedraw() +
    # scale_fill_gdocs() +
    facet_wrap(~par, scales='free') +
    theme(axis.text.x=element_blank()) +
    geom_boxplot(outlier.shape=NA) +
    ggtitle(species[s]) +
    xlab('Region') +
    ylab('Distribution of Median Parameter Estimates')
  # g
  plot(g)
}#next s
dev.off()

#Page PArs
pdf(file.path(dir.figs,'Pars_2_freeX.pdf'), height=6, width=9)
p <- 1
for(p in 1:length(pars)) {
  
  g <- ggplot(par.df[par.df$par==pars[p],], 
              aes(y=med, x=region, fill=region)) +
    theme_linedraw() +
    # scale_fill_gdocs() +
    facet_wrap(~species, scales='free_x') +
    theme(axis.text.x=element_blank()) +
    geom_boxplot(outlier.shape=NA) +
    geom_jitter(alpha=0.3) +
    ggtitle(pars[p]) +
    xlab('Region') +
    ylab('Distribution of Median Parameter Estimates')
  # g
  plot(g)
}#next s
dev.off()

pdf(file.path(dir.figs,'Pars_2_free.pdf'), height=6, width=9)
p <- 1
for(p in 1:length(pars)) {
  
  g <- ggplot(par.df[par.df$par==pars[p],], 
              aes(y=med, x=region, fill=region)) +
    theme_linedraw() +
    # scale_fill_gdocs() +
    facet_wrap(~species, scales='free') +
    theme(axis.text.x=element_blank()) +
    geom_boxplot(outlier.shape=NA) +
    geom_jitter(alpha=0.3) +
    ggtitle(pars[p]) +
    xlab('Region') +
    ylab('Distribution of Median Parameter Estimates')
  # g
  plot(g)
}#next s
dev.off()

#Grouped aCross Regions

pdf(file.path(dir.figs,'Pars_3.pdf'), height=6, width=9)

  g <- ggplot(par.df, 
              aes(y=med, x=species, fill=species)) +
    theme_linedraw() +
    scale_fill_gdocs() +
    facet_wrap(~par, scales='free', ncol=length(pars)) +
    theme(legend.position='top', axis.text.x=element_blank()) +
    geom_boxplot(outlier.shape=NA) +
    geom_jitter(alpha=0.2) +
    ggtitle('Parameter Comparison') +
    xlab('Species') +
    ylab('Distribution of Median Parameter Estimates')
  g
  plot(g)
  
  g <- ggplot(filter(par.df, par %in% c('alpha','beta','beta.std','sigma_oe')), 
              aes(y=med, x=species, fill=species)) +
    theme_linedraw() +
    scale_fill_gdocs() +
    facet_wrap(~par, scales='free', ncol=length(pars)) +
    theme(legend.position='top', axis.text.x=element_blank()) +
    geom_boxplot(outlier.shape=NA) +
    geom_jitter(alpha=0.2) +
    ggtitle('Parameter Comparison') +
    xlab('Species') +
    ylab('Distribution of Median Parameter Estimates')
  g
  plot(g)
  
  g <- ggplot(par.df, 
              aes(y=med, x=species, fill=species)) +
    theme_linedraw() +
    scale_fill_gdocs() +
    facet_wrap(~par, scales='free', ncol=length(pars)) +
    theme(legend.position='na', axis.text.x=element_blank()) +
    geom_violin(lwd=FALSE, alpha=0.5) +
    geom_boxplot(outlier.shape=NA, width=0.3) +
    # geom_jitter(alpha=0.2) +
    ggtitle('Parameter Comparison') +
    xlab('Species') +
    ylab('Distribution of Median Parameter Estimates')
  g
  plot(g)
  
dev.off()




