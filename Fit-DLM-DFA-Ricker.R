#==================================================================================================
#Project Name: COAST-WIDE SALMON RECRUITMENT ANALYSIS - Fit DLM-DFA Ricker
#Creator: Curry James Cunningham, NOAA/NMFS, ABL
#Date: 8.20.18
#
#Purpose: To fit a dynamic linerarized Ricker with common Processes
#
#
#==================================================================================================
#NOTES:
#  a) Clark and Frelund - More datasets available from the TotalAge.csv data file. 20 vs. 8
#  b) Brenner and Couture sockeye has Chignik Multiple times. Best to remove Chignik Watershed, and Chignik
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

#Define Workflow Paths ====================================================
# *Assumes you are working from the Coastwide_Salmon_Analysis R project
wd <- getwd()
dir.data <- file.path(wd,"Data")
dir.figs <- file.path(wd,"Figs","DLM-DFA-Ricker")
dir.output <- file.path(wd,"Output","DLM-DFA-Ricker")
dir.R <- file.path(wd,"R")

#Create
dir.create(dir.figs, recursive=TRUE)
dir.create(dir.output, recursive=TRUE)

#Source Necessary Functions =====================================
# source(file.path(dir.R, "plot-model-fit.R"))
# source(file.path(dir.R, "plot-coef-ts.R"))
# source(file.path(dir.R, "plot-other-pars.R"))

#Read SR Data ===================================================
dat <- read.csv(file.path(dir.data,"AK-WCoast-Salmon-SR.csv"), header=TRUE, stringsAsFactors=FALSE)

#Add rps and ln.rps
dat$rps <- dat$rec/dat$spawn
dat$ln.rps <- log(dat$rps)

#Read PDO/NPGO Data =============================================
# covar.dat <- read.csv(file.path(dir.data, "Covariates","covars.list.csv"))
# npgo.dat <- covar.dat %>% filter(Covar=='NPGO')
# pdo.dat <- covar.dat %>% filter(Covar=='PDO_Jeff')

covar.dat <- read.csv(file.path(dir.data,"Covariates","covar.dat.csv"))
(covars <- unique(covar.dat$name))

#Remind me the start of time series
min(covar.dat$Year[covar.dat$name=='pdo'])
min(covar.dat$Year[covar.dat$name=='npgo'])

#Extract Metadata ===================================================
species <- unique(dat$species)
n.species <- length(species)


dat[dat$species=='',]

length(which(is.infinite(dat$ln.rps)==TRUE))
length(which(is.na(dat$ln.rps)==TRUE))

length(which(dat$spawn==0))
length(which(dat$rec==0))

sort(unique(dat$rps))

dat[is.infinite(dat$ln.rps),]
dat[is.na(dat$ln.rps),]

#Specify Offsets for Covariates ================================
#From BroodYear
offset <- c(2,2,1,2)#rep(2, n.species)
offset.table <- cbind(species, offset)
print(offset.table)
write.csv(offset.table, file=file.path(dir.figs,"offset.table.csv"))

#Subset SR Data ==============================================
dat.2 <- dat %>% filter(!is.infinite(ln.rps), !is.na(ln.rps), !is.nan(ln.rps), broodYr>=1950, broodYr<=2010)
dim(dat.2)

#Fit STAN Models ==============================================
start <- date()
s <- 2
# for(s in 1:n.species) {
# for(s in 2:n.species) {  
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
  
  #Collect Data
  N <- nrow(dat.input)
  ln_rps <- dat.input$ln.rps
  spawn <- dat.input$spawn
  K <- 2
  covars <- array(dim=c(N,K))  #[N,K]
  #Collect covar data
  n <- 1
  for(n in 1:N) {
    covar.year <- dat.input$broodYr[n] + offset[s]
    covars[n,1] <- covar.dat$avg[covar.dat$name=='pdo' & covar.dat$Year==covar.year] #PDO
    covars[n,2] <- covar.dat$avg[covar.dat$name=='npgo' & covar.dat$Year==covar.year] #NPGO
  }#next n
  
  #Fit the model
  fit <- stan(file=file.path(dir.R,"DLM-Ricker2.stan"),
              model_name="DLM-Ricker",
              data=list("N"=N, "ln_rps"=ln_rps, "spawn"=spawn,
                        "K"=K, "covars"=covars),
              chains=3, iter=5e4, thin=5,
              # chains=3, iter=1e3, thin=1,
              cores=3, verbose=FALSE,
              seed=101)
  #Save Output
  saveRDS(fit, file=file.path(dir.output,paste0(temp.species,'-',temp.stock,'-fit.rds')))
  
  # fit <- readRDS(file=file.path(dir.output,paste0(temp.species,'-',temp.stock,'-fit.rds')))
  
  # Evaluate Output
  # fit
  
  # stan_trace(fit, pars=list('beta'))
  # # stan_par(fit)
  # ## extract samples as a list of arrays
  # pars <- rstan::extract(fit)
  # 
  # #Plot parameters over time
  # med.coef <- apply(pars$coef,c(2,3), median)
  # plot(med.coef[,1])
  # 
  # # rstantools::bayes_R2(fit)
  # # saveRDS(fit, file=file.path(dir.figs,"fit.rds"))
  # 
  # #Plots
  # # bayesplot::mcmc_areas(a2)
  # # bayesplot::mcmc_areas_ridges()
  # 
  # #Shinystan
  # # sso <- shinystan::as.shinystan(fit, model_name="DLM-Ricker")
  # # shinystan::launch_shinystan(sso)
  # 
  # #Plot the fit
  # plot_model_fit(fit=fit, dat.input=dat.input)
  # 
  # #Plot Coefficients
  # plot_coef_ts(fit=fit, dat.input=dat.input, temp.offset=offset[s])
  # 
  # #Other Parameters
  # plot_other_pars(fit=fit)
  
  
  #Unload Unnecessary DLLs
  # loaded_dlls = getLoadedDLLs()
  # loaded_dlls = loaded_dlls[stringr::str_detect(names(loaded_dlls), '^file')]
  # if(length(loaded_dlls) > 5) {
  #   dll <- head(loaded_dlls)[1]
  #   dyn.unload(as.character(dll[[1]][2]))
  # }
  
}#next p - stock
# }#next s - species

end <- date()




