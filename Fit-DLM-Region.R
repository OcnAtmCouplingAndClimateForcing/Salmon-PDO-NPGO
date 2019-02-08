#==================================================================================================
#Project Name: COAST-WIDE SALMON RECRUITMENT ANALYSIS - Fit DLM-Ricker
#Creator: Curry James Cunningham, NOAA/NMFS, ABL
#Date: 8.20.18
#
#Purpose: To fit a dynamic linear Ricker model with 
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
dir.figs <- file.path(wd,"Figs","DLM-Region")
dir.output <- file.path(wd,"Output","DLM-Region")
dir.R <- file.path(wd,"R")

#Create
dir.create(dir.figs, recursive=TRUE)
dir.create(dir.output, recursive=TRUE)

#Source Necessary Functions =====================================
# source(file.path(dir.R, "plot-model-fit.R"))
# source(file.path(dir.R, "plot-coef-ts.R"))
# source(file.path(dir.R, "plot-other-pars.R"))

#Read SR Data ===================================================
dat <- read.csv(file.path("data","AK-WCoast-Salmon-SR.csv"), header=TRUE, stringsAsFactors=TRUE)
#Add rps and ln.rps
dat$rps <- dat$rec/dat$spawn
dat$ln.rps <- log(dat$rps)

#Subset SR Data ==============================================
dat.2 <- dat %>% filter(!is.infinite(ln.rps), !is.na(ln.rps), broodYr>=1950, broodYr<=2010)

#Read PDO/NPGO Data =============================================
covar.dat <- read.csv(file.path(dir.data,"Covariates","covar.dat.csv"))

#Extract Metadata ===================================================
species <- unique(dat.2$species)
n.species <- length(species)

regions <- unique(dat.2$large.region)
n.regions <- length(regions)

#Specify Offsets for Covariates ================================
#From BroodYear
offset <- c(2,2,1,2,3)
offset.table <- cbind(species, offset)
write.csv(offset.table, file=file.path(dir.figs,"offset.table.csv"))

#Fit STAN Models ==============================================
start <- date()
s <- 1
# for(s in 1:n.species) {
# for(s in 2:n.species) {  
  print(paste('###### s',s,'of',n.species))
  temp.species <- species[s]
  dat.3 <- dat.2 %>% filter(species==temp.species)
  stocks <- unique(dat.3$stock)
  n.stocks <- length(stocks)
  
  stock.regions <- unique(dat.3$large.region)
  n.stock.regions <- length(stock.regions)
  
  #Create Data Objects
  S <- n.stocks
  N <- vector(length=n.stocks)
  R <- n.stock.regions #Number of Regions
  region <- vector(length=n.stocks)
  K <- 2 #Number of covariates PDO, NPGO
  # covars <- array(dim=c(n.stocks,100,K)) #We will start with a temporary length of 100 years then trim down to the max N
  PDO <- array(data=0, dim=c(n.stocks,100))
  NPGO <- array(data=0, dim=c(n.stocks,100))
  
  
  #Ricker Parameters
  ln_rps <- array(data=0, dim=c(n.stocks, 100))
  spawn <- array(data=0, dim=c(n.stocks, 100))

  p <- 1
  for(p in 1:n.stocks) {
  
    #Retreive Data =================
    temp.stock <- stocks[p]
    dat.input <- dat.3 %>% filter(stock==temp.stock) %>% arrange(broodYr)
  
    #Assign STAN Inputs ============
    N[p] <- nrow(dat.input)
    region[p] <- which(stock.regions==unique(dat.input$large.region))
  
    ln_rps[p, 1:N[p]] <- dat.input$ln.rps 
    spawn[p, 1:N[p]] <- dat.input$spawn
  
    #Assign Covars ===============
    n <- 1
    for(n in 1:N[p]) {
      covar.year <- dat.input$broodYr[n] + offset[s]
      # covars[p,n,1] <- covar.dat$avg[covar.dat$name=='pdo' & covar.dat$Year==covar.year] #PDO
      # covars[p,n,2] <- covar.dat$avg[covar.dat$name=='npgo' & covar.dat$Year==covar.year] #NPGO
      PDO[p,n] <- covar.dat$avg[covar.dat$name=='pdo' & covar.dat$Year==covar.year] #PDO
      NPGO[p,n] <- covar.dat$avg[covar.dat$name=='npgo' & covar.dat$Year==covar.year] #NPGO
    }#next n
  }#next p

  #Determine maximum length of covariates =====================
  maxN <- max(N)
  temp.regions <- regions[unique(region)]

  #Truncate STAN Input Objects ======================
  ln_rps <- ln_rps[,1:maxN]
  spawn <- spawn[,1:maxN]
  # covars <- covars[,1:maxN,]
  PDO <- PDO[,1:maxN]
  NPGO <- NPGO[,1:maxN]
  
  
  #Call STAN =======================================================
  fit <- stan(file=file.path(dir.R,"DLM-Region.stan"),
              model_name="DLM-Region",
              data=list("N"=N, "maxN"=maxN,
                        "ln_rps"=ln_rps, "spawn"=spawn,
                        "K"=K, 
                        #"covars"=covars,
                        "PDO"=PDO, "NPGO"=NPGO,
                        "S"=S,
                        "R"=R,
                        "region"=region),
              chains=3, iter=5e3, thin=5,
              # chains=3, iter=1e3, thin=1,
              cores=3, verbose=FALSE,
              seed=101)
  #Save Output
  saveRDS(fit, file=file.path(dir.output,paste0(temp.species,'-fit.rds')))
  saveRDS(stock.regions, file=file.path(dir.output,paste0(temp.species,'-stock.regions.rds')))
  
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
  
# }#next p - stock
}#next s - species

end <- date()




