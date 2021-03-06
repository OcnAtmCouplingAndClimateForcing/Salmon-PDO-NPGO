#==================================================================================================
#Project Name: COAST-WIDE SALMON RECRUITMENT ANALYSIS - Fit DLM-Ricker with Time-varying Components Shared Among Regions
#Creator: Curry James Cunningham, NOAA/NMFS, ABL
#Date: 8.20.18
#
#Purpose: To fit a dynamic linear Ricker model with 
#
#
#==================================================================================================
#NOTES:
#  a) Corrected problem with Chinook salmon - Lewis River had data sets for both Late-fall and Fall run 
#       with the same "stock" identifier. I updated these in the input .xlsx and .csv files.
# 
# TIMINGS: 
# [1] "n.iter: 30000"
# [1] "n.thin: 15"
# [1] "Sat Oct 26 13:52:39 2019"
# [1] "Sun Oct 27 15:11:44 2019"

# Change in AR coefficient
# [1] "n.iter: 30000"
# [1] "n.thin: 15"
# [1] "Sat Mar 21 19:29:08 2020"
# [1] "Sun Mar 22 11:32:00 2020"

# [1] "n.iter: 30000"
# > print(paste('n.thin:',n.thin))
# [1] "n.thin: 10"
# > print(start)
# [1] "Sun Mar 22 11:51:38 2020"
# > print(end)
# [1] "Mon Mar 23 07:32:44 2020"

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
require(reshape2)

#CONTROL SECTION ==========================================================
do.est <- FALSE

n.chains <- 4
n.iter <- 3e4#2e4
n.thin <- 10
(n.iter/n.thin*0.5)*n.chains

# Model Designations
# Original
# model_file <- "DLM-Region.stan"
# model_name <- "DLM-Region"

model_file <- "DLM-Region-ARerror.stan"
model_name <- "DLM-Region-ARerror-UPDATED"

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
dat <- read.csv(file.path("data","AK-WCoast-Salmon-SR.csv"), header=TRUE, stringsAsFactors=FALSE)
#Add rps and ln.rps
dat$rps <- dat$rec/dat$spawn
dat$ln.rps <- log(dat$rps)

# How many total datasets
# species.stock <- paste(dat$species, dat$stock)

# length(unique(species.stock))

# nums.tab <- dat %>% group_by(species) %>% summarize(n.pops=length(unique(stock)))
# View(nums.tab); sum(nums.tab$n.pops)
#Subset SR Data ==============================================
# dat.2 <- dat %>% filter(!is.infinite(ln.rps), !is.na(ln.rps), broodYr>=1950, broodYr<=2010)
# dat.2 <- dat %>% filter(!is.infinite(ln.rps), !is.na(ln.rps), !is.nan(ln.rps), rec>0, broodYr>=1950, broodYr<=2010)
dat.2 <- dat %>% filter(!is.infinite(ln.rps), !is.na(ln.rps), !is.nan(ln.rps), rec>0, broodYr>=1950, broodYr<=2010)

#Update Duplicate Stock Identifiers (mostly Chinook) ==========================
# Lewis River Fall-run and Late fall are the sticky wickets.

# dat.2$stock[dat.2$species=="Chinook" & dat.2$stock=="Lewis River",] <- paste(dat.2$run[dat.2$species=="Chinook" & dat.2$stock=="Lewis River",],"-", dat.2$stock[dat.2$species=="Chinook" & dat.2$stock=="Lewis River",])

# sum.id <- dat.2 %>% group_by(species, stock) %>% summarize(n.stockid=length(unique(stock.id)))
# View(sum.id)

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
for(s in 1:n.species) {
# for(s in c(1,3,4,5)) { #Currently Excluding Chinook
  print(paste('###### s',s,'of',n.species))
  temp.species <- species[s]
  dat.3 <- dat.2 %>% filter(species==temp.species)
  stocks <- unique(dat.3$stock)
  n.stocks <- length(stocks)
  
  #Required Attributs of the species
  stock.regions <- unique(dat.3$large.region)
  n.stock.regions <- length(stock.regions)
  
  stock.years <- min(dat.3$broodYr):max(dat.3$broodYr)
  n.stock.years <- length(stock.years)
  
  
  #Create Data Objects
  maxN <- n.stock.years
  S <- n.stocks
  N <- vector(length=n.stocks)
  R <- n.stock.regions #Number of Regions
  region <- vector(length=n.stocks)
  K <- 2 #Number of covariates PDO, NPGO
  # covars <- array(dim=c(n.stocks,100,K)) #We will start with a temporary length of 100 years then trim down to the max N
  PDO <- array(data=0, dim=c(n.stocks,maxN))
  NPGO <- array(data=0, dim=c(n.stocks,maxN))
  
  #Year Pointer - For Referencing Coefficient Values
  years <- array(data=0, dim=c(n.stocks,maxN))
  
  #Ricker Parameters
  ln_rps <- array(data=0, dim=c(n.stocks, maxN))
  spawn <- array(data=0, dim=c(n.stocks, maxN))

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
    
    years[p,1:N[p]] <- which(stock.years %in% dat.input$broodYr )
      
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
  # maxN <- max(N)
  temp.regions <- regions[unique(region)]

  #Truncate STAN Input Objects ======================
  # ln_rps <- ln_rps[,1:maxN]
  # spawn <- spawn[,1:maxN]
  # # covars <- covars[,1:maxN,]
  # PDO <- PDO[,1:maxN]
  # NPGO <- NPGO[,1:maxN]
  # 
  #Call STAN =======================================================
  if(do.est==TRUE) {
  fit <- stan(file=file.path(dir.R, model_file),
              model_name=model_name,
              data=list("N"=N, "maxN"=maxN,
                        "ln_rps"=ln_rps, "spawn"=spawn,
                        "K"=K, 
                        #"covars"=covars,
                        "PDO"=PDO, "NPGO"=NPGO,
                        "S"=S,
                        "R"=R,
                        "region"=region, "years"=years),
              chains=n.chains, iter=n.iter, thin=n.thin,
              cores=n.chains, verbose=FALSE,
              seed=101,
              # sample_prior=TRUE,
              control = list(adapt_delta = 0.99))
  
  #Save Output
  saveRDS(fit, file=file.path(dir.output,paste0(temp.species,'-fit.rds')))
  # Save Summary
  write.csv(summary(fit)$summary, file=file.path(dir.figs, "summary.csv"))
  }else {
    fit <- readRDS(file=file.path(dir.output,paste0(temp.species,'-fit.rds')))
  }
  #Save Extras
  extras <- NULL
  extras$stocks <- stocks
  extras$n.stocks <- n.stocks
  extras$stock.regions <- stock.regions
  extras$n.stock.regions <- n.stock.regions
  extras$stock.years <- stock.years
  extras$n.stock.years <- n.stock.years
  saveRDS(extras, file=file.path(dir.output,paste0(temp.species,'-extras.rds')))
  
  # fit <- readRDS(file=file.path(dir.output,paste0(temp.species,'-',temp.stock,'-fit.rds')))
  
  # Evaluate Output
  # fit
  
  # stan_trace(fit, pars=list('beta'))
  # # stan_par(fit)
  # ## extract samples as a list of arrays
  # pars <- rstan::extract(fit)
  # 
  # #Quickly Plot coefficients over time.
  # coefs.pdo <- apply(pars$coef_PDO, c(2,3), quantile, probs=c(0.025,0.25,0.5,0.75,0.975))
  # coefs.npgo <- apply(pars$coef_NPGO, c(2,3), quantile, probs=c(0.025,0.25,0.5,0.75,0.975))
  # 
  # dimnames(coefs.pdo)[[2]] <- stock.regions
  # dimnames(coefs.pdo)[[3]] <- stock.years
  # 
  # dimnames(coefs.npgo)[[2]] <- stock.regions
  # dimnames(coefs.npgo)[[3]] <- stock.years
  # 
  # list.coefs.pdo <- melt(coefs.pdo)
  # head(list.coefs.pdo)
  # names(list.coefs.pdo) <- c('quantile')
  # 
  # g <- ggplot(filter(list.coefs.pdo, Var1=='50%'), aes(x=Var3, y=value, color=factor(Var2))) +
  #        geom_line()
  # g
  
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
  # shinystan::launch_shinystan(fit)
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

print(paste('n.iter:',n.iter))
print(paste('n.thin:',n.thin))
print(start)
print(end)


