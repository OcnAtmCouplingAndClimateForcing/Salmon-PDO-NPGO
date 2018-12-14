#==================================================================================================
#Project Name: COAST-WIDE SALMON RECRUITMENT ANALYSIS - Fit MARSS-Ricker
#Creator: Curry James Cunningham, NOAA/NMFS, ABL
#Date: 8.20.18
#
#Purpose: To fit a dynamic linear Ricker model with MARSS
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
# require(rstan)
# require(rstantools)
# require(bayesplot)
# require(shinystan)
require(BEST)
require(MARSS)

#Define Workflow Paths ====================================================
# *Assumes you are working from the Coastwide_Salmon_Analysis R project
wd <- getwd()
dir.data <- file.path(wd,"Data")
dir.figs <- file.path(wd,"Figs","MARSS-Ricker")
dir.output <- file.path(wd,"Output","MARSS-Ricker")
dir.R <- file.path(wd,"R")

#Create
dir.create(dir.figs, recursive=TRUE)
dir.create(dir.output, recursive=TRUE)

#Source Necessary Functions =====================================
# source(file.path(dir.R, "plot-model-fit.R"))
# source(file.path(dir.R, "plot-coef-ts.R"))
# source(file.path(dir.R, "plot-other-pars.R"))
source(file.path(dir.R, "MARSS-Ricker2.R"))

#Read SR Data ===================================================
dat <- read.csv(file.path(dir.data,"full.dat.4.csv"), header=TRUE, stringsAsFactors=FALSE)

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
dat.2 <- dat %>% filter(!is.infinite(ln.rps), !is.na(ln.rps), broodYr>=1950, broodYr<=2010)
dim(dat.2)

#Create Output Object ==========================================
vect.n.stocks <- vector(length=n.species)
for(s in 1:n.species) { 
  temp.species <- species[s]
  dat.3 <- dat.2 %>% filter(species==temp.species)
  stocks <- unique(dat.3$stock)
  vect.n.stocks[s] <-  length(stocks)
}

out.list <- list(vector('list', length=vect.n.stocks[1]),
                 vector('list', length=vect.n.stocks[2]),
                 vector('list', length=vect.n.stocks[3]),
                 vector('list', length=vect.n.stocks[4]))

#Fit MARSS Models ==============================================
start <- date()
s <- 4
# for(s in 1:n.species) {
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
    # ln_rps <- dat.input$ln.rps
    # spawn <- dat.input$spawn
    K <- 2
    covars <- array(dim=c(N,K))  #[N,K]
    #Collect covar data
    n <- 1
    for(n in 1:N) {
      covar.year <- dat.input$broodYr[n] + offset[s]
      covars[n,1] <- covar.dat$avg[covar.dat$name=='pdo' & covar.dat$Year==covar.year] #PDO
      covars[n,2] <- covar.dat$avg[covar.dat$name=='npgo' & covar.dat$Year==covar.year] #NPGO
    }#next n
  
    #FIT MARSS MODEL
    out <- MARSS_Ricker2(dat.input=dat.input, covars=covars)
    saveRDS(out, file=file.path(dir.output,paste0(temp.species,'-',temp.stock,'-out.rds')))
 
  
  }#next p - stock
# }#next s - species

end <- date()




