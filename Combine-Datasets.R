#==================================================================================================
#Project Name: COAST-WIDE SALMON RECRUITMENT ANALYSIS - Combine Datasets
#Creator: Curry James Cunningham, NOAA/NMFS, ABL
#Date: 8.8.18
#
#Purpose: To combine various salmon recruitment datasets into a single data object.
#
#
#==================================================================================================
#NOTES:
#  a) Clark and Frelund - More datasets available from the TotalAge.csv data file. 20 vs. 8
#  b) Brenner and Couture sockeye has Chignik Multiple times. Best to remove Chignik Watershed, and Chignik
#==================================================================================================
require(tidyverse)
require(ggthemes)
require(cowplot)
require(viridis)

wd <- getwd()

#Define Workflow Paths ====================================================
# *Assumes you are working from the Coastwide_Salmon_Analysis R project
dir.data <- file.path(wd,"Data")
dir.figs <- file.path(wd,"Figs","Data")
dir.output <- file.path(wd,"Output")

#Create
dir.create(dir.figs, recursive=TRUE)

#CONTROL SECTION ==========================================================
height <- 8
width <- 10

#LOAD DATASETS ============================================================

#KNB: Clark and Brenner 2017 =================================
cb.dat <- read.csv(file.path(dir.data,"Clark and Brenner 2017/Sockeye_BroodTables.csv"), header=TRUE)

#Calculate total Recruitment
TotalRec <- apply(cb.dat[,(which(names(cb.dat)=='TotalEscapement')+1):ncol(cb.dat)], 1, sum, na.rm=TRUE)
cb.dat <- data.frame(cb.dat,TotalRec)
cb.dat.2 <- cb.dat %>% filter(UseFlag==1) %>%  select(Species, Stock, Region, Sub.Region, BroodYear, TotalEscapement, TotalRec) %>% arrange(Region)
names(cb.dat.2)[5:7] <- c('broodYr','spawn','rec')

g <- ggplot(cb.dat.2, aes(y=Stock, x=broodYr, cex=rec/spawn, fill=Region)) +
      theme_bw() + 
      geom_point(alpha=0.5, pch=21, color='black') +
      ggtitle('Sockeye', subtitle='Clark and Brenner 2017')
g
ggsave(file.path(dir.figs,"Clark-Brenner-2017.pdf"), plot=g, height=height, width=width, unit='in')


#KNB: Brenner and Couture 2017 ====================================
bc.dat <- read.csv(file.path(dir.data,"Brenner and Couture 2017","SalmonSpawnRecruit.csv"))
names(bc.dat)[c(2,7:8)] <- c('broodYr','spawn','rec')
head(bc.dat)

bc.species <- unique(bc.dat$Species)
pdf(file.path(dir.figs,"Brenner-Couture-2017.pdf"), height=height, width=width)
for(s in bc.species) {
g <- ggplot(filter(bc.dat, Species==s), aes(y=Stock, x=broodYr, cex=rec/spawn, fill=Region)) +
             theme_bw() +
             geom_point(alpha=0.5, pch=21, color='black') +
             ggtitle(s, subtitle='Brenner and Couture 2017')
plot(g)
}#next s
dev.off()

#Clark and Frelund =======================================
cf.dat.age <- read.csv(file.path(dir.data,"Clark and Frelund","ChinookBroodTables_FW_SW_Age.csv"), header=TRUE)

cf.dat.total <- read.csv(file.path(dir.data,"Clark and Frelund","ChinookBroodTables_TotalAge.csv"), header=TRUE)

#Compare
sort(unique(cf.dat.age$Stock))
sort(unique(cf.dat.total$Stock))# More stocks in the total age dataset

head(cf.dat.total)

cf.dat.total.2 <- cf.dat.total %>% filter(UseFlag==1) %>% select(Stock.ID, Species, Stock, Region, BroodYear, TotalEscapement, TotalRecruits)
names(cf.dat.total.2)[c(5,6,7)] <- c('broodYr','spawn','rec')

g <- ggplot(cf.dat.total.2, aes(y=Stock, x=broodYr, cex=rec/spawn, fill=Region)) +
       theme_bw() + 
       geom_point(alpha=0.5, pch=21, color='black') +
       ggtitle('Chinook', subtitle='Clark-Frelund')
g
ggsave(file.path(dir.figs,"Clark-Frelund.pdf"), plot=g, height=height, width=width, unit='in')

#Malick and Cox 2016 =======================================

mc.pink.dat <- read.csv(file.path(dir.data,"Malick and Cox 2016","pink_data.csv"), header=TRUE)
mc.pink.dat.2 <- mc.pink.dat %>% filter(use==1) %>% select(-use)
names(mc.pink.dat.2)[6:8] <- c('broodYr','spawn','rec')

mc.chum.dat <- read.csv(file.path(dir.data,"Malick and Cox 2016","chum_data.csv"), header=TRUE)
mc.chum.dat.2 <- mc.chum.dat %>% filter(use==1) %>% select(stock.id, species, stock, region, sub.region, brood.yr, spawners, recruits)
names(mc.chum.dat.2)[6:8] <- c('broodYr','spawn','rec')

#Combine 
mc.dat.2 <- rbind(mc.pink.dat.2, mc.chum.dat.2)

#Plot
mc.species <- unique(mc.dat.2$species)
pdf(file.path(dir.figs,"Malick-Cox-2016.pdf"), height=height, width=width)
for(s in mc.species) {
  g <- ggplot(filter(mc.dat.2, species==s), aes(y=stock, x=broodYr, cex=rec/spawn, fill=region)) +
    theme_bw() +
    geom_point(alpha=0.5, pch=21, color='black') +
    ggtitle(s, subtitle='Malick and Cox 2016')
  plot(g)
}#next s
dev.off()


#Bring Them in The Darkness and Bind Them ============================

cb.dat.2 #Sockeye - 48 stocks
bc.dat
cf.dat.total.2 #Use for Chinook
mc.dat.2 #Pink and Chum 

length(unique(cb.dat.2$Stock[c]))

length(unique(mc.dat.2$stock[mc.dat.2$species=='Chum']))

length(unique(cf.dat.total.2$Stock[cf.dat.total.2$Species=='Chinook']))

names(cb.dat.2)
names(cf.dat.total.2)
names(mc.dat.2)


#Combine
full.dat <- cb.dat.2
full.dat.2 <- data.frame(NA,full.dat)
names(full.dat.2) <- c('stock.id','species','stock','region','sub.region',
                       'broodYr','spawn','rec')
#Attach CF
names(full.dat.2)
names(cf.dat.total.2)
cf.add <- data.frame(cf.dat.total.2[1:4],NA,cf.dat.total.2[5:ncol(cf.dat.total.2)])
names(cf.add) <- c('stock.id','species','stock','region','sub.region',
                   'broodYr','spawn','rec')
full.dat.3 <- rbind(full.dat.2,cf.add)

#Attach MC
names(full.dat.3)
names(mc.dat.2)
full.dat.4 <- rbind(full.dat.3, mc.dat.2)

#Write Combined Output File - Without Brenner and Couture ===============
#Change name of Alsek/Klukshu to Alsek
# full.dat.4$stock[full.dat.4$stock=="Alsek/Klukshu"] <- 'Alsek'
# levels(full.dat.4$stock)=="Alsek/Klukshu"
write.csv(full.dat.4, file=file.path(dir.data,"full.dat.4.csv"), row.names=FALSE)


# Dynamic Linear Model
require(mgcv)

temp.dat <- cb.dat.2 %>% filter(broodYr>=1950 & rec>0)
stocks <- sort(unique(temp.dat$Stock))
n.stocks <- length(stocks)
#Add Log R/S
temp.dat$ln.rps <- log(temp.dat$rec/temp.dat$spawn)


#Get Covariate Data
covar.dat <- read.csv(file.path(dir.data, "Covariates","covars.list.csv"))
npgo.dat <- covar.dat %>% filter(Covar=='NPGO')
pdo.dat <- covar.dat %>% filter(Covar=='PDO_Jeff')

offset <- 2

#Attach covariates to SR data
temp.dat$npgo <- NA
temp.dat$pdo <- NA
temp.dat$regime <- NA

i <- 1

for(i in 1:nrow(temp.dat)) {
  temp.dat$npgo[i] <- npgo.dat$value[npgo.dat$Year==(temp.dat$broodYr[i] + offset)]
  temp.dat$pdo[i] <- pdo.dat$value[pdo.dat$Year==(temp.dat$broodYr[i] + offset)]
  #
  if(temp.dat$broodYr[i]<1978) {
    temp.dat$regime[i] <- "<1978"
  }else {
    if(temp.dat$broodYr[i]<1989) {
      temp.dat$regime[i] <- "1978-1989"
    }else {
      temp.dat$regime[i] <- ">1989"
    }
  }
}

r <- S*exp(a-bs)

data <- filter(temp.dat, Stock=='Bear')

mod <- lm(ln.rps ~ spawn + regime:pdo + regime:npgo, data=data)



