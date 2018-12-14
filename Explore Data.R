#==================================================================================================
#Project Name: COAST-WIDE SALMON RECRUITMENT ANALYSIS - Data Exploration
#Creator: Curry James Cunningham, NOAA/NMFS, ABL
#Date: 12.6.17
#
#Purpose: To generate some exploratory figures describing the available data
#
#
#==================================================================================================
#NOTES:
#  a) 
#
#==================================================================================================
require(tidyverse)
require(ggthemes)
require(cowplot)
require(viridis)

wd <- getwd()

#Other directories
plot.dir <- paste0(wd,"/Figs")
output.dir <- paste0(wd,"/Output")

#Load Sockeye data

sock.dat <- read.csv("Data/Clark and Brenner 2017/Sockeye_BroodTables.csv", header=TRUE)
sock.meta <- read.csv("Data/Clark and Brenner 2017/Sockeye_StockInfo.csv", header=TRUE)

#Calculate total Recruitment
TotalRec <- apply(sock.dat[,(which(names(sock.dat)=='TotalEscapement')+1):ncol(sock.dat)], 1, sum, na.rm=TRUE)
sock.dat <- data.frame(sock.dat,TotalRec)

#========================================
#Filter

#Determine which flag to uste
length(which(sock.dat$UseFlag==0))
length(which(sock.dat$UseFlag==1))

sock.dat.2 <- sock.dat %>% filter(UseFlag==1, !is.na(TotalEscapement), TotalRec>0)


#========================================
#Metadata
n.data <- nrow(sock.dat.2)

stocks <- unique(sock.dat.2$Stock)
n.stocks <- length(stocks)




#========================================
#Exploratory plots


g <- ggplot(sock.dat.2, aes(x=TotalEscapement/1e6, y=TotalRec/1e6, color=BroodYear)) +
       theme_gray() +
       geom_point() +
       scale_color_viridis() +
       facet_wrap(~Stock, ncol=8, scales='free') +
       stat_smooth()

g
ggsave(paste0(plot.dir,"/Sockeye Data_1.png"), plot=g, height=10, width=15, units='in', dpi=600)


# Calculate Ratio 3/2 ocean ============================
# sock.dat.3$meanOage <- sock.dat$


#=======================================================
#Function to fit Ricker Model
sim



# 
# a <- seq(from=0, to=100, length.out=1000)
# b <- log(a)
# 
# expected <- 20
# sigma <- 0.75
# 
# 
# 
# norm.like <- dnorm(b, log(expected), sigma)
# 
# 
# lognorm.like <- dlnorm(a, log(expected), sigma)
# 
# require(BEST)
# par(mfrow=c(2,1), oma=c(3,0,3,0), mar=c(0,4,0,1))
# plot(x=a, y=norm.like, type='h')
# axis(3)
# plot(x=a, y=lognorm.like, type='h')
# 
# 
