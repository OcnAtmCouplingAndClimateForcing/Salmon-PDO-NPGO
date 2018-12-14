#==================================================================================================
#Project Name: COAST-WIDE SALMON RECRUITMENT ANALYSIS - Get PDO/NPGO Data
#Creator: Curry James Cunningham, NOAA/NMFS, ABL
#Date: 8.15.18
#
#Purpose: To download latest PDO and NPGO data
#
#
#==================================================================================================
#NOTES:
#
#==================================================================================================
require(rpdo)
require(rsoi)
require(tidyverse)
require(dplyr)

#Define Workflow Paths ====================================================
# *Assumes you are working from the Coastwide_Salmon_Analysis R project
wd <- getwd()
dir.data <- file.path(wd,"Data")
dir.figs <- file.path(wd,"figs")
dir.output <- file.path(wd,"Output")

#Load PDO/NPGO =====================================
pdo.dat <- rpdo::download_pdo()
npgo.dat <- rsoi::download_npgo()

head(pdo.dat)
head(npgo.dat)

#Create Annual Averages
pdo.avg <- pdo.dat %>% group_by(Year) %>% summarize('avg'=mean(PDO)) %>% 
             mutate('std'= (avg-mean(avg))/sd(avg) )

npgo.avg <- npgo.dat %>% group_by(Year) %>% summarize('avg'=mean(NPGO)) %>% 
              mutate('std'= (avg-mean(avg))/sd(avg) )

pdo.avg.winter <- pdo.dat %>% filter(Month<5) %>% group_by(Year) %>% 
                    summarize('avg'=mean(PDO)) %>% 
                    mutate('std'= (avg-mean(avg))/sd(avg) )

npgo.avg.winter <- npgo.dat %>% filter(Month %in% c('Jan','Feb','Mar','Apl')) %>% group_by(Year) %>% 
                     summarize('avg'=mean(NPGO)) %>% 
                     mutate('std'= (avg-mean(avg))/sd(avg) )

#Plotting Timeseries ============================================
png(file.path(dir.figs,'PDO-NPGO.png'), height=7, width=8, units='in', res=500)

par(mfrow=c(2,1), oma=c(2,2,1,1), mar=c(2,2,2,0))
#Trial plotting
plot(std~Year, data=pdo.avg, type='l', col='blue',
     ylab='Standardized Value', main='Annual Mean')
lines(std~Year, data=npgo.avg, type='l', col='red')
legend('bottomleft', legend=c('PDO','NPGO'), text.col=c('blue','red'), bty='n', ncol=2)

#Plotting winter values
plot(std~Year, data=pdo.avg.winter, type='l', col='blue',
     ylab='Standardized Value',main='Jan-Mar Mean')
lines(std~Year, data=npgo.avg.winter, type='l', col='red')


mtext('Standardized Value', side=2, outer=TRUE, font=2, line=0.5)
mtext('Year', side=1, outer=TRUE, font=2, line=0.5)

dev.off()

#Combine and Save Output ====================================
pdo.avg$name <- 'pdo'
pdo.avg.winter$name <- 'pdo.w'
npgo.avg$name <- 'npgo'
npgo.avg.winter$name <- 'npgo.w'

covar.dat <- rbind(pdo.avg, pdo.avg.winter, npgo.avg, npgo.avg.winter)

write.csv(covar.dat, file=file.path(dir.data,"Covariates","covar.dat.csv"), row.names=FALSE)
