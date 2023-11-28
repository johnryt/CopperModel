# Import Data and require packages
library(plm)
library(stargazer)
library(lmtest)
library(lavaan)
library(readxl)
library(diptest)

# Set path, read data
wd <- ("C:/Users/ryter/Dropbox (MIT)/Group Research Folder_Olivetti/Displacement/00 Simulation/06 Module Integration/Data/scenario results/Multimodality")
setwd(wd)

# large50high <- read.csv('large50high.csv')
# large40high <- read.csv('large40high.csv')
# large30high <- read.csv('large30high.csv')
# large20high <- read.csv('large20high.csv')
# large10high <- read.csv('large10high.csv')
# largeall    <- read.csv('large0high.csv')
# 
# total50high <- read.csv('total50high.csv')
# total40high <- read.csv('total40high.csv')
# total30high <- read.csv('total30high.csv')
# total20high <- read.csv('total20high.csv')
# total10high <- read.csv('total10high.csv')
# totalall    <- read.csv('total0high.csv')

small50high <- read.csv('small50high.csv')
small40high <- read.csv('small40high.csv')
small30high <- read.csv('small30high.csv')
small20high <- read.csv('small20high.csv')
small10high <- read.csv('small10high.csv')
smallall    <- read.csv('small0high.csv')

cols <- colnames(small50high)[3:length(small50high)]
result50 = rep(0,length(cols))
result40 = rep(0,length(cols))
result30 = rep(0,length(cols))
result20 = rep(0,length(cols))
result10 = rep(0,length(cols))
sequ = seq(from=1,to=length(cols),by=1)
some_bimodality = rep(0,5)
significant_bimodality = rep(0,5)

for (i in sequ)
{
  result50[i] = dip.test(small50high[,cols[i]])$p.value
  result40[i] = dip.test(small40high[,cols[i]])$p.value
  result30[i] = dip.test(small30high[,cols[i]])$p.value
  result20[i] = dip.test(small20high[,cols[i]])$p.value
  result10[i] = dip.test(small10high[,cols[i]])$p.value
}
result50

cutoff = 0.1
some_bimodality[1] = length(result50[result50<=cutoff])/length(result50)*100
some_bimodality[2] = length(result40[result40<=cutoff])/length(result40)*100
some_bimodality[3] = length(result30[result30<=cutoff])/length(result30)*100
some_bimodality[4] = length(result20[result20<=cutoff])/length(result20)*100
some_bimodality[5] = length(result10[result10<=cutoff])/length(result10)*100

cutoff = 0.05
length(result50[result50<=cutoff])#/length(result)*100
significant_bimodality[1] = length(result50[result50<=cutoff])/length(result50)*100
significant_bimodality[2] = length(result40[result40<=cutoff])/length(result40)*100
significant_bimodality[3] = length(result30[result30<=cutoff])/length(result30)*100
significant_bimodality[4] = length(result20[result20<=cutoff])/length(result20)*100
significant_bimodality[5] = length(result10[result10<=cutoff])/length(result10)*100

print(some_bimodality)
print(significant_bimodality)




# cols[which(result<=cutoff)]


#x = c(0,1,2,3)
#y = c(0,0,0,0)
#y[x>1 & x<3] = 1
#x
#y
#m = as.data.frame(rbind(x,y))
#m[m[,'V2'] >0 & m[,'V2']<15,'V1'] = 1
#m
