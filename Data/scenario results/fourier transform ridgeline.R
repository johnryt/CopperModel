# Import Data and require packages
library(plm)
library(stargazer)
library(lmtest)
library(lavaan)
library(readxl)
library(ggridges)
install.packages('ggplot2')
library(ggplot2)

# Set path, read data
wd <- ("C:/Users/ryter/Dropbox (MIT)/Group Research Folder_Olivetti/Displacement/00 Simulation/06 Module Integration/Data/scenario results")
setwd(wd)
ft_rl <- read.csv("fourier transform ridgeline.csv")
ft_s <- read.csv('fourier transform ridgeline -50 opp.csv')[,c(1,2,3,8)]
ft_b <- read.csv('fourier transform ridgeline +50 opp.csv')[,c(1,2,3,8)]
ft_s$period = as.factor(ft_s$period)

ggplot(ft_s) +
  geom_density_ridges(mapping=aes(x=value,y=period,group=period,color=category,point_color=category,fill=category)) + 
  scale_y_discrete(expand = c(0,0)) + scale_x_continuous(expand = c(0,0)) +
  scale_fill_manual(values = c("#D55E0050", "#0072B250"), labels = c("female", "male")) 
  
