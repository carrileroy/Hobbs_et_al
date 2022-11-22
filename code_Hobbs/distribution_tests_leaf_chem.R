#Ridgeline plots and Z-tests for Hobbs et al. 2022

# Data in hobbs.csv, mass.csv, N.csv

# Update R and RStudio: 
install.packages("installr")
library(installr)
updateR()
# To update RStudio, go to "Help" and click "Check for Updates"
# To update packages, go to "Tools" and "Check for Package Updates"

if(!require(ggridges)){install.packages("ggridges")}
if(!require(ggplot2)){install.packages("ggplot2")}
if(!require(BSDA)){install.packages("BSDA")}
if(!require(twosamples)){install.packages("twosamples")}
if(!require(gridExtra)){install.packages("gridExtra")}

library(ggridges)
library(ggplot2)
library(BSDA)
library(twosamples)
library(gridExtra) 

#######
# Molecule Mass
#######
#z-test: requires data be in two separate vectors, and that you know the std dev for each (sigma): data = mass.csv

vec1 <- mass$F
vec2 <- mass$M

z.test(x=vec1, y=vec2, mu=0, sigma.x=178.81, sigma.y=147.12, alternative = "two.sided", conf.level=0.95)


#Monte Carlo Kolmogorov-Smirnov test: requires data be in two separate vectors: data = mass.csv

ks.test(vec1, vec2, alternative = "two.sided", exact = TRUE, simulate.p.value = TRUE, B = 2000)

#######
#% Nitrogen
#######
#z-test: requires data be in two separate vectors, and that you know the std dev for each (sigma): data = mass.csv

vec3 <- N$F
vec4 <- N$M

z.test(x=vec3, y=vec4, mu=0, sigma.x=6.604, sigma.y=7.522, alternative = "two.sided", conf.level=0.95)


#Monte Carlo Kolmogorov-Smirnov test: requires data be in two separate vectors: data = mass.csv

ks.test(vec3, vec4, alternative = "two.sided", exact = TRUE, simulate.p.value = TRUE, B = 2000)

#######
# Retention Time
#######
#z-test: requires data be in two separate vectors, and that you know the std dev for each (sigma): data = mass.csv

vec5 <- RT$F
vec6 <- RT$M

z.test(x=vec5, y=vec6, mu=0, sigma.x=5.86, sigma.y=4.11, alternative = "two.sided", conf.level=0.95)


#Monte Carlo Kolmogorov-Smirnov test: requires data be in two separate vectors: data = mass.csv

ks.test(vec5, vec6, alternative = "two.sided", exact = TRUE, simulate.p.value = TRUE, B = 2000)


########
# Ridgeline plots
########

#plots require data in one column: hobbs.csv

hobbs$Sex = factor(hobbs$Sex,
                   levels=unique(hobbs$Sex))

str(hobbs)

plotM <- ggplot(hobbs, aes(x = mass, y = Sex, fill = Sex)) +
  geom_density_ridges() +
  theme_ridges() + 
  theme(legend.position = "none") +
  ylab("Willow Sex") + xlab("Molecular Mass") +
  scale_fill_manual(values = c("gray", "gray50"))
plotM

ggsave(file="Fig4.pdf", plotM, width=18.2, height=14, units="cm", dpi=600)

plotRT <- ggplot(hobbs, aes(x = RT, y = Sex, fill = Sex)) +
  geom_density_ridges() +
  theme_ridges() + 
  theme(legend.position = "none") +
  ylab("Willow Sex") + xlab("Retention Time") +
  scale_fill_manual(values = c("gray", "gray50"))
plotRT

ggsave(file="Fig5.pdf", plotRT, width=18.2, height=14, units="cm", dpi=600)

plotN <- ggplot(hobbs, aes(x = X.N, y = Sex, fill = Sex)) +
  geom_density_ridges() +
  theme_ridges() + 
  theme(legend.position = "none") +
  ylab("Willow Sex") + xlab("% N") +
  scale_fill_manual(values = c("gray", "gray50"))
plotN

ggsave(file="Fig6.pdf", plotN, width=18.2, height=14, units="cm", dpi=600)

gM <- ggplotGrob(plotM)
gN <- ggplotGrob(plotN)
g <- arrangeGrob(gM, gN, ncol=2)
#
ggsave(file="Fig4-update.pdf", g, width=20, height=12, units="cm", dpi=800)