# Carri J. LeRoy, 8-10-22
#
# Two-way ANOVAs using permutation-Resampling!
#
# For Hobbs et al. litter paper: Initial Litter Chemistry, 2-way pANOVAs;
# Datasets = CHNTanninUpdate.csv (updated C:N to molar element ratios)
# Check to see if you need to install the lmPerm package and other packages needed for analysis and plotting
#
if(!require(agricolae)){install.packages("agricolae")}
if(!require(lmPerm)){install.packages("lmPerm")}
if(!require(psych)){install.packages("psych")}
if(!require(FSA)){install.packages("FSA")}
if(!require(ggplot2)){install.packages("ggplot2")}
if(!require(car)){install.packages("car")}
if(!require(multcompView)){install.packages("multcompView")}
if(!require(lsmeans)){install.packages("lsmeans")}
if(!require(tidyr)){install.packages("tidyr")}
if(!require(dplyr)){install.packages("dplyr")}
if(!require(magrittr)){install.packages("magrittr")}
if(!require(gridExtra)){install.packages("gridExtra")}
if(!require(tidyverse)){install.packages("tidyverse")}
if(!require(readxl)){install.packages("readxl")}
if(!require(devtools)){install.packages("devtools")}
if(!require(ggpubr)){install.packages("ggpubr")}
#
# Call the packages
library(lmPerm)           # permutative analyses
library(psych)            # stats for psych research (useful)
library(FSA)              # Fisheries stock assessment methods
library(multcompView)     # Visualizations of paired comparisons
library(lsmeans)          # least-squared means
library(tidyr)          	# data re-shaping
library(ggplot2)        	# plotting & data
library(dplyr)          	# data manipulation
library(magrittr)       	# pipe operator
library(gridExtra)     	  # provides side-by-side plotting
library(car)     		      # companion to applied regression
library(tidyverse)        # tidyverse
library(readxl)           # reads Excel files into R
library(devtools)
library(ggpubr)
library(agricolae)
#
#
# Make sure your grouping variables are classified as factors, not characters and not numbers
# This also orders your factors as they are in the datatable, otherwise R with alphabetize them. 
CHNTanninUpdate$Weevil = factor(CHNTanninUpdate$Weevil,
                           levels=unique(CHNTanninUpdate$Weevil))
CHNTanninUpdate$Sex = factor(CHNTanninUpdate$Sex,
                         levels=unique(CHNTanninUpdate$Sex))
CHNTanninUpdate$dummy = factor(CHNTanninUpdate$dummy,
                            levels=unique(CHNTanninUpdate$dummy))
#
# Set a random seed
set.seed(1431)
#
# Run a permutative 2-way ANOVA for %C (Yvar~Xvar1*Xvar2, data=YourDataSet, perm="Prob")
fit1 <- aovp(C~Weevil*Sex, data=CHNTanninUpdate, perm="Prob")
anova(fit1)
#
# Then if you want to get Tukey values, use this:
fit2 <- aovp(C ~ dummy, data = CHNTanninUpdate)
summary(fit2)

TukeyHSD(fit2)

HSD.test(fit2, "dummy", group=TRUE)
out1 <- HSD.test(fit2, "dummy", alpha = 0.05, group = TRUE, main = NULL, unbalanced = TRUE, console = TRUE)
#
# Repeat for %N
fit3 <- aovp(N~Weevil*Sex, data=CHNTanninUpdate, perm="Prob")
anova(fit3)
#
# Then if you want to get Tukey values, use this:
fit4 <- aovp(N ~ dummy, data = CHNTanninUpdate)
summary(fit4)

TukeyHSD(fit4)

HSD.test(fit4, "dummy", group=TRUE)
out1 <- HSD.test(fit4, "dummy", alpha = 0.05, group = TRUE, main = NULL, unbalanced = TRUE, console = TRUE)
#
#
# Repeat for C:N
fit5 <- aovp(CNmolar~Weevil*Sex, data=CHNTanninUpdate, perm="Prob")
anova(fit5)
#
# Then if you want to get Tukey values, use this:
fit6 <- aovp(CNmolar ~ dummy, data = CHNTanninUpdate)
summary(fit6)

TukeyHSD(fit6)

HSD.test(fit6, "dummy", group=TRUE)
out1 <- HSD.test(fit6, "dummy", alpha = 0.05, group = TRUE, main = NULL, unbalanced = TRUE, console = TRUE)
#
#
# Repeat for %CT
fit7 <- aovp(CTannin~Weevil*Sex, data=CHNTanninUpdate, perm="Prob")
anova(fit7)
#
# Then if you want to get Tukey values, use this:
fit8 <- aovp(CTannin ~ dummy, data = CHNTanninUpdate)
summary(fit8)

TukeyHSD(fit8)

HSD.test(fit8, "dummy", group=TRUE)
out1 <- HSD.test(fit8, "dummy", alpha = 0.05, group = TRUE, main = NULL, unbalanced = TRUE, console = TRUE)
#
#
#
#Now make a plot in grayscale - success!
#
#
#
plotN <- ggboxplot(CHNTanninUpdate, x = "Weevil", y = "N", fill = "Sex") +
  scale_fill_manual(values = c("gray50", "gray90")) +
  labs(tag = "a)", x = element_blank(), 
       y = "% N") +
  ylim(0, 3.1) +
  theme(legend.position = c(0.2, 0.9)) +
  theme(axis.text.x = element_blank()) +
  theme(text = element_text(size = 12))+
  theme(axis.line = element_line(size=0.5, colour = "black"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_blank(), panel.background = element_blank()) +
  theme(legend.title=element_blank()) +
  theme(legend.background = element_rect(fill="white", size=0.5)) +
  theme(legend.text = element_text(colour="black", size = 10)) +
  annotate(geom="text", size=2, x=0.8, y=0.9, label="a") +
  annotate(geom="text", size=2, x=1.2, y=1.0, label="a") +
  annotate(geom="text", size=2, x=1.8, y=3.1, label="b") +
  annotate(geom="text", size=2, x=2.2, y=2.8, label="b") +
  annotate(geom="text", size=2, x=1.6, y=0.8, label="F(3,42)= 178.28, p < 0.0001", hjust = 0) +
  annotate(geom="text", size=2, x=1.6, y=0.6, label="Sex: p = 0.8235", hjust = 0) +
  annotate(geom="text", size=2, x=1.6, y=0.4, label="Weevil: p < 0.0001", hjust = 0) +
  annotate(geom="text", size=2, x=1.6, y=0.2, label="Sex*Weevil: p = 0.1047", hjust = 0)
plotN
#
#
plotC <- ggboxplot(CHNTanninUpdate, x = "Weevil", y = "C", fill = "Sex") +
  scale_fill_manual(values = c("gray50", "gray90")) +
  labs(tag = "b)", x = element_blank(), y = "% C") +
  ylim(42, 49) +
  theme(legend.position ="none") +
  theme(axis.text.x = element_blank()) +
  theme(text = element_text(size = 12))+
  theme(axis.line = element_line(size=0.5, colour = "black"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_blank(), panel.background = element_blank()) +
  annotate(geom="text", size=2, x=0.6, y=44, label="F(3,42)= 1.82, p = 0.1582", hjust = 0) +
  annotate(geom="text", size=2, x=0.6, y=43.5, label="Sex: p = 0.5412", hjust = 0) +
  annotate(geom="text", size=2, x=0.6, y=43, label="Weevil: p = 0.1443", hjust = 0) +
  annotate(geom="text", size=2, x=0.6, y=42.5, label="Sex*Weevil: p = 0.9804", hjust = 0)
plotC
#
plotCN <- ggboxplot(CHNTanninUpdate, x = "Weevil", y = "CNmolar", fill = "Sex") +
  scale_fill_manual(values = c("gray50", "gray90")) +
  labs(tag = "c)", x = element_blank(), 
       y = "C : N") +
  ylim(0, 125) +
  theme(legend.position ="none") +
  theme(axis.text.x = element_blank()) +
  theme(text = element_text(size = 12))+
  theme(axis.line = element_line(size=0.5, colour = "black"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_blank(), panel.background = element_blank()) +
  theme(legend.title=element_blank()) +
  theme(legend.background = element_rect(fill="white", size=0.5)) +
  theme(legend.text = element_text(colour="black", size = 15)) +
  geom_text(size=2, x=1.8, y=45, label="a") +
  geom_text(size=2, x=2.2, y=47, label="a") +
  geom_text(size=2, x=0.8, y=118, label="b") +
  geom_text(size=2, x=1.2, y=124, label="b") +
  annotate(geom="text", size=2, x=0.6, y=40, label="F(3,42)= 406.7, p < 0.0001", hjust = 0) +
  annotate(geom="text", size=2, x=0.6, y=32, label="Sex: p = 0.3438", hjust = 0) +
  annotate(geom="text", size=2, x=0.6, y=24, label="Weevil: p < 0.0001", hjust = 0) +
  annotate(geom="text", size=2, x=0.6, y=16, label="Sex*Weevil: p = 0.0993", hjust = 0)
plotCN
#
#
plotCT <- ggboxplot(CHNTanninUpdate, x = "Weevil", y = "CTannin", fill = "Sex") +
  scale_fill_manual(values = c("gray50", "gray90")) +
  labs(tag = "d)", x = "Weevil", 
       y = "% CT") +
  ylim(0, 33) +
  theme(legend.position ="none") +
  theme(axis.text.x = element_text(angle=0)) +
  theme(text = element_text(size = 12))+
  theme(axis.line = element_line(size=0.5, colour = "black"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_blank(), panel.background = element_blank()) +
  theme(legend.title=element_blank()) +
  theme(legend.background = element_rect(fill="white", size=0.5)) +
  theme(legend.text = element_text(colour="black", size = 15)) +
  geom_text(size=2, x=0.8, y=32, label="b") +
  geom_text(size=2, x=1.2, y=31, label="b") +
  geom_text(size=2, x=1.8, y=23, label="a") +
  geom_text(size=2, x=2.2, y=16, label="a") +
  annotate(geom="text", size=2, x=0.6, y=10, label="F(3,42)= 93.3, p < 0.0001", hjust = 0) +
  annotate(geom="text", size=2, x=0.6, y=7.5, label="Sex: p = 0.2934", hjust = 0) +
  annotate(geom="text", size=2, x=0.6, y=5, label="Weevil: p < 0.0001", hjust = 0) +
  annotate(geom="text", size=2, x=0.6, y=2.5, label="Sex*Weevil: p = 0.3514", hjust = 0)
plotCT
# 
#
# Arrange plots A and B in one column:
gN <- ggplotGrob(plotN)
gC <- ggplotGrob(plotC)
gCN <- ggplotGrob(plotCN)
gCT <- ggplotGrob(plotCT)
g <- arrangeGrob(gN, gC, gCN, gCT, nrow=4)
# Warning message: Removed 7 rows containing non-finite values (stat_boxplot). This is okay - 7 missing (NA) values for CT
#
ggsave(file="Fig2-grays.pdf", g, width=8.4, height=20, units="cm", dpi=800)

