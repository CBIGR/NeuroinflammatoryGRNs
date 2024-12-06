##### Barplots for presecision and recall

# http://www.sthda.com/english/wiki/ggplot2-barplots-quick-start-guide-r-software-and-data-visualization 

library(tidyverse)
setwd("C:/Users/hannepu/OneDrive - UGent/Master 2/Master's dissertation/Analysis/")

KnockTF <- read.table("KnockTF/PR_all_KnockTF.txt")
OmniPath <- read.table("OmniPath/PR_Omnipath.txt")
Marbach <- read.table("OmniPath/PR_Marbach.txt")
UniBind <- read.table("UniBind/PR_all_UniBind.txt")

# Make dataframes for each plot
AD_bulk <- data.frame("Ground_truth" = c(OmniPath[1,], OmniPath[5,], Marbach[1,], Marbach[5,], KnockTF[1,], KnockTF[5,], UniBind[1,], UniBind[5,]),
                      "Source" = c("OmniPath", "OmniPath", "Marbach", "Marbach", "KnockTF", "KnockTF", "UniBind", "UniBind"), 
                      "Measure" = c("Precision", "Recall", "Precision", "Recall", "Precision", "Recall", "Precision", "Recall"))

MDD_bulk <- data.frame("Ground_truth" = c(OmniPath[2,], OmniPath[6,], Marbach[2,], Marbach[6,], KnockTF[2,], KnockTF[6,], UniBind[2,], UniBind[6,]),
                       "Source" = c("OmniPath", "OmniPath", "Marbach", "Marbach", "KnockTF", "KnockTF", "UniBind", "UniBind"), 
                       "Measure" = c("Precision", "Recall", "Precision", "Recall", "Precision", "Recall", "Precision", "Recall"))


AD_sc <- data.frame("Ground_truth" = c(OmniPath[3,], OmniPath[7,], Marbach[3,], Marbach[7,], KnockTF[3,], KnockTF[7,], UniBind[3,], UniBind[7,]),
                    "Source" = c("OmniPath", "OmniPath", "Marbach", "Marbach", "KnockTF", "KnockTF", "UniBind", "UniBind"), 
                    "Measure" = c("Precision", "Recall", "Precision", "Recall", "Precision", "Recall", "Precision", "Recall"))

MDD_sc <- data.frame("Ground_truth" = c(OmniPath[4,], OmniPath[8,], Marbach[3,], Marbach[7,], KnockTF[4,], KnockTF[8,], UniBind[4,], UniBind[8,]),
                     "Source" = c("OmniPath", "OmniPath", "Marbach", "Marbach", "KnockTF", "KnockTF", "UniBind", "UniBind"), 
                     "Measure" = c("Precision", "Recall", "Precision", "Recall", "Precision", "Recall", "Precision", "Recall"))

AD_bulk$Source <- as.factor(AD_bulk$Source)
AD_bulk$Measure <- as.factor(AD_bulk$Measure)
MDD_bulk$Source <- as.factor(MDD_bulk$Source)
MDD_bulk$Measure <- as.factor(MDD_bulk$Measure)
AD_sc$Source <- as.factor(AD_sc$Source)
AD_sc$Measure <- as.factor(AD_sc$Measure)
MDD_sc$Source <- as.factor(MDD_sc$Source)
MDD_sc$Measure <- as.factor(MDD_sc$Measure)

#### Plots ####
ggplot(AD_bulk, aes(x = Source, y = Ground_truth, fill = Measure)) +
  geom_bar(stat="identity", position = position_dodge()) +
  scale_fill_manual(values = c('#99DDFF', '#44BB99')) +
  theme_classic() +
  labs(title = "AD consensus network", x = "Source", y = "Value") +
  theme(text = element_text(size = 18)) +
  geom_text(aes(label = round(Ground_truth, 3)), vjust = -0.3, color = "black",
            position = position_dodge(0.9), size = 3.5)

ggplot(MDD_bulk, aes(x = Source, y = Ground_truth, fill = Measure)) +
  geom_bar(stat="identity", position= position_dodge()) +
  scale_fill_manual(values=c('#99DDFF', '#44BB99')) +
  theme_classic() +
  labs(title="MDD concensus network", x = "Source", y = "Value") +
  theme(text = element_text(size = 18)) +
  geom_text(aes(label = round(Ground_truth, 3)), vjust = -0.3, color = "black",
            position = position_dodge(0.9), size = 3.5)


ggplot(AD_sc, aes(x = Source, y = Ground_truth, fill = Measure)) +
  geom_bar(stat="identity", position= position_dodge()) +
  scale_fill_manual(values=c('#99DDFF', '#44BB99')) +
  theme_classic() +
  labs(title="AD single-cell network", x = "Source", y = "Value") +
  theme(text = element_text(size = 18)) +
  geom_text(aes(label = round(Ground_truth, 3)), vjust = -0.3, color = "black",
            position = position_dodge(0.9), size = 3.5)

ggplot(MDD_sc, aes(x = Source, y = Ground_truth, fill = Measure)) +
  geom_bar(stat="identity", position= position_dodge()) +
  scale_fill_manual(values=c('#99DDFF', '#44BB99')) +
  theme_classic() +
  labs(title="MDD single-cell network", x = "Source", y = "Value") +
  theme(text = element_text(size = 18)) +
  geom_text(aes(label = round(Ground_truth, 3)), vjust = -0.3, color = "black",
            position = position_dodge(0.9), size = 3.5)
