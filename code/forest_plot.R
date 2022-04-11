library(ggplot2)
library(data.table)
library(here)

#setwd("~/farhan/data")

%args = commandArgs(trailingOnly=TRUE)


# test if there is at least one argument: if not, return an error
#if (length(args)==0) {
#name = c("RV_noFilter", "RV_proteinDomain")
#} else {
#  name = args
#}



filename = here('farhan', 'data', 'proteinDomain', 'all_lr.csv')
dt <- read.table(filename, sep=",", header=TRUE)

# --------------------------------------------------------------
#exome vs protein functional domain

name = c("RV_noFilter", "RV_proteinDomain")
title = "Enrichment in Exome and Protein Functional Domain"
num <- length(name)

subset <- dt[dt$geneset %in% name,]
subset$geneset <- c("Exome", "Exome", "Exome", "Exome", "Exome", "Exome",
                    "Protein Functional Domain","Protein Functional Domain","Protein Functional Domain","Protein Functional Domain","Protein Functional Domain","Protein Functional Domain")
subset <- subset[c(6,1,2,3,4,5,12,7,8,9,10,11),]
subset$terms <- c("Synonymous", "Benign missense", "Missense", "Damaging Missense", "Missense MPC >= 2", "PTV",
                  "Synonymous", "Benign missense", "Missense", "Damaging Missense", "Missense MPC >= 2", "PTV")

ggplot(data = subset, aes(x=terms, y=coef, ymin=ci.lower, ymax=ci.upper)) +
  geom_pointrange(aes(col=terms)) +
  geom_hline(aes(fill=terms), yintercept = 1, linetype = 2) +
  #rename the axes
  xlab('Logistic regression') + ylab('Odds ratio and 95% CI') +
  geom_errorbar(aes(ymin=ci.lower, ymax=ci.upper, col=terms), width=0.5, cex=1) +
  #facet wrap the forest plots for the different models together, giving each its own row
  facet_wrap(~geneset, strip.position = "left", nrow=num, scales = "free_y", labeller = label_wrap_gen(width = 20, multi_line = TRUE)) +
  #flip the forest plot to be horizontal rather than vertical 
  coord_flip(clip = "off")+theme(axis.text.y = element_blank(), axis.ticks.y=element_blank()) +
  theme_classic() +
  theme(legend.position="bottom", legend.title=element_text(size=14, face = "bold"), legend.text=element_text(size=12), axis.title = element_text(size = 14, face = "bold"), 
        strip.text.x = element_blank(), strip.text= element_text(size=12), axis.text.y = element_blank(), axis.ticks.y = element_blank(), 
        axis.text.x = element_text(size=12, colour = "black"), 
        plot.margin = margin(1,3,1,1, "lines"), plot.title = element_text(size=18)) +
  #add the p-values as labels
  geom_text(aes(y=1.12, label=formatC(pvalue, format = "e", digits = 2)), size=4, hjust=0)+
  theme(panel.spacing.x = unit(5, "lines"))+
  ggtitle(paste0(title, collapse = ", "))+
  scale_color_manual(name="Variant class", breaks = c("Synonymous", "Benign missense", "Missense", "Damaging missense", "Missense MPC >= 2", "PTV"), 
                     values = c("#89B374", "#3c92e0", "#674EA7", "#f0c542", "#e88e2f", "#f0452e"))


# ---------------------------------------------------------

name = c("RV_allHigh", "RV_allMedium", "RV_allLow")
title = "Enrichment in Expression Levels within Exome"
num <- length(name)

subset <- dt[dt$geneset %in% name,]
subset$geneset <- c("High", "High", "High", "High", "High", "High", 
                    "Medium", "Medium","Medium","Medium","Medium","Medium",
                    "Low","Low","Low","Low","Low","Low")
subset <- subset[c(6,1,2,3,4,5,12,7,8,9,10,11,18,13,14,15,16,17),]
subset$terms <- c("Synonymous", "Benign missense", "Missense", "Damaging Missense", "Missense MPC >= 2", "PTV",
                  "Synonymous", "Benign missense", "Missense", "Damaging Missense", "Missense MPC >= 2", "PTV",
                  "Synonymous", "Benign missense", "Missense", "Damaging Missense", "Missense MPC >= 2", "PTV")

ggplot(data = subset, aes(x=terms, y=coef, ymin=ci.lower, ymax=ci.upper)) +
  geom_pointrange(aes(col=terms)) +
  geom_hline(aes(fill=terms), yintercept = 1, linetype = 2) +
  #rename the axes
  xlab('Logistic regression') + ylab('Odds ratio and 95% CI') +
  geom_errorbar(aes(ymin=ci.lower, ymax=ci.upper, col=terms), width=0.5, cex=1) +
  #facet wrap the forest plots for the different models together, giving each its own row
  facet_wrap(~geneset, strip.position = "left", nrow=num, scales = "free_y", labeller = label_wrap_gen(width = 20, multi_line = TRUE)) +
  #flip the forest plot to be horizontal rather than vertical 
  coord_flip(clip = "off")+theme(axis.text.y = element_blank(), axis.ticks.y=element_blank()) +
  theme_classic() +
  theme(legend.position="bottom", legend.title=element_text(size=14, face = "bold"), legend.text=element_text(size=12), axis.title = element_text(size = 14, face = "bold"), 
        strip.text.x = element_blank(), strip.text= element_text(size=12), axis.text.y = element_blank(), axis.ticks.y = element_blank(), 
        axis.text.x = element_text(size=12, colour = "black"), 
        plot.margin = margin(1,3,1,1, "lines"), plot.title = element_text(size=18)) +
  #add the p-values as labels
  geom_text(aes(y=1.12, label=formatC(pvalue, format = "e", digits = 2)), size=4, hjust=-1)+
  theme(panel.spacing.x = unit(5, "lines"))+
  ggtitle(paste0(title, collapse = ", "))+
  scale_color_manual(name="Variant class", breaks = c("Synonymous", "Benign missense", "Missense", "Damaging missense", "Missense MPC >= 2", "PTV"), 
                     values = c("#89B374", "#3c92e0", "#674EA7", "#f0c542", "#e88e2f", "#f0452e"))


# ----------------------------------------------------------
# plotfor brain expression levels in protein domain
name = c("RV_brainHigh", "RV_brainMedium", "RV_brainLow")
title = "Enrichment in Brain Expression Levels within Exome"
num <- length(name)

subset <- dt[dt$geneset %in% name,]
subset$geneset <- c("High", "High", "High", "High", "High", "High", 
                    "Medium", "Medium","Medium","Medium","Medium","Medium",
                    "Low","Low","Low","Low","Low","Low")
subset <- subset[c(6,1,2,3,4,5,12,7,8,9,10,11,18,13,14,15,16,17),]
subset$terms <- c("Synonymous", "Benign missense", "Missense", "Damaging Missense", "Missense MPC >= 2", "PTV",
                  "Synonymous", "Benign missense", "Missense", "Damaging Missense", "Missense MPC >= 2", "PTV",
                  "Synonymous", "Benign missense", "Missense", "Damaging Missense", "Missense MPC >= 2", "PTV")

ggplot(data = subset, aes(x=terms, y=coef, ymin=ci.lower, ymax=ci.upper)) +
  geom_pointrange(aes(col=terms)) +
  geom_hline(aes(fill=terms), yintercept = 1, linetype = 2) +
  #rename the axes
  xlab('Logistic regression') + ylab('Odds ratio and 95% CI') +
  geom_errorbar(aes(ymin=ci.lower, ymax=ci.upper, col=terms), width=0.5, cex=1) +
  #facet wrap the forest plots for the different models together, giving each its own row
  facet_wrap(~geneset, strip.position = "left", nrow=num, scales = "free_y", labeller = label_wrap_gen(width = 20, multi_line = TRUE)) +
  #flip the forest plot to be horizontal rather than vertical 
  coord_flip(clip = "off")+theme(axis.text.y = element_blank(), axis.ticks.y=element_blank()) +
  theme_classic() +
  theme(legend.position="bottom", legend.title=element_text(size=14, face = "bold"), legend.text=element_text(size=12), axis.title = element_text(size = 14, face = "bold"), 
        strip.text.x = element_blank(), strip.text= element_text(size=12), axis.text.y = element_blank(), axis.ticks.y = element_blank(), 
        axis.text.x = element_text(size=12, colour = "black"), 
        plot.margin = margin(1,3,1,1, "lines"), plot.title = element_text(size=18)) +
  #add the p-values as labels
  geom_text(aes(y=1.12, label=formatC(pvalue, format = "e", digits = 2)), size=4, hjust=-1)+
  theme(panel.spacing.x = unit(5, "lines"))+
  ggtitle(paste0(title, collapse = ", "))+
  scale_color_manual(name="Variant class", breaks = c("Synonymous", "Benign missense", "Missense", "Damaging missense", "Missense MPC >= 2", "PTV"), 
                     values = c("#89B374", "#3c92e0", "#674EA7", "#f0c542", "#e88e2f", "#f0452e"))

# -------------------------------------------------------------
# all expression levels in protein domain

name = c("RV_pd_allHigh", "RV_pd_allMedium", "RV_pd_allLow")
title = "Enrichment in Expression Levels within Protein Functional Domain"
num <- length(name)

subset <- dt[dt$geneset %in% name,]
subset$geneset <- c("High", "High", "High", "High", "High", "High", 
                    "Medium", "Medium","Medium","Medium","Medium","Medium",
                    "Low","Low","Low","Low","Low","Low")
subset <- subset[c(6,1,2,3,4,5,12,7,8,9,10,11,18,13,14,15,16,17),]
subset$terms <- c("Synonymous", "Benign missense", "Missense", "Damaging Missense", "Missense MPC >= 2", "PTV",
                  "Synonymous", "Benign missense", "Missense", "Damaging Missense", "Missense MPC >= 2", "PTV",
                  "Synonymous", "Benign missense", "Missense", "Damaging Missense", "Missense MPC >= 2", "PTV")

ggplot(data = subset, aes(x=terms, y=coef, ymin=ci.lower, ymax=ci.upper)) +
  geom_pointrange(aes(col=terms)) +
  geom_hline(aes(fill=terms), yintercept = 1, linetype = 2) +
  #rename the axes
  xlab('Logistic regression') + ylab('Odds ratio and 95% CI') +
  geom_errorbar(aes(ymin=ci.lower, ymax=ci.upper, col=terms), width=0.5, cex=1) +
  #facet wrap the forest plots for the different models together, giving each its own row
  facet_wrap(~geneset, strip.position = "left", nrow=num, scales = "free_y", labeller = label_wrap_gen(width = 20, multi_line = TRUE)) +
  #flip the forest plot to be horizontal rather than vertical 
  coord_flip(clip = "off")+theme(axis.text.y = element_blank(), axis.ticks.y=element_blank()) +
  theme_classic() +
  theme(legend.position="bottom", legend.title=element_text(size=14, face = "bold"), legend.text=element_text(size=12), axis.title = element_text(size = 14, face = "bold"), 
        strip.text.x = element_blank(), strip.text= element_text(size=12), axis.text.y = element_blank(), axis.ticks.y = element_blank(), 
        axis.text.x = element_text(size=12, colour = "black"), 
        plot.margin = margin(1,3,1,1, "lines"), plot.title = element_text(size=18)) +
  #add the p-values as labels
  geom_text(aes(y=1.12, label=formatC(pvalue, format = "e", digits = 2)), size=4, hjust=-6)+
  theme(panel.spacing.x = unit(5, "lines"))+
  ggtitle(paste0(title, collapse = ", "))+
  scale_color_manual(name="Variant class", breaks = c("Synonymous", "Benign missense", "Missense", "Damaging missense", "Missense MPC >= 2", "PTV"), 
                     values = c("#89B374", "#3c92e0", "#674EA7", "#f0c542", "#e88e2f", "#f0452e"))


# ----------------------------------------------------------
# plotfor brain expression levels in protein domain
name = c("RV_pd_brainHigh", "RV_pd_brainMedium", "RV_pd_brainLow")
title = "Enrichment in Brain Expression Levels within Functional Protein Domain"
num <- length(name)

subset <- dt[dt$geneset %in% name,]
subset$geneset <- c("High", "High", "High", "High", "High", "High", 
                    "Medium", "Medium","Medium","Medium","Medium","Medium",
                    "Low","Low","Low","Low","Low","Low")
subset <- subset[c(6,1,2,3,4,5,12,7,8,9,10,11,18,13,14,15,16,17),]
subset$terms <- c("Synonymous", "Benign missense", "Missense", "Damaging Missense", "Missense MPC >= 2", "PTV",
                  "Synonymous", "Benign missense", "Missense", "Damaging Missense", "Missense MPC >= 2", "PTV",
                  "Synonymous", "Benign missense", "Missense", "Damaging Missense", "Missense MPC >= 2", "PTV")

ggplot(data = subset, aes(x=terms, y=coef, ymin=ci.lower, ymax=ci.upper)) +
  geom_pointrange(aes(col=terms)) +
  geom_hline(aes(fill=terms), yintercept = 1, linetype = 2) +
  #rename the axes
  xlab('Logistic regression') + ylab('Odds ratio and 95% CI') +
  geom_errorbar(aes(ymin=ci.lower, ymax=ci.upper, col=terms), width=0.5, cex=1) +
  #facet wrap the forest plots for the different models together, giving each its own row
  facet_wrap(~geneset, strip.position = "left", nrow=num, scales = "free_y", labeller = label_wrap_gen(width = 20, multi_line = TRUE)) +
  #flip the forest plot to be horizontal rather than vertical 
  coord_flip(clip = "off")+theme(axis.text.y = element_blank(), axis.ticks.y=element_blank()) +
  theme_classic() +
  theme(legend.position="bottom", legend.title=element_text(size=14, face = "bold"), legend.text=element_text(size=12), axis.title = element_text(size = 14, face = "bold"), 
        strip.text.x = element_blank(), strip.text= element_text(size=12), axis.text.y = element_blank(), axis.ticks.y = element_blank(), 
        axis.text.x = element_text(size=12, colour = "black"), 
        plot.margin = margin(1,3,1,1, "lines"), plot.title = element_text(size=18)) +
  #add the p-values as labels
  geom_text(aes(y=1.12, label=formatC(pvalue, format = "e", digits = 2)), size=4, hjust=-6)+
  theme(panel.spacing.x = unit(5, "lines"))+
  ggtitle(paste0(title, collapse = ", "))+
  scale_color_manual(name="Variant class", breaks = c("Synonymous", "Benign missense", "Missense", "Damaging missense", "Missense MPC >= 2", "PTV"), 
                     values = c("#89B374", "#3c92e0", "#674EA7", "#f0c542", "#e88e2f", "#f0452e"))