library(ggplot2)
library(tidyr)
library(dplyr)
library(cowplot)

## FIRST TEST: No of ORFs

tablaoriginal <- read.csv("C:\\Users/SES089/OneDrive - University of Salford/Documents/Project_LES/VIRUS_RefSeq/viral_genomes.ORIGINAL.tsv", header = T, sep = "\t")
tablaprodigal <- read.csv("C:\\Users/SES089/OneDrive - University of Salford/Documents/Project_LES/VIRUS_RefSeq/viral_genomes.PRODIGAL.tsv", header = T, sep = "\t")
tablametaprodigal <- read.csv("C:\\Users/SES089/OneDrive - University of Salford/Documents/Project_LES/VIRUS_RefSeq/viral_genomes.METAPRODIGAL.tsv", header = T, sep = "\t")
tablaglimmer <- read.csv("C:\\Users/SES089/OneDrive - University of Salford/Documents/Project_LES/VIRUS_RefSeq/viral_genomes.GLIMMER.tsv", header = T, sep = "\t")
tablagenemarks <- read.csv("C:\\Users/SES089/OneDrive - University of Salford/Documents/Project_LES/VIRUS_RefSeq/viral_genomes.GENEMARKS.tsv", header = T, sep = "\t")
tablaphanotate <- read.csv("C:\\Users/SES089/OneDrive - University of Salford/Documents/Project_LES/VIRUS_RefSeq/viral_genomes.PHANOTATE.tsv", header = T, sep = "\t")
tablafraggenescan <- read.csv("C:\\Users/SES089/OneDrive - University of Salford/Documents/Project_LES/VIRUS_RefSeq/viral_genomes.FRAGGENESCAN.tsv", header = T, sep = "\t")
tablamga <- read.csv("C:\\Users/SES089/OneDrive - University of Salford/Documents/Project_LES/VIRUS_RefSeq/viral_genomes.MGA.tsv", header = T, sep = "\t")
tablaaugustussa <- read.csv("C:\\Users/SES089/OneDrive - University of Salford/Documents/Project_LES/VIRUS_RefSeq/viral_genomes.AUGUSTUS_SAUREUS.tsv", header = T, sep = "\t")
tablaaugustusec <- read.csv("C:\\Users/SES089/OneDrive - University of Salford/Documents/Project_LES/VIRUS_RefSeq/viral_genomes.AUGUSTUS_ECOLI.tsv", header = T, sep = "\t")
tablaaugustushs <- read.csv("C:\\Users/SES089/OneDrive - University of Salford/Documents/Project_LES/VIRUS_RefSeq/viral_genomes.AUGUSTUS_HUMAN.tsv", header = T, sep = "\t")
tablahost <- read.csv("C:\\Users/SES089/OneDrive - University of Salford/Documents/Project_LES/VIRUS_RefSeq/list_organisms.tsv", header = T, sep = "\t")
tablanucleic <- read.csv("C:\\Users/SES089/OneDrive - University of Salford/Documents/Project_LES/VIRUS_RefSeq/list_viralnucleic_2.csv", header = T, sep = "\t")

tabula <- tablaoriginal %>% 
  full_join(tablahost, by = c("ID" = "ID")) %>%
  full_join(tablanucleic, by = c("ID" = "ID")) %>%
  full_join(tablaprodigal, by = c("ID" = "ID")) %>% 
  full_join(tablametaprodigal, by = c("ID" = "ID")) %>% 
  full_join(tablaglimmer, by = c("ID" = "ID")) %>% 
  full_join(tablagenemarks, by = c("ID" = "ID")) %>% 
  full_join(tablaphanotate, by = c("ID" = "ID")) %>% 
  full_join(tablafraggenescan, by = c("ID" = "ID")) %>% 
  full_join(tablamga, by = c("ID" = "ID")) %>% 
  full_join(tablaaugustussa, by = c("ID" = "ID")) %>% 
  full_join(tablaaugustusec, by = c("ID" = "ID")) %>% 
  full_join(tablaaugustushs, by = c("ID" = "ID"))

prodigallm <- lm(CDSs_PRODIGAL ~ CDSs_EXPECTED + 0, data = tabula)
metaprodigallm <- lm(CDSs_METAPRODIGAL ~ CDSs_EXPECTED + 0, data = tabula)
glimmerlm <- lm(CDSs_GLIMMER ~ CDSs_EXPECTED + 0, data = tabula)
genemarkslm <- lm(CDSs_GENEMARKS ~ CDSs_EXPECTED + 0, data = tabula)
phanotatelm <- lm(CDSs_PHANOTATE ~ CDSs_EXPECTED + 0, data = tabula)
fraggenescanlm <- lm(CDSs_FRAGGENESCAN ~ CDSs_EXPECTED + 0, data = tabula)
mgalm <- lm(CDSs_MGA ~ CDSs_EXPECTED + 0, data = tabula)
augustuslm <- lm(CDSs_AUGUSTUS_SAUREUS ~ CDSs_EXPECTED + 0, data = tabula)
augustus2lm <- lm(CDSs_AUGUSTUS_ECOLI ~ CDSs_EXPECTED + 0, data = tabula)
augustus3lm <- lm(CDSs_AUGUSTUS_HUMAN ~ CDSs_EXPECTED + 0, data = tabula)

summary(prodigallm)
summary(metaprodigallm)
summary(glimmerlm)
summary(genemarkslm)
summary(phanotatelm)
summary(fraggenescanlm)
summary(mgalm)
summary(augustuslm)
summary(augustus2lm)
summary(augustus3lm)

tabula.long <- tabula %>%
  pivot_longer(cols=c(2,5:14), names_to = "Program", values_to = "No_CDS")

tabula.long.2 <- tabula %>%
  pivot_longer(cols=c(5:14), names_to = "Program", values_to = "No_CDS")
generallm <- lm(CDSs_EXPECTED ~ No_CDS * Program + 0, data = tabula.long.2)
anova(generallm)
summary(generallm)

# Eukaryotic viruses

tabula_euk <- tabula[grepl("Eukaryote", tabula$Host_domain),]
prodigallm_euk <- lm(CDSs_PRODIGAL ~ CDSs_EXPECTED + 0, data = tabula_euk)
metaprodigallm_euk <- lm(CDSs_METAPRODIGAL ~ CDSs_EXPECTED + 0, data = tabula_euk)
glimmerlm_euk <- lm(CDSs_GLIMMER ~ CDSs_EXPECTED + 0, data = tabula_euk)
genemarkslm_euk <- lm(CDSs_GENEMARKS ~ CDSs_EXPECTED + 0, data = tabula_euk)
phanotatelm_euk <- lm(CDSs_PHANOTATE ~ CDSs_EXPECTED + 0, data = tabula_euk)
fraggenescanlm_euk <- lm(CDSs_FRAGGENESCAN ~ CDSs_EXPECTED + 0, data = tabula_euk)
mgalm_euk <- lm(CDSs_MGA ~ CDSs_EXPECTED + 0, data = tabula_euk)
augustuslm_euk <- lm(CDSs_AUGUSTUS_SAUREUS ~ CDSs_EXPECTED + 0, data = tabula_euk)
augustus2lm_euk <- lm(CDSs_AUGUSTUS_ECOLI ~ CDSs_EXPECTED + 0, data = tabula_euk)
augustus3lm_euk <- lm(CDSs_AUGUSTUS_HUMAN ~ CDSs_EXPECTED + 0, data = tabula_euk)
summary(prodigallm_euk)
summary(metaprodigallm_euk)
summary(glimmerlm_euk)
summary(genemarkslm_euk)
summary(phanotatelm_euk)
summary(fraggenescanlm_euk)
summary(mgalm_euk)
summary(augustuslm_euk)
summary(augustus2lm_euk)
summary(augustus3lm_euk)

tabula_euk.long.2 <- tabula %>%
  pivot_longer(cols=c(5:14), names_to = "Program", values_to = "No_CDS")
generallm_euk <- lm(CDSs_EXPECTED ~ No_CDS * Program + 0, data = tabula_euk.long.2)
anova(generallm_euk)
summary(generallm_euk)

# Archaeal viruses

tabula_arc <- tabula[grepl("Archaea", tabula$Host_domain),]
prodigallm_arc <- lm(CDSs_PRODIGAL ~ CDSs_EXPECTED + 0, data = tabula_arc)
metaprodigallm_arc <- lm(CDSs_METAPRODIGAL ~ CDSs_EXPECTED + 0, data = tabula_arc)
glimmerlm_arc <- lm(CDSs_GLIMMER ~ CDSs_EXPECTED + 0, data = tabula_arc)
genemarkslm_arc <- lm(CDSs_GENEMARKS ~ CDSs_EXPECTED + 0, data = tabula_arc)
phanotatelm_arc <- lm(CDSs_PHANOTATE ~ CDSs_EXPECTED + 0, data = tabula_arc)
fraggenescanlm_arc <- lm(CDSs_FRAGGENESCAN ~ CDSs_EXPECTED + 0, data = tabula_arc)
mgalm_arc <- lm(CDSs_MGA ~ CDSs_EXPECTED + 0, data = tabula_arc)
augustuslm_arc <- lm(CDSs_AUGUSTUS_SAUREUS ~ CDSs_EXPECTED + 0, data = tabula_arc)
augustus2lm_arc <- lm(CDSs_AUGUSTUS_ECOLI ~ CDSs_EXPECTED + 0, data = tabula_arc)
augustus3lm_arc <- lm(CDSs_AUGUSTUS_HUMAN ~ CDSs_EXPECTED + 0, data = tabula_arc)
summary(prodigallm_arc)
summary(metaprodigallm_arc)
summary(glimmerlm_arc)
summary(genemarkslm_arc)
summary(phanotatelm_arc)
summary(fraggenescanlm_arc)
summary(mgalm_arc)
summary(augustuslm_arc)
summary(augustus2lm_arc)
summary(augustus3lm_arc)

# Bacteriophages

tabula_bac <- tabula[grepl("Bacteria", tabula$Host_domain),]
prodigallm_bac <- lm(CDSs_PRODIGAL ~ CDSs_EXPECTED + 0, data = tabula_bac)
metaprodigallm_bac <- lm(CDSs_METAPRODIGAL ~ CDSs_EXPECTED + 0, data = tabula_bac)
glimmerlm_bac <- lm(CDSs_GLIMMER ~ CDSs_EXPECTED + 0, data = tabula_bac)
genemarkslm_bac <- lm(CDSs_GENEMARKS ~ CDSs_EXPECTED + 0, data = tabula_bac)
phanotatelm_bac <- lm(CDSs_PHANOTATE ~ CDSs_EXPECTED + 0, data = tabula_bac)
fraggenescanlm_bac <- lm(CDSs_FRAGGENESCAN ~ CDSs_EXPECTED + 0, data = tabula_bac)
mgalm_bac <- lm(CDSs_MGA ~ CDSs_EXPECTED + 0, data = tabula_bac)
augustuslm_bac <- lm(CDSs_AUGUSTUS_SAUREUS ~ CDSs_EXPECTED + 0, data = tabula_bac)
augustus2lm_bac <- lm(CDSs_AUGUSTUS_ECOLI ~ CDSs_EXPECTED + 0, data = tabula_bac)
augustus3lm_bac <- lm(CDSs_AUGUSTUS_HUMAN ~ CDSs_EXPECTED + 0, data = tabula_bac)
summary(prodigallm_bac)
summary(metaprodigallm_bac)
summary(glimmerlm_bac)
summary(genemarkslm_bac)
summary(phanotatelm_bac)
summary(fraggenescanlm_bac)
summary(mgalm_bac)
summary(augustuslm_bac)
summary(augustus2lm_bac)
summary(augustus3lm_bac)

tabula_bac.long.2 <- tabula %>%
  pivot_longer(cols=c(5:14), names_to = "Program", values_to = "No_CDS")
generallm_bac <- lm(CDSs_EXPECTED ~ No_CDS * Program + 0, data = tabula_bac.long.2)
anova(generallm_bac)
summary(generallm_bac)

# DNA viruses

tabula_dna <- tabula[grepl("DNA", tabula$Nucleic_acid),]
prodigallm_dna <- lm(CDSs_PRODIGAL ~ CDSs_EXPECTED + 0, data = tabula_dna)
metaprodigallm_dna <- lm(CDSs_METAPRODIGAL ~ CDSs_EXPECTED + 0, data = tabula_dna)
glimmerlm_dna <- lm(CDSs_GLIMMER ~ CDSs_EXPECTED + 0, data = tabula_dna)
genemarkslm_dna <- lm(CDSs_GENEMARKS ~ CDSs_EXPECTED + 0, data = tabula_dna)
phanotatelm_dna <- lm(CDSs_PHANOTATE ~ CDSs_EXPECTED + 0, data = tabula_dna)
fraggenescanlm_dna <- lm(CDSs_FRAGGENESCAN ~ CDSs_EXPECTED + 0, data = tabula_dna)
mgalm_dna <- lm(CDSs_MGA ~ CDSs_EXPECTED + 0, data = tabula_dna)
augustuslm_dna <- lm(CDSs_AUGUSTUS_SAUREUS ~ CDSs_EXPECTED + 0, data = tabula_dna)
augustus2lm_dna <- lm(CDSs_AUGUSTUS_ECOLI ~ CDSs_EXPECTED + 0, data = tabula_dna)
augustus3lm_dna <- lm(CDSs_AUGUSTUS_HUMAN ~ CDSs_EXPECTED + 0, data = tabula_dna)
summary(prodigallm_dna)
summary(metaprodigallm_dna)
summary(glimmerlm_dna)
summary(genemarkslm_dna)
summary(phanotatelm_dna)
summary(fraggenescanlm_dna)
summary(mgalm_dna)
summary(augustuslm_dna)
summary(augustus2lm_dna)
summary(augustus3lm_dna)

tabula_dna.long.2 <- tabula %>%
  pivot_longer(cols=c(5:14), names_to = "Program", values_to = "No_CDS")
generallm_dna <- lm(CDSs_EXPECTED ~ No_CDS * Program + 0, data = tabula_dna.long.2)
anova(generallm_dna)
summary(generallm_dna)

# RNA viruses

tabula_rna <- tabula[grepl("RNA", tabula$Nucleic_acid),]
prodigallm_rna <- lm(CDSs_PRODIGAL ~ CDSs_EXPECTED + 0, data = tabula_rna)
metaprodigallm_rna <- lm(CDSs_METAPRODIGAL ~ CDSs_EXPECTED + 0, data = tabula_rna)
glimmerlm_rna <- lm(CDSs_GLIMMER ~ CDSs_EXPECTED + 0, data = tabula_rna)
genemarkslm_rna <- lm(CDSs_GENEMARKS ~ CDSs_EXPECTED + 0, data = tabula_rna)
phanotatelm_rna <- lm(CDSs_PHANOTATE ~ CDSs_EXPECTED + 0, data = tabula_rna)
fraggenescanlm_rna <- lm(CDSs_FRAGGENESCAN ~ CDSs_EXPECTED + 0, data = tabula_rna)
mgalm_rna <- lm(CDSs_MGA ~ CDSs_EXPECTED + 0, data = tabula_rna)
augustuslm_rna <- lm(CDSs_AUGUSTUS_SAUREUS ~ CDSs_EXPECTED + 0, data = tabula_rna)
augustus2lm_rna <- lm(CDSs_AUGUSTUS_ECOLI ~ CDSs_EXPECTED + 0, data = tabula_rna)
augustus3lm_rna <- lm(CDSs_AUGUSTUS_HUMAN ~ CDSs_EXPECTED + 0, data = tabula_rna)
summary(prodigallm_rna)
summary(metaprodigallm_rna)
summary(glimmerlm_rna)
summary(genemarkslm_rna)
summary(phanotatelm_rna)
summary(fraggenescanlm_rna)
summary(mgalm_rna)
summary(augustuslm_rna)
summary(augustus2lm_rna)
summary(augustus3lm_rna)

tabula_rna.long.2 <- tabula %>%
  pivot_longer(cols=c(5:14), names_to = "Program", values_to = "No_CDS")
generallm_rna <- lm(CDSs_EXPECTED ~ No_CDS * Program + 0, data = tabula_rna.long.2)
anova(generallm_rna)
summary(generallm_rna)

# Graph
tabula2 <- reshape2::melt(tabula, id.vars="CDSs_EXPECTED")
tabula2 <- tabula2[!grepl("ID", tabula2$variable),]
tabula2 <- tabula2[!grepl("Host_domain", tabula2$variable),]
tabula2 <- tabula2[!grepl("Nucleic_acid", tabula2$variable),]

genplot <- ggplot(tabula2, aes(x=CDSs_EXPECTED, y=as.numeric(value),
                               color = as.factor(variable))) +
  xlim(0,4000) + ylim(0,4000) + ylab("CDSs_OBSERVED") + geom_point() +
  geom_abline(intercept=0, size=1, slope=prodigallm$coefficients[1], color='firebrick') +
  geom_abline(intercept=0, size=1, slope=metaprodigallm$coefficients[1], color='orange') +
  geom_abline(intercept=0, size=1, slope=glimmerlm$coefficients[1], color='chartreuse3') +
  geom_abline(intercept=0, size=1, slope=genemarkslm$coefficients[1], color='springgreen4') +
  geom_abline(intercept=0, size=1, slope=phanotatelm$coefficients[1], color='cyan3') +
  geom_abline(intercept=0, size=1, slope=fraggenescanlm$coefficients[1], color='blue') +
  geom_abline(intercept=0, size=1, slope=mgalm$coefficients[1], color='deeppink') +
  geom_abline(intercept=0, size=1, slope=augustuslm$coefficients[1], color='purple') +
  geom_abline(intercept=0, size=1, slope=augustus2lm$coefficients[1], color='lightpink') +
  geom_abline(intercept=0, size=1, slope=augustus3lm$coefficients[1], color='goldenrod') +
  geom_abline(intercept = 0, slope=1, size=1.5, color = "black", linetype = "dashed") +
  theme(legend.position = "top", legend.title = element_blank()) + 
  scale_color_manual(values=c('firebrick', 'orange', 'chartreuse3', 'springgreen4', 'cyan3', 
    'blue', 'deeppink', 'purple', 'lightpink', 'goldenrod'),
                     labels=c("Prodigal", "Metaprodigal", "GLIMMER", "GeneMarkS", 
                      "PHANOTATE", "FragGeneScan", "MGA", "AUGUSTUS - S.aureus", 
                      "AUGUSTUS - E. coli", "AUGUSTUS - H. sapiens")) + 
  labs(x = "", y = "No. CDS (Program)")

tabula_arc.long.2 <- tabula %>%
  pivot_longer(cols=c(5:14), names_to = "Program", values_to = "No_CDS")
generallm_arc <- lm(CDSs_EXPECTED ~ No_CDS * Program + 0, data = tabula_arc.long.2)
anova(generallm_arc)
summary(generallm_arc)

tabula2_arc <- reshape2::melt(tabula_arc, id.vars="CDSs_EXPECTED")
tabula2_arc <- tabula2_arc[!grepl("ID", tabula2_arc$variable),]
tabula2_arc <- tabula2_arc[!grepl("Host_domain", tabula2_arc$variable),]
tabula2_arc <- tabula2_arc[!grepl("Nucleic_acid", tabula2_arc$variable),]

arcplot <- ggplot(tabula2_arc, aes(x=CDSs_EXPECTED, y=as.numeric(value),
                                   color = as.factor(variable))) +
  scale_color_manual(values=c('firebrick', 'orange', 'chartreuse3', 'springgreen4', 
    'cyan3', 'blue', 'deeppink', 'purple', 'lightpink', 'goldenrod')) +
  xlim(0,500) + ylim(0,500) + geom_point() +
  geom_abline(intercept=0, size=1, slope=prodigallm_arc$coefficients[1], color='firebrick') +
  geom_abline(intercept=0, size=1, slope=metaprodigallm_arc$coefficients[1], color='orange') +
  geom_abline(intercept=0, size=1, slope=glimmerlm_arc$coefficients[1], color='chartreuse3') +
  geom_abline(intercept=0, size=1, slope=genemarkslm_arc$coefficients[1], color='springgreen4') +
  geom_abline(intercept=0, size=1, slope=phanotatelm_arc$coefficients[1], color='cyan3') +
  geom_abline(intercept=0, size=1, slope=fraggenescanlm_arc$coefficients[1], color='blue') +
  geom_abline(intercept=0, size=1, slope=mgalm_arc$coefficients[1], color='deeppink') +
  geom_abline(intercept=0, size=1, slope=augustuslm_arc$coefficients[1], color='purple') +
  geom_abline(intercept=0, size=1, slope=augustus2lm_arc$coefficients[1], color='lightpink') +
  geom_abline(intercept=0, size=1, slope=augustus3lm_arc$coefficients[1], color='goldenrod') +
  geom_abline(intercept = 0, slope=1, size=1.5, color = "black", linetype = "dashed") +
  theme(legend.position = "none", legend.title = element_blank()) + 
  labs(x = "", y = "")

tabula2_bac <- reshape2::melt(tabula_bac, id.vars="CDSs_EXPECTED")
tabula2_bac <- tabula2_bac[!grepl("ID", tabula2_bac$variable),]
tabula2_bac <- tabula2_bac[!grepl("Host_domain", tabula2_bac$variable),]
tabula2_bac <- tabula2_bac[!grepl("Nucleic_acid", tabula2_bac$variable),]

bacplot <- ggplot(tabula2_bac, aes(x=CDSs_EXPECTED, y=as.numeric(value),
                                   color = as.factor(variable))) +
  scale_color_manual(values=c('firebrick', 'orange', 'chartreuse3', 'springgreen4', 
                              'cyan3', 'blue', 'deeppink', 'purple', 'lightpink', 'goldenrod')) +
  xlim(0,1000) + ylim(0,1000) + geom_point() +
  geom_abline(intercept=0, size = 1, slope=prodigallm_bac$coefficients[1], color='firebrick') +
  geom_abline(intercept=0, size = 1, slope=metaprodigallm_bac$coefficients[1], color='orange') +
  geom_abline(intercept=0, size = 1, slope=glimmerlm_bac$coefficients[1], color='chartreuse3') +
  geom_abline(intercept=0, size = 1, slope=genemarkslm_bac$coefficients[1], color='springgreen4') +
  geom_abline(intercept=0, size = 1, slope=phanotatelm_bac$coefficients[1], color='cyan3') +
  geom_abline(intercept=0, size = 1, slope=fraggenescanlm_bac$coefficients[1], color='blue') +
  geom_abline(intercept=0, size = 1, slope=mgalm_bac$coefficients[1], color='deeppink') +
  geom_abline(intercept=0, size = 1, slope=augustuslm_bac$coefficients[1], color='purple') +
  geom_abline(intercept=0, size = 1, slope=augustus2lm_bac$coefficients[1], color='lightpink') +
  geom_abline(intercept=0, size = 1, slope=augustus3lm_bac$coefficients[1], color='goldenrod') +
  geom_abline(intercept = 0, slope=1, size = 1.5, color = "black", linetype = "dashed") +
  theme(legend.position = "none", legend.title = element_blank()) + 
  labs(x = "", y = "No. CDS (Program)")

tabula2_euk <- reshape2::melt(tabula_euk, id.vars="CDSs_EXPECTED")
tabula2_euk <- tabula2_euk[!grepl("ID", tabula2_euk$variable),]
tabula2_euk <- tabula2_euk[!grepl("Host_domain", tabula2_euk$variable),]
tabula2_euk <- tabula2_euk[!grepl("Nucleic_acid", tabula2_euk$variable),]

eukplot <- ggplot(tabula2_euk, aes(x=CDSs_EXPECTED, y=as.numeric(value),
                                   color = as.factor(variable))) +
  scale_color_manual(values=c('firebrick', 'orange', 'chartreuse3', 'springgreen4', 
                              'cyan3', 'blue', 'deeppink', 'purple', 'lightpink', 'goldenrod')) +
  xlim(0,4000) + ylim(0,4000) + geom_point() +
  geom_abline(intercept=0, size = 1, slope=prodigallm_euk$coefficients[1], color='firebrick') +
  geom_abline(intercept=0, size = 1, slope=metaprodigallm_euk$coefficients[1], color='orange') +
  geom_abline(intercept=0, size = 1, slope=glimmerlm_euk$coefficients[1], color='chartreuse3') +
  geom_abline(intercept=0, size = 1, slope=genemarkslm_euk$coefficients[1], color='springgreen4') +
  geom_abline(intercept=0, size = 1, slope=phanotatelm_euk$coefficients[1], color='cyan3') +
  geom_abline(intercept=0, size = 1, slope=fraggenescanlm_euk$coefficients[1], color='blue') +
  geom_abline(intercept=0, size = 1, slope=mgalm_euk$coefficients[1], color='deeppink') +
  geom_abline(intercept=0, size = 1, slope=augustuslm_euk$coefficients[1], color='purple') +
  geom_abline(intercept=0, size = 1, slope=augustus2lm_euk$coefficients[1], color='lightpink') +
  geom_abline(intercept=0, size = 1, slope=augustus3lm_euk$coefficients[1], color='goldenrod') +
  geom_abline(intercept = 0, slope=1, size=1.5, color = "black", linetype = "dashed") +
  theme(legend.position = "none", legend.title = element_blank()) + 
  labs(x = "", y = "")

tabula2_dna <- reshape2::melt(tabula_dna, id.vars="CDSs_EXPECTED")
tabula2_dna <- tabula2_dna[!grepl("ID", tabula2_dna$variable),]
tabula2_dna <- tabula2_dna[!grepl("Host_domain", tabula2_dna$variable),]
tabula2_dna <- tabula2_dna[!grepl("Nucleic_acid", tabula2_dna$variable),]

dnaplot <- ggplot(tabula2_dna, aes(x=CDSs_EXPECTED, y=as.numeric(value),
                                   color = as.factor(variable))) +
  scale_color_manual(values=c('firebrick', 'orange', 'chartreuse3', 'springgreen4', 
                              'cyan3', 'blue', 'deeppink', 'purple', 'lightpink', 'goldenrod')) +
  xlim(0,4000) + ylim(0,4000) + geom_point() +
  geom_abline(intercept=0, size = 1, slope=prodigallm_dna$coefficients[1], color='firebrick') +
  geom_abline(intercept=0, size = 1, slope=metaprodigallm_dna$coefficients[1], color='orange') +
  geom_abline(intercept=0, size = 1, slope=glimmerlm_dna$coefficients[1], color='chartreuse3') +
  geom_abline(intercept=0, size = 1, slope=genemarkslm_dna$coefficients[1], color='springgreen4') +
  geom_abline(intercept=0, size = 1, slope=phanotatelm_dna$coefficients[1], color='cyan3') +
  geom_abline(intercept=0, size = 1, slope=fraggenescanlm_dna$coefficients[1], color='blue') +
  geom_abline(intercept=0, size = 1, slope=mgalm_dna$coefficients[1], color='deeppink') +
  geom_abline(intercept=0, size = 1, slope=augustuslm_dna$coefficients[1], color='purple') +
  geom_abline(intercept=0, size = 1, slope=augustus2lm_dna$coefficients[1], color='lightpink') +
  geom_abline(intercept=0, size = 1, slope=augustus3lm_dna$coefficients[1], color='goldenrod') +
  geom_abline(intercept = 0, slope=1, size = 1.5, color = "black", linetype = "dashed") +
  theme(legend.position = "none", legend.title = element_blank()) + 
  labs(x = "No. CDS (RefSeq)", y = "No. CDS (Program)")

tabula2_rna <- reshape2::melt(tabula_rna, id.vars="CDSs_EXPECTED")
tabula2_rna <- tabula2_rna[!grepl("ID", tabula2_rna$variable),]
tabula2_rna <- tabula2_rna[!grepl("Host_domain", tabula2_rna$variable),]
tabula2_rna <- tabula2_rna[!grepl("Nucleic_acid", tabula2_rna$variable),]

rnaplot <- ggplot(tabula2_rna, aes(x=CDSs_EXPECTED, y=as.numeric(value),
                                   color = as.factor(variable))) +
  scale_color_manual(values=c('firebrick', 'orange', 'chartreuse3', 'springgreen4', 
                              'cyan3', 'blue', 'deeppink', 'purple', 'lightpink', 'goldenrod')) +
  xlim(0,75) + ylim(0,75) + geom_point() +
  geom_abline(intercept=0, size = 1, slope=prodigallm_rna$coefficients[1], color='firebrick') +
  geom_abline(intercept=0, size = 1, slope=metaprodigallm_rna$coefficients[1], color='orange') +
  geom_abline(intercept=0, size = 1, slope=glimmerlm_rna$coefficients[1], color='chartreuse3') +
  geom_abline(intercept=0, size = 1, slope=genemarkslm_rna$coefficients[1], color='springgreen4') +
  geom_abline(intercept=0, size = 1, slope=phanotatelm_rna$coefficients[1], color='cyan3') +
  geom_abline(intercept=0, size = 1, slope=fraggenescanlm_rna$coefficients[1], color='blue') +
  geom_abline(intercept=0, size = 1, slope=mgalm_rna$coefficients[1], color='deeppink') +
  geom_abline(intercept=0, size = 1, slope=augustuslm_rna$coefficients[1], color='purple') +
  geom_abline(intercept=0, size = 1, slope=augustus2lm_rna$coefficients[1], color='lightpink') +
  geom_abline(intercept=0, size = 1, slope=augustus3lm_rna$coefficients[1], color='goldenrod') +
  geom_abline(intercept = 0, slope=1, size = 1.5, color = "black", linetype = "dashed") +
  theme(legend.position = "none", legend.title = element_blank()) + 
  labs(x = "No. CDS (RefSeq)", y = "")

legend <- get_legend(genplot)
svg("figure_1.svg", width = 7, height = 11)
topgrid <- plot_grid(genplot + theme(legend.position = "null"), arcplot, bacplot, 
  eukplot, dnaplot, rnaplot, labels = c("A", "B", "C", "D", "E", "F"), ncol = 2)
plot_grid(topgrid, legend, ncol = 1, rel_heights = c(3,.25))
dev.off()

## SECOND TEST: Coordinates

TPtabula_prodigal <- read.csv("C:\\Users/SES089/OneDrive - University of Salford/Documents/Project_LES/VIRUS_RefSeq/viral_genomes.PRODIGAL.TPtable.tsv", header = T, sep = "\t") %>% 
  full_join(tablahost, by = c("Genome" = "ID")) %>%
  full_join(tablanucleic, by = c("Genome" = "ID"))
TPtabula_metaprodigal <- read.csv("C:\\Users/SES089/OneDrive - University of Salford/Documents/Project_LES/VIRUS_RefSeq/viral_genomes.METAPRODIGAL.TPtable.tsv", header = T, sep = "\t") %>% 
  full_join(tablahost, by = c("Genome" = "ID")) %>%
  full_join(tablanucleic, by = c("Genome" = "ID"))
TPtabula_glimmer <- read.csv("C:\\Users/SES089/OneDrive - University of Salford/Documents/Project_LES/VIRUS_RefSeq/viral_genomes.GLIMMER.TPtable.tsv", header = T, sep = "\t") %>% 
  full_join(tablahost, by = c("Genome" = "ID")) %>%
  full_join(tablanucleic, by = c("Genome" = "ID"))
TPtabula_genemarks <- read.csv("C:\\Users/SES089/OneDrive - University of Salford/Documents/Project_LES/VIRUS_RefSeq/viral_genomes.GENEMARKS.TPtable.tsv", header = T, sep = "\t") %>% 
  full_join(tablahost, by = c("Genome" = "ID")) %>%
  full_join(tablanucleic, by = c("Genome" = "ID"))
TPtabula_phanotate <- read.csv("C:\\Users/SES089/OneDrive - University of Salford/Documents/Project_LES/VIRUS_RefSeq/viral_genomes.PHANOTATE.TPtable.tsv", header = T, sep = "\t") %>% 
  full_join(tablahost, by = c("Genome" = "ID")) %>%
  full_join(tablanucleic, by = c("Genome" = "ID"))
TPtabula_fraggenescan <- read.csv("C:\\Users/SES089/OneDrive - University of Salford/Documents/Project_LES/VIRUS_RefSeq/viral_genomes.FRAGGENESCAN.TPtable.tsv", header = T, sep = "\t") %>% 
  full_join(tablahost, by = c("Genome" = "ID")) %>%
  full_join(tablanucleic, by = c("Genome" = "ID"))
TPtabula_mga <- read.csv("C:\\Users/SES089/OneDrive - University of Salford/Documents/Project_LES/VIRUS_RefSeq/viral_genomes.MGA.TPtable.tsv", header = T, sep = "\t") %>% 
  full_join(tablahost, by = c("Genome" = "ID")) %>%
  full_join(tablanucleic, by = c("Genome" = "ID"))
TPtabula_augustus_sa <- read.csv("C:\\Users/SES089/OneDrive - University of Salford/Documents/Project_LES/VIRUS_RefSeq/viral_genomes.AUGUSTUS_SAUREUS.TPtable.tsv", header = T, sep = "\t") %>% 
  full_join(tablahost, by = c("Genome" = "ID")) %>%
  full_join(tablanucleic, by = c("Genome" = "ID"))
TPtabula_augustus_ec <- read.csv("C:\\Users/SES089/OneDrive - University of Salford/Documents/Project_LES/VIRUS_RefSeq/viral_genomes.AUGUSTUS_ECOLI.TPtable.tsv", header = T, sep = "\t") %>% 
  full_join(tablahost, by = c("Genome" = "ID")) %>%
  full_join(tablanucleic, by = c("Genome" = "ID"))
TPtabula_augustus_hs <- read.csv("C:\\Users/SES089/OneDrive - University of Salford/Documents/Project_LES/VIRUS_RefSeq/viral_genomes.AUGUSTUS_HUMAN.TPtable.tsv", header = T, sep = "\t") %>% 
  full_join(tablahost, by = c("Genome" = "ID")) %>%
  full_join(tablanucleic, by = c("Genome" = "ID"))

# General

sum_prodigal <- TPtabula_prodigal %>%
  summarize_if(is.numeric, sum, na.rm=TRUE)
sum_metaprodigal <- TPtabula_metaprodigal %>%
  summarize_if(is.numeric, sum, na.rm=TRUE)
sum_glimmer <- TPtabula_glimmer %>%
  summarize_if(is.numeric, sum, na.rm=TRUE)
sum_genemarks <- TPtabula_genemarks %>%
  summarize_if(is.numeric, sum, na.rm=TRUE)
sum_phanotate <- TPtabula_phanotate %>%
  summarize_if(is.numeric, sum, na.rm=TRUE)
sum_fraggenescan <- TPtabula_fraggenescan %>%
  summarize_if(is.numeric, sum, na.rm=TRUE)
sum_mga <- TPtabula_mga %>%
  summarize_if(is.numeric, sum, na.rm=TRUE)
sum_augustus_sa <- TPtabula_augustus_sa %>%
  summarize_if(is.numeric, sum, na.rm=TRUE)
sum_augustus_ec <- TPtabula_augustus_ec %>%
  summarize_if(is.numeric, sum, na.rm=TRUE)
sum_augustus_hs <- TPtabula_augustus_hs %>%
  summarize_if(is.numeric, sum, na.rm=TRUE)

prodigal_f1 = (2 * sum_prodigal$TP) / (2 * sum_prodigal$TP + sum_prodigal$FP + sum_prodigal$FN)
prodigal_precision = sum_prodigal$TP / (sum_prodigal$TP + sum_prodigal$FP)
prodigal_sensitivity = sum_prodigal$TP / (sum_prodigal$TP + sum_prodigal$FN)
prodigal_fdr = sum_prodigal$FP / (sum_prodigal$TP + sum_prodigal$FP)
prodigal_fnr = sum_prodigal$FN / (sum_prodigal$TP + sum_prodigal$FN)
prodigal_f1
prodigal_precision
prodigal_sensitivity
prodigal_fdr
prodigal_fnr

metaprodigal_f1 = (2 * sum_metaprodigal$TP) / (2 * sum_metaprodigal$TP + sum_metaprodigal$FP + sum_metaprodigal$FN)
metaprodigal_precision = sum_metaprodigal$TP / (sum_metaprodigal$TP + sum_metaprodigal$FP)
metaprodigal_sensitivity = sum_metaprodigal$TP / (sum_metaprodigal$TP + sum_metaprodigal$FN)
metaprodigal_fdr = sum_metaprodigal$FP / (sum_metaprodigal$TP + sum_metaprodigal$FP)
metaprodigal_fnr = sum_metaprodigal$FN / (sum_metaprodigal$TP + sum_metaprodigal$FN)
metaprodigal_f1
metaprodigal_precision
metaprodigal_sensitivity
metaprodigal_fdr
metaprodigal_fnr

glimmer_f1 = (2 * sum_glimmer$TP) / (2 * sum_glimmer$TP + sum_glimmer$FP + sum_glimmer$FN)
glimmer_precision = sum_glimmer$TP / (sum_glimmer$TP + sum_glimmer$FP)
glimmer_sensitivity = sum_glimmer$TP / (sum_glimmer$TP + sum_glimmer$FN)
glimmer_fdr = sum_glimmer$FP / (sum_glimmer$TP + sum_glimmer$FP)
glimmer_fnr = sum_glimmer$FN / (sum_glimmer$TP + sum_glimmer$FN)
glimmer_f1
glimmer_precision
glimmer_sensitivity
glimmer_fdr
glimmer_fnr

genemarks_f1 = (2 * sum_genemarks$TP) / (2 * sum_genemarks$TP + sum_genemarks$FP + sum_genemarks$FN)
genemarks_precision = sum_genemarks$TP / (sum_genemarks$TP + sum_genemarks$FP)
genemarks_sensitivity = sum_genemarks$TP / (sum_genemarks$TP + sum_genemarks$FN)
genemarks_fdr = sum_genemarks$FP / (sum_genemarks$TP + sum_genemarks$FP)
genemarks_fnr = sum_genemarks$FN / (sum_genemarks$TP + sum_genemarks$FN)
genemarks_f1
genemarks_precision
genemarks_sensitivity
genemarks_fdr
genemarks_fnr

phanotate_f1 = (2 * sum_phanotate$TP) / (2 * sum_phanotate$TP + sum_phanotate$FP + sum_phanotate$FN)
phanotate_precision = sum_phanotate$TP / (sum_phanotate$TP + sum_phanotate$FP)
phanotate_sensitivity = sum_phanotate$TP / (sum_phanotate$TP + sum_phanotate$FN)
phanotate_fdr = sum_phanotate$FP / (sum_phanotate$TP + sum_phanotate$FP)
phanotate_fnr = sum_phanotate$FN / (sum_phanotate$TP + sum_phanotate$FN)
phanotate_f1
phanotate_precision
phanotate_sensitivity
phanotate_fdr
phanotate_fnr

fraggenescan_f1 = (2 * sum_fraggenescan$TP) / (2 * sum_fraggenescan$TP + sum_fraggenescan$FP + sum_fraggenescan$FN)
fraggenescan_precision = sum_fraggenescan$TP / (sum_fraggenescan$TP + sum_fraggenescan$FP)
fraggenescan_sensitivity = sum_fraggenescan$TP / (sum_fraggenescan$TP + sum_fraggenescan$FN)
fraggenescan_fdr = sum_fraggenescan$FP / (sum_fraggenescan$TP + sum_fraggenescan$FP)
fraggenescan_fnr = sum_fraggenescan$FN / (sum_fraggenescan$TP + sum_fraggenescan$FN)
fraggenescan_f1
fraggenescan_precision
fraggenescan_sensitivity
fraggenescan_fdr
fraggenescan_fnr

mga_f1 = (2 * sum_mga$TP) / (2 * sum_mga$TP + sum_mga$FP + sum_mga$FN)
mga_precision = sum_mga$TP / (sum_mga$TP + sum_mga$FP)
mga_sensitivity = sum_mga$TP / (sum_mga$TP + sum_mga$FN)
mga_fdr = sum_mga$FP / (sum_mga$TP + sum_mga$FP)
mga_fnr = sum_mga$FN / (sum_mga$TP + sum_mga$FN)
mga_f1
mga_precision
mga_sensitivity
mga_fdr
mga_fnr

augustus_sa_f1 = (2 * sum_augustus_sa$TP) / (2 * sum_augustus_sa$TP + sum_augustus_sa$FP + sum_augustus_sa$FN)
augustus_sa_precision = sum_augustus_sa$TP / (sum_augustus_sa$TP + sum_augustus_sa$FP)
augustus_sa_sensitivity = sum_augustus_sa$TP / (sum_augustus_sa$TP + sum_augustus_sa$FN)
augustus_sa_fdr = sum_augustus_sa$FP / (sum_augustus_sa$TP + sum_augustus_sa$FP)
augustus_sa_fnr = sum_augustus_sa$FN / (sum_augustus_sa$TP + sum_augustus_sa$FN)
augustus_sa_f1
augustus_sa_precision
augustus_sa_sensitivity
augustus_sa_fdr
augustus_sa_fnr

augustus_ec_f1 = (2 * sum_augustus_ec$TP) / (2 * sum_augustus_ec$TP + sum_augustus_ec$FP + sum_augustus_ec$FN)
augustus_ec_precision = sum_augustus_ec$TP / (sum_augustus_ec$TP + sum_augustus_ec$FP)
augustus_ec_sensitivity = sum_augustus_ec$TP / (sum_augustus_ec$TP + sum_augustus_ec$FN)
augustus_ec_fdr = sum_augustus_ec$FP / (sum_augustus_ec$TP + sum_augustus_ec$FP)
augustus_ec_fnr = sum_augustus_ec$FN / (sum_augustus_ec$TP + sum_augustus_ec$FN)
augustus_ec_f1
augustus_ec_precision
augustus_ec_sensitivity
augustus_ec_fdr
augustus_ec_fnr

augustus_hs_f1 = (2 * sum_augustus_hs$TP) / (2 * sum_augustus_hs$TP + sum_augustus_hs$FP + sum_augustus_hs$FN)
augustus_hs_precision = sum_augustus_hs$TP / (sum_augustus_hs$TP + sum_augustus_hs$FP)
augustus_hs_sensitivity = sum_augustus_hs$TP / (sum_augustus_hs$TP + sum_augustus_hs$FN)
augustus_hs_fdr = sum_augustus_hs$FP / (sum_augustus_hs$TP + sum_augustus_hs$FP)
augustus_hs_fnr = sum_augustus_hs$FN / (sum_augustus_hs$TP + sum_augustus_hs$FN)
augustus_hs_f1
augustus_hs_precision
augustus_hs_sensitivity
augustus_hs_fdr
augustus_hs_fnr

# Eukaryotic viruses

sum_prodigal_euk <- TPtabula_prodigal %>%
  filter(Host_domain == "Eukaryote") %>%
  summarize_if(is.numeric, sum, na.rm=TRUE)
sum_metaprodigal_euk <- TPtabula_metaprodigal %>%
  filter(Host_domain == "Eukaryote") %>%
  summarize_if(is.numeric, sum, na.rm=TRUE)
sum_glimmer_euk <- TPtabula_glimmer %>%
  filter(Host_domain == "Eukaryote") %>%
  summarize_if(is.numeric, sum, na.rm=TRUE)
sum_genemarks_euk <- TPtabula_genemarks %>%
  filter(Host_domain == "Eukaryote") %>%
  summarize_if(is.numeric, sum, na.rm=TRUE)
sum_phanotate_euk <- TPtabula_phanotate %>%
  filter(Host_domain == "Eukaryote") %>%
  summarize_if(is.numeric, sum, na.rm=TRUE)
sum_fraggenescan_euk <- TPtabula_fraggenescan %>%
  filter(Host_domain == "Eukaryote") %>%
  summarize_if(is.numeric, sum, na.rm=TRUE)
sum_mga_euk <- TPtabula_mga %>%
  filter(Host_domain == "Eukaryote") %>%
  summarize_if(is.numeric, sum, na.rm=TRUE)
sum_augustus_sa_euk <- TPtabula_augustus_sa %>%
  filter(Host_domain == "Eukaryote") %>%
  summarize_if(is.numeric, sum, na.rm=TRUE)
sum_augustus_ec_euk <- TPtabula_augustus_ec %>%
  filter(Host_domain == "Eukaryote") %>%
  summarize_if(is.numeric, sum, na.rm=TRUE)
sum_augustus_hs_euk <- TPtabula_augustus_hs %>%
  filter(Host_domain == "Eukaryote") %>%
  summarize_if(is.numeric, sum, na.rm=TRUE)

prodigal_euk_f1 = (2 * sum_prodigal_euk$TP) / (2 * sum_prodigal_euk$TP + sum_prodigal_euk$FP + sum_prodigal_euk$FN)
prodigal_euk_precision = sum_prodigal_euk$TP / (sum_prodigal_euk$TP + sum_prodigal_euk$FP)
prodigal_euk_sensitivity = sum_prodigal_euk$TP / (sum_prodigal_euk$TP + sum_prodigal_euk$FN)
prodigal_euk_fdr = sum_prodigal_euk$FP / (sum_prodigal_euk$TP + sum_prodigal_euk$FP)
prodigal_euk_fnr = sum_prodigal_euk$FN / (sum_prodigal_euk$TP + sum_prodigal_euk$FN)
prodigal_euk_f1
prodigal_euk_precision
prodigal_euk_sensitivity
prodigal_euk_fdr
prodigal_euk_fnr

metaprodigal_euk_f1 = (2 * sum_metaprodigal_euk$TP) / (2 * sum_metaprodigal_euk$TP + sum_metaprodigal_euk$FP + sum_metaprodigal_euk$FN)
metaprodigal_euk_precision = sum_metaprodigal_euk$TP / (sum_metaprodigal_euk$TP + sum_metaprodigal_euk$FP)
metaprodigal_euk_sensitivity = sum_metaprodigal_euk$TP / (sum_metaprodigal_euk$TP + sum_metaprodigal_euk$FN)
metaprodigal_euk_fdr = sum_metaprodigal_euk$FP / (sum_metaprodigal_euk$TP + sum_metaprodigal_euk$FP)
metaprodigal_euk_fnr = sum_metaprodigal_euk$FN / (sum_metaprodigal_euk$TP + sum_metaprodigal_euk$FN)
metaprodigal_euk_f1
metaprodigal_euk_precision
metaprodigal_euk_sensitivity
metaprodigal_euk_fdr
metaprodigal_euk_fnr

glimmer_euk_f1 = (2 * sum_glimmer_euk$TP) / (2 * sum_glimmer_euk$TP + sum_glimmer_euk$FP + sum_glimmer_euk$FN)
glimmer_euk_precision = sum_glimmer_euk$TP / (sum_glimmer_euk$TP + sum_glimmer_euk$FP)
glimmer_euk_sensitivity = sum_glimmer_euk$TP / (sum_glimmer_euk$TP + sum_glimmer_euk$FN)
glimmer_euk_fdr = sum_glimmer_euk$FP / (sum_glimmer_euk$TP + sum_glimmer_euk$FP)
glimmer_euk_fnr = sum_glimmer_euk$FN / (sum_glimmer_euk$TP + sum_glimmer_euk$FN)
glimmer_euk_f1
glimmer_euk_precision
glimmer_euk_sensitivity
glimmer_euk_fdr
glimmer_euk_fnr

genemarks_euk_f1 = (2 * sum_genemarks_euk$TP) / (2 * sum_genemarks_euk$TP + sum_genemarks_euk$FP + sum_genemarks_euk$FN)
genemarks_euk_precision = sum_genemarks_euk$TP / (sum_genemarks_euk$TP + sum_genemarks_euk$FP)
genemarks_euk_sensitivity = sum_genemarks_euk$TP / (sum_genemarks_euk$TP + sum_genemarks_euk$FN)
genemarks_euk_fdr = sum_genemarks_euk$FP / (sum_genemarks_euk$TP + sum_genemarks_euk$FP)
genemarks_euk_fnr = sum_genemarks_euk$FN / (sum_genemarks_euk$TP + sum_genemarks_euk$FN)
genemarks_euk_f1
genemarks_euk_precision
genemarks_euk_sensitivity
genemarks_euk_fdr
genemarks_euk_fnr

phanotate_euk_f1 = (2 * sum_phanotate_euk$TP) / (2 * sum_phanotate_euk$TP + sum_phanotate_euk$FP + sum_phanotate_euk$FN)
phanotate_euk_precision = sum_phanotate_euk$TP / (sum_phanotate_euk$TP + sum_phanotate_euk$FP)
phanotate_euk_sensitivity = sum_phanotate_euk$TP / (sum_phanotate_euk$TP + sum_phanotate_euk$FN)
phanotate_euk_fdr = sum_phanotate_euk$FP / (sum_phanotate_euk$TP + sum_phanotate_euk$FP)
phanotate_euk_fnr = sum_phanotate_euk$FN / (sum_phanotate_euk$TP + sum_phanotate_euk$FN)
phanotate_euk_f1
phanotate_euk_precision
phanotate_euk_sensitivity
phanotate_euk_fdr
phanotate_euk_fnr

fraggenescan_euk_f1 = (2 * sum_fraggenescan_euk$TP) / (2 * sum_fraggenescan_euk$TP + sum_fraggenescan_euk$FP + sum_fraggenescan_euk$FN)
fraggenescan_euk_precision = sum_fraggenescan_euk$TP / (sum_fraggenescan_euk$TP + sum_fraggenescan_euk$FP)
fraggenescan_euk_sensitivity = sum_fraggenescan_euk$TP / (sum_fraggenescan_euk$TP + sum_fraggenescan_euk$FN)
fraggenescan_euk_fdr = sum_fraggenescan_euk$FP / (sum_fraggenescan_euk$TP + sum_fraggenescan_euk$FP)
fraggenescan_euk_fnr = sum_fraggenescan_euk$FN / (sum_fraggenescan_euk$TP + sum_fraggenescan_euk$FN)
fraggenescan_euk_f1
fraggenescan_euk_precision
fraggenescan_euk_sensitivity
fraggenescan_euk_fdr
fraggenescan_euk_fnr

mga_euk_f1 = (2 * sum_mga_euk$TP) / (2 * sum_mga_euk$TP + sum_mga_euk$FP + sum_mga_euk$FN)
mga_euk_precision = sum_mga_euk$TP / (sum_mga_euk$TP + sum_mga_euk$FP)
mga_euk_sensitivity = sum_mga_euk$TP / (sum_mga_euk$TP + sum_mga_euk$FN)
mga_euk_fdr = sum_mga_euk$FP / (sum_mga_euk$TP + sum_mga_euk$FP)
mga_euk_fnr = sum_mga_euk$FN / (sum_mga_euk$TP + sum_mga_euk$FN)
mga_euk_f1
mga_euk_precision
mga_euk_sensitivity
mga_euk_fdr
mga_euk_fnr

augustus_sa_euk_f1 = (2 * sum_augustus_sa_euk$TP) / (2 * sum_augustus_sa_euk$TP + sum_augustus_sa_euk$FP + sum_augustus_sa_euk$FN)
augustus_sa_euk_precision = sum_augustus_sa_euk$TP / (sum_augustus_sa_euk$TP + sum_augustus_sa_euk$FP)
augustus_sa_euk_sensitivity = sum_augustus_sa_euk$TP / (sum_augustus_sa_euk$TP + sum_augustus_sa_euk$FN)
augustus_sa_euk_fdr = sum_augustus_sa_euk$FP / (sum_augustus_sa_euk$TP + sum_augustus_sa_euk$FP)
augustus_sa_euk_fnr = sum_augustus_sa_euk$FN / (sum_augustus_sa_euk$TP + sum_augustus_sa_euk$FN)
augustus_sa_euk_f1
augustus_sa_euk_precision
augustus_sa_euk_sensitivity
augustus_sa_euk_fdr
augustus_sa_euk_fnr

augustus_ec_euk_f1 = (2 * sum_augustus_ec_euk$TP) / (2 * sum_augustus_ec_euk$TP + sum_augustus_ec_euk$FP + sum_augustus_ec_euk$FN)
augustus_ec_euk_precision = sum_augustus_ec_euk$TP / (sum_augustus_ec_euk$TP + sum_augustus_ec_euk$FP)
augustus_ec_euk_sensitivity = sum_augustus_ec_euk$TP / (sum_augustus_ec_euk$TP + sum_augustus_ec_euk$FN)
augustus_ec_euk_fdr = sum_augustus_ec_euk$FP / (sum_augustus_ec_euk$TP + sum_augustus_ec_euk$FP)
augustus_ec_euk_fnr = sum_augustus_ec_euk$FN / (sum_augustus_ec_euk$TP + sum_augustus_ec_euk$FN)
augustus_ec_euk_f1
augustus_ec_euk_precision
augustus_ec_euk_sensitivity
augustus_ec_euk_fdr
augustus_ec_euk_fnr

augustus_hs_euk_f1 = (2 * sum_augustus_hs_euk$TP) / (2 * sum_augustus_hs_euk$TP + sum_augustus_hs_euk$FP + sum_augustus_hs_euk$FN)
augustus_hs_euk_precision = sum_augustus_hs_euk$TP / (sum_augustus_hs_euk$TP + sum_augustus_hs_euk$FP)
augustus_hs_euk_sensitivity = sum_augustus_hs_euk$TP / (sum_augustus_hs_euk$TP + sum_augustus_hs_euk$FN)
augustus_hs_euk_fdr = sum_augustus_hs_euk$FP / (sum_augustus_hs_euk$TP + sum_augustus_hs_euk$FP)
augustus_hs_euk_fnr = sum_augustus_hs_euk$FN / (sum_augustus_hs_euk$TP + sum_augustus_hs_euk$FN)
augustus_hs_euk_f1
augustus_hs_euk_precision
augustus_hs_euk_sensitivity
augustus_hs_euk_fdr
augustus_hs_euk_fnr

# Archaeal viruses

sum_prodigal_arc <- TPtabula_prodigal %>%
  filter(Host_domain == "Archaea") %>%
  summarize_if(is.numeric, sum, na.rm=TRUE)
sum_metaprodigal_arc <- TPtabula_metaprodigal %>%
  filter(Host_domain == "Archaea") %>%
  summarize_if(is.numeric, sum, na.rm=TRUE)
sum_glimmer_arc <- TPtabula_glimmer %>%
  filter(Host_domain == "Archaea") %>%
  summarize_if(is.numeric, sum, na.rm=TRUE)
sum_genemarks_arc <- TPtabula_genemarks %>%
  filter(Host_domain == "Archaea") %>%
  summarize_if(is.numeric, sum, na.rm=TRUE)
sum_phanotate_arc <- TPtabula_phanotate %>%
  filter(Host_domain == "Archaea") %>%
  summarize_if(is.numeric, sum, na.rm=TRUE)
sum_fraggenescan_arc <- TPtabula_fraggenescan %>%
  filter(Host_domain == "Archaea") %>%
  summarize_if(is.numeric, sum, na.rm=TRUE)
sum_mga_arc <- TPtabula_mga %>%
  filter(Host_domain == "Archaea") %>%
  summarize_if(is.numeric, sum, na.rm=TRUE)
sum_augustus_sa_arc <- TPtabula_augustus_sa %>%
  filter(Host_domain == "Archaea") %>%
  summarize_if(is.numeric, sum, na.rm=TRUE)
sum_augustus_ec_arc <- TPtabula_augustus_ec %>%
  filter(Host_domain == "Archaea") %>%
  summarize_if(is.numeric, sum, na.rm=TRUE)
sum_augustus_hs_arc <- TPtabula_augustus_hs %>%
  filter(Host_domain == "Archaea") %>%
  summarize_if(is.numeric, sum, na.rm=TRUE)

prodigal_arc_f1 = (2 * sum_prodigal_arc$TP) / (2 * sum_prodigal_arc$TP + sum_prodigal_arc$FP + sum_prodigal_arc$FN)
prodigal_arc_precision = sum_prodigal_arc$TP / (sum_prodigal_arc$TP + sum_prodigal_arc$FP)
prodigal_arc_sensitivity = sum_prodigal_arc$TP / (sum_prodigal_arc$TP + sum_prodigal_arc$FN)
prodigal_arc_fdr = sum_prodigal_arc$FP / (sum_prodigal_arc$TP + sum_prodigal_arc$FP)
prodigal_arc_fnr = sum_prodigal_arc$FN / (sum_prodigal_arc$TP + sum_prodigal_arc$FN)
prodigal_arc_f1
prodigal_arc_precision
prodigal_arc_sensitivity
prodigal_arc_fdr
prodigal_arc_fnr

metaprodigal_arc_f1 = (2 * sum_metaprodigal_arc$TP) / (2 * sum_metaprodigal_arc$TP + sum_metaprodigal_arc$FP + sum_metaprodigal_arc$FN)
metaprodigal_arc_precision = sum_metaprodigal_arc$TP / (sum_metaprodigal_arc$TP + sum_metaprodigal_arc$FP)
metaprodigal_arc_sensitivity = sum_metaprodigal_arc$TP / (sum_metaprodigal_arc$TP + sum_metaprodigal_arc$FN)
metaprodigal_arc_fdr = sum_metaprodigal_arc$FP / (sum_metaprodigal_arc$TP + sum_metaprodigal_arc$FP)
metaprodigal_arc_fnr = sum_metaprodigal_arc$FN / (sum_metaprodigal_arc$TP + sum_metaprodigal_arc$FN)
metaprodigal_arc_f1
metaprodigal_arc_precision
metaprodigal_arc_sensitivity
metaprodigal_arc_fdr
metaprodigal_arc_fnr

glimmer_arc_f1 = (2 * sum_glimmer_arc$TP) / (2 * sum_glimmer_arc$TP + sum_glimmer_arc$FP + sum_glimmer_arc$FN)
glimmer_arc_precision = sum_glimmer_arc$TP / (sum_glimmer_arc$TP + sum_glimmer_arc$FP)
glimmer_arc_sensitivity = sum_glimmer_arc$TP / (sum_glimmer_arc$TP + sum_glimmer_arc$FN)
glimmer_arc_fdr = sum_glimmer_arc$FP / (sum_glimmer_arc$TP + sum_glimmer_arc$FP)
glimmer_arc_fnr = sum_glimmer_arc$FN / (sum_glimmer_arc$TP + sum_glimmer_arc$FN)
glimmer_arc_f1
glimmer_arc_precision
glimmer_arc_sensitivity
glimmer_arc_fdr
glimmer_arc_fnr

genemarks_arc_f1 = (2 * sum_genemarks_arc$TP) / (2 * sum_genemarks_arc$TP + sum_genemarks_arc$FP + sum_genemarks_arc$FN)
genemarks_arc_precision = sum_genemarks_arc$TP / (sum_genemarks_arc$TP + sum_genemarks_arc$FP)
genemarks_arc_sensitivity = sum_genemarks_arc$TP / (sum_genemarks_arc$TP + sum_genemarks_arc$FN)
genemarks_arc_fdr = sum_genemarks_arc$FP / (sum_genemarks_arc$TP + sum_genemarks_arc$FP)
genemarks_arc_fnr = sum_genemarks_arc$FN / (sum_genemarks_arc$TP + sum_genemarks_arc$FN)
genemarks_arc_f1
genemarks_arc_precision
genemarks_arc_sensitivity
genemarks_arc_fdr
genemarks_arc_fnr

phanotate_arc_f1 = (2 * sum_phanotate_arc$TP) / (2 * sum_phanotate_arc$TP + sum_phanotate_arc$FP + sum_phanotate_arc$FN)
phanotate_arc_precision = sum_phanotate_arc$TP / (sum_phanotate_arc$TP + sum_phanotate_arc$FP)
phanotate_arc_sensitivity = sum_phanotate_arc$TP / (sum_phanotate_arc$TP + sum_phanotate_arc$FN)
phanotate_arc_fdr = sum_phanotate_arc$FP / (sum_phanotate_arc$TP + sum_phanotate_arc$FP)
phanotate_arc_fnr = sum_phanotate_arc$FN / (sum_phanotate_arc$TP + sum_phanotate_arc$FN)
phanotate_arc_f1
phanotate_arc_precision
phanotate_arc_sensitivity
phanotate_arc_fdr
phanotate_arc_fnr

fraggenescan_arc_f1 = (2 * sum_fraggenescan_arc$TP) / (2 * sum_fraggenescan_arc$TP + sum_fraggenescan_arc$FP + sum_fraggenescan_arc$FN)
fraggenescan_arc_precision = sum_fraggenescan_arc$TP / (sum_fraggenescan_arc$TP + sum_fraggenescan_arc$FP)
fraggenescan_arc_sensitivity = sum_fraggenescan_arc$TP / (sum_fraggenescan_arc$TP + sum_fraggenescan_arc$FN)
fraggenescan_arc_fdr = sum_fraggenescan_arc$FP / (sum_fraggenescan_arc$TP + sum_fraggenescan_arc$FP)
fraggenescan_arc_fnr = sum_fraggenescan_arc$FN / (sum_fraggenescan_arc$TP + sum_fraggenescan_arc$FN)
fraggenescan_arc_f1
fraggenescan_arc_precision
fraggenescan_arc_sensitivity
fraggenescan_arc_fdr
fraggenescan_arc_fnr

mga_arc_f1 = (2 * sum_mga_arc$TP) / (2 * sum_mga_arc$TP + sum_mga_arc$FP + sum_mga_arc$FN)
mga_arc_precision = sum_mga_arc$TP / (sum_mga_arc$TP + sum_mga_arc$FP)
mga_arc_sensitivity = sum_mga_arc$TP / (sum_mga_arc$TP + sum_mga_arc$FN)
mga_arc_fdr = sum_mga_arc$FP / (sum_mga_arc$TP + sum_mga_arc$FP)
mga_arc_fnr = sum_mga_arc$FN / (sum_mga_arc$TP + sum_mga_arc$FN)
mga_arc_f1
mga_arc_precision
mga_arc_sensitivity
mga_arc_fdr
mga_arc_fnr

augustus_sa_arc_f1 = (2 * sum_augustus_sa_arc$TP) / (2 * sum_augustus_sa_arc$TP + sum_augustus_sa_arc$FP + sum_augustus_sa_arc$FN)
augustus_sa_arc_precision = sum_augustus_sa_arc$TP / (sum_augustus_sa_arc$TP + sum_augustus_sa_arc$FP)
augustus_sa_arc_sensitivity = sum_augustus_sa_arc$TP / (sum_augustus_sa_arc$TP + sum_augustus_sa_arc$FN)
augustus_sa_arc_fdr = sum_augustus_sa_arc$FP / (sum_augustus_sa_arc$TP + sum_augustus_sa_arc$FP)
augustus_sa_arc_fnr = sum_augustus_sa_arc$FN / (sum_augustus_sa_arc$TP + sum_augustus_sa_arc$FN)
augustus_sa_arc_f1
augustus_sa_arc_precision
augustus_sa_arc_sensitivity
augustus_sa_arc_fdr
augustus_sa_arc_fnr

augustus_ec_arc_f1 = (2 * sum_augustus_ec_arc$TP) / (2 * sum_augustus_ec_arc$TP + sum_augustus_ec_arc$FP + sum_augustus_ec_arc$FN)
augustus_ec_arc_precision = sum_augustus_ec_arc$TP / (sum_augustus_ec_arc$TP + sum_augustus_ec_arc$FP)
augustus_ec_arc_sensitivity = sum_augustus_ec_arc$TP / (sum_augustus_ec_arc$TP + sum_augustus_ec_arc$FN)
augustus_ec_arc_fdr = sum_augustus_ec_arc$FP / (sum_augustus_ec_arc$TP + sum_augustus_ec_arc$FP)
augustus_ec_arc_fnr = sum_augustus_ec_arc$FN / (sum_augustus_ec_arc$TP + sum_augustus_ec_arc$FN)
augustus_ec_arc_f1
augustus_ec_arc_precision
augustus_ec_arc_sensitivity
augustus_ec_arc_fdr
augustus_ec_arc_fnr

augustus_hs_arc_f1 = (2 * sum_augustus_hs_arc$TP) / (2 * sum_augustus_hs_arc$TP + sum_augustus_hs_arc$FP + sum_augustus_hs_arc$FN)
augustus_hs_arc_precision = sum_augustus_hs_arc$TP / (sum_augustus_hs_arc$TP + sum_augustus_hs_arc$FP)
augustus_hs_arc_sensitivity = sum_augustus_hs_arc$TP / (sum_augustus_hs_arc$TP + sum_augustus_hs_arc$FN)
augustus_hs_arc_fdr = sum_augustus_hs_arc$FP / (sum_augustus_hs_arc$TP + sum_augustus_hs_arc$FP)
augustus_hs_arc_fnr = sum_augustus_hs_arc$FN / (sum_augustus_hs_arc$TP + sum_augustus_hs_arc$FN)
augustus_hs_arc_f1
augustus_hs_arc_precision
augustus_hs_arc_sensitivity
augustus_hs_arc_fdr
augustus_hs_arc_fnr

# Bacteriophages

sum_prodigal_bac <- TPtabula_prodigal %>%
  filter(Host_domain == "Bacteria") %>%
  summarize_if(is.numeric, sum, na.rm=TRUE)
sum_metaprodigal_bac <- TPtabula_metaprodigal %>%
  filter(Host_domain == "Bacteria") %>%
  summarize_if(is.numeric, sum, na.rm=TRUE)
sum_glimmer_bac <- TPtabula_glimmer %>%
  filter(Host_domain == "Bacteria") %>%
  summarize_if(is.numeric, sum, na.rm=TRUE)
sum_genemarks_bac <- TPtabula_genemarks %>%
  filter(Host_domain == "Bacteria") %>%
  summarize_if(is.numeric, sum, na.rm=TRUE)
sum_phanotate_bac <- TPtabula_phanotate %>%
  filter(Host_domain == "Bacteria") %>%
  summarize_if(is.numeric, sum, na.rm=TRUE)
sum_fraggenescan_bac <- TPtabula_fraggenescan %>%
  filter(Host_domain == "Eukaryote") %>%
  summarize_if(is.numeric, sum, na.rm=TRUE)
sum_mga_bac <- TPtabula_mga %>%
  filter(Host_domain == "Bacteria") %>%
  summarize_if(is.numeric, sum, na.rm=TRUE)
sum_augustus_sa_bac <- TPtabula_augustus_sa %>%
  filter(Host_domain == "Bacteria") %>%
  summarize_if(is.numeric, sum, na.rm=TRUE)
sum_augustus_ec_bac <- TPtabula_augustus_ec %>%
  filter(Host_domain == "Bacteria") %>%
  summarize_if(is.numeric, sum, na.rm=TRUE)
sum_augustus_hs_bac <- TPtabula_augustus_hs %>%
  filter(Host_domain == "Bacteria") %>%
  summarize_if(is.numeric, sum, na.rm=TRUE)

prodigal_bac_f1 = (2 * sum_prodigal_bac$TP) / (2 * sum_prodigal_bac$TP + sum_prodigal_bac$FP + sum_prodigal_bac$FN)
prodigal_bac_precision = sum_prodigal_bac$TP / (sum_prodigal_bac$TP + sum_prodigal_bac$FP)
prodigal_bac_sensitivity = sum_prodigal_bac$TP / (sum_prodigal_bac$TP + sum_prodigal_bac$FN)
prodigal_bac_fdr = sum_prodigal_bac$FP / (sum_prodigal_bac$TP + sum_prodigal_bac$FP)
prodigal_bac_fnr = sum_prodigal_bac$FN / (sum_prodigal_bac$TP + sum_prodigal_bac$FN)
prodigal_bac_f1
prodigal_bac_precision
prodigal_bac_sensitivity
prodigal_bac_fdr
prodigal_bac_fnr

metaprodigal_bac_f1 = (2 * sum_metaprodigal_bac$TP) / (2 * sum_metaprodigal_bac$TP + sum_metaprodigal_bac$FP + sum_metaprodigal_bac$FN)
metaprodigal_bac_precision = sum_metaprodigal_bac$TP / (sum_metaprodigal_bac$TP + sum_metaprodigal_bac$FP)
metaprodigal_bac_sensitivity = sum_metaprodigal_bac$TP / (sum_metaprodigal_bac$TP + sum_metaprodigal_bac$FN)
metaprodigal_bac_fdr = sum_metaprodigal_bac$FP / (sum_metaprodigal_bac$TP + sum_metaprodigal_bac$FP)
metaprodigal_bac_fnr = sum_metaprodigal_bac$FN / (sum_metaprodigal_bac$TP + sum_metaprodigal_bac$FN)
metaprodigal_bac_f1
metaprodigal_bac_precision
metaprodigal_bac_sensitivity
metaprodigal_bac_fdr
metaprodigal_bac_fnr

glimmer_bac_f1 = (2 * sum_glimmer_bac$TP) / (2 * sum_glimmer_bac$TP + sum_glimmer_bac$FP + sum_glimmer_bac$FN)
glimmer_bac_precision = sum_glimmer_bac$TP / (sum_glimmer_bac$TP + sum_glimmer_bac$FP)
glimmer_bac_sensitivity = sum_glimmer_bac$TP / (sum_glimmer_bac$TP + sum_glimmer_bac$FN)
glimmer_bac_fdr = sum_glimmer_bac$FP / (sum_glimmer_bac$TP + sum_glimmer_bac$FP)
glimmer_bac_fnr = sum_glimmer_bac$FN / (sum_glimmer_bac$TP + sum_glimmer_bac$FN)
glimmer_bac_f1
glimmer_bac_precision
glimmer_bac_sensitivity
glimmer_bac_fdr
glimmer_bac_fnr

genemarks_bac_f1 = (2 * sum_genemarks_bac$TP) / (2 * sum_genemarks_bac$TP + sum_genemarks_bac$FP + sum_genemarks_bac$FN)
genemarks_bac_precision = sum_genemarks_bac$TP / (sum_genemarks_bac$TP + sum_genemarks_bac$FP)
genemarks_bac_sensitivity = sum_genemarks_bac$TP / (sum_genemarks_bac$TP + sum_genemarks_bac$FN)
genemarks_bac_fdr = sum_genemarks_bac$FP / (sum_genemarks_bac$TP + sum_genemarks_bac$FP)
genemarks_bac_fnr = sum_genemarks_bac$FN / (sum_genemarks_bac$TP + sum_genemarks_bac$FN)
genemarks_bac_f1
genemarks_bac_precision
genemarks_bac_sensitivity
genemarks_bac_fdr
genemarks_bac_fnr

phanotate_bac_f1 = (2 * sum_phanotate_bac$TP) / (2 * sum_phanotate_bac$TP + sum_phanotate_bac$FP + sum_phanotate_bac$FN)
phanotate_bac_precision = sum_phanotate_bac$TP / (sum_phanotate_bac$TP + sum_phanotate_bac$FP)
phanotate_bac_sensitivity = sum_phanotate_bac$TP / (sum_phanotate_bac$TP + sum_phanotate_bac$FN)
phanotate_bac_fdr = sum_phanotate_bac$FP / (sum_phanotate_bac$TP + sum_phanotate_bac$FP)
phanotate_bac_fnr = sum_phanotate_bac$FN / (sum_phanotate_bac$TP + sum_phanotate_bac$FN)
phanotate_bac_f1
phanotate_bac_precision
phanotate_bac_sensitivity
phanotate_bac_fdr
phanotate_bac_fnr

fraggenescan_bac_f1 = (2 * sum_fraggenescan_bac$TP) / (2 * sum_fraggenescan_bac$TP + sum_fraggenescan_bac$FP + sum_fraggenescan_bac$FN)
fraggenescan_bac_precision = sum_fraggenescan_bac$TP / (sum_fraggenescan_bac$TP + sum_fraggenescan_bac$FP)
fraggenescan_bac_sensitivity = sum_fraggenescan_bac$TP / (sum_fraggenescan_bac$TP + sum_fraggenescan_bac$FN)
fraggenescan_bac_fdr = sum_fraggenescan_bac$FP / (sum_fraggenescan_bac$TP + sum_fraggenescan_bac$FP)
fraggenescan_bac_fnr = sum_fraggenescan_bac$FN / (sum_fraggenescan_bac$TP + sum_fraggenescan_bac$FN)
fraggenescan_bac_f1
fraggenescan_bac_precision
fraggenescan_bac_sensitivity
fraggenescan_bac_fdr
fraggenescan_bac_fnr

mga_bac_f1 = (2 * sum_mga_bac$TP) / (2 * sum_mga_bac$TP + sum_mga_bac$FP + sum_mga_bac$FN)
mga_bac_precision = sum_mga_bac$TP / (sum_mga_bac$TP + sum_mga_bac$FP)
mga_bac_sensitivity = sum_mga_bac$TP / (sum_mga_bac$TP + sum_mga_bac$FN)
mga_bac_fdr = sum_mga_bac$FP / (sum_mga_bac$TP + sum_mga_bac$FP)
mga_bac_fnr = sum_mga_bac$FN / (sum_mga_bac$TP + sum_mga_bac$FN)
mga_bac_f1
mga_bac_precision
mga_bac_sensitivity
mga_bac_fdr
mga_bac_fnr

augustus_sa_bac_f1 = (2 * sum_augustus_sa_bac$TP) / (2 * sum_augustus_sa_bac$TP + sum_augustus_sa_bac$FP + sum_augustus_sa_bac$FN)
augustus_sa_bac_precision = sum_augustus_sa_bac$TP / (sum_augustus_sa_bac$TP + sum_augustus_sa_bac$FP)
augustus_sa_bac_sensitivity = sum_augustus_sa_bac$TP / (sum_augustus_sa_bac$TP + sum_augustus_sa_bac$FN)
augustus_sa_bac_fdr = sum_augustus_sa_bac$FP / (sum_augustus_sa_bac$TP + sum_augustus_sa_bac$FP)
augustus_sa_bac_fnr = sum_augustus_sa_bac$FN / (sum_augustus_sa_bac$TP + sum_augustus_sa_bac$FN)
augustus_sa_bac_f1
augustus_sa_bac_precision
augustus_sa_bac_sensitivity
augustus_sa_bac_fdr
augustus_sa_bac_fnr

augustus_ec_bac_f1 = (2 * sum_augustus_ec_bac$TP) / (2 * sum_augustus_ec_bac$TP + sum_augustus_ec_bac$FP + sum_augustus_ec_bac$FN)
augustus_ec_bac_precision = sum_augustus_ec_bac$TP / (sum_augustus_ec_bac$TP + sum_augustus_ec_bac$FP)
augustus_ec_bac_sensitivity = sum_augustus_ec_bac$TP / (sum_augustus_ec_bac$TP + sum_augustus_ec_bac$FN)
augustus_ec_bac_fdr = sum_augustus_ec_bac$FP / (sum_augustus_ec_bac$TP + sum_augustus_ec_bac$FP)
augustus_ec_bac_fnr = sum_augustus_ec_bac$FN / (sum_augustus_ec_bac$TP + sum_augustus_ec_bac$FN)
augustus_ec_bac_f1
augustus_ec_bac_precision
augustus_ec_bac_sensitivity
augustus_ec_bac_fdr
augustus_ec_bac_fnr

augustus_hs_bac_f1 = (2 * sum_augustus_hs_bac$TP) / (2 * sum_augustus_hs_bac$TP + sum_augustus_hs_bac$FP + sum_augustus_hs_bac$FN)
augustus_hs_bac_precision = sum_augustus_hs_bac$TP / (sum_augustus_hs_bac$TP + sum_augustus_hs_bac$FP)
augustus_hs_bac_sensitivity = sum_augustus_hs_bac$TP / (sum_augustus_hs_bac$TP + sum_augustus_hs_bac$FN)
augustus_hs_bac_fdr = sum_augustus_hs_bac$FP / (sum_augustus_hs_bac$TP + sum_augustus_hs_bac$FP)
augustus_hs_bac_fnr = sum_augustus_hs_bac$FN / (sum_augustus_hs_bac$TP + sum_augustus_hs_bac$FN)
augustus_hs_bac_f1
augustus_hs_bac_precision
augustus_hs_bac_sensitivity
augustus_hs_bac_fdr
augustus_hs_bac_fnr

# DNA viruses

sum_prodigal_dna <- TPtabula_prodigal %>%
  filter(Nucleic_acid == "DNA") %>%
  summarize_if(is.numeric, sum, na.rm=TRUE)
sum_metaprodigal_dna <- TPtabula_metaprodigal %>%
  filter(Nucleic_acid == "DNA") %>%
  summarize_if(is.numeric, sum, na.rm=TRUE)
sum_glimmer_dna <- TPtabula_glimmer %>%
  filter(Nucleic_acid == "DNA") %>%
  summarize_if(is.numeric, sum, na.rm=TRUE)
sum_genemarks_dna <- TPtabula_genemarks %>%
  filter(Nucleic_acid == "DNA") %>%
  summarize_if(is.numeric, sum, na.rm=TRUE)
sum_phanotate_dna <- TPtabula_phanotate %>%
  filter(Nucleic_acid == "DNA") %>%
  summarize_if(is.numeric, sum, na.rm=TRUE)
sum_fraggenescan_dna <- TPtabula_fraggenescan %>%
  filter(Nucleic_acid == "DNA") %>%
  summarize_if(is.numeric, sum, na.rm=TRUE)
sum_mga_dna <- TPtabula_mga %>%
  filter(Nucleic_acid == "DNA") %>%
  summarize_if(is.numeric, sum, na.rm=TRUE)
sum_augustus_sa_dna <- TPtabula_augustus_sa %>%
  filter(Nucleic_acid == "DNA") %>%
  summarize_if(is.numeric, sum, na.rm=TRUE)
sum_augustus_ec_dna <- TPtabula_augustus_ec %>%
  filter(Nucleic_acid == "DNA") %>%
  summarize_if(is.numeric, sum, na.rm=TRUE)
sum_augustus_hs_dna <- TPtabula_augustus_hs %>%
  filter(Nucleic_acid == "DNA") %>%
  summarize_if(is.numeric, sum, na.rm=TRUE)

prodigal_dna_f1 = (2 * sum_prodigal_dna$TP) / (2 * sum_prodigal_dna$TP + sum_prodigal_dna$FP + sum_prodigal_dna$FN)
prodigal_dna_precision = sum_prodigal_dna$TP / (sum_prodigal_dna$TP + sum_prodigal_dna$FP)
prodigal_dna_sensitivity = sum_prodigal_dna$TP / (sum_prodigal_dna$TP + sum_prodigal_dna$FN)
prodigal_dna_fdr = sum_prodigal_dna$FP / (sum_prodigal_dna$TP + sum_prodigal_dna$FP)
prodigal_dna_fnr = sum_prodigal_dna$FN / (sum_prodigal_dna$TP + sum_prodigal_dna$FN)
prodigal_dna_f1
prodigal_dna_precision
prodigal_dna_sensitivity
prodigal_dna_fdr
prodigal_dna_fnr

metaprodigal_dna_f1 = (2 * sum_metaprodigal_dna$TP) / (2 * sum_metaprodigal_dna$TP + sum_metaprodigal_dna$FP + sum_metaprodigal_dna$FN)
metaprodigal_dna_precision = sum_metaprodigal_dna$TP / (sum_metaprodigal_dna$TP + sum_metaprodigal_dna$FP)
metaprodigal_dna_sensitivity = sum_metaprodigal_dna$TP / (sum_metaprodigal_dna$TP + sum_metaprodigal_dna$FN)
metaprodigal_dna_fdr = sum_metaprodigal_dna$FP / (sum_metaprodigal_dna$TP + sum_metaprodigal_dna$FP)
metaprodigal_dna_fnr = sum_metaprodigal_dna$FN / (sum_metaprodigal_dna$TP + sum_metaprodigal_dna$FN)
metaprodigal_dna_f1
metaprodigal_dna_precision
metaprodigal_dna_sensitivity
metaprodigal_dna_fdr
metaprodigal_dna_fnr

glimmer_dna_f1 = (2 * sum_glimmer_dna$TP) / (2 * sum_glimmer_dna$TP + sum_glimmer_dna$FP + sum_glimmer_dna$FN)
glimmer_dna_precision = sum_glimmer_dna$TP / (sum_glimmer_dna$TP + sum_glimmer_dna$FP)
glimmer_dna_sensitivity = sum_glimmer_dna$TP / (sum_glimmer_dna$TP + sum_glimmer_dna$FN)
glimmer_dna_fdr = sum_glimmer_dna$FP / (sum_glimmer_dna$TP + sum_glimmer_dna$FP)
glimmer_dna_fnr = sum_glimmer_dna$FN / (sum_glimmer_dna$TP + sum_glimmer_dna$FN)
glimmer_dna_f1
glimmer_dna_precision
glimmer_dna_sensitivity
glimmer_dna_fdr
glimmer_dna_fnr

genemarks_dna_f1 = (2 * sum_genemarks_dna$TP) / (2 * sum_genemarks_dna$TP + sum_genemarks_dna$FP + sum_genemarks_dna$FN)
genemarks_dna_precision = sum_genemarks_dna$TP / (sum_genemarks_dna$TP + sum_genemarks_dna$FP)
genemarks_dna_sensitivity = sum_genemarks_dna$TP / (sum_genemarks_dna$TP + sum_genemarks_dna$FN)
genemarks_dna_fdr = sum_genemarks_dna$FP / (sum_genemarks_dna$TP + sum_genemarks_dna$FP)
genemarks_dna_fnr = sum_genemarks_dna$FN / (sum_genemarks_dna$TP + sum_genemarks_dna$FN)
genemarks_dna_f1
genemarks_dna_precision
genemarks_dna_sensitivity
genemarks_dna_fdr
genemarks_dna_fnr

phanotate_dna_f1 = (2 * sum_phanotate_dna$TP) / (2 * sum_phanotate_dna$TP + sum_phanotate_dna$FP + sum_phanotate_dna$FN)
phanotate_dna_precision = sum_phanotate_dna$TP / (sum_phanotate_dna$TP + sum_phanotate_dna$FP)
phanotate_dna_sensitivity = sum_phanotate_dna$TP / (sum_phanotate_dna$TP + sum_phanotate_dna$FN)
phanotate_dna_fdr = sum_phanotate_dna$FP / (sum_phanotate_dna$TP + sum_phanotate_dna$FP)
phanotate_dna_fnr = sum_phanotate_dna$FN / (sum_phanotate_dna$TP + sum_phanotate_dna$FN)
phanotate_dna_f1
phanotate_dna_precision
phanotate_dna_sensitivity
phanotate_dna_fdr
phanotate_dna_fnr

fraggenescan_dna_f1 = (2 * sum_fraggenescan_dna$TP) / (2 * sum_fraggenescan_dna$TP + sum_fraggenescan_dna$FP + sum_fraggenescan_dna$FN)
fraggenescan_dna_precision = sum_fraggenescan_dna$TP / (sum_fraggenescan_dna$TP + sum_fraggenescan_dna$FP)
fraggenescan_dna_sensitivity = sum_fraggenescan_dna$TP / (sum_fraggenescan_dna$TP + sum_fraggenescan_dna$FN)
fraggenescan_dna_fdr = sum_fraggenescan_dna$FP / (sum_fraggenescan_dna$TP + sum_fraggenescan_dna$FP)
fraggenescan_dna_fnr = sum_fraggenescan_dna$FN / (sum_fraggenescan_dna$TP + sum_fraggenescan_dna$FN)
fraggenescan_dna_f1
fraggenescan_dna_precision
fraggenescan_dna_sensitivity
fraggenescan_dna_fdr
fraggenescan_dna_fnr

mga_dna_f1 = (2 * sum_mga_dna$TP) / (2 * sum_mga_dna$TP + sum_mga_dna$FP + sum_mga_dna$FN)
mga_dna_precision = sum_mga_dna$TP / (sum_mga_dna$TP + sum_mga_dna$FP)
mga_dna_sensitivity = sum_mga_dna$TP / (sum_mga_dna$TP + sum_mga_dna$FN)
mga_dna_fdr = sum_mga_dna$FP / (sum_mga_dna$TP + sum_mga_dna$FP)
mga_dna_fnr = sum_mga_dna$FN / (sum_mga_dna$TP + sum_mga_dna$FN)
mga_dna_f1
mga_dna_precision
mga_dna_sensitivity
mga_dna_fdr
mga_dna_fnr

augustus_sa_dna_f1 = (2 * sum_augustus_sa_dna$TP) / (2 * sum_augustus_sa_dna$TP + sum_augustus_sa_dna$FP + sum_augustus_sa_dna$FN)
augustus_sa_dna_precision = sum_augustus_sa_dna$TP / (sum_augustus_sa_dna$TP + sum_augustus_sa_dna$FP)
augustus_sa_dna_sensitivity = sum_augustus_sa_dna$TP / (sum_augustus_sa_dna$TP + sum_augustus_sa_dna$FN)
augustus_sa_dna_fdr = sum_augustus_sa_dna$FP / (sum_augustus_sa_dna$TP + sum_augustus_sa_dna$FP)
augustus_sa_dna_fnr = sum_augustus_sa_dna$FN / (sum_augustus_sa_dna$TP + sum_augustus_sa_dna$FN)
augustus_sa_dna_f1
augustus_sa_dna_precision
augustus_sa_dna_sensitivity
augustus_sa_dna_fdr
augustus_sa_dna_fnr

augustus_ec_dna_f1 = (2 * sum_augustus_ec_dna$TP) / (2 * sum_augustus_ec_dna$TP + sum_augustus_ec_dna$FP + sum_augustus_ec_dna$FN)
augustus_ec_dna_precision = sum_augustus_ec_dna$TP / (sum_augustus_ec_dna$TP + sum_augustus_ec_dna$FP)
augustus_ec_dna_sensitivity = sum_augustus_ec_dna$TP / (sum_augustus_ec_dna$TP + sum_augustus_ec_dna$FN)
augustus_ec_dna_fdr = sum_augustus_ec_dna$FP / (sum_augustus_ec_dna$TP + sum_augustus_ec_dna$FP)
augustus_ec_dna_fnr = sum_augustus_ec_dna$FN / (sum_augustus_ec_dna$TP + sum_augustus_ec_dna$FN)
augustus_ec_dna_f1
augustus_ec_dna_precision
augustus_ec_dna_sensitivity
augustus_ec_dna_fdr
augustus_ec_dna_fnr

augustus_hs_dna_f1 = (2 * sum_augustus_hs_dna$TP) / (2 * sum_augustus_hs_dna$TP + sum_augustus_hs_dna$FP + sum_augustus_hs_dna$FN)
augustus_hs_dna_precision = sum_augustus_hs_dna$TP / (sum_augustus_hs_dna$TP + sum_augustus_hs_dna$FP)
augustus_hs_dna_sensitivity = sum_augustus_hs_dna$TP / (sum_augustus_hs_dna$TP + sum_augustus_hs_dna$FN)
augustus_hs_dna_fdr = sum_augustus_hs_dna$FP / (sum_augustus_hs_dna$TP + sum_augustus_hs_dna$FP)
augustus_hs_dna_fnr = sum_augustus_hs_dna$FN / (sum_augustus_hs_dna$TP + sum_augustus_hs_dna$FN)
augustus_hs_dna_f1
augustus_hs_dna_precision
augustus_hs_dna_sensitivity
augustus_hs_dna_fdr
augustus_hs_dna_fnr

# RNA viruses

sum_prodigal_rna <- TPtabula_prodigal %>%
  filter(Nucleic_acid == "RNA") %>%
  summarize_if(is.numeric, sum, na.rm=TRUE)
sum_metaprodigal_rna <- TPtabula_metaprodigal %>%
  filter(Nucleic_acid == "RNA") %>%
  summarize_if(is.numeric, sum, na.rm=TRUE)
sum_glimmer_rna <- TPtabula_glimmer %>%
  filter(Nucleic_acid == "RNA") %>%
  summarize_if(is.numeric, sum, na.rm=TRUE)
sum_genemarks_rna <- TPtabula_genemarks %>%
  filter(Nucleic_acid == "RNA") %>%
  summarize_if(is.numeric, sum, na.rm=TRUE)
sum_phanotate_rna <- TPtabula_phanotate %>%
  filter(Nucleic_acid == "RNA") %>%
  summarize_if(is.numeric, sum, na.rm=TRUE)
sum_fraggenescan_rna <- TPtabula_fraggenescan %>%
  filter(Nucleic_acid == "RNA") %>%
  summarize_if(is.numeric, sum, na.rm=TRUE)
sum_mga_rna <- TPtabula_mga %>%
  filter(Nucleic_acid == "RNA") %>%
  summarize_if(is.numeric, sum, na.rm=TRUE)
sum_augustus_sa_rna <- TPtabula_augustus_sa %>%
  filter(Nucleic_acid == "RNA") %>%
  summarize_if(is.numeric, sum, na.rm=TRUE)
sum_augustus_ec_rna <- TPtabula_augustus_ec %>%
  filter(Nucleic_acid == "RNA") %>%
  summarize_if(is.numeric, sum, na.rm=TRUE)
sum_augustus_hs_rna <- TPtabula_augustus_hs %>%
  filter(Nucleic_acid == "RNA") %>%
  summarize_if(is.numeric, sum, na.rm=TRUE)

prodigal_rna_f1 = (2 * sum_prodigal_rna$TP) / (2 * sum_prodigal_rna$TP + sum_prodigal_rna$FP + sum_prodigal_rna$FN)
prodigal_rna_precision = sum_prodigal_rna$TP / (sum_prodigal_rna$TP + sum_prodigal_rna$FP)
prodigal_rna_sensitivity = sum_prodigal_rna$TP / (sum_prodigal_rna$TP + sum_prodigal_rna$FN)
prodigal_rna_fdr = sum_prodigal_rna$FP / (sum_prodigal_rna$TP + sum_prodigal_rna$FP)
prodigal_rna_fnr = sum_prodigal_rna$FN / (sum_prodigal_rna$TP + sum_prodigal_rna$FN)
prodigal_rna_f1
prodigal_rna_precision
prodigal_rna_sensitivity
prodigal_rna_fdr
prodigal_rna_fnr

metaprodigal_rna_f1 = (2 * sum_metaprodigal_rna$TP) / (2 * sum_metaprodigal_rna$TP + sum_metaprodigal_rna$FP + sum_metaprodigal_rna$FN)
metaprodigal_rna_precision = sum_metaprodigal_rna$TP / (sum_metaprodigal_rna$TP + sum_metaprodigal_rna$FP)
metaprodigal_rna_sensitivity = sum_metaprodigal_rna$TP / (sum_metaprodigal_rna$TP + sum_metaprodigal_rna$FN)
metaprodigal_rna_fdr = sum_metaprodigal_rna$FP / (sum_metaprodigal_rna$TP + sum_metaprodigal_rna$FP)
metaprodigal_rna_fnr = sum_metaprodigal_rna$FN / (sum_metaprodigal_rna$TP + sum_metaprodigal_rna$FN)
metaprodigal_rna_f1
metaprodigal_rna_precision
metaprodigal_rna_sensitivity
metaprodigal_rna_fdr
metaprodigal_rna_fnr

glimmer_rna_f1 = (2 * sum_glimmer_rna$TP) / (2 * sum_glimmer_rna$TP + sum_glimmer_rna$FP + sum_glimmer_rna$FN)
glimmer_rna_precision = sum_glimmer_rna$TP / (sum_glimmer_rna$TP + sum_glimmer_rna$FP)
glimmer_rna_sensitivity = sum_glimmer_rna$TP / (sum_glimmer_rna$TP + sum_glimmer_rna$FN)
glimmer_rna_fdr = sum_glimmer_rna$FP / (sum_glimmer_rna$TP + sum_glimmer_rna$FP)
glimmer_rna_fnr = sum_glimmer_rna$FN / (sum_glimmer_rna$TP + sum_glimmer_rna$FN)
glimmer_rna_f1
glimmer_rna_precision
glimmer_rna_sensitivity
glimmer_rna_fdr
glimmer_rna_fnr

genemarks_rna_f1 = (2 * sum_genemarks_rna$TP) / (2 * sum_genemarks_rna$TP + sum_genemarks_rna$FP + sum_genemarks_rna$FN)
genemarks_rna_precision = sum_genemarks_rna$TP / (sum_genemarks_rna$TP + sum_genemarks_rna$FP)
genemarks_rna_sensitivity = sum_genemarks_rna$TP / (sum_genemarks_rna$TP + sum_genemarks_rna$FN)
genemarks_rna_fdr = sum_genemarks_rna$FP / (sum_genemarks_rna$TP + sum_genemarks_rna$FP)
genemarks_rna_fnr = sum_genemarks_rna$FN / (sum_genemarks_rna$TP + sum_genemarks_rna$FN)
genemarks_rna_f1
genemarks_rna_precision
genemarks_rna_sensitivity
genemarks_rna_fdr
genemarks_rna_fnr

phanotate_rna_f1 = (2 * sum_phanotate_rna$TP) / (2 * sum_phanotate_rna$TP + sum_phanotate_rna$FP + sum_phanotate_rna$FN)
phanotate_rna_precision = sum_phanotate_rna$TP / (sum_phanotate_rna$TP + sum_phanotate_rna$FP)
phanotate_rna_sensitivity = sum_phanotate_rna$TP / (sum_phanotate_rna$TP + sum_phanotate_rna$FN)
phanotate_rna_fdr = sum_phanotate_rna$FP / (sum_phanotate_rna$TP + sum_phanotate_rna$FP)
phanotate_rna_fnr = sum_phanotate_rna$FN / (sum_phanotate_rna$TP + sum_phanotate_rna$FN)
phanotate_rna_f1
phanotate_rna_precision
phanotate_rna_sensitivity
phanotate_rna_fdr
phanotate_rna_fnr

fraggenescan_rna_f1 = (2 * sum_fraggenescan_rna$TP) / (2 * sum_fraggenescan_rna$TP + sum_fraggenescan_rna$FP + sum_fraggenescan_rna$FN)
fraggenescan_rna_precision = sum_fraggenescan_rna$TP / (sum_fraggenescan_rna$TP + sum_fraggenescan_rna$FP)
fraggenescan_rna_sensitivity = sum_fraggenescan_rna$TP / (sum_fraggenescan_rna$TP + sum_fraggenescan_rna$FN)
fraggenescan_rna_fdr = sum_fraggenescan_rna$FP / (sum_fraggenescan_rna$TP + sum_fraggenescan_rna$FP)
fraggenescan_rna_fnr = sum_fraggenescan_rna$FN / (sum_fraggenescan_rna$TP + sum_fraggenescan_rna$FN)
fraggenescan_rna_f1
fraggenescan_rna_precision
fraggenescan_rna_sensitivity
fraggenescan_rna_fdr
fraggenescan_rna_fnr

mga_rna_f1 = (2 * sum_mga_rna$TP) / (2 * sum_mga_rna$TP + sum_mga_rna$FP + sum_mga_rna$FN)
mga_rna_precision = sum_mga_rna$TP / (sum_mga_rna$TP + sum_mga_rna$FP)
mga_rna_sensitivity = sum_mga_rna$TP / (sum_mga_rna$TP + sum_mga_rna$FN)
mga_rna_fdr = sum_mga_rna$FP / (sum_mga_rna$TP + sum_mga_rna$FP)
mga_rna_fnr = sum_mga_rna$FN / (sum_mga_rna$TP + sum_mga_rna$FN)
mga_rna_f1
mga_rna_precision
mga_rna_sensitivity
mga_rna_fdr
mga_rna_fnr

augustus_sa_rna_f1 = (2 * sum_augustus_sa_rna$TP) / (2 * sum_augustus_sa_rna$TP + sum_augustus_sa_rna$FP + sum_augustus_sa_rna$FN)
augustus_sa_rna_precision = sum_augustus_sa_rna$TP / (sum_augustus_sa_rna$TP + sum_augustus_sa_rna$FP)
augustus_sa_rna_sensitivity = sum_augustus_sa_rna$TP / (sum_augustus_sa_rna$TP + sum_augustus_sa_rna$FN)
augustus_sa_rna_fdr = sum_augustus_sa_rna$FP / (sum_augustus_sa_rna$TP + sum_augustus_sa_rna$FP)
augustus_sa_rna_fnr = sum_augustus_sa_rna$FN / (sum_augustus_sa_rna$TP + sum_augustus_sa_rna$FN)
augustus_sa_rna_f1
augustus_sa_rna_precision
augustus_sa_rna_sensitivity
augustus_sa_rna_fdr
augustus_sa_rna_fnr

augustus_ec_rna_f1 = (2 * sum_augustus_ec_rna$TP) / (2 * sum_augustus_ec_rna$TP + sum_augustus_ec_rna$FP + sum_augustus_ec_rna$FN)
augustus_ec_rna_precision = sum_augustus_ec_rna$TP / (sum_augustus_ec_rna$TP + sum_augustus_ec_rna$FP)
augustus_ec_rna_sensitivity = sum_augustus_ec_rna$TP / (sum_augustus_ec_rna$TP + sum_augustus_ec_rna$FN)
augustus_ec_rna_fdr = sum_augustus_ec_rna$FP / (sum_augustus_ec_rna$TP + sum_augustus_ec_rna$FP)
augustus_ec_rna_fnr = sum_augustus_ec_rna$FN / (sum_augustus_ec_rna$TP + sum_augustus_ec_rna$FN)
augustus_ec_rna_f1
augustus_ec_rna_precision
augustus_ec_rna_sensitivity
augustus_ec_rna_fdr
augustus_ec_rna_fnr

augustus_hs_rna_f1 = (2 * sum_augustus_hs_rna$TP) / (2 * sum_augustus_hs_rna$TP + sum_augustus_hs_rna$FP + sum_augustus_hs_rna$FN)
augustus_hs_rna_precision = sum_augustus_hs_rna$TP / (sum_augustus_hs_rna$TP + sum_augustus_hs_rna$FP)
augustus_hs_rna_sensitivity = sum_augustus_hs_rna$TP / (sum_augustus_hs_rna$TP + sum_augustus_hs_rna$FN)
augustus_hs_rna_fdr = sum_augustus_hs_rna$FP / (sum_augustus_hs_rna$TP + sum_augustus_hs_rna$FP)
augustus_hs_rna_fnr = sum_augustus_hs_rna$FN / (sum_augustus_hs_rna$TP + sum_augustus_hs_rna$FN)
augustus_hs_rna_f1
augustus_hs_rna_precision
augustus_hs_rna_sensitivity
augustus_hs_rna_fdr
augustus_hs_rna_fnr

# Graphs
all_precision <- c(prodigal_precision, metaprodigal_precision, glimmer_precision, genemarks_precision, 
                   phanotate_precision, fraggenescan_precision, mga_precision, augustus_sa_precision, augustus_ec_precision, 
                   augustus_hs_precision)
all_sensitivity <- c(prodigal_sensitivity, metaprodigal_sensitivity, glimmer_sensitivity, genemarks_sensitivity, 
                     phanotate_sensitivity, fraggenescan_sensitivity, mga_sensitivity, augustus_sa_sensitivity, augustus_ec_sensitivity, 
                     augustus_hs_sensitivity)
all_programs <- c("Prodigal", "Metaprodigal", "Glimmer", "GeneMarkS", "PHANOTATE", "FragGeneScan", "MGA", "AUGUSTUS_SA", 
                  "AUGUSTUS_EC", "AUGUSTUS_HS")
maintablestats <- tibble(programs = all_programs, precision = all_precision, sensitivity = all_sensitivity)
genplot2 <- ggplot(maintablestats, aes(x= precision, y= sensitivity, colour=as.factor(programs), label=programs)) +
  geom_text(aes(label=programs),hjust=0, vjust=0, size=5) + xlim(0,1) + ylim (0,1) +
  scale_color_manual(values=c('lightpink', 'goldenrod', 'purple', 'blue', 'springgreen4', 'chartreuse3', 'orange', 'deeppink', 'cyan3', 'firebrick'))

all_arc_precision <- c(prodigal_arc_precision, metaprodigal_arc_precision, glimmer_arc_precision, genemarks_arc_precision, 
                       phanotate_arc_precision, fraggenescan_arc_precision, mga_arc_precision, augustus_sa_arc_precision, 
                       augustus_ec_arc_precision, augustus_hs_arc_precision)
all_arc_sensitivity <- c(prodigal_arc_sensitivity, metaprodigal_arc_sensitivity, glimmer_arc_sensitivity, genemarks_arc_sensitivity, 
                         phanotate_arc_sensitivity, fraggenescan_arc_sensitivity, mga_arc_sensitivity, augustus_sa_arc_sensitivity, 
                         augustus_ec_arc_sensitivity, augustus_hs_arc_sensitivity)
all_arc_programs <- c("Prodigal", "Metaprodigal", "Glimmer", "GeneMarkS", "PHANOTATE", "FragGeneScan", "MGA", "AUGUSTUS_SA", 
                      "AUGUSTUS_EC", "AUGUSTUS_HS")
maintablestats_arc <- tibble(programs = all_arc_programs, precision = all_arc_precision, sensitivity = all_arc_sensitivity)
arcplot2 <- ggplot(maintablestats_arc, aes(x= precision, y= sensitivity, colour=as.factor(programs), label=programs)) +
  geom_text(aes(label=programs),hjust=0, vjust=0, size=5) + xlim(0,1) + ylim (0,1) +
  scale_color_manual(values=c('lightpink', 'goldenrod', 'purple', 'blue', 'springgreen4', 'chartreuse3', 'orange', 'deeppink', 'cyan3', 'firebrick'))

all_bac_precision <- c(prodigal_bac_precision, metaprodigal_bac_precision, glimmer_bac_precision, genemarks_bac_precision, 
                       phanotate_bac_precision, fraggenescan_bac_precision, mga_bac_precision, augustus_sa_bac_precision, 
                       augustus_ec_bac_precision, augustus_hs_bac_precision)
all_bac_sensitivity <- c(prodigal_bac_sensitivity, metaprodigal_bac_sensitivity, glimmer_bac_sensitivity, genemarks_bac_sensitivity, 
                         phanotate_bac_sensitivity, fraggenescan_bac_sensitivity, mga_bac_sensitivity, augustus_sa_bac_sensitivity, 
                         augustus_ec_bac_sensitivity, augustus_hs_bac_sensitivity)
all_bac_programs <- c("Prodigal", "Metaprodigal", "Glimmer", "GeneMarkS", "PHANOTATE", "FragGeneScan", "MGA", "AUGUSTUS_SA", 
                      "AUGUSTUS_EC", "AUGUSTUS_HS")
maintablestats_bac <- tibble(programs = all_bac_programs, precision = all_bac_precision, sensitivity = all_bac_sensitivity)
bacplot2 <- ggplot(maintablestats_bac, aes(x= precision, y= sensitivity, colour=as.factor(programs), label=programs)) +
  geom_text(aes(label=programs),hjust=0, vjust=0, size=5) + xlim(0,1) + ylim (0,1) +
  scale_color_manual(values=c('lightpink', 'goldenrod', 'purple', 'blue', 'springgreen4', 'chartreuse3', 'orange', 'deeppink', 'cyan3', 'firebrick'))

all_euk_precision <- c(prodigal_euk_precision, metaprodigal_euk_precision, glimmer_euk_precision, genemarks_euk_precision, 
                       phanotate_euk_precision, fraggenescan_euk_precision, mga_euk_precision, augustus_sa_euk_precision, 
                       augustus_ec_euk_precision, augustus_hs_euk_precision)
all_euk_sensitivity <- c(prodigal_euk_sensitivity, metaprodigal_euk_sensitivity, glimmer_euk_sensitivity, genemarks_euk_sensitivity, 
                         phanotate_euk_sensitivity, fraggenescan_euk_sensitivity, mga_euk_sensitivity, augustus_sa_euk_sensitivity, 
                         augustus_ec_euk_sensitivity, augustus_hs_euk_sensitivity)
all_euk_programs <- c("Prodigal", "Metaprodigal", "Glimmer", "GeneMarkS", "PHANOTATE", "FragGeneScan", "MGA", "AUGUSTUS_SA", 
                      "AUGUSTUS_EC", "AUGUSTUS_HS")
maintablestats_euk <- tibble(programs = all_euk_programs, precision = all_euk_precision, sensitivity = all_euk_sensitivity)
eukplot2 <- ggplot(maintablestats_euk, aes(x= precision, y= sensitivity, colour=as.factor(programs), label=programs)) +
  geom_text(aes(label=programs),hjust=0, vjust=0, size=5) + xlim(0,1) + ylim (0,1) +
  scale_color_manual(values=c('lightpink', 'goldenrod', 'purple', 'blue', 'springgreen4', 'chartreuse3', 'orange', 'deeppink', 'cyan3', 'firebrick'))

all_dna_precision <- c(prodigal_dna_precision, metaprodigal_dna_precision, glimmer_dna_precision, genemarks_dna_precision, 
                       phanotate_dna_precision, fraggenescan_dna_precision, mga_dna_precision, augustus_sa_dna_precision, 
                       augustus_ec_dna_precision, augustus_hs_dna_precision)
all_dna_sensitivity <- c(prodigal_dna_sensitivity, metaprodigal_dna_sensitivity, glimmer_dna_sensitivity, genemarks_dna_sensitivity, 
                         phanotate_dna_sensitivity, fraggenescan_dna_sensitivity, mga_dna_sensitivity, augustus_sa_dna_sensitivity, 
                         augustus_ec_dna_sensitivity, augustus_hs_dna_sensitivity)
all_dna_programs <- c("Prodigal", "Metaprodigal", "Glimmer", "GeneMarkS", "PHANOTATE", "FragGeneScan", "MGA", "AUGUSTUS_SA", 
                      "AUGUSTUS_EC", "AUGUSTUS_HS")
maintablestats_dna <- tibble(programs = all_dna_programs, precision = all_dna_precision, sensitivity = all_dna_sensitivity)
dnaplot2 <- ggplot(maintablestats_dna, aes(x= precision, y= sensitivity, colour=as.factor(programs), label=programs)) +
  geom_text(aes(label=programs),hjust=0, vjust=0, size=5) + xlim(0,1) + ylim (0,1) +
  scale_color_manual(values=c('lightpink', 'goldenrod', 'purple', 'blue', 'springgreen4', 'chartreuse3', 'orange', 'deeppink', 'cyan3', 'firebrick'))

all_rna_precision <- c(prodigal_rna_precision, metaprodigal_rna_precision, glimmer_rna_precision, genemarks_rna_precision, 
                       phanotate_rna_precision, fraggenescan_rna_precision, mga_rna_precision, augustus_sa_rna_precision, 
                       augustus_ec_rna_precision, augustus_hs_rna_precision)
all_rna_sensitivity <- c(prodigal_rna_sensitivity, metaprodigal_rna_sensitivity, glimmer_rna_sensitivity, genemarks_rna_sensitivity, 
                         phanotate_rna_sensitivity, fraggenescan_rna_sensitivity, mga_rna_sensitivity, augustus_sa_rna_sensitivity, 
                         augustus_ec_rna_sensitivity, augustus_hs_rna_sensitivity)
all_rna_programs <- c("Prodigal", "Metaprodigal", "Glimmer", "GeneMarkS", "PHANOTATE", "FragGeneScan", "MGA", "AUGUSTUS_SA", 
                      "AUGUSTUS_EC", "AUGUSTUS_HS")
maintablestats_rna <- tibble(programs = all_rna_programs, precision = all_rna_precision, sensitivity = all_rna_sensitivity)
rnaplot2 <- ggplot(maintablestats_rna, aes(x= precision, y= sensitivity, colour=as.factor(programs), label=programs)) +
  geom_text(aes(label=programs),hjust=0, vjust=0, size=5) + xlim(0,1) + ylim (0,1) +
  scale_color_manual(values=c('lightpink', 'goldenrod', 'purple', 'blue', 'springgreen4', 'chartreuse3', 'orange', 'deeppink', 'cyan3', 'firebrick'))

legend2 <- get_legend(genplot2)
svg("figure_2.svg", width = 10, height = 10)
plot_grid(genplot2 + theme(legend.position = "null"), arcplot2 + theme(legend.position = "null"), bacplot2 + theme(legend.position = "null"), 
                     eukplot2 + theme(legend.position = "null"), dnaplot2 + theme(legend.position = "null"), rnaplot2 + theme(legend.position = "null"), labels = c("A", "B", "C", "D", "E", "F"), ncol = 2)
#plot_grid(topgrid2, legend2, ncol = 1, rel_heights = c(3,.25))
dev.off()

