library(ggplot2)
library(tidyr)
library(dplyr)
library(cowplot)
library(emmeans)
library(FSA)

## FIRST TEST: No of ORFs

tablaoriginal <- read.csv("viral_genomes.ORIGINAL.tsv", header = T, sep = "\t")
tablaprodigal <- read.csv("viral_genomes.PRODIGAL.tsv", header = T, sep = "\t")
tablametaprodigal <- read.csv("viral_genomes.METAPRODIGAL.tsv", header = T, sep = "\t")
tablaglimmer <- read.csv("viral_genomes.GLIMMER.tsv", header = T, sep = "\t")
tablagenemarks <- read.csv("viral_genomes.GENEMARKS.tsv", header = T, sep = "\t")
tablaphanotate <- read.csv("viral_genomes.PHANOTATE.tsv", header = T, sep = "\t")
tablafraggenescan <- read.csv("viral_genomes.FRAGGENESCAN.tsv", header = T, sep = "\t")
tablamga <- read.csv("viral_genomes.MGA.tsv", header = T, sep = "\t")
tablaaugustussa <- read.csv("viral_genomes.AUGUSTUS_SAUREUS.tsv", header = T, sep = "\t")
tablaaugustusec <- read.csv("viral_genomes.AUGUSTUS_ECOLI.tsv", header = T, sep = "\t")
tablaaugustushs <- read.csv("viral_genomes.AUGUSTUS_HUMAN.tsv", header = T, sep = "\t")
tablahost <- read.csv("list_organisms.tsv", header = T, sep = "\t")
tablanucleic <- read.csv("list_viralnucleic_2.csv", header = T, sep = "\t")

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
pairs(emtrends(generallm, ~Program, var="No_CDS"))

# Eukaryotic viruses

tabula_euk <- tabula %>% filter(Host_domain == "Eukaryote")
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

tabula_euk.long.2 <- tabula_euk %>%
  pivot_longer(cols=c(5:14), names_to = "Program", values_to = "No_CDS")
generallm_euk <- lm(CDSs_EXPECTED ~ No_CDS * Program + 0, data = tabula_euk.long.2)
anova(generallm_euk)
summary(generallm_euk)
pairs(emtrends(generallm_euk, ~Program, var="No_CDS"))


# Archaeal viruses

tabula_arc <- tabula %>% filter(Host_domain == "Archaea")
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

tabula_arc.long.2 <- tabula_arc %>%
  pivot_longer(cols=c(5:14), names_to = "Program", values_to = "No_CDS")
generallm_arc <- lm(CDSs_EXPECTED ~ No_CDS * Program + 0, data = tabula_arc.long.2)
anova(generallm_arc)
summary(generallm_arc)
pairs(emtrends(generallm_arc, ~Program, var="No_CDS"))


# Bacteriophages

tabula_bac <- tabula %>% filter(Host_domain == "Bacteria")
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

tabula_bac.long.2 <- tabula_bac %>%
  pivot_longer(cols=c(5:14), names_to = "Program", values_to = "No_CDS")
generallm_bac <- lm(CDSs_EXPECTED ~ No_CDS * Program + 0, data = tabula_bac.long.2)
anova(generallm_bac)
summary(generallm_bac)
pairs(emtrends(generallm_bac, ~Program, var="No_CDS"))


# DNA viruses

tabula_dna <- tabula %>% filter(Nucleic_acid == "DNA")
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

tabula_dna.long.2 <- tabula_dna %>%
  pivot_longer(cols=c(5:14), names_to = "Program", values_to = "No_CDS")
generallm_dna <- lm(CDSs_EXPECTED ~ No_CDS * Program + 0, data = tabula_dna.long.2)
anova(generallm_dna)
summary(generallm_dna)
pairs(emtrends(generallm_dna, ~Program, var="No_CDS"))


# RNA viruses

tabula_rna <- tabula %>% filter(Nucleic_acid == "RNA")
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

tabula_rna.long.2 <- tabula_rna %>%
  pivot_longer(cols=c(5:14), names_to = "Program", values_to = "No_CDS")
generallm_rna <- lm(CDSs_EXPECTED ~ No_CDS * Program + 0, data = tabula_rna.long.2)
anova(generallm_rna)
summary(generallm_rna)
pairs(emtrends(generallm_rna, ~Program, var="No_CDS"))

## SECOND TEST: Coordinates

Basal_table <- read.csv("viral_genomes.PRODIGAL.TPtable.tsv", header = T, sep = "\t") %>% 
  rename(TP_Prodigal = TP, FP_Prodigal = FP, FN_Prodigal = FN) %>%
  full_join(read.csv("viral_genomes.METAPRODIGAL.TPtable.tsv", header = T, sep = "\t"), by = c("Genome" = "Genome")) %>%
  rename(TP_MetaProdigal = TP, FP_MetaProdigal = FP, FN_MetaProdigal = FN) %>%
  full_join(read.csv("viral_genomes.GLIMMER.TPtable.tsv", header = T, sep = "\t"), by = c("Genome" = "Genome")) %>%
  rename(TP_Glimmer = TP, FP_Glimmer = FP, FN_Glimmer = FN) %>%
  full_join(read.csv("viral_genomes.GENEMARKS.TPtable.tsv", header = T, sep = "\t"), by = c("Genome" = "Genome")) %>%
  rename(TP_GeneMarkS = TP, FP_GeneMarkS = FP, FN_GeneMarkS = FN) %>%
  full_join(read.csv("viral_genomes.PHANOTATE.TPtable.tsv", header = T, sep = "\t"), by = c("Genome" = "Genome")) %>%
  rename(TP_Phanotate = TP, FP_Phanotate = FP, FN_Phanotate = FN) %>%
  full_join(read.csv("viral_genomes.FRAGGENESCAN.TPtable.tsv", header = T, sep = "\t"), by = c("Genome" = "Genome")) %>%
  rename(TP_FragGeneScan = TP, FP_FragGeneScan = FP, FN_FragGeneScan = FN) %>%
  full_join(read.csv("viral_genomes.MGA.TPtable.tsv", header = T, sep = "\t"), by = c("Genome" = "Genome")) %>%
  rename(TP_MGA = TP, FP_MGA = FP, FN_MGA = FN) %>%
  full_join(read.csv("viral_genomes.AUGUSTUS_SAUREUS.TPtable.tsv", header = T, sep = "\t"), by = c("Genome" = "Genome")) %>%
  rename(TP_AUGUSTUS_SA = TP, FP_AUGUSTUS_SA = FP, FN_AUGUSTUS_SA = FN) %>%
  full_join(read.csv("viral_genomes.AUGUSTUS_ECOLI.TPtable.tsv", header = T, sep = "\t"), by = c("Genome" = "Genome")) %>%
  rename(TP_AUGUSTUS_EC = TP, FP_AUGUSTUS_EC = FP, FN_AUGUSTUS_EC = FN) %>%
  full_join(read.csv("viral_genomes.AUGUSTUS_HUMAN.TPtable.tsv", header = T, sep = "\t"), by = c("Genome" = "Genome")) %>%
  rename(TP_AUGUSTUS_HS = TP, FP_AUGUSTUS_HS = FP, FN_AUGUSTUS_HS = FN) %>%
  full_join(tablahost, by = c("Genome" = "ID")) %>%
  full_join(tablanucleic, by = c("Genome" = "ID"))

# General (Real, absolute values)

sum_all <-Basal_table %>%
  summarize_if(is.numeric, sum, na.rm=TRUE) %>%
  rowwise() %>%
  mutate(Prodigal_F1 = (2 * TP_Prodigal) / (2 * TP_Prodigal + FP_Prodigal + FN_Prodigal),
         Prodigal_Precision = TP_Prodigal / (TP_Prodigal + FP_Prodigal),
         Prodigal_Sensitivity = TP_Prodigal / (TP_Prodigal + FN_Prodigal),
         Prodigal_FDR = FP_Prodigal / (TP_Prodigal + FP_Prodigal),
         Prodigal_FNR = FN_Prodigal / (TP_Prodigal + FN_Prodigal),
         MetaProdigal_F1 = (2 * TP_MetaProdigal) / (2 * TP_MetaProdigal + FP_MetaProdigal + FN_MetaProdigal),
         MetaProdigal_Precision = TP_MetaProdigal / (TP_MetaProdigal + FP_MetaProdigal),
         MetaProdigal_Sensitivity = TP_MetaProdigal / (TP_MetaProdigal + FN_MetaProdigal),
         MetaProdigal_FDR = FP_MetaProdigal / (TP_MetaProdigal + FP_MetaProdigal),
         MetaProdigal_FNR = FN_MetaProdigal / (TP_MetaProdigal + FN_MetaProdigal),
         Glimmer_F1 = (2 * TP_Glimmer) / (2 * TP_Glimmer + FP_Glimmer + FN_Glimmer),
         Glimmer_Precision = TP_Glimmer / (TP_Glimmer + FP_Glimmer),
         Glimmer_Sensitivity = TP_Glimmer / (TP_Glimmer + FN_Glimmer),
         Glimmer_FDR = FP_Glimmer / (TP_Glimmer + FP_Glimmer),
         Glimmer_FNR = FN_Glimmer / (TP_Glimmer + FN_Glimmer),
         GeneMarkS_F1 = (2 * TP_GeneMarkS) / (2 * TP_GeneMarkS + FP_GeneMarkS + FN_GeneMarkS),
         GeneMarkS_Precision = TP_GeneMarkS / (TP_GeneMarkS + FP_GeneMarkS),
         GeneMarkS_Sensitivity = TP_GeneMarkS / (TP_GeneMarkS + FN_GeneMarkS),
         GeneMarkS_FDR = FP_GeneMarkS / (TP_GeneMarkS + FP_GeneMarkS),
         GeneMarkS_FNR = FN_GeneMarkS / (TP_GeneMarkS + FN_GeneMarkS),
         Phanotate_F1 = (2 * TP_Phanotate) / (2 * TP_Phanotate + FP_Phanotate + FN_Phanotate),
         Phanotate_Precision = TP_Phanotate / (TP_Phanotate + FP_Phanotate),
         Phanotate_Sensitivity = TP_Phanotate / (TP_Phanotate + FN_Phanotate),
         Phanotate_FDR = FP_Phanotate / (TP_Phanotate + FP_Phanotate),
         Phanotate_FNR = FN_Phanotate / (TP_Phanotate + FN_Phanotate),
         FragGeneScan_F1 = (2 * TP_FragGeneScan) / (2 * TP_FragGeneScan + FP_FragGeneScan + FN_FragGeneScan),
         FragGeneScan_Precision = TP_FragGeneScan / (TP_FragGeneScan + FP_FragGeneScan),
         FragGeneScan_Sensitivity = TP_FragGeneScan / (TP_FragGeneScan + FN_FragGeneScan),
         FragGeneScan_FDR = FP_FragGeneScan / (TP_FragGeneScan + FP_FragGeneScan),
         FragGeneScan_FNR = FN_FragGeneScan / (TP_FragGeneScan + FN_FragGeneScan),
         MGA_F1 = (2 * TP_MGA) / (2 * TP_MGA + FP_MGA + FN_MGA),
         MGA_Precision = TP_MGA / (TP_MGA + FP_MGA),
         MGA_Sensitivity = TP_MGA / (TP_MGA + FN_MGA),
         MGA_FDR = FP_MGA / (TP_MGA + FP_MGA),
         MGA_FNR = FN_MGA / (TP_MGA + FN_MGA),
         AUGUSTUS_SA_F1 = (2 * TP_AUGUSTUS_SA) / (2 * TP_AUGUSTUS_SA + FP_AUGUSTUS_SA + FN_AUGUSTUS_SA),
         AUGUSTUS_SA_Precision = TP_AUGUSTUS_SA / (TP_AUGUSTUS_SA + FP_AUGUSTUS_SA),
         AUGUSTUS_SA_Sensitivity = TP_AUGUSTUS_SA / (TP_AUGUSTUS_SA + FN_AUGUSTUS_SA),
         AUGUSTUS_SA_FDR = FP_AUGUSTUS_SA / (TP_AUGUSTUS_SA + FP_AUGUSTUS_SA),
         AUGUSTUS_SA_FNR = FN_AUGUSTUS_SA / (TP_AUGUSTUS_SA + FN_AUGUSTUS_SA),
         AUGUSTUS_EC_F1 = (2 * TP_AUGUSTUS_EC) / (2 * TP_AUGUSTUS_EC + FP_AUGUSTUS_EC + FN_AUGUSTUS_EC),
         AUGUSTUS_EC_Precision = TP_AUGUSTUS_EC / (TP_AUGUSTUS_EC + FP_AUGUSTUS_EC),
         AUGUSTUS_EC_Sensitivity = TP_AUGUSTUS_EC / (TP_AUGUSTUS_EC + FN_AUGUSTUS_EC),
         AUGUSTUS_EC_FDR = FP_AUGUSTUS_EC / (TP_AUGUSTUS_EC + FP_AUGUSTUS_EC),
         AUGUSTUS_EC_FNR = FN_AUGUSTUS_EC / (TP_AUGUSTUS_EC + FN_AUGUSTUS_EC),
         AUGUSTUS_HS_F1 = (2 * TP_AUGUSTUS_HS) / (2 * TP_AUGUSTUS_HS + FP_AUGUSTUS_HS + FN_AUGUSTUS_HS),
         AUGUSTUS_HS_Precision = TP_AUGUSTUS_HS / (TP_AUGUSTUS_HS + FP_AUGUSTUS_HS),
         AUGUSTUS_HS_Sensitivity = TP_AUGUSTUS_HS / (TP_AUGUSTUS_HS + FN_AUGUSTUS_HS),
         AUGUSTUS_HS_FDR = FP_AUGUSTUS_HS / (TP_AUGUSTUS_HS + FP_AUGUSTUS_HS),
         AUGUSTUS_HS_FNR = FN_AUGUSTUS_HS / (TP_AUGUSTUS_HS + FN_AUGUSTUS_HS))

# General (Subsapling)

listofdfs = list()
for(i in 1:1000){
  listofdfs[[i]] <- Basal_table %>%
    slice_sample(prop = 0.1)
}
dfofdfs <- tibble(d = 1:1000, data = listofdfs)

sum_all_subsamples <- list()
for(i in 1:1000){
  sum_all_subsamples[[i]] <- dfofdfs$data[[i]] %>%
    summarize_if(is.numeric, sum, na.rm=TRUE)
}

subsamples_sum_all <- bind_rows(sum_all_subsamples, .id = "subsample_number") %>%
  rowwise() %>%
  mutate(Prodigal_F1 = (2 * TP_Prodigal) / (2 * TP_Prodigal + FP_Prodigal + FN_Prodigal),
         Prodigal_Precision = TP_Prodigal / (TP_Prodigal + FP_Prodigal),
         Prodigal_Sensitivity = TP_Prodigal / (TP_Prodigal + FN_Prodigal),
         Prodigal_FDR = FP_Prodigal / (TP_Prodigal + FP_Prodigal),
         Prodigal_FNR = FN_Prodigal / (TP_Prodigal + FN_Prodigal),
         MetaProdigal_F1 = (2 * TP_MetaProdigal) / (2 * TP_MetaProdigal + FP_MetaProdigal + FN_MetaProdigal),
         MetaProdigal_Precision = TP_MetaProdigal / (TP_MetaProdigal + FP_MetaProdigal),
         MetaProdigal_Sensitivity = TP_MetaProdigal / (TP_MetaProdigal + FN_MetaProdigal),
         MetaProdigal_FDR = FP_MetaProdigal / (TP_MetaProdigal + FP_MetaProdigal),
         MetaProdigal_FNR = FN_MetaProdigal / (TP_MetaProdigal + FN_MetaProdigal),
         Glimmer_F1 = (2 * TP_Glimmer) / (2 * TP_Glimmer + FP_Glimmer + FN_Glimmer),
         Glimmer_Precision = TP_Glimmer / (TP_Glimmer + FP_Glimmer),
         Glimmer_Sensitivity = TP_Glimmer / (TP_Glimmer + FN_Glimmer),
         Glimmer_FDR = FP_Glimmer / (TP_Glimmer + FP_Glimmer),
         Glimmer_FNR = FN_Glimmer / (TP_Glimmer + FN_Glimmer),
         GeneMarkS_F1 = (2 * TP_GeneMarkS) / (2 * TP_GeneMarkS + FP_GeneMarkS + FN_GeneMarkS),
         GeneMarkS_Precision = TP_GeneMarkS / (TP_GeneMarkS + FP_GeneMarkS),
         GeneMarkS_Sensitivity = TP_GeneMarkS / (TP_GeneMarkS + FN_GeneMarkS),
         GeneMarkS_FDR = FP_GeneMarkS / (TP_GeneMarkS + FP_GeneMarkS),
         GeneMarkS_FNR = FN_GeneMarkS / (TP_GeneMarkS + FN_GeneMarkS),
         Phanotate_F1 = (2 * TP_Phanotate) / (2 * TP_Phanotate + FP_Phanotate + FN_Phanotate),
         Phanotate_Precision = TP_Phanotate / (TP_Phanotate + FP_Phanotate),
         Phanotate_Sensitivity = TP_Phanotate / (TP_Phanotate + FN_Phanotate),
         Phanotate_FDR = FP_Phanotate / (TP_Phanotate + FP_Phanotate),
         Phanotate_FNR = FN_Phanotate / (TP_Phanotate + FN_Phanotate),
         FragGeneScan_F1 = (2 * TP_FragGeneScan) / (2 * TP_FragGeneScan + FP_FragGeneScan + FN_FragGeneScan),
         FragGeneScan_Precision = TP_FragGeneScan / (TP_FragGeneScan + FP_FragGeneScan),
         FragGeneScan_Sensitivity = TP_FragGeneScan / (TP_FragGeneScan + FN_FragGeneScan),
         FragGeneScan_FDR = FP_FragGeneScan / (TP_FragGeneScan + FP_FragGeneScan),
         FragGeneScan_FNR = FN_FragGeneScan / (TP_FragGeneScan + FN_FragGeneScan),
         MGA_F1 = (2 * TP_MGA) / (2 * TP_MGA + FP_MGA + FN_MGA),
         MGA_Precision = TP_MGA / (TP_MGA + FP_MGA),
         MGA_Sensitivity = TP_MGA / (TP_MGA + FN_MGA),
         MGA_FDR = FP_MGA / (TP_MGA + FP_MGA),
         MGA_FNR = FN_MGA / (TP_MGA + FN_MGA),
         AUGUSTUSSA_F1 = (2 * TP_AUGUSTUS_SA) / (2 * TP_AUGUSTUS_SA + FP_AUGUSTUS_SA + FN_AUGUSTUS_SA),
         AUGUSTUSSA_Precision = TP_AUGUSTUS_SA / (TP_AUGUSTUS_SA + FP_AUGUSTUS_SA),
         AUGUSTUSSA_Sensitivity = TP_AUGUSTUS_SA / (TP_AUGUSTUS_SA + FN_AUGUSTUS_SA),
         AUGUSTUSSA_FDR = FP_AUGUSTUS_SA / (TP_AUGUSTUS_SA + FP_AUGUSTUS_SA),
         AUGUSTUSSA_FNR = FN_AUGUSTUS_SA / (TP_AUGUSTUS_SA + FN_AUGUSTUS_SA),
         AUGUSTUSEC_F1 = (2 * TP_AUGUSTUS_EC) / (2 * TP_AUGUSTUS_EC + FP_AUGUSTUS_EC + FN_AUGUSTUS_EC),
         AUGUSTUSEC_Precision = TP_AUGUSTUS_EC / (TP_AUGUSTUS_EC + FP_AUGUSTUS_EC),
         AUGUSTUSEC_Sensitivity = TP_AUGUSTUS_EC / (TP_AUGUSTUS_EC + FN_AUGUSTUS_EC),
         AUGUSTUSEC_FDR = FP_AUGUSTUS_EC / (TP_AUGUSTUS_EC + FP_AUGUSTUS_EC),
         AUGUSTUSEC_FNR = FN_AUGUSTUS_EC / (TP_AUGUSTUS_EC + FN_AUGUSTUS_EC),
         AUGUSTUSHS_F1 = (2 * TP_AUGUSTUS_HS) / (2 * TP_AUGUSTUS_HS + FP_AUGUSTUS_HS + FN_AUGUSTUS_HS),
         AUGUSTUSHS_Precision = TP_AUGUSTUS_HS / (TP_AUGUSTUS_HS + FP_AUGUSTUS_HS),
         AUGUSTUSHS_Sensitivity = TP_AUGUSTUS_HS / (TP_AUGUSTUS_HS + FN_AUGUSTUS_HS),
         AUGUSTUSHS_FDR = FP_AUGUSTUS_HS / (TP_AUGUSTUS_HS + FP_AUGUSTUS_HS),
         AUGUSTUSHS_FNR = FN_AUGUSTUS_HS / (TP_AUGUSTUS_HS + FN_AUGUSTUS_HS))

#Testing normality for the variables F1, Precision, Sensitivity, FDR and FNR of the different programs
listpvalues_shapiro = list()

for(i in 32:length(subsamples_sum_all)){
  listpvalues_shapiro[[i]] <- shapiro.test(subsamples_sum_all[[i]])$p.value
}

listpvalues_shapiro <- Filter(Negate(is.null), listpvalues_shapiro)
listpvalues_shapiro <- as.list(p.adjust(listpvalues_shapiro, method = "fdr"))
datapvalues_shapiro <- tibble(d = 1:length(listpvalues_shapiro),
                              variable = c("Prodigal_F1", "Prodigal_Precision", "Prodigal_Sensitivity", "Prodigal_FDR", "Prodigal_FNR",
                                           "MetaProdigal_F1", "MetaProdigal_Precision", "MetaProdigal_Sensitivity", "MetaProdigal_FDR", "MetaProdigal_FNR",
                                           "Glimmer_F1", "Glimmer_Precision", "Glimmer_Sensitivity", "Glimmer_FDR", "Glimmer_FNR",
                                           "GeneMarkS_F1", "GeneMarkS_Precision", "GeneMarkS_Sensitivity", "GeneMarkS_FDR", "GeneMarkS_FNR",
                                           "Phanotate_F1", "Phanotate_Precision", "Phanotate_Sensitivity", "Phanotate_FDR", "Phanotate_FNR", 
                                           "FragGeneScan_F1", "FragGeneScan_Precision", "FragGeneScan_Sensitivity", "FragGeneScan_FDR", "FragGeneScan_FNR", 
                                           "MGA_F1", "MGA_Precision", "MGA_Sensitivity", "MGA_FDR", "MGA_FNR",
                                           "AUGUSTUS_SA_F1", "AUGUSTUS_SA_Precision", "AUGUSTUS_SA_Sensitivity", "AUGUSTUS_SA_FDR", "AUGUSTUS_SA_FNR",
                                           "AUGUSTUS_EC_F1", "AUGUSTUS_EC_Precision", "AUGUSTUS_EC_Sensitivity", "AUGUSTUS_EC_FDR", "AUGUSTUS_EC_FNR",
                                           "AUGUSTUS_HS_F1", "AUGUSTUS_HS_Precision", "AUGUSTUS_HS_Sensitivity", "AUGUSTUS_HS_FDR", "AUGUSTUS_HS_FNR"),
                              p_value = listpvalues_shapiro) %>%
  select(-d) %>%
  pivot_wider(names_from = variable, values_from = p_value)

#All p-values (or near all) were significant, so there is no chance to assume normality in our data
#Thus, the only way to compare values among them is via Kruskal-Wallis followed by pairwise tests (Dunn test)

#Measuring the differences between F1 values
datastatsrep <- subsamples_sum_all %>%
  select(c(1, 32:length(subsamples_sum_all))) %>%
  gather(key = var_name, value = value, 2:51) %>%
  separate(col = var_name, sep = "_", into = c("Program", "Statistic")) %>%
  select(-subsample_number)

F1_datarep <- datastatsrep %>% filter(Statistic == "F1")
Precision_datarep <- datastatsrep %>% filter(Statistic == "Precision")
Sensitivity_datarep <- datastatsrep %>% filter(Statistic == "Sensitivity")
FDR_datarep <- datastatsrep %>% filter(Statistic == "FDR")
FNR_datarep <- datastatsrep %>% filter(Statistic == "FNR")

kruskal.test(x = F1_datarep$value, g = F1_datarep$Program)
dunnTest(value ~ Program, data = F1_datarep, method = "holm")

kruskal.test(x = Precision_datarep$value, g = Precision_datarep$Program)
dunnTest(value ~ Program, data = Precision_datarep, method = "holm")

kruskal.test(x = Sensitivity_datarep$value, g = Sensitivity_datarep$Program)
dunnTest(value ~ Program, data = Sensitivity_datarep, method = "holm")

kruskal.test(x = FDR_datarep$value, g = FDR_datarep$Program)
dunnTest(value ~ Program, data = FDR_datarep, method = "holm")

kruskal.test(x = FNR_datarep$value, g = FNR_datarep$Program)
dunnTest(value ~ Program, data = FNR_datarep, method = "holm")

# Eukaryote (Real)

sum_euk <- Basal_table %>%
  filter(Host_domain == "Eukaryote") %>%
  summarize_if(is.numeric, sum, na.rm=TRUE) %>%
  rowwise() %>%
  mutate(Prodigal_F1 = (2 * TP_Prodigal) / (2 * TP_Prodigal + FP_Prodigal + FN_Prodigal),
         Prodigal_Precision = TP_Prodigal / (TP_Prodigal + FP_Prodigal),
         Prodigal_Sensitivity = TP_Prodigal / (TP_Prodigal + FN_Prodigal),
         Prodigal_FDR = FP_Prodigal / (TP_Prodigal + FP_Prodigal),
         Prodigal_FNR = FN_Prodigal / (TP_Prodigal + FN_Prodigal),
         MetaProdigal_F1 = (2 * TP_MetaProdigal) / (2 * TP_MetaProdigal + FP_MetaProdigal + FN_MetaProdigal),
         MetaProdigal_Precision = TP_MetaProdigal / (TP_MetaProdigal + FP_MetaProdigal),
         MetaProdigal_Sensitivity = TP_MetaProdigal / (TP_MetaProdigal + FN_MetaProdigal),
         MetaProdigal_FDR = FP_MetaProdigal / (TP_MetaProdigal + FP_MetaProdigal),
         MetaProdigal_FNR = FN_MetaProdigal / (TP_MetaProdigal + FN_MetaProdigal),
         Glimmer_F1 = (2 * TP_Glimmer) / (2 * TP_Glimmer + FP_Glimmer + FN_Glimmer),
         Glimmer_Precision = TP_Glimmer / (TP_Glimmer + FP_Glimmer),
         Glimmer_Sensitivity = TP_Glimmer / (TP_Glimmer + FN_Glimmer),
         Glimmer_FDR = FP_Glimmer / (TP_Glimmer + FP_Glimmer),
         Glimmer_FNR = FN_Glimmer / (TP_Glimmer + FN_Glimmer),
         GeneMarkS_F1 = (2 * TP_GeneMarkS) / (2 * TP_GeneMarkS + FP_GeneMarkS + FN_GeneMarkS),
         GeneMarkS_Precision = TP_GeneMarkS / (TP_GeneMarkS + FP_GeneMarkS),
         GeneMarkS_Sensitivity = TP_GeneMarkS / (TP_GeneMarkS + FN_GeneMarkS),
         GeneMarkS_FDR = FP_GeneMarkS / (TP_GeneMarkS + FP_GeneMarkS),
         GeneMarkS_FNR = FN_GeneMarkS / (TP_GeneMarkS + FN_GeneMarkS),
         Phanotate_F1 = (2 * TP_Phanotate) / (2 * TP_Phanotate + FP_Phanotate + FN_Phanotate),
         Phanotate_Precision = TP_Phanotate / (TP_Phanotate + FP_Phanotate),
         Phanotate_Sensitivity = TP_Phanotate / (TP_Phanotate + FN_Phanotate),
         Phanotate_FDR = FP_Phanotate / (TP_Phanotate + FP_Phanotate),
         Phanotate_FNR = FN_Phanotate / (TP_Phanotate + FN_Phanotate),
         FragGeneScan_F1 = (2 * TP_FragGeneScan) / (2 * TP_FragGeneScan + FP_FragGeneScan + FN_FragGeneScan),
         FragGeneScan_Precision = TP_FragGeneScan / (TP_FragGeneScan + FP_FragGeneScan),
         FragGeneScan_Sensitivity = TP_FragGeneScan / (TP_FragGeneScan + FN_FragGeneScan),
         FragGeneScan_FDR = FP_FragGeneScan / (TP_FragGeneScan + FP_FragGeneScan),
         FragGeneScan_FNR = FN_FragGeneScan / (TP_FragGeneScan + FN_FragGeneScan),
         MGA_F1 = (2 * TP_MGA) / (2 * TP_MGA + FP_MGA + FN_MGA),
         MGA_Precision = TP_MGA / (TP_MGA + FP_MGA),
         MGA_Sensitivity = TP_MGA / (TP_MGA + FN_MGA),
         MGA_FDR = FP_MGA / (TP_MGA + FP_MGA),
         MGA_FNR = FN_MGA / (TP_MGA + FN_MGA),
         AUGUSTUS_SA_F1 = (2 * TP_AUGUSTUS_SA) / (2 * TP_AUGUSTUS_SA + FP_AUGUSTUS_SA + FN_AUGUSTUS_SA),
         AUGUSTUS_SA_Precision = TP_AUGUSTUS_SA / (TP_AUGUSTUS_SA + FP_AUGUSTUS_SA),
         AUGUSTUS_SA_Sensitivity = TP_AUGUSTUS_SA / (TP_AUGUSTUS_SA + FN_AUGUSTUS_SA),
         AUGUSTUS_SA_FDR = FP_AUGUSTUS_SA / (TP_AUGUSTUS_SA + FP_AUGUSTUS_SA),
         AUGUSTUS_SA_FNR = FN_AUGUSTUS_SA / (TP_AUGUSTUS_SA + FN_AUGUSTUS_SA),
         AUGUSTUS_EC_F1 = (2 * TP_AUGUSTUS_EC) / (2 * TP_AUGUSTUS_EC + FP_AUGUSTUS_EC + FN_AUGUSTUS_EC),
         AUGUSTUS_EC_Precision = TP_AUGUSTUS_EC / (TP_AUGUSTUS_EC + FP_AUGUSTUS_EC),
         AUGUSTUS_EC_Sensitivity = TP_AUGUSTUS_EC / (TP_AUGUSTUS_EC + FN_AUGUSTUS_EC),
         AUGUSTUS_EC_FDR = FP_AUGUSTUS_EC / (TP_AUGUSTUS_EC + FP_AUGUSTUS_EC),
         AUGUSTUS_EC_FNR = FN_AUGUSTUS_EC / (TP_AUGUSTUS_EC + FN_AUGUSTUS_EC),
         AUGUSTUS_HS_F1 = (2 * TP_AUGUSTUS_HS) / (2 * TP_AUGUSTUS_HS + FP_AUGUSTUS_HS + FN_AUGUSTUS_HS),
         AUGUSTUS_HS_Precision = TP_AUGUSTUS_HS / (TP_AUGUSTUS_HS + FP_AUGUSTUS_HS),
         AUGUSTUS_HS_Sensitivity = TP_AUGUSTUS_HS / (TP_AUGUSTUS_HS + FN_AUGUSTUS_HS),
         AUGUSTUS_HS_FDR = FP_AUGUSTUS_HS / (TP_AUGUSTUS_HS + FP_AUGUSTUS_HS),
         AUGUSTUS_HS_FNR = FN_AUGUSTUS_HS / (TP_AUGUSTUS_HS + FN_AUGUSTUS_HS))

# Eukaryotes (Subsapling)

listofdfs = list()
for(i in 1:1000){
  listofdfs[[i]] <- Basal_table %>%
    filter(Host_domain == "Eukaryote") %>%
    slice_sample(prop = 0.1)
}
dfofdfs <- tibble(d = 1:1000, data = listofdfs)

sum_euk_subsamples <- list()
for(i in 1:1000){
  sum_euk_subsamples[[i]] <- dfofdfs$data[[i]] %>%
    summarize_if(is.numeric, sum, na.rm=TRUE)
}

subsamples_sum_euk <- bind_rows(sum_euk_subsamples, .id = "subsample_number") %>%
  rowwise() %>%
  mutate(Prodigal_F1 = (2 * TP_Prodigal) / (2 * TP_Prodigal + FP_Prodigal + FN_Prodigal),
         Prodigal_Precision = TP_Prodigal / (TP_Prodigal + FP_Prodigal),
         Prodigal_Sensitivity = TP_Prodigal / (TP_Prodigal + FN_Prodigal),
         Prodigal_FDR = FP_Prodigal / (TP_Prodigal + FP_Prodigal),
         Prodigal_FNR = FN_Prodigal / (TP_Prodigal + FN_Prodigal),
         MetaProdigal_F1 = (2 * TP_MetaProdigal) / (2 * TP_MetaProdigal + FP_MetaProdigal + FN_MetaProdigal),
         MetaProdigal_Precision = TP_MetaProdigal / (TP_MetaProdigal + FP_MetaProdigal),
         MetaProdigal_Sensitivity = TP_MetaProdigal / (TP_MetaProdigal + FN_MetaProdigal),
         MetaProdigal_FDR = FP_MetaProdigal / (TP_MetaProdigal + FP_MetaProdigal),
         MetaProdigal_FNR = FN_MetaProdigal / (TP_MetaProdigal + FN_MetaProdigal),
         Glimmer_F1 = (2 * TP_Glimmer) / (2 * TP_Glimmer + FP_Glimmer + FN_Glimmer),
         Glimmer_Precision = TP_Glimmer / (TP_Glimmer + FP_Glimmer),
         Glimmer_Sensitivity = TP_Glimmer / (TP_Glimmer + FN_Glimmer),
         Glimmer_FDR = FP_Glimmer / (TP_Glimmer + FP_Glimmer),
         Glimmer_FNR = FN_Glimmer / (TP_Glimmer + FN_Glimmer),
         GeneMarkS_F1 = (2 * TP_GeneMarkS) / (2 * TP_GeneMarkS + FP_GeneMarkS + FN_GeneMarkS),
         GeneMarkS_Precision = TP_GeneMarkS / (TP_GeneMarkS + FP_GeneMarkS),
         GeneMarkS_Sensitivity = TP_GeneMarkS / (TP_GeneMarkS + FN_GeneMarkS),
         GeneMarkS_FDR = FP_GeneMarkS / (TP_GeneMarkS + FP_GeneMarkS),
         GeneMarkS_FNR = FN_GeneMarkS / (TP_GeneMarkS + FN_GeneMarkS),
         Phanotate_F1 = (2 * TP_Phanotate) / (2 * TP_Phanotate + FP_Phanotate + FN_Phanotate),
         Phanotate_Precision = TP_Phanotate / (TP_Phanotate + FP_Phanotate),
         Phanotate_Sensitivity = TP_Phanotate / (TP_Phanotate + FN_Phanotate),
         Phanotate_FDR = FP_Phanotate / (TP_Phanotate + FP_Phanotate),
         Phanotate_FNR = FN_Phanotate / (TP_Phanotate + FN_Phanotate),
         FragGeneScan_F1 = (2 * TP_FragGeneScan) / (2 * TP_FragGeneScan + FP_FragGeneScan + FN_FragGeneScan),
         FragGeneScan_Precision = TP_FragGeneScan / (TP_FragGeneScan + FP_FragGeneScan),
         FragGeneScan_Sensitivity = TP_FragGeneScan / (TP_FragGeneScan + FN_FragGeneScan),
         FragGeneScan_FDR = FP_FragGeneScan / (TP_FragGeneScan + FP_FragGeneScan),
         FragGeneScan_FNR = FN_FragGeneScan / (TP_FragGeneScan + FN_FragGeneScan),
         MGA_F1 = (2 * TP_MGA) / (2 * TP_MGA + FP_MGA + FN_MGA),
         MGA_Precision = TP_MGA / (TP_MGA + FP_MGA),
         MGA_Sensitivity = TP_MGA / (TP_MGA + FN_MGA),
         MGA_FDR = FP_MGA / (TP_MGA + FP_MGA),
         MGA_FNR = FN_MGA / (TP_MGA + FN_MGA),
         AUGUSTUSSA_F1 = (2 * TP_AUGUSTUS_SA) / (2 * TP_AUGUSTUS_SA + FP_AUGUSTUS_SA + FN_AUGUSTUS_SA),
         AUGUSTUSSA_Precision = TP_AUGUSTUS_SA / (TP_AUGUSTUS_SA + FP_AUGUSTUS_SA),
         AUGUSTUSSA_Sensitivity = TP_AUGUSTUS_SA / (TP_AUGUSTUS_SA + FN_AUGUSTUS_SA),
         AUGUSTUSSA_FDR = FP_AUGUSTUS_SA / (TP_AUGUSTUS_SA + FP_AUGUSTUS_SA),
         AUGUSTUSSA_FNR = FN_AUGUSTUS_SA / (TP_AUGUSTUS_SA + FN_AUGUSTUS_SA),
         AUGUSTUSEC_F1 = (2 * TP_AUGUSTUS_EC) / (2 * TP_AUGUSTUS_EC + FP_AUGUSTUS_EC + FN_AUGUSTUS_EC),
         AUGUSTUSEC_Precision = TP_AUGUSTUS_EC / (TP_AUGUSTUS_EC + FP_AUGUSTUS_EC),
         AUGUSTUSEC_Sensitivity = TP_AUGUSTUS_EC / (TP_AUGUSTUS_EC + FN_AUGUSTUS_EC),
         AUGUSTUSEC_FDR = FP_AUGUSTUS_EC / (TP_AUGUSTUS_EC + FP_AUGUSTUS_EC),
         AUGUSTUSEC_FNR = FN_AUGUSTUS_EC / (TP_AUGUSTUS_EC + FN_AUGUSTUS_EC),
         AUGUSTUSHS_F1 = (2 * TP_AUGUSTUS_HS) / (2 * TP_AUGUSTUS_HS + FP_AUGUSTUS_HS + FN_AUGUSTUS_HS),
         AUGUSTUSHS_Precision = TP_AUGUSTUS_HS / (TP_AUGUSTUS_HS + FP_AUGUSTUS_HS),
         AUGUSTUSHS_Sensitivity = TP_AUGUSTUS_HS / (TP_AUGUSTUS_HS + FN_AUGUSTUS_HS),
         AUGUSTUSHS_FDR = FP_AUGUSTUS_HS / (TP_AUGUSTUS_HS + FP_AUGUSTUS_HS),
         AUGUSTUSHS_FNR = FN_AUGUSTUS_HS / (TP_AUGUSTUS_HS + FN_AUGUSTUS_HS))

#Testing normality for the variables F1, Precision, Sensitivity, FDR and FNR of the different programs
listpvalues_shapiro_euk = list()

for(i in 32:length(subsamples_sum_euk)){
  listpvalues_shapiro_euk[[i]] <- shapiro.test(subsamples_sum_euk[[i]])$p.value
}

listpvalues_shapiro_euk <- Filter(Negate(is.null), listpvalues_shapiro_euk)
listpvalues_shapiro_euk <- as.list(p.adjust(listpvalues_shapiro_euk, method = "fdr"))
datapvalues_shapiro_euk <- tibble(d = 1:length(listpvalues_shapiro_euk),
                              variable = c("Prodigal_F1", "Prodigal_Precision", "Prodigal_Sensitivity", "Prodigal_FDR", "Prodigal_FNR",
                                           "MetaProdigal_F1", "MetaProdigal_Precision", "MetaProdigal_Sensitivity", "MetaProdigal_FDR", "MetaProdigal_FNR",
                                           "Glimmer_F1", "Glimmer_Precision", "Glimmer_Sensitivity", "Glimmer_FDR", "Glimmer_FNR",
                                           "GeneMarkS_F1", "GeneMarkS_Precision", "GeneMarkS_Sensitivity", "GeneMarkS_FDR", "GeneMarkS_FNR",
                                           "Phanotate_F1", "Phanotate_Precision", "Phanotate_Sensitivity", "Phanotate_FDR", "Phanotate_FNR", 
                                           "FragGeneScan_F1", "FragGeneScan_Precision", "FragGeneScan_Sensitivity", "FragGeneScan_FDR", "FragGeneScan_FNR", 
                                           "MGA_F1", "MGA_Precision", "MGA_Sensitivity", "MGA_FDR", "MGA_FNR",
                                           "AUGUSTUS_SA_F1", "AUGUSTUS_SA_Precision", "AUGUSTUS_SA_Sensitivity", "AUGUSTUS_SA_FDR", "AUGUSTUS_SA_FNR",
                                           "AUGUSTUS_EC_F1", "AUGUSTUS_EC_Precision", "AUGUSTUS_EC_Sensitivity", "AUGUSTUS_EC_FDR", "AUGUSTUS_EC_FNR",
                                           "AUGUSTUS_HS_F1", "AUGUSTUS_HS_Precision", "AUGUSTUS_HS_Sensitivity", "AUGUSTUS_HS_FDR", "AUGUSTUS_HS_FNR"),
                              p_value = listpvalues_shapiro_euk) %>%
  select(-d) %>%
  pivot_wider(names_from = variable, values_from = p_value)

#All p-values (or near all) were significant, so there is no chance to assume normality in our data
#Thus, the only way to compare values among them is via Kruskal-Wallis followed by pairwise tests (Dunn test)

#Measuring the differences between F1 values
datastatsrepeuk <- subsamples_sum_euk %>%
  select(c(1, 32:length(subsamples_sum_euk))) %>%
  gather(key = var_name, value = value, 2:51) %>%
  separate(col = var_name, sep = "_", into = c("Program", "Statistic")) %>%
  select(-subsample_number)

F1_datarep_euk <- datastatsrepeuk %>% filter(Statistic == "F1")
Precision_datarep_euk <- datastatsrepeuk %>% filter(Statistic == "Precision")
Sensitivity_datarep_euk <- datastatsrepeuk %>% filter(Statistic == "Sensitivity")
FDR_datarep_euk <- datastatsrepeuk %>% filter(Statistic == "FDR")
FNR_datarep_euk <- datastatsrepeuk %>% filter(Statistic == "FNR")

kruskal.test(x = F1_datarep_euk$value, g = F1_datarep_euk$Program)
dunnTest(value ~ Program, data = F1_datarep_euk, method = "holm")

kruskal.test(x = Precision_datarep_euk$value, g = Precision_datarep_euk$Program)
dunnTest(value ~ Program, data = Precision_datarep_euk, method = "holm")

kruskal.test(x = Sensitivity_datarep_euk$value, g = Sensitivity_datarep_euk$Program)
dunnTest(value ~ Program, data = Sensitivity_datarep_euk, method = "holm")

kruskal.test(x = FDR_datarep_euk$value, g = FDR_datarep_euk$Program)
dunnTest(value ~ Program, data = FDR_datarep_euk, method = "holm")

kruskal.test(x = FNR_datarep_euk$value, g = FNR_datarep_euk$Program)
dunnTest(value ~ Program, data = FNR_datarep_euk, method = "holm")

# Bacteria (Real)

sum_bac <- Basal_table %>%
  filter(Host_domain == "Bacteria") %>%
  summarize_if(is.numeric, sum, na.rm=TRUE) %>%
  rowwise() %>%
  mutate(Prodigal_F1 = (2 * TP_Prodigal) / (2 * TP_Prodigal + FP_Prodigal + FN_Prodigal),
         Prodigal_Precision = TP_Prodigal / (TP_Prodigal + FP_Prodigal),
         Prodigal_Sensitivity = TP_Prodigal / (TP_Prodigal + FN_Prodigal),
         Prodigal_FDR = FP_Prodigal / (TP_Prodigal + FP_Prodigal),
         Prodigal_FNR = FN_Prodigal / (TP_Prodigal + FN_Prodigal),
         MetaProdigal_F1 = (2 * TP_MetaProdigal) / (2 * TP_MetaProdigal + FP_MetaProdigal + FN_MetaProdigal),
         MetaProdigal_Precision = TP_MetaProdigal / (TP_MetaProdigal + FP_MetaProdigal),
         MetaProdigal_Sensitivity = TP_MetaProdigal / (TP_MetaProdigal + FN_MetaProdigal),
         MetaProdigal_FDR = FP_MetaProdigal / (TP_MetaProdigal + FP_MetaProdigal),
         MetaProdigal_FNR = FN_MetaProdigal / (TP_MetaProdigal + FN_MetaProdigal),
         Glimmer_F1 = (2 * TP_Glimmer) / (2 * TP_Glimmer + FP_Glimmer + FN_Glimmer),
         Glimmer_Precision = TP_Glimmer / (TP_Glimmer + FP_Glimmer),
         Glimmer_Sensitivity = TP_Glimmer / (TP_Glimmer + FN_Glimmer),
         Glimmer_FDR = FP_Glimmer / (TP_Glimmer + FP_Glimmer),
         Glimmer_FNR = FN_Glimmer / (TP_Glimmer + FN_Glimmer),
         GeneMarkS_F1 = (2 * TP_GeneMarkS) / (2 * TP_GeneMarkS + FP_GeneMarkS + FN_GeneMarkS),
         GeneMarkS_Precision = TP_GeneMarkS / (TP_GeneMarkS + FP_GeneMarkS),
         GeneMarkS_Sensitivity = TP_GeneMarkS / (TP_GeneMarkS + FN_GeneMarkS),
         GeneMarkS_FDR = FP_GeneMarkS / (TP_GeneMarkS + FP_GeneMarkS),
         GeneMarkS_FNR = FN_GeneMarkS / (TP_GeneMarkS + FN_GeneMarkS),
         Phanotate_F1 = (2 * TP_Phanotate) / (2 * TP_Phanotate + FP_Phanotate + FN_Phanotate),
         Phanotate_Precision = TP_Phanotate / (TP_Phanotate + FP_Phanotate),
         Phanotate_Sensitivity = TP_Phanotate / (TP_Phanotate + FN_Phanotate),
         Phanotate_FDR = FP_Phanotate / (TP_Phanotate + FP_Phanotate),
         Phanotate_FNR = FN_Phanotate / (TP_Phanotate + FN_Phanotate),
         FragGeneScan_F1 = (2 * TP_FragGeneScan) / (2 * TP_FragGeneScan + FP_FragGeneScan + FN_FragGeneScan),
         FragGeneScan_Precision = TP_FragGeneScan / (TP_FragGeneScan + FP_FragGeneScan),
         FragGeneScan_Sensitivity = TP_FragGeneScan / (TP_FragGeneScan + FN_FragGeneScan),
         FragGeneScan_FDR = FP_FragGeneScan / (TP_FragGeneScan + FP_FragGeneScan),
         FragGeneScan_FNR = FN_FragGeneScan / (TP_FragGeneScan + FN_FragGeneScan),
         MGA_F1 = (2 * TP_MGA) / (2 * TP_MGA + FP_MGA + FN_MGA),
         MGA_Precision = TP_MGA / (TP_MGA + FP_MGA),
         MGA_Sensitivity = TP_MGA / (TP_MGA + FN_MGA),
         MGA_FDR = FP_MGA / (TP_MGA + FP_MGA),
         MGA_FNR = FN_MGA / (TP_MGA + FN_MGA),
         AUGUSTUS_SA_F1 = (2 * TP_AUGUSTUS_SA) / (2 * TP_AUGUSTUS_SA + FP_AUGUSTUS_SA + FN_AUGUSTUS_SA),
         AUGUSTUS_SA_Precision = TP_AUGUSTUS_SA / (TP_AUGUSTUS_SA + FP_AUGUSTUS_SA),
         AUGUSTUS_SA_Sensitivity = TP_AUGUSTUS_SA / (TP_AUGUSTUS_SA + FN_AUGUSTUS_SA),
         AUGUSTUS_SA_FDR = FP_AUGUSTUS_SA / (TP_AUGUSTUS_SA + FP_AUGUSTUS_SA),
         AUGUSTUS_SA_FNR = FN_AUGUSTUS_SA / (TP_AUGUSTUS_SA + FN_AUGUSTUS_SA),
         AUGUSTUS_EC_F1 = (2 * TP_AUGUSTUS_EC) / (2 * TP_AUGUSTUS_EC + FP_AUGUSTUS_EC + FN_AUGUSTUS_EC),
         AUGUSTUS_EC_Precision = TP_AUGUSTUS_EC / (TP_AUGUSTUS_EC + FP_AUGUSTUS_EC),
         AUGUSTUS_EC_Sensitivity = TP_AUGUSTUS_EC / (TP_AUGUSTUS_EC + FN_AUGUSTUS_EC),
         AUGUSTUS_EC_FDR = FP_AUGUSTUS_EC / (TP_AUGUSTUS_EC + FP_AUGUSTUS_EC),
         AUGUSTUS_EC_FNR = FN_AUGUSTUS_EC / (TP_AUGUSTUS_EC + FN_AUGUSTUS_EC),
         AUGUSTUS_HS_F1 = (2 * TP_AUGUSTUS_HS) / (2 * TP_AUGUSTUS_HS + FP_AUGUSTUS_HS + FN_AUGUSTUS_HS),
         AUGUSTUS_HS_Precision = TP_AUGUSTUS_HS / (TP_AUGUSTUS_HS + FP_AUGUSTUS_HS),
         AUGUSTUS_HS_Sensitivity = TP_AUGUSTUS_HS / (TP_AUGUSTUS_HS + FN_AUGUSTUS_HS),
         AUGUSTUS_HS_FDR = FP_AUGUSTUS_HS / (TP_AUGUSTUS_HS + FP_AUGUSTUS_HS),
         AUGUSTUS_HS_FNR = FN_AUGUSTUS_HS / (TP_AUGUSTUS_HS + FN_AUGUSTUS_HS))

# Bacteria (Subsampling)

listofdfs = list()
for(i in 1:1000){
  listofdfs[[i]] <- Basal_table %>%
    filter(Host_domain == "Bacteria") %>%
    slice_sample(prop = 0.1)
}
dfofdfs <- tibble(d = 1:1000, data = listofdfs)

sum_bac_subsamples <- list()
for(i in 1:1000){
  sum_bac_subsamples[[i]] <- dfofdfs$data[[i]] %>%
    summarize_if(is.numeric, sum, na.rm=TRUE)
}

subsamples_sum_bac <- bind_rows(sum_bac_subsamples, .id = "subsample_number") %>%
  rowwise() %>%
  mutate(Prodigal_F1 = (2 * TP_Prodigal) / (2 * TP_Prodigal + FP_Prodigal + FN_Prodigal),
         Prodigal_Precision = TP_Prodigal / (TP_Prodigal + FP_Prodigal),
         Prodigal_Sensitivity = TP_Prodigal / (TP_Prodigal + FN_Prodigal),
         Prodigal_FDR = FP_Prodigal / (TP_Prodigal + FP_Prodigal),
         Prodigal_FNR = FN_Prodigal / (TP_Prodigal + FN_Prodigal),
         MetaProdigal_F1 = (2 * TP_MetaProdigal) / (2 * TP_MetaProdigal + FP_MetaProdigal + FN_MetaProdigal),
         MetaProdigal_Precision = TP_MetaProdigal / (TP_MetaProdigal + FP_MetaProdigal),
         MetaProdigal_Sensitivity = TP_MetaProdigal / (TP_MetaProdigal + FN_MetaProdigal),
         MetaProdigal_FDR = FP_MetaProdigal / (TP_MetaProdigal + FP_MetaProdigal),
         MetaProdigal_FNR = FN_MetaProdigal / (TP_MetaProdigal + FN_MetaProdigal),
         Glimmer_F1 = (2 * TP_Glimmer) / (2 * TP_Glimmer + FP_Glimmer + FN_Glimmer),
         Glimmer_Precision = TP_Glimmer / (TP_Glimmer + FP_Glimmer),
         Glimmer_Sensitivity = TP_Glimmer / (TP_Glimmer + FN_Glimmer),
         Glimmer_FDR = FP_Glimmer / (TP_Glimmer + FP_Glimmer),
         Glimmer_FNR = FN_Glimmer / (TP_Glimmer + FN_Glimmer),
         GeneMarkS_F1 = (2 * TP_GeneMarkS) / (2 * TP_GeneMarkS + FP_GeneMarkS + FN_GeneMarkS),
         GeneMarkS_Precision = TP_GeneMarkS / (TP_GeneMarkS + FP_GeneMarkS),
         GeneMarkS_Sensitivity = TP_GeneMarkS / (TP_GeneMarkS + FN_GeneMarkS),
         GeneMarkS_FDR = FP_GeneMarkS / (TP_GeneMarkS + FP_GeneMarkS),
         GeneMarkS_FNR = FN_GeneMarkS / (TP_GeneMarkS + FN_GeneMarkS),
         Phanotate_F1 = (2 * TP_Phanotate) / (2 * TP_Phanotate + FP_Phanotate + FN_Phanotate),
         Phanotate_Precision = TP_Phanotate / (TP_Phanotate + FP_Phanotate),
         Phanotate_Sensitivity = TP_Phanotate / (TP_Phanotate + FN_Phanotate),
         Phanotate_FDR = FP_Phanotate / (TP_Phanotate + FP_Phanotate),
         Phanotate_FNR = FN_Phanotate / (TP_Phanotate + FN_Phanotate),
         FragGeneScan_F1 = (2 * TP_FragGeneScan) / (2 * TP_FragGeneScan + FP_FragGeneScan + FN_FragGeneScan),
         FragGeneScan_Precision = TP_FragGeneScan / (TP_FragGeneScan + FP_FragGeneScan),
         FragGeneScan_Sensitivity = TP_FragGeneScan / (TP_FragGeneScan + FN_FragGeneScan),
         FragGeneScan_FDR = FP_FragGeneScan / (TP_FragGeneScan + FP_FragGeneScan),
         FragGeneScan_FNR = FN_FragGeneScan / (TP_FragGeneScan + FN_FragGeneScan),
         MGA_F1 = (2 * TP_MGA) / (2 * TP_MGA + FP_MGA + FN_MGA),
         MGA_Precision = TP_MGA / (TP_MGA + FP_MGA),
         MGA_Sensitivity = TP_MGA / (TP_MGA + FN_MGA),
         MGA_FDR = FP_MGA / (TP_MGA + FP_MGA),
         MGA_FNR = FN_MGA / (TP_MGA + FN_MGA),
         AUGUSTUSSA_F1 = (2 * TP_AUGUSTUS_SA) / (2 * TP_AUGUSTUS_SA + FP_AUGUSTUS_SA + FN_AUGUSTUS_SA),
         AUGUSTUSSA_Precision = TP_AUGUSTUS_SA / (TP_AUGUSTUS_SA + FP_AUGUSTUS_SA),
         AUGUSTUSSA_Sensitivity = TP_AUGUSTUS_SA / (TP_AUGUSTUS_SA + FN_AUGUSTUS_SA),
         AUGUSTUSSA_FDR = FP_AUGUSTUS_SA / (TP_AUGUSTUS_SA + FP_AUGUSTUS_SA),
         AUGUSTUSSA_FNR = FN_AUGUSTUS_SA / (TP_AUGUSTUS_SA + FN_AUGUSTUS_SA),
         AUGUSTUSEC_F1 = (2 * TP_AUGUSTUS_EC) / (2 * TP_AUGUSTUS_EC + FP_AUGUSTUS_EC + FN_AUGUSTUS_EC),
         AUGUSTUSEC_Precision = TP_AUGUSTUS_EC / (TP_AUGUSTUS_EC + FP_AUGUSTUS_EC),
         AUGUSTUSEC_Sensitivity = TP_AUGUSTUS_EC / (TP_AUGUSTUS_EC + FN_AUGUSTUS_EC),
         AUGUSTUSEC_FDR = FP_AUGUSTUS_EC / (TP_AUGUSTUS_EC + FP_AUGUSTUS_EC),
         AUGUSTUSEC_FNR = FN_AUGUSTUS_EC / (TP_AUGUSTUS_EC + FN_AUGUSTUS_EC),
         AUGUSTUSHS_F1 = (2 * TP_AUGUSTUS_HS) / (2 * TP_AUGUSTUS_HS + FP_AUGUSTUS_HS + FN_AUGUSTUS_HS),
         AUGUSTUSHS_Precision = TP_AUGUSTUS_HS / (TP_AUGUSTUS_HS + FP_AUGUSTUS_HS),
         AUGUSTUSHS_Sensitivity = TP_AUGUSTUS_HS / (TP_AUGUSTUS_HS + FN_AUGUSTUS_HS),
         AUGUSTUSHS_FDR = FP_AUGUSTUS_HS / (TP_AUGUSTUS_HS + FP_AUGUSTUS_HS),
         AUGUSTUSHS_FNR = FN_AUGUSTUS_HS / (TP_AUGUSTUS_HS + FN_AUGUSTUS_HS))

#Testing normality for the variables F1, Precision, Sensitivity, FDR and FNR of the different programs
listpvalues_shapiro_bac = list()

for(i in 32:length(subsamples_sum_bac)){
  listpvalues_shapiro_bac[[i]] <- shapiro.test(subsamples_sum_bac[[i]])$p.value
}

listpvalues_shapiro_bac <- Filter(Negate(is.null), listpvalues_shapiro_bac)
listpvalues_shapiro_bac <- as.list(p.adjust(listpvalues_shapiro_bac, method = "fdr"))
datapvalues_shapiro_bac <- tibble(d = 1:length(listpvalues_shapiro_bac),
                                  variable = c("Prodigal_F1", "Prodigal_Precision", "Prodigal_Sensitivity", "Prodigal_FDR", "Prodigal_FNR",
                                               "MetaProdigal_F1", "MetaProdigal_Precision", "MetaProdigal_Sensitivity", "MetaProdigal_FDR", "MetaProdigal_FNR",
                                               "Glimmer_F1", "Glimmer_Precision", "Glimmer_Sensitivity", "Glimmer_FDR", "Glimmer_FNR",
                                               "GeneMarkS_F1", "GeneMarkS_Precision", "GeneMarkS_Sensitivity", "GeneMarkS_FDR", "GeneMarkS_FNR",
                                               "Phanotate_F1", "Phanotate_Precision", "Phanotate_Sensitivity", "Phanotate_FDR", "Phanotate_FNR", 
                                               "FragGeneScan_F1", "FragGeneScan_Precision", "FragGeneScan_Sensitivity", "FragGeneScan_FDR", "FragGeneScan_FNR", 
                                               "MGA_F1", "MGA_Precision", "MGA_Sensitivity", "MGA_FDR", "MGA_FNR",
                                               "AUGUSTUS_SA_F1", "AUGUSTUS_SA_Precision", "AUGUSTUS_SA_Sensitivity", "AUGUSTUS_SA_FDR", "AUGUSTUS_SA_FNR",
                                               "AUGUSTUS_EC_F1", "AUGUSTUS_EC_Precision", "AUGUSTUS_EC_Sensitivity", "AUGUSTUS_EC_FDR", "AUGUSTUS_EC_FNR",
                                               "AUGUSTUS_HS_F1", "AUGUSTUS_HS_Precision", "AUGUSTUS_HS_Sensitivity", "AUGUSTUS_HS_FDR", "AUGUSTUS_HS_FNR"),
                                  p_value = listpvalues_shapiro_bac) %>%
  select(-d) %>%
  pivot_wider(names_from = variable, values_from = p_value)

#Some p-values were significant, so there is no chance to assume normality in our data
#Thus, the only way to compare values among them is via Kruskal-Wallis followed by pairwise tests (Dunn test)

#Measuring the differences between F1 values
datastatsrepbac <- subsamples_sum_bac %>%
  select(c(1, 32:length(subsamples_sum_bac))) %>%
  gather(key = var_name, value = value, 2:51) %>%
  separate(col = var_name, sep = "_", into = c("Program", "Statistic")) %>%
  select(-subsample_number)

F1_datarep_bac <- datastatsrepbac %>% filter(Statistic == "F1")
Precision_datarep_bac <- datastatsrepbac %>% filter(Statistic == "Precision")
Sensitivity_datarep_bac <- datastatsrepbac %>% filter(Statistic == "Sensitivity")
FDR_datarep_bac <- datastatsrepbac %>% filter(Statistic == "FDR")
FNR_datarep_bac <- datastatsrepbac %>% filter(Statistic == "FNR")

kruskal.test(x = F1_datarep_bac$value, g = F1_datarep_bac$Program)
dunnTest(value ~ Program, data = F1_datarep_bac, method = "holm")

kruskal.test(x = Precision_datarep_bac$value, g = Precision_datarep_bac$Program)
dunnTest(value ~ Program, data = Precision_datarep_bac, method = "holm")

kruskal.test(x = Sensitivity_datarep_bac$value, g = Sensitivity_datarep_bac$Program)
dunnTest(value ~ Program, data = Sensitivity_datarep_bac, method = "holm")

kruskal.test(x = FDR_datarep_bac$value, g = FDR_datarep_bac$Program)
dunnTest(value ~ Program, data = FDR_datarep_bac, method = "holm")

kruskal.test(x = FNR_datarep_bac$value, g = FNR_datarep_bac$Program)
dunnTest(value ~ Program, data = FNR_datarep_bac, method = "holm")

# Archaea (Real)

sum_arc <- Basal_table %>%
  filter(Host_domain == "Archaea") %>%
  summarize_if(is.numeric, sum, na.rm=TRUE) %>%
  rowwise() %>%
  mutate(Prodigal_F1 = (2 * TP_Prodigal) / (2 * TP_Prodigal + FP_Prodigal + FN_Prodigal),
         Prodigal_Precision = TP_Prodigal / (TP_Prodigal + FP_Prodigal),
         Prodigal_Sensitivity = TP_Prodigal / (TP_Prodigal + FN_Prodigal),
         Prodigal_FDR = FP_Prodigal / (TP_Prodigal + FP_Prodigal),
         Prodigal_FNR = FN_Prodigal / (TP_Prodigal + FN_Prodigal),
         MetaProdigal_F1 = (2 * TP_MetaProdigal) / (2 * TP_MetaProdigal + FP_MetaProdigal + FN_MetaProdigal),
         MetaProdigal_Precision = TP_MetaProdigal / (TP_MetaProdigal + FP_MetaProdigal),
         MetaProdigal_Sensitivity = TP_MetaProdigal / (TP_MetaProdigal + FN_MetaProdigal),
         MetaProdigal_FDR = FP_MetaProdigal / (TP_MetaProdigal + FP_MetaProdigal),
         MetaProdigal_FNR = FN_MetaProdigal / (TP_MetaProdigal + FN_MetaProdigal),
         Glimmer_F1 = (2 * TP_Glimmer) / (2 * TP_Glimmer + FP_Glimmer + FN_Glimmer),
         Glimmer_Precision = TP_Glimmer / (TP_Glimmer + FP_Glimmer),
         Glimmer_Sensitivity = TP_Glimmer / (TP_Glimmer + FN_Glimmer),
         Glimmer_FDR = FP_Glimmer / (TP_Glimmer + FP_Glimmer),
         Glimmer_FNR = FN_Glimmer / (TP_Glimmer + FN_Glimmer),
         GeneMarkS_F1 = (2 * TP_GeneMarkS) / (2 * TP_GeneMarkS + FP_GeneMarkS + FN_GeneMarkS),
         GeneMarkS_Precision = TP_GeneMarkS / (TP_GeneMarkS + FP_GeneMarkS),
         GeneMarkS_Sensitivity = TP_GeneMarkS / (TP_GeneMarkS + FN_GeneMarkS),
         GeneMarkS_FDR = FP_GeneMarkS / (TP_GeneMarkS + FP_GeneMarkS),
         GeneMarkS_FNR = FN_GeneMarkS / (TP_GeneMarkS + FN_GeneMarkS),
         Phanotate_F1 = (2 * TP_Phanotate) / (2 * TP_Phanotate + FP_Phanotate + FN_Phanotate),
         Phanotate_Precision = TP_Phanotate / (TP_Phanotate + FP_Phanotate),
         Phanotate_Sensitivity = TP_Phanotate / (TP_Phanotate + FN_Phanotate),
         Phanotate_FDR = FP_Phanotate / (TP_Phanotate + FP_Phanotate),
         Phanotate_FNR = FN_Phanotate / (TP_Phanotate + FN_Phanotate),
         FragGeneScan_F1 = (2 * TP_FragGeneScan) / (2 * TP_FragGeneScan + FP_FragGeneScan + FN_FragGeneScan),
         FragGeneScan_Precision = TP_FragGeneScan / (TP_FragGeneScan + FP_FragGeneScan),
         FragGeneScan_Sensitivity = TP_FragGeneScan / (TP_FragGeneScan + FN_FragGeneScan),
         FragGeneScan_FDR = FP_FragGeneScan / (TP_FragGeneScan + FP_FragGeneScan),
         FragGeneScan_FNR = FN_FragGeneScan / (TP_FragGeneScan + FN_FragGeneScan),
         MGA_F1 = (2 * TP_MGA) / (2 * TP_MGA + FP_MGA + FN_MGA),
         MGA_Precision = TP_MGA / (TP_MGA + FP_MGA),
         MGA_Sensitivity = TP_MGA / (TP_MGA + FN_MGA),
         MGA_FDR = FP_MGA / (TP_MGA + FP_MGA),
         MGA_FNR = FN_MGA / (TP_MGA + FN_MGA),
         AUGUSTUS_SA_F1 = (2 * TP_AUGUSTUS_SA) / (2 * TP_AUGUSTUS_SA + FP_AUGUSTUS_SA + FN_AUGUSTUS_SA),
         AUGUSTUS_SA_Precision = TP_AUGUSTUS_SA / (TP_AUGUSTUS_SA + FP_AUGUSTUS_SA),
         AUGUSTUS_SA_Sensitivity = TP_AUGUSTUS_SA / (TP_AUGUSTUS_SA + FN_AUGUSTUS_SA),
         AUGUSTUS_SA_FDR = FP_AUGUSTUS_SA / (TP_AUGUSTUS_SA + FP_AUGUSTUS_SA),
         AUGUSTUS_SA_FNR = FN_AUGUSTUS_SA / (TP_AUGUSTUS_SA + FN_AUGUSTUS_SA),
         AUGUSTUS_EC_F1 = (2 * TP_AUGUSTUS_EC) / (2 * TP_AUGUSTUS_EC + FP_AUGUSTUS_EC + FN_AUGUSTUS_EC),
         AUGUSTUS_EC_Precision = TP_AUGUSTUS_EC / (TP_AUGUSTUS_EC + FP_AUGUSTUS_EC),
         AUGUSTUS_EC_Sensitivity = TP_AUGUSTUS_EC / (TP_AUGUSTUS_EC + FN_AUGUSTUS_EC),
         AUGUSTUS_EC_FDR = FP_AUGUSTUS_EC / (TP_AUGUSTUS_EC + FP_AUGUSTUS_EC),
         AUGUSTUS_EC_FNR = FN_AUGUSTUS_EC / (TP_AUGUSTUS_EC + FN_AUGUSTUS_EC),
         AUGUSTUS_HS_F1 = (2 * TP_AUGUSTUS_HS) / (2 * TP_AUGUSTUS_HS + FP_AUGUSTUS_HS + FN_AUGUSTUS_HS),
         AUGUSTUS_HS_Precision = TP_AUGUSTUS_HS / (TP_AUGUSTUS_HS + FP_AUGUSTUS_HS),
         AUGUSTUS_HS_Sensitivity = TP_AUGUSTUS_HS / (TP_AUGUSTUS_HS + FN_AUGUSTUS_HS),
         AUGUSTUS_HS_FDR = FP_AUGUSTUS_HS / (TP_AUGUSTUS_HS + FP_AUGUSTUS_HS),
         AUGUSTUS_HS_FNR = FN_AUGUSTUS_HS / (TP_AUGUSTUS_HS + FN_AUGUSTUS_HS))

# Archaea (Subsampling)

listofdfs = list()
for(i in 1:1000){
  listofdfs[[i]] <- Basal_table %>%
    filter(Host_domain == "Archaea") %>%
    slice_sample(prop = 0.1)
}
dfofdfs <- tibble(d = 1:1000, data = listofdfs)

sum_arc_subsamples <- list()
for(i in 1:1000){
  sum_arc_subsamples[[i]] <- dfofdfs$data[[i]] %>%
    summarize_if(is.numeric, sum, na.rm=TRUE)
}

subsamples_sum_arc <- bind_rows(sum_arc_subsamples, .id = "subsample_number") %>%
  rowwise() %>%
  mutate(Prodigal_F1 = (2 * TP_Prodigal) / (2 * TP_Prodigal + FP_Prodigal + FN_Prodigal),
         Prodigal_Precision = TP_Prodigal / (TP_Prodigal + FP_Prodigal),
         Prodigal_Sensitivity = TP_Prodigal / (TP_Prodigal + FN_Prodigal),
         Prodigal_FDR = FP_Prodigal / (TP_Prodigal + FP_Prodigal),
         Prodigal_FNR = FN_Prodigal / (TP_Prodigal + FN_Prodigal),
         MetaProdigal_F1 = (2 * TP_MetaProdigal) / (2 * TP_MetaProdigal + FP_MetaProdigal + FN_MetaProdigal),
         MetaProdigal_Precision = TP_MetaProdigal / (TP_MetaProdigal + FP_MetaProdigal),
         MetaProdigal_Sensitivity = TP_MetaProdigal / (TP_MetaProdigal + FN_MetaProdigal),
         MetaProdigal_FDR = FP_MetaProdigal / (TP_MetaProdigal + FP_MetaProdigal),
         MetaProdigal_FNR = FN_MetaProdigal / (TP_MetaProdigal + FN_MetaProdigal),
         Glimmer_F1 = (2 * TP_Glimmer) / (2 * TP_Glimmer + FP_Glimmer + FN_Glimmer),
         Glimmer_Precision = TP_Glimmer / (TP_Glimmer + FP_Glimmer),
         Glimmer_Sensitivity = TP_Glimmer / (TP_Glimmer + FN_Glimmer),
         Glimmer_FDR = FP_Glimmer / (TP_Glimmer + FP_Glimmer),
         Glimmer_FNR = FN_Glimmer / (TP_Glimmer + FN_Glimmer),
         GeneMarkS_F1 = (2 * TP_GeneMarkS) / (2 * TP_GeneMarkS + FP_GeneMarkS + FN_GeneMarkS),
         GeneMarkS_Precision = TP_GeneMarkS / (TP_GeneMarkS + FP_GeneMarkS),
         GeneMarkS_Sensitivity = TP_GeneMarkS / (TP_GeneMarkS + FN_GeneMarkS),
         GeneMarkS_FDR = FP_GeneMarkS / (TP_GeneMarkS + FP_GeneMarkS),
         GeneMarkS_FNR = FN_GeneMarkS / (TP_GeneMarkS + FN_GeneMarkS),
         Phanotate_F1 = (2 * TP_Phanotate) / (2 * TP_Phanotate + FP_Phanotate + FN_Phanotate),
         Phanotate_Precision = TP_Phanotate / (TP_Phanotate + FP_Phanotate),
         Phanotate_Sensitivity = TP_Phanotate / (TP_Phanotate + FN_Phanotate),
         Phanotate_FDR = FP_Phanotate / (TP_Phanotate + FP_Phanotate),
         Phanotate_FNR = FN_Phanotate / (TP_Phanotate + FN_Phanotate),
         FragGeneScan_F1 = (2 * TP_FragGeneScan) / (2 * TP_FragGeneScan + FP_FragGeneScan + FN_FragGeneScan),
         FragGeneScan_Precision = TP_FragGeneScan / (TP_FragGeneScan + FP_FragGeneScan),
         FragGeneScan_Sensitivity = TP_FragGeneScan / (TP_FragGeneScan + FN_FragGeneScan),
         FragGeneScan_FDR = FP_FragGeneScan / (TP_FragGeneScan + FP_FragGeneScan),
         FragGeneScan_FNR = FN_FragGeneScan / (TP_FragGeneScan + FN_FragGeneScan),
         MGA_F1 = (2 * TP_MGA) / (2 * TP_MGA + FP_MGA + FN_MGA),
         MGA_Precision = TP_MGA / (TP_MGA + FP_MGA),
         MGA_Sensitivity = TP_MGA / (TP_MGA + FN_MGA),
         MGA_FDR = FP_MGA / (TP_MGA + FP_MGA),
         MGA_FNR = FN_MGA / (TP_MGA + FN_MGA),
         AUGUSTUSSA_F1 = (2 * TP_AUGUSTUS_SA) / (2 * TP_AUGUSTUS_SA + FP_AUGUSTUS_SA + FN_AUGUSTUS_SA),
         AUGUSTUSSA_Precision = TP_AUGUSTUS_SA / (TP_AUGUSTUS_SA + FP_AUGUSTUS_SA),
         AUGUSTUSSA_Sensitivity = TP_AUGUSTUS_SA / (TP_AUGUSTUS_SA + FN_AUGUSTUS_SA),
         AUGUSTUSSA_FDR = FP_AUGUSTUS_SA / (TP_AUGUSTUS_SA + FP_AUGUSTUS_SA),
         AUGUSTUSSA_FNR = FN_AUGUSTUS_SA / (TP_AUGUSTUS_SA + FN_AUGUSTUS_SA),
         AUGUSTUSEC_F1 = (2 * TP_AUGUSTUS_EC) / (2 * TP_AUGUSTUS_EC + FP_AUGUSTUS_EC + FN_AUGUSTUS_EC),
         AUGUSTUSEC_Precision = TP_AUGUSTUS_EC / (TP_AUGUSTUS_EC + FP_AUGUSTUS_EC),
         AUGUSTUSEC_Sensitivity = TP_AUGUSTUS_EC / (TP_AUGUSTUS_EC + FN_AUGUSTUS_EC),
         AUGUSTUSEC_FDR = FP_AUGUSTUS_EC / (TP_AUGUSTUS_EC + FP_AUGUSTUS_EC),
         AUGUSTUSEC_FNR = FN_AUGUSTUS_EC / (TP_AUGUSTUS_EC + FN_AUGUSTUS_EC),
         AUGUSTUSHS_F1 = (2 * TP_AUGUSTUS_HS) / (2 * TP_AUGUSTUS_HS + FP_AUGUSTUS_HS + FN_AUGUSTUS_HS),
         AUGUSTUSHS_Precision = TP_AUGUSTUS_HS / (TP_AUGUSTUS_HS + FP_AUGUSTUS_HS),
         AUGUSTUSHS_Sensitivity = TP_AUGUSTUS_HS / (TP_AUGUSTUS_HS + FN_AUGUSTUS_HS),
         AUGUSTUSHS_FDR = FP_AUGUSTUS_HS / (TP_AUGUSTUS_HS + FP_AUGUSTUS_HS),
         AUGUSTUSHS_FNR = FN_AUGUSTUS_HS / (TP_AUGUSTUS_HS + FN_AUGUSTUS_HS))

#Testing normality for the variables F1, Precision, Sensitivity, FDR and FNR of the different programs
listpvalues_shapiro_arc = list()

for(i in 32:length(subsamples_sum_arc)){
  listpvalues_shapiro_arc[[i]] <- shapiro.test(subsamples_sum_arc[[i]])$p.value
}

listpvalues_shapiro_arc <- Filter(Negate(is.null), listpvalues_shapiro_arc)
listpvalues_shapiro_arc <- as.list(p.adjust(listpvalues_shapiro_arc, method = "fdr"))
datapvalues_shapiro_arc <- tibble(d = 1:length(listpvalues_shapiro_arc),
                                  variable = c("Prodigal_F1", "Prodigal_Precision", "Prodigal_Sensitivity", "Prodigal_FDR", "Prodigal_FNR",
                                               "MetaProdigal_F1", "MetaProdigal_Precision", "MetaProdigal_Sensitivity", "MetaProdigal_FDR", "MetaProdigal_FNR",
                                               "Glimmer_F1", "Glimmer_Precision", "Glimmer_Sensitivity", "Glimmer_FDR", "Glimmer_FNR",
                                               "GeneMarkS_F1", "GeneMarkS_Precision", "GeneMarkS_Sensitivity", "GeneMarkS_FDR", "GeneMarkS_FNR",
                                               "Phanotate_F1", "Phanotate_Precision", "Phanotate_Sensitivity", "Phanotate_FDR", "Phanotate_FNR", 
                                               "FragGeneScan_F1", "FragGeneScan_Precision", "FragGeneScan_Sensitivity", "FragGeneScan_FDR", "FragGeneScan_FNR", 
                                               "MGA_F1", "MGA_Precision", "MGA_Sensitivity", "MGA_FDR", "MGA_FNR",
                                               "AUGUSTUS_SA_F1", "AUGUSTUS_SA_Precision", "AUGUSTUS_SA_Sensitivity", "AUGUSTUS_SA_FDR", "AUGUSTUS_SA_FNR",
                                               "AUGUSTUS_EC_F1", "AUGUSTUS_EC_Precision", "AUGUSTUS_EC_Sensitivity", "AUGUSTUS_EC_FDR", "AUGUSTUS_EC_FNR",
                                               "AUGUSTUS_HS_F1", "AUGUSTUS_HS_Precision", "AUGUSTUS_HS_Sensitivity", "AUGUSTUS_HS_FDR", "AUGUSTUS_HS_FNR"),
                                  p_value = listpvalues_shapiro_arc) %>%
  select(-d) %>%
  pivot_wider(names_from = variable, values_from = p_value)

#All p-values (or near all) were significant, so there is no chance to assume normality in our data
#Thus, the only way to compare values among them is via Kruskal-Wallis followed by pairwise tests (Dunn test)

#Measuring the differences between F1 values
datastatsreparc <- subsamples_sum_arc %>%
  select(c(1, 32:length(subsamples_sum_arc))) %>%
  gather(key = var_name, value = value, 2:51) %>%
  separate(col = var_name, sep = "_", into = c("Program", "Statistic")) %>%
  select(-subsample_number)

F1_datarep_arc <- datastatsreparc %>% filter(Statistic == "F1")
Precision_datarep_arc <- datastatsreparc %>% filter(Statistic == "Precision")
Sensitivity_datarep_arc <- datastatsreparc %>% filter(Statistic == "Sensitivity")
FDR_datarep_arc <- datastatsreparc %>% filter(Statistic == "FDR")
FNR_datarep_arc <- datastatsreparc %>% filter(Statistic == "FNR")

kruskal.test(x = F1_datarep_arc$value, g = F1_datarep_arc$Program)
dunnTest(value ~ Program, data = F1_datarep_arc, method = "holm")

kruskal.test(x = Precision_datarep_arc$value, g = Precision_datarep_arc$Program)
dunnTest(value ~ Program, data = Precision_datarep_arc, method = "holm")

kruskal.test(x = Sensitivity_datarep_arc$value, g = Sensitivity_datarep_arc$Program)
dunnTest(value ~ Program, data = Sensitivity_datarep_arc, method = "holm")

kruskal.test(x = FDR_datarep_arc$value, g = FDR_datarep_arc$Program)
dunnTest(value ~ Program, data = FDR_datarep_arc, method = "holm")

kruskal.test(x = FNR_datarep_arc$value, g = FNR_datarep_arc$Program)
dunnTest(value ~ Program, data = FNR_datarep_arc, method = "holm")

# DNA viruses (Real)

sum_dna <- Basal_table %>%
  filter(Nucleic_acid == "DNA") %>%
  summarize_if(is.numeric, sum, na.rm=TRUE) %>%
  rowwise() %>%
  mutate(Prodigal_F1 = (2 * TP_Prodigal) / (2 * TP_Prodigal + FP_Prodigal + FN_Prodigal),
         Prodigal_Precision = TP_Prodigal / (TP_Prodigal + FP_Prodigal),
         Prodigal_Sensitivity = TP_Prodigal / (TP_Prodigal + FN_Prodigal),
         Prodigal_FDR = FP_Prodigal / (TP_Prodigal + FP_Prodigal),
         Prodigal_FNR = FN_Prodigal / (TP_Prodigal + FN_Prodigal),
         MetaProdigal_F1 = (2 * TP_MetaProdigal) / (2 * TP_MetaProdigal + FP_MetaProdigal + FN_MetaProdigal),
         MetaProdigal_Precision = TP_MetaProdigal / (TP_MetaProdigal + FP_MetaProdigal),
         MetaProdigal_Sensitivity = TP_MetaProdigal / (TP_MetaProdigal + FN_MetaProdigal),
         MetaProdigal_FDR = FP_MetaProdigal / (TP_MetaProdigal + FP_MetaProdigal),
         MetaProdigal_FNR = FN_MetaProdigal / (TP_MetaProdigal + FN_MetaProdigal),
         Glimmer_F1 = (2 * TP_Glimmer) / (2 * TP_Glimmer + FP_Glimmer + FN_Glimmer),
         Glimmer_Precision = TP_Glimmer / (TP_Glimmer + FP_Glimmer),
         Glimmer_Sensitivity = TP_Glimmer / (TP_Glimmer + FN_Glimmer),
         Glimmer_FDR = FP_Glimmer / (TP_Glimmer + FP_Glimmer),
         Glimmer_FNR = FN_Glimmer / (TP_Glimmer + FN_Glimmer),
         GeneMarkS_F1 = (2 * TP_GeneMarkS) / (2 * TP_GeneMarkS + FP_GeneMarkS + FN_GeneMarkS),
         GeneMarkS_Precision = TP_GeneMarkS / (TP_GeneMarkS + FP_GeneMarkS),
         GeneMarkS_Sensitivity = TP_GeneMarkS / (TP_GeneMarkS + FN_GeneMarkS),
         GeneMarkS_FDR = FP_GeneMarkS / (TP_GeneMarkS + FP_GeneMarkS),
         GeneMarkS_FNR = FN_GeneMarkS / (TP_GeneMarkS + FN_GeneMarkS),
         Phanotate_F1 = (2 * TP_Phanotate) / (2 * TP_Phanotate + FP_Phanotate + FN_Phanotate),
         Phanotate_Precision = TP_Phanotate / (TP_Phanotate + FP_Phanotate),
         Phanotate_Sensitivity = TP_Phanotate / (TP_Phanotate + FN_Phanotate),
         Phanotate_FDR = FP_Phanotate / (TP_Phanotate + FP_Phanotate),
         Phanotate_FNR = FN_Phanotate / (TP_Phanotate + FN_Phanotate),
         FragGeneScan_F1 = (2 * TP_FragGeneScan) / (2 * TP_FragGeneScan + FP_FragGeneScan + FN_FragGeneScan),
         FragGeneScan_Precision = TP_FragGeneScan / (TP_FragGeneScan + FP_FragGeneScan),
         FragGeneScan_Sensitivity = TP_FragGeneScan / (TP_FragGeneScan + FN_FragGeneScan),
         FragGeneScan_FDR = FP_FragGeneScan / (TP_FragGeneScan + FP_FragGeneScan),
         FragGeneScan_FNR = FN_FragGeneScan / (TP_FragGeneScan + FN_FragGeneScan),
         MGA_F1 = (2 * TP_MGA) / (2 * TP_MGA + FP_MGA + FN_MGA),
         MGA_Precision = TP_MGA / (TP_MGA + FP_MGA),
         MGA_Sensitivity = TP_MGA / (TP_MGA + FN_MGA),
         MGA_FDR = FP_MGA / (TP_MGA + FP_MGA),
         MGA_FNR = FN_MGA / (TP_MGA + FN_MGA),
         AUGUSTUS_SA_F1 = (2 * TP_AUGUSTUS_SA) / (2 * TP_AUGUSTUS_SA + FP_AUGUSTUS_SA + FN_AUGUSTUS_SA),
         AUGUSTUS_SA_Precision = TP_AUGUSTUS_SA / (TP_AUGUSTUS_SA + FP_AUGUSTUS_SA),
         AUGUSTUS_SA_Sensitivity = TP_AUGUSTUS_SA / (TP_AUGUSTUS_SA + FN_AUGUSTUS_SA),
         AUGUSTUS_SA_FDR = FP_AUGUSTUS_SA / (TP_AUGUSTUS_SA + FP_AUGUSTUS_SA),
         AUGUSTUS_SA_FNR = FN_AUGUSTUS_SA / (TP_AUGUSTUS_SA + FN_AUGUSTUS_SA),
         AUGUSTUS_EC_F1 = (2 * TP_AUGUSTUS_EC) / (2 * TP_AUGUSTUS_EC + FP_AUGUSTUS_EC + FN_AUGUSTUS_EC),
         AUGUSTUS_EC_Precision = TP_AUGUSTUS_EC / (TP_AUGUSTUS_EC + FP_AUGUSTUS_EC),
         AUGUSTUS_EC_Sensitivity = TP_AUGUSTUS_EC / (TP_AUGUSTUS_EC + FN_AUGUSTUS_EC),
         AUGUSTUS_EC_FDR = FP_AUGUSTUS_EC / (TP_AUGUSTUS_EC + FP_AUGUSTUS_EC),
         AUGUSTUS_EC_FNR = FN_AUGUSTUS_EC / (TP_AUGUSTUS_EC + FN_AUGUSTUS_EC),
         AUGUSTUS_HS_F1 = (2 * TP_AUGUSTUS_HS) / (2 * TP_AUGUSTUS_HS + FP_AUGUSTUS_HS + FN_AUGUSTUS_HS),
         AUGUSTUS_HS_Precision = TP_AUGUSTUS_HS / (TP_AUGUSTUS_HS + FP_AUGUSTUS_HS),
         AUGUSTUS_HS_Sensitivity = TP_AUGUSTUS_HS / (TP_AUGUSTUS_HS + FN_AUGUSTUS_HS),
         AUGUSTUS_HS_FDR = FP_AUGUSTUS_HS / (TP_AUGUSTUS_HS + FP_AUGUSTUS_HS),
         AUGUSTUS_HS_FNR = FN_AUGUSTUS_HS / (TP_AUGUSTUS_HS + FN_AUGUSTUS_HS))

# DNA (Subsampling)

listofdfs = list()
for(i in 1:1000){
  listofdfs[[i]] <- Basal_table %>%
    filter(Nucleic_acid == "DNA") %>%
    slice_sample(prop = 0.1)
}
dfofdfs <- tibble(d = 1:1000, data = listofdfs)

sum_dna_subsamples <- list()
for(i in 1:1000){
  sum_dna_subsamples[[i]] <- dfofdfs$data[[i]] %>%
    summarize_if(is.numeric, sum, na.rm=TRUE)
}

subsamples_sum_dna <- bind_rows(sum_dna_subsamples, .id = "subsample_number") %>%
  rowwise() %>%
  mutate(Prodigal_F1 = (2 * TP_Prodigal) / (2 * TP_Prodigal + FP_Prodigal + FN_Prodigal),
         Prodigal_Precision = TP_Prodigal / (TP_Prodigal + FP_Prodigal),
         Prodigal_Sensitivity = TP_Prodigal / (TP_Prodigal + FN_Prodigal),
         Prodigal_FDR = FP_Prodigal / (TP_Prodigal + FP_Prodigal),
         Prodigal_FNR = FN_Prodigal / (TP_Prodigal + FN_Prodigal),
         MetaProdigal_F1 = (2 * TP_MetaProdigal) / (2 * TP_MetaProdigal + FP_MetaProdigal + FN_MetaProdigal),
         MetaProdigal_Precision = TP_MetaProdigal / (TP_MetaProdigal + FP_MetaProdigal),
         MetaProdigal_Sensitivity = TP_MetaProdigal / (TP_MetaProdigal + FN_MetaProdigal),
         MetaProdigal_FDR = FP_MetaProdigal / (TP_MetaProdigal + FP_MetaProdigal),
         MetaProdigal_FNR = FN_MetaProdigal / (TP_MetaProdigal + FN_MetaProdigal),
         Glimmer_F1 = (2 * TP_Glimmer) / (2 * TP_Glimmer + FP_Glimmer + FN_Glimmer),
         Glimmer_Precision = TP_Glimmer / (TP_Glimmer + FP_Glimmer),
         Glimmer_Sensitivity = TP_Glimmer / (TP_Glimmer + FN_Glimmer),
         Glimmer_FDR = FP_Glimmer / (TP_Glimmer + FP_Glimmer),
         Glimmer_FNR = FN_Glimmer / (TP_Glimmer + FN_Glimmer),
         GeneMarkS_F1 = (2 * TP_GeneMarkS) / (2 * TP_GeneMarkS + FP_GeneMarkS + FN_GeneMarkS),
         GeneMarkS_Precision = TP_GeneMarkS / (TP_GeneMarkS + FP_GeneMarkS),
         GeneMarkS_Sensitivity = TP_GeneMarkS / (TP_GeneMarkS + FN_GeneMarkS),
         GeneMarkS_FDR = FP_GeneMarkS / (TP_GeneMarkS + FP_GeneMarkS),
         GeneMarkS_FNR = FN_GeneMarkS / (TP_GeneMarkS + FN_GeneMarkS),
         Phanotate_F1 = (2 * TP_Phanotate) / (2 * TP_Phanotate + FP_Phanotate + FN_Phanotate),
         Phanotate_Precision = TP_Phanotate / (TP_Phanotate + FP_Phanotate),
         Phanotate_Sensitivity = TP_Phanotate / (TP_Phanotate + FN_Phanotate),
         Phanotate_FDR = FP_Phanotate / (TP_Phanotate + FP_Phanotate),
         Phanotate_FNR = FN_Phanotate / (TP_Phanotate + FN_Phanotate),
         FragGeneScan_F1 = (2 * TP_FragGeneScan) / (2 * TP_FragGeneScan + FP_FragGeneScan + FN_FragGeneScan),
         FragGeneScan_Precision = TP_FragGeneScan / (TP_FragGeneScan + FP_FragGeneScan),
         FragGeneScan_Sensitivity = TP_FragGeneScan / (TP_FragGeneScan + FN_FragGeneScan),
         FragGeneScan_FDR = FP_FragGeneScan / (TP_FragGeneScan + FP_FragGeneScan),
         FragGeneScan_FNR = FN_FragGeneScan / (TP_FragGeneScan + FN_FragGeneScan),
         MGA_F1 = (2 * TP_MGA) / (2 * TP_MGA + FP_MGA + FN_MGA),
         MGA_Precision = TP_MGA / (TP_MGA + FP_MGA),
         MGA_Sensitivity = TP_MGA / (TP_MGA + FN_MGA),
         MGA_FDR = FP_MGA / (TP_MGA + FP_MGA),
         MGA_FNR = FN_MGA / (TP_MGA + FN_MGA),
         AUGUSTUSSA_F1 = (2 * TP_AUGUSTUS_SA) / (2 * TP_AUGUSTUS_SA + FP_AUGUSTUS_SA + FN_AUGUSTUS_SA),
         AUGUSTUSSA_Precision = TP_AUGUSTUS_SA / (TP_AUGUSTUS_SA + FP_AUGUSTUS_SA),
         AUGUSTUSSA_Sensitivity = TP_AUGUSTUS_SA / (TP_AUGUSTUS_SA + FN_AUGUSTUS_SA),
         AUGUSTUSSA_FDR = FP_AUGUSTUS_SA / (TP_AUGUSTUS_SA + FP_AUGUSTUS_SA),
         AUGUSTUSSA_FNR = FN_AUGUSTUS_SA / (TP_AUGUSTUS_SA + FN_AUGUSTUS_SA),
         AUGUSTUSEC_F1 = (2 * TP_AUGUSTUS_EC) / (2 * TP_AUGUSTUS_EC + FP_AUGUSTUS_EC + FN_AUGUSTUS_EC),
         AUGUSTUSEC_Precision = TP_AUGUSTUS_EC / (TP_AUGUSTUS_EC + FP_AUGUSTUS_EC),
         AUGUSTUSEC_Sensitivity = TP_AUGUSTUS_EC / (TP_AUGUSTUS_EC + FN_AUGUSTUS_EC),
         AUGUSTUSEC_FDR = FP_AUGUSTUS_EC / (TP_AUGUSTUS_EC + FP_AUGUSTUS_EC),
         AUGUSTUSEC_FNR = FN_AUGUSTUS_EC / (TP_AUGUSTUS_EC + FN_AUGUSTUS_EC),
         AUGUSTUSHS_F1 = (2 * TP_AUGUSTUS_HS) / (2 * TP_AUGUSTUS_HS + FP_AUGUSTUS_HS + FN_AUGUSTUS_HS),
         AUGUSTUSHS_Precision = TP_AUGUSTUS_HS / (TP_AUGUSTUS_HS + FP_AUGUSTUS_HS),
         AUGUSTUSHS_Sensitivity = TP_AUGUSTUS_HS / (TP_AUGUSTUS_HS + FN_AUGUSTUS_HS),
         AUGUSTUSHS_FDR = FP_AUGUSTUS_HS / (TP_AUGUSTUS_HS + FP_AUGUSTUS_HS),
         AUGUSTUSHS_FNR = FN_AUGUSTUS_HS / (TP_AUGUSTUS_HS + FN_AUGUSTUS_HS))

#Testing normality for the variables F1, Precision, Sensitivity, FDR and FNR of the different programs
listpvalues_shapiro_dna = list()

for(i in 32:length(subsamples_sum_dna)){
  listpvalues_shapiro_dna[[i]] <- shapiro.test(subsamples_sum_dna[[i]])$p.value
}

listpvalues_shapiro_dna <- Filter(Negate(is.null), listpvalues_shapiro_dna)
listpvalues_shapiro_dna <- as.list(p.adjust(listpvalues_shapiro_dna, method = "fdr"))
datapvalues_shapiro_dna <- tibble(d = 1:length(listpvalues_shapiro_dna),
                                  variable = c("Prodigal_F1", "Prodigal_Precision", "Prodigal_Sensitivity", "Prodigal_FDR", "Prodigal_FNR",
                                               "MetaProdigal_F1", "MetaProdigal_Precision", "MetaProdigal_Sensitivity", "MetaProdigal_FDR", "MetaProdigal_FNR",
                                               "Glimmer_F1", "Glimmer_Precision", "Glimmer_Sensitivity", "Glimmer_FDR", "Glimmer_FNR",
                                               "GeneMarkS_F1", "GeneMarkS_Precision", "GeneMarkS_Sensitivity", "GeneMarkS_FDR", "GeneMarkS_FNR",
                                               "Phanotate_F1", "Phanotate_Precision", "Phanotate_Sensitivity", "Phanotate_FDR", "Phanotate_FNR", 
                                               "FragGeneScan_F1", "FragGeneScan_Precision", "FragGeneScan_Sensitivity", "FragGeneScan_FDR", "FragGeneScan_FNR", 
                                               "MGA_F1", "MGA_Precision", "MGA_Sensitivity", "MGA_FDR", "MGA_FNR",
                                               "AUGUSTUS_SA_F1", "AUGUSTUS_SA_Precision", "AUGUSTUS_SA_Sensitivity", "AUGUSTUS_SA_FDR", "AUGUSTUS_SA_FNR",
                                               "AUGUSTUS_EC_F1", "AUGUSTUS_EC_Precision", "AUGUSTUS_EC_Sensitivity", "AUGUSTUS_EC_FDR", "AUGUSTUS_EC_FNR",
                                               "AUGUSTUS_HS_F1", "AUGUSTUS_HS_Precision", "AUGUSTUS_HS_Sensitivity", "AUGUSTUS_HS_FDR", "AUGUSTUS_HS_FNR"),
                                  p_value = listpvalues_shapiro_dna) %>%
  select(-d) %>%
  pivot_wider(names_from = variable, values_from = p_value)

#All p-values (or near all) were significant, so there is no chance to assume normality in our data
#Thus, the only way to compare values among them is via Kruskal-Wallis followed by pairwise tests (Dunn test)

#Measuring the differences between F1 values
datastatsrepdna <- subsamples_sum_dna %>%
  select(c(1, 32:length(subsamples_sum_dna))) %>%
  gather(key = var_name, value = value, 2:51) %>%
  separate(col = var_name, sep = "_", into = c("Program", "Statistic")) %>%
  select(-subsample_number)

F1_datarep_dna <- datastatsrepdna %>% filter(Statistic == "F1")
Precision_datarep_dna <- datastatsrepdna %>% filter(Statistic == "Precision")
Sensitivity_datarep_dna <- datastatsrepdna %>% filter(Statistic == "Sensitivity")
FDR_datarep_dna <- datastatsrepdna %>% filter(Statistic == "FDR")
FNR_datarep_dna <- datastatsrepdna %>% filter(Statistic == "FNR")

kruskal.test(x = F1_datarep_dna$value, g = F1_datarep_dna$Program)
dunnTest(value ~ Program, data = F1_datarep_dna, method = "holm")

kruskal.test(x = Precision_datarep_dna$value, g = Precision_datarep_dna$Program)
dunnTest(value ~ Program, data = Precision_datarep_dna, method = "holm")

kruskal.test(x = Sensitivity_datarep_dna$value, g = Sensitivity_datarep_dna$Program)
dunnTest(value ~ Program, data = Sensitivity_datarep_dna, method = "holm")

kruskal.test(x = FDR_datarep_dna$value, g = FDR_datarep_dna$Program)
dunnTest(value ~ Program, data = FDR_datarep_dna, method = "holm")

kruskal.test(x = FNR_datarep_dna$value, g = FNR_datarep_dna$Program)
dunnTest(value ~ Program, data = FNR_datarep_dna, method = "holm")

# RNA viruses (Real)

sum_rna <- Basal_table %>%
  filter(Nucleic_acid == "RNA") %>%
  summarize_if(is.numeric, sum, na.rm=TRUE) %>%
  rowwise() %>%
  mutate(Prodigal_F1 = (2 * TP_Prodigal) / (2 * TP_Prodigal + FP_Prodigal + FN_Prodigal),
         Prodigal_Precision = TP_Prodigal / (TP_Prodigal + FP_Prodigal),
         Prodigal_Sensitivity = TP_Prodigal / (TP_Prodigal + FN_Prodigal),
         Prodigal_FDR = FP_Prodigal / (TP_Prodigal + FP_Prodigal),
         Prodigal_FNR = FN_Prodigal / (TP_Prodigal + FN_Prodigal),
         MetaProdigal_F1 = (2 * TP_MetaProdigal) / (2 * TP_MetaProdigal + FP_MetaProdigal + FN_MetaProdigal),
         MetaProdigal_Precision = TP_MetaProdigal / (TP_MetaProdigal + FP_MetaProdigal),
         MetaProdigal_Sensitivity = TP_MetaProdigal / (TP_MetaProdigal + FN_MetaProdigal),
         MetaProdigal_FDR = FP_MetaProdigal / (TP_MetaProdigal + FP_MetaProdigal),
         MetaProdigal_FNR = FN_MetaProdigal / (TP_MetaProdigal + FN_MetaProdigal),
         Glimmer_F1 = (2 * TP_Glimmer) / (2 * TP_Glimmer + FP_Glimmer + FN_Glimmer),
         Glimmer_Precision = TP_Glimmer / (TP_Glimmer + FP_Glimmer),
         Glimmer_Sensitivity = TP_Glimmer / (TP_Glimmer + FN_Glimmer),
         Glimmer_FDR = FP_Glimmer / (TP_Glimmer + FP_Glimmer),
         Glimmer_FNR = FN_Glimmer / (TP_Glimmer + FN_Glimmer),
         GeneMarkS_F1 = (2 * TP_GeneMarkS) / (2 * TP_GeneMarkS + FP_GeneMarkS + FN_GeneMarkS),
         GeneMarkS_Precision = TP_GeneMarkS / (TP_GeneMarkS + FP_GeneMarkS),
         GeneMarkS_Sensitivity = TP_GeneMarkS / (TP_GeneMarkS + FN_GeneMarkS),
         GeneMarkS_FDR = FP_GeneMarkS / (TP_GeneMarkS + FP_GeneMarkS),
         GeneMarkS_FNR = FN_GeneMarkS / (TP_GeneMarkS + FN_GeneMarkS),
         Phanotate_F1 = (2 * TP_Phanotate) / (2 * TP_Phanotate + FP_Phanotate + FN_Phanotate),
         Phanotate_Precision = TP_Phanotate / (TP_Phanotate + FP_Phanotate),
         Phanotate_Sensitivity = TP_Phanotate / (TP_Phanotate + FN_Phanotate),
         Phanotate_FDR = FP_Phanotate / (TP_Phanotate + FP_Phanotate),
         Phanotate_FNR = FN_Phanotate / (TP_Phanotate + FN_Phanotate),
         FragGeneScan_F1 = (2 * TP_FragGeneScan) / (2 * TP_FragGeneScan + FP_FragGeneScan + FN_FragGeneScan),
         FragGeneScan_Precision = TP_FragGeneScan / (TP_FragGeneScan + FP_FragGeneScan),
         FragGeneScan_Sensitivity = TP_FragGeneScan / (TP_FragGeneScan + FN_FragGeneScan),
         FragGeneScan_FDR = FP_FragGeneScan / (TP_FragGeneScan + FP_FragGeneScan),
         FragGeneScan_FNR = FN_FragGeneScan / (TP_FragGeneScan + FN_FragGeneScan),
         MGA_F1 = (2 * TP_MGA) / (2 * TP_MGA + FP_MGA + FN_MGA),
         MGA_Precision = TP_MGA / (TP_MGA + FP_MGA),
         MGA_Sensitivity = TP_MGA / (TP_MGA + FN_MGA),
         MGA_FDR = FP_MGA / (TP_MGA + FP_MGA),
         MGA_FNR = FN_MGA / (TP_MGA + FN_MGA),
         AUGUSTUS_SA_F1 = (2 * TP_AUGUSTUS_SA) / (2 * TP_AUGUSTUS_SA + FP_AUGUSTUS_SA + FN_AUGUSTUS_SA),
         AUGUSTUS_SA_Precision = TP_AUGUSTUS_SA / (TP_AUGUSTUS_SA + FP_AUGUSTUS_SA),
         AUGUSTUS_SA_Sensitivity = TP_AUGUSTUS_SA / (TP_AUGUSTUS_SA + FN_AUGUSTUS_SA),
         AUGUSTUS_SA_FDR = FP_AUGUSTUS_SA / (TP_AUGUSTUS_SA + FP_AUGUSTUS_SA),
         AUGUSTUS_SA_FNR = FN_AUGUSTUS_SA / (TP_AUGUSTUS_SA + FN_AUGUSTUS_SA),
         AUGUSTUS_EC_F1 = (2 * TP_AUGUSTUS_EC) / (2 * TP_AUGUSTUS_EC + FP_AUGUSTUS_EC + FN_AUGUSTUS_EC),
         AUGUSTUS_EC_Precision = TP_AUGUSTUS_EC / (TP_AUGUSTUS_EC + FP_AUGUSTUS_EC),
         AUGUSTUS_EC_Sensitivity = TP_AUGUSTUS_EC / (TP_AUGUSTUS_EC + FN_AUGUSTUS_EC),
         AUGUSTUS_EC_FDR = FP_AUGUSTUS_EC / (TP_AUGUSTUS_EC + FP_AUGUSTUS_EC),
         AUGUSTUS_EC_FNR = FN_AUGUSTUS_EC / (TP_AUGUSTUS_EC + FN_AUGUSTUS_EC),
         AUGUSTUS_HS_F1 = (2 * TP_AUGUSTUS_HS) / (2 * TP_AUGUSTUS_HS + FP_AUGUSTUS_HS + FN_AUGUSTUS_HS),
         AUGUSTUS_HS_Precision = TP_AUGUSTUS_HS / (TP_AUGUSTUS_HS + FP_AUGUSTUS_HS),
         AUGUSTUS_HS_Sensitivity = TP_AUGUSTUS_HS / (TP_AUGUSTUS_HS + FN_AUGUSTUS_HS),
         AUGUSTUS_HS_FDR = FP_AUGUSTUS_HS / (TP_AUGUSTUS_HS + FP_AUGUSTUS_HS),
         AUGUSTUS_HS_FNR = FN_AUGUSTUS_HS / (TP_AUGUSTUS_HS + FN_AUGUSTUS_HS))

# RNA (Subsampling)

listofdfs = list()
for(i in 1:1000){
  listofdfs[[i]] <- Basal_table %>%
    filter(Nucleic_acid == "RNA") %>%
    slice_sample(prop = 0.1)
}
dfofdfs <- tibble(d = 1:1000, data = listofdfs)

sum_rna_subsamples <- list()
for(i in 1:1000){
  sum_rna_subsamples[[i]] <- dfofdfs$data[[i]] %>%
    summarize_if(is.numeric, sum, na.rm=TRUE)
}

subsamples_sum_rna <- bind_rows(sum_rna_subsamples, .id = "subsample_number") %>%
  rowwise() %>%
  mutate(Prodigal_F1 = (2 * TP_Prodigal) / (2 * TP_Prodigal + FP_Prodigal + FN_Prodigal),
         Prodigal_Precision = TP_Prodigal / (TP_Prodigal + FP_Prodigal),
         Prodigal_Sensitivity = TP_Prodigal / (TP_Prodigal + FN_Prodigal),
         Prodigal_FDR = FP_Prodigal / (TP_Prodigal + FP_Prodigal),
         Prodigal_FNR = FN_Prodigal / (TP_Prodigal + FN_Prodigal),
         MetaProdigal_F1 = (2 * TP_MetaProdigal) / (2 * TP_MetaProdigal + FP_MetaProdigal + FN_MetaProdigal),
         MetaProdigal_Precision = TP_MetaProdigal / (TP_MetaProdigal + FP_MetaProdigal),
         MetaProdigal_Sensitivity = TP_MetaProdigal / (TP_MetaProdigal + FN_MetaProdigal),
         MetaProdigal_FDR = FP_MetaProdigal / (TP_MetaProdigal + FP_MetaProdigal),
         MetaProdigal_FNR = FN_MetaProdigal / (TP_MetaProdigal + FN_MetaProdigal),
         Glimmer_F1 = (2 * TP_Glimmer) / (2 * TP_Glimmer + FP_Glimmer + FN_Glimmer),
         Glimmer_Precision = TP_Glimmer / (TP_Glimmer + FP_Glimmer),
         Glimmer_Sensitivity = TP_Glimmer / (TP_Glimmer + FN_Glimmer),
         Glimmer_FDR = FP_Glimmer / (TP_Glimmer + FP_Glimmer),
         Glimmer_FNR = FN_Glimmer / (TP_Glimmer + FN_Glimmer),
         GeneMarkS_F1 = (2 * TP_GeneMarkS) / (2 * TP_GeneMarkS + FP_GeneMarkS + FN_GeneMarkS),
         GeneMarkS_Precision = TP_GeneMarkS / (TP_GeneMarkS + FP_GeneMarkS),
         GeneMarkS_Sensitivity = TP_GeneMarkS / (TP_GeneMarkS + FN_GeneMarkS),
         GeneMarkS_FDR = FP_GeneMarkS / (TP_GeneMarkS + FP_GeneMarkS),
         GeneMarkS_FNR = FN_GeneMarkS / (TP_GeneMarkS + FN_GeneMarkS),
         Phanotate_F1 = (2 * TP_Phanotate) / (2 * TP_Phanotate + FP_Phanotate + FN_Phanotate),
         Phanotate_Precision = TP_Phanotate / (TP_Phanotate + FP_Phanotate),
         Phanotate_Sensitivity = TP_Phanotate / (TP_Phanotate + FN_Phanotate),
         Phanotate_FDR = FP_Phanotate / (TP_Phanotate + FP_Phanotate),
         Phanotate_FNR = FN_Phanotate / (TP_Phanotate + FN_Phanotate),
         FragGeneScan_F1 = (2 * TP_FragGeneScan) / (2 * TP_FragGeneScan + FP_FragGeneScan + FN_FragGeneScan),
         FragGeneScan_Precision = TP_FragGeneScan / (TP_FragGeneScan + FP_FragGeneScan),
         FragGeneScan_Sensitivity = TP_FragGeneScan / (TP_FragGeneScan + FN_FragGeneScan),
         FragGeneScan_FDR = FP_FragGeneScan / (TP_FragGeneScan + FP_FragGeneScan),
         FragGeneScan_FNR = FN_FragGeneScan / (TP_FragGeneScan + FN_FragGeneScan),
         MGA_F1 = (2 * TP_MGA) / (2 * TP_MGA + FP_MGA + FN_MGA),
         MGA_Precision = TP_MGA / (TP_MGA + FP_MGA),
         MGA_Sensitivity = TP_MGA / (TP_MGA + FN_MGA),
         MGA_FDR = FP_MGA / (TP_MGA + FP_MGA),
         MGA_FNR = FN_MGA / (TP_MGA + FN_MGA),
         AUGUSTUSSA_F1 = (2 * TP_AUGUSTUS_SA) / (2 * TP_AUGUSTUS_SA + FP_AUGUSTUS_SA + FN_AUGUSTUS_SA),
         AUGUSTUSSA_Precision = TP_AUGUSTUS_SA / (TP_AUGUSTUS_SA + FP_AUGUSTUS_SA),
         AUGUSTUSSA_Sensitivity = TP_AUGUSTUS_SA / (TP_AUGUSTUS_SA + FN_AUGUSTUS_SA),
         AUGUSTUSSA_FDR = FP_AUGUSTUS_SA / (TP_AUGUSTUS_SA + FP_AUGUSTUS_SA),
         AUGUSTUSSA_FNR = FN_AUGUSTUS_SA / (TP_AUGUSTUS_SA + FN_AUGUSTUS_SA),
         AUGUSTUSEC_F1 = (2 * TP_AUGUSTUS_EC) / (2 * TP_AUGUSTUS_EC + FP_AUGUSTUS_EC + FN_AUGUSTUS_EC),
         AUGUSTUSEC_Precision = TP_AUGUSTUS_EC / (TP_AUGUSTUS_EC + FP_AUGUSTUS_EC),
         AUGUSTUSEC_Sensitivity = TP_AUGUSTUS_EC / (TP_AUGUSTUS_EC + FN_AUGUSTUS_EC),
         AUGUSTUSEC_FDR = FP_AUGUSTUS_EC / (TP_AUGUSTUS_EC + FP_AUGUSTUS_EC),
         AUGUSTUSEC_FNR = FN_AUGUSTUS_EC / (TP_AUGUSTUS_EC + FN_AUGUSTUS_EC),
         AUGUSTUSHS_F1 = (2 * TP_AUGUSTUS_HS) / (2 * TP_AUGUSTUS_HS + FP_AUGUSTUS_HS + FN_AUGUSTUS_HS),
         AUGUSTUSHS_Precision = TP_AUGUSTUS_HS / (TP_AUGUSTUS_HS + FP_AUGUSTUS_HS),
         AUGUSTUSHS_Sensitivity = TP_AUGUSTUS_HS / (TP_AUGUSTUS_HS + FN_AUGUSTUS_HS),
         AUGUSTUSHS_FDR = FP_AUGUSTUS_HS / (TP_AUGUSTUS_HS + FP_AUGUSTUS_HS),
         AUGUSTUSHS_FNR = FN_AUGUSTUS_HS / (TP_AUGUSTUS_HS + FN_AUGUSTUS_HS))

#Testing normality for the variables F1, Precision, Sensitivity, FDR and FNR of the different programs
listpvalues_shapiro_rna = list()

for(i in 32:length(subsamples_sum_rna)){
  listpvalues_shapiro_rna[[i]] <- shapiro.test(subsamples_sum_rna[[i]])$p.value
}

listpvalues_shapiro_rna <- Filter(Negate(is.null), listpvalues_shapiro_rna)
listpvalues_shapiro_rna <- as.list(p.adjust(listpvalues_shapiro_rna, method = "fdr"))
datapvalues_shapiro_rna <- tibble(d = 1:length(listpvalues_shapiro_rna),
                                  variable = c("Prodigal_F1", "Prodigal_Precision", "Prodigal_Sensitivity", "Prodigal_FDR", "Prodigal_FNR",
                                               "MetaProdigal_F1", "MetaProdigal_Precision", "MetaProdigal_Sensitivity", "MetaProdigal_FDR", "MetaProdigal_FNR",
                                               "Glimmer_F1", "Glimmer_Precision", "Glimmer_Sensitivity", "Glimmer_FDR", "Glimmer_FNR",
                                               "GeneMarkS_F1", "GeneMarkS_Precision", "GeneMarkS_Sensitivity", "GeneMarkS_FDR", "GeneMarkS_FNR",
                                               "Phanotate_F1", "Phanotate_Precision", "Phanotate_Sensitivity", "Phanotate_FDR", "Phanotate_FNR", 
                                               "FragGeneScan_F1", "FragGeneScan_Precision", "FragGeneScan_Sensitivity", "FragGeneScan_FDR", "FragGeneScan_FNR", 
                                               "MGA_F1", "MGA_Precision", "MGA_Sensitivity", "MGA_FDR", "MGA_FNR",
                                               "AUGUSTUS_SA_F1", "AUGUSTUS_SA_Precision", "AUGUSTUS_SA_Sensitivity", "AUGUSTUS_SA_FDR", "AUGUSTUS_SA_FNR",
                                               "AUGUSTUS_EC_F1", "AUGUSTUS_EC_Precision", "AUGUSTUS_EC_Sensitivity", "AUGUSTUS_EC_FDR", "AUGUSTUS_EC_FNR",
                                               "AUGUSTUS_HS_F1", "AUGUSTUS_HS_Precision", "AUGUSTUS_HS_Sensitivity", "AUGUSTUS_HS_FDR", "AUGUSTUS_HS_FNR"),
                                  p_value = listpvalues_shapiro_rna) %>%
  select(-d) %>%
  pivot_wider(names_from = variable, values_from = p_value)

#All p-values (or near all) were significant, so there is no chance to assume normality in our data
#Thus, the only way to compare values among them is via Kruskal-Wallis followed by pairwise tests (Dunn test)

#Measuring the differences between F1 values
datastatsreprna <- subsamples_sum_rna %>%
  select(c(1, 32:length(subsamples_sum_rna))) %>%
  gather(key = var_name, value = value, 2:51) %>%
  separate(col = var_name, sep = "_", into = c("Program", "Statistic")) %>%
  select(-subsample_number)

F1_datarep_rna <- datastatsreprna %>% filter(Statistic == "F1")
Precision_datarep_rna <- datastatsreprna %>% filter(Statistic == "Precision")
Sensitivity_datarep_rna <- datastatsreprna %>% filter(Statistic == "Sensitivity")
FDR_datarep_rna <- datastatsreprna %>% filter(Statistic == "FDR")
FNR_datarep_rna <- datastatsreprna %>% filter(Statistic == "FNR")

kruskal.test(x = F1_datarep_rna$value, g = F1_datarep_rna$Program)
dunnTest(value ~ Program, data = F1_datarep_rna, method = "holm")

kruskal.test(x = Precision_datarep_rna$value, g = Precision_datarep_rna$Program)
dunnTest(value ~ Program, data = Precision_datarep_rna, method = "holm")

kruskal.test(x = Sensitivity_datarep_rna$value, g = Sensitivity_datarep_rna$Program)
dunnTest(value ~ Program, data = Sensitivity_datarep_rna, method = "holm")

kruskal.test(x = FDR_datarep_rna$value, g = FDR_datarep_rna$Program)
dunnTest(value ~ Program, data = FDR_datarep_rna, method = "holm")

kruskal.test(x = FNR_datarep_rna$value, g = FNR_datarep_rna$Program)
dunnTest(value ~ Program, data = FNR_datarep_rna, method = "holm")





%>% 
  spread(key = names(subsamples_sum_all)[1],value = 'value') 




  
  


# dunnTest(flipper_length_mm ~ species,
#          data = dat,
#          method = "holm"
# )
# 
# subsamples_sum_all$...
# Prodigal_F1 = (2 * TP_Prodigal) / (2 * TP_Prodigal + FP_Prodigal + FN_Prodigal),
# Prodigal_Precision = TP_Prodigal / (TP_Prodigal + FP_Prodigal),
# Prodigal_Sensitivity = TP_Prodigal / (TP_Prodigal + FN_Prodigal),
# Prodigal_FDR = FP_Prodigal / (TP_Prodigal + FP_Prodigal),
# Prodigal_FNR = FN_Prodigal / (TP_Prodigal + FN_Prodigal),
# MetaProdigal_F1 = (2 * TP_MetaProdigal) / (2 * TP_MetaProdigal + FP_MetaProdigal + FN_MetaProdigal),
# MetaProdigal_Precision = TP_MetaProdigal / (TP_MetaProdigal + FP_MetaProdigal),
# MetaProdigal_Sensitivity = TP_MetaProdigal / (TP_MetaProdigal + FN_MetaProdigal),
# MetaProdigal_FDR = FP_MetaProdigal / (TP_MetaProdigal + FP_MetaProdigal),
# MetaProdigal_FNR = FN_MetaProdigal / (TP_MetaProdigal + FN_MetaProdigal),
# Glimmer_F1 = (2 * TP_Glimmer) / (2 * TP_Glimmer + FP_Glimmer + FN_Glimmer),
# Glimmer_Precision = TP_Glimmer / (TP_Glimmer + FP_Glimmer),
# Glimmer_Sensitivity = TP_Glimmer / (TP_Glimmer + FN_Glimmer),
# Glimmer_FDR = FP_Glimmer / (TP_Glimmer + FP_Glimmer),
# Glimmer_FNR = FN_Glimmer / (TP_Glimmer + FN_Glimmer),
# GeneMarkS_F1 = (2 * TP_GeneMarkS) / (2 * TP_GeneMarkS + FP_GeneMarkS + FN_GeneMarkS),
# GeneMarkS_Precision = TP_GeneMarkS / (TP_GeneMarkS + FP_GeneMarkS),
# GeneMarkS_Sensitivity = TP_GeneMarkS / (TP_GeneMarkS + FN_GeneMarkS),
# GeneMarkS_FDR = FP_GeneMarkS / (TP_GeneMarkS + FP_GeneMarkS),
# GeneMarkS_FNR = FN_GeneMarkS / (TP_GeneMarkS + FN_GeneMarkS),
# Phanotate_F1 = (2 * TP_Phanotate) / (2 * TP_Phanotate + FP_Phanotate + FN_Phanotate),
# Phanotate_Precision = TP_Phanotate / (TP_Phanotate + FP_Phanotate),
# Phanotate_Sensitivity = TP_Phanotate / (TP_Phanotate + FN_Phanotate),
# Phanotate_FDR = FP_Phanotate / (TP_Phanotate + FP_Phanotate),
# Phanotate_FNR = FN_Phanotate / (TP_Phanotate + FN_Phanotate),
# FragGeneScan_F1 = (2 * TP_FragGeneScan) / (2 * TP_FragGeneScan + FP_FragGeneScan + FN_FragGeneScan),
# FragGeneScan_Precision = TP_FragGeneScan / (TP_FragGeneScan + FP_FragGeneScan),
# FragGeneScan_Sensitivity = TP_FragGeneScan / (TP_FragGeneScan + FN_FragGeneScan),
# FragGeneScan_FDR = FP_FragGeneScan / (TP_FragGeneScan + FP_FragGeneScan),
# FragGeneScan_FNR = FN_FragGeneScan / (TP_FragGeneScan + FN_FragGeneScan),
# MGA_F1 = (2 * TP_MGA) / (2 * TP_MGA + FP_MGA + FN_MGA),
# MGA_Precision = TP_MGA / (TP_MGA + FP_MGA),
# MGA_Sensitivity = TP_MGA / (TP_MGA + FN_MGA),
# MGA_FDR = FP_MGA / (TP_MGA + FP_MGA),
# MGA_FNR = FN_MGA / (TP_MGA + FN_MGA),
# AUGUSTUS_SA_F1 = (2 * TP_AUGUSTUS_SA) / (2 * TP_AUGUSTUS_SA + FP_AUGUSTUS_SA + FN_AUGUSTUS_SA),
# AUGUSTUS_SA_Precision = TP_AUGUSTUS_SA / (TP_AUGUSTUS_SA + FP_AUGUSTUS_SA),
# AUGUSTUS_SA_Sensitivity = TP_AUGUSTUS_SA / (TP_AUGUSTUS_SA + FN_AUGUSTUS_SA),
# AUGUSTUS_SA_FDR = FP_AUGUSTUS_SA / (TP_AUGUSTUS_SA + FP_AUGUSTUS_SA),
# AUGUSTUS_SA_FNR = FN_AUGUSTUS_SA / (TP_AUGUSTUS_SA + FN_AUGUSTUS_SA),
# AUGUSTUS_EC_F1 = (2 * TP_AUGUSTUS_EC) / (2 * TP_AUGUSTUS_EC + FP_AUGUSTUS_EC + FN_AUGUSTUS_EC),
# AUGUSTUS_EC_Precision = TP_AUGUSTUS_EC / (TP_AUGUSTUS_EC + FP_AUGUSTUS_EC),
# AUGUSTUS_EC_Sensitivity = TP_AUGUSTUS_EC / (TP_AUGUSTUS_EC + FN_AUGUSTUS_EC),
# AUGUSTUS_EC_FDR = FP_AUGUSTUS_EC / (TP_AUGUSTUS_EC + FP_AUGUSTUS_EC),
# AUGUSTUS_EC_FNR = FN_AUGUSTUS_EC / (TP_AUGUSTUS_EC + FN_AUGUSTUS_EC),
# AUGUSTUS_HS_F1 = (2 * TP_AUGUSTUS_HS) / (2 * TP_AUGUSTUS_HS + FP_AUGUSTUS_HS + FN_AUGUSTUS_HS),
# AUGUSTUS_HS_Precision = TP_AUGUSTUS_HS / (TP_AUGUSTUS_HS + FP_AUGUSTUS_HS),
# AUGUSTUS_HS_Sensitivity = TP_AUGUSTUS_HS / (TP_AUGUSTUS_HS + FN_AUGUSTUS_HS),
# AUGUSTUS_HS_FDR = FP_AUGUSTUS_HS / (TP_AUGUSTUS_HS + FP_AUGUSTUS_HS),
# AUGUSTUS_HS_FNR = FN_AUGUSTUS_HS / (TP_AUGUSTUS_HS + FN_AUGUSTUS_HS)

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

