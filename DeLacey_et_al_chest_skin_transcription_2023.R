### DeLacey et al. 2023 
### Transcriptional changes capture sex differences in vascularization in a skin coloration signal in a wild primate

# Sections:
# 1) Chest redness in male and female geladas: Figures 2, S1, and S2
# 2) Sex differences in gene expression: Figures S5, S2, S3, S4, 3
# 3) Enrichment analysis vascularization: Figure 4, S5
# 4) Enrichment analysis androgen and estrogen expression: Figure S6


################################################################################
###    1) Chest redness in male and female geladas: Figure 2 and S1      ###
################################################################################

# prepare workspace

# clear workspace
rm(list = ls())

# load libraries
library(dplyr)
library(tidyr)
library(ggplot2)
library(lme4) 
library(lmerTest) 
library(ggpubr)

# set working data
setwd("[Insert file location here]")

# load data for Results subsection: "Chest redness in male and female geladas"
load("chest_red_nat_anes.RData")
# loads 4 data frames:
# red_nat: Figure 2a, stat 
# red_range_nat: Figure 2b, stat
# red_anes: Figure 2c, stat
# red_heat_change_anes: Figure S1

# prep transparent colors for plotting
makeTransparent <- function(black, alpha = 200){
  newColor <- col2rgb(black)
  apply(newColor, 2, function(curcoldata)			
  {rgb(red = curcoldata[1],
       green = curcoldata[2],
       blue = curcoldata[3],
       alpha = alpha,
       maxColorValue=  255)})
}
tBlack <- makeTransparent("black")
tPurple <- makeTransparent("purple4") # for females
tCyan <- makeTransparent("darkcyan") # for males

############################---------------------###############################
# Under natural conditions (not anesthetized), compare male and female chest redness
dim(red_nat) # 144 photos from 24 adult males and 13 adult females 

model1 <- lmer(rg ~ (1|id) + (1|camera_brand) + sex, data=red_nat)
summary(model1)
# Under natural conditions (i.e., not anesthetized), males and females overlap in redness substantially
# Males display only marginally redder chests than females (Beta=0.11, P=0.06). 

# remove objects to clean up work space
rm(model1)

############################---------------------###############################
# Under natural conditions (not anesthetized), compare male and female chest redness

range(red_nat$rg)

# Plot Figure 2a
fig2a <- ggplot(red_nat, aes(x=sex, y=rg, fill=sex)) +
  geom_boxplot(width=0.5) + 
  ggtitle("Natural Conditions: Redness") + 
  geom_dotplot(binaxis='y', stackdir='center', position=position_dodge(1), dotsize=1.5, binwidth=0.025) +
  scale_fill_manual(values = c(tPurple, tCyan)) +
  scale_y_continuous(name="Redness (Red/Green)",
                     limits= c(1.2,3.2),
                     breaks = seq(1.2,3.2,0.2)) +
  scale_x_discrete(labels=c("Female", "Male"),
                   name="") +
  theme(axis.text =element_text(size=18,
                               family = "Arial"),
        axis.title.y=element_text(size=18,
                                family = "Arial",
                                margin =margin(r=12)), 
        plot.title=element_text(size=18,
                                family="Arial"),
        panel.background = element_blank(),
        axis.line=element_line(color = "black"),
        legend.position = "none") + 
  annotate("text", x=1.5, y=3.2, label = "P=0.06", size=6)
fig2a 

############################---------------------###############################
# Under natural conditions (not anesthetized), compare male and female RANGE in chest redness (max - min for each individual)
dim(red_range_nat) # 23 adult males, 13 adult females

# Plot Figure 2b
fig2b <- ggplot(red_range_nat, aes(x=sex, y=range, fill=sex)) +
  geom_boxplot(width=0.5) +
  ggtitle("Natural Conditions: Range in Redness") +
  geom_dotplot(binaxis='y', stackdir='center', position=position_dodge(1), dotsize=1.0, binwidth=0.025) +
  scale_fill_manual(values = c(tPurple, tCyan)) +
  scale_y_continuous(name="Range in Redness (Max - Min)",
                     limits= c(0,1.4),
                     breaks = seq(0,1.4,0.2)) +
  scale_x_discrete(labels=c("Female", "Male"),
                   name="") +
  theme(axis.text=element_text(size=18,
                               family = "Arial"),
        axis.title.y=element_text(size=18,
                                  family = "Arial",
                                  margin =margin(r=12)), 
        plot.title=element_text(size=18,
                                family="Arial"),
        panel.background = element_blank(),
        axis.line=element_line(color = "black"),
        legend.position = "none") + 
  annotate("text", x=1.5, y=1.4, label = "P=0.007**", size=6)
fig2b

# Stats to accompany Figure 2b
model2 <- lm(range ~ sex*camera_brand, data=red_range_nat)
summary(model2)
# Males have a wider range in redness within individuals compared to females (Beta=0.63, P=0.007)

############################---------------------###############################
# While anesthetized, compare male and female redness
dim(red_anes) # 20 males, 18 females

# Plot Figure 2
fig2c <- ggplot(red_anes, aes(x=sex, y=rg, fill=sex)) +
  geom_boxplot(width=0.5) + 
  geom_dotplot(binaxis='y', stackdir='center', position=position_dodge(1), dotsize=0.75, binwidth =0.02) + 
  scale_fill_manual(values = c(tPurple, tCyan)) +
  scale_y_continuous(name="Redness (Red/Green)",
                     limits = c(1.0, 1.8), 
                     breaks = seq(1.0, 1.8, 0.2)) +
  scale_x_discrete(labels=c("Female", "Male"),
                   name="") + 
  ggtitle("Anesthetized Conditions: Redness") + 
  theme(axis.text=element_text(size=18,
                               family = "Arial"),
        axis.title.y=element_text(size=18,
                                  family = "Arial",
                                  margin =margin(r=12)),
        plot.title=element_text(size=18,
                                family="Arial"),
        panel.background = element_blank(),
        axis.line=element_line(color = "black"),
        legend.position = "none") +
  annotate("text", x=1.5, y=1.8, label = "P=0.40", size=6)
fig2c 

# Stats to accompany figure 2
model3 <-  lm(rg ~ camera_brand + sex, data=red_anes)
summary(model3)
# While under anesthesia, male geladas did not have redder chests than females (Beta=0.04, P=0.40)

# Make a multipaneled plot for Fig 2a-c
fig2 <- ggarrange(fig2a, fig2b, fig2c, 
                  nrow=1, ncol=3,
                  legend="none", 
                  common.legend = FALSE, 
                  labels = c("A", "B", "C"),
                  font.label = list(size=18, 
                                    color="black", 
                                    face="bold", 
                                    family='Arial', 
                                    vjust=1.5))
fig2 


# save plot
#setwd("[File location to save figures")]
# ggsave("fig2.jpeg", plot = fig2, width=18, height=8, units="in", dpi=600)
dev.off()

# remove objects to clean up work space
rm(red_anes, fig2, fig2a, fig2b, fig2c, model1,model2, model3)

############################---------------------###############################
# Supplementary Material: 
# Change in redness between baseline and application of a heat pack directly to the skin
# compared between males and females
# 13 males, 3 females

# Plot Figure S1
figs1 <- ggplot(red_heat_change_anes, aes(x=sex, y=change, fill=sex)) +
  geom_boxplot() +
  ggtitle("Anesthetized Conditions") +
  geom_dotplot(binaxis='y', stackdir='center', position=position_dodge(1), dotsize=1.1, binwidth =0.02) +
  scale_fill_manual(values = c(tPurple, tCyan)) + 
  scale_x_discrete(labels = c("Female", "Male"),
                   name = "") + 
  scale_y_continuous(name="Change in Redness: Heat-Baseline",
                     limits= c(-0.2,0.4),
                     breaks = seq(-0.2,0.4,0.2)) +
  theme(axis.text = element_text(size=24,
                                 family = 'Arial'),
        axis.title=element_text(size=24,
                                family = 'Arial'),
        panel.background = element_blank(),
        axis.line=element_line(color = "black"),
        plot.title=element_text(size=24,
                                family = 'Arial'),
        legend.position = "none")
figs1


# save plot
#setwd("[File location to save figures")]
# ggsave("figs1.jpeg", plot = figs1, width=8, height=6, units="in", dpi=600)
dev.off()

# remove objects to clean up work space
rm(red_heat_change_anes, figs1)

################################################################################
###     2) Sex Differences in Gene Expression: Figures S2, S3, S4, 3         ###
################################################################################

# clear workspace to begin new section
rm(list = ls())

# load libraries 
library(biomaRt)
library(cobs) 
library(doParallel)
library(edgeR)
library(FactoMineR)
library(limma)
library(reshape2)
library(stringr)
library(qvalue)
library(topGO)
library(venn)
library(foreach)
library(mosaic)
library(see)
library(ggpubr)
library(ggrepel)
library(tibble)
library(EMMREML)


# prep transparent colors
makeTransparent <- function(black, alpha = 200){
  newColor <- col2rgb(black)
  apply(newColor, 2, function(curcoldata)			
  {rgb(red = curcoldata[1],
       green = curcoldata[2],
       blue = curcoldata[3],
       alpha = alpha,
       maxColorValue=  255)})
}
tBlack <- makeTransparent("black") # all figures
tGreen <- makeTransparent("darkgreen") # for figure S2
tOrange <- makeTransparent("orange") # for figure S2
tPurple <- makeTransparent("purple4") # for adult females
tLightPurple <- makeTransparent("plum2") # for subadult females
tThistle <- makeTransparent("thistle1") # for FDR threshold of 5% for females
tCyan <- makeTransparent("darkcyan") # for adult males
tPTurq <- makeTransparent("paleturquoise") # for subadult males
tLightCyan <- makeTransparent("lightcyan") # for FDR threshold of 5% for males
tRed <- makeTransparent("darkred") # for figure 4

############################---------------------###############################
# Prepare RNA-Seq data for analysis

# set working data
setwd("/Users/patsydelacey/Documents/1. Research/3. RNA-Seq Project/1. MOLEC ECOL SUBMISSION/3. DRYAD")

## load data
load("PID_10120_mmul10.RData")
# loads R objects for RNA-Seq analyses
# 
# for section 2) Sex Differences in Gene Expression: Figures S2, S3, S4, 3 
  # genes: gene names for ensembl(V1) and VGNC(V2)
  # metadata: technical and biological variables from all samples 
  # meta_teeth: uploads teeth morphometric data needed for aging adults
  # PID_10120_mmul10: gene counts mapped to the annotated macaca mulatta genome
  # hb_rrna_genes: list of hemoglobin and ribosomal RNA genes from the annotated macaca mulatta genome
  # mmul10_gene_chr: used for annotation of chromosomal regions
#
# for section 4: 
  # AR: loads ophid results for Androgen Receptor protein
  # ESR1: loads ophid results for Estrogen Receptor Alpha protein
  # ESR2: loads ophid results for Estrogen Receptor Beta protein

## proportion of reads mapping to rRNA or HB genes
hist(apply(PID_10120_mmul10[hb_rrna_genes,],2,sum)/apply(PID_10120_mmul10,2,sum),25)
## remove rRNA or HB genes mapped reads
counts=PID_10120_mmul10[!rownames(PID_10120_mmul10)%in%hb_rrna_genes,]

## View LID by RQN and rna_conc
ggplot(metadata,aes(x=RQN,y=rna_conc,label=LID))+geom_point()+geom_label_repel()
summary(metadata$RQN)
# check whether LID_103739 & LID_103759 are less than the mean-2*SD
sd(metadata$RQN)
8.755 - (2*1.538461) # both fall below the 5.68 RQN cutoff

## remove samples with low RQNs (RNA qulaity): LID_103739 & LID_103759
new_counts=counts[,!colnames(counts)%in%c("LID_103739","LID_103759")]
new_meta <- dplyr::slice(metadata, 1, 3:21, 23:38)

# add teeth-based age categories
meta_teeth <- meta_teeth %>%
  dplyr::select(LID, age, age_cat_teeth)
# add to new_meta
new_meta <- left_join(new_meta, meta_teeth, by="LID")

# remove objects from work space
rm(counts, metadata, meta_teeth, PID_10120_mmul10)

# Info for Supplementary Material Table S1
tables1 <- new_meta %>%
  dplyr::group_by(Year, age_cat_teeth) %>% 
  summarise(N = n())
rm(tables1)

###------------------------------------###
# Read Count Normalization

# convert new counts to counts per million
cpm=apply(new_counts,2,function(x){x/sum(x)*1e6}) 
dim(cpm) # all = 22,508

## Divide male and female counts
# convert cpm from list to dataframe
cpm <- as.data.frame(cpm)
# Males only cpm thresholds
cpm_male <- cpm %>% dplyr::select(LID_103738, LID_103740,LID_103741,LID_103743,LID_103744,LID_103746,LID_103747,LID_103752,
                                  LID_103755,LID_103758,LID_103762,LID_103763,LID_103764,LID_103765,LID_103766,LID_103769,
                                  LID_103773,LID_103774,LID_103776,LID_103777)
# Males filter out cpm < 10
cpm_10_m <- cpm_male[apply(cpm_male,1,function(x){median(x) > 10}), ]
dim(cpm_10_m) # 10,011
# convert rowname to column (gene)
cpm_10_m <- cpm_10_m %>% 
  rownames_to_column() %>% 
  dplyr::rename(gene = rowname)
# Females only cpm thresholds
cpm_female <- cpm %>% dplyr::select(LID_103742,LID_103745,LID_103748,LID_103749,LID_103750,LID_103751,LID_103753,LID_103754,
                                    LID_103756,LID_103757,LID_103760,LID_103767,LID_103770,LID_103771,LID_103772,LID_103775)
# Females filter out cpm < 10
cpm_10_f <- cpm_female[apply(cpm_female,1,function(x){median(x) > 10}), ]
dim(cpm_10_f) #9,645
# convert rowname to column (gene)
cpm_10_f <- cpm_10_f %>% 
  rownames_to_column() %>% 
  dplyr::rename(gene = rowname)

# Filtered df for female and male to genes only
cpm_10_f_gene <- cpm_10_f %>% dplyr::select(gene)
cpm_10_m_gene <- cpm_10_m %>% dplyr::select(gene)

## Final decision to keep genes that are greater than 10cpm for either females OR males (union) 
# union of gene name only
cpm_10_genes<- dplyr::union(cpm_10_f_gene,cpm_10_m_gene)
dim(cpm_10_genes) # 10,226 genes total 

## get raw counts for the genes that passed the threshold
# make new_counts a data frame
new_counts <- as.data.frame(new_counts)
# change rownames to column
new_counts <- new_counts %>% 
  rownames_to_column() %>% 
  dplyr::rename(gene = rowname)
# get raw counts for only have the genes that passed threshold
filtered_gene_counts <- dplyr::semi_join(new_counts, cpm_10_genes, by = "gene")
# convert gene column to rownames
filtered_gene_counts <- filtered_gene_counts %>% 
  column_to_rownames("gene") 
# check that dimension is 10,226 genes by 36 individuals 
dim(filtered_gene_counts)

# remove data frames that are not needed from workspace 
rm(cpm, cpm_10_f, cpm_10_f_gene, 
   cpm_10_m, cpm_10_m_gene, 
   cpm_female, cpm_male, cpm_10_genes, hb_rrna_genes)

###------------------------------------###
# Correlation matrix of metadata variables

# transpose filtered gene counts
t_counts <- as.data.frame(as.matrix(t(filtered_gene_counts)))
# now LID is in the row names and genes are columns
t_counts$sum <- apply(t_counts, 1, sum)
# creates object of counts for LIDs across all genes
t_counts <- t_counts %>% 
  rownames_to_column() %>% 
  dplyr::rename(LID = rowname)
# change RID_date to a factor then numeric so it can be used for correlation plot
new_meta$RID_date_cor <- as.factor(new_meta$RID_date)
new_meta$RID_date_cor <- as.numeric(new_meta$RID_date_cor)

## join metadata with gene counts and format
t_norm_counts_meta <- left_join(t_counts, new_meta, by = "LID") %>%
  dplyr::select(c(sum, Sex, Year, RID_date_cor, rna_conc, RQN)) %>%
  mutate(Year = as.numeric(Year)) %>% 
  mutate(Sex = ifelse(Sex == "M", 1, 0))
# make correlation plot
library(corrplot)
library(RColorBrewer)
corrplot::corrplot(cor(t_norm_counts_meta), method = "number", type = "upper", col = brewer.pal(n = 8, name = "Blues"), tl.col = "black")
dev.off()

# remove objects from workspace
rm(t_norm_counts_meta, t_counts) 

############################---------------------###############################
# Modeling the effect of sex on gene expression

# Focus on adult individuals for this sexually selected trait
# Adults were determined by eruption of third molars

# create meta data df of adults
meta_adult <- new_meta %>%
  filter(age_cat_teeth != "subadult female") %>% 
  filter(age_cat_teeth != "subadult male")

# select those 28 adult individuals from filtered_gene_counts
LID_adult <-  meta_adult$LID
filtered_gene_counts_adult <- filtered_gene_counts %>% 
  dplyr::select(LID_adult)
dim(filtered_gene_counts_adult) #10,226 x 28 

# create design matrix for the base hypothesis including batch effects with meta_adult
d_meta <- model.matrix(~Sex + Year + RID_date_cor + rna_conc + RQN, data=meta_adult)
# voom matrix with counts and design matrix 
v_matrix <- voom(counts=filtered_gene_counts_adult, design=d_meta) 
# linear model of full model
lm_full <- eBayes(lmFit(v_matrix, d_meta)) # do not need resid matrix since batch effects are included
betas <- as.data.frame(lm_full$coefficients)
pvals <- as.data.frame(lm_full$p.value)
# view p-value distribution for sex
hist(pvals$SexM ,100,  xlab = "pvals of SexM (Male)", main = "Histogram of pvals SexM (Male)")
# strong effect when only adult males and adult females 

###------------------------------------###
# Run EMMREML for adults 

# create Z matrix (identity matrix)
Z_matrix <- diag(1, nrow = 28)
# set the parameters for parallel computing for lmm using EMMREML
ncores <- detectCores(logical = TRUE)  # find number of available cores
clus <- makeCluster(ncores)  # create a cluster object using available cores
registerDoParallel(cores = ncores)
clusterExport(clus,
              varlist = c("meta_adult", "Z_matrix"),
              envir = environment())  # initiate local cluster 
emma_gel_adult2 <- data.frame(t(parApply(clus, v_matrix, 1, function(y){ 
  library(EMMREML)
  emma <- emmreml(y = y,  # model each gene
                  X = model.matrix(~ Sex  + Year + RID_date_cor + rna_conc + RQN, meta_adult),  # design matrix
                  Z = as.matrix(Z_matrix),  # identity matrix
                  K = as.matrix(Z_matrix),  # identity matrix in place of kinship matrix
                  varbetahat = T, varuhat = T, PEVuhat = T, test = T)
  p <- emma$pvalbeta  # store p-values as object
  b <- emma$betahat  # store effect sizes as object
  v <- emma$varbetahat  # store variance in the effect sizes
  l <- emma$loglik # store loglikelihood for the model
  return(c(b[[2]], v[[2]], l, p[2, "none"]))   # return desired effect size and p-value for sex only
})))
colnames(emma_gel_adult2) <- c("beta", "var_beta", "loglik", "pval")  # rename columns


emma_gel_adult2 <- emma_gel_adult2 %>% 
  rownames_to_column() %>% 
  dplyr::rename(gene = rowname)

# calculate standard beta
emma_gel_adult2$std_beta=emma_gel_adult2$beta/sqrt(emma_gel_adult2$var_beta)
emma_gel_adult2 <- emma_gel_adult2 %>% dplyr::mutate(abs_std_beta = abs(std_beta))
# Find qvalue for 20% FDR and add to df
q <- qvalue(p = emma_gel_adult2$pval, fdr.level = 0.20, pfdr = FALSE)
q_emma_gel_adult2 <- as.data.frame(q$qvalues)
names(q_emma_gel_adult2)[1] <- "qvalues_FDR20" 
emma_gel_adult2 <- dplyr::bind_cols(emma_gel_adult2, q_emma_gel_adult2)
rm(q_emma_gel_adult2)

# double check that I calculated std_beta correctly 
qplot(abs(emma_gel_adult2$beta/sqrt(emma_gel_adult2$var_beta)),-log10(emma_gel_adult2$pval))
# perfect correlation between abs(std_beta) and p-values (plotted on a -lo10 scale on the yaxis)

# histogram EMMREML Model
hist(emma_gel_adult2$pval,100, main = "Histogram of p-values for the effect of sex \n in detectably expressed genes", 
     xlab = "p-values", col = tPurple, ylim = c(0,700), cex.main=1.5)
dev.off()

# genes under FDR 5%
emma_gel_adult2_FDR5 <- emma_gel_adult2 %>% dplyr::filter(qvalues_FDR20 < 0.05)
dim(emma_gel_adult2_FDR5) #410 genes
#rm(emma_gel_adult2_FDR5)
# genes under FDR 10%
emma_gel_adult2_FDR10 <- emma_gel_adult2 %>% dplyr::filter(qvalues_FDR20 < 0.1)
dim(emma_gel_adult2_FDR10) #639 genes
#rm(emma_gel_adult2_FDR10)
# how many genes fall under FDR 20% 
emma_gel_adult2_FDR20 <- emma_gel_adult2 %>% dplyr::filter(qvalues_FDR20 < 0.2)
dim(emma_gel_adult2_FDR20) # 1,077 genes pass 

# clean up workspace
rm(betas, clus, d_meta, lm_full, pvals, q, v_matrix, Z_matrix, ncores, LID_adult)

############################---------------------###############################
# Annotate selected model (emma_gel_adult2) for chromosomal locations

# join to emma_gel_adult2
emma_gel_adult2 <- dplyr::left_join(emma_gel_adult2, mmul10_gene_chr, by = "gene")
# join to the 1077 genes that pass a FDR 20% threshold
emma_gel_adult2_FDR20 <- dplyr::left_join(emma_gel_adult2_FDR20, mmul10_gene_chr, by="gene")
# clean workspace
rm(mmul10_gene_chr)

# See how many genes are map to the Y-chromosome
emma_gel_adult2_Y <- emma_gel_adult2 %>%
  filter(chr_scaff_name == "Y")
dim(emma_gel_adult2_Y ) # 14 genes map to Y- chromosome

# See how many genes map to the X-chromosome
emma_gel_adult2_X <- emma_gel_adult2 %>%
  filter(chr_scaff_name == "X")
dim(emma_gel_adult2_X) # 327 genes map to X chromosome

# remove Y-linked genes from subsequent analyses
emma_gel_adult2 <- emma_gel_adult2 %>%
  dplyr::filter(chr_scaff_name != "Y")
dim(emma_gel_adult2) #10,212 detectably expressed genes

# check how many Y-linked genes were in our subset of passing a FDR 20%
emma_gel_adult2_FDR20 <- emma_gel_adult2_FDR20 %>% 
  filter(chr_scaff_name != "Y")
dim(emma_gel_adult2_FDR20) # 1,068 genes left that pass after Y-linked genes are removed

length(emma_gel_adult2_FDR20$gene)/length(emma_gel_adult2$gene)*100
# 10.5% of the 10,212 detectably expressed genes exhibited significant 
# differential expression across males and females (n genes=1,068, FDR<20%)

############################---------------------###############################

# Volcano plot for Figure S5 

emma_gel_adult2_volcano <- emma_gel_adult2 %>% 
  mutate(diff_exp = case_when(
    std_beta > 3.0 ~ "FDR5_males", 
    std_beta < -3.0 ~ "FDR5_females",
    std_beta < 3.0 & std_beta > 2.65 ~ "FDR10_males", 
    std_beta > -3.0 & std_beta < -2.65 ~ "FDR10_females", 
    std_beta < 2.65 & std_beta > 2.208 ~ "FDR20_males", 
    std_beta > -2.65 & std_beta < -2.208 ~ "FDR20_females", 
    std_beta > -2.208 & std_beta < 2.208 ~ "not")) 

# order in a way I prefer
emma_gel_adult2_volcano$diff_exp <- factor(emma_gel_adult2_volcano$diff_exp, 
                                           levels = c("FDR5_females", "FDR10_females", "FDR20_females", 
                                                      "FDR5_males", "FDR10_males", "FDR20_males", "not"))

emma_gel_adult2_volcano <- emma_gel_adult2_volcano %>% 
  mutate(log_10 = -log10(pval))

# PLOT FIGURE S5
figs5 <- ggplot(data=emma_gel_adult2_volcano, aes(x=beta, y=-log10(pval), col=diff_exp)) + 
  geom_point(alpha = 0.4) + 
  scale_color_manual(values=c("purple4", "mediumpurple1", "plum", "darkslategray", "cyan4", "paleturquoise", "black"), 
                     labels=c("FDR < 5% females", "FDR < 10% females", "FDR < 20% females", "FDR < 5% males", "FDR < 10% males", "FDR < 20% males", "Not differentially expressed")) + 
  scale_x_continuous(name="logFC") + 
  theme_classic(base_size=20) + 
  geom_hline(yintercept=2.56, linetype = "dashed", color = "gray50", size = 0.75) + 
  geom_hline(yintercept=2.1, linetype = "dashed", color = "gray50", size = 0.75) + 
  geom_hline(yintercept=1.565, linetype = "dashed", color = "gray50", size = 0.75) + 
  theme(panel.background = element_blank(),
        axis.text=element_text(size=18,
                               family = 'Arial'),
        axis.title=element_text(size=18,
                                family = 'Arial'),
        axis.line=element_line(color = "black"), 
        legend.position = "none")
figs5

# save plot
#setwd("[File location to save figures")]
# ggsave("figs5.jpeg", plot = figs5, width=6, height=6, units="in", dpi=600)
dev.off()

rm(figs5, emma_gel_adult2_volcano)


############################---------------------###############################
# PCA Plots (Fig S2, S3, and S4)

# Y-linked removed but X-linked included

# remove Y-linked genes from filtered gene counts
filtered_gene_counts_noY <- filtered_gene_counts %>%
  tibble::rownames_to_column("gene")
# remove Y linked genes from filtered_gene_counts
filtered_gene_counts_noY <- filtered_gene_counts_noY %>%
  filter(gene %in% emma_gel_adult2$gene) 
dim(filtered_gene_counts_noY )

rownames(filtered_gene_counts_noY) <- filtered_gene_counts_noY[,1]
filtered_gene_counts_noY <- filtered_gene_counts_noY[,-1]

# voom normalize the raw gene counts
v_exp_noY=voom(filtered_gene_counts_noY)
pca_noY=prcomp(cor(v_exp_noY$E))
summary(pca_noY) # gives proportion of variance for each pc

# make a df of pca data to prepare for plotting
pca_df_noY <- as.data.frame(pca_noY$x)
# make LID a column
pca_df_noY <- pca_df_noY %>%
  tibble::rownames_to_column("LID")
# add year, sex, and age_cat_teeth from new_meta
new_meta_pca <- new_meta %>% 
  dplyr::select(LID, Year, Sex, age_cat_teeth)
# connect info by LID
pca_df_noY <- left_join(pca_df_noY, new_meta_pca, by="LID")
# make year categorical
pca_df_noY$Year <- as.character(pca_df_noY$Year)
# clean workspace
rm(new_meta_pca)

# PLOT FIGURE S2: PCA 1 and PCA 2 by YEAR 
figs2 <- ggplot(pca_df_noY, aes(x=PC1, y=PC2, col=Year)) +
  geom_point(size=5) + 
  scale_color_manual(values = c(tOrange, tGreen),
                     labels=c("2017", "2019")) +
  ggtitle("") +
  scale_x_continuous(name= "PC1 (78.2%)",
                     limits=c(-0.5,2), breaks=seq(-0.5,2,0.5)) +
  scale_y_continuous(name="PC2 (15.3%)",
                     limits =c(-0.5,0.4), breaks=seq(-0.5,0.3,0.2)) + 
  theme(axis.text=element_text(size=18,
                               family = 'Arial'),
        axis.title=element_text(size=18,
                                family = 'Arial'),
        panel.background = element_blank(),
        axis.line=element_line(color = "black"),
        plot.title=element_text(size=24,
                                family = 'Arial'),
        legend.title = element_blank(),
        legend.text=element_text(size=16,
                                 family = 'Arial'),
        legend.position = c(0.8, 0.3), 
        legend.box.background = element_rect(color="black"))
figs2

# save plot
#setwd("[File location to save figures")]
# ggsave("figs2.jpeg", plot = figs2, width=8, height=6, units="in", dpi=600)
dev.off()
rm(figs2)

# Stats to accompany figure S2 (body of paper in results section)
# Influence of year on PC1 and PC2
pc1_lm_year <- lm(PC1 ~ Year, data=pca_df_noY)
summary(pc1_lm_year) # Year is associated with PC1
# PC1: year Beta=0.39, P=0.02
pc2_lm_year <- lm(PC2 ~ Year, data=pca_df_noY)
summary(pc2_lm_year) # Year is also associated with PC2
# PC2: year Beta=-0.22, P=0.003

# scree plot to visualize how much each pc explains the variance
library(devtools)
library(factoextra)
get_eig(pca_noY)
fviz_eig(pca_noY)
dev.off()
# PC1 explains 78.2% of the variance
# PC2 explains 15.5% of the variance

###------------------------------------###
library(bbmle) # for AIC comparisons

# STATS FOR PC1 - SEX and YEAR
# lmm for PC1
pc_model1 <- lmer(PC1 ~ (1|Year) + Sex, pca_df_noY)
summary(pc_model1) # sex is significant along PC1 when controlling for Year 

# lm for PC1
pc_model2 <- lm(PC1 ~ Year + Sex, pca_df_noY)
summary(pc_model2) # Sex is significant Beta=-0.39,  P=0.01
# model is stronger is we include year in the model

# lm for PC1 with interaction 
pc_model3 <- lm(PC1 ~ Year*Sex, pca_df_noY)
summary(pc_model3) # interaction effect does not make model stronger 


AICtab(pc_model1, pc_model2, pc_model3, weights=T) #lm_mf is stronger (including year and sex)
rm(pc_model1, pc_model2, pc_model3)

# Figure out percent variance explained for each variable using our selected LMs
pc1_year <- lm(PC1 ~ Year, pca_df_noY)
summary(pc1_year) # 12.21% for only Year
pc1_lmsexyear <- lm(PC1 ~ Sex + Year, pca_df_noY)
summary(pc1_lmsexyear) # 25.35% of the variance is explained by year and sex (which means 13.14% is explained by sex)

###------------------------------------###
# Partial residual plot for PC1

# Residuals of LMER of effect of Sex while controlling for year (using model pc1_year above)
pca_df_noY$resids <- residuals(pc1_year)

# reorder variables in a way we think is meaningful for the data developmentally
pca_df_noY$age_cat_teeth <- factor(pca_df_noY$age_cat_teeth, levels = c("adult female", "subadult female","subadult male","adult male"))
summary(pca_df_noY$resids)
library(scales) # used to round on the y-axis


# PLOT FIGURE S3: PARTIAL RESIDUAL PLOT FOR PC1
figs3 <- ggplot(pca_df_noY, aes(x=age_cat_teeth, y=resids, fill=age_cat_teeth)) +
  geom_boxplot() + 
  labs(y= "Residuals of PC1 ~ Year", x= "Age Category") + 
  geom_dotplot(binaxis='y', stackdir='center', position=position_dodge(1), dotsize=1, binwidth =0.02) + 
  scale_fill_manual(values = c(tPurple, tLightPurple, tPTurq, tCyan), 
                    labels=c("Adult female", "Subadult female","Subadult male","Adult male")) + 
  scale_x_discrete(labels=c("Adult female", "Subadult female","Subadult male", "Adult male"), 
                   name = "") + 
  scale_y_continuous(limits = c(-0.6, 1.8), breaks = seq(-0.6,1.8,0.2), labels = scales::number_format(accuracy = 0.1)) +
  geom_hline(yintercept = 0, linetype=2, col=tBlack) + 
  ggtitle("") +
  theme(legend.position = "none",
        axis.text=element_text(size=18,
                               family = "Arial"),
        axis.title=element_text(size=24,
                                family = "Arial"),
        panel.background = element_blank(),
        axis.line=element_line(color = "black"),
        plot.title=element_text(size=24,
                                family = "Arial",
        ))
figs3

# save plot
#setwd("[File location to save figures")]
# ggsave("figs3.jpeg", plot = figs3, width=9, height=6, units="in", dpi=600)
dev.off()
rm(figs3)

###------------------------------------###
# Test whether results of the sex difference in gene expression are influenced by X-linked genes

# Model PCA again with both X and Y linked genes removed

# remove X-linked genes from filtered gene counts
dim(emma_gel_adult2_X) #327 genes to remove
filtered_gene_counts_noXY <- filtered_gene_counts_noY %>%
  tibble::rownames_to_column("gene") %>%
  filter(!gene %in% emma_gel_adult2_X$gene)
dim(filtered_gene_counts_noXY) #9,885 genes remain

rownames(filtered_gene_counts_noXY) <- filtered_gene_counts_noXY[,1]
filtered_gene_counts_noXY <- filtered_gene_counts_noXY[,-1]

# voom normalize the raw gene counts
v_exp_noXY=voom(filtered_gene_counts_noXY)
pca_noXY=prcomp(cor(v_exp_noXY$E))
summary(pca_noXY)

# create df for plotting
pca_df_noXY <- as.data.frame(pca_noXY$x)
# make LID a column
pca_df_noXY <- pca_df_noXY %>%
  tibble::rownames_to_column("LID")
# add year, sex, and age_cat_teeth from new_meta
new_meta_pca <- new_meta %>% 
  dplyr::select(LID, Year, Sex, age_cat_teeth)
# connect info by LID
pca_df_noXY <- left_join(pca_df_noXY, new_meta_pca, by="LID")
# clean workspace
rm(new_meta_pca)

# reorder variables
pca_df_noY$Sex <- factor(pca_df_noXY$Sex, levels = c("F","M"))
pca_df_noXY$Sex <- factor(pca_df_noXY$Sex, levels = c("F","M"))

# PLOT FIGURE S4A: Y-Linked Removed for Sex by PC1 and PC2
figs4a <- ggplot(pca_df_noY, aes(x=PC1, y=PC2, col=Sex)) +
  geom_point(size=5) + 
  scale_color_manual(values = c(tPurple, tCyan),
                     labels=c("Female", "Male")) +
  ggtitle("PCA Results: Y-Linked Removed") +
  scale_x_continuous(name= "PC1 (78.2%)",
                     limits=c(-0.5,2), breaks=seq(-0.5,2,0.5)) +
  scale_y_continuous(name="PC2 (15.3%)",
                     limits =c(-0.5,0.4), breaks=seq(-0.5,0.3,0.2)) + 
  theme(axis.text=element_text(size=24,
                               family = 'Arial'),
        axis.title=element_text(size=24,
                                family = 'Arial'),
        panel.background = element_blank(),
        axis.line=element_line(color = "black"),
        plot.title=element_text(size=24,
                                family = 'Arial'),
        legend.title = element_blank(),
        legend.text=element_text(size=20,
                                 family = 'Arial'),
        legend.position = c(0.8, 0.3), 
        legend.box.background = element_rect(color="black"))
figs4a  

# PLOT FIGURE S4B: Autosome Only for Sex by PC1 and PC2
# call plot
figs4b <- ggplot(pca_df_noXY, aes(x=PC1, y=PC2, col=Sex)) +
  geom_point(size=5) + 
  scale_color_manual(values = c(tPurple, tCyan),
                     labels=c("Female","Male")) +
  ggtitle("PCA Results: Autosomes Only") +
  scale_x_continuous(name= "PC1 (78.2%)",
                     limits=c(-0.5,2), breaks=seq(-0.5,2,0.5)) +
  scale_y_continuous(name="PC2 (15.3%)",
                     limits =c(-0.5,0.4), breaks=seq(-0.5,0.3,0.2)) + 
  theme(axis.text=element_text(size=24,
                               family = 'Arial'),
        axis.title=element_text(size=24,
                                family = 'Arial'),
        panel.background = element_blank(),
        axis.line=element_line(color = "black"),
        plot.title=element_text(size=24,
                                family = 'Arial'),
        legend.title = element_blank(),
        legend.text=element_text(size=20,
                                 family = 'Arial'),
        legend.position = c(0.8, 0.3), 
        legend.box.background = element_rect(color="black"))
figs4b

# Make multipanel plot for the supplementary material figure of PCAs with and without X chromosome inclusion
figs4 <- ggarrange(figs4a, figs4b, 
                       nrow =1, ncol =2,
                       common.legend = F,
                       labels = c("A", "B"),
                       font.label = list(size = 24,
                                         color = "black",
                                         face = "bold",
                                         family = 'Arial',
                                         vjust=1.5))
figs4

# save plot
#setwd("[File location to save figures")]
# ggsave("figs4.jpeg", plot = figs4, width=16, height=8, units="in", dpi=600)
dev.off()
rm(figs4a, figs4b, figs4)

# clean workspace
rm(pc1_lm_year, pc1_lmsexyear, pc1_year, pca_df_noXY, pca_df_noY, v_exp_noXY, 
   v_exp_noY, pca_noXY, pca_noY, pc2_lm_year)


############################---------------------###############################
# Average standardized sex-bias gene expression level for each individual

# Use significantly differentially expressed genes from our model
dim(emma_gel_adult2_FDR20) # 1068 left after removing y-linked genes

###---AVERAGE MALE EXPRESSION
# select genes from emmremml that pass FDR 20% threshold that are male bias 
male_bias_adult2 <- emma_gel_adult2_FDR20 %>% dplyr::filter(std_beta > 0)
dim(male_bias_adult2) # 867 genes 

# voom normalize the raw gene counts
v_exp=voom(filtered_gene_counts_noY) 
v_exp_df <- as.data.frame(v_exp) 
# format count data
v_exp_df <- v_exp_df %>% 
  rownames_to_column() %>% 
  dplyr::rename(gene = rowname) 
exp_male_genes <- v_exp_df %>% dplyr::filter(gene %in% male_bias_adult2$gene)
# convert gene column to rowname
exp_male_genes <- exp_male_genes  %>%
  column_to_rownames("gene") 
dim(exp_male_genes)
# transpose 
exp_male_genes <- t(exp_male_genes)
dim(exp_male_genes)
# 36 x 867

# apply scale to rows 
scale_male_bias <- apply(exp_male_genes,2,scale)
# average across columns 
avg_male_exp <- apply(scale_male_bias,1,mean)

# create data frame for avg_male_exp
avg_male_exp <- as.data.frame(as.matrix(avg_male_exp))
dim(avg_male_exp) 

avg_male_exp <- bind_cols(avg_male_exp, new_meta)
avg_male_exp <- avg_male_exp %>% dplyr::rename(avg_exp = V1)

###---AVERAGE FEMALE EXPRESSION
# select genes from emma_gel_heavy that pass FDR 20% threshold that are female bias 
female_bias_genes <- emma_gel_adult2_FDR20 %>% dplyr::filter(std_beta < 0) 
dim(female_bias_genes) # 201 genes

exp_female_genes <- v_exp_df %>% dplyr::filter(gene %in% female_bias_genes$gene)
# convert gene column to rowname
exp_female_genes <- exp_female_genes  %>%
  column_to_rownames("gene") 
dim(exp_female_genes)
# transpose 
exp_female_genes <- t(exp_female_genes)
dim(exp_female_genes)

# apply scale to rows 
scale_female_bias <- apply(exp_female_genes,2,scale)
# average across columns 
avg_female_exp <- apply(scale_female_bias,1,mean)

# create data frame for avg_female_exp
avg_female_exp <- as.data.frame(as.matrix(avg_female_exp))
dim(avg_female_exp)

avg_female_exp <- bind_cols(avg_female_exp, new_meta)
avg_female_exp <- avg_female_exp %>% dplyr::rename(avg_exp = V1)

###---COMBINED MALE-LIKE AND FEMALE-LIKE EXPRESSION
# make combined figure where positive numbers = males; negative numbers = females

# combine scale_male_bias and scale_female_bias genes then average
# flip the sign of female bias genes
scale_female_bias_flip <- scale_female_bias*-1
# combine male and female genes
scale_male_bias_df <- as.data.frame(scale_male_bias)
scale_female_bias_flip_df <- as.data.frame(scale_female_bias_flip)

scale_both_sex <- bind_cols(scale_male_bias_df, scale_female_bias_flip_df)
dim(scale_both_sex) # 36 x 1,068
# convert back to matrix
scale_both_sex <- as.matrix(scale_both_sex)

# average across columns 
avg_combined_exp <- apply(scale_both_sex,1,mean)

# create data frame for avg_combined_exp
avg_combined_exp <- as.data.frame(as.matrix(avg_combined_exp))
dim(avg_combined_exp)

avg_combined_exp <- bind_cols(avg_combined_exp, new_meta)
avg_combined_exp <- avg_combined_exp %>% dplyr::rename(avg_exp = V1)

# prep graph - put in the order we think is showing the developmental trajectory
avg_combined_exp$age_cat_teeth <- factor(avg_combined_exp$age_cat_teeth, levels=c("adult female","subadult female", "subadult male", "adult male"))

# PLOT FIGURE 3: Boxplot for mean expression level for each age_cat_teeth category
fig3 <- ggplot(avg_combined_exp, aes(x=age_cat_teeth, y=avg_exp, fill=age_cat_teeth)) +
  geom_boxplot()
fig3 <- fig3  + labs(y= "Mean Expression Level", x= "Age Category") + 
  geom_dotplot(binaxis='y', stackdir='center', position=position_dodge(1), dotsize=1.5, binwidth =0.035) + 
  scale_fill_manual(values = c(tPurple, tLightPurple, tPTurq, tCyan), 
                    labels=c("Adult female", "Subadult female","Subadult male","Adult male")) + 
  scale_x_discrete(labels=c("Adult female", "Subadult female","Subadult male", "Adult male"),
                   name="") + 
  scale_y_continuous(limits = c(-2, 1.2), breaks = seq(-2,1,0.5)) +
  geom_hline(yintercept = 0, linetype=2, col=tBlack) + 
  ggtitle("") +
  theme(legend.position = "none",
        axis.text=element_text(size=24,
                               family = "Arial"),
        axis.title.y=element_text(size=24,
                                family = "Arial",
                                margin=margin(r=12)),
        panel.background = element_blank(),
        axis.line=element_line(color = "black"),
        plot.title=element_text(size=24,
                                family = "Arial",
        ))
fig3

# ADD SIGNIFICANCE BARS
# prep for adding significance
stat_test_fig3 <- compare_means(avg_exp ~ age_cat_teeth, data = avg_combined_exp)

fig3<- fig3 +
  geom_bracket(aes(xmin="subadult male",
                   xmax="adult male",
                   label="P=0.003**"),
               data=stat_test_fig3,
               y.position=0.95,
               inherit.aes = FALSE,
               label.size=6) +
  geom_bracket(aes(xmin="adult female", 
                   xmax="subadult male",
                   label="P=0.04*"),
               data=stat_test_fig3,
               y.position=0.8,
               inherit.aes = FALSE, 
               label.size=6) +
  geom_bracket(aes(xmin="adult female", 
                   xmax="adult male",
                   label="P=7.95e-8***"),
                   data=stat_test_fig3,
                   y.position=1.2,
                   inherit.aes = FALSE,
               label.size=6) +
  geom_bracket(aes(xmin="adult female", 
                   xmax="subadult female",
                   label="P=0.02*"),
               data=stat_test_fig3,
               y.position=0.55,
               inherit.aes = FALSE,
               label.size=6)
fig3

# save plot
#setwd("[File location to save figures")]
# ggsave("fig3.jpeg", plot = fig3, width=12, height=8, units="in", dpi=600)
dev.off()
rm(fig3)

# Run stats to accompany fig3 
# Adult males compared to all 
avg_combined_exp$age_cat_teeth <- factor(avg_combined_exp$age_cat_teeth, levels = c("adult male", "adult female", "subadult male", "subadult female"))
lm_avgexp_m <- lm(avg_exp ~ age_cat_teeth, avg_combined_exp)
summary(lm_avgexp_m)
# adult males and subadult males differ (Beta=-0.66, P=0.003)
# Adult males and subadult females do not differ (small subadult female sample size)
# adult males and adult females differ (Beta=-1.11, P=7.95x10^-8***)

# Subadult males compared to all 
avg_combined_exp$age_cat_teeth <- factor(avg_combined_exp$age_cat_teeth, levels=c("subadult male", "subadult female", "adult female", "adult male"))
lm_avgexp_sa <- lm(avg_exp ~ age_cat_teeth, avg_combined_exp)
summary(lm_avgexp_sa)
# subadult males differed from adult males (Beta=-0.67, P=0.003) and adult females (Beta=0.45, P=0.04)

# Adult females compared to all
# run stats of mean, rescaled gene expression
avg_combined_exp$age_cat_teeth <- factor(avg_combined_exp$age_cat_teeth, levels = c("adult female", "adult male", "subadult male", "subadult female"))
lm_avgexp_f <- lm(avg_exp ~ age_cat_teeth, avg_combined_exp)
summary(lm_avgexp_f)
# adult females and subadult males differ (Beta=0.45, P=0.04)
# Adult females and males differ (Beta=1.11, P<0.001)
# Adult females and subadult females differ (Beta=0.81,P=0.02)

# clean workspace
rm(avg_combined_exp, avg_female_exp, avg_male_exp, exp_female_genes, exp_male_genes)
rm(female_bias_genes, lm_avgexp_f, lm_avgexp_m, lm_avgexp_sa, male_bias_adult2, 
   scale_both_sex, scale_female_bias, scale_female_bias_flip, scale_female_bias_flip_df)
rm(scale_male_bias, scale_male_bias_df, stat_test_fig3)

rm(emma_gel_adult2_FDR20, emma_gel_adult2_X, emma_gel_adult2_Y, 
   filtered_gene_counts, filtered_gene_counts_adult,
   filtered_gene_counts_noXY,
   genes, meta_adult, new_counts)


################################################################################
###        3) Enrichment analysis vascularization: Figure 4, S5              ###
################################################################################
# Gene Ontology Analysis

# PREPARE FOR GO ANALYSIS
# Get gene information ----------------
# Create a biomart object with ensembl mmul genes
mmul <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                         dataset = paste0("mmulatta", "_gene_ensembl"),
                         host = "https://www.ensembl.org") #MART object from BioMart; query biodatabases


# Create a vector of gene names from expression matrix
ensembl_gene_names <- emma_gel_adult2$gene
# Create a biomart object of GO terms associated with genes from study (takes a while to run)
mmul_GO <- biomaRt::getBM(attributes = c('ensembl_gene_id', 
                                         'go_id',
                                         'name_1006'), 
                          filters = 'ensembl_gene_id',
                          values = ensembl_gene_names, 
                          mart = mmul,
                          uniqueRows = TRUE) 
# Create a list of 
gene_ID2GO <- lapply(unique(mmul_GO$ensembl_gene_id),
                     function(x){sort(mmul_GO[mmul_GO$ensembl_gene_id == x,
                                              'go_id'])})
names(gene_ID2GO) <- unique(mmul_GO$ensembl_gene_id)
#save
# setwd("[select file location to save"])
save(mmul, mmul_GO, ensembl_gene_names, file = "chestDEsex_tm3_GO_emmageladult2.RData")

# can skip above and load this in subsequent runs 
load("chestDEsex_tm3_GO_emmageladult2.RData")

############################---------------------###############################
# Enrichment Analysis: Angiogenesis related genes

# Connect mmul_GO categories to all genes in emma_gel2
mmul_GO <- mmul_GO %>% 
  dplyr::rename(gene = ensembl_gene_id)
# join mmul_GO to all genes from emma_gel_adult2
emma_gel_adult2_GO <- left_join(emma_gel_adult2, mmul_GO, by = "gene") 

# Search for angiogenesis
angiogenesis <- grep("angiogenesis", emma_gel_adult2_GO$name_1006, perl=T, value=F)
angiogenesis_GO <- emma_gel_adult2_GO %>% dplyr::slice(angiogenesis)

# remove terms unrelated to angiogenesis (deal with lymph system)
# remove "lymphangiogenesis"
angiogenesis_GO <- angiogenesis_GO %>% dplyr::filter(!grepl("lymphangiogenesis", name_1006))

angiogenesis_GO_unique <- angiogenesis_GO  %>% dplyr::select(1:8)
angiogenesis_GO_unique  <- distinct(angiogenesis_GO_unique)
dim(angiogenesis_GO_unique) # 195 genes total
angiogenesis_GO_unique$Category <- "Angiogenesis gene"

NOT_angiogenesis_GO <- emma_gel_adult2 %>%
  dplyr::filter(!gene %in% angiogenesis_GO_unique$gene)
dim(NOT_angiogenesis_GO) # 10,017 genes 
NOT_angiogenesis_GO <- NOT_angiogenesis_GO %>% dplyr::select(1:8)
NOT_angiogenesis_GO$Category <- "All other genes"

# bind rows for graph 
density_df_angio <- dplyr::bind_rows(angiogenesis_GO_unique,NOT_angiogenesis_GO)

# find median for each distribution
median(NOT_angiogenesis_GO$std_beta)
median(angiogenesis_GO_unique$std_beta)
# find range
summary(density_df_angio$std_beta) 

density_df_angio$Category <- factor(density_df_angio$Category , levels = c("Angiogenesis gene", "All other genes"))

# Stats to accompany Figure 4
ks.test(NOT_angiogenesis_GO$std_beta,angiogenesis_GO_unique$std_beta, alternative="greater")
# Genes related to angiogenesis are more highly expressed in males 
# K-S Test: D=0.19, P=7.26x10^7

# prep for plotting
library(ggbeeswarm)
# Put df in format for plotting
angio_GO_violindf <- angiogenesis_GO %>%
  dplyr::select(std_beta, go_id, name_1006)
# reorder by std_beta
angio_GO_violindf$name_1006 = with(angio_GO_violindf, reorder(name_1006, std_beta, median))

# PLOT FIGURE 4
fig4 <- ggplot(angio_GO_violindf, aes(x=name_1006, y=std_beta)) +
  geom_beeswarm(cex=1, col=tBlack) +
  coord_flip() + 
  ylab("Standardized effect size of sex") + 
  xlab("Gene Ontology Category") + 
  scale_y_continuous(limits=c(-6,6), breaks=seq(-6,6,2)) +
  geom_hline(yintercept=0, linetype="dotted", color=tBlack,size=1) +
  stat_summary(fun = "mean", 
               geom = "crossbar",
               width =0.5,
               col = tBlack) +
  scale_x_discrete(labels = c("Angiogenesis involved in coronary \nvascular morphogenesis",
                              "Regulation of sprouting angiogenesis",
                              "Regulation of cell migration involved \nin sprouting angiogenesis",
                              "Blood vessel endothelial cell proliferation \ninvolved in sprouting angiogenesis",
                              "Regulation of angiogenesis",
                              "Sprouting angiogenesis",
                              "Positive regulation of sprouting angiogenesis",
                              "Intussusceptive angiogenesis",
                              "Positive regulation of angiogenesis",
                              "Negative regulation of angiogenesis",
                              "Negative regulation of sprouting angiogenesis",
                              "Angiogenesis involved in \nwound healing",
                              "Angiogenesis",
                              "Cell migration involved in \nsprouting angiogenesis",
                              "Positive regulation of cell migration \ninvolved in sprouting angiogenesis",
                              "Negative regulation of blood vessel endothelial cell \nproliferation involved in sprouting angiogenesis",
                              "Negative regulation of cell migration \ninvolved in sprouting angiogenesis",
                              "Positive regulation of blood vessel endothelial cell \nproliferation involved in sprouting angiogenesis")) + 
  geom_hline(yintercept=0.4952145, linetype="solid", color=tRed,size=1) + # add vertical line for all other genes median std beta
  theme(axis.text.y = element_text(size=16,family = "Arial"),
        axis.text.x = element_text(size=24, family="Arial"),
        axis.title=element_text(size=24,family = "Arial"),
        panel.background = element_blank(),
        axis.line=element_line(color = "black"),
        legend.title = element_blank(),
        legend.position = "none") 
fig4
# save plot
#setwd("[File location to save figures")]
# ggsave("fig4.jpeg", plot = fig4, width=16, height=12, units="in", dpi=600)
dev.off()
rm(fig4)

# clean workspace
rm(angio_GO_violindf, angiogenesis, angiogenesis_GO, angiogenesis_GO_unique, NOT_angiogenesis_GO,density_df_angio)

############################---------------------###############################
# Enrichment Analysis: Blood flow/bood pressure related genes

# Search for blood terms
blood <- grep("blood", emma_gel_adult2_GO$name_1006, perl=T, value=F)
blood_GO <- emma_gel_adult2_GO %>% dplyr::slice(blood)
# remove blood-brain terms
blood_GO <- blood_GO %>% dplyr::filter(!grepl("blood-brain", name_1006))
# remove blood-nerve
blood_GO <- blood_GO %>% dplyr::filter(!grepl("blood-nerve", name_1006))
# remove coagulation
blood_GO <- blood_GO %>% dplyr::filter(!grepl("coagulation", name_1006))
# remove placenta
blood_GO <- blood_GO %>% dplyr::filter(!grepl("placenta", name_1006))
# remove phosphatidylserine exposure on blood platelet
blood_GO <- blood_GO %>% dplyr::filter(name_1006 != "phosphatidylserine exposure on blood platelet")
# remove "establishment of blood-retinal barrier"
blood_GO <- blood_GO %>% dplyr::filter(!grepl("establishment of blood-retinal barrier", name_1006))
# remove "blood microparticle"
blood_GO <- blood_GO %>% dplyr::filter(!grepl("blood microparticle", name_1006))
# remove "blood microparticle"
blood_GO <- blood_GO %>% dplyr::filter(!grepl("positive regulation of blood microparticle formation", name_1006))
# sift through 
dim(blood_GO)

# Genes in the blood pressure/blood vessel maintenance category
blood_GO_unique <- blood_GO %>% dplyr::select(1:8)
blood_GO_unique <- distinct(blood_GO_unique)
dim(blood_GO_unique) # 148 genes total
blood_GO_unique$Category <- "Blood pressure/blood vessel maintenace gene"

# All other genes not in the blood pressure/blood vessel maintenance category
NOT_blood_GO <- emma_gel_adult2 %>%
  dplyr::filter(!gene %in% blood_GO_unique$gene)
dim(NOT_blood_GO ) # 10,064 genes 
NOT_blood_GO <- NOT_blood_GO %>% dplyr::select(1:8)
NOT_blood_GO$Category <- "All other genes"

# bind rows for graph 
density_df_blood <- dplyr::bind_rows(blood_GO_unique,NOT_blood_GO)

# Stats to accompany Figure S5
ks.test(NOT_blood_GO$std_beta,blood_GO_unique$std_beta, alternative="greater")
# Genes related to angiogenesis are more highly expressed in males 
# K-S Test: D=0.19, P=2.59x10^5

# find median for each distribution
median(NOT_blood_GO$std_beta)
median(blood_GO_unique$std_beta)
# find range
summary(density_df_blood$std_beta) 

# create df for fig s6
blood_GO_violindf <- blood_GO %>%
  dplyr::select(std_beta, go_id, name_1006)

# see how many categories there are
blood_GO_violindf_unique <- blood_GO_violindf %>%
  dplyr::select(name_1006)
blood_GO_violindf_unique  <- distinct(blood_GO_violindf_unique)
dim(blood_GO_violindf_unique) # 43 total
rm(blood_GO_violindf_unique)

# find which categories are the most represented
blood_GO_violindf_sum <- blood_GO_violindf %>%
  dplyr::group_by(name_1006) %>%
  summarise(N= n())
# select those with N=5 or greater
blood_GO_violindf_sum <- blood_GO_violindf_sum %>%
  filter(N>4) # left with 13 GO categories

# filter violindf to only include these 
blood_GO_violindf <- blood_GO_violindf %>%
  filter(name_1006 %in% blood_GO_violindf_sum$name_1006)

# reorder variable by std_beta
blood_GO_violindf$name_1006 = with(blood_GO_violindf, reorder(name_1006, std_beta, median))

median(NOT_blood_GO$std_beta)

# PLOT FIGURE S6
figs6 <- ggplot(blood_GO_violindf, aes(x=name_1006, y=std_beta)) +
  geom_beeswarm(cex=1, col=tBlack) +
  geom_jitter(position = position_jitter(seed =1, width = 0.2)) +
  coord_flip() +
  ylab("Standardized effect size of sex") + 
  xlab("Gene Ontology Category") + 
  scale_y_continuous(limits=c(-6,6), breaks=seq(-6,6,2)) +
  geom_hline(yintercept=0, linetype="dotted", color=tBlack,size=1) +
  stat_summary(fun = "mean", 
               geom = "crossbar",
               width =0.5,
               col = tBlack) +
  geom_hline(yintercept=0.4983315, linetype="solid", color=tRed,size=1) + # add vertical line for all other genes median std beta
  theme(axis.text.y = element_text(size=16,family = "Arial"),
        axis.text.x = element_text(size=24, family="Arial"),
        axis.title=element_text(size=24,family = "Arial"),
        panel.background = element_blank(),
        axis.line=element_line(color = "black"),
        legend.title = element_blank(),
        legend.position = "none") 
figs6

# save plot
#setwd("[File location to save figures")]
# ggsave("figs6.jpeg", plot = figs6, width=18, height=10, units="in", dpi=600)
dev.off()
rm(figs6)

# clean workspace
rm(blood_GO, blood_GO_unique, blood_GO_violindf, blood_GO_violindf_sum,
   density_df_blood, NOT_blood_GO)

################################################################################
###    4) Enrichment analysis androgen and estrogen expression: Figure S6    ###
################################################################################
# Enrichment Analysis for genes involved in protein-protein interaction networks for ESR1, ESR2, and AR

# Use Interlocus Interaction Database 
# name results
emremml_results = emma_gel_adult2
# order by ascending std_beta
emremml_results = emremml_results[order(emremml_results$std_beta),]
# rename gene to match other dfs
emremml_results <- emremml_results %>% dplyr::rename("ensembl_gene_id" = "gene")

colnames(emremml_results)

# name ophid proteins and their descriptions
ophid.proteins = c('P10275','P03372','Q92731')
ophid.names = c('AR','ESR1','ESR2')
ophid.descriptions =c('360 androgen receptor PPIs (OPHID)', '758 estrogen receptor alpha PPIs (OPHID)', '506 estrogen receptor beta PPIs (OPHID)')
names(ophid.names) = names(ophid.descriptions) = ophid.proteins


hsap = biomaRt::useEnsembl(biomart = 'ENSEMBL_MART_ENSEMBL',dataset='hsapiens_gene_ensembl', mirror="useast")
mmul = biomaRt::useEnsembl(biomart = 'ENSEMBL_MART_ENSEMBL',dataset='mmulatta_gene_ensembl', mirror="useast")


mmul.orthology = getBM(
  attributes = c('ensembl_gene_id','external_gene_name','hsapiens_homolog_ensembl_gene','hsapiens_homolog_orthology_type','hsapiens_homolog_associated_gene_name'),
  mart = mmul
)

hsap.uniprot = getBM(
  attributes = c('ensembl_gene_id','uniprotswissprot'),
  mart = hsap
)

colnames(emremml_results)
colnames(mmul.orthology)
# find one to one orthologs 
emremml.orthologs = merge(
  subset(mmul.orthology,hsapiens_homolog_orthology_type %in% 'ortholog_one2one'),
  emremml_results,
  by='ensembl_gene_id',all.x=FALSE,all.y=FALSE)
dim(emremml.orthologs) # 8390

### PPI FOR ESR1 ENRICHMENT ANALYSIS ###
# filter hsap.uniprot to only include those included in ESR1$UniProt2
hsap_ESR1 <- hsap.uniprot %>%
  filter(uniprotswissprot %in% ESR1$UniProt2)
dim(hsap_ESR1) #820 genes

emremml.orthologs_ESR1 <- emremml.orthologs %>% 
  filter(hsapiens_homolog_ensembl_gene %in% hsap_ESR1$ensembl_gene_id)
dim(emremml.orthologs_ESR1) #553 genes
emremml.orthologs_ESR1$Category <- "ESR1 PPI gene"

emremml.orthologs_notESR1 <- emremml.orthologs %>%
  filter(!ensembl_gene_id %in% emremml.orthologs_ESR1$ensembl_gene_id)
dim(emremml.orthologs_notESR1) # 7837 genes
emremml.orthologs_notESR1$Category <- "All other genes"

# Stats accompanying ESR1 PPI Enrichment 
ks.test(emremml.orthologs_notESR1$std_beta,emremml.orthologs_ESR1$std_beta, alternative="greater")
# Males did not have increased expression related to PPI for ESR1: D=0.05, P=0.07 


### PPI FOR ESR2 ENRICHMENT ANALYSIS ###
# filter hsap.uniprot to only include those included in ESR2$UniProt2
hsap_ESR2 <- hsap.uniprot %>%
  filter(uniprotswissprot %in% ESR2$UniProt2)
dim(hsap_ESR2) # 548 genes

emremml.orthologs_ESR2 <- emremml.orthologs %>% 
  filter(hsapiens_homolog_ensembl_gene %in% hsap_ESR2$ensembl_gene_id)
dim(emremml.orthologs_ESR2) # 345 genes
emremml.orthologs_ESR2$Category <- "ESR2 PPI gene"

emremml.orthologs_notESR2 <- emremml.orthologs %>% 
  filter(!ensembl_gene_id %in% emremml.orthologs_ESR2$ensembl_gene_id)
dim(emremml.orthologs_notESR2) # 8045 genes
emremml.orthologs_notESR2$Category <- "All other genes"

# Stats accompanying ESR2 PPI Enrichment 
ks.test(emremml.orthologs_notESR2$std_beta,emremml.orthologs_ESR2$std_beta, alternative="greater")
#  Males did not have increased expression related to PPI for ESR1: D=0.05, P=0.17


### PPI FOR AR ENRICHMENT ANALYSIS ###
# filter hsap.uniprot to only include those included in ESR2$UniProt2
hsap_AR <- hsap.uniprot %>%
  filter(uniprotswissprot %in% AR$UniProt2)
dim(hsap_AR) # 390 genes

emremml.orthologs_AR <- emremml.orthologs %>% 
  filter(hsapiens_homolog_ensembl_gene %in% hsap_AR$ensembl_gene_id)
dim(emremml.orthologs_AR) # 240 genes
emremml.orthologs_AR$Category <- "AR PPI gene"

emremml.orthologs_notAR <- emremml.orthologs %>% 
  filter(! ensembl_gene_id %in% emremml.orthologs_AR$ensembl_gene_id)
dim(emremml.orthologs_notAR) # 8150
emremml.orthologs_notAR$Category <- "All other genes"

#  Stats accompanying AR PPI Enrichment 
ks.test(emremml.orthologs_notAR$std_beta,emremml.orthologs_AR$std_beta, alternative="greater")
# Males did not have increased expression related to PPI for AR: D=0.04, P=0.46


############################---------------------###############################
# Average standardized ESR1, ESR2, and AR gene expression level for each individual 
# correlated with anesthetized redness values
# Included in supplement (Figure S7)

# load red_anes_rna - list of individual's baseline anesthetized redness (rg) value
# Note: we only have both redness data and RNA-Seq data for N=18 individuals
load("red_anes_rna.RData")

############################
#         ESR1             #
############################
# Genes associated with ESR1
dim(emremml.orthologs_ESR1)

# use voom normalized gene counts 
exp_esr1_genes <- v_exp_df %>% dplyr::filter(gene %in% emremml.orthologs_ESR1$ensembl_gene_id)
# convert gene column to rowname
exp_esr1_genes <- exp_esr1_genes %>% 
  column_to_rownames("gene")
dim(exp_esr1_genes) 
# transpose 
exp_esr1_genes <- t(exp_esr1_genes)
dim(exp_esr1_genes) 

# apply scale to rows 
scale_esr1 <- apply(exp_esr1_genes,2,scale)
# average across columns
avg_esr1_exp <- apply(scale_esr1,1,mean)

# create data frame for avg_esr1_exp
avg_esr1_exp <- as.data.frame(as.matrix(avg_esr1_exp))
dim(avg_esr1_exp) 
# format df and combine with new_meta
avg_esr1_exp <- bind_cols(avg_esr1_exp, new_meta)
avg_esr1_exp <- avg_esr1_exp %>% dplyr::rename(avg_exp = V1)
# add in redness data
avg_esr1_exp_rg <- left_join(avg_esr1_exp, red_anes_rna, by = "LID")
avg_esr1_exp_rg <- avg_esr1_exp_rg %>% 
  filter(!is.na(rg))

summary(avg_esr1_exp_rg$rg)
summary(avg_esr1_exp_rg$avg_exp)

# plot average expression of ESR1 genes
esr1_exp_rg <- ggplot(avg_esr1_exp_rg, aes(x=rg, y=avg_exp, col=Sex)) +
  geom_point(size=5) + 
  ggtitle("ER-Alpha") + 
  geom_smooth(method=lm) + 
  scale_color_manual(values=c(tPurple, tCyan),
                     labels=c("Female", "Male")) +
  scale_x_continuous(name= "Redness (Red/Green)",
                     limits = c(1.1,1.8), breaks=seq(1.1,1.8,0.1)) +
  scale_y_continuous(name="Mean Expression Level", 
                     limits = c(-2,1), breaks=seq(-2,1,0.5)) +
  theme(axis.text=element_text(size=18,
                               family = 'Arial'),
        axis.title=element_text(size=18,
                                family = 'Arial'), 
        panel.background = element_blank(),
        axis.line=element_line(color = "black"),
        plot.title=element_text(size=30,
                                family = 'Arial'),
        legend.title = element_blank(),
        legend.text=element_text(size=18,
                                 family = 'Arial'),
        legend.position = "top", 
        legend.box.background = element_rect(color="black")) 
esr1_exp_rg

# stats within females
avg_esr1_exp_rg_f <- avg_esr1_exp_rg[avg_esr1_exp_rg$Sex == "F", ]
lm_esr1_f <- lm(avg_exp ~ rg, data=avg_esr1_exp_rg_f)
summary(lm_esr1_f)
# Beta=3.797, p=0.271

# stats within males
avg_esr1_exp_rg_m <- avg_esr1_exp_rg[avg_esr1_exp_rg$Sex == "M", ]
lm_esr1_m <- lm(avg_exp ~ rg, data=avg_esr1_exp_rg_m)
summary(lm_esr1_m)
# Beta=0.47, P=0.23

# Stats of expression by sex (does not include redness data)
lm_esr1_exp_sex <- lm(avg_exp ~ Sex, data = avg_esr1_exp_rg)
summary(lm_esr1_exp_sex)
# Males have higher average expression of ER-Alpha genes compared to females 
# Beta = 0.59, P=0.04

############################
#         ESR2             #
############################
# Genes associated with ESR2
dim(emremml.orthologs_ESR2) 

# use voom normmalized gene counts 
exp_esr2_genes <- v_exp_df %>% dplyr::filter(gene %in% emremml.orthologs_ESR2$ensembl_gene_id)
# convert gene column to rowname
exp_esr2_genes <- exp_esr2_genes %>% 
  column_to_rownames("gene")
dim(exp_esr2_genes) 
# transpose 
exp_esr2_genes <- t(exp_esr2_genes)
dim(exp_esr2_genes) 

# apply scale to rows 
scale_esr2 <- apply(exp_esr2_genes,2,scale)
# average across columns
avg_esr2_exp <- apply(scale_esr2,1,mean)

# create data frame for avg_esr1_exp
avg_esr2_exp <- as.data.frame(as.matrix(avg_esr2_exp))
dim(avg_esr2_exp) 
# format df and combine with new_meta
avg_esr2_exp <- bind_cols(avg_esr2_exp, new_meta)
avg_esr2_exp <- avg_esr2_exp %>% dplyr::rename(avg_exp = V1)
# add in redness data
avg_esr2_exp_rg <- left_join(avg_esr2_exp, red_anes_rna, by = "LID")
avg_esr2_exp_rg <- avg_esr2_exp_rg %>% 
  filter(!is.na(rg))

summary(avg_esr2_exp_rg$rg)
summary(avg_esr2_exp_rg$avg_exp)

# plot average expression of ESR2 genes
esr2_exp_rg <- ggplot(avg_esr2_exp_rg, aes(x=rg, y=avg_exp, col=Sex)) +
  geom_point(size=5) + 
  ggtitle("ER-Beta") + 
  geom_smooth(method=lm) + 
  scale_color_manual(values=c(tPurple, tCyan),
                     labels=c("Female", "Male")) +
  scale_x_continuous(name= "Redness (Red/Green)",
                     limits = c(1.1,1.8), breaks=seq(1.1,1.8,0.1)) +
  scale_y_continuous(name="Mean Expression Level", 
                     limits = c(-2,1), breaks=seq(-2,1,0.5)) +
  theme(axis.text=element_text(size=18,
                               family = 'Arial'),
        axis.title=element_text(size=18,
                                family = 'Arial'),
        panel.background = element_blank(),
        axis.line=element_line(color = "black"),
        plot.title=element_text(size=30,
                                family = 'Arial'),
        legend.title = element_blank(),
        legend.text=element_text(size=18,
                                 family = 'Arial'),
        legend.position = "top", 
        legend.box.background = element_rect(color="black")) 
esr2_exp_rg

# stats within females
avg_esr2_exp_rg_f <- avg_esr2_exp_rg[avg_esr2_exp_rg$Sex == "F", ]
lm_esr2_f <- lm(avg_exp ~ rg, data=avg_esr2_exp_rg_f)
summary(lm_esr2_f)
# Beta=3.47, p=0.288

# stats within males
avg_esr2_exp_rg_m <- avg_esr2_exp_rg[avg_esr2_exp_rg$Sex == "M", ]
lm_esr2_m <- lm(avg_exp ~ rg, data=avg_esr2_exp_rg_m)
summary(lm_esr2_m)
# Beta=0.46, P=0.25

# Stats of expression by sex (does not include redness data)
lm_esr2_exp_sex <- lm(avg_exp ~ Sex, data = avg_esr2_exp_rg)
summary(lm_esr2_exp_sex)
# Males have higher average expression of ER-Beta genes compared to females 
# Beta = 0.58, P=0.04

############################
#           AR             #
############################
# Genes associated with AR
dim(emremml.orthologs_AR) 

# use voom normmalized gene counts 
exp_ar_genes <- v_exp_df %>% dplyr::filter(gene %in% emremml.orthologs_AR$ensembl_gene_id)
# convert gene column to rowname
exp_ar_genes <- exp_ar_genes %>% 
  column_to_rownames("gene")
dim(exp_ar_genes) 
# transpose 
exp_ar_genes <- t(exp_ar_genes)
dim(exp_ar_genes) 

# apply scale to rows 
scale_ar <- apply(exp_ar_genes,2,scale)
# average across columns
avg_ar_exp <- apply(scale_ar,1,mean)

# create data frame for avg_esr1_exp
avg_ar_exp <- as.data.frame(as.matrix(avg_ar_exp))
dim(avg_ar_exp) 
# format df and combine with new_meta
avg_ar_exp <- bind_cols(avg_ar_exp, new_meta)
avg_ar_exp <- avg_ar_exp %>% dplyr::rename(avg_exp = V1)
# add in redness data
avg_ar_exp_rg <- left_join(avg_ar_exp, red_anes_rna, by = "LID")
avg_ar_exp_rg <- avg_ar_exp_rg %>% 
  filter(!is.na(rg))

summary(avg_ar_exp_rg$avg_exp)

# plot average expression of AR genes
ar_exp_rg <- ggplot(avg_ar_exp_rg, aes(x=rg, y=avg_exp, col=Sex)) +
  geom_point(size=5) + 
  ggtitle("AR") + 
  geom_smooth(method=lm) + 
  scale_color_manual(values=c(tPurple, tCyan),
                     labels=c("Female", "Male")) +
  scale_x_continuous(name= "Redness (Red/Green)",
                     limits = c(1.1,1.8), breaks=seq(1.1,1.8,0.1)) +
  scale_y_continuous(name="Mean Expression Level", 
                     limits = c(-2,1), breaks=seq(-2,1,0.5)) +
  theme(axis.text=element_text(size=18,
                               family = 'Arial'),
        axis.title=element_text(size=18,
                                family = 'Arial'),
        panel.background = element_blank(),
        axis.line=element_line(color = "black"),
        plot.title=element_text(size=30,
                                family = 'Arial'),
        legend.title = element_blank(),
        legend.text=element_text(size=18,
                                 family = 'Arial'),
        legend.position = "top", 
        legend.box.background = element_rect(color="black")) 
ar_exp_rg

figs7 <- ggarrange(esr1_exp_rg, esr2_exp_rg, ar_exp_rg,
          nrow=1, ncol=3, 
          labels = c("A", "B", "C"),
          common.legend = FALSE,
          font.label = list(size= 24,
                            color = "black",
                            face = "bold",
                            family = 'Arial',
                            vjust=1.5))
figs7 

# setwd
# setwd("[file location for saving figures")
# ggsave("figs7.jpeg", plot = figs7, width=24, height=10, units="in", dpi=600)


# stats within females
avg_ar_exp_rg_f <- avg_ar_exp_rg[avg_ar_exp_rg$Sex == "F", ]
lm_ar_f <- lm(avg_exp ~ rg, data=avg_ar_exp_rg_f)
summary(lm_ar_f)
# Beta=2.91, p=0.374

# stats within males
avg_ar_exp_rg_m <- avg_ar_exp_rg[avg_ar_exp_rg$Sex == "M", ]
lm_ar_m <- lm(avg_exp ~ rg, data=avg_ar_exp_rg_m)
summary(lm_ar_m)
# Beta=0.47, P=0.325

# Stats of expression by sex (does not include redness data)
lm_ar_exp_sex <- lm(avg_exp ~ Sex, data = avg_ar_exp_rg)
summary(lm_ar_exp_sex)
# Males gave marginally higer average expression of AR genes compared to females
# Beta = 0.53, P=0.053

# clear workspace
rm(list = ls())

################################################################################
###                                END                                       ###
################################################################################
