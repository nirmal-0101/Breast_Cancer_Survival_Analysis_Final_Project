


# Loading required libraries
library(tidyverse)
library(survival)
library(survminer)
library(ggplot2)
library(pheatmap)
library(TCGAbiolinks)

### 1. LOADING RNA, CNV and Clinical data

# Loading RNA data
gene_expression_data <- read.delim("C:/Users/nirma/Documents/FINAL RESEARCH INTERNSHIP PROJECT/Datasets/HiSeqV2",row.names = 1, check.names = FALSE)

head(gene_expression_data[,1:5])

# Loading CNV data
data_cnv <- read.delim("C:/Users/nirma/Documents/FINAL RESEARCH INTERNSHIP PROJECT/Datasets/Gistic2_CopyNumber_Gistic2_all_thresholded.by_genes",row.names =1, check.names = FALSE)


head(data_cnv[,1:5])

# transposing to make samples as rows as genes as columns
data_cnv <- t(data_cnv)
head(data_cnv[,1:5])

# loading clinical data
data_survival <- read.delim("C:/Users/nirma/Documents/FINAL RESEARCH INTERNSHIP PROJECT/Datasets/patient survival data.txt",row.names = 1)

head(data_survival)

head(gene_expression_data[, 1:5])


# transpose to make samples as rows and genes as columns
gene_expression_data <- t(gene_expression_data)

head(gene_expression_data[, 1:5])

# loading metadata
metadata <- GDCquery_clinic(project = "TCGA-BRCA",type = "clinical")

head(metadata)


### 2. PREPROCESSING DATA
# Checking for missing values
sum(is.na(gene_expression_data))
sum(is.na(data_cnv))
sum(is.na(data_survival))

colSums(is.na(gene_expression_data))
colSums(is.na(data_cnv))
colSums(is.na(data_survival))

# Cleaning the survival data

# keeping the OS and OS.time columns 
data_survival <- data_survival[,c("OS","OS.time")]

# drop samples with NA in OS.time
data_survival <- data_survival[!is.na(data_survival$OS.time), ]

head(data_survival)

# checking the values in OS column of survival data
unique(data_survival$OS)


# Matching samples
common_samples <- Reduce(intersect,list(
  rownames(gene_expression_data),
  rownames(data_cnv),
  rownames(data_survival)
))

# subsetting each dataset

gene_expression_data <- gene_expression_data[common_samples, ]
data_cnv <- data_cnv[common_samples, ]
data_survival <- data_survival[common_samples, ]

head(gene_expression_data[,1:5])
head(data_cnv[,1:5])
head(data_survival)


# Filter top 5000 most variable genes from the expression
variance_genes <- apply(gene_expression_data,2,var)
top_genes <- names(sort(variance_genes,decreasing = TRUE))[1:5000]
gene_expression_data <- gene_expression_data[, top_genes]

head(top_genes)

# scale expression data
gene_expression_data <- scale(gene_expression_data)

head(gene_expression_data[,1:5])



              ###### SINGLE-OMICS DATA ANALYSIS #######



# 1. RNA-ONLY SURVIVAL MODEL

# matching and aligning
gene_expression_common_samples <- intersect(rownames(gene_expression_data),rownames(data_survival))


gene_expression_data_for_survival <- gene_expression_data[gene_expression_common_samples,]    

survival_gene_expression <- data_survival[gene_expression_common_samples, ]  

head(gene_expression_data_for_survival[,1:5])


# Performing PCA for dimensionality reduction
pca_gene_expression <- prcomp(gene_expression_data_for_survival,center = TRUE,scale. = TRUE)

summary(pca_gene_expression)  


top_pca <- as.data.frame(pca_gene_expression$x[,1:7])  

head(top_pca)  


# Adding survival data
top_pca$OS <- survival_gene_expression$OS

top_pca$OS.time <- survival_gene_expression$OS.time  



# Performing Cox Regression and Survival Plot

survival_object_gene_expression <- Surv(time = top_pca$OS.time,
                                        event = top_pca$OS)


# using PC1 to build a cox model
cox_gene_expression <- coxph(survival_object_gene_expression ~ PC1,data = top_pca)

summary(cox_gene_expression)  


# grouping samples
top_pca$Group <- ifelse(top_pca$PC1 > median(top_pca$PC1),"High","Low")

# Creating Kaplan-Meier Plot
km_plot_gene_expression <- survfit(Surv(OS.time,OS) ~ Group, data = top_pca)


# visualising KM plot
ggsurvplot(km_plot_gene_expression,
           data = top_pca,
           pval = TRUE,
           title = "Survival based on PC1 (Gene expression)",
           risk.table = TRUE,
           legend.title = "PC1 group")





###### SURVIVAL ANALYSIS USING CNV DATA ALONE ####

pca_cnv <- prcomp(data_cnv,center = TRUE,scale. = TRUE)

summary(pca_cnv)  

# getting PC1 scores
pc1_cnv <- pca_cnv$x[,1]

head(pc1_cnv)


# Merging PC1 with survival data
survival_data_cnv <- data.frame(
  PC1 = pc1_cnv,
  OS = data_survival$OS,
  OS.time = data_survival$OS.time
)


# creating survival object

survival_object_cnv <- Surv(time = survival_data_cnv$OS.time,event = survival_data_cnv$OS)


# Cox regression on PC1
cox_cnv <- coxph(survival_object_cnv ~ PC1,data = survival_data_cnv)

summary(cox_cnv)


# grouping samples based on PC1
survival_data_cnv$Group <- ifelse(survival_data_cnv$PC1 > median(survival_data_cnv$PC1),"High","Low")

# creating KM plot
km_plot_cnv <- survfit(Surv(OS.time,OS) ~ Group,data = survival_data_cnv)

# visualising plot
ggsurvplot(
  km_plot_cnv,
  data = survival_data_cnv,
  pval = TRUE,
  risk.table = TRUE,
  title = "Survival based on PC1 (CNV)",
  legend.title = "PC1 group"
)




                ##### MULTI-OMICS FACTOR ANALYSIS (MOFA2) ####


library(data.table)
library(MOFA2)


# Preparing the data for MOFA2

# MOFA requires samples to be stored in columns and feature in rows
gene_expression_matrix <- t(as.matrix(gene_expression_data))
cnv_matrix <- t(as.matrix(data_cnv))

head(gene_expression_matrix[,1:5])

head(cnv_matrix[,1:5])


mofa_data <- list(
  "RNA" = gene_expression_matrix,
  "CNV" = cnv_matrix
  
)


# Creating MOFA object
MOFAobject <- create_mofa(mofa_data)



# Visualise the structure of the data
plot_data_overview(MOFAobject)



# Defining training options
data_opts <- get_default_data_options(MOFAobject)
head(data_opts)


# Defining model options
model_opts <- get_default_model_options(MOFAobject)
head(model_opts)



# Defining train options
train_opts <- get_default_training_options(MOFAobject)
head(train_opts)



# Preparing MOFA object
MOFAobject <- prepare_mofa(
  object = MOFAobject,
  data_options = data_opts,
  model_options = model_opts,
  training_options = train_opts
)


# Train the MOFA model
outfile = file.path(getwd(),"mofa_model_new.hdf5")
MOFAobject.trained <- run_mofa(MOFAobject,outfile,use_basilisk = TRUE)






######## DOWNSTREAM ANALYSIS ########

#loading libraries
library(ggplot2)
library(MOFA2)


# Loading trained model
MOFAobject.trained <- load_model("C:/Users/nirma/Documents/FINAL RESEARCH INTERNSHIP PROJECT/mofa_model_new.hdf5")

# getting sample metadata

samples_metadata <- samples_metadata(MOFAobject.trained)


# matching sample IDs
common_samples_metadata <- intersect(samples_metadata$sample,rownames(data_survival))

samples_metadata <- samples_metadata[samples_metadata$sample %in% common_samples_metadata, ]

data_survival <- data_survival[common_samples_metadata, ]


# Adding OS and OS.time
samples_metadata$OS <- data_survival$OS
samples_metadata$OS.time <- data_survival$OS.time


# updating MOFA object with this modified metadata
samples_metadata(MOFAobject.trained) <- samples_metadata


# Overall data overview
plot_data_overview(MOFAobject.trained)


# visualisation of MOFA factor scores by survival
plot_factor(MOFAobject.trained,factors = 1,color_by = "OS")

# Plotting variance explained by each factor
plot_variance_explained(MOFAobject.trained)

plot_variance_explained(MOFAobject.trained,x = "group",y = "factor",plot_total = T)




install.packages("GGally")

# visualisation of combinations of factors
plot_factors(MOFAobject.trained,
             factors = 1:3,
             color_by = "OS")

plot_weights(MOFAobject.trained,
             view = "CNV",
             factor = 2,
             nfeatures = 10,     # Number of features to highlight
             scale = T,          # Scale weights from -1 to 1
             abs = F             # Take the absolute value?
)

# visualisation of feature weights
plot_top_weights(MOFAobject.trained,
                 view = "RNA",
                 factor = 1,
                 nfeatures = 10
)
plot_top_weights(MOFAobject.trained,
                 view = "CNV",
                 factor = 1,
                 nfeatures = 10
)





plot_data_scatter(MOFAobject.trained,
                  view = "RNA",         # view of interest
                  factor = 1,             # factor of interest
                  features = 5,           # number of features to plot (they are selected by weight)
                  add_lm = TRUE,          # add linear regression
                  color_by = "OS"
)

plot_data_scatter(MOFAobject.trained,
                  view = "CNV",         # view of interest
                  factor = 1,             # factor of interest
                  features = 5,           # number of features to plot (they are selected by weight)
                  add_lm = TRUE,          # add linear regression
                  color_by = "OS"
)




###### SURVIVAL ANALYSIS USING COX REGRESSION ########

library(survival)
library(survminer)


# getting MOFA factors
data_factors <- get_factors(MOFAobject.trained,factors = "all",as.data.frame = TRUE)

head(data_factors)
test <- data_factors %>% dplyr::filter(factor == "Factor2")


#reshaping data factors
library(tidyr)
data_factors_wide <- pivot_wider(data_factors,names_from = factor,values_from = value)

head(data_factors_wide)


# merging with survival data

MOFA_metadata <- samples_metadata(MOFAobject.trained)

# merging using sample ID
cox_data <- merge(data_factors_wide,MOFA_metadata,by.x = "sample",by.y = "sample")
head(cox_data[,1:10])


  
  # Running COX regression on each factor
  
  # creating a survival object
  object_survival <- Surv(time = cox_data$OS.time,event = cox_data$OS)
  
  
  factors <- names(cox_data)[grep("^Factor", names(cox_data))]
  
  cox_regression_result <- setNames(
    lapply(factors, function(fac) {
      formula <- as.formula(paste("object_survival ~", fac))
      coxph(formula, data = cox_data)
    }),
    factors
  )
  
  # summarising p-values and hazard ratios
  cox_results_summary <- data.frame(
    Factor = factors,
    HR = sapply(cox_regression_result,function(x) exp(coef(x))),
    pvalue = sapply(cox_regression_result,function(x) summary(x)$coefficients[,"Pr(>|z|)"])
    
  )
  
  print(cox_results_summary)

  
  # filtering significant factors
  significant_factors <- cox_results_summary[cox_results_summary$pvalue < 0.05, ]
  print(significant_factors)
  
  
  # Plotting kaplan-meir curve
  
  cox_data$Group <- ifelse(cox_data$Factor2 > median(cox_data$Factor2),"High","Low")

  
  # Plot
  km_plot <- survfit(Surv(OS.time,OS) ~ Group,data = cox_data)
  
  ggsurvplot(km_plot,
             data = cox_data,
             pval = TRUE,
             title = "Survival Plot of MOFA Factor 2",
             risk.table = TRUE,
             legend.title = "Factor2 Group")
             

  
  
  
  
 
  
  #### TOP PREDICTIVE GENES FOR MOFA FACTOR 2 ####
  
  # extracting weights
  weights_factor2_RNA_MOFA <- get_weights(MOFAobject.trained,view = "RNA",factor = 2,as.data.frame = FALSE)
  str(weights_factor2_RNA_MOFA)
  
  weights_matrix_RNA <- weights_factor2_RNA_MOFA$RNA
  
  weights_df_RNA <- data.frame(
    Gene = rownames(weights_matrix_RNA),
    Weight = as.numeric(weights_matrix_RNA),
    stringsAsFactors = FALSE
  )
  
  # sort by absolute weight
  top_genes_factor2_RNA_MOFA <- weights_df_RNA[order(abs(weights_df_RNA$Weight),decreasing = TRUE),][1:10, ]
  
  # printing top genes
  print(top_genes_factor2_RNA_MOFA)
  
  #saving top genes
  write.csv(top_genes_factor2_RNA_MOFA,"top_factor2_RNA_MOFA.csv",row.names = FALSE)
  
  
  # extracting weights for CNV
  weights_factor2_CNV_MOFA <- get_weights(MOFAobject.trained,view = "CNV",factor = 2,as.data.frame = FALSE)
  str(weights_factor2_CNV_MOFA)
  
  weights_matrix_CNV <- weights_factor2_CNV_MOFA$CNV
  
  weights_df_CNV <- data.frame(
    Gene = rownames(weights_matrix_CNV),
    Weight = as.numeric(weights_matrix_CNV),
    stringsAsFactors = FALSE
  )
  
  # sort by absolute weight
  top_genes_factor2_CNV_MOFA <- weights_df_CNV[order(abs(weights_df_CNV$Weight),decreasing = TRUE),][1:10, ]
  
  # printing top genes
  print(top_genes_factor2_CNV_MOFA)
  
  #saving top genes
  write.csv(top_genes_factor2_CNV_MOFA,"top_factor2_CNV_MOFA.csv",row.names = FALSE)
  
  
  
  
          ########### SURVIVAL ANALYSIS USING MIXOMICS ##################
  
  library(mixOmics)

  
  # Preparing input data 
  data_diablo = list(
    RNA = gene_expression_data,
    CNV = data_cnv)

  # checking the dimensions
  lapply(data_diablo,dim)
  
  
  # Preparing the Y variable 
  
  
  # creating binary group for survival
  Y <- ifelse(data_survival$OS.time > median(data_survival$OS.time),"long survival","short survival")

  Y <- as.factor(Y)

  
  # Initital Analysis
  list.keepX = c(25, 25) 
  list.keepY = c(25, 25)  

  
  # Generating a pairwise PLS model
  pls <- spls(data_diablo[["RNA"]],data_diablo[["CNV"]],
              keepX = list.keepX,keepY = list.keepY)

  
  # plotting features of PLS
  plotVar(pls, cutoff = 0.5, title = "RNA vs CNV", 
          legend = c("RNA", "CNV"), 
          var.names = FALSE, style = 'graphics', 
          pch = c(16, 17), cex = c(2,2), 
          col = c('darkorchid', 'lightgreen'))

  
  # calculate correlation of RNA and CNV
  cor(pls$variates$X, pls$variates$Y)  
  
  
  # Initial DIABLO MODEL
  design = matrix(0.1, ncol = length(data_diablo), nrow = length(data_diablo), 
                  dimnames = list(names(data_diablo), names(data_diablo)))
  diag(design) = 0 # set diagonal to 0s
  
  design
  
  
  # form basic DIABLO model
  basic.diablo.model = block.splsda(X = data_diablo, Y = Y, ncomp = 5, design = design) 
  
  
  # Tuning the number of components
  
  # run component number tuning with repeated CV
  perf.diablo = perf(basic.diablo.model, validation = 'Mfold', 
                     folds = 10, nrepeat = 10) 
  
  # plot output of tuning
  plot(perf.diablo) 

  
  # set the optimal ncomp value
  ncomp = perf.diablo$choice.ncomp$WeightedVote["Overall.BER", "centroids.dist"] 

  
  # show the optimal choice for ncomp for each dist metric
  perf.diablo$choice.ncomp$WeightedVote 

  
  # set grid of values for each component to test
  list.keepX = list(RNA = 25,CNV = 25)
  ncomp = 3

  library(BiocParallel)
  

  # Building the final DIABLO model 
  final.diablo.model <- block.splsda(
    X = data_diablo,
    Y = Y,
    ncomp = ncomp,
    keepX = list.keepX,
    design = design
  )
  # design matrix for the final model
  final.diablo.model$design 
  
  # the features selected to form the first component
  selectVar(final.diablo.model, block = 'RNA', comp = 1)$RNA$name   

  
  # Sample plots
  plotDiablo(final.diablo.model, ncomp = 1)

  plotIndiv(final.diablo.model, ind.names = FALSE, legend = TRUE, 
            title = 'DIABLO Sample Plots')  

  plotArrow(final.diablo.model, ind.names = FALSE, legend = TRUE, 
            title = 'DIABLO')  

  
  
  # Variable plots
  plotVar(final.diablo.model, var.names = FALSE, 
          style = 'graphics', legend = TRUE,
          pch = c(16, 17), cex = c(2,2), 
          col = c('brown1', 'lightgreen'))
  
  
  network(final.diablo.model, blocks = c(1,2),
          color.node = c('brown1', 'lightgreen'), cutoff = 0.4)
  
  
 
  plotLoadings(final.diablo.model, comp = 2, contrib = 'max', method = 'median')
  
  cimDiablo(final.diablo.model)
  
  # saving the model
  saveRDS(final.diablo.model, file = "final_diablo_model.rds")
  
  # Loading trained diablomodel
  final.diablo.model<- readRDS("C:/Users/nirma/Documents/FINAL RESEARCH INTERNSHIP PROJECT/final_diablo_model.rds")
  
  
  
  
  # Extracting diablo latent components
  components_diablo <- data.frame(
    component1_RNA = final.diablo.model$variates$RNA[,1],
    component1_CNV = final.diablo.model$variates$CNV[,1],
    sample = rownames(final.diablo.model$variates$RNA)
  )

  # merging with survival data
  head(data_survival)
  data_survival$sample <- rownames(data_survival)
  diablo_survival_data <- merge(components_diablo,data_survival,by = "sample")

  head(diablo_survival_data)
  
  
  # Cox regression on DIABLO Components
  cox_RNA <- coxph(Surv(OS.time,OS) ~ component1_RNA,data = diablo_survival_data)
  summary(cox_RNA)
  
  cox_CNV <- coxph(Surv(OS.time,OS) ~ component1_CNV,data = diablo_survival_data)

  summary(cox_CNV)  
  
  
  
  # Kaplan-Meier curve for DIABLO
  diablo_survival_data$Group <- ifelse(diablo_survival_data$component1_RNA > median(diablo_survival_data$component1_RNA),"High","Low")

  KM_plot_diablo <- survfit(Surv(OS.time,OS) ~ Group,data = diablo_survival_data)  

  print(KM_plot_diablo)
  
  windows()
  x11()
  # plotting KM curve
  ggsurvplot(KM_plot_diablo,
             data = diablo_survival_data,
             pval = TRUE,
             title = "Survival analysis DIABLO Component1
                      (RNA)",
             risk.table = TRUE,
             legend.title = "Group")
  

  
  # Creating a table of top features with weights
  RNA_weights_diablo <- selectVar(final.diablo.model,block = "RNA",comp = 1)$RNA

  CNV_weights_diablo <- selectVar(final.diablo.model,block = "CNV",comp =1)$CNV
  
  # creating a dataframe with top_features for RNA and weights
  top_features_diablo_RNA <- data.frame(
    Gene = RNA_weights_diablo$name,
    Weight = RNA_weights_diablo$value)

  head(top_features_diablo_RNA)  
  
  # creating a dataframe with top_features for CNV and weights
  
  top_features_diablo_CNV <- data.frame(
    Gene = CNV_weights_diablo$name,
    Weight = CNV_weights_diablo$value)

    
  head(top_features_diablo_CNV) 

  
  # saving top features to csv
  write.csv(top_features_diablo_RNA,"top_RNA_features_DIABLO.csv",row.names = FALSE)

  write.csv(top_features_diablo_CNV,"top_CNV_features_DIABLO.csv",row.names = FALSE)
  
  
  ######### COMPARE RNA TOP GENES and CNV top features from MOFA and DIABLO ##########
  
  
  #loading top genes files
  mofa_rna_top <- read.csv("C:/Users/nirma/Documents/FINAL RESEARCH INTERNSHIP PROJECT/top_factor2_RNA_MOFA.csv")

  diablo_rna_top <- read.csv("C:/Users/nirma/Documents/FINAL RESEARCH INTERNSHIP PROJECT/top_RNA_features_DIABLO.csv")

  head(mofa_rna_top)
  head(diablo_rna_top)
  
  # compare between genes
  common_genes_RNA <- intersect(mofa_rna_top$Gene,diablo_rna_top$Gene)
  genes_mofa_unique_rna <- setdiff(mofa_rna_top$Gene,diablo_rna_top$Gene)
  genes_diablo_unique_rna <- setdiff(diablo_rna_top$Gene,mofa_rna_top$Gene)

  cat("Overalapping Genes found:\n");print(common_genes_RNA)  

  cat("\n Genes unique to MOFA (RNA) only:\n");print(genes_mofa_unique_rna)
  cat("\n Genes unique to DIABLO (RNA) only:\n"); print(genes_diablo_unique_rna)  
  
  
  
  #### PATHWAY ENRICHMENT ANALYSIS ####
  
  genes_mofa_pathwayenrichment <- gsub("_RNA","",top_genes_factor2_RNA_MOFA$Gene)
  genes_diablo_pathwayenrichment <- gsub("_RNA","",top_features_diablo_RNA$Gene)

  head(genes_mofa_pathwayenrichment)
  head(genes_diablo_pathwayenrichment)

  
  genes_mofa_pathwayenrichment    
  genes_diablo_pathwayenrichment
  
  
  # Install required packages
  BiocManager::install("clusterProfiler")
  BiocManager::install("org.Hs.eg.db")
  install.packages("enrichplot")
  
  # Load libraries
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(enrichplot)
  
  
  # Convert MOFA gene symbols to Entrez
  genes_mofa_entrez <- bitr(genes_mofa_pathwayenrichment,
                            fromType = "SYMBOL",
                            toType = "ENTREZID",
                            OrgDb = org.Hs.eg.db)
  
  # Convert DIABLO gene symbols
  genes_diablo_entrez <- bitr(genes_diablo_pathwayenrichment,
                              fromType = "SYMBOL",
                              toType = "ENTREZID",
                              OrgDb = org.Hs.eg.db)

  
  
  # MOFA genes enrichment
  go_mofa <- enrichGO(gene = genes_mofa_entrez$ENTREZID,
                       OrgDb = org.Hs.eg.db,
                       ont = "BP",
                       pAdjustMethod = "BH",
                       pvalueCutoff = 0.05,
                       readable = TRUE)
  
  # DIABLO genes enrichment
  go_diablo <- enrichGO(gene = genes_diablo_entrez$ENTREZID,
                         OrgDb = org.Hs.eg.db,
                         ont = "BP",
                         pAdjustMethod = "BH",
                         pvalueCutoff = 0.05,
                         readable = TRUE)
  
  
  # Dotplots
  dotplot(ego_mofa, showCategory = 10, title = "MOFA RNA Enrichment")
  dotplot(ego_diablo, showCategory = 10, title = "DIABLO RNA Enrichment")
  
  # Barplots
  barplot(ego_mofa, showCategory = 10, title = "MOFA")
  barplot(ego_diablo, showCategory = 10, title = "DIABLO")
  
  
  head(as.data.frame(ego_diablo))
  
  
  # Summarising top enriched pathways
  head(as.data.frame(go_mofa)[,c("Description","p.adjust")],10)
  head(as.data.frame(go_diablo)[,c("Description","p.adjust")],10)
  
  
  
  
  ######## COMPARING C-INDEXES ACROSS THE MODELS ##########
  
  BiocManager::install("survcomp")

  library(survcomp)  

  
  # Calculating C-index for MOFA (factor 2)
  mofa_c_index <-  concordance.index(
    x = cox_data$Factor2,
    surv.time = cox_data$OS.time,
    surv.event = cox_data$OS,
    method = "noether"
  )

  mofa_c_index$c.index
  
  
  # calculating c-index for DIABLO component 1
  diablo_c_index <- concordance.index(
    x = diablo_survival_data$component1_RNA,
    surv.time = diablo_survival_data$OS.time,
    surv.event = diablo_survival_data$OS,
    method = "noether"
  )

  
  diablo_c_index$c.index  
  
  
  # Calculating C-index for RNA PCA
  rna_pca_c_index <- concordance.index(
    x = top_pca$PC1,
    surv.time = top_pca$OS.time,
    surv.event = top_pca$OS,
    method = "noether"
  )

  
  rna_pca_c_index$c.index  

  
  # Calculating C-index for CNV PCA
  
  cnv_pca_c_index <- concordance.index(
    x = survival_data_cnv$PC1,
    surv.time = survival_data_cnv$OS.time,
    surv.event = survival_data_cnv$OS,
    method = "noether"
  )

  cnv_pca_c_index$c.index  

  
  
  # Summarising all the c-index results in a table
  all_c_index <- data.frame(
    Model = c("MOFA (Factor 2)","DIABLO (Component 1 RNA)","Single-Omics RNA","Single-Omics CNV"),
    C_Index = c( mofa_c_index$c.index,  diablo_c_index$c.index, rna_pca_c_index$c.index,cnv_pca_c_index$c.index)
  
  )

  print(all_c_index)  
  
  
  
`             ######## Creating multifactor cox model MOFA ###########
 
  
  
   cox_multifactor <- coxph(Surv(OS.time, OS) ~ Factor2 + Factor11 + Factor13, data = cox_data)
  summary(cox_multifactor)
  

  c_index_multifactor<- summary(cox_multifactor)$concordance[1]

  all_c_index <- rbind(
    all_c_index,
    data.frame(
      Model = "(MOFA Multi-factor)",
      C_Index = c_index_multifactor
    )
  )  

  
  all_c_index  
  
  cox_data$Risk_score <-stats::predict(cox_multifactor)
  cox_data$Risk_Group <- ifelse(cox_data$Risk_score > median(cox_data$Risk_score),"High Risk","Low Risk")
  
  # Creating KM plot
  KM_multifactor <- survfit(Surv(OS.time,OS) ~ Risk_Group, data = cox_data)
  
  # Plotting KM curve
  ggsurvplot(KM_multifactor,
             data = cox_data,
             pval = TRUE,
             title = "Survival Analysis(MOFA Multi-factor)",
             risk.table = TRUE,
             legend.title = "Risk Group")

  
  
  # visualisation of C-index comparison
  library(ggplot2)
  
  ggplot(all_c_index, aes(x = reorder(Model, C_Index), y = C_Index, fill = Model)) +
    geom_bar(stat = "identity", width = 0.6) +
    coord_flip() +
    theme_minimal() +
    labs(title = "C-index Comparison of 
         Survival Models",
         x = "Model",
         y = "C-index") +
    theme(legend.position = "none")
  
  
  

  
  
 

 
  
  
  
  
  