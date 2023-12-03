if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("hgu133plus2.db")
BiocManager::install("hgu133plus2cdf")
BiocManager::install("affy")
BiocManager::install("annotate")
BiocManager::install("enrichR")
BiocManager::install("rDGIdb")
BiocManager::install("EnsDb.Hsapiens.v86")
# ------------------------------------------------

library(affy)
library(hgu133plus2cdf)
library(hgu133plus2.db) # Opens library with annotation data.
library(help=hgu133plus2.db) # Shows availability and syntax for annotation data.
library(annotate)
library(greenbrown)

# INPUT

input_pval = 0.05
input_upfc = 1.
input_downfc = 1.
page_rank_limit = 100
input_fdr = 0.05

# Read normal sample files in CEL format
rawDataNormal <- ReadAffy(celfile.path = "nsamples/")

# Read tumor sample files in CEL format
rawDataTumor <- ReadAffy(celfile.path = "tsamples/")

# Create normalized and background corrected expression values using the #RMA method. The generated data are stored as ExpressionSet class in the 'eset' #object. For large data sets use the more memory efficient justRMA() function.
normalizedDataNormal <- rma(rawDataNormal)
normalizedDataTumor <- rma(rawDataTumor)

nprobeIDs <- featureNames(normalizedDataNormal)
tprobeIDs <- featureNames(normalizedDataTumor)

nonlyExp <- exprs(normalizedDataNormal)
tonlyExp <- exprs(normalizedDataTumor)

# Perform mapping between Manufacturer Identifiers and Gene Symbols
geneAnnot <- data.frame(SYMBOL=sapply(contents(hgu133plus2ENTREZID), paste, collapse=", "))

# Create a dataframe from the onlyexp data
nonlyExpDF <- data.frame(nonlyExp)
tonlyExpDF <- data.frame(tonlyExp)

# Bind two dataframes
onlyExpDF <- cbind(nonlyExpDF, tonlyExpDF)

# Add the annotations to the dataframe
onlyExpDF <- cbind(onlyExpDF, annot = as.vector(geneAnnot[match(nprobeIDs, rownames(geneAnnot)), 1]))

# Drop the missing columns
onlyExpDF <- subset(onlyExpDF, onlyExpDF["annot"] != "NA")

# Grouping the values which has the same annotation and taking mean
onlyExpDF_agg <- aggregate(onlyExpDF[, 1:(dim(onlyExpDF)[2] - 1)], list(onlyExpDF$annot), mean)

# Differential expression analysis
# Compute fold change: difference between sample means (in logarithmic scale)
if(dim(nonlyExpDF)[2] > 1){
  fcS <- rowMeans(onlyExpDF_agg[,(2 + dim(nonlyExpDF)[2]):(1 + dim(nonlyExpDF)[2] + dim(tonlyExpDF)[2])]) - rowMeans(onlyExpDF_agg[,2:(1 + dim(nonlyExpDF)[2])])
}else{
  fcS <- onlyExpDF_agg[,3] - onlyExpDF_agg[,2]
}

#p.val = apply(onlyExpDF_agg, 1, function(x) { t.test(as.numeric(x[(2 + dim(nonlyExpDF)[2]):(1 + dim(nonlyExpDF)[2] + dim(tonlyExpDF)[2])]), as.numeric(x[2:(1 + dim(nonlyExpDF)[2])]))$p.value } )


#BURAYI SORDUM
p.val = apply(onlyExpDF_agg, 1, function(x) { 
  if(AllEqual(x[2:(1 + dim(nonlyExpDF)[2] + dim(tonlyExpDF)[2])])){return(1)}
  else if(AllEqual(x[2:(1 + dim(nonlyExpDF)[2])]) && AllEqual(x[(2 + dim(nonlyExpDF)[2]):(1 + dim(nonlyExpDF)[2] + dim(tonlyExpDF)[2])])){return(0.0000001)} 
  else {return(t.test(as.numeric(x[(2 + dim(nonlyExpDF)[2]):(1 + dim(nonlyExpDF)[2] + dim(tonlyExpDF)[2])]), as.numeric(x[2:(1 + dim(nonlyExpDF)[2])]))$p.value) }} )


fdr.pval = p.adjust(p.val, method="fdr")
analyzed_data <- cbind(onlyExpDF_agg, fcS, p.val, fdr.pval)
analyzed_data <- analyzed_data[, c("Group.1", "fcS", "p.val", "fdr.pval")]

up_reg_table <- analyzed_data[analyzed_data["fcS"] > input_upfc, ]
down_reg_table <- analyzed_data[analyzed_data["fcS"] < input_downfc, ]

up_reg_table <- up_reg_table[order(-up_reg_table$fcS),]
down_reg_table <- down_reg_table[order(down_reg_table$fcS),]

# GRAPH

edge_list <- read.table("edge_list.txt", sep=" ", header=T)

edge_list <- subset(edge_list, !(is.na(protCol1_sep.Enterezid)) & !(is.na(protCol2_sep.Enterezid)) )

library(igraph)

graph <- graph.data.frame(edge_list)

E(graph)$weigth <- get.data.frame(graph)$ppl_df2.combined_score

vertices <- as.array(V(graph))

personalized_vector <- c()

for (i in 1:length(vertices)){
  personalized_vector <- c( personalized_vector, ifelse(vertices[i] %in% analyzed_data$Group.1, analyzed_data[which(vertices[i] == analyzed_data$Group.1),]$p.val, Inf ))
}

personalized_vector <- 1 / personalized_vector

p_sum <- sum(personalized_vector)

personalized_vector <- personalized_vector / p_sum

page_ranks <- page_rank(graph, directed = F, personalized = personalized_vector)

page_ranks_vector <- as.matrix(page_ranks[["vector"]])

# Sort values
page_ranks_vector <- page_ranks_vector[order(page_ranks_vector[,1],decreasing=TRUE),,drop=FALSE] 

top_genes <- rownames(page_ranks_vector)[1:page_rank_limit]

in_data <- analyzed_data[ which(analyzed_data$Group.1 %in% top_genes), ]

in_data_pos <- subset(in_data, fcS > 0)
in_data_neg <- subset(in_data, fcS < 0)

selected_genes_pos <- in_data_pos[ (abs(in_data_pos$fcS) >= input_upfc & in_data_pos$fdr.pval <= input_fdr), ]
selected_genes_neg <- in_data_neg[ (abs(in_data_neg$fcS) >= input_downfc & in_data_neg$fdr.pval <= input_fdr), ]

library(rDGIdb)
library(EnsDb.Hsapiens.v86)

#Drug analysis for fcS<0 genes


selected_genes_neg$symbol <- mapIds(EnsDb.Hsapiens.v86,
                                       keys = selected_genes_neg$Group.1,
                                       keytype = 'ENTREZID',
                                       column = 'SYMBOL')

neg_gen_df<-selected_genes_neg$symbol

neg_gen_df_result <- queryDGIdb(neg_gen_df,
                                sourceDatabases = c("CGI","COSMIC","ChemblInteractions", "CIViC", "CancerCommons", "ClearityFoundationBiomarkers", "ClearityFoundationClinicalTrial","DTC", "DoCM", "DrugBank","FDA", "JAX-CKB", "MyCancerGenomeClinicalTrial", "OncoKB", "TALC", "TTD", "GuideToPharmacology","MyCancerGenome", "NCI","PharmGKB","TEND","TdgClinicalTrial"),
                                interactionTypes = c("activator","agonist","chaperone","cofactor","inducer","partial agonist", "positive modulator","stimulator","vaccine"))

resultSummary(neg_gen_df_result)


selected_genes_pos$symbol <- mapIds(EnsDb.Hsapiens.v86,
                                       keys = selected_genes_pos$Group.1,
                                       keytype = 'ENTREZID',
                                       column = 'SYMBOL')

pos_gen_df<-selected_genes_pos$symbol

pos_gen_df_result <- queryDGIdb(pos_gen_df,
                                sourceDatabases = c("CGI","COSMIC","ChemblInteractions", "CIViC", "CancerCommons", "ClearityFoundationBiomarkers", "ClearityFoundationClinicalTrial","DTC", "DoCM", "DrugBank","FDA", "JAX-CKB", "MyCancerGenomeClinicalTrial", "OncoKB", "TALC", "TTD", "GuideToPharmacology","MyCancerGenome", "NCI","PharmGKB","TEND","TdgClinicalTrial"),
                                interactionTypes = c("antagonist","antibody","antisense oligonucleotide","blocker","cleavage","inhibitor", "inhibitory allosteric modulator","inverse agonist","negative modulator","partial antagonist", "suppressor"))

resultSummary(pos_gen_df_result)
write.table(pos_gen_df_result, file="test.txt", sep="\t")