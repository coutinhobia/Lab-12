library(limma)
library(edgeR)
library(tidyverse)
library(ggplot2)
library(RColorBrewer) #need colors to make heatmaps
library(gplots) #the heatmap2 function in this package is a primary tool for making heatmaps

# Introduction ----
# The goal of this short script to it make sure you are able to read the count data and the study design for the COVID hackdash into your R environment

# Begin by reading in study design that includes all info for both human and ferret samples
targets <- read_tsv("covid_metadata.txt")
# then read in the human covid data and convert to a matrix with gene symbols as rownames
human_covid_data <- read_tsv("GSE147507_RawReadCounts_Human.tsv")
human_covid_data <- as.matrix(column_to_rownames(human_covid_data, "...1"))


# Now proceed with your exploration and analysis of the data!

# our hypothesis: The transcription profile of A549 cells infected with SARS CoV2 will 
#show differences in immune and cytokine responses compared to the transcription profile of 
#A549-Ace2 cells infected with SARS CoV2.
 
# first we want to select just the groups we want from our dataset.
 
 
test_data<-targets %>% dplyr::filter(sample=="Series5_A549_SARS-CoV-2_1"|
                                       sample=="Series5_A549_SARS-CoV-2_2"|
                                       sample=="Series5_A549_SARS-CoV-2_3"|
                                       sample=="Series16_A549-ACE2_SARS-CoV-2_1"|
                                       sample=="Series16_A549-ACE2_SARS-CoV-2_2"|
                                       sample=="Series16_A549-ACE2_SARS-CoV-2_3")
 
expression_data <- human_covid_data[,test_data$sample]  ## this selects data where the columns match the terms in the test_data character vector
 
#Creating DGE List
 
myDGEList <- DGEList(expression_data)
log2.cpm <- cpm(myDGEList, log=TRUE)
log2.cpm.df <- as_tibble(log2.cpm, rownames = "geneID")
sampleLabels <- targets$sample
 
#Filtering and normalizing the DGE list
 
log2.cpm.df <- as_tibble(log2.cpm, rownames = "geneID")
colnames(log2.cpm.df) <- c("geneID", sampleLabels)
 
 
#This is not strictly needed
log2.cpm.df.pivot <- pivot_longer(log2.cpm.df, # dataframe to be pivoted
                                  cols = !"geneID", # column names to be stored as a SINGLE variable
                                  names_to = "samples", # name of that new variable (column)
                                  values_to = "expression") # name of new variable (column) storing all the values (data)
 
cpm <- cpm(myDGEList)
keepers <- rowSums(cpm>1)>=3 #user defined
myDGEList.filtered <- myDGEList[keepers,]
 
log2.cpm.filtered <- cpm(myDGEList.filtered, log=TRUE)
log2.cpm.filtered.df <- as_tibble(log2.cpm.filtered, rownames = "geneID")
colnames(log2.cpm.filtered.df) <- c("geneID", sampleLabels)
 
myDGEList.filtered.norm <- calcNormFactors(myDGEList.filtered, method = "TMM")
log2.cpm.filtered.norm <- cpm(myDGEList.filtered.norm, log=TRUE)

 
## R did not like the levels in the dataset - called them syntactically invalid. I suspect it doesn't like the dash in "A549-ACE2"
## So I am renaming the levels using the trick Camila taught us. 
 
"%ni%" <- Negate("%in%")
test_data <- test_data %>%
  mutate(cell_line = case_when(
    group %in% "A549-ACE2" ~ "A549Ace2",
    group %in% "A549" ~ "A549",
  ))

# Set up your design matrix ----
group <- factor(test_data$cell_line)
design <- model.matrix(~0 + group)
colnames(design) <- levels(group)


v.DEGList.filtered.norm <- voom(myDGEList.filtered.norm, design, plot = FALSE)
# fit a linear model to your data
fit <- lmFit(v.DEGList.filtered.norm, design)


contrast.matrix <- makeContrasts(celltype = A549Ace2 - A549,
                                 levels=design)
 
#Now trying the contrast matrix again. 
##You MAY have to go back and re-run all of the above steps so that the name change is incorporated into all relevant parts of the script.
contrast.matrix <- makeContrasts(celltype = A549Ace2 - A549,
                                 levels=design)
 
fits <- contrasts.fit(fit, contrast.matrix)
ebFit <- eBayes(fits)
 
#Not strictly necessary, but here is the table of top 50 hits
myTopHits <- topTable(ebFit, adjust ="BH", coef=1, number=50, sort.by="logFC")
myTopHits.df <- myTopHits %>%
  as_tibble(rownames = "geneID")
 
#Now we need to include all of the genes, not just the top 50 for the volcano plot (I just put 40,000 to be bigger than the number of genes)
myTopHits <- topTable(ebFit, adjust ="BH", coef=1, number=40000, sort.by="logFC")
myTopHits.df <- myTopHits %>%
  as_tibble(rownames = "geneID")
 
#Creating a volcano plot with all of the same parameters (reference lines) as Dan had in his graphs
vplot <- ggplot(myTopHits) +
  aes(y=-log10(adj.P.Val), x=logFC,) +
  geom_point(size=2) +
  geom_hline(yintercept = -log10(0.01), linetype="longdash", colour="grey", size=1) +
  geom_vline(xintercept = 1, linetype="longdash", colour="#BE684D", size=1) +
  geom_vline(xintercept = -1, linetype="longdash", colour="#2C467A", size=1) +
  #annotate("rect", xmin = 1, xmax = 12, ymin = -log10(0.01), ymax = 7.5, alpha=.2, fill="#BE684D") +
  #annotate("rect", xmin = -1, xmax = -12, ymin = -log10(0.01), ymax = 7.5, alpha=.2, fill="#2C467A") +
  labs(title="Volcano plot",
       subtitle = "Cell types A549 (without Ace receptor) vs A549Ace2 (with Ace2 receptor)",
       caption="Note: Genes that are more highly expressed in A549-Ace2 cells are found on the right side of the volcano") +
  theme_bw()
vplot


#heatmap creation

#clustering

results <- decideTests(ebFit, method="global", adjust.method="BH", p.value=0.01, lfc=2)
colnames(v.DEGList.filtered.norm$E) <- sampleLabels
diffGenes <- v.DEGList.filtered.norm$E[results[,1] !=0,]
dim(diffGenes)

myheatcolors <- rev(brewer.pal(name="RdBu", n=11))
clustRows <- hclust(as.dist(1-cor(t(diffGenes), method="pearson")), method="complete") #cluster rows by pearson correlation
clustColumns <- hclust(as.dist(1-cor(diffGenes, method="spearman")), method="complete")
module.assign <- cutree(clustRows, k=2)
module.color <- rainbow(length(unique(module.assign)), start=0.1, end=0.9) 
module.color <- module.color[as.vector(module.assign)] 

#actual heatmap creation 

heatmap.2(diffGenes, 
          Rowv=as.dendrogram(clustRows), 
          Colv=as.dendrogram(clustColumns),
          RowSideColors=module.color,
          col=myheatcolors, scale='row', labRow=NA,
          density.info="none", trace="none",  
          cexRow=1, cexCol=1, margins=c(20, 30)) 
dev.off()

names(module.color) <- names(module.assign) 

module.assign.df <- as_tibble(as.list(module.assign))
module.assign.pivot <- pivot_longer(module.assign.df, # dataframe to be pivoted
                                    cols = 1:1702, # column names to be stored as a SINGLE variable
                                    names_to = "geneID", # name of that new variable (column)
                                    values_to = "module") # name of new variable (column) storing all the values (data)

module.assign.pivot <- module.assign.pivot %>%
  mutate(moduleColor = case_when(
    module == 1 ~ "#FF9900",
    module == 2 ~ "#FF0099"))


ggplot(module.assign.pivot) +
  aes(module) +
  geom_bar(aes(fill=moduleColor)) +
  theme_bw()

#Choose cluster 1 of interest because it has more genes
modulePick <- 2 #use 'c()' to grab more than one cluster from the heatmap.  e.g., c(1,2)
#now we pull out the genes from this module using a fancy subsetting operation on a named vector
myModule <- diffGenes[names(module.assign[module.assign %in% modulePick]),] 
hrsub <- hclust(as.dist(1-cor(t(myModule), method="pearson")), method="complete") 

# Create heatmap for chosen sub-cluster.
heatmap.2(myModule, 
          Rowv=as.dendrogram(hrsub), 
          Colv=NA, 
          labRow = NA,
          col=rev(myheatcolors), scale="row", 
          density.info="none", trace="none", 
          RowSideColors=module.color[module.assign%in%modulePick], margins=c(20, 30)) 

# Export modules for downstream analysis ----
#prints out genes in the order you see them in the cluster
moduleSymbols <- tibble(geneID = rev(hrsub$labels[hrsub$order]))
moduleData <- diffGenes[moduleSymbols$geneID,]
moduleData.df <- as_tibble(moduleData, rownames = "geneSymbol")
write_tsv(moduleData.df,"module_upRegulated.tsv")


