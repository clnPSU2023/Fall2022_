### Microbiome Analysis of Maize Leaf and Seed Endophytes ###
### Theo Newbold - 2022 to 2023 - Pennsylvania State University 

#### Phyloseq Pipeline ####

#Download and Install Necessary Packages 
install.packages("rlang") #Phyloseq requires specific version of rlang; need to download most current version to properly use phyloseq downstream 
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c("phyloseq","Biostrings","RDPutils","ggplot2","dplyr", "ggpubr", "vegan", "emmeans", "plyr"))

#Load Packages
library(phyloseq)
library(Biostrings)
library(RDPutils) # not available for R version of, 4.3.2
library(ggplot2)
library(dplyr)
library(ggpubr)
library(vegan)
library(emmeans)
library(plyr)
install.packages("tidyr")
library(tidyr)
BiocManager::install("microbiome")
library(microbiome)

#Import Metadata or Mapping File
Metadata <- read.csv("~/Desktop/Phyllo_Micro/MrDNA_Pipeline/MetaData_NegConRemove_23Jun23.csv", header = TRUE, row.names = 1, stringsAsFactors = F, na.strings = c(".", "-")) 

#Verify Metadata or Mapping file structure.
str(Metadata)
head(Metadata, 10)

# Changing stage, management style, plot type
# and plot to factors makes them easier to work with in the future
Metadata$Stage <- as.factor(Metadata$Stage)
Metadata$Man_Style <- as.factor(Metadata$Man_Style)
Metadata$Plot_Type <- as.factor(Metadata$Plot_Type)
Metadata$Plot <- as.factor(Metadata$Plot)

# Bring in your taxa and OTU tables

# Read in taxa table
# the end part (row.names = 1) is very important. It removes column one from the actual data 
taxaITS <- read.csv("~/Desktop/Phyllo_Micro/MrDNA_Pipeline/Taxa_Table.csv", row.names = 1)
str(taxaITS)

# As it is, the taxa table is a data frame but it needs to be a matrix in order to merge
taxaITS <- as.matrix(taxaITS) # makes it a matrix
class(taxaITS) # check class
# its a matrix, hooray! 

# Read in OTU Table
# Again, the end part (row.names = 1) is very important. It removes column one from the actual data 
# This is important here because the first column in not numeric and the rest of the table needs to be
OTUITS <-read.csv("~/Desktop/Phyllo_Micro/MrDNA_Pipeline/OTU.csv",row.names=1 )

#you will see that it reads the counts in as integers, you must change this to numeric
str(OTUITS)
#find out the number of columns 
ncol(OTUITS)

#create a subset of the data OTU[1:2676] then use lapply to convert it to 
#numeric and then assign it back to those columns
OTUITS[1:85] <- lapply(OTUITS[1:85], as.numeric)
str(OTUITS)

# transpose OTU table
otus.t <- t(OTUITS)

#now it should work to create a merged phyloseq object 
psITS <- phyloseq(otu_table(otus.t, taxa_are_rows=FALSE), tax_table(taxaITS))

#can create a phyloseq object of each individual thing (otu and taxa)
#can then merge these 2 phyloseq objects into one 
otuITS <- otu_table(otus.t, taxa_are_rows = FALSE)
taxITS <- tax_table(taxaITS)

psITS <- merge_phyloseq(otuITS, taxITS)

#once that works you can do it with your metadata
meta <- sample_data(Metadata)
psITS.bac <- merge_phyloseq(psITS, meta)

#Create a Phyloseq object. Here we merged the otu table, taxonomy table, and Metadata/Mapping File.
psITS.bac <- merge_phyloseq(otu_table(otus.t,taxa_are_rows = FALSE), tax_table(taxaITS),sample_data(Metadata))

#Duplicate the phyloseq object to preserve original, this is important! 

psITS.bac1 <- psITS.bac #working phyloseq object 

#Order samples by the number of reads or counts or library size or depth of coverage
ordered(sample_sums(psITS.bac1))

#Save RDS and Read RDS
psITS.bac1
saveRDS(psITS.bac1, "MR_DNA_Fungi_Data.rds")

#Read .RDS By using <-, we are restoring the RDS under the same name we used for our Phyloseq object.
psITS.bac <- readRDS("~/Desktop/Phyllo_Micro/2023_05_04_041323SCits_Raw_Data_UDI/demux/MR_DNA_Fungi_Data.rds") # error appeared but looks to have worked anyway???
psITS.bac

### Now, let's start Filtering, Subsetting, and Performing Richness, Ordination, and Multivariate analysis ###

#Verify the prevalence of each taxon in other words, the number of samples in which a taxon appears at least once.
# Remember you didn't do any of the filtering so you are still using ps.bac, you never made a ps.bac.filt
#Prevalence of each feature and store as data.frame

prevdfITS <- apply(X = otu_table(psITS.bac),
                   MARGIN = ifelse(taxa_are_rows(psITS.bac), yes = 1, no = 2),
                   FUN = function(x){sum(x > 0)})

#Add taxonomy and total read counts to this data.frame
prevdfITS <- data.frame(Prevalence = prevdfITS,
                        TotalAbundance = taxa_sums(psITS.bac),
                        tax_table(psITS.bac))

#Average prevalence and total read counts of the features in Genera
plyr::ddply(prevdfITS, "Genus", function(df1){cbind(mean(df1$Prevalence),sum(df1$Prevalence))})

#Plot the Prevalence vs Total Abundance or Counts

ggplot(prevdfITS, aes(TotalAbundance, Prevalence / nsamples(psITS.bac),color=Genus)) +
  geom_hline(yintercept = 0.05, alpha = 0.5, linetype = 2) +  geom_point(size = 2, alpha = 0.7) + #Include a guess for parameter for prevalence treshold, here is 5%
  scale_x_log10() +  xlab("Total Abundance") + ylab("Prevalence [Frac. Samples]") +
  facet_wrap(~Genus) + theme(legend.position="none")

#Define a prevalence treshold based on the plot, here we used 5% of total samples
prevalenceThresholdITS <- 0.0000000001*nsamples(psITS.bac)

#### Plotting - Stacked Bar Plot ####

#Convert to simple probplots or relative abundance. 
ps.bac.normalizedITS <- transform_sample_counts(psITS.bac, function (x) {x/ sum(x)}) ###counts are divided by the total library size. 

#to make sure it worked
View(ps.bac.normalizedITS)

#Agglomerate taxa at the Genus level and convert to long format for plotting.
ps.bac.aggloITS <- tax_glom(ps.bac.normalizedITS, taxrank = "Genus") # for genus level

ps.bac.aggloITS_SP <- tax_glom(ps.bac.normalizedITS, taxrank = "Species") # for species level


ps.bac.agglo.longITS <- ps.bac.aggloITS %>% #this is a function of the dplyr package not R base
  tax_glom(taxrank = "Genus") %>%                     
  psmelt() %>%                                         
  arrange(Genus)

# Condense low abundance taxa into an "Other" category for the barplot; you use the ps.bac.agglo
# don't think you ever use the ps.bac.agglo


#Genera 
dataITS <- psmelt(ps.bac.aggloITS)
dataITS$Genus <- as.character(dataITS$Genus)
dataITS$Genus[dataITS$Abundance <0.01] <- "1% abund."
mediansITS <- ddply(dataITS, ~Genus, function(x) c(mediuan=median(x$Abundance)))
remainderITS <- mediansITS[mediansITS$mediuan <= 0.01,]$Genus
dataITS[dataITS$Genus %in% remainderITS,]$Genus <- "Genera < 1% abundance"
dataITS$Genus[dataITS$Abundance < 0.01] <- "Genera <1% abundance"

#Species 
dataITS_SP <- psmelt(ps.bac.aggloITS_SP)
dataITS_SP$Genus <- as.character(dataITS_SP$Genus)
dataITS_SP$Genus[dataITS_SP$Abundance <0.01] <- "10% abund."
mediansITS_SP <- ddply(dataITS_SP, ~Genus, function(x) c(mediuan=median(x$Abundance)))
remainderITS_SP <- mediansITS_SP[mediansITS_SP$mediuan <= 0.01,]$Genus
dataITS_SP[dataITS_SP$Genus %in% remainderITS_SP,]$Genus <- "Genera < 1% abundance"
dataITS_SP$Genus[dataITS_SP$Abundance < 0.01] <- "Genera <1% abundance"

#Create a stacked bar plot to observe community composition at the Genus Level.
colors <- c("#CBD588", "#5F7FC7", "orange","#DA5724", "#508578", "#CD9BCD",
            "#AD6F3B", "#673770","#D14285", "#652926", "#C84248","#8569D5", 
            "#5E738F","#D1A33D", "#8A7C64", "#599861","#999999", "#E69F00",
            "#56B4E9", "darkblue", "darkgoldenrod1", "darkseagreen", "darkorchid", 
            "darkolivegreen1", "lightskyblue", "darkgreen", "deeppink", "khaki2", "firebrick", 
            "brown1", "darkorange1", "cyan1", "royalblue4", "darksalmon", "darkblue",
            "royalblue4", "dodgerblue3", "steelblue1", "lightskyblue", "darkseagreen", 
            "darkgoldenrod1", "darkseagreen", "darkorchid", "darkolivegreen1", "brown1", 
            "darkorange1", "cyan1", "darkgrey", "aquamarine", "azure3", "bisque2", "cadetblue", "burlywood",
            "blue1", "azure4", "coral4", "darkorange1", "chocolate", "chartreuse3", "cornsilk4")

#this creates a very straight forward bar plot; but its combining the two time points
#if you change the following line: geom_bar(stat = "identity", position="fill") to geom_bar(stat = "identity", position="stack") it will plot total abundance instead of relative
#all you did was change the word fill to stack 

#Genera 

stackedITS <- GenusITS_Stacked_Bar_Plot <- ggplot(data = dataITS, aes(x = Stage, y = Abundance, fill = Genus)) +
  geom_bar(stat = "identity", position="fill") + #fill vs stack 
  theme_grey() + 
  theme (legend.position = "bottom") +
  scale_fill_manual(values = colors) +
  guides(fill = guide_legend(reverse = TRUE, keywidth = 0.6, keyheight = 0.6)) +
  ylab("Relative Abundance") +
  xlab("Developmental Stage") + theme(axis.text.x = element_text(angle =90, hjust = 1, vjust = 0.5)) +
  ggtitle("Genera of Cumulative Samples Across the Season") 
print(stackedITS)


# Subset by Management Style
org.sample <- filter(dataITS, Man_Style == "ORG")

con.sample <- filter(dataITS, Man_Style == "CON")


#Genera 
quartz(title="PCoA Endophyte Communities of Maize")
par(mfrow=c(2,1))

#organic
stackedITS_org <- GenusITS_Stacked_Bar_Plot <- ggplot(data = org.sample, aes(x = Stage, y = Abundance, fill = Genus)) +
  geom_bar(stat = "identity", position="fill") + #fill vs stack 
  theme_grey() + 
  theme (legend.position = "bottom") +
  scale_fill_manual(values = colors) +
  guides(fill = guide_legend(reverse = TRUE, keywidth = 0.6, keyheight = 0.6)) +
  ylab("Relative Abundance") +
  xlab("Developmental Stage") + theme(axis.text.x = element_text(angle =90, hjust = 1, vjust = 0.5)) +
  ggtitle("Genera of Organic") 
print(stackedITS_org)


#conventional
stackedITS_con <- GenusITS_Stacked_Bar_Plot <- ggplot(data = con.sample, aes(x = Stage, y = Abundance, fill = Genus)) +
  geom_bar(stat = "identity", position="fill") + #fill vs stack 
  theme_grey() + 
  theme (legend.position = "bottom") +
  scale_fill_manual(values = colors) +
  guides(fill = guide_legend(reverse = TRUE, keywidth = 0.6, keyheight = 0.6)) +
  ylab("Relative Abundance") +
  xlab("Developmental Stage") + theme(axis.text.x = element_text(angle =90, hjust = 1, vjust = 0.5)) +
  ggtitle("Genera of Conventional") 
print(stackedITS_con)

#Trying to combine plots into one figure
library(ggpubr)
ggpubr::ggarrange(stackedITS_org, stackedITS_con, ncol=2, nrow=1, common.legend = T, legend="bottom")


# Abundance for fusarium species only
fus.sample <- filter(dataITS_SP, Genus == "fusarium")

stackedITS_fu <- GenusITS_Stacked_Bar_Plot <- ggplot(data = fus.sample, aes(x = Stage, y = Abundance, fill = Species)) +
  geom_bar(stat = "identity", position="fill") + #fill vs stack 
  theme_grey() + 
  theme (legend.position = "bottom") +
  scale_fill_manual(values = colors) +
  guides(fill = guide_legend(reverse = TRUE, keywidth = 0.6, keyheight = 0.6)) +
  ylab("Relative Abundance") +
  xlab("Plot Type") + theme(axis.text.x = element_text(angle =90, hjust = 1, vjust = 0.5)) +
  ggtitle("Fusarium Species Identified Through ITS") 
print(stackedITS_fu)



### Top Taxa, See PERMANOVA code for subsetting

#Species 

#I care about developmental stage, so changing x to reflect that
stackedITS_SP <- SpeciesITS_Stacked_Bar_Plot <- ggplot(data = dataITS_SP, aes(x = Stage, y = Abundance, fill = Species)) +
  geom_bar(stat = "identity", position="fill") + #fill vs stack 
  theme_grey() + 
  theme (legend.position = "bottom") +
  scale_fill_manual(values = colors) +
  guides(fill = guide_legend(reverse = TRUE, keywidth = 0.6, keyheight = 0.6)) +
  ylab("Relative Abundance") +
  xlab("Leaf Sample") + theme(axis.text.x = element_text(angle =90, hjust = 1, vjust = 0.5)) +
  ggtitle("B") 
print(stackedITS_SP)

#### Diversity analysis ####

## 1. Alpha Diversity
library(ggpubr)
library(reshape2)
alpha_diversityITS <- estimate_richness(psITS.bac , measures= c("Observed", "Shannon", "InvSimpson"))

alpha_comparisonITS <- cbind(alpha_diversityITS, sample_data(psITS.bac )) 

melt_plotITS <- melt(alpha_comparisonITS) 
alpha_diversityITS

# future interest in trying hill diversity to address uneven sample sizes 
melt_plotITS


melt_plot_ShannonITS <- subset(melt_plotITS, variable =="Shannon") 
melt_plot_observedITS <- subset(melt_plotITS, variable =="Observed")
melt_plot_InvSimpsonITS <- subset(melt_plotITS, variable =="InvSimpson")

#Plotting Shannon, Observed, and InvSimpson 
aplot_alphaITS <- ggplot(data = melt_plot_ShannonITS, aes(x = factor(Stage), y=value, fill=Stage)) +
  geom_boxplot() +
  theme_minimal() + 
  xlab("Life Stage") + theme(axis.text.x = element_text(angle =90, hjust = 1, vjust = 0.5)) +
  scale_fill_brewer(palette="Dark2") + 
  ggtitle("Shannon Diversity") +
  stat_compare_means(label = "p.signif") # wilcox test (non-parametric test) for significance
aplot_alphaITS
ggsave("Shannon_DiversityITS.jpg", plot = aplot_alphaITS, width = 10, height=7, units= "in", dpi =1000)

bplot_alphaITS <- ggplot(data = melt_plot_observedITS, aes(x = factor(Stage), y=value, fill=Stage)) +
  geom_boxplot() +
  theme_minimal() + 
  xlab("Life Stage") + theme(axis.text.x = element_text(angle =90, hjust = 1, vjust = 0.5)) +
  scale_fill_brewer(palette="Dark2") + 
  ggtitle("Observed Diversity") +
  stat_compare_means(label = "p.signif") # wilcox test (non-parametric test) for significance
bplot_alphaITS
ggsave("Observed_DiversityITS .jpg", plot = bplot_alphaITS, width = 10, height=7, units= "in", dpi =1000)

cplot_alphaITS <- ggplot(data = melt_plot_InvSimpsonITS, aes(x = factor(Stage), y=value, fill=Stage)) +
  geom_boxplot() +
  theme_minimal() + 
  xlab("Life Stage") + theme(axis.text.x = element_text(angle =90, hjust = 1, vjust = 0.5)) +
  scale_fill_brewer(palette="Dark2") + 
  ggtitle("Inverted Simpson Diversity") +
  stat_compare_means(label = "p.signif") # wilcox test (non-parametric test) for significance
cplot_alphaITS
ggsave("InvSimpson_DiversityITS .jpg", plot = cplot_alphaITS, width = 10, height=7, units= "in", dpi =1000)
 

#### Multivariate Community Analysis (PCoA, PERMANOVA, etc.) ####

###PCoA

#First, we fix sample_ID in sample data and then we ordinate at the samples level using PCoA and Bray-Curtis distances.  
#sample_data is phyloseq language not the name of anything you created
sample_data(ps.bac.normalizedITS)$Stage <- factor(sample_data(ps.bac.normalizedITS)$Stage, levels = c("ES", "R3", "R5", "SC", "Seed", "V5", "V12", "VT"))

#Unconstrained Ordination and plot ordination.
PCoAITS <- ordinate (ps.bac.normalizedITS, method =  "PCoA", distance = "bray")

PCoA_PlotITS <- plot_ordination(physeq = ps.bac.normalizedITS, ordination = PCoAITS, color = "Stage", title = "Endophyte Communities of Maize Over Time"
) + scale_fill_manual (values = colors) + geom_point(aes(color = Stage), alpha = 0.7, size =4) + theme_bw()

print(PCoA_PlotITS)
PCoA_PlotITS + stat_ellipse(aes(lty=factor(Stage)), type = "t", level = 0.95) + scale_linetype_manual(values=c(1,2,3,4,5,6,7,8,9)) + theme_bw()

# Repeat for plot to find outliers 
PCoA_PlotITS <- plot_ordination(physeq = ps.bac.normalizedITS, ordination = PCoAITS, color = "Plot", title = "Endophyte Communities of Maize Over Time by Plot"
) + scale_fill_manual (values = colors) + geom_point(aes(color = Plot), alpha = 0.7, size =4) + theme_bw()

print(PCoA_PlotITS)
PCoA_PlotITS + stat_ellipse(aes(lty=factor(Stage)), type = "t", level = 0.95) + scale_linetype_manual(values=c(1,2,3,4,5,6,7,8,9)) + theme_bw()

##Repeat for organic and conventional plots
sample_data(ps.bac.normalizedITS)$Man_Style <- factor (sample_data(ps.bac.normalizedITS)$Man_Style, levels = c("CON", "ORG"))

#Unconstrained Ordination and plot ordination.
PCoAITS <- ordinate (ps.bac.normalizedITS, method =  "PCoA", distance = "bray")

PCoA_PlotITS <- plot_ordination(physeq = ps.bac.normalizedITS, ordination = PCoAITS, color = "Man_Style",
) + scale_fill_manual (values = colors) + geom_point(aes(color = Man_Style), alpha = 0.7, size =4) + theme_bw()

print(PCoA_PlotITS)
PCoA_PlotITS + stat_ellipse(aes(lty=factor(Man_Style)), type = "t", level = 0.95) + scale_linetype_manual(values=c(1,2,3,4,5,6,7,8,9)) + theme_bw()

# Note to senf: This code stops where Andrews "Microbiome_Analysis-4.R" code begins the PERMANOVA analysis

# We want to know how much of each taxanomic level we have for our results
# to do this we will try the following'
get_taxa_unique(psITS.bac, taxonomic.rank = "Genus")


