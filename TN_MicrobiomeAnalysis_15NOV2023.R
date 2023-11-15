### Propject Title: Understanding the selective capability of maize developmental stage on fungal endophyte communities 
### Author: Chelsea Newbold - cln68@psu.edu
### Course: ENT 535 FALL 2023 - Penn State

### This code was written and run on a Windows computer

# Download and Install Necessary Packages (this is not all packages, other packages noted where relvent)

# TN- Installed on PC, 24 July 2023
install.packages("rlang") #Phyloseq requires specific version of rlang; need to download most current version to properly use phyloseq downstream 
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c("phyloseq","Biostrings","RDPutils","ggplot2","dplyr", "ggpubr", "vegan", "emmeans", "plyr"))

#Load Packages
library(phyloseq)
library(Biostrings)
library(ggplot2)
library(dplyr)
library(ggpubr)
library(vegan)
library(emmeans)
library(plyr)
# install.packages("tidyr")
library(tidyr)
# BiocManager::install("microbiome")
library(microbiome)
library(vegan)
library(tidyselect)
install.packages("viridis") # for relative abundance plots mostly
library(viridis)
library(reshape2)

#### Bring in OTU Table, Taxonomy Table and Metadata table from Github ####

# Import Metadata File
# if you have a hard time finding the path, use your terminal and drag/drop the file
meta_d <- read.csv("https://raw.githubusercontent.com/clnPSU2023/Fall2022_EndoMaizeTime/main/MetaData_NegConRemove_23Jun23.csv",row.names = 1) 

#Metadata structure.
str(meta_d)
head(meta_d, 10)

#changing to factor
meta_d$Stage <- as.factor(meta_d$Stage)
meta_d$Man_Style <- as.factor(meta_d$Man_Style)
meta_d$Plot_Type <- as.factor(meta_d$Plot_Type)
meta_d$Plot <- as.factor(meta_d$Plot)
meta_d$Num_Stage <- as.factor(meta_d$Num_Stage)

str(meta_d)

#Taxa and OTU tables
#read in taxa table
#the end part (row.names = 1) is very important. It removes column one from the actual data 
taxa_d <- read.csv("https://raw.githubusercontent.com/clnPSU2023/Fall2022_EndoMaizeTime/main/Taxa_Table.csv", row.names = 1)
str(taxa_d)

#as it is, the taxa table is a dataframe but it needs to be a matrix in order to merge into a phyloseq object
taxaITS <- as.matrix(taxa_d)
class(taxaITS)

#read in OTU Table
#again, the end part (row.names = 1) is very important. It removes column one from the actual data 
#this is important here because the first column in not numeric and the rest of the table needs to be
OTU_d <-read.csv("https://raw.githubusercontent.com/clnPSU2023/Fall2022_EndoMaizeTime/main/OTU.csv",row.names=1 )

#you will see that it reads the counts in as integers, you must change this to numeric
str(OTU_d)
#find out the number of columns 
ncol(OTU_d)

#create a subset of the data OTU[1:85] then use lapply to convert it to 
#numeric and then assign it back to those columns
OTU_d[1:85] <- lapply(OTU_d[1:85], as.numeric)
str(OTU_d)

# transpose OTU table
otus.t <- t(OTU_d)

# Now create a merged phyloseq object for all three tables; meta, taxa, and otu
psITS <- phyloseq(otu_table(otus.t, taxa_are_rows=FALSE), tax_table(taxaITS)) # first just taxa and otu
meta <- sample_data(meta_d) # sample level data, prepare for pairing in phyloseq 
psITS.bac <- merge_phyloseq(psITS, meta)# now add the metadata

#Duplicate the phyloseq object to preserve original
psITS.bac1 <- psITS.bac #working phyloseq object

#Order samples by the number of reads
ordered(sample_sums(psITS.bac1)) # this is helpful when discussing data before processing 

### Note: Save where you want the RDS object saved on your computer, see code for saving here! 
#Save RDS and Read RDS
# psITS.bac1
# saveRDS(psITS.bac1, "MR_DNA_Fungi_Data_order.rds")

# Bring in RDS if saved and working on at a later date, see code here! 
# Read .RDS By using <-, we are restoring the RDS under the same name we used for our Phyloseq object.
# psITS.bac <- readRDS("C:/Users/16262/OneDrive/Desktop/MrDNA_Pipeline/MR_DNA_Fungi_Data_order.rds") # error appeared but looks to have worked anyway???


#### Filtering, Subsetting, and Performing Richness, Ordination, and Multivariate analyses ####

# Remember our basic questions:
# How do the communities change over time?
# How do the communities differ between conventional and organic sites? 
# How do the communities change between cover crop treatments? 

#Verify the prevalence of each taxon, in other words, the number of samples in which a taxon appears at least once.

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

#Plot the Prevalence vs Total Abundance 
geom_col() # must run this line of code to get this to work on Windows PC

ggplot(prevdfITS, aes(TotalAbundance, Prevalence / nsamples(psITS.bac),color=Genus)) +
  geom_hline(yintercept = 0.05, alpha = 0.5, linetype = 2) +  geom_point(size = 2, alpha = 0.7) + #Include a guess for parameter for prevalence treshold, here is 5%
  scale_x_log10() +  xlab("Total Abundance") + ylab("Prevalence [Frac. Samples]") +
  facet_wrap(~Genus) + theme(legend.position="none")
# By plotting the abundace you get a good idea for which genera are most abundant before moving forward with visualization

### Plotting - Stacked Bar Plot

#Convert to simple probability plots or relative abundance. 
ps.bac.normalizedITS <- transform_sample_counts(psITS.bac, function (x) {x/ sum(x)}) ###counts are divided by the total library size. 

#to make sure it worked
View(ps.bac.normalizedITS)

#Agglomerate taxa at the Genus and Species level and convert to long format for plotting.
ps.bac.aggloITS <- tax_glom(ps.bac.normalizedITS, taxrank = "Genus")

ps.bac.aggloITS_SP <- tax_glom(ps.bac.normalizedITS, taxrank = "Species")


ps.bac.agglo.longITS <- ps.bac.aggloITS %>% #this is a function of the dplyr package not R base
  tax_glom(taxrank = "Genus") %>%                     
  psmelt() %>%                                         
  arrange(Genus)

#Condense low abundance taxa into an "Other" category for the barplot; you use the ps.bac.aggloITS

#Genera 

dataITS <- psmelt(ps.bac.aggloITS) # combines meta and taxa data into one  dataframe
dataITS$Genus <- as.character(dataITS$Genus) # make genus a character
dataITS$Genus[dataITS$Abundance <0.01] <- "1% abund." # low abudnance reads are those <1%
mediansITS <- ddply(dataITS, ~Genus, function(x) c(mediuan=median(x$Abundance)))
remainderITS <- mediansITS[mediansITS$mediuan <= 0.01,]$Genus
dataITS[dataITS$Genus %in% remainderITS,]$Genus <- "Genera < 1% abundance"
dataITS$Genus[dataITS$Abundance < 0.01] <- "Genera <1% abundance"

#Species 

dataITS_SP <- psmelt(ps.bac.aggloITS_SP) # combines meta and taxa data into one  dataframe
dataITS_SP$Species <- as.character(dataITS_SP$Species) # make genus a character
dataITS_SP$Species[dataITS_SP$Abundance <0.01] <- "1% abund"  # low abudnance reads are those <1%
mediansITS_SP <- ddply(dataITS_SP, ~Species, function(x) c(mediuan=median(x$Abundance)))
remainderITS_SP <- mediansITS_SP[mediansITS_SP$mediuan <= 0.01,]$Species
dataITS_SP[dataITS_SP$Species %in% remainderITS_SP,]$Species <- "Species < 1% abundance"
dataITS_SP$Species[dataITS_SP$Abundance < 0.01] <- "Species <1% abundance"


### Species - Relative Abundance ###

# We need a lot of colors for microbiome data visualization, especially colors that can be seen distinctly, viridis is great for this,
# or at least better than other options I have used
#install.packages("viridis") # if not done already
#library(viridis) # if not brought in earlier

# realtive abundance for total species for conventional and organic combined
stackedITS_sp <- GenusITS_Stacked_Bar_Plot <- ggplot(data = dataITS_SP, aes(x = Man_Style, y = Abundance, fill = Species)) +
  geom_bar(stat = "identity", position="fill") + #fill vs stack 
  theme_grey() + 
  theme (legend.position = "bottom") +
  scale_color_viridis() + # This has been the best option so far!!!
  guides(fill = guide_legend(reverse = TRUE, keywidth = 0.6, keyheight = 0.6)) +
  ylab("Relative Abundance") +
  xlab("Management Style") + theme(axis.text.x = element_text(angle =90, hjust = 1, vjust = 0.5)) +
  ggtitle("Species of Organic and Covnentional Sites") 
print(stackedITS_sp)

# realtive abudnance for species over time, total conventional and orgnaic
dataITS_names <- dataITS_SP %>%
  mutate(Num_Stage = recode(Num_Stage, 
                            "1" = "Planting - Seed",
                            "2" = "Harvest - Seed",
                            "3" = "V5.",
                            "4" = "V12", "5" = "VT", "6" = "R3", "7" = "R5", "8" = "Scenescent")) # change names

stackedITS_sp_name <- GenusITS_Stacked_Bar_Plot <- ggplot(data = dataITS_names, aes(x = Num_Stage, y = Abundance, fill = Species)) +
  geom_bar(stat = "identity", position="fill") + #fill vs stack 
  theme_grey() + 
  theme (legend.position = "bottom") +
  scale_color_viridis() +
  guides(fill = guide_legend(reverse = TRUE, keywidth = 0.6, keyheight = 0.6)) +
  ylab("Relative Abundance") +
  xlab("Developmental Stage") + theme(axis.text.x = element_text(angle =90, hjust = 1, vjust = 0.5)) +
  ggtitle("Species of Cumulative Samples Across the Season") 
print(stackedITS_sp_name)

# Subset by management styles
sp.org.sample <- filter(dataITS_names, Man_Style == "ORG")

sp.con.sample <- filter(dataITS_names, Man_Style == "CON")

# create relative abundance plots for species in organic and conventional plots seperately

#organic
stackedITS_org_sp <- GenusITS_Stacked_Bar_Plot <- ggplot(data = sp.org.sample, aes(x = Num_Stage, y = Abundance, fill = Species)) +
  geom_bar(stat = "identity", position="fill") + #fill vs stack 
  theme_grey() + 
  theme (legend.position = "bottom") +
  scale_color_viridis()  +
  guides(fill = guide_legend(reverse = TRUE, keywidth = 0.6, keyheight = 0.6)) +
  ylab("Relative Abundance") +
  xlab("Developmental Stage") + theme(axis.text.x = element_text(angle =90, hjust = 1, vjust = 0.5)) +
  ggtitle("Topp 30 Species from Organic Samples") 
print(stackedITS_org_sp)

#conventional
stackedITS_con_sp <- GenusITS_Stacked_Bar_Plot <- ggplot(data = sp.con.sample, aes(x = Num_Stage, y = Abundance, fill = Species)) +
  geom_bar(stat = "identity", position="fill") + #fill vs stack 
  theme_grey() + 
  theme (legend.position = "bottom") +
  scale_color_viridis() +
  guides(fill = guide_legend(reverse = TRUE, keywidth = 0.6, keyheight = 0.6)) +
  ylab("Relative Abundance") +
  xlab("Developmental Stage") + theme(axis.text.x = element_text(angle =90, hjust = 1, vjust = 0.5)) +
  ggtitle("Top 30 Species from Conventional Samples") 
print(stackedITS_con_sp)


### Genera - relative Abundance ###

# This will order the maize developmental stages as 1 through 8, with 1 and 2 being planting and harvest
# and the 3 through 8 showing V5, V12, VT, R3, R5 and SC respectively 
dataITS_shorter <- dataITS %>%
  mutate(Num_Stage = recode(Num_Stage, 
                            "1" = "Planting - Seed",
                            "2" = "Harvest - Seed",
                            "3" = "V5.",
                            "4" = "V12", "5" = "VT", "6" = "R3", "7" = "R5", "8" = "Scenescent")) # change names


# run relative abundance stacked bar plot for JUST organic and conventional with no developmental stage
stackedITS_TT <- GenusITS_Stacked_Bar_Plot <- ggplot(data = dataITS_shorter, aes(x = Man_Style, y = Abundance, fill = Genus)) +
  geom_bar(stat = "identity", position="fill") + #fill vs stack 
  theme_grey() + 
  theme (legend.position = "bottom") +
  scale_color_viridis() +
  guides(fill = guide_legend(reverse = TRUE, keywidth = 0.6, keyheight = 0.6)) +
  ylab("Relative Abundance") +
  xlab("Management Style") + theme(axis.text.x = element_text(angle =90, hjust = 1, vjust = 0.5)) +
  ggtitle("Genera of Organic and Covnentional Sites") 
print(stackedITS_TT)

## now do this over developmental stages

stackedITS <- GenusITS_Stacked_Bar_Plot <- ggplot(data = dataITS_shorter, aes(x = Num_Stage, y = Abundance, fill = Genus)) +
  geom_bar(stat = "identity", position="fill") + #fill vs stack 
  theme_grey() + 
  theme (legend.position = "bottom") +
  scale_color_viridis()  +
  guides(fill = guide_legend(reverse = TRUE, keywidth = 0.6, keyheight = 0.6)) +
  ylab("Relative Abundance") +
  xlab("Developmental Stage") + theme(axis.text.x = element_text(angle =90, hjust = 1, vjust = 0.5)) +
  ggtitle("Genera of Cumulative Samples Across the Season") 
print(stackedITS)

# Subset by Management Style
org.sample <- filter(dataITS_shorter, Man_Style == "ORG")

con.sample <- filter(dataITS_shorter, Man_Style == "CON")

# create relative abundance plots for genera in organic and conventional plots seperately

#organic
stackedITS_org <- GenusITS_Stacked_Bar_Plot <- ggplot(data = org.sample, aes(x = Num_Stage, y = Abundance, fill = Genus)) +
  geom_bar(stat = "identity", position="fill") + #fill vs stack 
  theme_grey() + 
  theme (legend.position = "bottom") +
  scale_color_viridis()  +
  guides(fill = guide_legend(reverse = TRUE, keywidth = 0.6, keyheight = 0.6)) +
  ylab("Relative Abundance") +
  xlab("Developmental Stage") + theme(axis.text.x = element_text(angle =90, hjust = 1, vjust = 0.5)) +
  ggtitle("Top 20 Genera from Organic Samples") 
print(stackedITS_org)

#conventional
stackedITS_con <- GenusITS_Stacked_Bar_Plot <- ggplot(data = con.sample, aes(x = Num_Stage, y = Abundance, fill = Genus)) +
  geom_bar(stat = "identity", position="fill") + #fill vs stack 
  theme_grey() + 
  theme (legend.position = "bottom") +
  scale_color_viridis()  +
  guides(fill = guide_legend(reverse = TRUE, keywidth = 0.6, keyheight = 0.6)) +
  ylab("Relative Abundance") +
  xlab("Developmental Stage") + theme(axis.text.x = element_text(angle =90, hjust = 1, vjust = 0.5)) +
  ggtitle("Top 20 Genera from Conventional Samples") 
print(stackedITS_con)

# Run relative abundance plots by cover crop for organic data (there are no cover crops in the coventional)

# if not done already subset organic data, see line 220
# org.sample <- filter(dataITS_shorter, Man_Style == "ORG")

# filter by cover crop and fallow
f.sample <- filter(org.sample, Plot_Type == "Fallow")
c.sample <- filter(org.sample, Plot_Type == "Clover")
r.sample <- filter(org.sample, Plot_Type == "Radish")
t.sample <- filter(org.sample, Plot_Type == "Triticale")

f_plot <- GenusITS_Stacked_Bar_Plot <- ggplot(data = f.sample, aes(x = Num_Stage, y = Abundance, fill = Genus)) +
  geom_bar(stat = "identity", position="fill") + #fill vs stack 
  theme_grey() + 
  theme (legend.position = "bottom") +
  scale_color_viridis()  +
  guides(fill = guide_legend(reverse = TRUE, keywidth = 0.6, keyheight = 0.6)) +
  ylab("Relative Abundance") +
  xlab("Developmental Stage") + theme(axis.text.x = element_text(angle =90, hjust = 1, vjust = 0.5)) +
  ggtitle("Top 20 Genera from Fallow Plot Samples") 
print(f_plot)

c_plot <- GenusITS_Stacked_Bar_Plot <- ggplot(data = c.sample, aes(x = Num_Stage, y = Abundance, fill = Genus)) +
  geom_bar(stat = "identity", position="fill") + #fill vs stack 
  theme_grey() + 
  theme (legend.position = "bottom") +
  scale_color_viridis()  +
  guides(fill = guide_legend(reverse = TRUE, keywidth = 0.6, keyheight = 0.6)) +
  ylab("Relative Abundance") +
  xlab("Developmental Stage") + theme(axis.text.x = element_text(angle =90, hjust = 1, vjust = 0.5)) +
  ggtitle("Top 20 Genera from Clover Plot Samples") 
print(c_plot)

r_plot <- GenusITS_Stacked_Bar_Plot <- ggplot(data = r.sample, aes(x = Num_Stage, y = Abundance, fill = Genus)) +
  geom_bar(stat = "identity", position="fill") + #fill vs stack 
  theme_grey() + 
  theme (legend.position = "bottom") +
  scale_color_viridis()  +
  guides(fill = guide_legend(reverse = TRUE, keywidth = 0.6, keyheight = 0.6)) +
  ylab("Relative Abundance") +
  xlab("Developmental Stage") + theme(axis.text.x = element_text(angle =90, hjust = 1, vjust = 0.5)) +
  ggtitle("Top 20 Genera from Radish Plot Samples") 
print(r_plot)

t_plot <- GenusITS_Stacked_Bar_Plot <- ggplot(data = t.sample, aes(x = Num_Stage, y = Abundance, fill = Genus)) +
  geom_bar(stat = "identity", position="fill") + #fill vs stack 
  theme_grey() + 
  theme (legend.position = "bottom") +
  scale_color_viridis()  +
  guides(fill = guide_legend(reverse = TRUE, keywidth = 0.6, keyheight = 0.6)) +
  ylab("Relative Abundance") +
  xlab("Developmental Stage") + theme(axis.text.x = element_text(angle =90, hjust = 1, vjust = 0.5)) +
  ggtitle("Top 20 Genera from Triticale Plot Samples") 
print(t_plot)

# To create a combined plot with shared legend and shared axes add "facet_grid(rows = vars(Plot_Type)) +"

org.sample_CC <- org.sample %>%
  mutate(Plot_Type = recode(Plot_Type, 
                            "Fallow" = "Fallow",
                            "Clover" = "Clover",
                            "Triticale" = "Triticale",
                            "Radish" = "Radish", "ORG" = "Organic Seed"))

stackedITS_org_stack <- GenusITS_Stacked_Bar_Plot <- ggplot(data = org.sample_CC, aes(x = Num_Stage, y = Abundance, fill = Genus)) +
  geom_bar(stat = "identity", position="fill") + #fill vs stack 
  theme_grey() + 
  theme (legend.position = "bottom") +
  scale_color_viridis() +
  guides(fill = guide_legend(reverse = TRUE, keywidth = 0.6, keyheight = 0.6)) +
  facet_grid(rows = vars(Plot_Type)) +
  ylab("Relative Abundance") +
  xlab("Developmental Stage") + theme(axis.text.x = element_text(angle =90, hjust = 1, vjust = 0.5)) +
  ggtitle("Relative Abundance of Top 20 Genera in Organic Samples Across Cover Crops") 
print(stackedITS_org_stack)

### Diversity analyses ###

## 1. Alpha Diversity Metrics observed, shannon, and simpson

alpha_diversityITS <- estimate_richness(psITS.bac , measures= c("Observed", "Shannon", "InvSimpson"))

alpha_comparisonITS <- cbind(alpha_diversityITS, sample_data(psITS.bac )) # add alspha diversity 

melt_plotITS <- melt(alpha_comparisonITS)

alpha_diversityITS # view diversit values 
melt_plotITS # see plot

# reassign names so that they appear in order, but note for some reason this does not always work
melt_plotITS_shorter <- melt_plotITS %>%
  mutate(Num_Stage = recode(Stage, 
                            "Seed" = "Planting - Seed",
                            "ES" = "Harvest - Seed",
                            "V5" = "V5",
                            "V12" = "V12", "VT" = "VT", "R3" = "R3", "R5" = "R5", "SC" = "Senescent"))

# seperate out the diversity metrics 
melt_plot_ShannonITS <- subset(melt_plotITS_shorter, variable =="Shannon") 
melt_plot_observedITS <- subset(melt_plotITS_shorter, variable =="Observed")
melt_plot_InvSimpsonITS <- subset(melt_plotITS_shorter, variable =="InvSimpson")

## Note: export to your PC location of choice, see code here! 
# export to pc
# basic code -> write.csv(df, file = "export_csv")
# write.csv(melt_plotITS_shorter, file = "C:/Users/16262/OneDrive/Desktop/TN_DiversityMetrics_2022.csv")

# Plotting Shannon, Observed, and InvSimpson

# start by comparing Shannon's between organic and conventional over all time points
aplot_alphaITS_TT <- ggplot(data = melt_plot_ShannonITS, aes(x = factor(Man_Style), y=value, fill=Man_Style)) +
  geom_boxplot() +
  theme_minimal() + 
  xlab("Man_Style") + theme(axis.text.x = element_text(angle =90, hjust = 1, vjust = 0.5)) +
  scale_fill_brewer(palette="Dark2") + 
  ggtitle("Shannon Diversity of Organic and Conventional Sites") +
  stat_compare_means(label = "p.signif") # wilcox test (non-parametric test) for significance
aplot_alphaITS_TT

bplot_alphaITS_TT <- ggplot(data = melt_plot_observedITS, aes(x = factor(Man_Style), y=value, fill=Man_Style)) +
  geom_boxplot() +
  theme_minimal() + 
  xlab("Management Style") + theme(axis.text.x = element_text(angle =90, hjust = 1, vjust = 0.5)) +
  scale_fill_brewer(palette="Dark2") + 
  ggtitle("Observed Diversity of Organic and Conventional Sites") +
  stat_compare_means(label = "p.signif") # wilcox test (non-parametric test) for significance
bplot_alphaITS_TT

cplot_alphaITS_TT <- ggplot(data = melt_plot_InvSimpsonITS <- subset(melt_plotITS_shorter, variable =="InvSimpson")
, aes(x = factor(Man_Style), y=value, fill=Man_Style)) +
  geom_boxplot() +
  theme_minimal() + 
  xlab("Management Style") + theme(axis.text.x = element_text(angle =90, hjust = 1, vjust = 0.5)) +
  scale_fill_brewer(palette="Dark2") + 
  ggtitle("Simpson Diversity of Organic and Conventional Sites") +
  stat_compare_means(label = "p.signif") # wilcox test (non-parametric test) for significance
cplot_alphaITS_TT

# run diversity metrics separately for organic and conventional

org.div <- filter(melt_plotITS_shorter, Man_Style == "ORG")
con.div <- filter(melt_plotITS_shorter, Man_Style == "CON")

# For organic

org_melt_plot_ShannonITS <- subset(org.div, variable =="Shannon") 
org_melt_plot_observedITS <- subset(org.div, variable =="Observed")
org_melt_plot_InvSimpsonITS <- subset(org.div, variable =="InvSimpson")

# alpha diversity for organic 
org.div <- filter(melt_plotITS_shorter, Man_Style == "ORG")
con.div <- filter(melt_plotITS_shorter, Man_Style == "CON")

# run diversity metrics 
org_melt_plot_ShannonITS <- subset(org.div, variable =="Shannon") 
org_melt_plot_observedITS <- subset(org.div, variable =="Observed")
org_melt_plot_InvSimpsonITS <- subset(org.div, variable =="InvSimpson")

#Plotting Shannon, Observed, and InvSimpson for organic; a, b, and c resepctively 
aplot_alphaITS_org <- ggplot(data = org_melt_plot_ShannonITS, aes(x = factor(Num_Stage), y=value, fill=Num_Stage)) +
  geom_boxplot() +
  theme_minimal() + 
  xlab("Life Stage") + theme(axis.text.x = element_text(angle =90, hjust = 1, vjust = 0.5)) +
  scale_fill_brewer(palette="Dark2") + 
  ggtitle("Shannon Diversity - Organic") +
  stat_compare_means(label = "p.signif") # wilcox test (non-parametric test) for significance
aplot_alphaITS_org

bplot_alphaITS_org <- ggplot(data = org_melt_plot_observedITS, aes(x = factor(Num_Stage), y=value, fill=Num_Stage)) +
  geom_boxplot() +
  theme_minimal() + 
  xlab("Life Stage") + theme(axis.text.x = element_text(angle =90, hjust = 1, vjust = 0.5)) +
  scale_fill_brewer(palette="Dark2") + 
  ggtitle("Observed Diversity - Organic") +
  stat_compare_means(label = "p.signif") # wilcox test (non-parametric test) for significance
bplot_alphaITS_org

cplot_alphaITS_org <- ggplot(data = org_melt_plot_InvSimpsonITS, aes(x = factor(Num_Stage), y=value, fill=Num_Stage)) +
  geom_boxplot() +
  theme_minimal() + 
  xlab("Life Stage") + theme(axis.text.x = element_text(angle =90, hjust = 1, vjust = 0.5)) +
  scale_fill_brewer(palette="Dark2") + 
  ggtitle("Inverted Simpson Diversity - Organic") +
  stat_compare_means(label = "p.signif") # wilcox test (non-parametric test) for significance
cplot_alphaITS_org

#Do the same for Conventional

con_melt_plot_ShannonITS <- subset(con.div, variable =="Shannon") 
con_melt_plot_observedITS <- subset(con.div, variable =="Observed")
con_melt_plot_InvSimpsonITS <- subset(con.div, variable =="InvSimpson")

#Plotting Shannon, Observed, and Simpsons for Conventional
aplot_alphaITS_con <- ggplot(data = con_melt_plot_ShannonITS, aes(x = factor(Num_Stage), y=value, fill=Num_Stage)) +
  geom_boxplot() +
  theme_minimal() + 
  xlab("Life Stage") + theme(axis.text.x = element_text(angle =90, hjust = 1, vjust = 0.5)) +
  scale_fill_brewer(palette="Dark2") + 
  ggtitle("Shannon Diversity - Conventional") +
  stat_compare_means(label = "p.signif") # wilcox test (non-parametric test) for significance
aplot_alphaITS_con

bplot_alphaITS_con <- ggplot(data = con_melt_plot_observedITS, aes(x = factor(Num_Stage), y=value, fill=Num_Stage)) +
  geom_boxplot() +
  theme_minimal() + 
  xlab("Life Stage") + theme(axis.text.x = element_text(angle =90, hjust = 1, vjust = 0.5)) +
  scale_fill_brewer(palette="Dark2") + 
  ggtitle("Observed Diversity - Conventional") +
  stat_compare_means(label = "p.signif") # wilcox test (non-parametric test) for significance
bplot_alphaITS_con

cplot_alphaITS_con <- ggplot(data = con_melt_plot_InvSimpsonITS, aes(x = factor(Num_Stage), y=value, fill=Num_Stage)) +
  geom_boxplot() +
  theme_minimal() + 
  xlab("Life Stage") + theme(axis.text.x = element_text(angle =90, hjust = 1, vjust = 0.5)) +
  scale_fill_brewer(palette="Dark2") + 
  ggtitle("Inverted Simpson Diversity - Conventional") +
  stat_compare_means(label = "p.signif") # wilcox test (non-parametric test) for significance
cplot_alphaITS_con

#### Multivariate Analyses ####

### PCoA - nondimentional scaline, aobserve dissimilarity between stages and management ###

# remember that this step requires ggplot and works with phyloseq object

sample_data(ps.bac.normalizedITS)$Stage <- factor(sample_data(ps.bac.normalizedITS)$Stage, levels = c("ES", "R3", "R5", "SC", "Seed", "V5", "V12", "VT"))


#Run an unconstrained Ordination and plot the ordination
PCoAITS <- ordinate (ps.bac.normalizedITS, method =  "PCoA", distance = "bray")

PCoA_PlotITS <- plot_ordination(physeq = ps.bac.normalizedITS, 
                                ordination = PCoAITS, 
                                color = "Stage", 
                                title = "Endophyte Communities of Maize Over Time") + # not needed for final paper
  scale_fill_manual (values = colors) + 
  geom_point(aes(color = Stage), alpha = 0.7, size =4) + 
  theme_bw()

PCoA_PlotITS + 
  stat_ellipse(aes(lty=factor(Stage)), type = "t", level = 0.95, show.legend = F) + 
  scale_linetype_manual(values=c(1,2,3,4,5,6,7,8, 9, 10)) + theme_grey() + # adds 95% confidence elipses without legend
  xlab("PCoA 1 (40.3%)") + ylab("PCoA 2 (14.1%)") # correct labels


### PERMANOVA ###

## Ratinoal:
# the distrobution of abudnance of ASVs/ taxa is skewed because we have a lot of 
# low abudnance ASVs
ggplot(data = dataITS_SP, 
       aes(x = Abundance)) +
  geom_histogram(aes(fill = Stage))
# The permanova is a nonparametric test, ie it odes not worry about distrobution 
# the test allows me to compare groups based on distance from the centroid
# for my data this is best accomplished using the dissimilarity matrics, such as bray-curtis 

# original tutorial I used is HERE -> https://github.com/joey711/phyloseq/issues/1046

ps.prop_allvend_ctrl <- transform_sample_counts(psITS.bac1, function(otu) {otu/sum(otu)})
set.seed(134)

# calculating bray curtis distance matrix
PCoAITS_bray <- phyloseq::distance(ps.prop_allvend_ctrl, method = "bray")

# make sample data
sampledf <- data.frame(sample_data(psITS.bac1))

# run adonis, pertmanova for development stage and management style
perm_stage <- adonis2(PCoAITS_bray ~ Num_Stage + Man_Style, data = sampledf)
perm_stage # view  output
# significant effect of developmental stage and management style on community dissimilarity/ composition
# num_stage p val = 0.001
# man_style p val = 0.004

# run adonis for both stage and management variables
set.seed(134) #reproducible results
perm_total <-  adonis2(PCoAITS_bray ~ Man_Style * Num_Stage, data = sampledf)
perm_total # view  output
# same result as perm_stage but the interaction term is also signficant, p cal = 0.001

#PERMANOVA dispersion (Code from Jared!)

# if not already open, open these libraries 
library(vegan)
library(tidyselect)
library(ggplot2)

end.res.adisper_tot <- betadisper(d = PCoAITS_bray, group = sampledf$Stage, type = 'centroid')
end.res.adisper_tot$distances

anova(end.res.adisper_tot)
TukeyHSD(end.res.adisper_tot) # this is our post hoc comparison across variables 
# see written results for specific comaprisons
# in sumamry there are no trends of biologica relevence
# except, senescent leaves are very different froom the rest of the developmental stages 

boxplot(end.res.adisper_tot)
plot(end.res.adisper_tot)

# I also want to run permanova for the organic seperately looking at cover crop treatments, n=2

# subset organic data
organic_df <- filter(dataITS, Man_Style == "ORG")
# subset data for abundances over 0%
organic_df_zro <- filter(organic_df, Abundance > 0)

# bray curtis 
org.dist<-vegdist(organic_df_zro$Abundance, method='bray')
org.dist

# start with organic by stage before cover crops
# permanova for org 

set.seed(134) #reproducible results
org.div<-adonis2(org.dist~organic_df_zro$Stage, data=organic_df_zro, permutations=9999)
org.div
#PERMANOVA dispersion
org.res.adisper <- betadisper(d = org.dist, group = organic_df_zro$Stage, type = 'centroid')

org.res.adisper$distances

anova(org.res.adisper)
TukeyHSD(org.res.adisper) # post hok comparisons 
# see result sectino for comparison results

boxplot(org.res.adisper)
plot(org.res.adisper)

# cover crops next

# permanova for cover crops, column called plot_type
set.seed(134) #reproducible results
org.div.ccc<-adonis2(org.dist~organic_df_zro$Plot_Type, data=organic_df_zro, permutations=9999)
org.div.ccc

# PERMANOVA dispersion
org.res.adisper <- betadisper(d = org.dist, group = organic_df$Stage, type = 'centroid')

org.res.adisper$distances

anova(org.res.adisper)
TukeyHSD(org.res.adisper)# this is our post hoc comparison across variables 
# see result sectino for comparison results

boxplot(org.res.adisper)
plot(org.res.adisper)

# in summary there is absolutely no differences in community composition between cover crop treatments 

# permanova for conventional samples looking at developmental stage 
# use the subset data for con now
conv_df <- filter(dataITS, Man_Style == "CON")
# subset data for abundances over 0%
conv_df_zro <- filter(conv_df, Abundance > 0)

# bray curtis dissimilarity for con
con.dist<-vegdist(conv_df_zro$Abundance, method='bray')
con.dist

set.seed(134) #reproducible results

con.div<-adonis2(con.dist~conv_df_zro$Stage, data=conv_df_zro, permutations=9999) # this will take a while
con.div

#PERMANOVA dispersion
con.res.adisper <- betadisper(d = con.dist, group = conv_df_zro$Stage, type = 'centroid')

con.res.adisper$distances

anova(con.res.adisper)
TukeyHSD(con.res.adisper) # post hoc comparisons
# see result sectino for comparison results

# visualize the model
boxplot(con.res.adisper)
plot(con.res.adisper)

# This ends data visualization and analyses 
# Thank you! 