
#
##
###
####
##### Loading necessary packages from your global enviroment 
####
###
##
#



library("phyloseq")
library("metagenomeSeq")
library("plyr")
library("dplyr")
library("ggplot2")
library("data.table")
library("tidyr")
library("forcats")
library("vegan")
library("pheatmap")
library("dendsort")
library("ggrepel")
library("devtools")
library("pairwiseAdonis")
library("DAtest")


#
##
###
####
##### Loading AMR files 
####
###
##
#



# Load in the metadata file
sample_metadata <- read.table('./AMR_metadata_portal.csv', header=T, row.names=1, sep=',')

# MEGARes has a separate annotation file that contains this information as well and facilitates the aggregation of counts.
# For example, we might want to sum all counts to gene accessions that confer resistance to drugs in the Tetracycline class.
annotations <- read.table('./megares_full_annotations_v2.00.csv', header=T, row.names=1, sep=',')

# Load resistome MEGARes count matrix                                       
amr <- read.table('./Small_Results_WithKraken0.1/ResistomeResults/AMR_analytic_matrix.csv', header=T, row.names=1, sep=',')

# Load the kraken2 count matrix                                        
kraken_microbiome <- read.table('./Small_Results_WithKraken0.1/KrakenResults/kraken_analytic_matrix.csv', header=T, row.names=1, sep=',')


#
##
###
####
##### Exploring metadata file 
####
###
##
#



# Let's view the metadata file we loaded into R
View(sample_metadata)

# View structure of the data
str(sample_metadata)

# We also know the names for the columns and they can be accessed like this:
sample_metadata$Sample.ID

# We can now check the values of our Sample.ID column
length(sample_metadata$Sample.ID)

# Now, we want to know many host species we have in the column, "Host"
# Just using length() gives us a total count
length(sample_metadata$Host)

# But, we can specify that we only want the unique variables
unique(sample_metadata$Host)
length(unique(sample_metadata$Host))



#
##
###
####
##### Exploring annotation file for MEGARes gene accessions
####
###
##
#



# Explore the first 25 lines of the annotation object 
head(annotations, n= 25L)

# Get all the row names for the "amr" object and place it in a new object called "amr_headers"
amr_headers <- row.names(amr)

# Try extracting all of the data from annotations that matches "amr_headers"
subset_annotations <- annotations[amr_headers,]

# Now, we we can use the "subset_annotations" object and see how many unique AMR classes and mechanisms we found
unique(subset_annotations$class)

# Notice that the output is numbered in brackets, but it seems different than the amount of "Levels" reported
unique(subset_annotations$mechanism)

# Notice that the output is numbered in brackets, but it seems different than the amount of "Levels" reported
unique(subset_annotations$type)

# Notice that the output is numbered in brackets, but it seems different than the amount of "Levels" reported
unique(subset_annotations$group)

# We can then use "length" to count the number of unique "class" values (add the column different columns to see their unique lengths)
length(unique(subset_annotations$class))



#
##
###
####
##### Creating the MEGARes resistome (shotgun reads) phyloseq object 
####
###
##
#



# We can convert our amr count object to the otu_table format required for phyloseq
amr <- otu_table(amr, taxa_are_rows = TRUE)

# We can now merge these objects to make a phyloseq object
amr.ps <- merge_phyloseq(amr, tax_table(as.matrix(annotations)), sample_data(sample_metadata))

# Estimating richness and diversity using the easy-to-use function estimate_richness()
amr_shotgun_diversity_values <- estimate_richness(amr.ps)

# We can now use the plot_bar functions to explore the resistome counts
plot_bar(amr.ps)




#
##
###
####
##### Creating the kraken2 microbiome (shotgun reads) phyloseq object
####
###
##
#



# Convert to format that phyloseq likes with otu_table()                                      
kraken_microbiome <- otu_table(kraken_microbiome, taxa_are_rows = TRUE)

# This part is a little more complicated, but you should be able to adopt this code for other similar uses.
# Notice that we are specifying the function after the comma, so in the columns.
# The main function we are using is "tstrsplit()". If you look at the documentation, you'll see that 
# "This is a convenient wrapper function to split a column using strsplit and assign the transposed result to individual columns"
# Here, we specify that we are splitting the "id" column by the "|" character
# To the left of the ":=" we specify the column names for the values split from the id column
kraken_taxonomy <- data.table(id=rownames(kraken_microbiome))
kraken_taxonomy[, c('domain',
                    'kingdom',
                    'phylum',
                    'class',
                    'order',
                    'family',
                    'genus',
                    'species') := tstrsplit(id, '|', type.convert = TRUE, fixed = TRUE)]

# Convert to data.frame
kraken_taxonomy <- as.data.frame(kraken_taxonomy)
# Use the id variable to rename the row.names
row.names(kraken_taxonomy) <- kraken_taxonomy$id
# Remove the "id" column
kraken_taxonomy <- within(kraken_taxonomy, rm(id))

# Create kraken phyloseq object
kraken_microbiome.ps <- merge_phyloseq(kraken_microbiome, tax_table(as.matrix(kraken_taxonomy)), sample_data(sample_metadata))

# Estimating richness and diversity using the easy-to-use function estimate_richness()
microbiome_shotgun_diversity_values <- estimate_richness(kraken_microbiome.ps)

# We can now use the plot_bar functions to explore the microbiome counts
plot_bar(kraken_microbiome.ps)

# Explore the microbiome phyloseq object
kraken_microbiome.ps
# We can access the individual components like this:
# Look at this site for a few more examples: https://joey711.github.io/phyloseq/import-data.html
otu_table(kraken_microbiome.ps)
tax_table(kraken_microbiome.ps)
sample_data(kraken_microbiome.ps)

# Here are some other functions that provide information about the phyloseq object
# The sample names
sample_names(kraken_microbiome.ps)
# Taxonomic rank names
rank_names(kraken_microbiome.ps)
# Variables in the metadata
sample_variables(kraken_microbiome.ps)
# Easy sums of counts by sample
sample_sums(kraken_microbiome.ps)



#
##
###
####
##### Agglomerate ASV counts to different taxonomic levels
####
###
##
#



# Using tax_glom(), we can easily aggregate counts to different levels (taxonomic or AMR annotation levels)
phylum.ps <- tax_glom(kraken_microbiome.ps, "phylum")
phylum.ps

# We can get sample counts at the phylum level.
sum(sample_sums(phylum.ps))

# By summing counts at each taxonomic level, we can observe the reduction in sample counts in lower taxonomic levels.
# Note the stop sign on the right hand side of the Console panel, when you see this it means R is running a command.
# You can click on the stop sign to stop the run.
species.ps <- tax_glom(kraken_microbiome.ps, "species") # Depending on your computer and file size, this can take a couple of seconds to a few minutes

# We can get sample counts at the species level.
sum(sample_sums(species.ps))

# We can calculate the percentage of counts at the species level out of classified at the phylum level.
sum(sample_sums(species.ps))  /  sum(sample_sums(phylum.ps)) * 100



#
##
###
####
##### Plots! plots! plots!
####
###
##


#####################
#  Sequencing depth #
#####################

# Make a column "Sample" to join the metadata file with the diversity indices
sample_metadata$Sample <- row.names(sample_data(sample_metadata))

# Plot sequencing depth 
ggplot(sample_metadata, aes(x = Sample , y = Raw.Paired.Reads)) + 
  geom_bar(stat="identity", position = position_dodge()) +
  labs(title = "Metagenomic sequencing depth", x = "Sample", y = "Raw Paired Reads") + 
  theme_classic() +
  theme(axis.text.x = element_text(size = 10, angle = 90, vjust = 0.5, hjust=1), axis.text.y = element_text(size = 10))

####################################
#  Microbiome diversity AD indices #
####################################

# Alpha-diversity indices are a common measurement to summarize the composition of the microbiome and resistome.
# "Richness" or "Observed" simply describes the number of unique taxa identified in a sample
# On the other hand, we use "Shannon's index" or "Inverse Simpsons's index" to describe "evenness", or
# distribution of counts among taxa in a sample

# Find more information about diversity indices here: http://www.jmb.or.kr/submission/Journal/027/JMB027-12-02_FDOC_2.pdf

# Estimating richness and diversity using the easy-to-use function estimate_richness()
microbiome_diversity_values <- estimate_richness(kraken_microbiome.ps)

# We can do the same for the phylum counts
estimate_richness(phylum.ps) 

# We can easily plot these values using the plot_richness() function, we'll just pick 3 commonly used alpha diversity indices
plot_richness(phylum.ps, x = "Host", color = "Host", measures = c("Observed", "Shannon", "InvSimpson")) 

# We can then modify this figure using ggplot2 functions like geom_boxplot() by using the "+" sign 
plot_richness(phylum.ps, x = "Host", color = "Host", measures = c("Observed", "Shannon", "InvSimpson")) + 
  geom_boxplot()

##########################
#  Microbiome diversity  #
##########################

# Estimate richness and diversity 
microbiome_diversity_values <- estimate_richness(kraken_microbiome.ps)
microbiome_diversity_values$Sample <- row.names(microbiome_diversity_values)

# With phyloseq objects have to convert to matrix first, then to dataframe
sample_metadata <- as.data.frame(as(sample_data(kraken_microbiome.ps), "matrix"))

# Make a column "Sample" to join the metadata file with the diversity indices
sample_metadata$Sample <- row.names(sample_data(sample_metadata))

# Now, we use left_join() to add the sample_metadata to the combined_diversity_values object
microbiome_diversity_values <- left_join(microbiome_diversity_values,sample_metadata, by = "Sample")

ggplot(microbiome_diversity_values, aes(x = Host, y = Observed, color = Host)) +
  geom_boxplot() +
  geom_point() +
  labs(title = "Unique features by Host", x = "Host", y = "Observed features") + 
  theme_classic()

# Same analysis by phylum (raw data)
phylum_microbiome.ps <- tax_glom(kraken_microbiome.ps, "phylum")
phylum_diversity_values <- estimate_richness(phylum_microbiome.ps)
phylum_diversity_values$Sample <- row.names(phylum_diversity_values)

phylum_diversity_values <- left_join(phylum_diversity_values,sample_metadata, by = "Sample")

ggplot(phylum_diversity_values, aes(x = Host, y = Observed, color = Host)) +
  geom_boxplot() +
  geom_point() +
  labs(title = "Unique features by Host at the phylum level", x = "Host", y = "Observed features") + 
  theme_classic()

##########################
#  Resistome diversity  #
##########################

# Estimate richness and diversity 
resistome_diversity_values <- estimate_richness(amr.ps)
resistome_diversity_values$Sample <- row.names(resistome_diversity_values)

# With phyloseq objects have to convert to matrix first, then to dataframe
sample_metadata <- as.data.frame(as(sample_data(amr.ps), "matrix"))

# Make a column "Sample" to join the metadata file with the diversity indices
sample_metadata$Sample <- row.names(sample_data(sample_metadata))

# Now, we use left_join() to add the sample_metadata to the combined_diversity_values object
resistome_diversity_values <- left_join(resistome_diversity_values,sample_metadata, by = "Sample")

ggplot(resistome_diversity_values, aes(x = Host, y = Observed, color = Host)) +
  geom_boxplot() +
  geom_point() +
  labs(title = "Unique features by Host", x = "Host", y = "Observed features") + 
  theme_classic()

##########################
# Microbiome composition #
##########################

# First let's filter and normalize our counts with CSS normalization and aggregate counts at the phylum level
filtered_microbiome.ps = filter_taxa(kraken_microbiome.ps, function(x) sum(x) > 5, TRUE)
filtered_microbiome.metaseq <- phyloseq_to_metagenomeSeq(filtered_microbiome.ps)
cumNorm(filtered_microbiome.metaseq)
CSS_microbiome_counts <- MRcounts(filtered_microbiome.metaseq, norm = TRUE)
CSS_normalized_microbiome.ps <- merge_phyloseq(otu_table(CSS_microbiome_counts, taxa_are_rows = TRUE),sample_data(kraken_microbiome.ps),tax_table(kraken_microbiome.ps))
CSS_normalized_phylum_microbiome.ps <- tax_glom(CSS_normalized_microbiome.ps, "phylum")

# Plot the microbiome compostion at the phylum level
plot_bar(CSS_normalized_phylum_microbiome.ps, fill = "phylum") + 
  facet_wrap(~ Host, scales = "free_x") +
  theme_classic()

# Convert OTU abundances to relative abundances
rel_phylum_microbiome.ps <- transform_sample_counts(CSS_normalized_phylum_microbiome.ps, function(x) x / sum(x) )

# We can plot these results, notice the y-axis relative abundance plots
plot_bar(rel_phylum_microbiome.ps, fill= "phylum")

# only keep OTUs present in greater than 0.5% of all OTUs across all samples.
minTotRelAbun = 0.005
x = taxa_sums(CSS_normalized_phylum_microbiome.ps)
keepTaxa = (x / sum(x)) > minTotRelAbun
length(keepTaxa[keepTaxa ==TRUE])
pruned_phylum_microbiome.ps = prune_taxa(keepTaxa, CSS_normalized_phylum_microbiome.ps)
plot_bar(pruned_phylum_microbiome.ps, fill = "phylum")

rel_pruned_phylum_microbiome.ps <- transform_sample_counts(pruned_phylum_microbiome.ps, function(x) x / sum(x) )

plot_bar(rel_pruned_phylum_microbiome.ps, fill = "phylum")

# We can calculate distance measures to cluster samples
mat_cluster_cols <- hclust(dist(t(otu_table(CSS_normalized_phylum_microbiome.ps))))

# We can also change the clustering of our samples using the function below
sort_hclust <- function(...) as.hclust(dendsort(as.dendrogram(...)))
sorted_mat_cluster_cols <- sort_hclust(mat_cluster_cols)

# Now let's cluster the microbiome features
mat_cluster_rows <- sort_hclust(hclust(dist(otu_table(CSS_normalized_phylum_microbiome.ps))))

# The features in the tax_table have unique IDs for each taxa, but those names
# don't get updated to the aggregated taxa labels. Instead, we'll pull out that information
# and use if for plotting below
phylum_taxa_labels <- as.data.frame(as(tax_table(CSS_normalized_phylum_microbiome.ps),"matrix"))
phylum_taxa_labels$phylum

# Change the taxa names for your phyloseq object
taxa_names(CSS_normalized_phylum_microbiome.ps) <- phylum_taxa_labels$phylum

# Here's an example of plotting a heatmap with the "pheatmap" package
# Notice we added a pseudocount, prior to log normalization.
# Without the pseudocount, CSS values < 1 are converted to negative values and cause an error
# Below are just a few examples of flags we can use, check out the pheatmap() documentation
# for more options.
pheatmap( log(otu_table(CSS_normalized_phylum_microbiome.ps) + 1), cluster_cols = sorted_mat_cluster_cols,
          cluster_rows = mat_cluster_rows, 
          drop_levels = TRUE , fontsize = 10, treeheight_row = 10)

##########################
# Resistome composition #
##########################

# Now for the resistome let's filter and normalize our counts with CSS normalization and aggregate counts at the class level
filtered_resistome.ps = filter_taxa(amr.ps, function(x) sum(x) > 5, TRUE)
filtered_resistome.metaseq <- phyloseq_to_metagenomeSeq(filtered_resistome.ps)
cumNorm(filtered_resistome.metaseq)
CSS_resistome_counts <- MRcounts(filtered_resistome.metaseq, norm = TRUE)
CSS_normalized_resistome.ps <- merge_phyloseq(otu_table(CSS_resistome_counts, taxa_are_rows = TRUE),sample_data(amr.ps),tax_table(amr.ps))
CSS_normalized_class_resistome.ps <- tax_glom(CSS_normalized_resistome.ps, "class")

# Plot the resitome composition at the class level
plot_bar(CSS_normalized_class_resistome.ps, fill = "class") + 
  facet_wrap(~ Host, scales = "free_x") +
  theme_classic() +
  theme(axis.text.x = element_text(size = 7), axis.text.y = element_text(size = 10))

# Convert OTU abundances to relative abundances
rel_phylum_resistome.ps <- transform_sample_counts(CSS_normalized_class_resistome.ps, function(x) x / sum(x) )

# We can plot these results, notice the y-axis relative abundance plots
plot_bar(rel_phylum_resistome.ps, fill= "class")

# only keep OTUs present in greater than 0.5% of all OTUs across all samples.
minTotRelAbun = 0.005
x = taxa_sums(CSS_normalized_class_resistome.ps)
keepTaxa = (x / sum(x)) > minTotRelAbun
length(keepTaxa[keepTaxa ==TRUE])
pruned_class_resistome.ps = prune_taxa(keepTaxa, CSS_normalized_class_resistome.ps)
plot_bar(pruned_class_resistome.ps, fill = "class")

rel_pruned_class_resistome.ps <- transform_sample_counts(pruned_class_resistome.ps, function(x) x / sum(x) )

plot_bar(rel_pruned_class_resistome.ps, fill = "class")

# We can calculate distance measures to cluster samples
mat_cluster_cols <- hclust(dist(t(otu_table(CSS_normalized_class_resistome.ps))))

# We can also change the clustering of our samples using the function belowz
sort_hclust <- function(...) as.hclust(dendsort(as.dendrogram(...)))
sorted_mat_cluster_cols <- sort_hclust(mat_cluster_cols)

# Now let's cluster the microbiome features
mat_cluster_rows <- sort_hclust(hclust(dist(otu_table(CSS_normalized_class_resistome.ps))))

class_taxa_labels <- as.data.frame(as(tax_table(CSS_normalized_class_resistome.ps),"matrix"))
class_taxa_labels$class

# Change the taxa names for your phyloseq object
taxa_names(CSS_normalized_class_resistome.ps) <- class_taxa_labels$class

# Plot a heatmap with a pheatmap like the microbiome example above
pheatmap( log(otu_table(CSS_normalized_class_resistome.ps) + 1), cluster_cols = sorted_mat_cluster_cols,
          cluster_rows = mat_cluster_rows, 
          drop_levels = TRUE , fontsize = 10, treeheight_row = 10)

###########################
#  Microbiome ordination  #
###########################

ordination_phylum_bray <- ordinate(CSS_normalized_phylum_microbiome.ps, method = "NMDS" , distance="bray")

# We can use "plot_ordination()" to plot the distance matrix
# We specify that we want to compare "samples" and color the points by the "Host" metadata variable
plot_ordination(CSS_normalized_phylum_microbiome.ps, ordination_phylum_bray, type = "samples",color = "Host")
group_variable = get_variable(CSS_normalized_phylum_microbiome.ps,"Host")
anosim(distance(CSS_normalized_phylum_microbiome.ps, "bray"), group_variable)

###########################
#  Resistome ordination  #
###########################

ordination_class_bray <- ordinate(CSS_normalized_class_resistome.ps, method = "NMDS" , distance="bray")

# We can use "plot_ordination()" to plot the distance matrix
# We specify that we want to compare "samples" and color the points by the "Host" metadata variable
plot_ordination(CSS_normalized_class_resistome.ps, ordination_class_bray, type = "samples",color = "Host")
group_variable = get_variable(CSS_normalized_class_resistome.ps,"Host")
anosim(distance(CSS_normalized_class_resistome.ps, "bray"), group_variable)

# Look at your Global Environment panel (usually top right)
# Unless you are saving the Workspace with R objects you created previously, this panel should be empty.
# If it isn't empty, we can remove all of those R objects 

# Let's check for any variables in your workspace
ls() # If you don't have any variables created, this will show "character(0)"

# To delete them all, we have to use the 
rm(list=ls())  # If your environment is empty already, this won't show anything


