  #! /usr/bin/Rscript --vanilla

  args <- commandArgs(trailingOnly = TRUE)

  options(warn = - 1)

  suppressMessages(library(dplyr))
  suppressMessages(library(tidyr))
  suppressMessages(library(ggplot2))
  suppressMessages(library(forcats))
  suppressMessages(library(readr))

  data = read.table(file = file.path(args[1], "final_genes.tsv"), sep = '\t', header = TRUE)
  genes = read.table(file = file.path(args[1], "geneNumbers.tsv"), sep = '\t',
                     header = FALSE, col.names = c("species_name", "total_genes"))

  # Rename species name not to include accession numbers
  data$species_name <- gsub("-.*", "", data$species_name)

  # Create count summary stats
  StatsData <- data %>%
   group_by(species_name) %>%
   summarise(unique_genes = n())

   # Merge dataframes
   allGenes = merge(StatsData, genes, by="species_name", all = T)

   # Remove any unexpected species name created due to eventual trailing spaces
   allGenes <- allGenes %>%
     filter(!is.na(total_genes))

   # Change all NA values to 0
   allGenes[is.na(allGenes)] <- 0

   # Calculate the number of non-unique genes
   allGenes <- allGenes %>%
     mutate(nonUnique_genes = total_genes - unique_genes)

   # Sort dataframe by decreasing number of unique genes
   allGenes <- allGenes[order(allGenes$unique_genes, decreasing = T),]

   # Reorder columns in dataframe
   allGenes <- allGenes[, c(1,2,4,3)]

   # Pivot longer to create column for barplot labels
   allGenes_label <- allGenes %>%
     subset(select = -total_genes) %>%
     # Reorder species name based on number of unique genes
     mutate(species_name = fct_reorder(species_name, unique_genes)) %>%
     pivot_longer(cols = c(unique_genes,nonUnique_genes),
                  names_to = "type", values_to = "counts")

  # Count the number of files in the CDS folder. This is the number of species
  # ran as input for the pipeline
  CDSnums <- length(list.files(args[2], pattern = "*fna"))

  # Define functions

  # Function to plot number of genes (for smaller databases)
plotting <- function(DB,ySize,gSize,cSize, xText, xTitle) {
  plot <- DB %>%
    ggplot(aes(y=species_name, x=counts, fill=type)) +
    geom_bar(stat="identity") +
    geom_col(colour = "black", size = cSize) +
    geom_label(data = DB %>% filter(type == "unique_genes"), aes(label=counts), colour = "black",
               x = -250, show.legend = F, label.size = 0.15,
               position = position_stack(),
               size = gSize,
               hjust = 0.5) +
    scale_x_continuous(sec.axis = dup_axis(), limits=c(-250, max(DB$counts)+1500)) +
    scale_fill_brewer(palette = "Pastel1",
                      breaks = c("unique_genes", "nonUnique_genes"),
                      labels = c("Unique", "Non-unique")) +
    labs(x="Number of taxon-specific genes") +
    theme_classic() +
    theme(
      axis.title.y = element_blank(),
      axis.ticks.length = unit(0.15, "cm"),
      axis.text.x=element_text(size=xText, margin = margin(r = 4)),
      axis.title.x = element_text(size = 16),
      axis.text.y=element_text(size=ySize, margin = margin(r = 2)),
      legend.justification = "center", legend.position = "bottom",
      legend.title = element_blank(),
      plot.margin = unit(c(0.1,0.5,0.5,0.5), "cm"), # second position is for the right most y-axis, third position is for bottom x-axis
      panel.grid.major.x = element_line(),
      legend.text = element_text(size=13),
      #legend.title = element_text(size=14),
      legend.background = element_blank(),
      legend.box.background = element_rect(colour="black", size = 0.5)
    )

}


# Function to save ggplot
saving_plot <- function(Width, Height) {
  plotname=file.path(args[1], "UniqueGenesPlot.pdf")
  ggsave(plotname, width = Width, height = Height, units = "in")
}

# plot and save based on number of CDS files (i.e. species) to plot on the
# y-axis

if ( CDSnums > 0 & CDSnums < 20 )  {

  plotting(allGenes_label, 14, 4, 0.5, 14, 16)
  saving_plot(10,8)

} else if ( CDSnums >= 20 & CDSnums < 35 ) {

  plotting(allGenes_label, 14, 3, 0.5, 14, 16)
  saving_plot(10,11)

} else if ( CDSnums >= 35 & CDSnums < 50 ) {

  plotting(allGenes_label, 13.5, 2.5, 0.5, 14, 16)
  saving_plot(10,13)


} else if ( CDSnums >= 50 & CDSnums < 60 ) {

  plotting(allGenes_label, 13, 2.5, 0.3, 14, 16)
  saving_plot(10,15)


} else if ( CDSnums >= 60 & CDSnums < 75 ) {

  plotting(allGenes_label, 12, 2, 0.3, 14, 16)
  saving_plot(10,15)

} else if ( CDSnums >=75  & CDSnums < 85 ) {

  plotting(allGenes_label, 11, 1.77, 0.3, 14, 16)
  saving_plot(10,17)

} else if ( CDSnums >= 85 & CDSnums < 95 ) {

  plotting(allGenes_label, 10.5, 1.75, 0.3, 14, 16)
  saving_plot(10,17.5)

} else if ( CDSnums >= 95 & CDSnums < 105 ) {

  plotting(allGenes_label, 10, 1.65, 0.4, 13, 15)
  saving_plot(10,19)

} else if ( CDSnums >= 105 & CDSnums < 125 ) {

  plotting(allGenes_label, 9.5, 1.50, 0.4, 13, 15)
  saving_plot(10,20)

} else if ( CDSnums >=125 & CDSnums < 135 ) {

  plotting(allGenes_label, 9, 1.40, 0.4, 13, 15)
  saving_plot(10,22)

} else if ( CDSnums >=135 & CDSnums <= 150 ) {

  plotting(allGenes_label, 9, 1.40, 0.3, 13, 15)
  saving_plot(10,24)

} else {

  plotting(allGenes_label, 8, 1.5, 0.25, 13, 15)
  saving_plot(10,30)

}

# Save dataframe as tsv for stats on pipeline output
filename=file.path(args[1], "Stats_pipelineOutput.tsv")
write_tsv(allGenes, file=filename, col_names = T)