#### this script uses R to modify a curated-curated sequence database (DictDB) into a training fasta file appropriate for assignTaxnomy() in DADA2
#### written by Jonathan Lin @ UC Davis

library(stringr)
library(dplyr)

#set working directory
setwd("/Users/jonathanylin/Desktop/DictDB")

#load fasta file
fasta <- read.table("DictDB_V3.fasta")

#load taxonomy file. Make sure that a ">" is added before each taxa identifier before laoding
taxonomy <- read.table("DictDB_V3_modified.taxonomy")

#create a new dataframe containing the same total number of lines (OTUs) as taxonomy file
combined <- select(taxonomy, V1)
combined$V1 <- NA
combined$V2 <- NA

#iterates through taxonomy and fasta files; pastes matching taxa string with ">" in new dataframe
#pastes matching sequence for each taxa string to second column (same row) in new dataframe
#takes ~1 hour to search through 55,000+ lines
for (i in taxonomy$V1) {
  for (ii in fasta$V1) {
    if (i == ii) {
      combined$V1[which(taxonomy$V1 == i)] = paste(">Bacteria;", taxonomy$V2[which(taxonomy$V1 == i)], sep = "")
      R = which(fasta$V1 == ii) + 1
      combined$V2[which(taxonomy$V1 == i)] = paste(fasta$V1[R], sep ="")
    }
  }
}

#reformat combined dataframe such that sequences follow taxa strings
#code taken from:https://stackoverflow.com/questions/23374100/convert-table-into-fasta-in-r
combined_fasta <- do.call(rbind, lapply(seq(nrow(combined)), function(i) t(combined[i, ])))

#export converted dataframe as fasta file
write.table(combined_fasta, "DictDB_DADA2_FINAL.fasta", row.names = FALSE, col.names = FALSE, quote = FALSE)

#file is now ready to use as custom training set for the assignTaxonomy() step in DADA2