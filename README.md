## Formatting a custom 16S sequence databases for use in DADA2
This script uses R to modify a custom-curated sequence database (DictDB) into a training fasta file appropriate for assignTaxnomy in DADA2

I will be using a termite-specific curated database called [DictDB](https://www.sciencedirect.com/science/article/pii/S0723202015001162) to classify and assign taxonomy for my termite 16S reads. This database improves taxonomy assignments for termite gut 16S rRNA gene sequences compared to the standard training sets from SILVA or RDP.

Download DictDB (version 3.0) [here](https://drive.google.com/file/d/0B5NJB7U2TQKNVEZneUhFODZFdmc/view).

Unzip and move the 2 files (fasta and taxonomy file) into a new folder named "DictDB" and navigate to that directory.

#### View and inspect the DictDB fasta and taxonomy files:
```{bash}
head -6 DictDB_V3.fasta

>UltBac92
AGAUUUUGCUUUUGGCUCAGAUUGAACGCUGGCGGCAGGCCUAACACAUGCUUGUCGAACGGUAACAGGAAACAG...
>UltBac89
AGAGUUUCUUCCUGGCUCAGAUUGAACGCUGGCGGCAGGCCUAACACAUGCAAGUCGAACGGUAACAGGAAGCAG...
>UltBac90
UGAACGCUGGCGGCAGGCCUAACACAUGCUUUGUCGAACGGUAACAGGAAGAAGCUUGCUUCUUUGCUGACGAGU...

head -6 DictDB_V3.taxonomy

UltBac92 Proteobacteria;Gammaproteobacteria_1;Enterobacteriales;Enterobacteriaceae;Escherichia-Shigella;
UltBac89 Proteobacteria;Gammaproteobacteria_1;Enterobacteriales;Enterobacteriaceae;Escherichia-Shigella;
UltBac90 Proteobacteria;Gammaproteobacteria_1;Enterobacteriales;Enterobacteriaceae;Escherichia-Shigella;
UltBac42 Proteobacteria;Gammaproteobacteria_1;Enterobacteriales;Enterobacteriaceae;Escherichia-Shigella;
UltBac44 Proteobacteria;Gammaproteobacteria_1;Enterobacteriales;Enterobacteriaceae;Escherichia-Shigella;
UltBac55 Proteobacteria;Gammaproteobacteria_1;Enterobacteriales;Enterobacteriaceae;Escherichia-Shigella;
```


The fasta file contains the taxa id followed by the sequence after a line. The taxonomy file contains the same taxa id with its associated taxonomy string all on the same line, in the same order as the taxa ids provided in the fasta file.

Thus, the fasta file should have exactly twice the number of lines as the taxonomy file:
```{bash}
wc -l DictDB_V3.fasta
wc -l DictDB_V3.taxonomy

111464 DictDB_V3.fasta
55732 DictDB_V3.taxonomy
```


This checks out.

#### To make sure the taxa id strings are matching in both files, add a ">" to the front of every taxa id in the taxonomy file:
```{bash}
sed -e 's/^/>/' DictDB_V3.taxonomy > DictDB_V3_modified.taxonomy

head -6 DictDB_V3_modified.taxonomy

>UltBac92    Proteobacteria;Gammaproteobacteria_1;Enterobacteriales;Enterobacteriaceae;Escherichia-Shigella;
>UltBac89    Proteobacteria;Gammaproteobacteria_1;Enterobacteriales;Enterobacteriaceae;Escherichia-Shigella;
>UltBac90    Proteobacteria;Gammaproteobacteria_1;Enterobacteriales;Enterobacteriaceae;Escherichia-Shigella;
>UltBac42    Proteobacteria;Gammaproteobacteria_1;Enterobacteriales;Enterobacteriaceae;Escherichia-Shigella;
>UltBac44    Proteobacteria;Gammaproteobacteria_1;Enterobacteriales;Enterobacteriaceae;Escherichia-Shigella;
>UltBac55    Proteobacteria;Gammaproteobacteria_1;Enterobacteriales;Enterobacteriaceae;Escherichia-Shigella;
```

With the ">" in front of both taxa id strings, this will allow us to to do string matching. It also preserves the fasta format for later on.

#### Load libraries in R, set directory, and read in fasta and taxonomy files as dataframes:
```{r}
library(stringr)
library(dplyr)

setwd("/Users/jonathanylin/Desktop/DictDB")

fasta <- read.table("DictDB_V3.fasta")
taxonomy <- read.table("DictDB_V3_modified.taxonomy")
```

#### View and inspect loaded files:
```{r}
head(fasta)

          V1
1 >UltBac92
2 AUCCCUAGCUGGUCUGAGAGGAUGACCAGCCACACUGGAA...
3 >UltBac89
4 AGAGUUUCUUCCUGGCUCAGAUUGAACGCUGGCGGCAGGC...
5 >UltBac90
6 UGAACGCUGGCGGCAGGCCUAACACAUGCUUUGUCGAACG...

head(taxonomy)

          V1
1 >UltBac92
2 >UltBac89
3 >UltBac90
4 >UltBac42
5 >UltBac44
6 >UltBac55
                                                                                               V2
1 Proteobacteria;Gammaproteobacteria_1;Enterobacteriales;Enterobacteriaceae;Escherichia-Shigella;
2 Proteobacteria;Gammaproteobacteria_1;Enterobacteriales;Enterobacteriaceae;Escherichia-Shigella;
3 Proteobacteria;Gammaproteobacteria_1;Enterobacteriales;Enterobacteriaceae;Escherichia-Shigella;
4 Proteobacteria;Gammaproteobacteria_1;Enterobacteriales;Enterobacteriaceae;Escherichia-Shigella;
5 Proteobacteria;Gammaproteobacteria_1;Enterobacteriales;Enterobacteriaceae;Escherichia-Shigella;
6 Proteobacteria;Gammaproteobacteria_1;Enterobacteriales;Enterobacteriaceae;Escherichia-Shigella;
```

Note that in the dataframe of the fasta file, the taxa id and accompanying sequence are in the same column, but separated in different rows. In the taxonomy dataframe, the taxa id and matching taxonomy string are on the same rows across two different columns.

The assignTaxonomy() command in DADA2 expects a fasta training file in the following format:

```
>Kingdom;Phylum;Class;Family;Genus;Species
sequence
>Kingdom;Phylum;Class;Fmaily;Genus;
sequence
```
See [here](https://benjjneb.github.io/dada2/training.html) for more information.
We'll need to search through the taxonomy and fasta dataframes and pull out the matching taxa ids with the taxonomy strings.

#### Create a new dataframe with the same number of rows and columns as the taxonomy file:
```{r}
combined <- select(taxonomy, V1)
combined$V1 <- NA
combined$V2 <- NA
```

#### Use this loop:
```{r}
for (i in taxonomy$V1) {
  for (ii in fasta$V1) {
    if (i == ii) {
      combined$V1[which(taxonomy$V1 == i)] = paste(">", taxonomy$V2[which(taxonomy$V1 == i)], sep = "")
      R = which(fasta$V1 == ii) + 1
      combined$V2[which(taxonomy$V1 == i)] = paste(fasta$V1[R], sep = "")
    }
  }
}
```
This:
- Iterates through the taxonomy and fasta files and pastes the matching taxa id, with a ">", into the new dataframe "combined"
- Pastes the matching sequence for each taxa string to a second column (same row) in the dataframe "combined"
- Takes ~ 1 hour to iterate through 55,000+ lines

#### View combined dataframe:
```{r}
head(combined)

                                                                                                V1
1 >Proteobacteria;Gammaproteobacteria_1;Enterobacteriales;Enterobacteriaceae;Escherichia-Shigella;
2 >Proteobacteria;Gammaproteobacteria_1;Enterobacteriales;Enterobacteriaceae;Escherichia-Shigella;
3 >Proteobacteria;Gammaproteobacteria_1;Enterobacteriales;Enterobacteriaceae;Escherichia-Shigella;
4 >Proteobacteria;Gammaproteobacteria_1;Enterobacteriales;Enterobacteriaceae;Escherichia-Shigella;
5 >Proteobacteria;Gammaproteobacteria_1;Enterobacteriales;Enterobacteriaceae;Escherichia-Shigella;
6 >Proteobacteria;Gammaproteobacteria_1;Enterobacteriales;Enterobacteriaceae;Escherichia-Shigella;
                                                                                                V2
1 AGAUUUUGCUUUUGGCUCAGAUUGAACGCUGGCGGCAGGCCUAACACAUGCUUGUCGAACGGUAACAGGAAACAGCUUGCUGUUUCGCUGACGAGU...
2 AGAGUUUCUUCCUGGCUCAGAUUGAACGCUGGCGGCAGGCCUAACACAUGCAAGUCGAACGGUAACAGGAAGCAGCUUGCUGCUUCGCUGACGAGU...
3 UGAACGCUGGCGGCAGGCCUAACACAUGCUUUGUCGAACGGUAACAGGAAGAAGCUUGCUUCUUUGCUGACGAGUGGCGGACGGGUGAGUAAUGUC...
4 UUGAACGCUGGCGGCAGGCCUAACACAUGUUUGCCGAACGGGAACAGGAAGAAGCUUGCUUCUUUGCUGACGAGUGGCGGGCGGGUGAGUAAUGUC...
5 AUUGAACGCUGGCGGCAGGCCUAACACAUGCAAGUCGAACGGUAACAGGAAGAAGCUUGCUUCUUUGCUGACGAGUGGCGGACGGGUGAGUAAUGU...
6 AGAGUUUGAUCCUGGCUCAGAUUGAACGCUGGCGGCAGGCCUAACACAUGCAAGUCGAACGGUAACAGGAAGAAGCUUGCUUCUUUGCUGACGAGU...

```

The taxonomy strings, now each with a ">", should be in a column. The matching sequence for each taxonomy string should also be in the next column.

Now we need to combine these two columns into a fasta-ready format, where the taxonomy string is followed with the sequence:
The function below was taken from [here](https://stackoverflow.com/questions/23374100/convert-table-into-fasta-in-r).
```{r}
combined_fasta <- do.call(rbind, lapply(seq(nrow(combined)), function(i) t(combined[i, ])))

head(combined_fasta)

    1                    
V1 ">Proteobacteria;Gammaproteobacteria_1;Enterobacteriales;Enterobacteriaceae;Escherichia-Shigella;"    
V2 "AGAUUUUGCUUUUGGCUCAGAUUGAACGCUGGCGGCAGGCCUAACACAUGCUUGUCGAACGGUAACAGGAAACAGCUUGCUGUUUCGCUGACGAGUG..."
V1 ">Proteobacteria;Gammaproteobacteria_1;Enterobacteriales;Enterobacteriaceae;Escherichia-Shigella;" 
V2 "AGAGUUUCUUCCUGGCUCAGAUUGAACGCUGGCGGCAGGCCUAACACAUGCAAGUCGAACGGUAACAGGAAGCAGCUUGCUGCUUCGCUGACGAGUG..."
V1 ">Proteobacteria;Gammaproteobacteria_1;Enterobacteriales;Enterobacteriaceae;Escherichia-Shigella;"
V2 "UGAACGCUGGCGGCAGGCCUAACACAUGCUUUGUCGAACGGUAACAGGAAGAAGCUUGCUUCUUUGCUGACGAGUGGCGGACGGGUGAGUAAUGUCU..."
```

You can see that the taxa strings and matching sequences are now all in the same columns. Now we can export the reformatted dataframe as a fasta file:
```{r}
write.table(combined_fasta, "DictDB_DADA2_FINAL.fasta", row.names = FALSE, col.names = FALSE, quote = FALSE)
```

#### Inspect the final product:
```
head -6 DictDB_DADA2_FINAL.fasta

>Proteobacteria;Gammaproteobacteria_1;Enterobacteriales;Enterobacteriaceae;Escherichia-Shigella;
AGAUUUUGCUUUUGGCUCAGAUUGAACGCUGGCGGCAGGCCUAACACAUGCUUGUCGAACGGUAACAGGAAACAGCUUGCUGUUUCGCUGACGAGUG...
>Proteobacteria;Gammaproteobacteria_1;Enterobacteriales;Enterobacteriaceae;Escherichia-Shigella;
AGAGUUUCUUCCUGGCUCAGAUUGAACGCUGGCGGCAGGCCUAACACAUGCAAGUCGAACGGUAACAGGAAGCAGCUUGCUGCUUCGCUGACGAGUG...
>Proteobacteria;Gammaproteobacteria_1;Enterobacteriales;Enterobacteriaceae;Escherichia-Shigella;
UGAACGCUGGCGGCAGGCCUAACACAUGCUUUGUCGAACGGUAACAGGAAGAAGCUUGCUUCUUUGCUGACGAGUGGCGGACGGGUGAGUAAUGUCU...
```

This file is now ready to use as a custom training set for the assignTaxonomy() step in DADA2!
