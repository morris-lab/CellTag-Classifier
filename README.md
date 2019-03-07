# CellTag-Classifier

### Description
This repository contains source functions required for CellTag classification and brief tutorial using a species mixing dataset that involves two distinct CellTags. For detailed explanation, please refer to the paper: Guo et. al., CellTag Indexing: a genetic barcode-based multiplexing tool for single-cell technologies, bioRxiv, 2018 (new version under review).

In brief, we introduce a new multiplexing approach using CellTag - a lentiviral-based barcode system. This system allows CellTags to be transribed and captured using single-cell RNA-sequencing, further enable us to use CellTag as the identifiers for each cell. A step in this process is to identify the true signal among different CellTags; thus, further recover the identities of cells in the single-cell transcriptome. Here, we share the source code and demonstrate the usage and flow of this classification approach on a human-mouse species-mixing dataset with two CellTags tags (Mouse Tag - GTTGGCTA, Human Tag - TGCTATAT).

### Method and General Tutorial
#### Before we start
Clone the Github repository
```
git clone https://github.com/morris-lab/CellTag-Classifier.git
cd CellTag-Classifier
```
Open a new R terminal. The source code can be loaded.
```r
setwd("path/to/CellTag-Classifier")
source(src/Classifier_Source.R)
```
#### Initial overview of the data
Single-cell transcriptome are collected and processed using 10x and CellRanger. CellTag count matrix can be extracted via previous demonstrated workflow - https://github.com/morris-lab/CellTagWorkflow or R package - CloneHunter:https://github.com/morris-lab/CloneHunter (Biddy et. al., Naure, 2018). Here the example count data file is located at data/2-tag species mixing before collapsing.csv.
```r
# load celltag counts (before Starcode collapse)
count.celltag <- read.csv("data/2-tag species mixing before collapsing.csv", header = TRUE, row.names=1, stringsAsFactors = F)
# check all CellTags before filtering
all.celltags <- rowSums(count.celltag)
hist(all.celltags, breaks = 50)
```
#### Starcode CellTag sequence analysis and collapse
CellTags are then extracted and compared. CellTags with similar sequences are collapsed to the centroid CellTag. For more information, please refer to starcode software - https://github.com/gui11aume/starcode. Briefly, starcode clusters DNA sequences based on the Levenshtein distances between each pair of sequences, from which we collapse similar CellTag sequences to correct for potential errors occurred during single-cell RNA-sequencing process. Default maximum distance from starcode was used to cluster the CellTags. Here the example starcode output result is located at data/2-tag species mixing starcode result.txt.
```r
# load starcode collapsed result
starcode.out <- read.table("data/2-tag species mixing starcode result.txt", sep = "\t", stringsAsFactors=FALSE)

# create a list to store consensus celltags for each centroid celltag
consensus.lookup <- list()
for (i in 1:nrow(starcode.out)) { # repeat for each row
  centroid <- 0
  # split consensus celltags separated by comma and populate list
  centroid <- strsplit(starcode.out[i,3], ",")
  names(centroid) <- starcode.out[i,1]
  consensus.lookup[[i]] <- split(unname(centroid), names(centroid))
}

# create new matrix where each count.celltag column is renamed to its centroid CellTag
count.celltag.collapsed <- as.data.frame(count.celltag)
for (j in 1:ncol(count.celltag.collapsed)) {
  n <- grep(colnames(count.celltag.collapsed[j]), consensus.lookup)
  colnames(count.celltag.collapsed)[j] <- names(consensus.lookup[[n]])
}

# collapse count.celltag by consensus CellTags
count.celltag.collapsed <- t(count.celltag.collapsed)
count.celltag.collapsed <- by(count.celltag.collapsed, INDICES=row.names(count.celltag.collapsed), FUN=colSums)
count.celltag.collapsed <- as.data.frame(do.call(cbind,count.celltag.collapsed))
```

#### Normalization
The CellTag count matrices are log-normalized (log2-based) as the following. We first normalize the total reads number and UMI counts across cells.
<p align="center">
  <img src="/equation/Normalization.png" height="72" width="450">
</p>

```r
# Log the count matrix
count.norm.expr <- log2(count.celltag.collapsed+1)
ct.dge.norm <- t(count.norm.expr)
rownames(ct.dge.norm) <- colnames(count.celltag.collapsed)
colnames(ct.dge.norm) <- rownames(count.celltag.collapsed)
# Normalize
norm.ct.dge <- normalize.function(ct.dge.norm)
```

We identify the top two most abundant CellTags, which match with what we knew in the experimental design. We extract those two for further analysis.
```r
count.norm.expr.t <- t(norm.ct.dge)
sort(colSums(count.norm.expr.t))

# pull out dge by 2 most abundant celltags
ct.dge <- count.norm.expr.t[,c("TGCTATAT", "GTTGGCTA")]
collapsed.orig.count <- count.celltag.collapsed[,c("TGCTATAT", "GTTGGCTA")]
```

we perform another normalization between different CellTags to make them more comparable to each other.
```r
norm.ct.dge.2 <- normalize.function(ct.dge)
```

#### Perform permutation sampling and dynamic CellTag detection
Assuming balanced loading of cell number from two groups without overloading on the 10x machine, we expect a relatively low multiplet rate, suggesting most cells should have one or the other CellTag. In such case, each CellTag expression across all cells should have ~50% zero and the remaining as significantly expressed. This allows us to assume that CellTags' expression across cells would share similar density functions. Under this assumption, we will examine the significance of each CellTag expression, i.e. what is the likelihood of occurrence of each expression value.

For each CellTag expression, we compute the density function *D* of its expression across all cells. Each cell's expression for each tag is then examined via a modified permutation test. In brief, for expression of a CellTag in a cell to be tested, we draw 1,000 samples from the density functions and calculate the proportion of samples that are greater than or equal to expression level at test. For instance, consider expression (*C<sub>ij</sub>*) of CellTag *j* in Cell *i*, we draw 1,000 sample *S* from the density of CellTag *j*, *D<sub>j</sub>*. The proportion was computed as shown below. This process is iterated for at least 50 times to make sure that the samples are representative of the overall density.
<p align="center">
  <img src="/equation/permutation.png" height="72" width="140">
</p>

This function takes the normalized matrix, sample size and iteration number as inputs. It outputs the result in a list format. For each item in the list, there are two parts with the first part = Cell barcode and second part = data frame that contains proportion information from each CellTag and each iteration (sample) as below.

||CellTag 1|CellTag 2|\<More CellTags\>|CellTag N|
|:-------------:|:-------------:|:-------------:|:-------------:|:-------------:|
|1|*P<sub>11</sub>*|*P<sub>12</sub>*|...|*P<sub>1N</sub>*|
|2|*P<sub>21</sub>*|*P<sub>22</sub>*|...|*P<sub>2N</sub>*|
|...|...|...|...|...|
|50|*P<sub>501</sub>*|*P<sub>502</sub>*|...|*P<sub>50N</sub>*|

```r
perc.ls <- dynamic.celltag.detection(norm.ct.dge.2)
```

