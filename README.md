# CellTag-Classifier

### Description
This repository contains source functions required for CellTag classification and brief tutorial using a species mixing dataset that involves two distinct CellTags. For detailed explanation, please refer to the paper: Guo et. al., CellTag Indexing: a genetic barcode-based multiplexing tool for single-cell technologies, bioRxiv, 2018 (new version under review).

In brief, we introduce a new multiplexing approach using CellTag - a lentiviral-based barcode system. This system allows CellTags to be transribed and captured using single-cell RNA-sequencing, further enable us to use CellTag as the identifiers for each cell. A step in this process is to identify the true signal among different CellTags; thus, further recover the identities of cells in the single-cell transcriptome. Here, we share the source code and demonstrate the usage and flow of this classification approach on a human-mouse species-mixing dataset with two CellTags tags.

### Method and General Tutorial
0. Clone the Github repository
```
git clone 
```

1. Single-cell transcriptome are collected and processed using 10x and CellRanger. CellTag count matrix can be extracted via previous demonstrated workflow - https://github.com/morris-lab/CellTagWorkflow or R package - CloneHunter:https://github.com/morris-lab/CloneHunter (Biddy et. al., Naure, 2018). Here the example count data file is located at data/2-tag species mixing before collapsing.csv.

Load in the matrix
```r

```

2. The CellTag count matrices are log-normalized (log2-based) as the following. We first normalize the total reads number and UMI counts across cells.
![](/equation/Normalization.png)

3. CellTags are then extracted and compared. CellTags with similar sequences are collapsed to the centroid CellTag. For more information, please refer to starcode software - https://github.com/gui11aume/starcode. Briefly, starcode clusters DNA sequences based on the Levenshtein distances between each pair of sequences, from which we collapse similar CellTag sequences to correct for potential errors occurred during single-cell RNA-sequencing process. Default maximum distance from starcode was used to cluster the CellTags.

4. With the collapsed CellTag matrix, we perform another normalization between different CellTags to make them more comparable to each other.

5. Having the coll


