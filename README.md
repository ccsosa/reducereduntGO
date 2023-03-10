## Description:

A full list of libraries needed is included below.

**Dependencies:** `R (>= 4.0.0)`

**Imports:** `base, utils,DBI, GO.db`


## Usage

The reducereduntGO.R function uses the results of enrichment analysis and remove redundant terms for biological process. Please load the function before try this this example. This function does:

- Get the children terms of each GO term using the function GOBPCHILDREN and stack them in a data frame
- Check if the children terms are in the enrichment output (no filtered by FDR) and add the FDR values and add the parent term as well to get a table with both, the children and parents
- Use the min with which.min the GO term id with the lowest FDR value to only conserve children terms with FDR values lower than the parent term
- Save the GO term id in a slot of a list with the same rows as the enrichment output
- Filter by the unique list of GO terms obtained.

## Usage
- Pending tasks: implement function to use molecular function and cellular compartment (DONE) 
LOAD reducereduntGOALL.R and add the mode. It can be "BP", "MF" or "CC". Don´t combine them!
``` r
x <- reducereduntGOALL(df = x,GO.ID="term_id",FDR_col="p_value",mode="MF")
```
- An example using GOCompare is provided in the R file Cancer_hallmark_reduce_terms.R

## Requirements
- To use this function you need a dataframe with the enrichment results. I recommend only use with FDR values thrseshold (e.g. FDR <0.05) to avoid excess of computational time. The functions parameters are:
- `df` the results dataframe which includes p values and a column feature where the categories (genelists are named hereif it is the case). The df looks like this:

feature | go_term_name | term_id | p_value
------------ | ------------- | ------------- | -------------
AID | Response to stress | GO:0006950  | 0.01
DCE | Response to external biotic stimulus | GO:0043207 | 0.04
AID | Regulation of cell size | GO:0008361  | 0.005

- `GO.ID`is the column name of the GO term id (e.g. GO:00001)
- `FDR_col.ID`is the column name with the FDR values

```r
###Example
require(gprofiler2);require(stringr);require(GOCompare)

url_file = "https://raw.githubusercontent.com/ccsosa/R_Examples/master/Hallmarks_of_Cancer_AT.csv"
x <- read.csv(url_file)
x[,1] <- NULL
CH <- c("AID","AIM","DCE","ERI","EGS","GIM","IA","RCD","SPS","TPI")


x_Hsap <- lapply(seq_len(length(CH)), function(i){
  x_unique <- unique(na.omit(x[,i]))
  x_unique <- x_unique[which(x_unique!="")]
  x_unique <- as.list(x_unique)
  return(x_unique)
})

names(x_Hsap) <- CH

#Using as background the unique genes for the ten CH.
GOterm_field <- "term_name"
x_s <-  gprofiler2::gost(query = x_Hsap,
                         organism = "hsapiens", ordered_query = FALSE,
                         multi_query = FALSE, significant = TRUE, exclude_iea = FALSE,
                         measure_underrepresentation = FALSE, evcodes = FALSE,
                         user_threshold = 0.05, correction_method = "g_SCS",
                         domain_scope = "annotated", custom_bg = unique(unlist(x_Hsap)),
                         numeric_ns = "", sources = "GO:BP", as_short_link = FALSE)

colnames(x_s$result)[1] <- "feature"

#Check number of enriched terms per category
tapply(x_s$result$feature,x_s$result$feature,length)

#removing redundant GO0 terms
x_s1 <- lapply(seq_len(length(CH)),function(i){
  x <- reducereduntGO(df=x_s$result[which(x_s$result$feature==CH[[i]]),],
                        GO.ID="term_id",FDR_col="p_value")
  return(x)
})

```

## Author
Main:Chrystian C. Sosa

## License
GNU GENERAL PUBLIC LICENSE Version 3
