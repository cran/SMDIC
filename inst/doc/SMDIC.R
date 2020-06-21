## ---- include = FALSE----------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup---------------------------------------------------------------
library(SMDIC)

## ------------------------------------------------------------------------
knitr::include_graphics("../inst/workflow.jpg")

## ----echo = T, results = 'hide'------------------------------------------
library(SMDIC)
library(GSVA)
exp.example<-GetExampleData("exp.example") # obtain example expression data 
cellmatrix.example<-exp2cell(exp.example,method="ssGSEA") #Cell abundance matrix,method must be one of "xCell","ssGSEA" and "CIBERSORT".


## ------------------------------------------------------------------------
#view first six rows and six colmns of cell infiltration score matrix.
cellmatrix.example[1:6,1:6]

## ------------------------------------------------------------------------

maf <- system.file("extdata","example.maf.gz",package = "SMDIC") #get path of the mutation annotation file.
mutmatrix.example<-maf2matrix(maf) 
head(mutmatrix.example)[,1:6]

## ----import, results = "hide"--------------------------------------------
#prepare data for following analysis.
cellmatrix<-GetExampleData("cellmatrix") # obtain example result from real rasult: cell abundance matrix from real data.
mutmatrix<-GetExampleData("mutmatrix")# select mutmatrix example result from real result: a binary mutations matrix
#mutcell<-mutcorcell(cell = cellmatrix,mutmatrix = mutmatrix,fisher.adjust = TRUE) ## perform the function `mutcorcell`.

## ------------------------------------------------------------------------
mutcell<-GetExampleData("mutcell") #get the result of the `mutcorcell` function
#view first ten rows and six colmns of mutcell matrix.
mutcell$mut_cell[1:6,1:6]
#mutcell$mut_cell_p
#mutcell$mut_cell_fdr
#mutcell$mut_cell_cellresponses

## ----echo=TRUE-----------------------------------------------------------
summary<-mutcellsummary(mutcell =mutcell,mutmatrix = mutmatrix,cellmatrix = cellmatrix)# The summary have four columns.The first column are gene names,the second column are the cells driven by the gene,the third column are the number of cells driven by the gene,the fourth column are mutation rates of gene.
head(summary)

## ------------------------------------------------------------------------
gene2cellsummary(gene="TP53",method="xCell",mutcell = mutcell) #a matrix shows the short name, full name, pvalue, fdr, sigtype of the cells driven by a somatic mutation

## ----fig.height=6, fig.width=8-------------------------------------------
library(pheatmap)
heatmapcell(gene = "TP53",mutcell = mutcell,cellmatrix = cellmatrix,mutmatrix = mutmatrix)

## ----echo=TRUE-----------------------------------------------------------
#file<-"dir" #dir must be an absolute path or the name  relatived to the current working directory.
#mutoncoplot(maffile = file,mutcell.summary = summary,cellnumcuoff =0)
#mutinteractions(maffile = file,mutcell.summary = summary,cellnumcuoff =0)

## ------------------------------------------------------------------------
knitr::include_graphics("../inst/plotwaterfall.jpeg")
knitr::include_graphics("../inst/plotinteractions.jpeg")

## ------------------------------------------------------------------------
mutcell<-GetExampleData("mutcell") # The result of `mutcorcell` function.
cellmatrix<-GetExampleData("cellmatrix") # Cell abundance matrix
surv<-GetExampleData("surv") # The survival data
survcell(gene ="TP53",mutcell=mutcell,cellmatrix=cellmatrix,surv=surv)

