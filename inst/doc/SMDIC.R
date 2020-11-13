## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
require(SMDIC)

## ------------------------------------------------------------------------
#Flow diagram of SMDIC.
knitr::include_graphics("../inst/workflow.jpg")

## ----echo = T, results = 'hide'------------------------------------------
library(SMDIC)
#get breast cancer gene expression profile.
exp.example<-GetExampleData("exp.example")

# perform the exp2cell method. The method must be one of "xCell","ssGSEA" and "CIBERSORT".
cellmatrix.example<-exp2cell(exp.example,method="ssGSEA")


## ------------------------------------------------------------------------
#get the result of the exp2cell function
#view the first six rows and six columns of the cell abundance matrix.
head(cellmatrix.example)

## ------------------------------------------------------------------------
# get the path of the mutation annotation file.
maf <- system.file("extdata","example.maf.gz",package = "SMDIC") 

# perform the maf2matrix method.
mutmatrix.example<-maf2matrix(maffile = maf,percent = 0.01) 

#get the result of the exp2cell function
#view the first six rows and six columns of the binary mutations matrix
head(mutmatrix.example)[1:6,1:6]

## ----import, results = "hide"--------------------------------------------
# get breast cancer cell abundance matrix, which can be the result of exp2cell function.
cellmatrix<-GetExampleData("cellmatrix") 

# get breast cancer binary mutations matrix, which can be the result of maf2matrix function.
mutmatrix<-GetExampleData("mutmatrix")

# perform the function `mutcorcell`.
mutcell<-mutcorcell(cellmatrix= cellmatrix,mutmatrix = mutmatrix,fisher.adjust = TRUE) 

#get the result of the `mutcorcell` function
mutcell<-GetExampleData("mutcell")

# the binary numerical matrix which shows the immune cells driven by somatic mutant gene.
mutcell$mut_cell[1:6,1:6]

#the numerical matrix which shows the pvalue of the immune cells driven by a somatic mutant gene
#mutcell$mut_cell_p

#the numerical matrix which show the fdr of the immune cells driven by somatic mutant gene
#mutcell$mut_cell_fdr

#the character matrix which shows the cell responses of the immune cells driven by a somatic mutant gene."up" means up-regulation, "down" means down-regulation, and "0" means no significant adjustment relationship
#mutcell$mut_cell_cellresponses

## ----echo=TRUE-----------------------------------------------------------
# perform the function mutcellsummary
summary<-mutcellsummary(mutcell =mutcell,mutmatrix = mutmatrix,cellmatrix = cellmatrix)

# get the result of the mutcellsummary function
head(summary)

## ------------------------------------------------------------------------
# perform the function gene2cellsummary
gene2cellsummary(gene="TP53",method="xCell",mutcell = mutcell) 

## ----fig.height=6, fig.width=8-------------------------------------------
# load dependent package.
require(pheatmap)

# plot significant up-regulation or down-regulation cells heat map specific for breast cancer
heatmapcell(gene = "TP53",mutcell = mutcell,cellmatrix = cellmatrix,mutmatrix = mutmatrix)

## ----echo=TRUE-----------------------------------------------------------
#maf<-"dir" 
#tips: dir is the name of the mutation annotation file (MAF) format data. It must be an absolute path or the name relative to the current working directory.

#plot the waterfall for mutation genes which drive immune cells
#plotwaterfall(maffile = maf,mutcell.summary = summary,cellnumcuoff =4)

#plot the co-occurrence and mutual exclusivity plots for mutation genes that drive immune cells.
#plotCoocMutex(maffile = maf,mutcell.summary = summary,cellnumcuoff =4)

#view the result of the plotwaterfall function
knitr::include_graphics("../inst/plotwaterfall.jpeg")

#view the result of the plotCoocMutex function
knitr::include_graphics("../inst/plotCoocMutex.jpeg")

## ------------------------------------------------------------------------
# get the result of `mutcorcell` function.
mutcell<-GetExampleData("mutcell") 

# get the result of `exp2cell` function.
cellmatrix<-GetExampleData("cellmatrix") 

#get the survival data, the first column is the sample name, the second column is the survival time, and the third is the survival event.
surv<-GetExampleData("surv")

#draw Kaplan–Meier survival curves
survcell(gene ="TP53",mutcell=mutcell,cellmatrix=cellmatrix,surv=surv,palette = c("#E7B800", "#2E9FDF"))

