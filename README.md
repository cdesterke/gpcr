# gpcr
R package to perform a GPCR network enrichment starting from a human gene list


### package installation
```r
library(devtools)
install_github("cdesterke/gpcr")
```

### description

Starting from an imput human gene list obtained by omics experiments this gpcr package allows to perform a GPCR family enrichment.

- In a first step as compared to human GPCR database, lastly updated in december 2022: https://www.guidetopharmacology.org/download.jsp. This initial step of enrichment performs with "gpcrcalc" function allows to output a result table with significance of each GPCR family found in your gene list. 

- In second step you could output your enrichment results in barplot for "enrichment scores with "plotes" function and for negative log10 p-values with "plotnlp" function. 

- In last step you can performed a GPCR network with Louvain communitie detection: fonction "gpcrnet".

## code for processing analysis

```r
##define a human gene list
custom<-c("TAS1R1","TAS2R3","TAS2R4","PTGDR","PTGDR2","PTGER1","PTGER2","P2RY4","P2RY6","P2RY11","P2RY12","P2RY13","P2RY14","CXCR4","CXCR6","CXCL11","CXCR2","OPN1LW","OPN1MW","OPN1SW","RHO","OPN3","OPN4","OPN5")
## calculates the GPC family enrichments
res<-gpcrcalc(custom)
res
```
