# gpcr
R package to perform a GPCR network enrichment starting from a human gene list


### package installation
```r
library(devtools)
install_github("cdesterke/gpcr")
```

### description

Starting from an input human gene list obtained by omics experiments this gpcr package allows to perform a GPCR family enrichment.

- In a first step as compared to human GPCR database, lastly updated in december 2022: https://www.guidetopharmacology.org/download.jsp. This initial step of enrichment performs with "gpcrcalc" function allows to output a result table with significance of each GPCR family found in your gene list. 

- In second step you could output your enrichment results in barplot for "enrichment scores with "plotes" function and for negative log10 p-values with "plotnlp" function. 

- In last step you can performed a GPCR network with Louvain communitie detection: fonction "gpcrnet". Significant families during the enrichment were represented with red edges and not significant families with gray edges. 

## code for processing analysis

### compute enrichment

```r
##define a human gene list
custom<-c("TAS1R1","TAS2R3","TAS2R4","PTGDR","PTGDR2","PTGER1","PTGER2","P2RY4","P2RY6","P2RY11","P2RY12","P2RY13","P2RY14","CXCR4","CXCR6","CXCL11","CXCR2","OPN1LW","OPN1MW","OPN1SW","RHO","OPN3","OPN4","OPN5")
## calculates the GPCR family enrichments
res<-gpcrcalc(custom)
res
```


![res](https://github.com/cdesterke/gpcr/blob/main/res.png)

### barplot of enrichment scores

```r
plotes(res)

```
![es](https://github.com/cdesterke/gpcr/blob/main/es.png)


### barplot of negative log10 p-values

```r
plotnlp(res)

```
![nlp](https://github.com/cdesterke/gpcr/blob/main/nlp.png)


## network of GPCR enriched
default parameters

```r
gpcrnet(custom,res,layout=layout_as_star,cex=1,distance=2)

```
![default](https://github.com/cdesterke/gpcr/blob/main/netdefault.png)

### custom the network with some parameters

- cex parameter: change size of the vertex (nodes) label

- distance parameter: change the distance between the vertex and its label

- layout parameter: change the design of the network and have several options such as: (layout_as_star, layout_components, layout_in_circle, layout_nicely,layout_on_grid,
layout_on_sphere, layout_randomly, layout_with_dh, layout_with_drl, layout_with_fr, layout_with_gem,
layout_with_graphopt, layout_with_kk, layout_with_lgl, layout_with_mds)

```r
gpcrnet(custom,res,layout=layout_nicely,cex=0.8,distance=1.5)

```
![nicely](https://github.com/cdesterke/gpcr/blob/main/nicely.png)
