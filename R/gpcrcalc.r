#' @title gpcrcalc

#' Calculates GPCR family enrichment on the input human gene list
#' @param custom a vector of human gene symbol containing possible GPCRs
#' @usage custom<-c("TAS1R1","TAS2R3","TAS2R4","PTGDR","PTGDR2","PTGER1","PTGER2","P2RY4","P2RY6","P2RY11","P2RY12","P2RY13","P2RY14","CXCR4","CXCR6","CXCL11","CXCR2","OPN1LW","OPN1MW","OPN1SW","RHO","OPN3","OPN4","OPN5")
#' @usage res<-gpcrcalc(custom)
#' @examples custom<-c("TAS1R1","TAS2R3","TAS2R4","PTGDR","PTGDR2","PTGER1","PTGER2","P2RY4","P2RY6","P2RY11","P2RY12","P2RY13","P2RY14","CXCR4","CXCR6","CXCL11","CXCR2","OPN1LW","OPN1MW","OPN1SW","RHO","OPN3","OPN4","OPN5")
#' @examples res<-gpcrcalc(custom)
#' @examples res



gpcrcalc<-function(custom){

	# load necessary packages

	data(genesplit)

	export <- names(genesplit)
	export <- as.data.frame(export)

	for(i in 1:length(genesplit)){
  		export$intersect[i]<-length(intersect(custom,genesplit[[i]]))}

	export$input<-length(custom)

	for(i in 1:length(genesplit)){
  		export$geneset[i]<-length(genesplit[[i]])}

	export$totaldb<-length(unique(unlist(genesplit)))


	export<-export[(export$intersect != "0"),]
	df<-export[with(export,order(-intersect)),]

	row.names(df)<-df$export
	df$export<-NULL

	res1 <- NULL
	for (i in 1:nrow(df)){
  		table <- matrix(c(df[i,1], df[i,2], df[i,3], df[i,4]), ncol = 2, byrow = TRUE)
  		o <- fisher.test(table, alternative="two.sided")$estimate
  		# save all odds in a vector
  		res1 <- c(res1,o)
	}

	df$ES <- res1
	

	res2 <- NULL
	for (i in 1:nrow(df)){
  		table <- matrix(c(df[i,1], df[i,2], df[i,3], df[i,4]), ncol = 2, byrow = TRUE)
  		p <- fisher.test(table, alternative="two.sided")$p.value
  		# save all p values in a vector
  	res2 <- c(res2,p)
	}

	df$pvalues <- res2


	df$qvalues<-p.adjust(df$pvalues,method="fdr")


	df<-df[with(df,order(pvalues)),]
	###df<-df[(df$pvalues <= 0.05),] #filtration not necessary

	df$family<-row.names(df)
	df$NLP= -log(df$pvalues,10)
	return(df)
}