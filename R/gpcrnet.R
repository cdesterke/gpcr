#' @title gpcrnet

#' built the network of GPCR enriched
#' @param res a dataframe of results obtained with gpcrcalc function
#' @param custom a vector of human gene symbol containing possible GPCRs
#' @param layout string of layout type for the network (layout_as_star, layout_components, layout_in_circle, layout_nicely,layout_on_grid, layout_on_sphere, layout_randomly, layout_with_dh, layout_with_drl, layout_with_fr, layout_with_gem, layout_with_graphopt, layout_with_kk, layout_with_lgl, layout_with_mds)
#' @param cex size of the vertex labels
#' @param distance distance between vertex and its label
#' @usage custom<-c("TAS1R1","TAS2R3","TAS2R4","PTGDR","PTGDR2","PTGER1","PTGER2","P2RY4","P2RY6","P2RY11","P2RY12","P2RY13","P2RY14","CXCR4","CXCR6","CXCL11","CXCR2","OPN1LW","OPN1MW","OPN1SW","RHO","OPN3","OPN4","OPN5")
#' @usage res<-gpcrcalc(custom)
#' @usage gpcrnet(custom,res,layout=layout_as_star,cex=1,distance=2)
#' @examples custom<-c("TAS1R1","TAS2R3","TAS2R4","PTGDR","PTGDR2","PTGER1","PTGER2","P2RY4","P2RY6","P2RY11","P2RY12","P2RY13","P2RY14","CXCR4","CXCR6","CXCL11","CXCR2","OPN1LW","OPN1MW","OPN1SW","RHO","OPN3","OPN4","OPN5")
#' @examples res<-gpcrcalc(custom)
#' @examples gpcrnet(custom,res,layout=layout_as_star,cex=1,distance=2)


gpcrnet<-function(custom,res,layout=layout_as_star,cex=1,distance=2){

	
	## load necessary packages
		if(!require(plyr)){
    install.packages("plyr")
    library(plyr)}
		
		if(!require(dplyr)){
    install.packages("dplyr")
    library(dplyr)}

		if(!require(igraph)){
    install.packages("igraph")
    library(igraph)}
	
	## load necessary gpcr database present in the package
	data(genesplit)
	
	## intercept gpcr database with input vector
	
	inter<-rep(NA,length(genesplit))
	inter<-as.list(inter)
	for(i in 1:length(genesplit)){
		if (length(intersect(custom,genesplit[[i]]))>0){
		inter[[i]]<-intersect(custom,genesplit[[i]])}
	}

	## remove empty families after intersection
	names(inter)<-names(genesplit)
	inter<-lapply(inter, function(z){ z[!is.na(z) & z != "NA"]})
	inter<-inter[lapply(inter,length)>0]

	## row.names list intersection with enriched families

	sig=row.names(res)
	intersig<-inter[sig]

	## transform list in dataframe

	data <- ldply (intersig, data.frame)
	colnames(data)<-c("family","gene")


	data%>%inner_join(res,by="family")->final
	final%>%arrange(desc(NLP))%>%
		mutate(significance=ifelse(pvalues<=0.05,"YES","no"))->final


	## prepare results data for network analysis
	final%>%select(family,gene,NLP,significance)->new	
	colnames(new)<-c("from","to","weight","significance")
	
	## replace significance data by distinct color for edges
	new%>%mutate(significance = replace(significance,significance=="YES","red"))%>%
		mutate(significance = replace(significance,significance=="no","gray30"))->new
		
	## build the network
	g <- graph_from_data_frame(new, directed = FALSE)
  
	## performed louvain communities detection
	lc <- cluster_louvain(g)
	membership(lc)
	communities(lc)

    ## plot the network with communities and significance
	plot(lc,g,edge.width=new$weight*1,edge.color=new$significance,vertex.color="black",
	vertex.label.cex=cex,vertex.label.dist=distance,vertex.label.family="sans",layout=layout)

}