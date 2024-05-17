#' @title plotnlp

#' Draw the barplot of negative log10 p-values for GPCR families
#' @param df a dataframe of results obtained with gpcrcalc function
#' @usage custom<-c("TAS1R1","TAS2R3","TAS2R4","PTGDR","PTGDR2","PTGER1","PTGER2","P2RY4","P2RY6","P2RY11","P2RY12","P2RY13","P2RY14","CXCR4","CXCR6","CXCL11","CXCR2","OPN1LW","OPN1MW","OPN1SW","RHO","OPN3","OPN4","OPN5")
#' @usage res<-gpcrcalc(custom)
#' @usage plotnlp(res)
#' @examples custom<-c("TAS1R1","TAS2R3","TAS2R4","PTGDR","PTGDR2","PTGER1","PTGER2","P2RY4","P2RY6","P2RY11","P2RY12","P2RY13","P2RY14","CXCR4","CXCR6","CXCL11","CXCR2","OPN1LW","OPN1MW","OPN1SW","RHO","OPN3","OPN4","OPN5")
#' @examples res<-gpcrcalc(custom)
#' @examples plotnlp(res)


plotnlp<-function(res){

		## load necessary packages
		  	if(!require(ggplot2)){
    	install.packages("ggplot2")
    	library(ggplot2)}
		
			if(!require(pals)){
    	install.packages("pals")
    	library(pals)}
		
		## perform the barplot
		p=ggplot(data=res,aes(x=reorder(family,ES),y=NLP,fill=family))+geom_bar(stat="identity")+
			#ylim(0,max(res$NLP)+1)+
			coord_flip()+
			theme_minimal()+
			xlab("GPCR Genesets")+
			ylab("Negative log10 of p-values")+
			scale_fill_manual(values=cols25())+
			geom_text(aes(label=round(NLP,2)),hjust=0, vjust=0.5,color="navy",position= position_dodge(0),size=5,angle=0)+
			ggtitle("") +theme(text = element_text(size = 16))+theme(legend.position="none")

		return(p)
}
