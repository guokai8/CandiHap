#' @title plot snp with gene
#' @importFrom trackViewer lolliplot
#' @importFrom trackViewer dandelion.plot
#' @importFrom GenomicRanges mcols
#' @importFrom GenomicRanges `mcols<-`
#' @importFrom VennDetail setcolor
#' @importFrom grid grid.text
#' @importFrom grid gpar
#' @param hap SeqHap object include a GRanges, GRangesList
#' @param gene character gene
#' @param hapname haplotype want to display
#' @param layout a character to indiate layout method (lolliplot or
#' dandelion.plot)
#' @param mutateOnly only show the alt is different with the REF or not
#' @param snp.colors color for the snp, random colors will be used if NULL
#' @param angle a numeric indicate the angle for the labels
#' @param side use both top and bottom to display the SNP (TRUE/FALSE)
#' @param type character. Could be circle, pie, pin, pie.stack or flag.
#' @param rescale rescale the x-axis or not.
#' @param snp.alpha alpha value for color
#' @param jitter jitter the position of nodes or labels.
#' @param random use random height to display SNP or not
#' @param snp.lwd The line width for displaying SNP, a _positive_ number.
#' @param label.cex label size
#' @param upstream Defines the number of nucleotides toward the 5'
#' @param downstream Defines the number of nucleotides toward the 3'
#' @export
#' @return figure
#' @author Kai Guo
snplot <- function(hap, gene, hapname = NULL,layout = "lolliplot",
                   mutateOnly = FALSE,
                   snp.colors = NULL, angle = 90, side = FALSE,
                   type="circle",rescale = FALSE,snp.alpha = 0.7,
                   jitter ="node",
                   random = TRUE,snp.lwd = 0.5,label.cex = 0.8,
                   upstream = NULL, downstream = NULL){
    set.seed(123)
    cols <- hap@cols
    mcol <-hap@mcol
    grs <- hap@gr[[gene]]
    gr <- grs
    mcols(gr)<-NULL
    times <- c(length(.TX_TYPES), length(c(.EXON_TYPES,.CDS_TYPES)),1,
               2,1,1)
    gcol <- rep(c("black","#1F78B4","brown","#DFA32D","brown","darkcyan"),times = times)
    names(gcol) <- c(.TX_TYPES, .EXON_TYPES,
                     .CDS_TYPES, .UTR_TYPES,.STOP_CODON_TYPES)
    hei <- rep(c(0.001,0.04,0.02,0.02,0.02,0.02), times = times)
    names(hei) <- c(.TX_TYPES, .EXON_TYPES,
                    .CDS_TYPES, .STOP_CODON_TYPES,.UTR_TYPES)
    gr$height<-hei[as.character(grs$type)]
    gr$fill<-gcol[as.character(grs$type)]
    gr$featureLayerID<-unlist(mcols(grs)[,mcol])
    names(gr)<-grs$type
    if(is.null(hapname)){
        haps <- hap@haplotype[[gene]][,c("REF","ALT")]
        names(haps)<-paste(haps$REF,haps$ALT,sep="->")
    }else{
        haps <- hap@haplotype[[gene]][,c("REF",hapname)]
        if(isTRUE(mutateOnly)){
            snpx <- .check_same_or_not(haps)
          if(sum(snpx)==0){
              print('All ALT is same as REF')
          }else{
              haps <- haps[snpx,]
          }
        }
        names(haps)<-paste(haps$REF,sub('[\\/\\|].*','',mcols(haps)[,hapname]),sep="->")
    }
    if(isTRUE(side) & layout == "lolliplot"){
        haps$SNPsideID <- sample(c("top", "bottom"),
                            length(haps),
                            replace=TRUE)
    }
    if(is.null(snp.colors)){
        haps$color <- setcolor(length(haps))
    }else{
        haps$color <- rep(snp.colors,length(haps))
    }

    haps$alpha <- snp.alpha
    if(is.null(upstream)){
        upstream <- hap@upstream
    }
    if(is.null(downstream)){
        downstream <- hap@downstream
    }
    haps$label.parameter.rot <- angle
    haps$lwd <- snp.lwd
    if(isTRUE(random) & layout == "lolliplot"){
        haps$score <- runif(length(haps))*length(haps)
    }

    if(type == "flag"){
        haps$label <- names(haps)
        names(haps) <- NULL
        haps$label.rot <- angle/5
    }
    gr.exp <- c(expandRange(range(gr),upstream,downstream),gr)
    gr.exp$height[1]<-gr$height[1]/10
    if("mRNA"%in%cols){
        names(gr.exp)[names(gr.exp) == "mRNA"] <- "genome"
        names(gr.exp)[1]<-"genome"
    }else if("transcript" %in% cols){
        names(gr.exp)[names(gr.exp) == "transcript"] <- "genome"
        names(gr.exp)[1]<-"gene"
    }else{
        names(gr.exp)[1]<-"genome"
    }
   # gr.exp$color[1]<-"black"
    gr.exp$fill[1]<-"black"
    if(layout=="dandelion"){
        maxgaps <- tile(gr.exp, n = 10)[[1]]
        dandelion.plot(haps,gr.exp,type=type,jitter = jitter,cex=label.cex,
                ylab="",rescale = rescale, xaxis = FALSE,maxgaps=maxgaps,
                yaxis = FALSE)
    }else{
        lolliplot(haps,gr.exp,type=type,jitter = jitter,cex=label.cex,
                ylab="",rescale = rescale, xaxis = FALSE,
                yaxis = FALSE)
    }
    if(isTRUE(as.character(strand(gr.exp))[1]=="+")){
        label <- paste0(gene," (5'->3')")
    }else{
        label <- paste0(gene," (3'<-5')")
    }
    grid.text(label, x=.20, y=.9, just="top",
              gp=gpar(cex=1, fontface="bold"))
}

#' @title boxplot for haplotype testing based on phenotype
#' @importFrom dplyr left_join
#' @importFrom dplyr mutate
#' @importFrom stats aov na.omit TukeyHSD
#' @importFrom tidyr separate
#' @importFrom broom tidy
#' @importFrom VennDetail setcolor
#' @importFrom scales scientific
#' @importFrom ggpubr ggboxplot
#' @importFrom ggpubr stat_pvalue_manual
#' @importFrom ggplot2 theme_light
#' @importFrom ggplot2 scale_fill_manual
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 xlab
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 element_text
#' @param hap SeqHap object include a GRanges, GRangesList
#' @param gene gene name
#' @param feature phenotype name
#' @export
#' @author Kai Guo
snboxplot <- function(hap, gene, feature){
    haps <- hap@haplotype[[gene]]
    pheno <- hap@pheno[,c("sample",feature)]
    pheno <- na.omit(pheno)
    haplist <- hap@haplist[[gene]]
    #hapn <- strsplit(as.character(haplist[, 1]), ";")
    hapl <- rep(names(haplist), times = unlist(lapply(haplist, length)))
    hapx <- data.frame(sample=unlist(haplist),group=hapl)
    pheno <- left_join(pheno, hapx, by=c("sample"="sample"))
    sel <- names(which(table(pheno$group)>=3))
    ###select haplotype with three phenotype values
    pheno <- subset(pheno,group%in%sel)
    colnames(pheno)[2]<-"feature"
    ptest <- tidy(TukeyHSD(aov(feature~group,data=pheno)))%>%
        separate(contrast,c("group1","group2"),sep="-")
    ptest <- subset(ptest,adj.p.value < 0.05)
    if(nrow(ptest)<1){
        stop("No significant found!")
    }
    ran <- nrow(ptest)
    ptest$adj.p.value<-scientific(ptest$adj.p.value,digits=2)
    unin <- length(unique(c(ptest$group1,ptest$group2)))
    y.position <- max(pheno$feature)+seq(1,ran*max(pheno$feature)/10,
                                    max(pheno$feature)/10)
    p <- ggboxplot(pheno,x="group",y="feature",fill="group",add="jitter",
              add.params = list(size=1,alpha=0.7))+xlab("")+
        scale_fill_manual(values=setcolor(unin))+theme_light(base_size=15)+
        theme(axis.text.x=element_text(angle=90))+
    stat_pvalue_manual(ptest,label = "p={adj.p.value}",y.position = y.position)
    print(p)
}
#'
#' plot haplotype with network
#' @importFrom pegas as.igraph.haploNet
#' @importFrom pegas haploNet
#' @importFrom igraph E
#' @importFrom igraph E<-
#' @importFrom igraph V
#' @importFrom igraph V<-
#' @importFrom igraph graph.edgelist
#' @importFrom GGally ggnet2
#' @importFrom magrittr %>%
#' @importFrom ggrepel geom_text_repel
#' @importFrom dplyr group_by summarise
#' @importFrom tidyr gather
#' @export
hapnet <- function(hap, gene, feature = NULL, freq = TRUE,
                   cex=1, altlinks=FALSE,high="red",low="skyblue",log=TRUE,
                   edge.label.size = 3,node.label.size=4,
                   node.alpha=0.75,...){
  require(igraph)
  happ<-hap@hap[[gene]]
  hapnet <- haploNet(happ)
    #hapnet <- hap@hapnet[[gene]]
  #  happ <- hap@hap[[gene]]
    hapnames <- hap@hapnames[[gene]]
    haplist <- hap@haplist[[gene]]
    pheno <- hap@pheno
    colnames(pheno)[1]<-"sample"
    attr(hapnet,"labels") <- hapnames
#    if(isTRUE(default)){
#        if(isTRUE(altlinks)){
#            par(mar=c(1,1,1,1))
#            plot(hapnet,cex=cex,...)
#        }else{
#            par(mar=c(1,1,1,1))
#            plot(hapnet,cex=cex,threshold=0,...)
#        }
 #   }else{
    alter <- attr(hapnet,"alter.links")
    #get links weight
    links <- as.data.frame(hapnet[,1:4])
    freq <- attr(hapnet,"freq")
    names(freq)<-hapnames
    g <- as.igraph.haploNet(hapnet,altlinks = altlinks)
    if(isTRUE(altlinks)){
        links <- rbind(links,alter)
    }
    E(g)$label<-links$step
    cols <- .color_scale(low, high)
    if(isTRUE(freq)){
        px <- freq
    }else{
       if(!is.null(feature)){
           pp <- pheno[, c("sample",feature)]
       }else{
           pp <- pheno
       }
           pp <- pp%>%gather(id,val,-sample)%>%
               group_by(sample)%>%summarise(mu=mean(val,na.rm=T))
           pp <- as.data.frame(pp)
           rownames(pp)<-pp$sample
           px <- unlist(lapply(haplist, function(x)mean(pp[x,"mu"],na.rm=T)))
           px[is.na(px)] <- 1E-10
    }
    px <- px[V(g)$name]
    if(isTRUE(log)){
        size <- log2(freq[V(g)$name])
    }
    V(g)$size <- size
    V(g)$color <- cols[sapply(px, .getIdx, min(px), max(px))]
    ggnet2(g,edge.label = E(g)$label,node.size = V(g)$size,
           legend.position = "none",
           node.color = V(g)$color,node.alpha = node.alpha)+
    geom_text_repel(label=V(g)$name)
#    }
}


