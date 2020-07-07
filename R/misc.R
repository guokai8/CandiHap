
GFF3_COLNAMES <- c("source","type","phase","ID","Name","Parent","GeneID")
GFF2_COLNAMES<-c("source","type","phase",
                 "gene_id","gene_name","transcript_id","transcript_name")

#' copy from GenomicFeatures
.GENE_TYPES <- c("gene", "pseudogene", "transposable_element_gene")

.TX_TYPES <- c("transcript", "pseudogenic_transcript", "primary_transcript",
               "mRNA", "ncRNA", "rRNA", "snoRNA", "snRNA", "tRNA", "tmRNA",
               "miRNA", "miRNA_primary_transcript",
               "RNase_P_RNA", "RNase_MRP_RNA", "SRP_RNA", "misc_RNA",
               "antisense_RNA", "antisense",
               "lnc_RNA", "antisense_lncRNA", "transcript_region",
               "pseudogenic_tRNA", "scRNA", "guide_RNA", "telomerase_RNA",
               "vault_RNA", "Y_RNA")

.EXON_TYPES <- c("exon", "pseudogenic_exon", "coding_exon",
                 "five_prime_coding_exon", "interior_coding_exon",
                 "three_prime_coding_exon", "exon_of_single_exon_gene",
                 "interior_exon", "noncoding_exon",
                 "five_prime_noncoding_exon", "three_prime_noncoding_exon")
.CDS_TYPES <- c("CDS", "transposable_element_CDS",
                "CDS_predicted", "edited_CDS")
.UTR_TYPES <- c("five_prime_UTR", "three_prime_UTR","UTR")

.STOP_CODON_TYPES <- c("start_codon","stop_codon")

### cory from GenomicFeatures
GFF_FEATURE_TYPES <- c(.GENE_TYPES, .TX_TYPES, .EXON_TYPES,
                       .CDS_TYPES, .STOP_CODON_TYPES,.UTR_TYPES)
GENE_FEATURE <- c(.EXON_TYPES,.CDS_TYPES,.UTR_TYPES,.STOP_CODON_TYPES)
#' detect file format for import gff or gff3 file
#' modified based on GenomicFeatures .detect_file_format function
#' @importFrom tools file_ext
#' @importFrom tools file_path_sans_ext
#' @importFrom rtracklayer FileForFormat
.detect_file <- function(file)
{
    if (isSingleString(file)) {
        file2 <- try(FileForFormat(file), silent=TRUE)
        if (inherits(file2, "try-error"))
            return(file_ext(file))
        file <- file2
    }
    if (is(file, "RTLFile")) {
        if (is(file, "GFF3File"))
            return("gff3")
        if (is(file, "GTFFile"))
            return("gtf")
        desc <- rtracklayer:::resourceDescription(file)
        if (is(file, "CompressedFile")) {
            ## Mimic what import,CompressedFile,missing,missing method does.
            desc <- file_path_sans_ext(desc)
        }
        type <- file_ext(desc)
        return(type)
    }
    stop(wmsg("Invalid 'file'. Must be a path to a file, or an URL, ",
              "or a GFF3File or GTFFile object."))
}
#' expand range from start and end
#' @importFrom BiocGenerics start
#' @importFrom BiocGenerics end
#' @importFrom BiocGenerics which
#' @importFrom GenomicRanges `start<-`
#' @importFrom GenomicRanges `end<-`
#' @param gr GRanges object
#' @param upstream flanking length for upstream
#' @param downstream flanking length for downstream
#' @export
#' @author Kai Guo
expandRange = function(gr, upstream=1000, downstream=1000) {
    min = start(range(gr))
    max = end(range(gr))
    strand_is_minus = strand(gr) == "-"
    on_plus = which(!strand_is_minus)
    on_minus = which(strand_is_minus)
    start(gr)[on_plus & start(gr) == min] = start(gr)[on_plus &
                        start(gr) == min] - upstream
    end(gr)[on_minus & end(gr) == max] = end(gr)[on_minus &
                        end(gr) == max] + upstream
    end(gr)[on_plus & end(gr) == max] = end(gr)[on_plus & end(gr) == max] +
        downstream
    start(gr)[on_minus & start(gr) == min] = start(gr)[on_minus &
                        start(gr) == min] - downstream
    gr
}

#' unlist and add names
#' @param list a list include character vector
#' @author Kai Guo
.unlist.name <- function(list){
    list.name <- rep(names(list),times=unlist(lapply(list,length)))
    list <- unlist(list)
    names(list) <- list.name
    list
}
#' change mcols data colnames
#' @importFrom GenomicRanges GRanges
#' @importFrom S4Vectors `elementMetadata<-`
#' @importFrom S4Vectors elementMetadata
#' @param gr GRanges object
#' @param names new column names
#' @author Kai Guo
.colnames.mcol <- function(gr,names){
    #### need to be change
    names(elementMetadata(gr))[3:ncol(elementMetadata(gr))] <- names
    gr
}
#' @title anova test for phenotype based on haplotype
#' @importFrom stats aov na.omit TukeyHSD
#' @importFrom dplyr left_join
#' @importFrom broom tidy
#' @importFrom magrittr %>%
#' @importFrom scales scientific
#' @importFrom tidyr separate
#' @param haps SeqHap object
#' @param gene gene name
#' @param feature pheaotype name
#' @export
#' @author Kai Guo
#'
aov.pheno.test <- function(haps, gene, feature){
    haps <- hap@haplotype[[gene]]
    pheno <- hap@pheno[,c("sample",feature)]
    pheno <- na.omit(pheno)
    haplist <- hap@haplist[[gene]]
    hapl <- rep(names(haplist), times = unlist(lapply(haplist, length)))
    hapx <- data.frame(sample=unlist(haplist),group=hapl)
    pheno <- left_join(pheno, hapx, by=c("sample"="sample"))
    sel <- names(which(table(pheno$group)>=3))
    ###select haplotype with three phenotype values
    pheno <- subset(pheno,group%in%sel)
    colnames(pheno)[2]<-"feature"
    ptest <- tidy(TukeyHSD(aov(feature~group,data=pheno)))%>%
        separate(comparison,c("group1","group2"),sep="-")
    ptest$adj.p.value<-scientific(ptest$scientific,digits=2)
    ptest
}


#' @title extract sequence form SeqHap object
#' @param object SeqHap object
#' @param gene gene
#' @export
#' @author Kai Guo
setMethod("seqs", signature = (object="SeqHap"), function(object,gene=NULL){
    res <- object@seqs[[gene]]
    name <- names(object@haplist[[gene]])
    names(res)<-name
    res
})

#' @title extract final results from SeqHap object
#' @param object SeqHap object
#' @param gene gene
#' @importFrom GenomicRanges as.data.frame
#' @export
#' @author Kai Guo
setMethod("results", signature = (object="SeqHap"), function(object,gene=NULL){
    if(!is.null(gene)){
        res <- object@haplotype[[gene]]
    }else{
        res <- unlist(object@haplotype)
    }
    res <- as.data.frame(res)
    return(res)
})

#' @title read hmp files and combine snp together
#' @importFrom readr read_delim
#' @param file the name of the file which the data are to be read from.
#' @param sep the field separator character.
#' @param comment character: a character vector of length one containing a
#' single character or an empty string.
#' @return data.frame
#' @author Kai Guo
#' @export
read_pheno <- function(file, sep = "\t"){
    pheno <- read_delim(file, delim = sep)
    colnames(pheno)[1] <- 'sample'
    pheno
}

#' @title get name of list
#' @param name gene name
#' @param x haplist
#' @param y name
#' @param hapname haplotype name
.getlist <- function(name,x,y,hapname){
    tx <- x[[name]]
    names(tx) <- paste0(hapname,1:length(tx))
    ty <- y[[name]]
    tmp <- lapply(tx,function(x)ty[x])
    return(tmp)
}

###
.color_scale <- function(c1="pink", c2="red") { #modified from DOSE
    pal <- colorRampPalette(c(c1, c2))
    colors <- pal(200)
    return(colors)
}
.getIdx <- function(v, MIN, MAX) { #modified from DOSE
    intervals <- seq(MIN, MAX, length.out=200)
    max(which(intervals <= v))
}

#' check if REF same as hap
#' @importFrom magrittr %>%
#' @importFrom dplyr mutate
#' @importFrom GenomicRanges as.data.frame
#' @importFrom GenomicRanges mcols
#' @importFrom GenomicRanges `mcols<-`
#' @param gr GRanges object
.check_same_or_not <- function(gr){
    gr <- mcols(gr)
    gr <- as.data.frame(gr)
    gr$alt <- sub('[\\/\\|].*','',gr[,2])
    gr$cond <- gr$alt!=gr$REF
    return(gr$cond)
}
#' @title read vcf files
#' @importFrom pegas read.vcf
#' @importFrom pegas VCFloci
#' @importFrom GenomicRanges makeGRangesFromDataFrame
#' @importFrom GenomicRanges mcols
#' @importFrom GenomicRanges `mcols<-`
#' @param file A ‘VcfFile’ (synonymous with ‘TabixFile’) instance or
#' character() name of the VCF file to be processed..
#' @export
#' @return GRanges object
#' @author Kai Guo
read_vcf <- function(file) {
    info <- VCFloci(file, quiet = TRUE)
    snp <- makeGRangesFromDataFrame(info, start.field = 'POS',
                end.field = 'POS', keep.extra.columns = TRUE)
    snp <- snp[, c('REF', 'ALT')]
    vcf <- read.vcf(file, from = 1,to = nrow(info), quiet = TRUE)
    mcols(snp)<-cbind(mcols(snp), t(as.data.frame(vcf)))
    snp
}
