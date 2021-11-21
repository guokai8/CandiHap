
#' import gff3 or gtf file
#' @title read GFF file into Granges object
#' @importFrom rtracklayer import
#' @importFrom GenomicRanges GRanges
#' @importFrom S4Vectors isSingleString
#' @param file filename for gtf or gff3
#' @param format character indicate the file format(gtf/gff3)
#' @return GRanges object
#' @author Kai Guo
#' @export
#'
importGFF <- function(file,format=c("auto","gtf","gff3")){
    format <- match.arg(format)
    if (format == "auto") {
        format <- .detect_file(file)
        if (!(format %in% c("gff3", "gff", "gtf")))
            stop(wmsg("Cannot detect whether 'file' is a GFF3 or GTF file. ",
                      "Please use the 'format' argument to specify the ",
                      "format (\"gff3\" or \"gtf\")."))
    }
    if (format == "gff3") {
        colnames <- GFF3_COLNAMES
    }
    else if (format == "gtf") {
        colnames <- GFF2_COLNAMES
    }
    else {
        colnames <- union(GFF3_COLNAMES, GFF2_COLNAMES)
    }
    message("Import genomic features from the file as a GRanges object ...\n ",
            appendLF = FALSE)
    gr <- import(file, format = format, colnames = colnames,
                 feature.type = GFF_FEATURE_TYPES)
    gr
}

#'
#' @title filter GRanges and add flanking regions
#' @importFrom dplyr filter
#' @importFrom magrittr %>%
#' @importFrom GenomicRanges GRanges
#' @importFrom GenomicRanges GRangesList
#' @importFrom S4Vectors elementMetadata
#' @importFrom rlang `!!`
#' @importFrom rlang sym
#' @param gr GRanges object
#' @param gene gene name (set to NULL if specify the chromosome name and region)
#' @param chr chromosome name
#' @param region a vector include the start and the end location
#' @param exon a logical value indicating whether use exon information or not
#' @param cds a logical value indicating whether use cds information or not
#' @param gene_name a logical value indicating whether use
#'                  gene_name or not(for gtf)
#' @param format If not missing, should be one of “gff”,“gff3” or “gtf”.
#' @return GRangesInfo object include only select genes information
#' @export
#' @author Kai Guo
preGranges <- function(gr, gene = NULL,chr=NULL,region=NULL, exon = TRUE, cds = FALSE,
                       gene_name = FALSE,format = NULL)
    {
    #gr <- .importGFF(file, type = "auto")
    if(is.null(format)){
        if("Parent" %in% names(mcols(gr))){
            format="gff3"
            }else{
                format="gtf"
            }
    }
    if(format == "gff3"){
        gene.col = "Parent"
    }else{
        if(isTRUE(gene_name)){
            gene.col = "gene_name"
        }else{
            gene.col = "gene_id"
        }
    }
    if(isTRUE(cds)){
        cols <- c("mRNA","transcript","CDS","three_prime_UTR","five_prime_UTR",
                  "start_codon","stop_codon")
    }else{
        cols <- c("mRNA","transcript","exon","three_prime_UTR","five_prime_UTR",
                  "start_codon","stop_codon")
    }
    cols <- intersect(cols,as.vector(gr$type))
    if(!is.null(chr)&!is.null(region)){
        gx <- GRanges(seqnames = chr, strand = "*", ranges = IRanges(start=region[1],end=region[2]))
        chr_r <- subsetByOverlaps(gff, gx, ignore.strand=TRUE)
        gene <- unique(sub('\\.\\d+','',unique(unlist(elementMetadata(chr_r)[,gene.col]))))
    }
    gr <- sapply(gene, function(x)gr %>% as.data.frame() %>%
                     filter(grepl(!!x,!!sym(gene.col))) %>%
                     filter(type %in% cols)%>%GRanges())%>%GRangesList()
    res <- new("GRangesInfo",
               gr = gr,
               gene = gene,
               format = format,
               mcol = gene.col,
               cols = cols)
    res
}

#' @title read hmp files and combine snp together
#' @importFrom readr read_delim
#' @importFrom GenomicRanges makeGRangesFromDataFrame
#' @param file the name of the file which the data are to be read from.
#' @param sep the field separator character.
#' @param comment character: a character vector of length one containing a
#' single character or an empty string.
#' @return data.frame
#' @author Kai Guo
#' @export
read_hmp <- function(file, sep = "\t", comment = "&"){
    hmp <- read_delim(file, delim = sep, comment = comment)
    colnames(hmp)[1]<-"seqnames"
    hmp$start <- hmp$POS
    hmp$end <- hmp$POS
    hmp <- makeGRangesFromDataFrame(hmp,keep.extra.columns = T)
    hmp
}

#' @title find overlap between target genes and hmp
#' @importFrom GenomicRanges GRangesList
#' @importFrom GenomicRanges mcols
#' @importFrom IRanges subsetByOverlaps
#' @importFrom GenomicRanges strand
#' @param gr GrangesInfo object or GrangeList
#' @param hmp data.frame include hapmap information
#' @param upstream Vectors of non-NA non-negative integers.
#' Defines the number of nucleotides toward the 5'
#' @param downstream Vectors of non-NA non-negative integers.
#' Defines the number toward the 3'end,relative to the transcription start site
#' @export
#' @return OverlapInfo object
#' @author Kai Guo
findover <- function(gr, hmp, upstream = 2000, downstream = 500){
    if(is(gr,"GRangesInfo")){
        gro <- gr@gr
        gene <- gr@gene
        format <- gr@format
        mcol <- gr@mcol
        cols <- gr@cols
    }else{
        gro <- gr
        gene <- NULL
        format <- NULL
        mcol <- NULL
        cols <- NULL
    }
    overlap <- lapply(gro, function(x)subsetByOverlaps(hmp, expandRange(x,
                            upstream, downstream)))
    # remove genes without snps
    overlap <- overlap[unlist(lapply(overlap,function(x)nrow(mcols(x))!=0))]
    gene <- intersect(gene, names(overlap))
    overlap <- overlap %>% GRangesList()
    res <- new("OverlapInfo",
               gr         = gro,
               overlap    = overlap,
               gene       = gene,
               format     = format,
               mcol       = mcol,
               upstream   = upstream,
               downstream = downstream,
               cols       = cols
               )
    res
}

#' @title extract fasta from genotype file and find haplotype
#' @importFrom readr read_delim
#' @importFrom pegas haplotype
#' @importFrom muscle muscle
#' @importFrom ape as.DNAbin
#' @importFrom Biostrings DNAStringSet
#' @importFrom GenomicRanges GRangesList
#' @importFrom BiocParallel bplapply
#' @importFrom BiocParallel MulticoreParam
#' @importFrom BiocParallel SnowParam
#' @param pheno phenotype data
#' @param grob GRanges of overlap results
#' @param hapname haplotype name
#' @param cpus number of cores
#' @export
#' @return SeqHap object
#' @author Kai Guo
snp2hap <- function(pheno,grob,hapname='haplotype',cpus=1){
    gr <-grob@overlap
    colnames(pheno)[1] <- 'sample'
    #BPPARAM = MulticoreParam(cpus)
    if (.Platform$OS.type=="windows") {
        BPPARAM <- BiocParallel::SnowParam(cpus,type="SOCK")
        } else {
        BPPARAM <- BiocParallel::MulticoreParam(cpus)
    }
    sequence <- bplapply(gr, function(x)as.data.frame(mcols(x)[,pheno$sample]),BPPARAM=BPPARAM)
    sequence <- bplapply(sequence, function(x)
        DNAStringSet(.unlist.name(sapply(x,function(y)snp2fasta(y)))),BPPARAM=BPPARAM)
    seqname <- lapply(sequence,function(x)names(x))
    ### check the snp sequence length
    cond <- names(which(unlist(lapply(sequence, function(x)max(width(x))>5))))
    sequence <- sequence[cond]
    if(length(sequence) == 0){
        stop("Gene got fewer snp !")
    }
    ss <- bplapply(sequence, function(x)as.DNAbin(muscle(x, quiet =TRUE)),BPPARAM=BPPARAM)
    haps <- bplapply(ss,function(x)haplotype(x,strict=T),BPPARAM=BPPARAM)
    hapind <- bplapply(haps, function(x)unlist(bplapply(attr(x,"index"),'[',1)),BPPARAM=BPPARAM)
    hapind <- sapply(names(hapind), function(x)seqname[[x]][hapind[[x]]],simplify = F)
    hapnames <- bplapply(haps, function(x)paste0(hapname,1:length(rownames(x))),BPPARAM=BPPARAM)
    hapfreq <-lapply(haps, function(x)as.integer(summary(x)))
    hapx <- lapply(haps, function(x)attr(x,"index"))
    haplist <- sapply(names(hapx), function(x).getlist(x,hapx,seqname,hapname),simplify = F)
    hapgr <- sapply(names(hapind), function(x)
        unlist(gr[x][, c("REF","ALT", hapind[[x]])],
               use.names = F))
    hap <- sapply(names(hapgr), function(x).colnames.mcol(hapgr[[x]],
                                                          hapnames[[x]]))
    hap <- GRangesList(hap)
    sequence <- sapply(names(hapind), function(x)
        sequence[[x]][hapind[[x]]])
    res <- new("SeqHap",
               gr = grob@gr,
               pheno = pheno,
               haplotype = hap,
               hap = haps,
               seqs = sequence,
               hapnames = hapnames,
               hapfreq = hapfreq,
               haplist = haplist,
               upstream   = grob@upstream,
               downstream  = grob@downstream,
               cols = grob@cols,
               mcol = grob@mcol,
               gene = names(gr))
    res
}

#' @title snp data frame to vector
#' @param snp data.frame with snp information
#' @export
#' @return sequence file
#' @author Kai Guo
snp2fasta <- function(snp){
    snp <- do.call(rbind,strsplit(snp, '[\\/\\|]'))
    cond <- all.equal(paste0(snp[, 1], collapse = ""),
                    paste0(snp[, 2], collapse = ""))
    if(isTRUE(cond)){
        snp <- paste0(snp[, 1], collapse = "")
    }else{
        snp <- c(paste0(snp[, 1], collapse = ""),
                paste0(snp[, 2], collapse = ""))
    }
    snp
}

#' @title prepare GWAS data
#' @importFrom readr read_delim
#' @importFrom GenomicRanges makeGRangesFromDataFrame
#' @param file the name of the file which the data are to be read from.
#' @param sep the field separator character.
#' @param comment character: a character vector of length one containing a
#' single character or an empty string.
#' @return GRanges
#' @author Kai Guo
#' @export
read_data <- function(file, sep = "\t", comment = "&"){
    dat <- read_delim(file, delim = sep, comment = comment)
    dat <- as.data.frame(dat)
    dd <- dat[,1,drop=F]
    colnames(dd)[1]<-"seqnames"
    dd$start <- dat[,2]
    dd$end <- dat[,2]
    dd <- cbind(dd,dat[,3:ncol(dat)])
    dd <- makeGRangesFromDataFrame(dd,keep.extra.columns = T)
    dd
}
