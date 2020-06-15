##' Class "SeqHap"
##' This class represents the result of overlap between snp and target gene
##'
##'
##' @name SeqHap-class
##'
##' @docType class
##' @slot gr GRangeList obejct
##' @slot pheno data frame of phenotype data
##' @slot haplotype list of haplotype object
##' @slot gene vector of character of genes
##' @exportClass SeqHap
##' @author Kai Guo
##' @keywords classes
setClass("SeqHap",
         representation=representation(
             gr         = "GRangesList",
             pheno         = "data.frame",
             haplotype   = "GRangesList",
             hap = "list",
             seqs = "list",
             hapnames = "list",
             hapfreq = "list",
             hapnet = "list",
             haplist = "list",
             upstream   = "numeric",
             downstream  = "numeric",
             cols = "character",
             mcol = "character",
             gene  = "character"
         )
)
##

##' Class "GRangesInfo"
##' This class represents the extract gene infomation from
##' @slot gr GrangesList of gff or gff3 files with filter genes
##' @slot gene vector of character of genes
##' @slot gene Gene IDs
##' @slot format a character gff/gtf or gff3
##' @slot mcol a character specify the gene id column
##' @slot cols a vector include gene structure information
##' @exportClass GRangesInfo
##' @author Kai Guo
##' @keywords classes
setClass("GRangesInfo",
         representation=representation(
             gr         = "GRangesList",
             gene       = "character",
             format     = "character",
             mcol       = "character",
             cols       = "character"
         )
)

##
##' Class "OverlapInfo"
##' This class represents the extract gene infomation from
##' @slot gr GrangesList of gff or gff3 files with filter genes
##' @slot overlap GRangesList
##' @slot gene vector of character of genes
##' @slot gene Gene IDs
##' @slot format a character gff/gtf or gff3
##' @slot mcol a character specify the gene id column
##' @slot cols a vector include gene structure information
##' @exportClass OverlapInfo
##' @author Kai Guo
##' @keywords classes
setClass("OverlapInfo",
         representation=representation(
             gr         = "GRangesList",
             overlap    = "GRangesList",
             gene       = "character",
             format     = "character",
             mcol       = "character",
             upstream   = "numeric",
             downstream  = "numeric",
             cols       = "character"
         )
)
