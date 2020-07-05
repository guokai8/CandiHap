# CandiHap
<a href="https://travis-ci.org/guokai8/CandiHap"><img src="https://travis-ci.org/guokai8/CandiHap.svg" alt="Build status"></a>
[![](https://img.shields.io/badge/devel%20version-0.0.9-green.svg)](https://github.com/guokai8/CandiHap)
## Install
```
library(devtools)
install_github("guokai8/CandiHap")
```
## Usage
```
library(CandiHap)
data(snp)
gff <- snp$gff
# gff <- importGFF("test.gff3",format="gff3")
gr<-preGranges(gff,gene="Si9g49910",cds = T)
hmp <- snp$hmp
ovl <- findover(gr,hmp)
# hmp <- read.hmp("haplotypes.hmp")
pheno <- snp$pheno
# pheno <- read.pheno("Phenotype.txt",sep="\t")
hap <- snp2hap(pheno,ovl)
snplot(hap,gene="Si9g49910",side=T)
snplot(hap,gene="Si9g49910",side=T,random = F,hapname="haplotype3")
snboxplot(hap,gene="Si9g49910",feature = "test")
hapnet(hap,gene="Si9g49910",feature = "test")
# plot gene track with snp
snptrack(snp$gff,dat=snp$dat,id="Parent")
## id is the gene name you want to display, in gff3 file should be 'Parent'
#show snp locate in gene only
snptrack(snp$gff,dat=snp$dat,id="Parent",geneOnly=T)
# exon only
snptrack(snp$gff,dat=snp$dat,id="Parent",exonOnly=T)
#define different color
snptrack(snp$gff,dat=snp$dat,id="Parent",low='green',high='pink')
#use r^2 stand for color
snptrack(snp$gff,dat=snp$dat,id="Parent",color='r2')
#change name, etc ...
?snptrack
```
