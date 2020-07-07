# CandiHap
<a href="https://travis-ci.org/guokai8/CandiHap"><img src="https://travis-ci.org/guokai8/CandiHap.svg" alt="Build status"></a>
[![](https://img.shields.io/badge/devel%20version-0.0.12-green.svg)](https://github.com/guokai8/CandiHap)
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
# hmp <- read_hmp("haplotypes.hmp")
# or if you use vcf file
# hmp <- read_vcf('test.vcf')
pheno <- snp$pheno
# pheno <- read_pheno("Phenotype.txt",sep="\t")
hap <- snp2hap(pheno,ovl)
## want to extract results
results(hap, gene="Si9g49910")
# make figures
snplot(hap,gene="Si9g49910",side=T)
snplot(hap,gene="Si9g49910",side=T,random = F,hapname="haplotype3")
snplot(hap,gene="Si9g49910",side=T,random = F,hapname="haplotype15",mutateOnly=TRUE)
## want to know more
?snplot
snboxplot(hap,gene="Si9g49910",feature = "test")
hapnet(hap,gene="Si9g49910",feature = "test")
# plot gene track with snp
dat=snp$dat
# or read your data
# notice that the dat should have the chromosome name in the first column, position in the second column and the values in the following columns 
# dat <- read_data(file,sep="\t")
snptrack(gff,dat=dat,id="Parent")
## id is the gene name you want to display, in gff3 file should be 'Parent'
#show snp locate in gene only
snptrack(gff,dat=dat,id="Parent",geneOnly=T)
# exon only
snptrack(gff,dat=dat,id="Parent",exonOnly=T)
#define different color
snptrack(gff,dat=dat,id="Parent",low='green',high='pink')
#use r^2 stand for color
snptrack(gff,dat=dat,id="Parent",color='r2')
## show some genes
snptrack(gff,dat=dat,id="Parent",gene='Si9g49910')
#change name, etc ...
?snptrack
```
