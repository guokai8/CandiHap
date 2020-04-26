# CandiHap
<a href="https://travis-ci.org/guokai8/CandiHap"><img src="https://travis-ci.org/guokai8/CandiHap.svg" alt="Build status"></a>
[![](https://img.shields.io/badge/devel%20version-0.0.6-green.svg)](https://github.com/guokai8/CandiHap)
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
hap <- snp2hap(pheno,grob=ovl)
# hap <- snp2hap(file="Phenotype.txt",grob=ovl)
snplot(hap,gene="Si9g49910",side=T)
snplot(hap,gene="Si9g49910",side=T,random = F,hapname="haplotype3")
snboxplot(hap,gene="Si9g49910",feature = "test")
```
