# CandiHap
<a href="https://travis-ci.org/guokai8/CandiHap"><img src="https://travis-ci.org/guokai8/CandiHap.svg" alt="Build status"></a>

## Install
```
library(devtools)
install_github("guokai8/CandiHap")
```
## Usage
```
library(Candihap)
gff <- importGFF("test.gff3",format="gff3")
gr<-preGranges(gff,gene="Si9g49910",cds = T)
hmp <- read.hmp("haplotypes.hmp")
ovl <- findover(gr,hmp)
hap <- snp2hap("Phenotype.txt",ovl)
snplot(hap,gene="Si9g49910",side=T)
snplot(hap,gene="Si9g49910",side=T,random = F,hapname="haplotype3")
snboxplot(hap,gene="Si9g49910",feature = "test")
```
