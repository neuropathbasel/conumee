UNDER CONSTRUCTION
================

conumeeSummary package with additional baseline correction
================

Enhanced copy-number variation analysis using Illumina DNA methylation arrays, forked from https://github.com/dstichel/conumee.

In this version of the conumee package an additional step for baseline-correction based on BAF of SNPs is included, which allows for automated scoring of copy number alterations for genomic segments. The package can be installed using <br />

```
library(devtools)
library(pastecs)
install_github("neuropathbasel/conumeeSummary")
library(conumeeSummary)

# Calling of copy number alterations for one sample can be done using
RGset <- read.metharray.exp(targets=...)
Mset <- preprocessIllumina(RGset, bg.correct=TRUE, normalize="controls")        
BAFsnps<-as.data.frame(getSnpBeta(RGset))
query.data <- CNV.load(Mset,BAFsnps)
x <- CNV.fit(query.data, CNV.data.convert(refEPIC.data@intensity,data.frame()), annoEPIC)
x <- CNV.bin(x)
x <- CNV.detail(x)
x <- CNV.segment(x)
x <- CNV.adjustbaseline(x,"BAF")
x <- CNV.evaluation(x)
CNV.write(x, what = "segments", file = paste(path, act.idat, ".BAF.segments.seg", sep=""))
```

Please see the vignette of the original conumee package (hovestadt/conumee) for more details on the method.
