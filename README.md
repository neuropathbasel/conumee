UNDER CONSTRUCTION
================

conumeeSummary package with additional baseline correction
================

Enhanced copy-number variation analysis using Illumina DNA methylation arrays.

In this version of the conumee package an additional step for baseline-correction based on BAF of SNPs is included, which allows for automated scoring of copy number alterations for genomic segments. The package can be installed using <br />
library(devtools) <br />
library(pastecs) <br />
install_github("neuropathbasel/conumeeSummary") <br />

Calling of copy number alterations for one sample can be done using <br />
RGset <- read.metharray.exp(targets=...) <br />
Mset <- preprocessIllumina(RGset, bg.correct=TRUE, normalize="controls") <br />             
BAFsnps<-as.data.frame(getSnpBeta(RGset)) <br />
query.data <- CNV.load(Mset,BAFsnps) <br />
x <- CNV.fit(query.data, CNV.data.convert(refEPIC.data@intensity,data.frame()), annoEPIC) <br />
x <- CNV.bin(x) <br />
x <- CNV.detail(x) <br />
x <- CNV.segment(x) <br />
x <- CNV.adjustbaseline(x,"BAF") <br />
x <- CNV.evaluation(x) <br />
CNV.write(x, what = "segments", file = paste(path, act.idat, ".BAF.segments.seg", sep="")) <br />
 
Please see the vignette of the original conumee package (hovestadt/conumee) for more details on the method.
