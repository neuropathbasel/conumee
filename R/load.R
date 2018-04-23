##### LOADING methods #####

#' CNV.load
#' @description Prepare combined intensities from various input objects.
#' @param input Object of MethylSet class (minfi package), data.frame class, matrix class or numeric class.
#' @param names Vector specifying sample names. If not supplied, colnames are used. For MethylSet input, the first column of pData(input) matching 'name' (grep) is used.
#' @param ... Additional parameters (\code{CNV.load} generic, currently not used).
#' @return \code{CNV.data} object.
#' @details This method gathers combined intensities of the Methylated and Unmethylated signals for all supplied probes. Probe IDs must be supplied as row names or in a seperate column named `ID_REF` or `TargetID`.
#' If column names match 'intensity', only those columns are used. Else, if column names match 'signal' or 'methylated', only those columns are used. Otherwise, all columns are used.
#' @examples
#' library(minfiData)
#' d <- CNV.load(MsetEx)
#' d
#' @author Volker Hovestadt \email{conumee@@hovestadt.bio}
#' @export
setGeneric("CNV.load", function(input, BAFsnps, ...) {
    standardGeneric("CNV.load")
})

#' @rdname CNV.load
setMethod("CNV.load", signature(input = "MethylSet", BAFsnps="data.frame"), function(input, BAFsnps, names = NULL) {
    message("Loading data")
    object <- new("CNV.data")
    methy.ba <- minfi::getMeth(input)
    unmethy.ba <- minfi::getUnmeth(input)
    object@intensity <- as.data.frame(methy.ba + unmethy.ba)

    data(BAFsnps_positions)
    BAFsnps.names<-rownames(BAFsnps)
    #BAFsnps.names<-rownames(input)[grepl("rs",rownames(input))]
    
    #object@BAFsnps<-data.frame(probe=BAFsnps.names,BAF=methy.ba[rownames(betas)%in%BAFsnps.names] / (methy.ba[rownames(betas)%in%BAFsnps.names] +unmethy.ba[rownames(betas)%in%BAFsnps.names] +100))
    object@BAFsnps<-data.frame(probe=BAFsnps_positions$probe,chrom=BAFsnps_positions$chrom,start=BAFsnps_positions$start,end=BAFsnps_positions$end,BAF=BAFsnps[match(BAFsnps_positions$probe,BAFsnps.names),] ) 
    
    input.names <- grep("Name", setdiff(colnames(minfi::pData(input)), 
        c("Basename", "filenames")), ignore.case = TRUE)
    if (length(input.names) > 0) 
        names(object) <- minfi::pData(input)[, grep("name", setdiff(colnames(minfi::pData(input)), 
            c("Basename", "filenames")), ignore.case = TRUE)[1]]
    if (!is.null(names)) 
        names(object) <- names
    
    
    #Method for gender prediction (ex minfi and MNP)
    MNPgetSex <- function(object = NULL, cutoff = -2){
      #.myisGenomicOrStop(object)
      #if(is(object, "GenomicMethylSet"))
      #  CN <- getCN(object)
      #if(is(object, "GenomicRatioSet"))
      CN <- getCN(object)
      ## FIXME: add test for logarithmic scale or non-log scale
      xIndex <- which(as.logical(seqnames(object) == "chrX"))
      yIndex <- which(as.logical(seqnames(object) == "chrY"))
      out <- .mygetSex(CN = CN, xIndex = xIndex,
                       yIndex = yIndex, cutoff = cutoff)
      return(out)
    }
    .mygetSex <- function(CN = NULL, xIndex = NULL, yIndex = NULL, cutoff=-2) {
      if(is.null(CN) | is.null(xIndex) | is.null(yIndex))
        stop("must provide CN, xIndex, and yIndex")
      ## FIXME: does not handle only females or only males
      ## this ought to be handled by the 'centers' (see below) being too close together
      xMed <- matrixStats::colMedians(CN, rows = xIndex, na.rm=TRUE)
      yMed <- matrixStats::colMedians(CN, rows = yIndex, na.rm=TRUE)
      dd <- yMed - xMed
      sex0 <- ifelse(dd < cutoff, "F", "M")
      df <- DataFrame(xMed = xMed, yMed = yMed, predictedSex = sex0)
      rownames(df) <- colnames(CN)
      df
    }
    
    #Determine gender
    object@gender<-MNPgetSex(minfi::mapToGenome(input))[1,3]
    
    object <- CNV.check(object)
    
    return(object)
})

#' @rdname CNV.load
setMethod("CNV.load", signature(input = "data.frame", BAFsnps="data.frame"), function(input, BAFsnps,
    names = NULL) {
    object <- new("CNV.data")
    object@date <- date()
    object@BAFsnps<-data.frame()
    
    if (any(grepl("TargetID", colnames(input)))) 
        rownames(input) <- input[, "TargetID"]
    if (any(grepl("ID_REF", colnames(input)))) 
        rownames(input) <- input[, "ID_REF"]
    if (is.null(rownames(input))) 
        stop("intensities not given for all probes.")
    
    if (any(grepl("intensity", colnames(input), ignore.case = TRUE))) {
        input.i <- grep("intensity", colnames(input), ignore.case = TRUE)
        input.n <- sapply(strsplit(colnames(input), "\\.", colnames(input))[input.i], 
            function(x) paste(x[!grepl("intensity", x, ignore.case = TRUE)], 
                collapse = "."))
        object@intensity <- as.data.frame(input[, input.i])
        colnames(object@intensity) <- make.names(input.n, unique = TRUE)
    } else if (any(grepl("signal", colnames(input), ignore.case = TRUE))) {
        input.i <- grep("signal", colnames(input), ignore.case = TRUE)
        input.n <- sapply(strsplit(colnames(input), "\\.", colnames(input))[input.i], 
            function(x) paste(x[!grepl("signal|methylated", x, ignore.case = TRUE)], 
                collapse = "."))
        if (!all(input.n[seq(1, length(input.i), 2)] == input.n[seq(2, 
            length(input.i), 2)])) 
            stop("names of both signal columns do not match.")
        object@intensity <- as.data.frame(input[, input.i[seq(1, length(input.i), 
            2)]] + input[, input.i[seq(2, length(input.i), 2)]])
        colnames(object@intensity) <- make.names(input.n[seq(1, length(input.i), 
            2)], unique = TRUE)

    } else {
        object@intensity <- as.data.frame(input)
    }
    if (!is.null(names)) 
        names(object) <- names
    
    object <- CNV.check(object)
    
    return(object)
})

#' @rdname CNV.load
setMethod("CNV.load", signature(input = "matrix", BAFsnps="data.frame"), function(input, BAFsnps, names = NULL) {
    CNV.load(as.data.frame(input), anno, names)
})

#' @rdname CNV.load
setMethod("CNV.load", signature(input = "numeric", BAFsnps="data.frame"), function(input, BAFsnps, names = NULL) {
    object <- new("CNV.data")
    
    if (is.null(names(input))) 
        stop("intensities not given for all probes.")
    object@intensity <- data.frame(sampleid = input)
    if (!is.null(names)) 
        names(object) <- names
    
    object <- CNV.check(object)
    
    return(object)
})


#' CNV.data.convert
#' @description Convert old CNV.data to new CNV.data
#' @param intensity
#' @param BAFsnps
#' @return \code{CNV.data} object.
#' @details This method converts old CNV.data to new CNV.data.
#' @author Damian Stichel \email{d.stichel@@dkfz.de}
#' @author Daniel Schrimpf \email{d.schrimpf@@dkfz.de}
#' @export
setGeneric("CNV.data.convert", function(intensity, BAFsnps, ...) {
    standardGeneric("CNV.data.convert")
})

#' @rdname CNV.data.convert
setMethod("CNV.data.convert", signature(intensity = "data.frame"), function(intensity) {
    object <- new("CNV.data")
    object@intensity<-intensity
    object@BAFsnps<-BAFsnps    

    object <- CNV.check(object)
    
    return(object)
})


#' CNV.check
#' @description Check intensity values.
#' @param object \code{CNV.data} object.
#' @return \code{CNV.data} object.
#' @details This method checks if intensities are positive and not NA. If not, they are set to 1. Warnings are given if intensities are abnormally high or low (> 50000 or < 5000, respectively).
#' @author Volker Hovestadt \email{conumee@@hovestadt.bio}
setGeneric("CNV.check", function(object) {
    standardGeneric("CNV.check")
})

#' @rdname CNV.check
setMethod("CNV.check", signature(object = "CNV.data"), function(object) {
    if (any(is.na(object@intensity))) {
        warning("some intensities are NA, now set to 1.")
        object@intensity[is.na(object@intensity)] <- 1
    }
    if (any(object@intensity < 0)) 
        warning("some intensities are smaller than 0, now set to 1.")
    object@intensity[object@intensity < 1] <- 1
    if (any(colMeans(object@intensity) < 5000)) 
        warning("intensities are abnormally low (< 5000).")
    if (any(colMeans(object@intensity) > 50000)) 
        warning("intensities are abnormally high (> 50000).")
    
    return(object)
})


#' read.450k.url
#' @description Read IDAT files from the web.
#' @param url URL of the directory in which the IDAT files are located.
#' @param idat Vector of IDAT names. \code{url} and \code{idat} default to the TCGA example described in the vignette.
#' @return \code{RGChannelSet} object.
#' @details This method downloads the provided list of IDAT files to a temporary folder (using the \code{RCurl} package). It then uses the `read.450k.exp` method of the `minfi` package.
#' @examples
#' RGsetTCGA <- read.450k.url()
#' @author Volker Hovestadt \email{conumee@@hovestadt.bio}
#' @export
read.450k.url <- function(url = NULL, idat = NULL) {
    if (is.null(url)) 
        url <- "https://tcga-data.nci.nih.gov/tcgafiles/ftp_auth/distro_ftpusers/anonymous/tumor/brca/cgcc/jhu-usc.edu/humanmethylation450/methylation/jhu-usc.edu_BRCA.HumanMethylation450.Level_1.8.8.0/"
    if (is.null(idat)) 
        idat <- c("6042324037_R05C02", "6042324037_R06C01")
    tmp <- paste0(tempdir(), .Platform$file.sep)
    if (!grepl("/$", url)) 
        url <- paste0(url, "/")
    if (!any(grepl("_Grn.idat", idat))) 
        idat <- unlist(lapply(idat, paste0, c("_Grn.idat", "_Red.idat")))
    for (i in idat) .curl(url = paste0(url, i), file = paste0(tmp, i))
    idatRG <- read.450k.exp(base = tmp)
    for (i in idat) file.remove(paste0(tmp, i))
    return(idatRG)
}

.curl <- function(url, file, verbose = TRUE) {
    if (.Platform$OS.type == "unix") {
        if (!RCurl::url.exists(url)) 
            stop("url does not exist.")
    } else {
        if (!RCurl::url.exists(url, .opts = list(ssl.verifypeer = FALSE))) 
            stop("url does not exist.")
    }
    if (verbose) 
        message("downloading ", tail(strsplit(url, "/")[[1]], 1), appendLF = FALSE)
    f <- RCurl::CFILE(file, mode = "wb")
    if (.Platform$OS.type == "unix") {
        r <- RCurl::curlPerform(url = url, writedata = f@ref, noprogress = TRUE)
    } else {
        r <- RCurl::curlPerform(url = url, writedata = f@ref, noprogress = TRUE, 
            .opts = list(ssl.verifypeer = FALSE))
    }
    RCurl::close(f)
    if (verbose) 
        message(" - done.")
} 
