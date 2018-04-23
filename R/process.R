##### PROCESSING methods #####

#' CNV.fit
#' @description Normalize query sample intensities by fitting intensities to reference set using a linear regression model.
#' @param query \code{CNV.data} object of query sample (single sample).
#' @param ref \code{CNV.data} object of reference set.
#' @param anno \code{CNV.anno} object. Use \code{CNV.create_anno} do create.
#' @param name character. Optional parameter to set query sample name.
#' @param intercept logical. Should intercept be considered? Defaults to \code{TRUE}.
#' @param ... Additional parameters (\code{CNV.fit} generic, currently not used).
#' @return \code{CNV.analysis} object.
#' @details The log2 ratio of query intensities versus a linear combination of reference set intensities that best reflects query intensities is calculated (as determined by linear regression). The annotations provided to \code{CNV.fit} are saved within the returned \code{CNV.analysis} object and used for subsequent analysis steps.
#' @examples
#' # prepare
#' library(minfiData)
#' data(MsetEx)
#' d <- CNV.load(MsetEx)
#' data(detail_regions)
#' anno <- CNV.create_anno(detail_regions = detail_regions)
#' 
#' # create object
#' x <- CNV.fit(query = d['GroupB_1'], ref = d[c('GroupA_1', 'GroupA_2', 'GroupA_3')], anno)
#' 
#' # modify object
#' #x <- CNV.bin(x)
#' #x <- CNV.detail(x)
#' #x <- CNV.segment(x)
#' 
#' # general information
#' x
#' show(x)
#' 
#' # coefficients of linear regression
#' coef(x)
#' 
#' # show or replace sample name
#' names(x)
#' names(x) <- 'Sample 1'
#' @author Volker Hovestadt \email{conumee@@hovestadt.bio}
#' @export
setGeneric("CNV.fit", function(query, ref, anno, ...) {
    standardGeneric("CNV.fit")
})

#' @rdname CNV.fit
setMethod("CNV.fit", signature(query = "CNV.data", ref = "CNV.data", anno = "CNV.anno"), 
    function(query, ref, anno, name = NULL, intercept = TRUE) {
        if (ncol(query@intensity) == 0) 
            stop("query intensities unavailable, run CNV.load")
        if (ncol(ref@intensity) == 0) 
            stop("reference set intensities unavailable, run CNV.load")
        
        if (ncol(query@intensity) != 1) 
            stop("query contains more than one sample.")
        if (ncol(ref@intensity) == 1) 
            warning("reference set contains only a single sample. use more samples for better results.")
        
        p <- names(anno@probes)  # ordered by location
        if (!all(is.element(p, rownames(query@intensity)))) 
            #stop("query intensities not given for all probes.")
            p<-p[-which(!is.element(p, rownames(query@intensity)))]
        if (!all(is.element(p, rownames(ref@intensity)))) 
            stop("reference set intensities not given for all probes.")
        
        object <- new("CNV.analysis")
        object@date <- date()
        object@fit$args <- list(intercept = intercept)
        
        if (!is.null(name)) {
            names(object) <- name
        } else {
            names(object) <- colnames(query@intensity)
        }
        object@anno <- anno
        
        if (max(cor(query@intensity[p, ], ref@intensity[p, ])[1, ]) > 0.99) 
           #stop("query sample seems to also be in the reference set. cannot fit against itself.")             #DS 14.09.16
         print("Warning: query sample seems to also be in the reference set.")
           ref@intensity<-ref@intensity[,-which.max(cor(query@intensity[p, ], ref@intensity[p, ])[1, ])]
        if (intercept) {
            ref.fit <- lm(y ~ ., data = data.frame(y = query@intensity[p, 
                1], X = ref@intensity[p, ]))
        } else {
            ref.fit <- lm(y ~ . - 1, data = data.frame(y = query@intensity[p, 
                1], X = ref@intensity[p, ]))
        }
        object@fit$coef <- ref.fit$coefficients
        
        ref.predict <- predict(ref.fit)
        ref.predict[ref.predict < 1] <- 1
        
        object@fit$ratio <- log2(query@intensity[p, 1]/ref.predict[p])
        object@fit$noise <- sqrt(mean((object@fit$ratio[-1] - object@fit$ratio[-length(object@fit$ratio)])^2, 
            na.rm = TRUE))
        
	object@BAFsnps<-query@BAFsnps
	object@gender<-query@gender

        return(object)
    })


#' CNV.bin
#' @description Combine single probe intensitiy values into predefined bins.
#' @param object \code{CNV.analysis} object.
#' @param ... Additional parameters (\code{CNV.bin} generic, currently not used).
#' @return \code{CNV.analysis} object.
#' @details The median intensity per bin is calculated. Bins are defined using \code{CNV.create_anno}. A value by which all probe and bin intensity values are shifted in subsequent analysis steps is calculated by minimizing the median absolute deviation from all bins to zero (ideally shifting the copy-neutral state to 0).
#' @examples
#' # prepare
#' library(minfiData)
#' data(MsetEx)
#' d <- CNV.load(MsetEx)
#' data(detail_regions)
#' anno <- CNV.create_anno(detail_regions = detail_regions)
#' 
#' # create object
#' x <- CNV.fit(query = d['GroupB_1'], ref = d[c('GroupA_1', 'GroupA_2', 'GroupA_3')], anno)
#' 
#' # modify object
#' x <- CNV.bin(x)
#' #x <- CNV.detail(x)
#' #x <- CNV.segment(x)
#' 
#' # general information
#' x
#' show(x)
#' 
#' # coefficients of linear regression
#' coef(x)
#' 
#' # show or replace sample name
#' names(x)
#' names(x) <- 'Sample 1'
#' @author Volker Hovestadt \email{conumee@@hovestadt.bio}
#' @export
setGeneric("CNV.bin", function(object, ...) {
    standardGeneric("CNV.bin")
})

#' @rdname CNV.bin
setMethod("CNV.bin", signature(object = "CNV.analysis"), function(object) {
    if (length(object@fit) == 0) 
        stop("fit unavailable, run CNV.fit")
    
    o1 <- as.matrix(findOverlaps(query = object@anno@bins, subject = object@anno@probes))
    o2 <- data.frame(bin = names(object@anno@bins)[o1[, "queryHits"]], 
        probe = names(object@anno@probes)[o1[, "subjectHits"]], stringsAsFactors = FALSE)
    
    object@bin$ratio <- sapply(split(object@fit$ratio[o2[, "probe"]], o2[, 
        "bin"]), median, na.rm = TRUE)[names(object@anno@bins)]
   # object@bin$shift <- optim(0, function(s) median(abs(object@bin$ratio - 
   #     s), na.rm = TRUE), method = "Brent", lower = -100, upper = 100)$par
   object@bin$shift <-0
    
    return(object)
})


#' CNV.detail
#' @description Combine single probe values within detail regions.
#' @param object \code{CNV.analysis} object.
#' @param ... Additional parameters (\code{CNV.detail} generic, currently not used).
#' @return \code{CNV.analysis} object.
#' @details The median intensity per detail region is calculated. Detail regions are defined using \code{CNV.create_anno(detail_bed=)}
#' @examples
#' # prepare
#' library(minfiData)
#' data(MsetEx)
#' d <- CNV.load(MsetEx)
#' data(detail_regions)
#' anno <- CNV.create_anno(detail_regions = detail_regions)
#' 
#' # create object
#' x <- CNV.fit(query = d['GroupB_1'], ref = d[c('GroupA_1', 'GroupA_2', 'GroupA_3')], anno)
#' 
#' # modify object
#' x <- CNV.bin(x)
#' x <- CNV.detail(x)
#' #x <- CNV.segment(x)
#' 
#' # general information
#' x
#' show(x)
#' 
#' # coefficients of linear regression
#' coef(x)
#' 
#' # show or replace sample name
#' names(x)
#' names(x) <- 'Sample 1'
#' @author Volker Hovestadt \email{conumee@@hovestadt.bio}
#' @export
setGeneric("CNV.detail", function(object, ...) {
    standardGeneric("CNV.detail")
})

#' @rdname CNV.detail
setMethod("CNV.detail", signature(object = "CNV.analysis"), function(object) {
    if (length(object@fit) == 0) 
        stop("fit unavailable, run CNV.fit")
    # if(length(object@bin) == 0) stop('bin unavailable, run CNV.bin')
    
    if (length(object@anno@detail) == 0) {
        message("no detail regions provided, define using CNV.create_anno")
    } else {
        d1 <- as.matrix(findOverlaps(query = object@anno@detail, subject = object@anno@probes))
        d2 <- data.frame(detail = values(object@anno@detail)$name[d1[, 
            "queryHits"]], probe = names(object@anno@probes[d1[, "subjectHits"]]), 
            stringsAsFactors = FALSE)
        
        object@detail$ratio <- sapply(split(object@fit$ratio[d2[, "probe"]], 
            d2[, "detail"]), median, na.rm = TRUE)[values(object@anno@detail)$name]
        object@detail$probes <- table(d2[, 1])[values(object@anno@detail)$name]
    }
    return(object)
})


#' @import DNAcopy
NULL

#' CNV.segment
#' @description Segment bin values (wrapper of \code{DNAcopy} package).
#' @param object \code{CNV.analysis} object.
#' @param alpha See details. Defaults to 0.001.
#' @param nperm See details. Defaults to 50000.
#' @param min.width See details. Defaults to 5.
#' @param undo.splits See details. Defaults to 'sdundo'.
#' @param undo.SD See details. Defaults to 2.2.
#' @param verbose See details. Defaults to 0.
#' @param ... Additional parameters supplied to the \code{segment} method of the \code{DNAcopy} package.
#' @return \code{CNV.analysis} object.
#' @details This method is a wrapper of the CNA, segment, segments.summary and segments.p methods of the DNAcopy package. Please refer to the respective man pages for more detailed information. The default parameters of \code{CNV.segment} override some of the default parameters of segment and are optimized for 450k data CNV analysis. 
#' @examples
#' # prepare
#' library(minfiData)
#' data(MsetEx)
#' d <- CNV.load(MsetEx)
#' data(detail_regions)
#' anno <- CNV.create_anno(detail_regions = detail_regions)
#' 
#' # create object
#' x <- CNV.fit(query = d['GroupB_1'], ref = d[c('GroupA_1', 'GroupA_2', 'GroupA_3')], anno)
#' 
#' # modify object
#' x <- CNV.bin(x)
#' x <- CNV.detail(x)
#' x <- CNV.segment(x)
#' 
#' # general information
#' x
#' show(x)
#' 
#' # coefficients of linear regression
#' coef(x)
#' 
#' # show or replace sample name
#' names(x)
#' names(x) <- 'Sample 1'
#' @author Volker Hovestadt \email{conumee@@hovestadt.bio}
#' @author Damian Stichel \email{d.stichel@@dkfz.de}
#' @export
setGeneric("CNV.segment", function(object, ...) {
    standardGeneric("CNV.segment")
})

setMethod("CNV.segment", signature(object = "CNV.analysis"), function(object, 
    alpha = 0.001, nperm = 50000, min.width = 5, undo.splits = "sdundo", 
    undo.SD = 2.2, verbose = 0, ...) {
    # if(length(object@fit) == 0) stop('fit unavailable, run CNV.fit')
    if (length(object@bin) == 0) 
        stop("bin unavailable, run CNV.bin")
    # if(length(object@detail) == 0) stop('bin unavailable, run
    # CNV.detail')
    
    a1 <- formals()
    a2 <- as.list(match.call())[-1]
    object@seg$args <- as.list(sapply(setdiff(unique(names(c(a1, a2))), 
        c("object", "verbose")), function(an) if (is.element(an, names(a2))) 
        a2[[an]] else a1[[an]], simplify = FALSE))
    
	    data(chromosome_positions)
    
    x1 <- DNAcopy::CNA(genomdat = object@bin$ratio[names(object@anno@bins)], 
        chrom = as.vector(seqnames(object@anno@bins)), maploc = values(object@anno@bins)$midpoint, 
        data.type = "logratio", sampleid = "sampleid")

		x1$chrom<-unlist(lapply(1:nrow(x1),function(i){ifelse(x1$maploc[i] < chromosome_positions$centromere[which(x1$chrom[i]==chromosome_positions$chromosome)],paste0(x1$chrom[i],"p"),paste0(x1$chrom[i],"q"))}))

		#x1$chrom[x1$chrom=="chr1"&x1$maploc>121535434]<-"chr1q"
		#x1$chrom[x1$chrom=="chr9"&x1$maploc<47367679]<-"chr9p"
		#x1$chrom[x1$chrom=="chr9"&x1$maploc>47367679]<-"chr9q"
		#x1$chrom[x1$chrom=="chr16"&x1$maploc<35335801]<-"chr16p"
		#x1$chrom[x1$chrom=="chr16"&x1$maploc>35335801]<-"chr16q"
    x2 <- DNAcopy::segment(x = x1, verbose = verbose, min.width = min.width, 
        nperm = nperm, alpha = alpha, undo.splits = undo.splits, undo.SD = undo.SD, 
        ...)
    object@seg$summary <- DNAcopy::segments.summary(x2)
    object@seg$summary$chrom <- as.vector(object@seg$summary$chrom)  # DNAcopy will factor chrom names. is there another way? 
    object@seg$p <- DNAcopy::segments.p(x2)
    object@seg$p$chrom <- as.vector(object@seg$p$chrom)
    
    object@seg$summary$chrom<-gsub("p","",object@seg$summary$chrom)
    object@seg$summary$chrom<-gsub("q","",object@seg$summary$chrom)
    
    return(object)
})



#' @author Damian Stichel \email{d.stichel@dkfz.de}
#' @author Daniel Schrimpf \email{d.schrimpf@dkfz.de}
#' @export
setGeneric("CNV.adjustbaseline", function(object, ...) {
  standardGeneric("CNV.adjustbaseline")
})
setMethod("CNV.adjustbaseline", signature(object = "CNV.analysis"), function(object, method.baseline.correction="BAF"){                                                                       
  
  #out<-NA
  #dnp.gr<-NA

  if( !exists("method.baseline.correction")){
    print("No method for baseline correction defined, using BAF (default). Options are 'MAD','MAXDENS','BAF'.")
     method.baseline.correction<-"BAF"
  }
  
  if( sum(is.na(object@BAFsnps$BAF))>50 & method.baseline.correction=="BAF"   ){
    print("Intensity for SNP positions not available. Using 'MAD' as method for baseline estimation.")
    method.baseline.correction<-"MAD"
  }

  if (!method.baseline.correction%in%c("MAD","BAF","MAXDENS")){
    print("Selected method not available, using BAF (default). Options are 'MAD','MAXDENS','BAF'.")
    method.baseline.correction<-"BAF"
  }
  
  object@method.baseline.correction<-method.baseline.correction

    if (method.baseline.correction=="MAD"){
	object@bin$shift<- optim(0, function(s) median(abs(object@bin$ratio - 
        s), na.rm = TRUE), method = "Brent", lower = -100, upper = 100)$par
    
  }else{
    x.histo <-hist(object@bin$ratio,breaks=seq(-20,20,by=0.01),xlim=c(-1,1),prob=TRUE,col="grey",plot=FALSE,warn.unused=FALSE)     #,plot=FALSE
    x.density <- density(object@bin$ratio,n=1024,bw=0.025)
    
    if (method.baseline.correction=="MAXDENS"){
	object@bin$shift<- x.density$x[which.max(x.density$y)]
    }else  if (method.baseline.correction=="BAF"){
	
      #Make it a time-series
      ts_y<-ts(x.density$y)
      #Compute turning points
      tp<-turnpoints(ts_y)
      
      out<-data.frame(pos=which(tp$peaks==TRUE | tp$pits==TRUE),tppos=tp$tppos,proba=tp$proba,info=tp$info)
      out$peak<-tp$peaks[out$pos]
      out$density<-x.density$y[out$tppos]
      #Keep only peaks, no pits
      out<-out[out$peak==TRUE,]
      #Keep only peaks with density higher than cutoff
      out<-out[out$density>=0.1,]
      
      dnp.gr<-GRanges(object@BAFsnps$chrom, IRanges(object@BAFsnps$start, object@BAFsnps$end), names=object@BAFsnps$probe,BAF=object@BAFsnps$BAF)
      
      seg.gr <- GRanges(object@seg$summary$chrom, IRanges(object@seg$summary$loc.start, object@seg$summary$loc.end), score=object@seg$summary$seg.median)
      df <-data.frame(xpos=sapply(as.list(findOverlaps(dnp.gr, seg.gr)), function(o) median(values(seg.gr)$score[o])),BAF=dnp.gr$BAF)
      df<-na.omit(df)
 
      #Associate rs-probes to intensity levels
      diff <- unlist(lapply(1:(length(x.density$x[out$tppos])-1),function(x) x.density$x[out$tppos][x+1]-x.density$x[out$tppos][x]))
      
      out$level <- c(1:length(x.density$x[out$tppos]))
      out$xpos  <- x.density$x[out$tppos]
      if (nrow(out)>5){out<-out[order(out$density,decreasing=TRUE),][1:5,]
      out<-out[order(out$xpos,decreasing=FALSE),] 
      out$level <-c(1:nrow(out))}
      
      border <- unlist(lapply(2:length(x.density$x[out$tppos]),function(x) x.density$x[out$tppos][x]-diff[x-1]/2))
      border <- c(out$xpos[1]-diff[1]/2,border,out$xpos[nrow(out)]+diff[length(diff)]/2)
      if (!is.na(border[1])){
        df$level<-unlist(lapply(1:nrow(df), function(x) {
          if( border[1]<df$xpos[x]&df$xpos[x]<border[length(border)] ){max(which(df$xpos[x] > border))}else{NA}
        }))
      } else df$level<-1
      df<-df[complete.cases(df),]
      
      ###### Get rid of levels, which have less than 6 SNPs assigned   ######
      levels.to.delete<-c()
      for (j in 1:nrow(out)){
        if (length(table(df$level)[names(table(df$level))==j])>0 && (table(df$level)[names(table(df$level))==j]>5) ){
          #print(j)
        }else {  
          levels.to.delete<-c(levels.to.delete,j)}
      }
      if(!is.null(levels.to.delete)&length(levels.to.delete)<nrow(out)){
        df$level[which(df$level%in%levels.to.delete)]<-NA
        out<-out[-levels.to.delete,]}
      
      ###Transform BAFs to be between 0 and 0.5
      df$BAF.transf[df$BAF<=0.5]<-df$BAF[df$BAF<=0.5]
      df$BAF.transf[df$BAF>0.5]<-abs(df$BAF[df$BAF>0.5]-1)
      
      #### Split with respect to intensity levels #####
      df.split<-split(df,df$level)
      
      ##### Compute BAF-Score #####  
      if (nrow(out)==1){out$score<-length(which(df$BAF>=0.4 & df$BAF<=0.6))/length(df$BAF)
      out$score2<-mean(abs(df$BAF[!is.na(df$BAF)]-0.5))
      out$score.abs<-length(which(df$BAF>=0.4 & df$BAF<=0.6))
      }else {out$score <- unlist(lapply(df.split,function(x){length(which(x$BAF>=0.4 & x$BAF<=0.6))/length(x$BAF)})) 
      out$score2 <- unlist(lapply(df.split,function(x){mean(abs(x$BAF[!is.na(x$BAF)]-0.5))}))
      out$score.abs<-unlist(lapply(df.split,function(x){length(which(x$BAF>=0.4 & x$BAF<=0.6))}))
      }
      
      ##### k-means Clustering ####  
      out$CMP<-unlist(lapply(df.split,function(x){
        clust<- kmeans(x$BAF.transf,centers=2)
        sort(clust$centers)[2]        #midpoint of upper cluster
      })) 
      out$level.candidate<-ifelse(out$CMP>0.4,2,ifelse(out$CMP<0.3,1,3))
      out$level.candidate[which(out$score.abs<3)]<-0            #  NEGLECT LEVELS WITH score.abs<2
      
      
      if ( length(which(out$level.candidate==2))==1 ) { shift<-out$xpos[which(out$level.candidate==2)]    
      } else if ( length(which(out$level.candidate==2))==0 ) {    ###if no candidate for level 2, use old method
        print("Warning: No candidate for level 2 detected")          
        if (out$score[1]>0.15&out$score.abs[1]>2){
          if(nrow(out)>1 &out$score.abs[1]<8 & out$score.abs[2]>8){
            shift<-out$xpos[2]
          } else {shift<-out$xpos[1]}
        }else{
          if(nrow(out)>2 &out$score.abs[2]<8 & out$score.abs[3]>8){
            shift<-out$xpos[3]
          } else{shift<-out$xpos[min(2,nrow(out))]}
        }  
      } else{
        print("Warning: More than 1 candidate for level 2 detected")
        shift<-min(out$xpos[which(out$level.candidate==2)])
      }
      
      object@bin$shift<-shift
        
	#object@bin$out<-out
	#object@bin$dnp.gr<-dnp.gr
    }else print("ERROR: Selected method is not available.")
    
    object@bin$ratio<-log2((2^object@bin$ratio)+(1-2^object@bin$shift))
    object@fit$ratio<-log2((2^object@fit$ratio)+(1-2^object@bin$shift))
    object@bin$ratio[is.nan(object@bin$ratio)]<- -1.2
    object@fit$ratio[is.nan(object@fit$ratio)]<- -1.2
    object@detail$ratio<-log2((2^object@detail$ratio)+(1-2^object@bin$shift))
    object@detail$ratio[is.nan(object@detail$ratio)]<- -1.2

    object@seg$summary$seg.mean<-log2((2^object@seg$summary$seg.mean)+(1-2^object@bin$shift))
    object@seg$summary$seg.sd<-log2((2^object@seg$summary$seg.sd)+(1-2^object@bin$shift))
    object@seg$summary$seg.median<-log2((2^object@seg$summary$seg.median)+(1-2^object@bin$shift))
    object@seg$summary$seg.mad<-log2((2^object@seg$summary$seg.mad)+(1-2^object@bin$shift))
    object@seg$summary$seg.mean[is.nan(object@seg$summary$seg.mean)]<- -1.2
    object@seg$summary$seg.sd[is.nan(object@seg$summary$seg.sd)]<- -1.2    
    object@seg$summary$seg.median[is.nan(object@seg$summary$seg.median)]<- -1.2
    object@seg$summary$seg.mad[is.nan(object@seg$summary$seg.mad)]<- -1.2
    
    }
  return(object)
})



#' @author Damian Stichel \email{d.stichel@dkfz.de}
#' @author Daniel Schrimpf \email{d.schrimpf@dkfz.de}
#' @export
setGeneric("CNV.evaluation", function(object, ...) {
  standardGeneric("CNV.evaluation")
})
setMethod("CNV.evaluation", signature(object = "CNV.analysis"), function(object){                                                                       
  
  #Read information, which cpg belongs to which chromosome arm
  data(cpg_chromosome_arms)
  cpg_chromosome_armsM<-cpg_chromosome_arms
  cpg_chromosome_armsF<-cpg_chromosome_arms[-which(cpg_chromosome_arms$Chromosome=="chrY"),]
  
  #Read information for chromosome_positions
  data(chromosome_positions)
  rownames(chromosome_positions)<-c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX","chrY")
  #chromosome_positions$chromosome <-rownames(chromosome_positions)
  chromosome_positionsM<-chromosome_positions
  chromosome_positionsF<-chromosome_positions[1:(nrow(chromosome_positions)-1),]  #without Y-chrom
  
      
      if(object@gender == "F") {
        chromosome_positions<-chromosome_positionsF       #get rid of Y for females, keep it for m and nd
        cpg_chromosome_arms<-cpg_chromosome_armsF
      }else {
        chromosome_positions<-chromosome_positionsM 
        cpg_chromosome_arms<-cpg_chromosome_armsM
      }
      
        #Data for probes
        if (length(object@fit) == 0) 
          stop("fit unavailable, run CNV.fit")
          # data <- data.frame(Chromosome = as.vector(seqnames(object@anno@probes)), 
           #             Start = start(object@anno@probes) - 1, 
            #            End = end(object@anno@probes), 
             #           Feature = names(object@anno@probes),
              #          Value = round(object@fit$ratio, 3), row.names = NULL)
           data <- data.frame(Chromosome = as.vector(seqnames(object@anno@probes)[names(object@anno@probes)%in%names(object@fit$ratio)]), 
                              Start = start(object@anno@probes)[names(object@anno@probes)%in%names(object@fit$ratio)] - 1, 
                              End = end(object@anno@probes)[names(object@anno@probes)%in%names(object@fit$ratio)], 
                              Feature = names(object@anno@probes)[names(object@anno@probes)%in%names(object@fit$ratio)],
                              Value = round(object@fit$ratio, 3), row.names = NULL)
        #colnames(data) <- sub("Value", object@name, colnames(data))
  
        colnames(data)<- sub("Value", "Intensity", colnames(data))
        data <- merge(data,cpg_chromosome_arms[,4:5],by="Feature")
        
        #Split with respect to chromosome arms
        data.pq <-split(data,data$pq)

        #Compute mean Intensity, sd, median                                      #oder log2(mean(as.numeric(2^x$Intensity),na.rm=TRUE))
        out <-data.frame(Chromosome=names(data.pq),meanIntensity=unlist(lapply(data.pq,function(x){mean(as.numeric(x$Intensity),na.rm=TRUE)})) )
        out$sd <-unlist(lapply(data.pq,function(x){sd(2^as.numeric(x$Intensity),na.rm=TRUE)}))
        out$medianIntensity <-unlist(lapply(data.pq,function(x){median(as.numeric(x$Intensity),na.rm=TRUE)}))

        #max(as.numeric(data.pq$chr2p$Intensity))
        #out$medianIntensity <- out$medianIntensity-1

        #Sort out
        tmp <-gsub("chr","",out$Chromosome)
        tmp <-gsub("p","",tmp)
        tmp <-gsub("q","",tmp)
        tmp <- as.numeric(tmp)
        out <- out[order(tmp),]

        #Compute, if gain/loss/balanced
        out$alteration <- ifelse(out$medianIntensity >=0.1,"gain",ifelse(out$medianIntensity <=-0.1,"loss","balanced"))
        #out$type <- ifelse(out$medianIntensity <=-0.6,"homo-del",ifelse(out$medianIntensity <=-0.1,"loss",ifelse(out$medianIntensity >=0.4,"amp",ifelse(out$medianIntensity >=0.1,"gain","balanced"))))

        object@arms$summary<-out

        #GENES OF INTEREST
        if (length(object@detail) == 0)
          stop("detail unavailable, run CNV.bin")
           #     det <- data.frame(chr = as.vector(seqnames(object@anno@detail)),
            #            start = start(object@anno@detail), end = end(object@anno@detail),
             #           name = names(object@detail$probes), sample = object@name, probes = object@detail$probes,
              #          value = round(object@detail$ratio, 3), row.names = NULL)

        object@detail$alteration <- ifelse(object@detail$ratio >=0.1,"gain",ifelse(object@detail$ratio <=-0.1,"loss","balanced"))
                                        # <- ifelse(det$value <=-0.4,"homo-del",ifelse(det$value <=-0.1,"loss",ifelse(det$value >=0.4,"amp",ifelse(det$value >=0.1,"gain","balanced"))))
       
        
        
        
              
        #############  Evaluate alterations for segments ##############
        
        #countsegments <- as.data.frame(table(segments$chrom))
        #rownames(countsegments)<-countsegments$Var1
        
        # test <- apply(data[1:1000,],1,function(x){
        #   tmp <- segments[grepl(paste("^",x["Chromosome"],"$", sep=""),segments$chrom),]
        #   select <-  tmp[tmp$loc.start<=as.numeric(x["Start"]) & tmp$loc.end>=as.numeric(x["Start"]), ]
        #   paste(select$chrom, ".", select$loc.start, ".", select$loc.end, sep="")
        # })
        # data$segment <- test
        
        if (length(object@seg) == 0) 
          stop("seg unavailable, run CNV.bin")
        # seg format, last numeric column is used in igv
        segments <- data.frame(ID = object@name, chrom = object@seg$summary$chrom, 
                        loc.start = object@seg$summary$loc.start, loc.end = object@seg$summary$loc.end, 
                        num.mark = object@seg$summary$num.mark, bstat = object@seg$p$bstat, 
                        pval = object@seg$p$pval, seg.mean = round(object@seg$summary$seg.mean
                                 , 3), seg.median = round(object@seg$summary$seg.median, 3), row.names = NULL)
        segments$segname <-paste(segments$chrom, ".", segments$loc.start, ".", segments$loc.end,sep="")
        
        object@seg$summary$alteration<- ifelse(object@seg$summary$seg.median >=0.1,"gain",ifelse(object@seg$summary$seg.median <=-0.1,"loss","balanced"))
          # ... <- ifelse(out2$medianIntensity <=-0.4,"homo-del",ifelse(out2$medianIntensity <=-0.1,"loss",ifelse(out2$medianIntensity >=0.4,"amp",ifelse(out2$medianIntensity >=0.1,"gain","balanced"))))
        
        object@seg$summary$chrom.arm<-  ifelse({as.numeric(object@seg$summary$loc.start)<as.numeric(chromosome_positions[object@seg$summary$chrom,2])&as.numeric(object@seg$summary$loc.end)<as.numeric(chromosome_positions[object@seg$summary$chrom,2])}, paste(object@seg$summary$chrom,"p",sep=""),    ifelse({as.numeric(object@seg$summary$loc.start)>as.numeric(chromosome_positions[object@seg$summary$chrom,2])&as.numeric(object@seg$summary$loc.end)>as.numeric(chromosome_positions[object@seg$summary$chrom,2])},paste(object@seg$summary$chrom,"q",sep=""), object@seg$summary$chrom))
        object@seg$summary$segment.length <- (as.numeric(object@seg$summary$loc.end)-as.numeric(object@seg$summary$loc.start))   #/1000000

        #Keep Y only for male:
        if (!object@gender=="M"&length(which(object@seg$summary$chrom %in% c("chrYp","chrYq","chrY")))>0){object@seg$summary<-object@seg$summary[-which(object@seg$summary$chrom %in% c("chrYp","chrYq","chrY")),]
                                                                                                          object@seg$p<-object@seg$p[-which(object@seg$p$chrom %in% c("chrYp","chrYq","chrY")),]}
        
         #CHROMOTHRIPSIS:
        if (length(which(table(object@seg$summary$chrom.arm)>9))>0){
          vec<-which(table(object@seg$summary$chrom.arm)>9)
          for (i in 1:length(vec)){
            newrow<-object@seg$summary[which(object@seg$summar$chrom.arm==names(vec)[i])[1],]
            newrow$loc.start<- NaN
            newrow$loc.end<- NaN
            newrow$num.mark<-NaN
            newrow$seg.mean<-NaN
            newrow$seg.median<-NaN
            newrow$seg.mad<-NaN
            newrow$alteration<-"chromothripsis"
            newrow$segment.length<-NaN
            object@seg$summary<-rbind(object@seg$summary,newrow)        
            object@seg$p<-rbind(object@seg$p,rep(NA,10))  
            #if Chromothripsis is found, change result for chr. arms:
            object@arms$summary$alteration[which(object@arms$summary$Chromosome==names(vec[i]))]<-"chromothripsis"   
          }
        }
        
         #summary of total number of alteration
         object@alteration.summary<-data.frame(object@name,alt_arm_tot=length(which(object@arms$summary$alteration=="loss"))+length(which(object@arms$summary$alteration=="gain")),alt_arm_gain=length(which(object@arms$summary$alteration=="gain")),alt_arm_loss=length(which(object@arms$summary$alteration=="loss")),
                         alt_seg_tot=length(which(object@seg$summary$alteration=="loss"))+length(which(object@seg$summary$alteration=="gain")),alt_seg_gain=length(which(object@seg$summary$alteration=="gain")),alt_seg_loss=length(which(object@seg$summary$alteration=="loss")),
                         alt_leng_tot=sum(object@seg$summary$segment.length[object@seg$summary$alteration=="gain"])+sum(object@seg$summary$segment.length[object@seg$summary$alteration=="loss"]),     #+sum(object@seg$summary$alteration$Segment_length[object@seg$summary$alteration$alteration=="amp"])+sum(object@seg$summary$alteration$Segment_length[object@seg$summary$alteration$alteration=="homo-del"]),
                         alt_leng_gain=sum(object@seg$summary$segment.length[object@seg$summary$alteration=="gain"]),#+sum(object@seg$summary$segment.length[object@seg$summary$alteration=="amp"]),
                         alt_leng_loss=sum(object@seg$summary$segment.length[object@seg$summary$alteration=="loss"])#+sum(object@seg$summary$segment.length[object@seg$summary$alteration=="homo-del"])
                         )
    
         seg.pq<-split(object@seg$summary,object@seg$summary$chrom.arm)
         object@arms$summary <-data.frame(object@name,Chromosome=names(seg.pq),alteration=unlist(lapply(seg.pq,function(y){
           if (length(unique(y$alteration))>1 ){
             #alt<-sort(unique(y$alteration))
             perc<-c(sum(y$segment.length[y$alteration=="balanced"])/sum(y$segment.length),
                      sum(y$segment.length[y$alteration=="gain"])/sum(y$segment.length),
                      sum(y$segment.length[y$alteration=="loss"])/sum(y$segment.length))
             perc<-round(perc*100,digits=0)
             #paste0("Segmental.changes:.",paste(sort(unique(y$alteration)),collapse="_"))
             paste0("Segmental.changes:","bal(", perc[1],"%)_gain(",perc[2],"%)_loss(",perc[3],"%)")
             }
           else{paste(sort(unique(y$alteration)),collapse="_")}
           #mean(as.numeric(x$Intensity),na.rm=TRUE)
           })) )
         
         
         # #Compute mean Intensity, sd, median                                      #oder log2(mean(as.numeric(2^x$Intensity),na.rm=TRUE))
         # out <-data.frame(Chromosome=names(data.pq),meanIntensity=unlist(lapply(data.pq,function(x){mean(as.numeric(x$Intensity),na.rm=TRUE)})) )
         # out$sd <-unlist(lapply(data.pq,function(x){sd(2^as.numeric(x$Intensity),na.rm=TRUE)}))
         # out$medianIntensity <-unlist(lapply(data.pq,function(x){median(as.numeric(x$Intensity),na.rm=TRUE)}))
         # 
         # #max(as.numeric(data.pq$chr2p$Intensity))
         # #out$medianIntensity <- out$medianIntensity-1
         # 
         # #Sort out
         # tmp <-gsub("chr","",out$Chromosome)
         # tmp <-gsub("p","",tmp)
         # tmp <-gsub("q","",tmp)
         # tmp <- as.numeric(tmp)
         # out <- out[order(tmp),]
         # 
         # #Compute, if gain/loss/balanced
         # out$alteration <- ifelse(out$medianIntensity >=0.1,"gain",ifelse(out$medianIntensity <=-0.1,"loss","balanced"))
         # #out$type <- ifelse(out$medianIntensity <=-0.6,"homo-del",ifelse(out$medianIntensity <=-0.1,"loss",ifelse(out$medianIntensity >=0.4,"amp",ifelse(out$medianIntensity >=0.1,"gain","balanced"))))
         # 
         # object@arms$summary<-out
         # 
         # #GENES OF INTEREST
         # if (length(object@detail) == 0) 
         #   stop("detail unavailable, run CNV.bin")
         #    #     det <- data.frame(chr = as.vector(seqnames(object@anno@detail)), 
         #     #            start = start(object@anno@detail), end = end(object@anno@detail), 
         #      #           name = names(object@detail$probes), sample = object@name, probes = object@detail$probes, 
         #       #          value = round(object@detail$ratio, 3), row.names = NULL)
         #         
         
         
         
return(object)
})




