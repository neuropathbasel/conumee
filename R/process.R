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
            stop("query intensities not given for all probes.")
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
#' @export
setGeneric("CNV.segment", function(object, ...) {
    standardGeneric("CNV.segment")
})

#' @rdname CNV.segment
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
    
    x1 <- DNAcopy::CNA(genomdat = object@bin$ratio[names(object@anno@bins)], 
        chrom = as.vector(seqnames(object@anno@bins)), maploc = values(object@anno@bins)$midpoint, 
        data.type = "logratio", sampleid = "sampleid")
    x2 <- DNAcopy::segment(x = x1, verbose = verbose, min.width = min.width, 
        nperm = nperm, alpha = alpha, undo.splits = undo.splits, undo.SD = undo.SD, 
        ...)
    object@seg$summary <- DNAcopy::segments.summary(x2)
    object@seg$summary$chrom <- as.vector(object@seg$summary$chrom)  # DNAcopy will factor chrom names. is there another way? 
    object@seg$p <- DNAcopy::segments.p(x2)
    object@seg$p$chrom <- as.vector(object@seg$p$chrom)
    
    return(object)
}) 


setGeneric("CNV.segment_DS", function(object, ...) {
    standardGeneric("CNV.segment_DS")
})

setMethod("CNV.segment_DS", signature(object = "CNV.analysis"), function(object, 
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
    

		chromosome_positions<-read.xlsx("/home/damian/sequencing/Projekte/conumee/data/chromosome_positions.xlsx")
		chromosome_positions$chromosome<-paste0("chr",chromosome_positions$chromosome)

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
    
    return(object)
})




setGeneric("CNV.adjustbaseline", function(object, method.baseline.correction, ...) {
  standardGeneric("CNV.adjustbaseline")
})
setMethod("CNV.adjustbaseline", signature(object = "CNV.analysis"), function(object, method.baseline.correction="BAF"){                                                                       
  
  if (!method.baseline.correction%in%c("MAD","BAF","MAXDENS")){
    print("Selected method not available, using BAF (default). Options are 'MAD','MAXDENS','BAF'.")
    method.baseline.correction<-"BAF"
  }
  
  object@method.baseline.correction<-method.baseline.correction

    if (method.baseline.correction=="MAD"){
	object@bin$shift<- optim(0, function(s) median(abs(object@bin$ratio - 
        s), na.rm = TRUE), method = "Brent", lower = -100, upper = 100)$par
    
  }else{
    x.histo <-hist(object@bin$ratio,breaks=seq(-20,20,by=0.01),xlim=c(-1,1),prob=TRUE,col="grey")     #,plot=FALSE
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
      
      ### Compute BAF (only 65 probes)
      #BAFs=data.frame(V4=rownames(betas)[which(rownames(betas) %in% pos$V4)],BAF=betas[which(rownames(betas) %in% pos$V4)])
      #BAFs<-as.data.frame(getSnpBeta(RGset))


      dnp.gr<-GRanges(object@BAFsnps$chrom, IRanges(object@BAFsnps$start, object@BAFsnps$end), names=object@BAFsnps$probe,BAF=object@BAFsnps$BAF)
      
      #### Merge BAFs to dnp.df, make sure to keep correct order   
     # dnp.df<-object@BAFsnps

      #for (ind in 1:nrow(dnp.df)){dnp.df$order[ind]<-ind}
      #dnp.df<-merge(dnp.df,BAFs,by="V4",all.x=TRUE,sort=FALSE)                 
      #dnp.df<-dnp.df[order(dnp.df$order),]
      
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
          print(j)
        }else {  
          print(paste0("Bad level: ", j))
          levels.to.delete<-c(levels.to.delete,j)}
      }
      if(!is.null(levels.to.delete)&length(levels.to.delete)<nrow(out)){
        df$level[which(df$level%in%levels.to.delete)]<-NA
        out<-out[-levels.to.delete,]}
      
      #if(isEmpty(which(df$level==1))) {df<-rbind(df,c(NA,NA,NA,1,1))}
      #if(max(df$level)>=2 & isEmpty(which(df$level==2))) {df<-rbind(df,c(NA,NA,NA,2,2))}
      #if(max(df$level)>=3 & isEmpty(which(df$level==3))) {df<-rbind(df,c(NA,NA,NA,3,3))}
      #if(max(df$level)>=4 & isEmpty(which(df$level==4))) {df<-rbind(df,c(NA,NA,NA,4,4))}
      
      ###Transform BAFs to be between 0 and 0.5
      df$BAF.transf[df$BAF<0.5]<-df$BAF[df$BAF<0.5]
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
        sort(clust$centers)[2]        #ClusterMittelPunkt von "oberem" Cluster ausgeben
      })) 
      out$level.candidate<-ifelse(out$CMP>0.4,2,ifelse(out$CMP<0.3,1,3))
      out$level.candidate[which(out$score.abs<3)]<-0            #  NEGLECT LEVELS WITH score.abs<2
      
      
      if ( length(which(out$level.candidate==2))==1 ) { shift_v3<-out$xpos[which(out$level.candidate==2)]    
      } else if ( length(which(out$level.candidate==2))==0 ) {    ###wenn kein Candidat für Level 2, dann nach alter Methode
        warning<-paste0(warning,"no candidate for level 2 detected")          
        if (out$score[1]>0.15&out$score.abs[1]>2){
          if(nrow(out)>1 &out$score.abs[1]<8 & out$score.abs[2]>8){
            shift_v3<-out$xpos[2]
          } else {shift_v3<-out$xpos[1]}
        }else{
          if(nrow(out)>2 &out$score.abs[2]<8 & out$score.abs[3]>8){
            shift_v3<-out$xpos[3]
          } else{shift_v3<-out$xpos[min(2,nrow(out))]}
        }  
      } else{
        warning<-paste0(warning,"more than 1 candidate for level 2 detected")
        shift_v3<-min(out$xpos[which(out$level.candidate==2)])
      }
      
      object@bin$shift<-shift_v3
        
	object@bin$out<-out
	object@bin$dnp.df<-dnp.df
	object@bin$dnp.gr<-dnp.gr
    }else print("ERROR: Selected method is not available.")
    }
  return(object)
})

