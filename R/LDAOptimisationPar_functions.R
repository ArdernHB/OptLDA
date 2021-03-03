




#' LDA correct cross-validation with equal sample size using parallel processing
#'
#' This function takes a matrix of variables of known group membership and returns
#' the results of a leave-one-out correct cross validation identification for each
#' specimen to provide a correct cross-validation percentage. This function is for
#' application to a set number of variables (or PCs), see the function
#' \code{LDACVStepwisePar} for stepwise testing of each set of
#' consecutive PCs.
#'
#' The function is primarily for use with resampling unequal groups to equal sample
#' size a set number of times. This process is carried out with parrallel processing.
#' If shape data is used the function offers the option to carry out a new GPA and
#' subsequent PCA for each resampling exercise (set ShapeGPA=TRUE). If raw data
#' is used the function carries out a new PCA with each resampling
#' exercise. In both these instances where a fresh PCA is carried
#' out the function will call on the input value in PClim to determine the number
#' of PCs to use.
#'
#' Note that the method in this function is to remove data leakage from unequal sample
#' size, but that the leave-one-out method here does not train the identification model
#' separately for each specimen. Instead it is assumed that the single specimen being
#' left out in the cross validation procedure contributes so little data leakage as to
#' be unimportant. However, if this is a concern please use the k-fold method in
#' \code{LDACVManyPar} and \code{LDACVManyStepwisePar} and set k to the sample size of
#' the smallest sample size.
#'
#' @param DiscriminationData A matrix of numeric data or an array for shape data. If shape data is used please set ShapeGPA=TRUE. A shape data array is expected to include landmarks as rows, landmark dimensions (either 2 or 3) as columns and specimens as slices. If a matrix is used, please ensure no columns have erroneously been loaded as factors.
#' @param GroupMembership A vector of group classifications, either as factors or as characters.
#' @param EqualIter an integer determining the number of times the data will be resampled to equal sample size. Default is set to 100.
#' @param SampleSize the sample size to be used. If set to NA the resampling exercise will automatically resample groups to the sample size of the smallest group.
#' @param Verbose logical (either TRUE or FALSE, default set to FALSE) to determine whether to return the identifications of each round of resampling (Verbose=TRUE) or to just return the correct-cross-validation percentages of each resampling exercise (Verbose=FALSE).
#' @param ShapeGPA logical (either TRUE or FALSE, default set to FALSE) to indicate whether the data is shape data and if so each subset of individuals from the resampling procedure will be processed through a Generalised Procrustes Analysis (GPA) and subsequent Principal Component Analysis (PCA) using the \code{\link[Morpho]{procSym}} function in the package \code{Morpho}.
#' @param Sliding if ShapeGPA is TRUE and the shape data has sliding landmarks then a list should be provided here with vectors of landmarks to be slid (see function description of \code{\link[Morpho]{procSym}}).
#' @param SlidingLMindex a vector of the landmark numbers that represent sliding landmarks.
#' @param SizeShape logical (either TRUE or FALSE, default set to FALSE). If ShapeGPA is TRUE and you wish to analyses form (i.e. shape+log size) then this argument will be passed to the \code{\link[Morpho]{procSym}} and the balanced LDA will be applied to form PCs (see function description of \code{\link[Morpho]{procSym}}).
#' @param PClim integer determining the number of PCs to be used in the case that ShapeGPA or PCA are set to TRUE. Default is arbitrarily set to 10.
#' @return Returns a matrix of the leave-one-out classifications for all the specimens along with their known classification.
#'
#'
#' @author Ardern Hulme-Beaman
#' @import MASS
#' @import Morpho
#' @import doParallel
#' @import parallel
#' @import foreach
#' @export



LDACVPar <- function(DiscriminationData, GroupMembership, EqualIter=100, SampleSize=NA, Verbose=FALSE, ShapeGPA=FALSE, Sliding=NULL, SlidingLMindex=NULL, SizeShape=FALSE, PClim=10){
  #DiscriminationData=BlackRatGPA$PCscores
  #GroupMembership=Groups
  #EqualIter=100
  #PClim=3
  #ShapeGPA=TRUE; Verbose=FALSE; Sliding=NA;
  #SizeShape=FALSE; PCA=FALSE

  #DiscriminationData = M1GPA$PCscores[RemPos,1:PCdim]
  #GroupMembership = Grouping[RemPos]
  #ShapeGPA = FALSE
  #PClim = PCdim

  if (ShapeGPA==TRUE & length(dim(DiscriminationData))==2){
    stop('DiscriminationData is a matrix, but ShapeGPA set to TRUE, if the data is shape data then it must be in array format, with LMs as rows, dimensions as columns and specimens as slices')
  }


  if (ShapeGPA==TRUE){
    if (dim(DiscriminationData)[3]!=length(GroupMembership)){
      stop('Number of specimens in DiscriminationData does not appear to match the number of speicmens listed in GroupMembership')
    }
  } else {
    if (dim(DiscriminationData)[1]!=length(GroupMembership)){
      stop('Number of specimens in DiscriminationData does not appear to match the number of speicmens listed in GroupMembership')
    }
  }



  chr2nu <- function(X){
    as.numeric(as.character(X))
  }

  ParOutput <- function(PreviousResults, ResultList){

    NewResults <- PreviousResults
    if (length(PreviousResults)==2){
      NewResults$Global.CVP <- rbind(PreviousResults$Global.CVP, ResultList$Global.CVP)
      NewResults$Grp.CVP <- rbind(PreviousResults$Grp.CVP, ResultList$Grp.CVP)

    } else {
      NewResults$Global.CVP <- rbind(PreviousResults$Global.CVP, ResultList$Global.CVP)
      NewResults$Grp.CVP <- cbind(PreviousResults$Grp.CVP, ResultList$Grp.CVP)
      NewResults$CV.IDs <- cbind(PreviousResults$CV.IDs, ResultList$CV.IDs)
      NewResults$True.IDs <- NewResults$True.IDs

    }

    return(NewResults)
  }

  BalancedGrps <- function(GroupMembership, GroupSize=SampleSize){
    #DiscriminationData=PairwiseShapePCAmat; GroupMembership=chr(Groups[GrpPos])

    if (is.na(GroupSize)){
      minSampSize <- min(table(as.character(GroupMembership)))
    } else {
      minSampSize <- GroupSize
    }

    sampleindex <- 1:length(GroupMembership)
    RandSamp <- stats::aggregate(x = sampleindex, by=list(factor(GroupMembership)), sample, size=minSampSize)
    Index <- c(t(RandSamp$x))
    GroupMem <- c(t(stats::aggregate(RandSamp$Group.1, list(RandSamp$Group.1), rep, minSampSize)$x))
    return(list(IndexedLocations=Index, Newfactors=GroupMem))

  }


  ParEqualIter <- function(DiscriminationData, GrpMem, ShapeGPA, Sliding, SlidingLMindex, ParVerbose=Verbose, SizeShape, PClim, SampleSize){
    #DiscriminationData=DiscriminationData; GrpMem=GroupMembership; ParTieBreaker='Report'; ParVerbose=FALSE
    #GrpMem=Groups; PClim=3; SampleSize=NA
    BalancingGrps <- BalancedGrps(GrpMem, SampleSize)


    if (ShapeGPA==TRUE){
      BalDataShape <- DiscriminationData[,,BalancingGrps$IndexedLocations]
      invisible(utils::capture.output(BalData <- Morpho::procSym(BalDataShape, sizeshape = SizeShape, outlines = Sliding, SMvector = SlidingLMindex)$PCscores))
      LDAres <- MASS::lda(x = BalData[,1:PClim], grouping=BalancingGrps$Newfactors, CV=TRUE)
    } else {
      BalData <- stats::prcomp(x = DiscriminationData[BalancingGrps$IndexedLocations,])
      LDAres <- MASS::lda(x = BalData$x[,1:PClim], grouping=BalancingGrps$Newfactors, CV=TRUE)
    }

    CVP <- sum(LDAres$class==BalancingGrps$Newfactors)/length(BalancingGrps$Newfactors)
    CVPbyGrp <- diag(table(LDAres$class, BalancingGrps$Newfactors))/table(BalancingGrps$Newfactors)

    if (ParVerbose==TRUE){
      Res <- list(Global.CVP=CVP, Grp.CVP=CVPbyGrp, CV.IDs=as.character(LDAres$class), True.IDs=as.character(BalancingGrps$Newfactors))
      return(Res)
    } else {
      return(list(Global.CVP=CVP, Grp.CVP=CVPbyGrp))
    }
  }


  cores <- parallel::detectCores()
  clust <- parallel::makeCluster(cores[1]-1)
  doParallel::registerDoParallel(clust)

  a <- 1
  ParResults <- foreach::foreach(a = 1:EqualIter, .combine = ParOutput) %dopar% {
    ParEqualIter(DiscriminationData, GrpMem=GroupMembership, ShapeGPA=ShapeGPA, Sliding=Sliding, SlidingLMindex=SlidingLMindex, ParVerbose=Verbose, PClim=PClim, SizeShape = SizeShape, SampleSize=SampleSize)
  }

  parallel::stopCluster(clust)


  ParResults$Grp.CVP

  if (Verbose==TRUE){
    names(ParResults) <- c('Iteration.Summaries', 'Verbose.Output')
    graphics::hist(ParResults$Iteration.Summaries)
    return(ParResults, main= paste('Correct cross-validation % for', PClim, 'PCs', sep = ' '), xlab = 'CCV%')
  } else {
    graphics::hist(ParResults[[1]], main= paste('Correct cross-validation % for', PClim, 'PCs', sep = ' '), xlab = 'CCV%')
    return(ParResults[[1]])
  }

}


#' Stepwise LDA correct cross-validation with equal sample size using parallel processing
#'
#' This function takes a matrix of variables of known group membership and returns
#' the results of a leave-one-out correct cross validation identification for each
#' specimen to provide a correct cross-validation percentage. This function is for
#' for stepwise testing of each set of consecutive PCs with resampling to equal sample
#' size at each incremental increase.
#'
#' The function is primarily for use with resampling unequal groups to equal sample
#' size a set number of times. This process is carried out with parallel processing.
#' If shape data is used the function offers the option to carry out a new GPA and
#' subsequent PCA for each resampling exercise (set ShapeGPA=TRUE). If raw data
#' is used the function carries out a new PCA with each resampling
#' exercise. In both these instances where a fresh PCA is carried
#' out the function will call on the input value in PClim to determine the number
#' of PCs to use.
#'
#' Note that the method in this function is to remove data leakage from unequal sample
#' size, but that the leave-one-out method here does not train the identification model
#' separately for each specimen. Instead it is assumed that the single specimen being
#' left out in the cross validation procedure contributes so little data leakage as to
#' be unimportant. However, if this is a concern please use the k-fold method in
#' \code{LDACVManyPar} and \code{LDACVManyStepwisePar} and set k to the sample size of
#' the smallest sample size.
#'
#' @param PlotResults logical (either TRUE or FALSE, default set to TRUE) to determine whether to plot the results of the stepwise discriminant analyses.
#' @param CombinePlots logical set to FALSE to indicate whether the results of this analysis should be plotted on a previous result that's currently open (e.g. for randomised analyses).
#' @inheritParams LDACVPar
#' @return Returns a matrix of the leave-one-out classifications for all the specimens along with their known classification.
#'
#'
#' @author Ardern Hulme-Beaman
#' @import MASS
#' @import Morpho
#' @import doParallel
#' @import parallel
#' @import foreach
#' @export


LDACVStepwisePar <- function(DiscriminationData, GroupMembership, EqualIter=100, SampleSize=NA, Verbose=FALSE, ShapeGPA=FALSE, Sliding=NULL, SlidingLMindex=NULL, SizeShape=FALSE, PClim=10, PlotResults=TRUE, CombinePlots=FALSE){

  #DiscriminationData = LMDataArray; GroupMembership = Groups; EqualIter=1000; SampleSize=NA; Verbose=FALSE; ShapeGPA=TRUE; Sliding=NA; SizeShape=FALSE; PClim=20; PlotResults=TRUE; CombinePlots=FALSE
  #DiscriminationData = Shape


  if (ShapeGPA==TRUE & length(dim(DiscriminationData))==2){
    stop('DiscriminationData is a matrix, but ShapeGPA set to TRUE, if the data is shape data then it must be in array format, with LMs as rows, dimensions as columns and specimens as slices')
  }


  if (ShapeGPA==TRUE){
    if (dim(DiscriminationData)[3]!=length(GroupMembership)){
      stop('Number of specimens in DiscriminationData does not appear to match the number of speicmens listed in GroupMembership')
    }
  } else {
    if (dim(DiscriminationData)[1]!=length(GroupMembership)){
      stop('Number of specimens in DiscriminationData does not appear to match the number of speicmens listed in GroupMembership')
    }
  }



  GroupMembership <- as.character(GroupMembership)

  MinSamp <- min(table(GroupMembership))

  if (PClim>MinSamp){
    warning('PClim is set to higher than the smallest sample size.')
  }


  #PClim=10
  #DiscriminationData=RatGPA$PCscores
  #DiscriminationData=RatM1data$LMs
  #GroupMembership=Groups; K=10; Equal=TRUE; EqualIter=100
  #SampleSize=NA; Verbose=TRUE
  #ShapeGPA=TRUE; Sliding=NA; SizeShape=FALSE; PCA=FALSE


  chr2nu <- function(X){
    as.numeric(as.character(X))
  }

  ParOutput <- function(PreviousResults, ResultList){

    NewResults <- PreviousResults
    if (length(PreviousResults)==1){
      NewResults$CVP <- rbind(PreviousResults$CVP, ResultList$CVP)
    } else {
      NewResults$CVP <- rbind(PreviousResults$CVP, ResultList$CVP)
      NewResults$CV.IDs <- abind::abind(PreviousResults$CV.IDs, ResultList$CV.IDs, along = 3)
      NewResults$True.IDs <- NewResults$True.IDs
      #Needs adding/fixing
      #NewResults$PCLoadings <- abind::abind(PreviousResults$PCLoadings, ResultList$PCLoadings, along = 3)

    }

    return(NewResults)
  }


  BalancedGrps <- function(GroupMembership, GroupSize=SampleSize){
    #DiscriminationData=PairwiseShapePCAmat; GroupMembership=chr(Groups[GrpPos])

    if (is.na(GroupSize)){
      minSampSize <- min(table(as.character(GroupMembership)))
    } else {
      minSampSize <- GroupSize
    }

    sampleindex <- 1:length(GroupMembership)
    RandSamp <- stats::aggregate(x = sampleindex, by=list(factor(GroupMembership)), sample, size=minSampSize)
    Index <- c(t(RandSamp$x))
    GroupMem <- c(t(stats::aggregate(RandSamp$Group.1, list(RandSamp$Group.1), rep, minSampSize)$x))
    return(list(IndexedLocations=Index, Newfactors=GroupMem))

  }


  ParEqualIterStepwise <- function(DiscriminationData, GrpMem, ShapeGPA, Sliding, SlidingLMindex, ParVerbose=Verbose, SizeShape, PClim, SampleSize){
    #DiscriminationData=DiscriminationData; GrpMem=GroupMembership
    #GrpMem=Groups
    BalancingGrps <- BalancedGrps(GrpMem, SampleSize)

    if (ShapeGPA==TRUE){
      BalDataShape <- DiscriminationData[,,BalancingGrps$IndexedLocations]
      invisible(utils::capture.output(BalRes <- Morpho::procSym(BalDataShape, sizeshape = SizeShape, outlines = Sliding, SMvector = SlidingLMindex)))
      BalData <- BalRes$PCscores
    } else {
      BalRes <- stats::prcomp(x = DiscriminationData[BalancingGrps$IndexedLocations,])
      BalData <- BalRes$x
    }

    CVPres <- rep(NA, PClim-1)
    CVPbyGrp <- matrix(NA, nrow = length(unique(BalancingGrps$Newfactors)), ncol = PClim-1, dimnames = list(unique(BalancingGrps$Newfactors), NULL))
    CVIDmat <- matrix(NA, nrow = length(BalancingGrps$Newfactors), ncol = PClim-1, dimnames = list(BalancingGrps$Newfactors, NULL))

    #consider getting rid of this TrueIDmat
    TrueIDmat <-  matrix(NA, nrow = length(BalancingGrps$Newfactors), ncol = PClim-1, dimnames = list(BalancingGrps$Newfactors, NULL))
    for(i in 2:PClim){
      LDAres <- MASS::lda(x = BalData[,1:i], grouping=BalancingGrps$Newfactors, CV=TRUE)
      CVPres[i-1] <- sum(LDAres$class==BalancingGrps$Newfactors)/length(BalancingGrps$Newfactors)
      CVPbyGrp[,i-1] <- diag(table(LDAres$class, BalancingGrps$Newfactors))/table(BalancingGrps$Newfactors)

      CVIDmat[,i-1] <- as.character(LDAres$class)
      TrueIDmat[,i-1] <- as.character(BalancingGrps$Newfactors)
    }




    if (ParVerbose==TRUE){
      Res <- list(CVP=CVPres, CV.IDs=CVIDmat, True.IDs=TrueIDmat)
      return(Res)
    } else {
      return(list(CVP=CVPres))
    }
  }


  cores <- parallel::detectCores()
  clust <- parallel::makeCluster(cores[1]-1)
  doParallel::registerDoParallel(clust)

  a <- 1
  ParResults <- foreach::foreach(a = 1:EqualIter, .combine = ParOutput) %dopar% {
    ParEqualIterStepwise(DiscriminationData, GrpMem=GroupMembership, ShapeGPA=ShapeGPA, Sliding=Sliding, SlidingLMindex=SlidingLMindex, ParVerbose=Verbose, PClim=PClim, SizeShape = SizeShape, SampleSize = SampleSize)
  }

  #dim(ParResults$CVP)

  parallel::stopCluster(clust)




  if (PlotResults==TRUE){

    graphics::plot(y = colMeans(ParResults$CVP), x = 2:PClim, type = 'n', ylim = c(10,105), ylab = 'CCV %', xlab = 'PCs')
    graphics::abline(h = seq(from = 20, to = 100, by = 10), v = seq(from = 2, to = PClim, by =2), lty = '1919')

    CCVRange <- apply(ParResults$CVP, MARGIN = 2, FUN = stats::quantile, probs = c(.05, .95))
    graphics::polygon(x = c(2:(PClim), (PClim):2),
                      y = c(CCVRange[1,], CCVRange[2,(PClim-1):1])*100,
                      col = transpar('darkblue', alpha = 25),
                      border = NA)

    graphics::lines(y = colMeans(ParResults$CVP*100), x = 2:PClim, col='darkblue', lwd=3)

    #graphics::legend('bottomright', legend = c('Weighted', 'Unweighted'), col = c('darkblue', 'lightblue'), lty=1, lwd=3, bty = 'o')
  }

  if (Verbose==TRUE){
    names(ParResults) <- c('Iteration.Summaries', 'Verbose.Output.Iteration.CV.IDs')
    return(ParResults)
  } else {
    return(ParResults[[1]])
  }

}



