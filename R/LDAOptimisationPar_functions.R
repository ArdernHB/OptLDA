




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
#' is used the funciton offers the option to carry out a new PCA with each resampling
#' exercise (set PCA=TRUE). In both these instances where a fresh PCA is carried
#' out the function will call on the input value in PClim to determine the number
#' of PCs to use. If a dataset that does not require additional PCA (e.g. a matrix
#' of PC scores) is examined then both the arguments ShapeGPA and PCA can be set
#' to FALSE and the resampling procedure will be carried out just on matrix provided.
#'
#' @param Data A matrix of numeric data or an array for shape data. If shape data is used please set ShapeGPA=TRUE. A shape data array is expected to include landmarks as rows, landmark dimensions (either 2 or 3) as columns and specimens as slices. If a matrix is used, please ensure no columns have erroneously been loaded as factors.
#' @param GroupMembership A vector of group classifications, either as factors or as characters.
#' @param EqualIter an integer determining the number of times the data will be resampled to equal sample size. Default is set to 100.
#' @param SampleSize the sample size to be used. If set to NA the resampling exercise will automatically resample groups to the sample size of the smallest group.
#' @param Verbose logical (either TRUE or FALSE, default set to FALSE) to determine whether to return the identifications of each round of resampling (Verbose=TRUE) or to just return the correct-cross-validation percentages of each resampling exercise (Verbose=FALSE).
#' @param ShapeGPA logical (either TRUE or FALSE, default set to FALSE) to indicate whether the data is shape data and if so each subset of individuals from the resampling procedure will be processed through a Generalised Procrustes Analysis (GPA) and subsequent Principal Component Analysis (PCA) using the \code{\link[Morpho]{procSym}} function in the package \code{Morpho}.
#' @param Sliding if ShapeGPA is TRUE and the shape data has sliding landmarks then a list should be provided here of landmarks to be slid (see function description of \code{\link[Morpho]{procSym}})
#' @param SizeShape logical (either TRUE or FALSE, default set to FALSE). If ShapeGPA is TRUE and you wish to analyses form (i.e. shape+log size) then this argument will be passed to the \code{\link[Morpho]{procSym}} and the balanced LDA will be applied to form PCs (see function description of \code{\link[Morpho]{procSym}}).
#' @param PCA logical (either TRUE or FALSE, default set to FALSE). If the data is raw data and requires a PCA first then a PCA is carried out using the \code{\link[stats]{prcomp}} function of the \code{stats} package.
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



LDACVPar <- function(Data, GroupMembership, EqualIter=100, SampleSize=NA, Verbose=FALSE, ShapeGPA=FALSE, Sliding=NA, SizeShape=FALSE, PCA=FALSE, PClim=10){
  #PClim=10
  #Data=RatGPA$PCscores
  #Data=RatM1data$LMs
  #GroupMembership=Groups; K=10; Equal=TRUE; EqualIter=100
  #SampleSize=NA; Verbose=TRUE
  #ShapeGPA=TRUE; Sliding=NA; SizeShape=FALSE; PCA=FALSE


  chr2nu <- function(X){
    as.numeric(as.character(X))
  }

  ParOutput <- function(PreviousResults, ResultList){

    NewResults <- PreviousResults
    for (i in 1:length(ResultList)){
      NewResults[[i]] <- rbind(PreviousResults[[i]], ResultList[[i]])
    }

    return(NewResults)
  }


  BalancedGrps <- function(GroupMembership, GroupSize=SampleSize){
    #Data=PairwiseShapePCAmat; GroupMembership=chr(Groups[GrpPos])

    if (is.na(GroupSize)){
      minSampSize <- min(table(GroupMembership))
    } else {
      minSampSize <- GroupSize
    }

    sampleindex <- 1:length(GroupMembership)
    RandSamp <- stats::aggregate(x = sampleindex, by=list(factor(GroupMembership)), sample, size=minSampSize)
    Index <- c(t(RandSamp$x))
    GroupMem <- c(t(stats::aggregate(RandSamp$Group.1, list(RandSamp$Group.1), rep, minSampSize)$x))
    return(list(IndexedLocations=Index, Newfactors=GroupMem))

  }


  ParEqualIter <- function(Data, GrpMem, ShapeGPA, Sliding, PCA, ParVerbose=Verbose, SizeShape, PClim){
    #DistData=PCAmat; GrpMem=GroupMembership; ParK=K; ParTieBreaker='Report'; ParVerbose=FALSE
    #GrpMem=Groups
    BalancingGrps <- BalancedGrps(GrpMem, SampleSize)


    if (ShapeGPA==TRUE){
      BalDataShape <- Data[,,BalancingGrps$IndexedLocations]
      if (is.list(Sliding)){
        BalData <- suppressMessages(Morpho::procSym(BalDataShape, sizeshape = SizeShape, outlines = Sliding)$PCscores)
        LDAres <- MASS::lda(x = BalData[,1:PClim], grouping=BalancingGrps$Newfactors, CV=TRUE)
      } else if (is.na(Sliding)){
        BalData <- suppressMessages(Morpho::procSym(BalDataShape, sizeshape = SizeShape)$PCscores)
        LDAres <- MASS::lda(x = BalData[,1:PClim], grouping=BalancingGrps$Newfactors, CV=TRUE)
      }
    } else if (PCA==TRUE){
      BalData <- stats::prcomp(x = Data[BalancingGrps$IndexedLocations,])
      LDAres <- MASS::lda(x = BalData[,1:PClim], grouping=BalancingGrps$Newfactors, CV=TRUE)
    } else {
      BalData <- Data[BalancingGrps$IndexedLocations,]
      LDAres <- MASS::lda(x = BalData, grouping=BalancingGrps$Newfactors, CV=TRUE)
    }

    CVP <- sum(LDAres$class==BalancingGrps$Newfactors)/length(BalancingGrps$Newfactors)

    if (ParVerbose==TRUE){
      Res <- list(CVP=CVP, IDs=as.character(LDAres$class))
      return(Res)
    } else {
      return(list(CVP))
    }
  }


  cores <- parallel::detectCores()
  clust <- parallel::makeCluster(cores[1]-1)
  doParallel::registerDoParallel(clust)

  a <- 1
  ParResults <- foreach::foreach(a = 1:EqualIter, .combine = ParOutput) %dopar% {
    ParEqualIter(Data, GrpMem=GroupMembership, ShapeGPA=ShapeGPA, Sliding=Sliding, PCA=PCA, ParVerbose=Verbose, PClim=PClim, SizeShape = SizeShape)
  }

  parallel::stopCluster(clust)


  if (Verbose==TRUE){
    names(ParResults) <- c('Iteration.Summaries', 'Verbose.Output')
    return(ParResults)
  } else {
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
#' size a set number of times. This process is carried out with parrallel processing.
#' If shape data is used the function offers the option to carry out a new GPA and
#' subsequent PCA for each resampling exercise (set ShapeGPA=TRUE). If raw data
#' is used the funciton offers the option to carry out a new PCA with each resampling
#' exercise (set PCA=TRUE). In both these instances where a fresh PCA is carried
#' out the function will call on the input value in PClim to determine the number
#' of PCs to use. If a dataset that does not require additional PCA (e.g. a matrix
#' of PC scores) is examined then both the arguments ShapeGPA and PCA can be set
#' to FALSE and the resampling procedure will be carried out just on matrix provided.
#'
#' @param PlotResults logical (either TRUE or FALSE, default set to TRUE) to determine whether to plot the results of the stepwise discriminant analyses.
#' @inheritParams LDACVPar
#' @return Returns a matrix of the leave-one-out classifications for all the specimens along with their known classification.
#'
#'
#' @keywords shape distances
#' @keywords Geometric morphometric distances
#' @author Ardern Hulme-Beaman
#' @import MASS
#' @import Morpho
#' @import doParallel
#' @import parallel
#' @import foreach
#' @export


LDACVStepwisePar <- function(Data, GroupMembership, EqualIter=100, SampleSize=NA, Verbose=FALSE, ShapeGPA=FALSE, Sliding=NA, SizeShape=FALSE, PClim=10, PlotResults=TRUE){

  #DistMat = VoleDistMat; GroupMembership = VoleGrps; Kmax = 10; TieBreaker = 'Remove'; PlotResults = TRUE; EqualIter = 100
  #Kmax=20
  #DistMat=ProcDTableRes; GroupMembership=Groups; K=10; Equal=TRUE; EqualIter=100
  #Weighted=TRUE; TieBreaker='Report'
  #PrintProg=TRUE
  #Verbose=TRUE; Equal=TRUE

  MinSamp <- min(table(as.character(GroupMembership)))

  if (PClim>MinSamp){
    warning('PClim is set to higher than the smallest sample size.')
  }


  #PClim=10
  #Data=RatGPA$PCscores
  #Data=RatM1data$LMs
  #GroupMembership=Groups; K=10; Equal=TRUE; EqualIter=100
  #SampleSize=NA; Verbose=TRUE
  #ShapeGPA=TRUE; Sliding=NA; SizeShape=FALSE; PCA=FALSE


  chr2nu <- function(X){
    as.numeric(as.character(X))
  }

  ParOutput <- function(PreviousResults, ResultList){

    NewResults <- PreviousResults
    for (i in 1:length(ResultList)){
      NewResults[[i]] <- rbind(PreviousResults[[i]], ResultList[[i]])
    }

    return(NewResults)
  }


  BalancedGrps <- function(GroupMembership, GroupSize=SampleSize){
    #Data=PairwiseShapePCAmat; GroupMembership=chr(Groups[GrpPos])

    if (is.na(GroupSize)){
      minSampSize <- min(table(GroupMembership))
    } else {
      minSampSize <- GroupSize
    }

    sampleindex <- 1:length(GroupMembership)
    RandSamp <- stats::aggregate(x = sampleindex, by=list(factor(GroupMembership)), sample, size=minSampSize)
    Index <- c(t(RandSamp$x))
    GroupMem <- c(t(stats::aggregate(RandSamp$Group.1, list(RandSamp$Group.1), rep, minSampSize)$x))
    return(list(IndexedLocations=Index, Newfactors=GroupMem))

  }


  ParEqualIterStepwise <- function(Data, GrpMem, ShapeGPA, Sliding, ParVerbose=Verbose, SizeShape, PClim){
    #DistData=PCAmat; GrpMem=GroupMembership; ParK=K; ParTieBreaker='Report'; ParVerbose=FALSE
    #GrpMem=Groups
    BalancingGrps <- BalancedGrps(GrpMem, SampleSize)

    CVPres <- rep(NA, PClim-1)
    CVIDmat <- matrix(NA, nrow = length(BalancingGrps$Newfactors), ncol = PClim-1, dimnames = list(BalancingGrps$Newfactors, NULL))
    for(i in 2:PClim){
      #i=2
      if (ShapeGPA==TRUE){
        BalDataShape <- Data[,,BalancingGrps$IndexedLocations]
        if (is.list(Sliding)){
          BalData <- suppressMessages(Morpho::procSym(BalDataShape, sizeshape = SizeShape, outlines = Sliding)$PCscores)
          LDAres <- MASS::lda(x = BalData[,1:i], grouping=BalancingGrps$Newfactors, CV=TRUE)
        } else if (is.na(Sliding)){
          BalData <- suppressMessages(Morpho::procSym(BalDataShape, sizeshape = SizeShape)$PCscores)
          LDAres <- MASS::lda(x = BalData[,1:i], grouping=BalancingGrps$Newfactors, CV=TRUE)
        }
      } else {
        BalData <- Data[BalancingGrps$IndexedLocations,]
        LDAres <- MASS::lda(x = BalData[,1:i], grouping=BalancingGrps$Newfactors, CV=TRUE)
      }
      CVPres[i-1] <- sum(LDAres$class==BalancingGrps$Newfactors)/length(BalancingGrps$Newfactors)
      CVIDmat[,i-1] <- as.character(LDAres$class)


    }




    if (ParVerbose==TRUE){
      Res <- list(CVP=CVPres, CV.IDs=CVIDmat)
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
    ParEqualIterStepwise(Data, GrpMem=GroupMembership, ShapeGPA=ShapeGPA, Sliding=Sliding, ParVerbose=Verbose, PClim=PClim, SizeShape = SizeShape)
  }



  parallel::stopCluster(clust)




  if (PlotResults==TRUE){

    graphics::plot(y = colMeans(ParResults$CVP), x = 2:PClim, type = 'n', ylim = c(10,105), ylab = 'CCV %', xlab = 'K')
    graphics::abline(h = seq(from = 20, to = 100, by = 10), v = seq(from = 2, to = PClim, by =2), lty = '1919')

    CCVRange <- apply(ParResults$CVP, MARGIN = 2, FUN = stats::quantile, probs = c(.05, .95))
    graphics::polygon(x = c(1:PClim, PClim:1),
                      y = c(CCVRange[1,], CCVRange[2,PClim:1])*100,
                      col = transpar('darkblue', alpha = 25),
                      border = NA)

    graphics::lines(y = colMeans(ParResults$CVP*100), x = 1:PClim, col='darkblue', lwd=3)

    graphics::legend('bottomright', legend = c('Weighted', 'Unweighted'), col = c('darkblue', 'lightblue'), lty=1, lwd=3, bty = 'o')
  }

  if (Verbose==TRUE){
    names(ParResults) <- c('Iteration.Summaries', 'Verbose.Output.Iteration.CV.IDs')
    return(ParResults)
  } else {
    return(ParResults[[1]])
  }

}


