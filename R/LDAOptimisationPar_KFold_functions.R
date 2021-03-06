




#' LDA leave many out cross-validation with equal sample size using parallel processing
#'
#' This function takes a matrix of variables of known group membership and returns
#' the results of a K-fold correct cross validation identification for each
#' specimen to provide a correct cross-validation percentage. This function is for
#' application to a set number of variables (or PCs), see the function
#' \code{LDACVManyStepwisePar} for stepwise testing of each set of consecutive PCs.
#'
#' The function is primarily for use with resampling unequal groups to equal sample
#' size a set number of times. This follows Evin et al 2013, but the addition of
#' K-fold resampling. This process is carried out with parallel processing.
#' If shape data is used the function offers the option to carry out a new GPA and
#' subsequent PCA for each resampling exercise (set ShapeGPA=TRUE). If raw data
#' is used the function carries out a new PCA with each resampling
#' exercise. In both these instances where a fresh PCA is carried
#' out the function will call on the input value in PClim to determine the number
#' of PCs to use.
#'
#' @param KFold value default set to 5
#' @param TestTraining default set to FALSE, indicates whether an ANOVA and MANOVA should be carried at each step and iteration.
#' @inheritParams LDACVStepwisePar
#' @return Returns a matrix of the leave-one-out classifications for all the specimens along with their known classification.
#'
#' @section Citations:
#'
#' Evin, A., Cucchi, T., Cardini, A., Vidarsdottir, U.S., Larson, G. and Dobney, K., 2013. The long and winding road: identifying pig domestication through molar size and shape. Journal of Archaeological Science, 40(1), pp.735-743.
#'
#' @author Ardern Hulme-Beaman
#' @import MASS
#' @import Morpho
#' @import doParallel
#' @import parallel
#' @import foreach
#' @export


LDACVManyPar <- function(DiscriminationData, GroupMembership, EqualIter=100, KFold=5, SampleSize=NA, Verbose=FALSE, ShapeGPA=FALSE, Sliding=NULL, SlidingLMindex=NULL, Bending=TRUE, SizeShape=FALSE, PClim=10, TestTraining=FALSE){
  #DiscriminationData=Scores
  #DiscriminationData=Shape
  #GroupMembership=Groups
  #EqualIter=10
  #PClim=3
  #ShapeGPA=TRUE; Verbose=FALSE; Sliding=NA;
  #SizeShape=FALSE; PCA=FALSE
  #KFold=2

  #DiscriminationData = BlackRatGPA$rotated
  #GroupMembership=Groups
  #EqualIter=1000
  #PClim=3; ShapeGPA = TRUE; KFold = 5


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
    for (i in 1:length(ResultList)){
      NewResults[[i]] <- rbind(PreviousResults[[i]], ResultList[[i]])
    }

    return(NewResults)
  }


  BalancedGrps <- function(GroupMembership, GroupSize=SampleSize){
    #DiscriminationData=PairwiseShapePCAmat
    #GroupMembership=chr(Groups[GrpPos])

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

  Array2Mat <- function(Array){
    Matrix <- matrix(NA, nrow = dim(Array)[3], ncol = length(c(t(Array[,,1]))))
    for (i in 1:dim(Array)[3]){
      #i <- 1
      Matrix[i,] <- c(t(Array[,,i]))
    }
    return(Matrix)
  }

  ParEqualIterKFold <- function(DiscriminationData, GrpMem, ShapeGPA, Sliding, SlidingLMindex, Bending, SizeShape, PClim, KFold, SampleSize, TestTraining){
    #DiscriminationData=DiscriminationData; GrpMem=GroupMembership; ParTieBreaker='Report'; ParVerbose=FALSE
    #GrpMem=Groups; PClim=3; SampleSize=NA
    #DiscriminationData; GrpMem=GroupMembership; ShapeGPA=ShapeGPA; Sliding=Sliding; PClim=PClim; SizeShape = SizeShape
    #SlidingLMindex=NULL
    BalancingGrps <- BalancedGrps(GrpMem, SampleSize)
    SampSizeFoo <- min(table(BalancingGrps$Newfactors))


    #This object below for the folding exercise is the lenght of the minimum sample size
    Folding <- rep(1:KFold, ceiling(SampSizeFoo/KFold))[1:SampSizeFoo]

    #As we're cbinding the new factors with the folding, the folding will be
    #duplicated equally until it matches the length of the factors
    #these then ensures that the now balanced sample sizes are subjected to the
    #K-fold proceedure equally
    NewFactFoldIndex <- cbind(BalancingGrps$Newfactors, Folding)

    CVPres <- rep(NA, KFold)
    GrpFoldRes <- matrix(NA, nrow = KFold, ncol = length(unique(NewFactFoldIndex[,1])), dimnames = list(paste('F', 1:KFold, sep=''), unique(NewFactFoldIndex[,1])))
    #CVIDmat <- matrix(NA, nrow = length(BalancingGrps$Newfactors), ncol = PClim-1, dimnames = list(BalancingGrps$Newfactors, NULL))
    if (TestTraining==TRUE){
      MANOVAmat <- matrix(NA, nrow = KFold, ncol = 6, dimnames = list(paste('F', 1:KFold, sep=''), 1:6))
      ANOVAmat <- matrix(NA, nrow = KFold, ncol = 5, dimnames = list(paste('F', 1:KFold, sep=''), 1:5))
    }


    if (ShapeGPA==TRUE){
      BalDataShape <- DiscriminationData[,,BalancingGrps$IndexedLocations]
    } else {
      BalData <- DiscriminationData[BalancingGrps$IndexedLocations,]
    }

    for (Kf in 1:KFold){
      #Kf <- 1
      FoldPos <- which(NewFactFoldIndex[,2]==Kf)

      #Figure out sliding issue
      # Note this logical statement starts in different places for the shape data versus non shape data
      # this is because the shape data needs two steps of transformation and each needs to take into
      # consideration the k folding
      if (ShapeGPA==TRUE){
        invisible(utils::capture.output(BalData <- Morpho::procSym(BalDataShape[,,-FoldPos], sizeshape = SizeShape, outlines = Sliding, SMvector = SlidingLMindex, bending = Bending)))
        #Spec=29
        #plot(BalDataShape[,,Spec], type='n', asp=1); text(BalDataShape[,1,Spec], BalDataShape[,2,Spec], labels=1:dim(BalDataShape[,,1])[1], cex=.5)
        #Sliding=list(c(1:21), c(21:50), c(50:63), c(63:74), c(74:105,1)); SlidingLMindex=c(1:105)[-c(1,21,50,63,74)]
        BalPCA <- stats::prcomp(Array2Mat(BalData$orpdata))


        if (!is.null(Sliding)){
          PreAlignBal <- suppressMessages(Morpho::align2procSym(BalData, BalDataShape[,,FoldPos]))
          BalTestSlid <- array(NA, dim = dim(PreAlignBal))
          for (kspec in 1:length(FoldPos)){
            invisible(utils::capture.output(BalTestSlid[,,kspec] <- Morpho::relaxLM(reference=BalData$mshape, lm=PreAlignBal[,,kspec], outlines = Sliding, SMvector = SlidingLMindex, bending=Bending)))
          }

          BalTest <- suppressMessages(Morpho::align2procSym(BalData, BalTestSlid))
        } else {
          BalTest <- suppressMessages(Morpho::align2procSym(BalData, BalDataShape[,,FoldPos]))

        }


        BalTestPCA <- stats::predict(BalPCA, newdata = Array2Mat(BalTest))
      } else {
        BalPCA <- stats::prcomp(x = BalData[-FoldPos,])$x

        BalTestPCA <- stats::predict(BalPCA, newdata = BalData[FoldPos,])
      }

      if (TestTraining==TRUE){
        MANOVARes <- summary(stats::manova(BalPCA$x[,1:PClim] ~ BalancingGrps$Newfactors[-FoldPos]))
        MANOVAmat[Kf,] <- MANOVARes$stats[1,]
        ANOVARes <- summary(stats::aov(BalPCA$x[,PClim]~BalancingGrps$Newfactors[-FoldPos]))
        ANOVAmat[Kf,] <- do.call(c, c(ANOVARes[[1]][1,]))
      }

      LDAres <- MASS::lda(x = BalPCA$x[,1:PClim], grouping=BalancingGrps$Newfactors[-FoldPos], CV=FALSE)
      LDAmany <- stats::predict(LDAres, newdata = BalTestPCA[,1:PClim])

      CVPres[Kf] <- sum(LDAmany$class==BalancingGrps$Newfactors[FoldPos])/length(LDAmany$class)

      GrpFoldRes[Kf,] <- diag(table(LDAmany$class,BalancingGrps$Newfactors[FoldPos]))/table(BalancingGrps$Newfactors[FoldPos])

    }

    if (TestTraining==TRUE){
      colnames(MANOVAmat) <- colnames(MANOVARes$stats)
      colnames(ANOVAmat) <- colnames(ANOVARes[[1]])
      return(list(GlobalCVP=CVPres, GrpCVP=GrpFoldRes, MANOVA=MANOVAmat, ANOVA=ANOVAmat))
    } else {
      return(list(GlobalCVP=CVPres, GrpCVP=GrpFoldRes))
    }

  }


  cores <- parallel::detectCores()
  clust <- parallel::makeCluster(cores[1]-1)
  doParallel::registerDoParallel(clust)

  a <- 1
  ParResults <- foreach::foreach(a = 1:EqualIter, .combine = ParOutput) %dopar% {
    ParEqualIterKFold(DiscriminationData, GrpMem=GroupMembership, ShapeGPA=ShapeGPA, Sliding=Sliding, SlidingLMindex=SlidingLMindex, Bending=Bending,  PClim=PClim, SizeShape = SizeShape, SampleSize = SampleSize, KFold = KFold, TestTraining = TestTraining)
  }

  parallel::stopCluster(clust)

  #colMeans(ParResults[[1]])

  #ParResults$GlobalCVP <- apply(ParResults$GlobalCVP, MARGIN = 1, FUN = 'median')

  if (Verbose==TRUE){
    names(ParResults) <- c('Iteration.Summaries', 'Verbose.Output')
    graphics::hist(apply(ParResults$Iteration.Summaries, 1, 'median'), main= paste('Correct cross-validation % for', PClim, 'PCs', sep = ' '), xlab = 'CCV%')
    return(ParResults)
  } else {
    graphics::hist(apply(ParResults[[1]], 1, 'median'), main= paste('Correct cross-validation % for', PClim, 'PCs', sep = ' '), xlab = 'CCV%')
    return(ParResults[[1]])
  }

}


#' Stepwise LDA correct cross-validation with equal sample size using parallel processing
#'
#' This function takes a matrix of variables of known group membership and returns
#' the results of a leave K out correct cross validation identification for each
#' specimen to provide a correct cross-validation percentage. This function is for
#' for stepwise testing of each set of consecutive PCs with resampling to equal sample
#' size at each incremental increase, following Baylac and Friess 2005.
#'
#' The function is primarily for use with resampling unequal groups to equal sample
#' size a set number of times. This follows Evin et al 2013, but the addition of
#' K-fold resampling. This process is carried out with parallel processing.
#' If shape data is used the function offers the option to carry out a new GPA and
#' subsequent PCA for each resampling exercise (set ShapeGPA=TRUE). If raw data
#' is used the function offers the option to carry out a new PCA with each resampling
#' exercise. In both these instances where a fresh PCA is carried out the function
#' will call on the input value in PClim to determine the number of PCs to use.
#'
#' @param PlotResults logical (either TRUE or FALSE, default set to TRUE) to determine whether to plot the results of the stepwise discriminant analyses.
#' @inheritParams LDACVManyPar
#' @return Returns a matrix of the leave-one-out classifications for all the specimens along with their known classification.
#'
#' @section Citations:
#'
#' Evin, A., Cucchi, T., Cardini, A., Vidarsdottir, U.S., Larson, G. and Dobney, K., 2013. The long and winding road: identifying pig domestication through molar size and shape. Journal of Archaeological Science, 40(1), pp.735-743.
#'
#' Baylac and Friess, 2005 Fourier descriptors, procrustes superimposition and data dimensionality: an example of cranial shape analysis in modern human populations D.E. Slice (Ed.), Modern Morphometrics in Physical Anthropology, Springer-Verlag, New York (2005), pp. 145-165
#'
#' @author Ardern Hulme-Beaman
#' @import MASS
#' @import Morpho
#' @import doParallel
#' @import parallel
#' @import foreach
#' @import abind
#' @export


LDACVManyStepwisePar <- function(DiscriminationData, GroupMembership, EqualIter=100, KFold=5, SampleSize=NA, Verbose=FALSE, ShapeGPA=FALSE, Sliding=NULL, SlidingLMindex=NULL, Bending=TRUE, SizeShape=FALSE, PClim=10, PlotResults=TRUE, TestTraining=FALSE){

  #DiscriminationData=BlackRatM1data$LMArray; GroupMembership = SpeciesGrps; PClim = 10; ShapeGPA = TRUE

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
    for (i in 1:length(ResultList)){
      NewResults[[i]] <- abind::abind(PreviousResults[[i]], ResultList[[i]], along = 1)
      #if (length(dim(PreviousResults[[i]]))==2){
      #  NewResults[[i]] <- rbind(PreviousResults[[i]], ResultList[[i]])
      #} else (length(dim(PreviousResults[[i]]))==3){
      #  NewResults[[i]] <- abind(PreviousResults[[i]], ResultList[[i]], along = 1)
      #}
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


  Array2Mat <- function(Array){
    Matrix <- matrix(NA, nrow = dim(Array)[3], ncol = length(c(t(Array[,,1]))))
    for (i in 1:dim(Array)[3]){
      #i <- 1
      Matrix[i,] <- c(t(Array[,,i]))
    }
    return(Matrix)
  }


  ParEqualIterStepwiseKFold <- function(DiscriminationData, GrpMem, ShapeGPA, Sliding, SlidingLMindex, Bending, KFold, ParVerbose=Verbose, SizeShape, PClim, TestTraining){
    #DiscriminationData=DiscriminationData; GrpMem=GroupMembership
    #GrpMem=Groups
    BalancingGrps <- BalancedGrps(GrpMem, SampleSize)
    SampSizeFoo <- min(table(BalancingGrps$Newfactors))


    #This object below for the folding exercise is the lenght of the minimum sample size
    Folding <- rep(1:KFold, ceiling(SampSizeFoo/KFold))[1:SampSizeFoo]

    #As we're cbinding the new factors with the folding, the folding will be
    #duplicated equally until it matches the length of the factors
    #these then ensures that the now balanced sample sizes are subjected to the
    #K-fold proceedure equally
    NewFactFoldIndex <- cbind(BalancingGrps$Newfactors, Folding)


    CVPres <- matrix(NA, nrow = KFold, ncol = PClim-1, dimnames = list(paste('F', 1:KFold, sep=''), paste('PC1-', 1:(PClim-1), sep='')))
    GrpFoldRes <- array(NA, dim = c(KFold, length(unique(NewFactFoldIndex[,1])), PClim-1), dimnames = list(paste('F', 1:KFold, sep=''), unique(NewFactFoldIndex[,1]), paste('PC1-', 1:(PClim-1), sep='')))
    #CVIDmat <- matrix(NA, nrow = length(BalancingGrps$Newfactors), ncol = PClim-1, dimnames = list(BalancingGrps$Newfactors, NULL))
    if (TestTraining==TRUE){
      MANOVAmat <- array(NA, dim = c(KFold, 6, PClim-1), dimnames = list(paste('F', 1:KFold, sep=''), 1:6, paste('PC1-', 1:(PClim-1), sep='')))
      ANOVAmat <- array(NA, dim = c(KFold, 5, PClim-1), dimnames = list(paste('F', 1:KFold, sep=''), 1:5, paste('PC', 1:(PClim-1), sep='')))
    }


    for(Kf in 1:KFold){
      #Kf=1

      #i=2
      if (ShapeGPA==TRUE){
        BalDataShape <- DiscriminationData[,,BalancingGrps$IndexedLocations]
      } else {
        BalData <- DiscriminationData[BalancingGrps$IndexedLocations,]
      }

      FoldPos <- which(NewFactFoldIndex[,2]==Kf)

      #Figure out sliding issue
      # Note this logical statement starts in different places for the shape data versus non shape data
      # this is because the shape data needs two steps of transformation and each needs to take into
      # consideration the k folding
      if (ShapeGPA==TRUE){
        invisible(utils::capture.output(BalData <- Morpho::procSym(BalDataShape[,,-FoldPos], sizeshape = SizeShape, outlines = Sliding, SMvector = SlidingLMindex, bending = Bending)))
        BalPCA <- stats::prcomp(Array2Mat(BalData$orpdata))

        #invisible(utils::capture.output(BalData2 <- Morpho::procSym(BalDataShape[,,-FoldPos], sizeshape = SizeShape, outlines = NULL, SMvector = NULL)))
        #BalPCA2 <- stats::prcomp(Array2Mat(BalData2$orpdata))


        if (!is.null(Sliding)){
          PreAlignBal <- suppressMessages(Morpho::align2procSym(BalData, BalDataShape[,,FoldPos]))
          BalTestSlid <- array(NA, dim = dim(PreAlignBal))
          for (kspec in 1:length(FoldPos)){
            invisible(utils::capture.output(BalTestSlid[,,kspec] <- Morpho::relaxLM(reference=BalData$mshape, lm=PreAlignBal[,,kspec], outlines = Sliding, SMvector = SlidingLMindex, bending=Bending)))
          }

          BalTest <- suppressMessages(Morpho::align2procSym(BalData, BalTestSlid))
        } else {
          BalTest <- suppressMessages(Morpho::align2procSym(BalData, BalDataShape[,,FoldPos]))

        }
        #dist(rbind(c(t(BalTest[,,1])),c(t(BalData$mshape))))
        #dist(rbind(c(t(BalTest2[,,1])),c(t(BalData$mshape))))
        #dist(rbind(c(t(BalDataShape[,,FoldPos[1]])),c(t(BalData$mshape))))
        #plot(BalTest[,,1], asp=1)
        #points(BalDataShape[,1,FoldPos[1]], BalDataShape[,2,FoldPos[1]], col='red')
        #points(PreAlignBal[,1,1], PreAlignBal[,2,1], col='red', pch=16)
        #points(BalData$mshape[,1], BalData$mshape[,2], col='blue', pch=16)



        BalTestPCA <- stats::predict(BalPCA, newdata = Array2Mat(BalTest))
        #BalTestPCA2 <- stats::predict(BalPCA, newdata = Array2Mat(BalTest2))
        #plot(BalPCA$x)
        #points(BalTestPCA[,1], BalTestPCA[,2], col='red')
        #points(BalPCA2$x[,1], BalPCA2$x[,2], col='black', pch=16)
        #plot(BalPCA2$x)
        #points(BalTestPCA2[,1], BalTestPCA2[,2], col='red', pch=16)
        #points(BalTestPCA2[1,1], BalTestPCA2[1,2], col='yellow', pch=16)
        #points(BalTestPCA[1,1], BalTestPCA[1,2], col='yellow')
      } else {
        BalPCA <- stats::prcomp(x = BalData[-FoldPos,])$x

        BalTestPCA <- stats::predict(BalPCA, newdata = BalData[FoldPos,])
      }


      for(i in 2:PClim){

        if (TestTraining==TRUE){
          MANOVARes <- summary(stats::manova(BalPCA$x[,1:i] ~ BalancingGrps$Newfactors[-FoldPos]))
          MANOVAmat[Kf,,i-1] <- MANOVARes$stats[1,]
          ANOVARes <- summary(stats::aov(BalPCA$x[,i]~BalancingGrps$Newfactors[-FoldPos]))
          ANOVAmat[Kf,,i-1] <- do.call(c, c(ANOVARes[[1]][1,]))
        }

        LDAres <- MASS::lda(x = BalPCA$x[,1:i], grouping=BalancingGrps$Newfactors[-FoldPos], CV=FALSE)
        LDAmany <- stats::predict(LDAres, newdata = BalTestPCA[,1:i])


        CVPres[Kf,i-1] <- sum(LDAmany$class==BalancingGrps$Newfactors[FoldPos])/length(LDAmany$class)

        GrpFoldRes[Kf,,i-1] <- diag(table(LDAmany$class,BalancingGrps$Newfactors[FoldPos]))/table(BalancingGrps$Newfactors[FoldPos])


      }

    }




    if (ParVerbose==TRUE){
      if (TestTraining==TRUE){
        colnames(MANOVAmat) <- colnames(MANOVARes$stats)
        colnames(ANOVAmat) <- colnames(ANOVARes[[1]])
        Res <- list(CVP=CVPres, CVPbyGrp=GrpFoldRes, MANOVA=MANOVAmat, ANOVA=ANOVAmat)
      } else {
        Res <- list(CVP=CVPres, CVPbyGrp=GrpFoldRes)
      }
      return(Res)
    } else {
      return(list(CVP=CVPres))
    }
  }


  cores <- parallel::detectCores()
  clust <- parallel::makeCluster(cores[1]-1)
  doParallel::registerDoParallel(clust)

  #NEED TIME TRIAL
  a <- 1
  ParResults <- foreach::foreach(a = 1:EqualIter, .combine = ParOutput, .packages = 'abind') %dopar% {
    ParEqualIterStepwiseKFold(DiscriminationData, GrpMem=GroupMembership, ShapeGPA=ShapeGPA, Sliding=Sliding, SlidingLMindex=SlidingLMindex, Bending=Bending, KFold=KFold, ParVerbose=Verbose, PClim=PClim, SizeShape = SizeShape, TestTraining=TestTraining)
  }



  parallel::stopCluster(clust)




  if (PlotResults==TRUE){


    FoldingSummary <- function(X, KFold){
      colMeans(t(matrix(X, nrow = length(X)/KFold, ncol = KFold, byrow = TRUE)))#, 2, median))
    }

    FoldCVP <- apply(ParResults$CVP, 2, FoldingSummary, KFold)

    graphics::plot(y = colMeans(FoldCVP), x = 2:PClim, type = 'n', ylim = c(10,105), ylab = 'CCV %', xlab = 'PCs')
    graphics::abline(h = seq(from = 20, to = 100, by = 10), v = seq(from = 2, to = PClim, by =2), lty = '1919')

    CCVRange <- apply(FoldCVP, MARGIN = 2, FUN = stats::quantile, probs = c(.05, .95))
    graphics::polygon(x = c(2:(PClim), (PClim):2),
                      y = c(CCVRange[1,], CCVRange[2,(PClim-1):1])*100,
                      col = transpar('darkblue', alpha = 25),
                      border = NA)

    graphics::lines(y = colMeans(FoldCVP*100), x = 2:PClim, col='darkblue', lwd=3)

    #graphics::legend('bottomright', legend = c('Weighted', 'Unweighted'), col = c('darkblue', 'lightblue'), lty=1, lwd=3, bty = 'o')
  }

  if (Verbose==TRUE){
    names(ParResults) <- c('Iteration.Summaries', 'Verbose.Output.Iteration.CV.IDs')
    return(ParResults)
  } else {
    return(ParResults[[1]])
  }

}



